#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Predict runtime, CPU and memory for a dispel4py workflow at a given scale.

An analytical pipeline model, not a task scheduler. dispel4py runs every PE
instance concurrently as its own process and streams items between them, so the
workflow does not finish when a critical path finishes -- it finishes when the
slowest *stage* has drained the stream. That makes the model:

    items_p  = items (stream length) x selectivity_p
                                       how many items reach this PE at all
    stage_p  = ceil(items_p / k_p)     how long this PE stays busy: its
               x secs_per_item_p       *busiest* instance's share of the items
    makespan = max_p stage_p + fill    the slowest PE, plus the part of the
                                       pipeline it cannot overlap

**Why the slowest PE and not the sum.** Every PE instance is its own process and
they all run at once, so stages do not queue behind each other. Once the pipeline
is full the workflow is busy exactly as long as its slowest stage needs to drain
the stream; every other stage finishes early and waits.

**Why the ceiling.** Items are whole things, handed out one at a time. Five items
across four instances go 2/1/1/1, so the stage lasts as long as the instance
holding two, not the 1.25 items an even split suggests. In the traced runs one
instance really did process zero items while another took the extra.

**What fill is.** The pipeline starts empty: before the bottleneck can begin, the
first item has to reach it, and after the bottleneck finishes, the last item must
still traverse whatever follows. That is one item's service time at every stage
*except the bottleneck itself*.

CPU is modelled separately rather than derived from runtime: a PE blocked on an
LLM call burns ~10% of a core for ten seconds, while a trivial PE burns ~100%
for microseconds. Using runtime as a proxy would be wrong by an order of
magnitude in both directions.

Memory is a per-process interpreter baseline plus what each PE adds on top of
it. The instrumentation records max observed process RSS and warns never to sum
it across PEs or instances, so the total here sums baselines and attributable
footprints, never raw RSS.
"""

import collections
import json
import logging
import math
import pathlib
from typing import Any, Dict, Optional

import networkx as nx

from . import resource_stats

logger = logging.getLogger(__name__)


def instance_counts(instance: pathlib.Path) -> Dict[str, int]:
    """How many instances of each PE a WfFormat instance holds."""
    spec = json.loads(pathlib.Path(instance).read_text())["workflow"]["specification"]
    return dict(collections.Counter(
        resource_stats.pe_type(t["id"]) for t in spec["tasks"]))


def pe_graph(instance: pathlib.Path) -> nx.DiGraph:
    """The PE-level dataflow graph behind a WfFormat instance.

    Tasks are per PE instance; collapsing them onto their PE recovers the
    abstract graph dispel4py was built from, which is what decides whether two
    stages run in sequence or side by side.
    """
    spec = json.loads(pathlib.Path(instance).read_text())["workflow"]["specification"]
    graph = nx.DiGraph()
    for task in spec["tasks"]:
        child = resource_stats.pe_type(task["id"])
        graph.add_node(child)
        for parent in task["parents"]:
            source = resource_stats.pe_type(parent)
            if source != child:      # replicas of one PE are not a dependency
                graph.add_edge(source, child)
    return graph


def _longest_before(graph: nx.DiGraph, target: str, cost) -> float:
    """Cost of the slowest path into ``target``, excluding ``target`` itself."""
    best: Dict[str, float] = {}
    for node in nx.topological_sort(graph):
        best[node] = max((best[p] + cost(p) for p in graph.predecessors(node)),
                         default=0.0)
    return best.get(target, 0.0)


def _longest_after(graph: nx.DiGraph, source: str, cost) -> float:
    """Cost of the slowest path out of ``source``, excluding ``source`` itself."""
    best: Dict[str, float] = {}
    for node in reversed(list(nx.topological_sort(graph))):
        best[node] = max((best[s] + cost(s) for s in graph.successors(node)),
                         default=0.0)
    return best.get(source, 0.0)


def fill_drain(graph: Optional[nx.DiGraph], bottleneck: str, cost, pes) -> tuple:
    """Time the bottleneck cannot overlap: filling the pipeline, then draining it.

    With the graph, this is the slowest path *into* the bottleneck plus the
    slowest path *out of* it. Branches that run beside the bottleneck are
    excluded, because their work overlaps the bottleneck rather than queueing
    behind it, and where several branches converge only the slowest is charged.

    Without a graph -- a caller who passed bare instance counts -- every other
    stage is summed instead. That is right only for a chain, and overstates any
    workflow with a branch that bypasses the bottleneck.

    :return: (seconds, how it was derived)
    """
    if graph is not None and bottleneck in graph:
        if not nx.is_directed_acyclic_graph(graph):
            logger.warning("the PE graph has a cycle, so there is no longest "
                           "path; summing every other stage instead")
        else:
            return (_longest_before(graph, bottleneck, cost)
                    + _longest_after(graph, bottleneck, cost)), "critical path"
    if graph is None:
        logger.warning("no dataflow graph given, so every non-bottleneck stage "
                       "is charged to fill/drain. That is exact for a chain and "
                       "pessimistic for a branching workflow -- pass a WfFormat "
                       "instance instead of bare counts to use the real paths.")
    return sum(cost(p) for p in pes if p != bottleneck), "sum of other stages"


def simulate(instance,
             stats,
             items: float = None,
             processes_per_instance: int = 1,
             graph: nx.DiGraph = None) -> Dict[str, Any]:
    """Predict cost for a workflow shape at a stream size.

    :param instance: a WfFormat instance path, or {pe_id: instance count}. Only
        a path carries the dataflow graph; bare counts fall back to a chain
        assumption for the fill/drain term (see `fill_drain`).
    :param graph: the PE-level dataflow graph, when `instance` is bare counts.
    :param stats: `resource_stats.learn` output, or a path to saved JSON.
    :param items: items to push through. Defaults to the widest observed run.
    :param processes_per_instance: 1 under timed_multi/timed_mpi, where every
        instance is its own process. Pass 0 for timed_simple, where all PEs
        share one process and only a single baseline is paid.
    :return: {"makespan_secs", "bottleneck", "cpu", "memory", "pes"}
    """
    if isinstance(stats, (str, pathlib.Path)):
        stats = resource_stats.load(stats)
    if isinstance(instance, (str, pathlib.Path)):
        counts = instance_counts(instance)
        graph = graph if graph is not None else pe_graph(instance)
    else:
        counts = dict(instance)
        if graph is None and stats.get("edges"):
            # the dataflow graph the statistics were learned from
            graph = nx.DiGraph(tuple(e) for e in stats["edges"])
    items = items if items is not None else stats["reference_items"]
    baseline = stats["baseline_rss_bytes"]

    pes: Dict[str, Any] = {}
    unknown = []
    for pe_id, k in sorted(counts.items()):
        pe = stats["pes"].get(pe_id)
        if pe is None:
            unknown.append(pe_id)
            continue
        pe_items = items * pe["selectivity"]  # the fraction of the stream that reaches this PE
        work = pe_items * pe["secs_per_item"] # the total time spent on this PE across all instances
        busiest = math.ceil(pe_items / k) if k else 0 # gets the PE that runs the most items (when not equal)
        stage = busiest * pe["secs_per_item"]
        # the same stage if every item ran as fast, or as slow, as the fastest
        # and slowest runs on record
        stage_lo = busiest * pe.get("secs_per_item_min", pe["secs_per_item"])
        stage_hi = busiest * pe.get("secs_per_item_max", pe["secs_per_item"])
        attributable = pe["rss_attributable_bytes"] or 0.0 # amount of memory this PE is responsible for, above the baseline
        pes[pe_id] = {
            "instances": k, # how many instances of this PE are running
            "items": pe_items,
            "work_secs": work,
            "stage_secs": stage,
            "stage_secs_range": (stage_lo, stage_hi),
            "items_on_busiest": busiest,
            "even_split_secs": work / k if k else 0.0,
            "cpu_secs": pe_items * pe["cpu_secs_per_item"],
            "cpu_percent": pe["cpu_percent"],
            "peak_cores": k * pe["cpu_percent"] / 100.0,
            "rss_per_instance_bytes": baseline * processes_per_instance + attributable,
            "rss_bytes": k * (baseline * processes_per_instance + attributable),
            "rss_attributable_bytes": attributable,
            "measured_memory": pe["rss_attributable_bytes"] is not None,
        }

    if unknown:
        logger.warning("no statistics for %s; excluded from the prediction. "
                       "The traces the stats came from did not contain these "
                       "PEs -- check that the instance and the traces are the "
                       "same workflow.", ", ".join(unknown))
    if not pes:
        raise SystemExit("none of the instance's PEs appear in the statistics")

    bottleneck = max(pes, key=lambda p: pes[p]["stage_secs"])
    # The pipeline starts empty and ends draining; the bottleneck's own per-item
    # cost is already inside stage_secs, so it is never charged again here.
    cost = lambda p: stats["pes"][p]["secs_per_item"] if p in stats["pes"] else 0.0
    fill, fill_basis = fill_drain(graph, bottleneck, cost, pes) 
    makespan = pes[bottleneck]["stage_secs"] + fill
    makespan_range = (pes[bottleneck]["stage_secs_range"][0] + fill,
                      pes[bottleneck]["stage_secs_range"][1] + fill)

    cpu_secs = sum(p["cpu_secs"] for p in pes.values())
    memory = sum(p["rss_bytes"] for p in pes.values())
    if processes_per_instance == 0:
        memory = baseline + sum(p["rss_attributable_bytes"] for p in pes.values())

    return {
        "items": items,
        "makespan_secs": makespan,
        "makespan_range_secs": makespan_range,
        "fill_secs": fill,
        "fill_basis": fill_basis,
        "bottleneck": bottleneck,
        "unknown_pes": unknown,
        "cpu": {
            "core_seconds": cpu_secs,
            "mean_cores": cpu_secs / makespan if makespan else 0.0,
            "peak_cores": sum(p["peak_cores"] for p in pes.values()),
            "utilization": (cpu_secs / (makespan * sum(p["instances"] for p in pes.values()))
                            if makespan else 0.0),
        },
        "memory": {
            "total_bytes": memory,
            "baseline_bytes": baseline,
            "attributable_bytes": sum(p["rss_attributable_bytes"] * p["instances"]
                                      for p in pes.values()),
            "processes": sum(p["instances"] for p in pes.values()) * processes_per_instance,
        },
        "pes": pes,
    }


def report(result: Dict[str, Any]) -> str:
    """Render a prediction as a table."""
    lines = []
    lines.append(f"{result['items']:.0f} items")
    lines.append("")
    lines.append(f"{'PE':28} {'inst':>4} {'items':>8} {'stage s':>10} "
                 f"{'cpu s':>9} {'cores':>6} {'mem MB':>9}")
    lines.append("-" * 80)
    for pe_id, pe in sorted(result["pes"].items(),
                            key=lambda kv: -kv[1]["stage_secs"]):
        flag = "" if pe["measured_memory"] else " *"
        lines.append(f"{pe_id:28} {pe['instances']:>4} {pe['items']:>8.0f} "
                     f"{pe['stage_secs']:>10.3f} {pe['cpu_secs']:>9.3f} "
                     f"{pe['peak_cores']:>6.2f} {pe['rss_bytes']/1e6:>9.1f}{flag}")
    if any(not pe["measured_memory"] for pe in result["pes"].values()):
        lines.append("  * memory not attributable from the traces; baseline only")

    cpu, mem = result["cpu"], result["memory"]
    lines.append("")
    lo, hi = result["makespan_range_secs"]
    lines.append(f"makespan      {result['makespan_secs']:>12.3f} s   "
                 f"[{lo:.1f} - {hi:.1f} across observed runs]")
    lines.append(f"{'':14}{'':>12}     "
                 f"bottleneck {result['bottleneck']} "
                 f"({result['pes'][result['bottleneck']]['stage_secs']:.3f} s)"
                 f" + fill/drain {result['fill_secs']:.3f} s "
                 f"({result['fill_basis']})")
    lines.append(f"cpu           {cpu['core_seconds']:>12.3f} core-s  "
                 f"mean {cpu['mean_cores']:.2f} cores, "
                 f"peak {cpu['peak_cores']:.2f}, "
                 f"util {cpu['utilization']*100:.1f}%")
    lines.append(f"memory        {mem['total_bytes']/1e6:>12.1f} MB      "
                 f"({mem['processes']} processes x "
                 f"{mem['baseline_bytes']/1e6:.0f} MB baseline + "
                 f"{mem['attributable_bytes']/1e6:.1f} MB attributable)")
    return "\n".join(lines)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("instance", type=pathlib.Path,
                        help="a WfFormat instance (real or synthetic)")
    parser.add_argument("-s", "--stats", required=True, type=pathlib.Path,
                        help="resource statistics JSON (see resource_stats)")
    parser.add_argument("-i", "--items", type=float, default=None,
                        help="items to push through (default: the widest run)")
    parser.add_argument("--shared-process", action="store_true",
                        help="timed_simple: every PE shares one process")
    args = parser.parse_args()

    logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")
    result = simulate(args.instance, args.stats, args.items,
                      processes_per_instance=0 if args.shared_process else 1)
    print(report(result))


if __name__ == "__main__":
    main()

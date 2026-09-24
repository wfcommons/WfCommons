#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Learn per-PE time, CPU and memory costs from dispel4py monitoring traces.

What the monitoring CSVs record, and how it is modelled here:

time
    ``total_secs`` in the per-run summary is the sum over ranks: the aggregate
    service time a PE spent on the whole stream, independent of how many
    instances it ran in. That makes it the quantity to scale by input size and
    divide by instance count, which is what the simulator does.

cpu
    ``total_cpu_secs`` is CPU time, not wall time, and the two diverge sharply:
    an LLM-calling PE blocks on the network at ~10% CPU while a trivial PE runs
    at ~100%. CPU is therefore learned on its own, never derived from runtime.

memory
    ``rss_*_bytes`` is the whole process's RSS, which is mostly interpreter and
    libraries. The floor across all instances is taken as a per-process
    baseline; what a PE adds on top of it is its attributable footprint. Under
    ``timed_simple`` every PE shares one process, so per-PE RSS there is the
    same number repeated -- only the baseline is trustworthy from such a run.
"""

import csv
import json
import logging
import pathlib
import re
import statistics
from typing import Any, Dict, List

import networkx as nx

logger = logging.getLogger(__name__)

_RUN_ID_RE = re.compile(r"_run(?P<run_id>[^.]+)\.csv$")


# A PE rewritten for parallel execution is the same PE: the parallel mapping's
# ParallelLLMSensorAgentPE4 and the simple mapping's LLMSensorAgentPE4 do the
# same work, so their measurements pool. Without this they stay separate, and a
# workflow grown from the simple graph never sees what the parallel runs
# measured -- including their far better item counts.
_ALIAS_PREFIXES = ("Parallel",)


def canonical_pe(pe_id: str) -> str:
    """The PE's name with a mapping-specific prefix removed."""
    for prefix in _ALIAS_PREFIXES:
        if pe_id.startswith(prefix) and len(pe_id) > len(prefix):
            return pe_id[len(prefix):]
    return pe_id


def pe_type(task_id: str) -> str:
    """The PE a task id belongs to.

    Converted instances number tasks ``<pe>_ID0000001`` and generated ones
    ``<pe>_00000001``; splitting on the last underscore handles both, and PE
    names containing underscores survive it. The result is canonical, so a task
    from a parallel mapping resolves to the same PE as its simple counterpart.
    """
    return canonical_pe(task_id.rsplit("_", 1)[0])


def _instance_rows(monitoring_dir: pathlib.Path, prefix: str = "monitor") -> List[dict]:
    """Every per-instance row of every run in a monitoring directory.

    All runs are read, not just the newest: repeat runs of the same
    configuration are the only evidence of how much a cost varies, and each run
    is aggregated separately anyway.
    """
    matches = sorted(monitoring_dir.glob(f"{prefix}_instances_run*.csv"))
    if not matches:
        raise FileNotFoundError(
            f"no {prefix}_instances_run*.csv in {monitoring_dir}")
    if len(matches) > 1:
        logger.info("%s holds %d runs; pooling them", monitoring_dir.name, len(matches))
    rows = []
    for match in matches:
        with match.open(newline="", encoding="utf-8") as handle:
            rows.extend(csv.DictReader(handle))
    return rows


def _number(row: dict, column: str, default: float = 0.0) -> float:
    value = row.get(column)
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _abstract_edges(monitoring_dir: pathlib.Path, prefix: str = "monitor") -> list:
    """PE-level edges from a run's abstract shape, with names canonicalised."""
    edges = []
    for path in sorted(monitoring_dir.glob(f"{prefix}_shape_run*.json")):
        shape = json.loads(path.read_text())
        for edge in shape.get("edges", []):
            source, target = canonical_pe(edge["from"]), canonical_pe(edge["to"])
            if source != target:
                edges.append((source, target))
    return edges


def split_feedback(edges) -> tuple:
    """Separate a PE graph into its acyclic part and its feedback edges.

    dispel4py allows loops; WfChef and the critical-path model both need a DAG.
    Rather than dropping the loop, the edges that close it are kept aside, so a
    consumer can see that a loop exists and how much it costs. Every edge inside
    a strongly connected component of more than one node closes some cycle,
    including a PE wired back to itself.

    :return: (acyclic edges, feedback edges)
    """
    graph = nx.DiGraph(edges)
    graph.add_edges_from(edges)
    looped = {node for component in nx.strongly_connected_components(graph)
              if len(component) > 1 for node in component}
    acyclic, feedback = [], []
    for source, target in sorted(set(edges)):
        if source == target or (source in looped and target in looped):
            feedback.append((source, target))
        else:
            acyclic.append((source, target))
    # Inside a component only the edges that close a cycle need to go; keep as
    # many as stay acyclic, so the path model still sees the loop body's order.
    kept = nx.DiGraph(acyclic)
    for source, target in list(feedback):
        kept.add_edge(source, target)
        if nx.is_directed_acyclic_graph(kept):
            acyclic.append((source, target))
            feedback.remove((source, target))
        else:
            kept.remove_edge(source, target)
    return sorted(acyclic), sorted(feedback)


def learn(monitoring_dirs, prefix: str = "monitor") -> Dict[str, Any]:
    """Build per-PE cost statistics from one or more monitoring directories.

    Runs are pooled per PE. Where the same PE appears in several runs its
    aggregate work is averaged, so a mapping that happened to process more items
    does not outweigh a smaller one.

    :return: {"baseline_rss_bytes", "reference_items", "pes": {pe_id: {...}}}
    """
    per_pe: Dict[str, List[dict]] = {}
    floors, runs = [], []
    # The dataflow graph, so a prediction made from bare instance counts still
    # knows which stages run beside each other and which queue behind.
    edges = set()
    seen = set()   # (run_id, instance_id): directories are often copies

    for monitoring_dir in (pathlib.Path(d) for d in monitoring_dirs):
        rows = _instance_rows(monitoring_dir, prefix)
        fresh = [r for r in rows
                 if (r.get("run_id"), r.get("instance_id")) not in seen]
        if not fresh:
            logger.info("%s holds run %s, already read; skipping",
                        monitoring_dir.name, rows[0].get("run_id"))
            continue
        if len(fresh) != len(rows):
            logger.warning("%s: %d of %d instances already seen in another "
                           "directory", monitoring_dir.name,
                           len(rows) - len(fresh), len(rows))
        seen.update((r.get("run_id"), r.get("instance_id")) for r in fresh)
        rows = fresh
        runs.append(monitoring_dir.name)
        edges.update(_abstract_edges(monitoring_dir, prefix))
        # One process per instance means per-instance RSS is that PE's own
        # process; a shared pid means the run cannot attribute memory at all.
        pids = {row.get("process_ids") for row in rows}
        shared = len(pids) == 1 and len(rows) > 1
        if shared:
            logger.info("%s runs every PE in one process; using it for the "
                        "memory baseline only", monitoring_dir.name)
        for row in rows:
            floors.append(_number(row, "rss_min_bytes"))
            canonical = canonical_pe(row["pe_id"])
            per_pe.setdefault(canonical, []).append(
                {**row, "_shared": shared, "_pe_id": row["pe_id"]})

    if not per_pe:
        raise ValueError("no instance rows found in any monitoring directory")

    baseline = min(f for f in floors if f > 0) if any(floors) else 0.0

    acyclic_edges, feedback_edges = split_feedback(edges)
    # A PE inside a loop processes each item several times, so its call count is
    # not the stream width. The width comes from the PEs that see every item
    # exactly once; only if every PE is looped is there nothing better to use.
    looped = {pe for edge in feedback_edges for pe in edge}
    if feedback_edges:
        logger.warning("feedback loop(s) in the dataflow: %s. The graph is kept "
                       "acyclic and the extra passes show up as an iteration "
                       "count per PE.",
                       ", ".join(f"{a}->{b}" for a, b in feedback_edges))

    # Reference stream width per run: how many items the run pushed through.
    reference = {}
    for pe_id, pe_rows in per_pe.items():
        if pe_id in looped and len(looped) < len(per_pe):
            continue
        for row in pe_rows:
            run = row.get("run_id", "")
            reference[run] = max(reference.get(run, 0.0),
                                 sum(_number(r, "total_count") for r in pe_rows
                                     if r.get("run_id") == run))

    pes: Dict[str, Any] = {}
    for pe_id, rows in sorted(per_pe.items()):
        by_run: Dict[str, List[dict]] = {}
        for row in rows:
            by_run.setdefault(row.get("run_id", ""), []).append(row)

        per_run, total_work, total_cpu, total_items = {}, 0.0, 0.0, 0.0
        for run, run_rows in by_run.items():
            work = sum(_number(r, "total_secs") for r in run_rows)
            cpu = sum(_number(r, "total_cpu_secs") for r in run_rows)
            items = sum(_number(r, "total_count") for r in run_rows)
            total_work += work
            total_cpu += cpu
            total_items += items
            per_run[run] = {
                "instances": len(run_rows),
                "items": items,
                "work_secs": work,
                "cpu_secs": cpu,
                # share of the run's stream this PE saw; data-dependent for a
                # filtering PE, so it is kept per run rather than averaged away
                "selectivity": items / reference[run] if reference.get(run) else 0.0,
                "secs_per_item": work / items if items else 0.0,
            }

        # An instance that processed no items has no measurement window, so its
        # cpu/memory columns are blank. It contributes nothing to the totals,
        # but averaging it in as a zero would understate every rate -- with five
        # items over four instances, one instance is always idle.
        busy = [r for r in rows if _number(r, "total_count") > 0]
        idle = len(rows) - len(busy)
        if idle:
            logger.info("%s: %d instance(s) processed nothing; excluded from "
                        "the cpu and memory averages", pe_id, idle)

        attributable = [
            max(0.0, _number(r, "rss_mean_endpoint_bytes") - baseline)
            for r in busy
            if not r["_shared"] and (r.get("rss_mean_endpoint_bytes") or "").strip()
        ]
        percents = [_number(r, "cpu_percent") for r in busy
                    if (r.get("cpu_percent") or "").strip()]
        # Prefer the widest run: selectivity measured on more items is the one
        # worth extrapolating from.
        widest = max(per_run, key=lambda run: reference.get(run, 0.0))

        # Repeat runs of one configuration are the only measure of how much a
        # cost varies. An LLM-calling PE can differ by more than 1.5x between
        # identical runs, which bounds how precise any prediction can be, so the
        # range is carried rather than collapsed into the mean.
        observed = [r["secs_per_item"] for r in per_run.values() if r["items"]]

        pes[pe_id] = {
            # per-item costs are the primitives: independent of how many items a
            # run happened to push and of how many instances ran
            "secs_per_item": total_work / total_items if total_items else 0.0,
            "secs_per_item_min": min(observed) if observed else 0.0,
            "secs_per_item_max": max(observed) if observed else 0.0,
            "runs_observed": len(observed),
            "cpu_secs_per_item": total_cpu / total_items if total_items else 0.0,
            "cpu_percent": statistics.mean(percents) if percents else 0.0,
            "selectivity": per_run[widest]["selectivity"],
            # >1 means the PE is revisited: a loop feeds each item back for
            # another pass. The cost model needs no special case, because every
            # pass was a real process() call and is counted as an item.
            "iterations": (per_run[widest]["selectivity"]
                           if pe_id in looped else 1.0),
            "in_feedback_loop": pe_id in looped,
            # memory a PE adds on top of the interpreter baseline. Process RSS is
            # never summed across PEs or instances -- the instrumentation records
            # max observed RSS, so a sum would count the baseline many times.
            "rss_attributable_bytes": (statistics.mean(attributable)
                                       if attributable else None),
            "rss_delta_mean_bytes": statistics.mean(
                [_number(r, "rss_delta_mean_bytes") for r in busy] or [0.0]),
            "idle_instances": idle,
            "aliases": sorted({r["_pe_id"] for r in rows}),
            "runs": per_run,
        }

    return {"baseline_rss_bytes": baseline,
            "reference_items": max(reference.values()) if reference else 0.0,
            "runs": runs,
            "edges": acyclic_edges,
            "feedback_edges": feedback_edges,
            "pes": pes}


def save(stats: Dict[str, Any], path: pathlib.Path) -> pathlib.Path:
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(stats, indent=2))
    return path


def load(path: pathlib.Path) -> Dict[str, Any]:
    return json.loads(pathlib.Path(path).read_text())


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("monitoring_dirs", nargs="+", type=pathlib.Path)
    parser.add_argument("-o", "--out", type=pathlib.Path, default=None,
                        help="write the stats as JSON")
    parser.add_argument("--prefix", default="monitor")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    stats = learn(args.monitoring_dirs, args.prefix)
    print(f"baseline RSS: {stats['baseline_rss_bytes']/1e6:.1f} MB   "
          f"widest run: {stats['reference_items']:.0f} items")
    print(f"{'PE':28} {'sel':>6} {'s/item':>10} {'cpu s/item':>11} "
          f"{'cpu%':>6} {'mem MB':>8}")
    for pe_id, pe in stats["pes"].items():
        mem = pe["rss_attributable_bytes"]
        print(f"{pe_id:28} {pe['selectivity']:>6.2f} {pe['secs_per_item']:>10.6f} "
              f"{pe['cpu_secs_per_item']:>11.6f} {pe['cpu_percent']:>6.1f} "
              f"{'n/a' if mem is None else f'{mem/1e6:>8.1f}'}")
    if args.out:
        print("wrote", save(stats, args.out))


if __name__ == "__main__":
    main()

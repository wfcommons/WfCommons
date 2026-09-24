#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 3: write synthetic instances from the installed recipe, and check them.

Each requested size is an upper bound: WfChef grows the base graph one
microstructure at a time, so 100 typically lands a little short. Anything that
is not a realisable dispel4py workflow is refused rather than left on disk.

    python -m wfcommons.wfstream.generate_workflows 100 50 75
"""

import collections
import json
import pathlib

import networkx as nx

from wfcommons import WorkflowGenerator

from .config import GENERATE_SIZES, RECIPE_NAME, SYNTHETIC_DIR, VERIFY_SIZES
from .resource_stats import pe_type
from .streaming_recipe import load_recipe, simple_base_graph, streaming_recipe


def _graph(spec: dict) -> nx.DiGraph:
    """The instance's task graph, keyed by task id."""
    graph = nx.DiGraph()
    for task in spec["tasks"]:
        graph.add_node(task["id"])
        for parent in task["parents"]:
            graph.add_edge(parent, task["id"])
    return graph


def _levels(spec: dict) -> list:
    """Widths of the real topological levels, source first."""
    return [len(gen) for gen in nx.topological_generations(_graph(spec))]


def source_instances(spec: dict) -> dict:
    """Instances per source PE: the PEs nothing feeds.

    Found from the graph rather than by name, so this holds for any workflow
    instead of only ones whose reader happens to be called ``read0``.
    """
    sources = collections.Counter(
        pe_type(task["id"]) for task in spec["tasks"] if not task["parents"])
    return dict(sources)


def check_instance(path: pathlib.Path) -> dict:
    """Validate a written instance; returns a report with an `errors` list."""
    spec = json.loads(path.read_text())["workflow"]["specification"]
    metrics, widths = spec["metrics"], _levels(spec)
    counts = collections.Counter(pe_type(t["id"]) for t in spec["tasks"])
    sources = source_instances(spec)
    components = nx.number_weakly_connected_components(_graph(spec))

    errors = []
    if not (metrics["numberOfLevels"] == len(widths)
            and metrics["minimumWidth"] == min(widths)
            and metrics["maximumWidth"] == max(widths)):
        errors.append(f"metrics {metrics} disagree with widths {widths}")
    # A dispel4py source PE always runs in exactly one process. More than one
    # instance of one means the requested size was too small for a base graph
    # with parallelism, so WfChef replicated the whole pipeline, reader and all.
    multiplied = {pe: n for pe, n in sources.items() if n != 1}
    if multiplied:
        errors.append(f"source PE(s) replicated: {multiplied}; a dispel4py "
                      f"source always runs in exactly one process")
    # A streaming workflow is one connected dataflow; several components mean
    # the graph was copied rather than widened.
    if components != 1:
        errors.append(f"{components} disconnected components; a dispel4py "
                      f"workflow is one dataflow")
    return {"metrics": metrics, "widths": widths, "counts": counts,
            "sources": sources, "components": components, "errors": errors}


def stamp_resources(workflow, stats) -> None:
    """Attach measured CPU and memory to a generated workflow's tasks.

    WfCommons' generator leaves `avg_cpu` and `memory` unset, so WfFormat's
    avgCPU and memoryInBytes come out null. Filling them from the traces is what
    lets a consumer of the instance see per-task resource use without also
    holding the statistics file.
    """
    from . import resource_stats as _rs
    if isinstance(stats, (str, pathlib.Path)):
        stats = _rs.load(stats)
    baseline = stats["baseline_rss_bytes"]
    missing = set()
    for task in workflow.tasks.values():
        kind = _rs.pe_type(task.task_id)
        pe = stats["pes"].get(kind)
        if pe is None:
            missing.add(kind)
            continue
        task.avg_cpu = pe["cpu_percent"]
        task.memory = int(baseline + (pe["rss_attributable_bytes"] or 0.0))
    if missing:
        print(f"      no statistics for {', '.join(sorted(missing))}; "
              f"left without cpu/memory")


def generate(sizes=None,
             output_dir: pathlib.Path = None,
             name: str = RECIPE_NAME,
             grow_from: str = None,
             stats=None) -> list:
    """Write one synthetic instance per requested task count.

    :param grow_from: the base graph to replicate PE instances into, named as
        ``<trace dir>-<task count>``. Defaults to the recipe's simple run, found
        by `simple_base_graph`: the one pipeline with a single instance per PE,
        so growing it adds instances rather than multiplying a shape that
        already has some. That name differs per workflow, hence the lookup.
    :param stats: `resource_stats` output (or a path to it). When given, each
        generated task carries the measured CPU and memory for its PE.
    """
    sizes = tuple(sizes) if sizes else GENERATE_SIZES
    output_dir = output_dir or SYNTHETIC_DIR
    output_dir.mkdir(parents=True, exist_ok=True)
    recipe_class = streaming_recipe(name, grow_from)

    written = []
    for num_tasks in sizes:
        workflow = WorkflowGenerator(recipe_class(num_tasks=num_tasks)).build_workflow()
        workflow.name = f"{name}-{num_tasks}"

        # WfChef names every task after its *type* while ids stay unique, but
        # wfcommons builds the graph behind the metrics block keyed by name and
        # draws its edges from parents, which are ids. When the two disagree the
        # edges land on auto-created nodes and the metrics describe a graph that
        # does not exist. Real traces have name == id; make these match too.
        for task in workflow.tasks.values():
            task.name = task.task_id

        if stats is not None:
            stamp_resources(workflow, stats)

        path = output_dir / f"{workflow.name}.json"
        workflow.write_json(path)

        report = check_instance(path)
        print(f"  {path.name}: {report['metrics']['numberOfTasks']} tasks, "
              f"{report['metrics']['numberOfFiles']} streams, "
              f"widths={report['widths']}")
        print("    " + ", ".join(f"{k}={v}" for k, v in sorted(report["counts"].items())))
        print(f"    sources: {report['sources']}, "
              f"components: {report['components']}")
        if report["errors"]:
            path.unlink()   # never leave a bad instance for something to train on
            raise SystemExit(f"{workflow.name}: " + "; ".join(report["errors"]) +
                             "\n  try a larger size, or add a trace whose "
                             "parallelism matches the scale you want "
                             "(the bad instance was not kept)")
        written.append(path)
    return written


def verify(sizes=None, name: str = RECIPE_NAME) -> None:
    """Build a few workflows from the *stock* recipe and report their shape.

    A sanity check on the installed recipe itself, before the streaming scaling
    rule is applied.
    """
    recipe_class = load_recipe(name)
    for num_tasks in (sizes or VERIFY_SIZES):
        workflow = WorkflowGenerator(recipe_class(num_tasks=num_tasks)).build_workflow()
        graph = nx.DiGraph((p, t) for t in workflow.tasks
                           for p in workflow.predecessors(t))
        graph.add_nodes_from(workflow.tasks)
        components = nx.number_weakly_connected_components(graph)
        sources = collections.Counter(
            pe_type(t) for t in workflow.tasks if graph.in_degree(t) == 0)
        files = {f.file_id for t in workflow.tasks.values()
                 for f in t.input_files + t.output_files}
        ok = all(n == 1 for n in sources.values()) and components == 1
        print(f"  num_tasks={num_tasks:4}: {len(workflow.tasks):4} tasks, "
              f"sources={dict(sources)}, components={components}, "
              f"streams={len(files)}  [{'ok' if ok else 'CHECK'}]")


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sizes", nargs="*", type=int,
                        help=f"task counts to generate (default: {GENERATE_SIZES})")
    parser.add_argument("-o", "--out", type=pathlib.Path, default=SYNTHETIC_DIR)
    parser.add_argument("-n", "--name", default=RECIPE_NAME)
    parser.add_argument("-g", "--grow-from", default=None,
                        help="base graph to grow (default: the recipe's simple run)")
    parser.add_argument("--verify", action="store_true",
                        help="also check the stock recipe at VERIFY_SIZES first")
    args = parser.parse_args()

    if args.verify:
        print("verifying the installed recipe")
        verify(name=args.name)
    print(f"generating -> {args.out}")
    generate(args.sizes, args.out, args.name, args.grow_from)


if __name__ == "__main__":
    main()

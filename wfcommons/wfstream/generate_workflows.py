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


def check_instance(path: pathlib.Path) -> dict:
    """Validate a written instance; returns a report with an `errors` list."""
    spec = json.loads(path.read_text())["workflow"]["specification"]
    metrics, widths = spec["metrics"], _levels(spec)
    counts = collections.Counter(t["id"].rsplit("_", 1)[0] for t in spec["tasks"])
    readers = counts.get("read0", 0)
    components = nx.number_weakly_connected_components(_graph(spec))

    errors = []
    if not (metrics["numberOfLevels"] == len(widths)
            and metrics["minimumWidth"] == min(widths)
            and metrics["maximumWidth"] == max(widths)):
        errors.append(f"metrics {metrics} disagree with widths {widths}")
    # A dispel4py source PE always runs in exactly one process, and a pipeline is
    # one connected graph. More than one of either means the requested size was
    # too small for a base graph with parallelism, so WfChef replicated the whole
    # pipeline, reader included.
    if readers != 1 or components != 1:
        errors.append(f"{readers} reader(s) in {components} component(s); a "
                      f"dispel4py workflow has exactly one of each")
    return {"metrics": metrics, "widths": widths, "counts": counts,
            "readers": readers, "components": components, "errors": errors}


def generate(sizes=None,
             output_dir: pathlib.Path = None,
             name: str = RECIPE_NAME,
             grow_from: str = None) -> list:
    """Write one synthetic instance per requested task count.

    :param grow_from: the base graph to replicate PE instances into, named as
        ``<trace dir>-<task count>``. Defaults to the recipe's simple run, found
        by `simple_base_graph`: the one pipeline with a single instance per PE,
        so growing it adds instances rather than multiplying a shape that
        already has some. That name differs per workflow, hence the lookup.
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

        path = output_dir / f"{workflow.name}.json"
        workflow.write_json(path)

        report = check_instance(path)
        print(f"  {path.name}: {report['metrics']['numberOfTasks']} tasks, "
              f"{report['metrics']['numberOfFiles']} streams, "
              f"widths={report['widths']}")
        print("    " + ", ".join(f"{k}={v}" for k, v in sorted(report["counts"].items())))
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
        counts = collections.Counter(t.rsplit("_", 1)[0] for t in workflow.tasks)
        graph = nx.DiGraph((p, t) for t in workflow.tasks
                           for p in workflow.predecessors(t))
        components = nx.number_weakly_connected_components(graph)
        files = {f.file_id for t in workflow.tasks.values()
                 for f in t.input_files + t.output_files}
        status = "ok" if counts.get("read0") == 1 and components == 1 else "CHECK"
        print(f"  num_tasks={num_tasks:4}: {len(workflow.tasks):4} tasks, "
              f"read0={counts.get('read0')}, components={components}, "
              f"streams={len(files)}  [{status}]")


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

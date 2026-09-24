#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""The two flows a registry drives, as callable entry points.

The registry (Laminar's) decides which case applies and runs the workflow
itself; wfstream is handed the monitoring directories those runs produced.

  new workflow    -> `on_new_workflow`: convert, cook, install, generate an
                     instance at the size the user asked for, then simulate.
  new size, known
  workflow        -> `on_new_size_run`: convert and store the new run, re-cook
                     the recipe, and stop. The next usage does the simulating.
"""

import datetime
import json
import logging
import pathlib

from . import (build_recipe, convert_traces, generate_workflows,
               resource_stats, simulate as simulate_module, update_traces)
from .config import BUILD_DIR, RECIPE_NAME, SYNTHETIC_DIR, WFFORMAT_DIR

logger = logging.getLogger(__name__)


def _task_counts(instances) -> list:
    """Task count of each converted instance."""
    return [len(json.loads(p.read_text())["workflow"]["specification"]["tasks"])
            for p in instances]


def check_against(stats: dict, monitoring_dir, name: str = None) -> list:
    """Predict each run in a monitoring directory and score it against itself.

    Called with statistics learned *before* the run existed, this measures how
    well the model extrapolates rather than how well it fits.
    """
    graph = None
    if stats.get("edges"):
        import networkx as nx
        graph = nx.DiGraph(tuple(e) for e in stats["edges"])

    records = []
    for run_id, observed in resource_stats.measured(monitoring_dir).items():
        try:
            predicted = simulate_module.simulate(
                observed["shape"], stats, items=observed["items"], graph=graph)
        except SystemExit as error:      # no overlap between run and statistics
            logger.warning("cannot score run %s: %s", run_id, error)
            continue
        record = simulate_module.compare(predicted, observed)
        record.update({"name": name, "run_id": run_id,
                       "source": str(monitoring_dir),
                       "checked_at": datetime.datetime.now().astimezone().isoformat()})
        records.append(record)
        ratio = record["bottleneck_ratio"]
        if ratio:
            logger.info("run %s: predicted %s at %.2fx the measured time",
                        run_id, record["bottleneck"], ratio)
    return records


def on_new_workflow(trace_dirs,
                    num_tasks: int,
                    items: float = None,
                    name: str = RECIPE_NAME,
                    wfformat_dir: pathlib.Path = None,
                    build_dir: pathlib.Path = None,
                    synthetic_dir: pathlib.Path = None,
                    input_files=None,
                    output_files=None,
                    grow_from: str = None,
                    stats_path: pathlib.Path = None,
                    register: bool = False,
                    simulate: bool = True) -> dict:
    """Bring a workflow the registry has not seen before up to a prediction.

    :param trace_dirs: monitoring directories from the registry's quick runs --
        the same workflow at a few scales (simple, and multi at 2x and 3x that
        size). One scale alone yields no microstructures, so the recipe cannot
        grow the graph and generation fails.
    :param num_tasks: the size the user asked about.
    :param items: how many items to push through the workflow when predicting.
        Defaults to the widest traced run.
    :param grow_from: the base graph generation replicates into. Defaults to
        the recipe's simple run, found by name, so a caller does not have to
        know it.
    :param stats_path: where to write the learned CPU/memory statistics;
        defaults to resource_stats.json beside the corpus.
    :param register: pip-install the cooked recipe so it becomes a WfChef entry
        point, visible to ``wfchef ls`` and to `get_recipe` in a fresh
        environment. Off by default because it shells out to pip; the faster
        path copies the data into the already-installed recipe instead.
    :param simulate: predict runtime, CPU and memory for the generated instance.
    :return: {"instances", "recipe", "synthetic", "stats", "simulation"}
    """
    wfformat_dir = wfformat_dir or WFFORMAT_DIR
    build_dir = build_dir or BUILD_DIR
    synthetic_dir = synthetic_dir or SYNTHETIC_DIR

    # The corpus is rebuilt from scratch: these traces define the new workflow.
    instances = convert_traces.convert_traces(
        trace_dirs, wfformat_dir, clean=True,
        input_files=input_files, output_files=output_files)
    counts = _task_counts(instances)
    if len(set(counts)) < 2:
        logger.warning(
            "all %d trace(s) are the same size %s; WfChef finds microstructures "
            "by comparing instances of different sizes, so the recipe will not "
            "be able to scale", len(instances), counts)

    cooked = build_recipe.cook(wfformat_dir, build_dir, name)
    if register:
        # pip-installs the cooked package, so the recipe becomes a WfChef entry
        # point and `wfchef ls` shows it like any built-in one
        build_recipe.register(build_dir)
    else:
        build_recipe.install(cooked, name=name)

    # Time, CPU and memory per PE, learned straight from the monitoring CSVs.
    # WfChef's own statistics cover runtime only, so this is kept alongside.
    stats = resource_stats.learn(trace_dirs)
    stats_path = stats_path or wfformat_dir.parent / "resource_stats.json"
    resource_stats.save(stats, stats_path)

    synthetic = generate_workflows.generate([num_tasks], synthetic_dir, name,
                                            grow_from, stats)[0]

    simulation, prediction_path, summary_path = None, None, None
    if simulate:
        simulation = simulate_module.simulate(synthetic, stats, items=items)
        # A prediction that is only returned is lost when the caller exits, so
        # there is nothing to hold the eventual real run against.
        prediction_path = simulate_module.save(
            simulation, synthetic.with_suffix(".prediction.json"), name=name,
            instance=synthetic, stats=stats_path, num_tasks=num_tasks)
        summary_path = synthetic.with_suffix(".summary.txt")
        summary_path.write_text(simulate_module.summary(simulation, name) + "\n")

    return {"instances": instances, "recipe": cooked, "synthetic": synthetic,
            "stats": stats_path, "simulation": simulation,
            "prediction": prediction_path, "summary": summary_path}


def on_new_size_run(trace_dirs,
                    name: str = RECIPE_NAME,
                    wfformat_dir: pathlib.Path = None,
                    build_dir: pathlib.Path = None,
                    input_files=None,
                    output_files=None,
                    stats_path: pathlib.Path = None,
                    log_path: pathlib.Path = None,
                    force: bool = False) -> dict:
    """Fold a real run at a new size into a known workflow's recipe.

    Nothing is generated or simulated here -- the point is that the recipe is
    better the next time someone asks. Re-cooking is skipped when the run was
    already in the corpus.

    :param log_path: where predicted-vs-measured records are appended;
        defaults to accuracy.jsonl beside the corpus.
    :return: {"added", "skipped", "recipe", "stats", "accuracy"}
    """
    wfformat_dir = wfformat_dir or WFFORMAT_DIR
    build_dir = build_dir or BUILD_DIR

    added, skipped = update_traces.merge_traces(
        trace_dirs, wfformat_dir, force=force,
        input_files=input_files, output_files=output_files)
    if not added:
        logger.info("no new runs to store (%s already in the corpus); "
                    "leaving the recipe alone", ", ".join(skipped))
        return {"added": [], "skipped": skipped, "recipe": None,
                "stats": None, "accuracy": []}

    # Before the new run is folded in, the model has never seen it: predicting
    # it now and checking the result is a genuine held-out test, and the only
    # accuracy record that accumulates by itself.
    stats_path = stats_path or wfformat_dir.parent / "resource_stats.json"
    checks = []
    if stats_path.exists():
        previous = resource_stats.load(stats_path)
        for monitoring_dir in trace_dirs:
            checks.extend(check_against(previous, monitoring_dir, name))
        if checks:
            log_path = log_path or wfformat_dir.parent / "accuracy.jsonl"
            with pathlib.Path(log_path).open("a") as handle:
                for check in checks:
                    handle.write(json.dumps(check, default=str) + "\n")

    cooked = build_recipe.cook(wfformat_dir, build_dir, name)
    build_recipe.install(cooked, name=name)

    # A real run at a new size sharpens the cost model as well as the recipe.
    stats = resource_stats.learn(trace_dirs)
    resource_stats.save(stats, stats_path)
    return {"added": added, "skipped": skipped, "recipe": cooked,
            "stats": stats_path, "accuracy": checks}

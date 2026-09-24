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

import json
import logging
import pathlib

from . import (build_recipe, convert_traces, generate_workflows,
               resource_stats, update_traces)
from .config import BUILD_DIR, RECIPE_NAME, SYNTHETIC_DIR, WFFORMAT_DIR

logger = logging.getLogger(__name__)


def _task_counts(instances) -> list:
    """Task count of each converted instance."""
    return [len(json.loads(p.read_text())["workflow"]["specification"]["tasks"])
            for p in instances]


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
    build_recipe.install(cooked)

    # Time, CPU and memory per PE, learned straight from the monitoring CSVs.
    # WfChef's own statistics cover runtime only, so this is kept alongside.
    stats = resource_stats.learn(trace_dirs)
    stats_path = stats_path or wfformat_dir.parent / "resource_stats.json"
    resource_stats.save(stats, stats_path)

    synthetic = generate_workflows.generate([num_tasks], synthetic_dir, name,
                                            grow_from, stats)[0]

    simulation = None
    if simulate:
        from .simulate import simulate as run_simulation
        simulation = run_simulation(synthetic, stats, items=items)

    return {"instances": instances, "recipe": cooked, "synthetic": synthetic,
            "stats": stats_path, "simulation": simulation}


def on_new_size_run(trace_dirs,
                    name: str = RECIPE_NAME,
                    wfformat_dir: pathlib.Path = None,
                    build_dir: pathlib.Path = None,
                    input_files=None,
                    output_files=None,
                    stats_path: pathlib.Path = None,
                    force: bool = False) -> dict:
    """Fold a real run at a new size into a known workflow's recipe.

    Nothing is generated or simulated here -- the point is that the recipe is
    better the next time someone asks. Re-cooking is skipped when the run was
    already in the corpus.

    :return: {"added", "skipped", "recipe", "stats"}
    """
    wfformat_dir = wfformat_dir or WFFORMAT_DIR
    build_dir = build_dir or BUILD_DIR

    added, skipped = update_traces.merge_traces(
        trace_dirs, wfformat_dir, force=force,
        input_files=input_files, output_files=output_files)
    if not added:
        logger.info("no new runs to store (%s already in the corpus); "
                    "leaving the recipe alone", ", ".join(skipped))
        return {"added": [], "skipped": skipped, "recipe": None, "stats": None}

    cooked = build_recipe.cook(wfformat_dir, build_dir, name)
    build_recipe.install(cooked)

    # A real run at a new size sharpens the cost model as well as the recipe.
    stats = resource_stats.learn(trace_dirs)
    stats_path = stats_path or wfformat_dir.parent / "resource_stats.json"
    resource_stats.save(stats, stats_path)
    return {"added": added, "skipped": skipped, "recipe": cooked,
            "stats": stats_path}

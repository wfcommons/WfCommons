#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Paths and per-trace settings shared by the wfstream scripts.

Everything installation-specific lives here so the step scripts stay generic.
"""

import pathlib

import wfcommons

# --- locations ---------------------------------------------------------------
TRACES = pathlib.Path("/home/taina/dispel4py_agentic_ai_traces")
WFFORMAT_DIR = TRACES / "wfformat"      # converted instances
SYNTHETIC_DIR = TRACES / "synthetic"    # generated instances
BUILD_DIR = pathlib.Path(__file__).resolve().parent / "build" / "climate"

RECIPE_NAME = "climate"
RECIPES_ROOT = pathlib.Path(wfcommons.__file__).parent / "wfchef" / "recipes"


def recipe_dir(name: str = None) -> pathlib.Path:
    """Where a cooked recipe is installed inside the wfcommons package."""
    return RECIPES_ROOT / f"wfchef_recipe_{name or RECIPE_NAME}"


RECIPE_DIR = recipe_dir(RECIPE_NAME)

# --- traces ------------------------------------------------------------------
# dispel4py receives the input path as a root input and hardcodes the output
# path, so neither reaches the trace; name them per monitoring directory here.
TRACE_DIRS = ["monitoring_simple", "monitoring_multi_16", "monitoring_multi_32"]
INPUT_FILES = {
    "monitoring_simple": ["sensor_data_agentic.json"],
    "monitoring_multi_16": ["sensor_data_parallel_100.json"],
    "monitoring_multi_32": ["sensor_data_parallel_100.json"],
}
OUTPUT_FILES = {
    "monitoring_simple": ["agentic_sensor_results.jsonl"],
    "monitoring_multi_16": ["agentic_parallel_results.jsonl"],
    "monitoring_multi_32": ["agentic_parallel_results.jsonl"],
}

# --- generation --------------------------------------------------------------
# The base graph to grow is not configured: it is the simple run, found by
# `streaming_recipe.simple_base_graph`, so this works for any workflow.
GENERATE_SIZES = (100, 50, 75)      # synthetic instances to write, by task count
VERIFY_SIZES = (50, 100, 200)       # sizes checked after installing the recipe


def trace_dirs(names=None):
    """Resolve monitoring directory names to paths (defaults to TRACE_DIRS)."""
    return [TRACES / name for name in (names or TRACE_DIRS)]

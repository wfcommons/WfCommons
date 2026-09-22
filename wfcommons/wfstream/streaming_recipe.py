#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""How a streaming workflow scales: more PE instances, not more pipelines.

WfChef's own duplication can only use microstructures discovered *within* the
chosen base graph. A one-instance-per-PE trace has none, so it falls back to
copying the whole graph, readers included -- wrong for dispel4py, where scale
means more instances of a PE.
"""

import importlib
import itertools
import logging
import pickle

import networkx as nx

from wfcommons.wfchef.duplicate import duplicate_nodes

from .config import RECIPE_NAME, recipe_dir

logger = logging.getLogger(__name__)

FICTITIOUS = ("SRC", "DST")   # WfChef bookends; not tasks


def load_recipe(name: str = RECIPE_NAME):
    """Import the installed WfChef recipe class for ``name``."""
    module = importlib.import_module(
        f"wfcommons.wfchef.recipes.wfchef_recipe_{name}.recipe")
    return getattr(module, f"{name.capitalize()}Recipe")


def _tasks(graph) -> list:
    return [n for n in graph if n not in FICTITIOUS]


def base_graphs(name: str = RECIPE_NAME) -> dict:
    """Every base graph in a recipe, by name."""
    root = recipe_dir(name) / "microstructures"
    if not root.is_dir():
        raise SystemExit(f"no installed recipe for {name!r} at {recipe_dir(name)}")
    graphs = {}
    for directory in sorted(root.iterdir()):
        pickled = directory / "base_graph.pickle"
        if directory.is_dir() and directory.name != "metric" and pickled.exists():
            graphs[directory.name] = pickle.loads(pickled.read_bytes())
    if not graphs:
        raise SystemExit(f"recipe {name!r} has no base graphs under {root}")
    return graphs


def simple_base_graph(name: str = RECIPE_NAME) -> str:
    """Name the recipe's simple run -- the base graph generation grows from.

    The simple run is the one with a single instance per PE: one plain pipeline,
    the only shape where replicating a PE adds parallelism rather than
    multiplying parallelism that is already there. It is found rather than
    configured, so this works for any registered workflow.

    Falls back to the smallest graph, with a warning, when no run is
    unparallelised -- growing that one will not produce clean instance counts.
    """
    graphs = base_graphs(name)
    by_size = sorted(graphs.items(), key=lambda kv: len(_tasks(kv[1])))
    for graph_name, graph in by_size:
        nodes = _tasks(graph)
        if len({graph.nodes[n]["type"] for n in nodes}) == len(nodes):
            return graph_name
    smallest = by_size[0][0]
    logger.warning(
        "no run of %r has one instance per PE, so there is no simple graph to "
        "grow; using the smallest (%s). Register a single-process run of the "
        "workflow to get clean scaling.", name, smallest)
    return smallest


def grow_pipeline(num_tasks: int,
                  grow_from: str = None,
                  name: str = RECIPE_NAME) -> nx.DiGraph:
    """Grow the base pipeline to ``num_tasks`` by replicating its PE instances.

    Every non-source PE is replicable, which is dispel4py's rule: any PE can take
    ``numprocesses > 1`` except a source, which always runs in one process. Types
    grow round-robin so instance counts stay balanced, and a replica inherits its
    original's wiring -- leaving consecutive stages fully connected, as the
    default shuffle routing makes them. Deterministic for a given ``num_tasks``.

    :param grow_from: base graph to grow, as ``<trace dir>-<task count>``.
        Defaults to the recipe's simple run, found by `simple_base_graph`.
    """
    graphs = base_graphs(name)
    if grow_from is None:
        grow_from = simple_base_graph(name)
    if grow_from not in graphs:
        raise SystemExit(f"no base graph {grow_from!r} in the {name} recipe; "
                         f"available: {', '.join(sorted(graphs))}")
    graph = graphs[grow_from]

    if num_tasks < len(_tasks(graph)):
        raise SystemExit(f"cannot build {num_tasks} tasks from {grow_from}, which "
                         f"already has {len(_tasks(graph))}")

    source_types = {graph.nodes[n]["type"] for n in graph.successors("SRC")}
    growable = sorted({graph.nodes[n]["type"] for n in _tasks(graph)} - source_types)
    for pe_type in itertools.cycle(growable):
        if len(_tasks(graph)) >= num_tasks:
            break
        # replicate an original rather than a copy, so each new instance inherits
        # the full set of neighbours and the stages stay all-to-all
        original = next(n for n in _tasks(graph)
                        if graph.nodes[n]["type"] == pe_type
                        and "duplicate_of" not in graph.nodes[n])
        duplicate_nodes(graph, {original})
    return graph


def streaming_recipe(name: str = RECIPE_NAME, grow_from: str = None):
    """The installed recipe, with `grow_pipeline` in place of its own scaling."""
    base_class = load_recipe(name)

    class StreamingRecipe(base_class):
        def generate_nx_graph(self):
            return grow_pipeline(self.num_tasks, grow_from, name)

    StreamingRecipe.__name__ = f"Streaming{base_class.__name__}"
    return StreamingRecipe

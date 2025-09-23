#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021-2024 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import json
import math
import networkx as nx
import numpy as np
import pathlib
import pickle

from hashlib import sha256
from itertools import chain, combinations
from typing import Union, Set, Optional, Dict, List
from uuid import uuid4

from .utils import create_graph, combine_hashes, annotate, draw

this_dir = pathlib.Path(__file__).resolve().parent


def comb(n: int, k: int) -> float:
    """
    Calculates the combination of two integers.

    :param n: number.
    :type n: int.
    :param k: number.
    :type k: int.
    
    :return: combination of two integers.
    :rtype: float.
    """
    return math.factorial(n) / (math.factorial(k) * math.factorial(n - k))


class ImbalancedMicrostructureError(Exception):
    pass


def get_children(graph: nx.DiGraph, node: str) -> List[str]:
    """
    Gets the children of a node.

    :param graph: graph that contains the node.
    :type graph: networkX DiGraph.
    :param node: a node.
    :type node: str.
    
    :return: list of the node's children.
    :rtype: List[str].
    """
    return [child for _, child in graph.out_edges(node)]


def get_parents(graph: nx.DiGraph, node: str) -> List[str]:
    """
    Gets the parents of a node.

    :param graph: graph that contains the node.
    :type graph: networkX DiGraph.
    :param node: a node.
    :type node: str.
    
    :return: list of the node's parents.
    :rtype: List[str].
    """
    return [parent for parent, _ in graph.in_edges(node)]


def get_relatives(graph: nx.DiGraph, node: str) -> Set[str]:
    """
    Gets all node's relatives (children and parents).

    :param graph: graph that contains the node.
    :type graph: networkX DiGraph.
    :param node: a node.
    :type node: str.
    
    :return: set of node's relative.
    :rtype: Set[str].
    """
    return set(chain(get_children(graph, node), get_parents(graph, node)))


def find_microstructure(graph: nx.DiGraph, n1: str, n2: str):
    """
    Detects a pattern (microstructure).

    :param graph: graph.
    :type graph: networkX DiGraph.
    :param n1:  a node in graph.
    :type n1: str.
    :param n2:  a different node in graph.
    :type n2: str.
    
    :return: sets of n1 related nodes, n2 related nodes, the nodes in common between n1 and n2 and 
            all the nodes involved in the process.
    :rtype: Set[str], Set[str], Set[str], Set[str].

    """
    n1_new_friends = get_relatives(graph, n1)
    n2_new_friends = get_relatives(graph, n2)

    n1_friends = {n1}
    n2_friends = {n2}
    common_friends = set()
    all_friends = {n1, n2}
    n1_new_friends = {n1}
    n2_new_friends = {n2}

    while n1_new_friends or n2_new_friends:
        if not n1_new_friends or not n2_new_friends:
            raise ImbalancedMicrostructureError()

        n1_new_friends = set.union(set(), *[get_relatives(graph, friend) for friend in n1_new_friends]) - all_friends
        n2_new_friends = set.union(set(), *[get_relatives(graph, friend) for friend in n2_new_friends]) - all_friends
        common_friends.update(n1_new_friends.intersection(n2_new_friends))
        all_friends.update(n1_new_friends.union(n2_new_friends))

        n1_new_friends -= common_friends
        n2_new_friends -= common_friends

        n1_friends.update(n1_new_friends)
        n2_friends.update(n2_new_friends)

    return n1_friends, n2_friends, common_friends, all_friends


def find_microstructures(graph: nx.DiGraph, verbose: bool = False):
    """
    Detects the patterns (microstructures) that are used for replication and graph expansion.

    :param graph: graph.
    :type graph: networkX DiGraph.
    :param verbose: if set, prints status messages.
    :type verbose: networkX DiGraph.

    :return: patterns (microstructures)
    :rtype: Set[str].
    """

    if verbose:
        print("Sorting nodes by type hash and parent")
    nodes_by_type_hash: Dict[str, Set[str]] = {}
    for node in graph.nodes:
        for child in get_children(graph, node):
            th = graph.nodes[child]["type_hash"]
            nodes_by_type_hash.setdefault((node, th), set())
            nodes_by_type_hash[(node, th)].add(child)

    if verbose:
        print("Finding microstructures for typehashes")
    microstructures = {}
    for (th, _), nodes in nodes_by_type_hash.items():
        if len(nodes) < 2:
            continue
        # if verbose:
        #     print(f"{th:5}: {len(nodes)} siblings ({comb(len(nodes), 2)} combos)")
        for n1, n2 in combinations(nodes, r=2):
            try:
                ms1, ms2, _, _ = find_microstructure(graph, n1, n2)
                if len(ms1) > len(ms2):
                    ms1, ms2 = ms2, ms1
                    n1, n2 = n2, n1
                key = combine_hashes(*{graph.nodes[node]["type_hash"] for node in ms1})
                microstructures.setdefault(key, set())
                microstructures[key].update({frozenset(ms1), frozenset(ms2)})
            except ImbalancedMicrostructureError:
                continue

    return microstructures


def sort_graphs(workflow_path: Union[pathlib.Path],
                verbose: bool = False) -> List[nx.DiGraph]:
    """
    Sort graphs in crescent order of number of tasks.

    :param workflow_path: path to the JSON instances.
    :type workflow_path: pathlib.Path.
    :param verbose: if set, prints status messages.
    :type verbose: networkX DiGraph.

    :return: sorted graphs
    :rtype: List[networkX.DiGraph].
    """
    if verbose:
        print(f"Working on {workflow_path}")
    graphs = []

    for path in workflow_path.glob("*.json"):
        graph = create_graph(path)
        annotate(graph)
        graph.graph["name"] = path.stem
        graphs.append(graph)

    if not graphs:
        raise ValueError(f"No graphs found in {workflow_path}")

    if verbose:
        print("Constructed graphs")

    sorted_graphs = sorted(graphs, key=lambda graph: len(graph.nodes))
    return sorted_graphs


def save_microstructures(workflow_path: Union[pathlib.Path],
                         savedir: pathlib.Path,
                         verbose: bool = False,
                         img_type: Optional[str] = 'png',
                         cutoff: int = 4000,
                         highlight_all_instances: bool = False) -> List[nx.DiGraph]:
    summary = {
        "frequencies": {},
        "base_graphs": {}
    }

    for graph in sort_graphs(workflow_path, verbose):
        if graph.order() > cutoff:
            print(f'This and the next workflows have more than {cutoff} tasks')
            break

        if verbose:
            print(f"Running for {graph.name}")
        g_savedir = savedir.joinpath(graph.name)
        g_savedir.mkdir(exist_ok=True, parents=True)

        base_graph_path = g_savedir.joinpath("base_graph.pickle")
        with open(base_graph_path, 'wb') as pf:
            pickle.dump(graph, pf, pickle.HIGHEST_PROTOCOL)
        summary["base_graphs"][graph.name] = {
            "size": graph.size(),
            "order": graph.order()
        }

        if img_type:
            base_graph_image_path = g_savedir.joinpath(f"base_graph")
            if verbose:
                print(f"Drawing base graph to {base_graph_image_path}")
            draw(graph, close=True, legend=False, extension=img_type, save=str(base_graph_image_path))

        if verbose:
            print("Finding microstructures")

        microstructures = find_microstructures(graph, verbose=verbose)
        mdatas = {}
        for _, (ms_hash, instances) in enumerate(microstructures.items()):
            ms_name = f"microstructure_{ms_hash}"

            summary["frequencies"].setdefault(ms_name, [])
            summary["frequencies"][ms_name].append((graph.order(), len(instances)))
            mdatas[ms_name] = {
                "name": ms_name,
                "nodes": list(map(list, instances)),
                "frequency": len(instances),
            }
            if img_type:
                print(f"Drawing {ms_name}")
                draw(
                    graph,
                    subgraph=list(instances)[0] if not highlight_all_instances else set.union(*instances),
                    with_labels=False,
                    extension=img_type,
                    save=str(g_savedir.joinpath(ms_name)),
                    close=True
                )

        if verbose:
            print()

        g_savedir.joinpath("microstructures").with_suffix(".json").write_text(json.dumps(mdatas, indent=2))

    savedir.joinpath("summary").with_suffix(".json").write_text(json.dumps(summary, indent=2))

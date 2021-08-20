#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import pathlib
import json
import pickle
import networkx as nx
from typing import Set, List, Union, Dict
from uuid import uuid4
import numpy as np
import random

this_dir = pathlib.Path(__file__).resolve().parent


class NoMicrostructuresError(Exception):
    pass


def duplicate_nodes(graph: nx.DiGraph, nodes: Set[str]) -> Dict:
    """
    Replicates nodes of a graph.

    :param graph: graph used to replicate and attach new nodes.
    :type graph: networkX DiGraph
    :param nodes: nodes to be replicated. 
    :type nodes: Set[str]. 
    
    :return: the new nodes replicated.
    :rtype: Dict[str].
    """
    new_nodes = {}
    for node in nodes:
        new_node = f"{node}_{uuid4()}"
        graph.add_node(new_node, **graph.nodes[node])
        nx.set_node_attributes(graph, {new_node: node}, "duplicate_of")
        new_nodes[node] = new_node

    for node, new_node in new_nodes.items():
        for parent, _ in graph.in_edges(node):
            if parent in new_nodes:
                graph.add_edge(new_nodes[parent], new_node)
            else:
                graph.add_edge(parent, new_node)

        for _, child in graph.out_edges(node):
            if child in new_nodes:
                graph.add_edge(new_node, new_nodes[child])
            else:
                graph.add_edge(new_node, child)

    return new_nodes


def duplicate(path: pathlib.Path,
              base: Union[str, pathlib.Path],
              num_nodes: int) -> nx.DiGraph:
    """
    Attaches replicated nodes to base graph.

    :param path: path to the summary JSON file.
    :type path: pathlib.Path.
    :param base: name (for samples available in WfCommons) or path to the specific 
                 graph to be used as base (if not set WfChef chooses the best fitting one). 
    :type base: str or pathlib.Path.
    :param num_nodes: total amount of nodes desired in the synthetic instance.
    :type num_nodes: int.

    :return: graph with the desired number of tasks.
    :rtype: networkX DiGraph.
    """
    summary = json.loads(path.joinpath("summary.json").read_text())
   
    if base:
        base_path = pathlib.Path(base)
        if not base_path.is_absolute():
            base_path = path.joinpath(base_path)
    else:
        base_path = path.joinpath(min(summary["base_graphs"].keys(), key=lambda k: summary["base_graphs"][k]["order"]))
    
        
    graph = pickle.loads(base_path.joinpath("base_graph.pickle").read_bytes())
    if num_nodes < graph.order():
        raise ValueError(
            f"Cannot create synthentic graph with {num_nodes} nodes from base graph with {graph.order()} nodes")

    all_microstructures = json.loads(base_path.joinpath("microstructures.json").read_text())
    microstructures, freqs = map(list, zip(*[(ms, ms["frequency"]) for ms_hash, ms in all_microstructures.items()]))

    p: List[float] = (np.array(freqs) / np.sum(freqs)).tolist()
    while graph.order() < num_nodes and microstructures:
        i = random.choice(range(len(microstructures)))
        ms = microstructures[i]
        while ms["nodes"]:
            j = random.choice(range(len(ms["nodes"])))
            structure = ms["nodes"][j]
            if graph.order() + len(structure) > num_nodes:
                del ms["nodes"][j]
            else:
                break

        if not ms["nodes"]:  # delete microstructure
            del microstructures[i]
            del p[i]
            continue

        duplicate_nodes(graph, structure)

    return graph

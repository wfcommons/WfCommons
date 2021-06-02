#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import networkx as nx
import pathlib
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches
from typing import Iterable, Union, Set, Optional, Tuple, Hashable
import json
from hashlib import sha256

this_dir = pathlib.Path(__file__).resolve().parent


def string_hash(obj: Hashable) -> str:
    return sha256(str(obj).encode("utf-8")).hexdigest()


def type_hash(_type: str, parent_types: Iterable[str]) -> str:
    return string_hash((_type, sorted(set(parent_types))))


def combine_hashes(*hashes: str) -> str:
    return string_hash(sorted(hashes))


def create_graph(path: Union[str, pathlib.Path]) -> nx.DiGraph:
    """
    Creates a networkX DiGraph from a JSON file in the WfFormat.

    :param path: name (for samples available in WfCommons) or the path to graphs JSON.
    :type path: str or pathlib.Path.
    
    :return: graph.
    :rtype: networkX DiGraph.
    """
    path = pathlib.Path(path)
    with path.open() as fp:
        content = json.load(fp)

        graph = nx.DiGraph()

        # Add src/dst nodes
        graph.add_node("SRC", label="SRC", type="SRC", id="SRC")
        graph.add_node("DST", label="DST", type="DST", id="DST")

        id_count = 0

        for job in content['workflow']['jobs']:

            # specific for epigenomics -- have to think about how to do it in general
            if "genome-dax" in content['name']:
                _type, *_ = job['name'].split('_')
                graph.add_node(job['name'], label=_type, type=_type, id=str(id_count))
                id_count += 1
            else:
                try:
                    _type, _id = job['name'].split('_ID')
                except ValueError:
                    _type, _id = job['name'].split('_0')
                graph.add_node(job['name'], label=_type, type=_type, id=_id)

            # for job in content['workflow']['jobs']:
            for parent in job['parents']:
                graph.add_edge(parent, job['name'])

        for node in graph.nodes:

            if node in ["SRC", "DST"]:
                continue
            if graph.in_degree(node) <= 0:
                graph.add_edge("SRC", node)
            if graph.out_degree(node) <= 0:
                graph.add_edge(node, "DST")

        return graph


def annotate(g: nx.DiGraph) -> None:
    """
    Annotates a networkX DiGraph with metadata such as the tasks top-down type hash, 
    bottom-up type hash, and type-hash.

    :param path: name (for samples available in WfCommons) or the path to graphs JSON.
    :type path: str or pathlib.Path.
    
    :return: annotated graph.
    :rtype: networkX DiGraph.
    """

    visited = set()
    queue = [(node, 1) for node in g.nodes if g.in_degree(node) <= 0]
    while queue:
        cur, level = queue.pop(0)
        g.nodes[cur]["level"] = level
        g.nodes[cur]["label"] = g.nodes[cur]["id"]
        parent_ths = [
            g.nodes[p]["top_down_type_hash"]
            for p, _ in g.in_edges(cur)
        ]
        g.nodes[cur]["top_down_type_hash"] = type_hash(g.nodes[cur]["type"], parent_ths)

        visited.add(cur)
        queue.extend([
            (child, level + 1) for _, child in g.out_edges(cur)
            if child not in visited and
               {sib for sib, _ in g.in_edges(child)}.issubset(visited)
        ])

    # REVERSE 
    visited = set()
    queue = [node for node in g.nodes if g.out_degree(node) <= 0]
    while queue:
        cur = queue.pop(0)
        parent_ths = [
            g.nodes[p]["bottom_up_type_hash"]
            for _, p in g.out_edges(cur)
        ]
        g.nodes[cur]["bottom_up_type_hash"] = type_hash(g.nodes[cur]["type"], parent_ths)
        g.nodes[cur]["type_hash"] = combine_hashes(g.nodes[cur]["top_down_type_hash"],
                                                   g.nodes[cur]["bottom_up_type_hash"])

        visited.add(cur)
        queue.extend([
            child for child, _ in g.in_edges(cur)
            if child not in visited and
               {sib for _, sib in g.out_edges(child)}.issubset(visited)
        ])

    return g


def draw(g: nx.DiGraph,
         extension: Optional[str] = 'png',
         with_labels: bool = False,
         ax: Optional[plt.Axes] = None,
         show: bool = False,
         save: Optional[Union[pathlib.Path, str]] = None,
         close: bool = False,
         legend: bool = False,
         node_size: int = 1000,
         linewidths: int = 5,
         subgraph: Set[str] = set()) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plots a netwrokX DiGraph.

    :param g: graph to be plotted.
    :type g: networkX DiGraph.
    :type extension: extension of the output file.
    :param extension: str.
    :param with_labels: if set, it prints the task types over their nodes.
    :type with_labels: bool.
    :param ax: plot axes.
    :type ax: plt.Axes.
    :param show: if set, displays the plot on screen.
    :type show: bool.
    :param save: path to directory to save the plot.
    :type save: pathlib.Path.
    :param close: if set, automatically closes window that displays plot.
    :type close: bool.
    :param legend: if set, displays legend of the plot.
    :type legend: bool.
    :param node_size: size of the nodes (circles) in the plot.
    :type node_size: int.
    :param linewidths: thickness of the edges in the plot.
    :type linewidths: int.
    :param subgraph: nodes that were added by replication and will be colored green.
    :type subgraph: Set[str].

    

    :return: the figure and the axis used.
    :rtype:  Tuple[plt.Figure, plt.Axes].
    """
    fig: plt.Figure
    ax: plt.Axes
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 10))
    else:
        fig = ax.get_figure()

    node_border_colors = {}
    if isinstance(subgraph, dict):
        for color, nodes in subgraph.items():
            for node in nodes:
                node_border_colors[node] = color
    else:
        for node in subgraph:
            node_border_colors[node] = "green"

    pos = nx.nx_agraph.pygraphviz_layout(g, prog='dot')
    type_set = sorted({g.nodes[node]["type"] for node in g.nodes})  # not type-hash
    types = {
        t: i for i, t in enumerate(type_set)
    }
    node_color = [types[g.nodes[node]["type"]] for node in g.nodes]  # not type-hash
    for node in g.nodes:
        if node in subgraph:
            g.nodes[node]["node_shape"] = "s"
        else:
            g.nodes[node]["node_shape"] = "c"
    edgecolors = [node_border_colors.get(node, "white") for node in g.nodes]
    edge_color = [
        node_border_colors.get(src) if node_border_colors.get(src, -1) == node_border_colors.get(dst, 1) else "black"
        for src, dst in g.edges
    ]
    cmap = cm.get_cmap('rainbow', len(type_set))
    nx.draw(g, pos, node_size=node_size, node_color=node_color, edgecolors=edgecolors, edge_color=edge_color,
            linewidths=linewidths, cmap=cmap, ax=ax, with_labels=with_labels)
    color_lines = [mpatches.Patch(color=cmap(types[t]), label=t) for t in type_set]

    if legend:
        legend = ax.legend(handles=color_lines, loc='lower right')

    if show:
        plt.show()

    if save is not None:
        if extension == 'dot':
            nx.drawing.nx_agraph.write_dot(g, save)
        else:
            fig.savefig(f'{save}.{extension}')

    if close:
        plt.close(fig)

    return fig, ax

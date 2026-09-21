#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Render a WfFormat instance as a runnable dispel4py workflow script.

The inverse of ``dispel_to_wfformat``, and deliberately not an LLM translation:
the mapping is mechanical once you know where each piece comes from.

  * **shape of the graph** -- from a dispel4py monitoring trace, not from the
    instance. WfFormat records one task per PE instance and names streams after
    the *producing* port only, so the consumer-side port name (dispel4py's
    ``to_connection``) is not recoverable from an instance alone; WfChef's
    synthetic instances also wire branch ports to whichever consumer the sampling
    picked. The trace's abstract shape has both ends of every connection and does
    not change with scale, so it is the authority on wiring.
  * **scale** -- from the instance: one PE class per task category, with
    ``numprocesses`` set to how many tasks of that category the instance holds.
  * **cost** -- from the instance: each stage burns its recorded runtime, divided
    across the items that flow through it.

The generated PEs are timing stand-ins. They carry the workflow's structure and
cost so it can be run under dispel4py's mappings, not its science: no sensor
readings are parsed and no model is called.

Usage:

    python -m wfcommons.wfchef.wfformat_to_dispel4py climate-100.json \\
        --shape monitoring_simple -o run_climate_100.py
"""

from __future__ import annotations

import argparse
import collections
import json
import logging
import pathlib
import re
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

_TRAILING_INDEX = re.compile(r"(?P<name>.*?)(?P<index>\d+)$")

# WfChef bookends every graph with fictitious SRC and DST nodes. They are not
# tasks and have no dispel4py counterpart -- a workflow's real source is its
# ProducerPE -- so anything that reaches here through a graph rather than through
# an instance's task list has to be stripped, or they become phantom PEs.
FICTITIOUS = frozenset({"SRC", "DST"})


def _pe_name_and_index(pe_id: str) -> Tuple[str, int]:
    """
    Split a dispel4py pe_id into the PE's name and its position in the graph.

    ``WorkflowNode`` builds ids as ``pe.name + str(node_counter)``, so the
    category ``NormalizeDataPE1`` is a PE named ``NormalizeDataPE`` that was the
    second node added. Reproducing both makes the regenerated workflow trace back
    to the same ids.
    """
    match = _TRAILING_INDEX.match(pe_id)
    if not match or not match.group("name"):
        return pe_id, 0
    return match.group("name"), int(match.group("index"))


def read_shape(monitoring_dir: pathlib.Path, prefix: str = "monitor") -> Dict[str, Any]:
    """Read the abstract (PE-level) shape of a monitoring run."""
    shapes = sorted(monitoring_dir.glob(f"{prefix}_shape_run*.json"))
    if not shapes:
        raise FileNotFoundError(f"no {prefix}_shape_run*.json in {monitoring_dir}")
    if len(shapes) > 1:
        logger.warning("%d shapes in %s, using %s", len(shapes), monitoring_dir, shapes[-1].name)
    shape = json.loads(shapes[-1].read_text())

    # Drop SRC/DST before anything reads the ports: an edge SRC -> read0 would
    # otherwise give the source PE an input connection and disguise it as a
    # middle stage (ProducerPE becomes IterativePE, and the workflow no longer
    # has a source at all).
    dropped = [node for node in shape["nodes"] if node in FICTITIOUS]
    if dropped:
        logger.info("dropping fictitious node(s) %s from the shape", ", ".join(dropped))
    shape["nodes"] = [node for node in shape["nodes"] if node not in FICTITIOUS]
    shape["edges"] = [edge for edge in shape["edges"]
                      if edge["from"] not in FICTITIOUS and edge["to"] not in FICTITIOUS]
    return shape


def summarize_instance(instance_path: pathlib.Path) -> Dict[str, Dict[str, Any]]:
    """Count tasks per category and average their runtimes."""
    doc = json.loads(instance_path.read_text())
    spec = doc["workflow"]["specification"]
    runtimes = {t["id"]: t.get("runtimeInSeconds", 0.0)
                for t in doc["workflow"].get("execution", {}).get("tasks", [])}

    per_type: Dict[str, List[float]] = collections.defaultdict(list)
    for task in spec["tasks"]:
        category = task.get("category") or task["id"].rsplit("_", 1)[0]
        if category in FICTITIOUS or task["id"] in FICTITIOUS:
            logger.info("skipping fictitious %s node", task["id"])
            continue
        per_type[category].append(runtimes.get(task["id"], 0.0))

    return {category: {"count": len(values),
                       "runtime": sum(values) / len(values) if values else 0.0}
            for category, values in per_type.items()}


def _ports(shape: Dict[str, Any]) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """Input and output connection names of every PE, from the shape's edges."""
    inputs: Dict[str, List[str]] = collections.defaultdict(list)
    outputs: Dict[str, List[str]] = collections.defaultdict(list)
    for edge in shape["edges"]:
        if edge["to_connection"] not in inputs[edge["to"]]:
            inputs[edge["to"]].append(edge["to_connection"])
        if edge["from_connection"] not in outputs[edge["from"]]:
            outputs[edge["from"]].append(edge["from_connection"])
    return inputs, outputs


def _identifier(name: str) -> str:
    """A Python-safe variable name for a PE."""
    cleaned = re.sub(r"\W", "_", name)
    return f"pe_{cleaned}" if cleaned[0].isdigit() else cleaned


def render(instance_path: pathlib.Path,
           monitoring_dir: pathlib.Path,
           items: int = 100,
           prefix: str = "monitor") -> str:
    """Render the dispel4py script for one WfFormat instance."""
    shape = read_shape(monitoring_dir, prefix)
    summary = summarize_instance(instance_path)
    inputs, outputs = _ports(shape)

    unknown = sorted(set(summary) - set(shape["nodes"]))
    if unknown:
        raise ValueError(
            f"{instance_path.name} has categories the shape does not describe: "
            f"{unknown}. Point --shape at a trace of the same workflow."
        )

    # Order the PEs the way dispel4py numbered them, so ids line up on a re-trace.
    pes = sorted(shape["nodes"], key=lambda pe_id: _pe_name_and_index(pe_id)[1])
    sources = [pe for pe in pes if not inputs.get(pe)]

    total_processes = sum(1 if pe in sources else summary.get(pe, {}).get("count", 1)
                          for pe in pes)

    lines: List[str] = []
    add = lines.append

    add('#!/usr/bin/env python')
    add('"""')
    add(f'dispel4py workflow generated from {instance_path.name}.')
    add('')
    add('Timing stand-ins: every PE burns the runtime the instance recorded for it,')
    add('divided across the items that pass through, so the workflow reproduces the')
    add("source workflow's structure and cost without its science.")
    add('')
    add('Run it:')
    add('')
    add(f'    dispel4py timed_multi {{this file}} -i 1 -n {total_processes} \\')
    add(f'        --timing-dir monitoring_out --timing-prefix monitor')
    add('')
    add(f'-n must be at least {total_processes}: one process per source PE plus every')
    add("other PE's numprocesses. Use the `simple` mapping to debug -- the multi")
    add('mappings hide worker tracebacks and hang instead of reporting them.')
    add('"""')
    add('')
    add('import time')
    add('')
    add('from dispel4py.base import ConsumerPE, GenericPE, IterativePE, ProducerPE')
    add('from dispel4py.workflow_graph import WorkflowGraph')
    add('')
    add(f'ITEMS = {items}  # items the source emits per iteration')
    add('')

    # --- PE classes ----------------------------------------------------------
    class_of: Dict[str, str] = {}
    for pe_id in pes:
        name, _ = _pe_name_and_index(pe_id)
        stats = summary.get(pe_id, {"count": 1, "runtime": 0.0})
        in_ports, out_ports = inputs.get(pe_id, []), outputs.get(pe_id, [])
        # The recorded runtime is one instance's whole share of the run; spread it
        # over the items that instance will see.
        per_item = stats["runtime"] * stats["count"] / items if items else 0.0
        class_of[pe_id] = name

        if not in_ports:                                    # source
            add(f'class {name}(ProducerPE):')
            add(f'    """Source: emits ITEMS items on "{out_ports[0]}"."""')
            add('')
            add('    def _process(self, inputs):')
            add(f'        for index in range(ITEMS):')
            add(f'            time.sleep({per_item:.8f})')
            add(f'            self.write("{out_ports[0]}", {{"item": index}})')
        elif not out_ports and in_ports == ["input"]:       # sink on the default port
            add(f'class {name}(ConsumerPE):')
            add('    """Sink: consumes items on "input"."""')
            add('')
            add('    def _process(self, data):')
            add(f'        time.sleep({per_item:.8f})')
        elif not out_ports:                                 # sink on named ports
            add(f'class {name}(GenericPE):')
            add(f'    """Sink: consumes items on {in_ports}."""')
            add('')
            add('    def __init__(self):')
            add('        super().__init__()')
            for port in in_ports:
                add(f'        self._add_input("{port}")')
            add('')
            add('    def process(self, inputs):')
            add(f'        time.sleep({per_item:.8f})')
        elif in_ports == ["input"] and out_ports == ["output"]:
            add(f'class {name}(IterativePE):')
            add('    """One input, one output."""')
            add('')
            add('    def _process(self, data):')
            add(f'        time.sleep({per_item:.8f})')
            add('        return data')
        else:                                               # branch or join
            add(f'class {name}(GenericPE):')
            add(f'    """Ports from the trace: in={in_ports}, out={out_ports}."""')
            add('')
            add('    def __init__(self):')
            add('        super().__init__()')
            for port in in_ports:
                add(f'        self._add_input("{port}")')
            for port in out_ports:
                add(f'        self._add_output("{port}")')
            if len(out_ports) > 1:
                add('        self._next_port = 0')
            add('')
            add('    def process(self, inputs):')
            add(f'        time.sleep({per_item:.8f})')
            add('        data = next(iter(inputs.values()))')
            if len(out_ports) > 1:
                add(f'        ports = {out_ports!r}')
                add('        # the real PE routes by content; alternate so every')
                add('        # downstream branch keeps receiving work')
                add('        port = ports[self._next_port % len(ports)]')
                add('        self._next_port += 1')
                add('        self.write(port, data)')
            else:
                add(f'        self.write("{out_ports[0]}", data)')
        add('')

    # --- graph ---------------------------------------------------------------
    for pe_id in pes:
        add(f'{_identifier(pe_id)} = {class_of[pe_id]}()')
    add('')
    add('# A PE id is its name plus the order it joined the graph, and that id is what')
    add('# the monitoring trace records -- so name the PEs without their index and add')
    add('# them in index order to reproduce the ids in the source instance.')
    for pe_id in pes:
        name, _ = _pe_name_and_index(pe_id)
        add(f'{_identifier(pe_id)}.name = "{name}"')
    add('')
    for pe_id in pes:
        count = 1 if pe_id in sources else summary.get(pe_id, {}).get("count", 1)
        comment = "  # a source PE always runs in one process" if pe_id in sources else ""
        add(f'{_identifier(pe_id)}.numprocesses = {count}{comment}')
    add('')
    add('graph = WorkflowGraph()')
    add('')
    add('# added before connecting, so each pe.id matches the source instance')
    for pe_id in pes:
        add(f'graph.add({_identifier(pe_id)})')
    add('')
    for edge in shape["edges"]:
        add(f'graph.connect({_identifier(edge["from"])}, "{edge["from_connection"]}", '
            f'{_identifier(edge["to"])}, "{edge["to_connection"]}")')
    add('')
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("instance", type=pathlib.Path,
                        help="WfFormat instance to render (real or synthetic)")
    parser.add_argument("--shape", required=True, type=pathlib.Path,
                        help="monitoring directory of the same workflow, for the "
                             "PE wiring and its connection names")
    parser.add_argument("-o", "--out", type=pathlib.Path, default=None,
                        help="file to write (default: <instance stem>.py)")
    parser.add_argument("-i", "--items", type=int, default=100,
                        help="items the source emits per iteration (default 100)")
    parser.add_argument("--prefix", default="monitor",
                        help="dispel4py --timing-prefix used for the trace")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format="%(levelname)s %(message)s")

    out = args.out or args.instance.with_suffix(".py")
    out.write_text(render(args.instance, args.shape, items=args.items, prefix=args.prefix))
    print(out)


if __name__ == "__main__":
    main()

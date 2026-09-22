#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert dispel4py monitoring traces (the ``monitor_*`` artifacts produced by the
``timed_simple`` / ``timed_multi`` / ``timed_mpi`` mappings) into WfFormat
workflow instances consumable by WfCommons / WfChef.

Granularity: one WfFormat task per PE *instance* (``pe_id@rank``).

dispel4py is a streaming system: PEs exchange in-memory data over named
connections rather than files. The first PE still reads a real file and the last one 
still writes one. Those paths never reach the trace (they arrive as dispel4py root 
inputs or module constants), so name them with 
``--input-file`` / ``--output-file``; their sizes are read off disk.

Inputs read from a monitoring directory:
  monitor_concrete_shape_run<id>.json  instance-level DAG (nodes + edges)
  monitor_shape_run<id>.json           abstract PE-level DAG (fallback for edges)
  monitor_instances_run<id>.csv        per-instance runtimes
"""

from __future__ import annotations

import csv
import json
import logging
import pathlib
import re
from typing import Any, Dict, List, Optional, Tuple

from wfcommons.common.file import File
from wfcommons.common.task import Task, TaskType
from wfcommons.common.workflow import Workflow

logger = logging.getLogger(__name__)

# Connection names to assume when a trace predates the *_connection fields.
_DEFAULT_OUT_CONNECTION = "output"
_DEFAULT_IN_CONNECTION = "input"

# An instance edge, carrying the connection (stream) names at both ends:
# (source instance, destination instance, source port, destination port).
InstanceEdge = Tuple[str, str, str, str]

# dispel4py writes run ids as a compact ISO-ish stamp, e.g. 20260916T183613549171Z
_RUN_ID_RE = re.compile(r"_run(?P<run_id>[^.]+)\.(?:json|csv|png)$")


def _find_run_ids(monitoring_dir: pathlib.Path, prefix: str) -> List[str]:
    """Return every run id present in a monitoring directory, newest last."""
    run_ids = set()
    for path in monitoring_dir.glob(f"{prefix}_instances_run*.csv"):
        match = _RUN_ID_RE.search(path.name)
        if match:
            run_ids.add(match.group("run_id"))
    return sorted(run_ids)


def _one(monitoring_dir: pathlib.Path, pattern: str) -> Optional[pathlib.Path]:
    """Return the single file matching a glob, or None."""
    matches = sorted(monitoring_dir.glob(pattern))
    if not matches:
        return None
    if len(matches) > 1:
        logger.warning("multiple matches for %s, using %s", pattern, matches[0].name)
    return matches[0]


def _read_instances(path: pathlib.Path) -> Dict[str, Dict[str, Any]]:
    """
    Parse monitor_instances_run<id>.csv.

    Columns: pe_id, rank, instance_id, total_count, total_secs, avg_secs,
             min_secs, p50_secs, p95_secs, max_secs
    """
    rows: Dict[str, Dict[str, Any]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            instance_id = row.get("instance_id") or f"{row['pe_id']}@{row['rank']}"
            rows[instance_id] = {
                "pe_id": row["pe_id"],
                "rank": row["rank"],
                "instance_id": instance_id,
                "total_count": int(float(row.get("total_count") or 0)),
                "total_secs": float(row.get("total_secs") or 0.0),
                "avg_secs": float(row.get("avg_secs") or 0.0),
                "max_secs": float(row.get("max_secs") or 0.0),
            }
    return rows


def _connections(edge: Dict[str, Any]) -> Tuple[str, str]:
    """Return an edge's (source port, destination port) connection names."""
    return (edge.get("from_connection") or _DEFAULT_OUT_CONNECTION,
            edge.get("to_connection") or _DEFAULT_IN_CONNECTION)


def _concrete_edges(concrete: Dict[str, Any]) -> List[InstanceEdge]:
    """Read instance-level edges, with their port names, off the concrete shape."""
    edges = []
    for edge in concrete.get("edges", []):
        from_connection, to_connection = _connections(edge)
        edges.append((edge["from"], edge["to"], from_connection, to_connection))
    return sorted(set(edges))


def _expand_abstract_edges(
    abstract: Dict[str, Any],
    instances_by_pe: Dict[str, List[str]],
) -> List[InstanceEdge]:
    """
    Fall back to the abstract shape: connect every instance of the source PE to
    every instance of the destination PE (dispel4py's default all-to-all
    grouping). Used only when the concrete shape carries no edges. Port names
    are PE-level, so they survive the expansion unchanged.
    """
    edges = []
    for edge in abstract.get("edges", []):
        from_connection, to_connection = _connections(edge)
        for src in instances_by_pe.get(edge["from"], []):
            for dst in instances_by_pe.get(edge["to"], []):
                edges.append((src, dst, from_connection, to_connection))
    return sorted(set(edges))


def _data_file(name: str, monitoring_dir: pathlib.Path) -> File:
    """
    Resolve a real on-disk file the workflow reads or writes.

    dispel4py never records these in the trace -- the path reaches the PE as a
    root input or a module constant -- so the caller names them. Sizes are read
    off disk, looking beside the monitoring directory as well, which is where
    dispel4py leaves them when it runs from the workflow's directory.
    """
    given = pathlib.Path(name)
    for candidate in (given, monitoring_dir / given, monitoring_dir.parent / given):
        if candidate.is_file():
            return File(file_id=given.name, size=candidate.stat().st_size)
    logger.warning("%s not found on disk; recording it with size 0", name)
    return File(file_id=given.name, size=0)


def _stream_id(instance_id: str, connection: str) -> str:
    """
    Name the in-memory stream a PE instance writes to one of its output ports.

    WfFormat file ids must match ``^[0-9a-zA-Z-_./:#]*$``, which excludes the
    "@" dispel4py uses between PE and rank, so ``read0@0`` becomes ``read0:0``.
    """
    return f"{instance_id.replace('@', ':')}.{connection}"


def _build_streams(edges: List[InstanceEdge]) -> Dict[str, Dict[str, Any]]:
    """
    Map each stream id to its single producing instance and its consumers.

    One stream per (producer, output port): a fan-out across the destination
    PE's ranks is one stream with several consumers, not one stream per rank.
    Keeping a single producer per id also matches what WfCommons' translators
    assume -- they key their file maps by file id, so a second producer would
    silently overwrite the first.
    """
    streams: Dict[str, Dict[str, Any]] = {}
    for src, dst, from_connection, _ in edges:
        stream = streams.setdefault(
            _stream_id(src, from_connection), {"producer": src, "consumers": []}
        )
        if dst not in stream["consumers"]:
            stream["consumers"].append(dst)
    return streams


def _break_cycles(edges: List[Tuple[str, str]], order: List[str]) -> List[Tuple[str, str]]:
    """
    WfFormat requires a DAG, but dispel4py graphs may contain feedback loops.
    Drop back-edges relative to the recorded topological order (or, absent one,
    relative to first-seen order) and report what was removed.
    """
    position = {node: index for index, node in enumerate(order)}
    if not position:
        seen: List[str] = []
        for src, dst in edges:
            for node in (src, dst):
                if node not in seen:
                    seen.append(node)
        position = {node: index for index, node in enumerate(seen)}

    kept, dropped = [], []
    for src, dst in edges:
        if position.get(src, 0) < position.get(dst, 0):
            kept.append((src, dst))
        else:
            dropped.append((src, dst))
    if dropped:
        logger.warning(
            "dropped %d cyclic edge(s) to keep the instance graph acyclic: %s",
            len(dropped),
            dropped,
        )
    return kept


def build_workflow(
    monitoring_dir: pathlib.Path | str,
    workflow_name: Optional[str] = None,
    prefix: str = "monitor",
    run_id: Optional[str] = None,
    mapping: Optional[str] = None,
    input_files: Optional[List[str]] = None,
    output_files: Optional[List[str]] = None,
) -> Workflow:
    """
    Build a WfFormat Workflow from one dispel4py monitoring run.

    :param monitoring_dir: directory holding the monitor_* artifacts.
    :param workflow_name: name for the instance (defaults to the directory name).
    :param prefix: dispel4py --timing-prefix used for the run (default "monitor").
    :param run_id: which run to convert; defaults to the only/latest one present.
    :param mapping: dispel4py mapping that produced the trace, recorded as the
                    runtime system version (e.g. "timed_multi").
    :param input_files: real files the workflow reads (e.g.
                        "sensor_data_parallel_100.json"), attached to every
                        source task. The trace does not record them: the path
                        reaches the reading PE as a dispel4py root input.
    :param output_files: real files the workflow writes (e.g.
                         "agentic_parallel_results.jsonl"), attached to every
                         sink task.
    """
    monitoring_dir = pathlib.Path(monitoring_dir)
    if not monitoring_dir.is_dir():
        raise NotADirectoryError(f"not a directory: {monitoring_dir}")

    if run_id is None:
        run_ids = _find_run_ids(monitoring_dir, prefix)
        if not run_ids:
            raise FileNotFoundError(
                f"no {prefix}_instances_run*.csv found in {monitoring_dir}"
            )
        if len(run_ids) > 1:
            logger.warning("found %d runs, converting the latest (%s)", len(run_ids), run_ids[-1])
        run_id = run_ids[-1]

    instances_csv = monitoring_dir / f"{prefix}_instances_run{run_id}.csv"
    if not instances_csv.exists():
        raise FileNotFoundError(instances_csv)
    instances = _read_instances(instances_csv)
    if not instances:
        raise ValueError(f"{instances_csv.name} contains no instance rows")

    concrete_path = _one(monitoring_dir, f"{prefix}_concrete_shape_run{run_id}.json")
    abstract_path = _one(monitoring_dir, f"{prefix}_shape_run{run_id}.json")

    concrete = json.loads(concrete_path.read_text()) if concrete_path else {}
    abstract = json.loads(abstract_path.read_text()) if abstract_path else {}

    # --- nodes -------------------------------------------------------------
    # Prefer the concrete shape's node table; fall back to whatever the
    # instances CSV recorded (always instance-level either way).
    node_rows = concrete.get("nodes") or [
        {"instance_id": k, "pe_id": v["pe_id"], "rank": v["rank"]}
        for k, v in sorted(instances.items())
    ]

    instances_by_pe: Dict[str, List[str]] = {}
    for node in node_rows:
        instances_by_pe.setdefault(node["pe_id"], []).append(node["instance_id"])

    # --- edges -------------------------------------------------------------
    instance_edges = _concrete_edges(concrete)
    if instance_edges:
        edge_source = "concrete shape"
    else:
        instance_edges = _expand_abstract_edges(abstract, instances_by_pe)
        edge_source = "abstract shape (expanded across ranks)"
    logger.info("instance edges derived from %s", edge_source)

    edges = _break_cycles(
        sorted({(src, dst) for src, dst, _, _ in instance_edges}),
        concrete.get("topological_order", []),
    )
    # Streams follow the dependencies that survived cycle breaking, so a dropped
    # back-edge takes its stream with it.
    kept = set(edges)
    streams = _build_streams(
        [edge for edge in instance_edges if (edge[0], edge[1]) in kept]
    )
    logger.info("%d in-memory stream(s) named from connection names", len(streams))

    # Streams are shared objects: one File per id, so the producer and every
    # consumer reference the same entry in the workflow's file table.
    stream_files = {
        stream_id: File(file_id=stream_id, size=0) for stream_id in streams
    }
    outputs_of: Dict[str, List[File]] = {}
    inputs_of: Dict[str, List[File]] = {}
    for stream_id, stream in sorted(streams.items()):
        outputs_of.setdefault(stream["producer"], []).append(stream_files[stream_id])
        for consumer in stream["consumers"]:
            inputs_of.setdefault(consumer, []).append(stream_files[stream_id])

    # --- real files --------------------------------------------------------
    # Streaming is only the middle of the workflow: the first PE still reads a
    # file off disk and the last one still writes one. Sources and sinks are
    # whatever the streams left unconnected, computed before the real files are
    # attached so they do not mask each other.
    all_instances = [node["instance_id"] for node in node_rows]
    sources = [i for i in all_instances if not inputs_of.get(i)]
    sinks = [i for i in all_instances if not outputs_of.get(i)]
    for name in input_files or []:
        data_file = _data_file(name, monitoring_dir)
        logger.info("%s (%d bytes) read by %s",
                    data_file.file_id, data_file.size, ", ".join(sources))
        for instance_id in sources:
            inputs_of.setdefault(instance_id, []).append(data_file)
    for name in output_files or []:
        data_file = _data_file(name, monitoring_dir)
        logger.info("%s (%d bytes) written by %s",
                    data_file.file_id, data_file.size, ", ".join(sinks))
        for instance_id in sinks:
            outputs_of.setdefault(instance_id, []).append(data_file)

    # --- workflow ----------------------------------------------------------
    workflow = Workflow(
        name=workflow_name or monitoring_dir.name,
        description=(
            f"dispel4py execution trace (run {run_id}) converted to WfFormat; "
            f"one task per PE instance"
        ),
        runtime_system_name="dispel4py",
        runtime_system_version=mapping or "unknown",
        runtime_system_url="https://github.com/StreamingFlow/d4py",
    )

    # Stable, deterministic numbering so repeated conversions are diffable.
    ordered_ids = [n["instance_id"] for n in sorted(node_rows, key=lambda n: n["instance_id"])]
    # wfchef parses the task type out of the id via id.split("_ID"), so the id
    # must be "<pe_id>_ID<n>" for the PE to be recognised as the task type.
    task_ids = {
        instance_id: f"{_pe_of(instance_id, node_rows)}_ID{index:07d}"
        for index, instance_id in enumerate(ordered_ids)
    }

    makespan = 0.0
    for instance_id in ordered_ids:
        stats = instances.get(instance_id)
        if stats is None:
            logger.warning("%s appears in the shape but not in %s; runtime set to 0",
                           instance_id, instances_csv.name)
            stats = {"total_secs": 0.0, "total_count": 0, "avg_secs": 0.0, "max_secs": 0.0}
        pe_id = _pe_of(instance_id, node_rows)
        runtime = stats["total_secs"]
        makespan = max(makespan, runtime)
        # NOTE: name must equal task_id. WfFormat records dependencies as task
        # ids, but wfchef's create_graph_from_json_object keys its nodes by
        # task["name"] and then draws edges from those id-valued parents --
        # any mismatch creates attribute-less phantom nodes and annotate()
        # dies with KeyError: 'id'. The dispel4py instance id is preserved as
        # the command argument instead.
        workflow.add_task(
            Task(
                name=task_ids[instance_id],
                task_id=task_ids[instance_id],
                runtime=runtime,
                cores=1.0,
                category=pe_id,
                program=pe_id,
                args=[instance_id],
                task_type=TaskType.COMPUTE,
                input_files=inputs_of.get(instance_id, []),
                output_files=outputs_of.get(instance_id, []),
            )
        )

    for src, dst in edges:
        if src in task_ids and dst in task_ids:
            workflow.add_dependency(task_ids[src], task_ids[dst])

    workflow.makespan = makespan
    return workflow


def _pe_of(instance_id: str, node_rows: List[Dict[str, Any]]) -> str:
    for node in node_rows:
        if node["instance_id"] == instance_id:
            return node["pe_id"]
    return instance_id.split("@")[0]


def _files_for(
    spec: Optional[List[str] | Dict[str, List[str]]],
    monitoring_dir: pathlib.Path,
) -> List[str]:
    """
    Pick one directory's real input (or output) files out of ``spec``.

    A plain list applies to every directory being converted; a dict keyed by
    directory name or path lets runs of the same workflow over different data
    declare their own (monitoring_simple reads sensor_data_agentic.json while
    monitoring_multi reads sensor_data_parallel_100.json).
    """
    if not spec:
        return []
    if isinstance(spec, dict):
        for key in (monitoring_dir.name, str(monitoring_dir)):
            if key in spec:
                return spec[key]
        return []
    return spec


def convert(
    monitoring_dirs: List[pathlib.Path | str],
    output_dir: pathlib.Path | str,
    workflow_name: Optional[str] = None,
    prefix: str = "monitor",
    input_files: Optional[List[str] | Dict[str, List[str]]] = None,
    output_files: Optional[List[str] | Dict[str, List[str]]] = None,
) -> List[pathlib.Path]:
    """
    Convert one or more monitoring directories into WfFormat JSON files.

    :param input_files: real files the workflows read, either as a list applied
                        to every directory or as a {directory: [files]} dict.
    :param output_files: real files the workflows write, same forms.
    """
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    written = []
    for monitoring_dir in monitoring_dirs:
        monitoring_dir = pathlib.Path(monitoring_dir)
        mapping = monitoring_dir.name.replace("monitoring_", "") or None
        workflow = build_workflow(
            monitoring_dir,
            workflow_name=workflow_name or monitoring_dir.name,
            prefix=prefix,
            mapping=mapping,
            input_files=_files_for(input_files, monitoring_dir),
            output_files=_files_for(output_files, monitoring_dir),
        )
        # wfchef groups instances by size, so encode the task count in the name.
        out = output_dir / f"{workflow.name}-{len(workflow.tasks)}.json"
        workflow.write_json(out)
        written.append(out)
        logger.info("wrote %s (%d tasks, %d edges)",
                    out.name, len(workflow.tasks), workflow.number_of_edges())
    return written


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("monitoring_dirs", nargs="+", type=pathlib.Path,
                        help="dispel4py monitoring directories (e.g. timings/ or monitoring_multi/)")
    parser.add_argument("-o", "--out", required=True, type=pathlib.Path,
                        help="directory to write WfFormat JSON instances into")
    parser.add_argument("-n", "--name", default=None,
                        help="workflow name (defaults to each directory's name)")
    parser.add_argument("--prefix", default="monitor",
                        help="dispel4py --timing-prefix used for the run")
    parser.add_argument("-i", "--input-file", action="append", dest="input_files",
                        metavar="FILE",
                        help="real file the workflow reads, attached to every source "
                             "task (repeatable). Not recorded in the trace, so it has "
                             "to be named here, e.g. sensor_data_parallel_100.json")
    parser.add_argument("--output-file", action="append", dest="output_files",
                        metavar="FILE",
                        help="real file the workflow writes, attached to every sink "
                             "task (repeatable), e.g. agentic_parallel_results.jsonl")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s %(message)s",
    )
    for path in convert(args.monitoring_dirs, args.out, args.name, args.prefix,
                        input_files=args.input_files,
                        output_files=args.output_files):
        print(path)


if __name__ == "__main__":
    main()
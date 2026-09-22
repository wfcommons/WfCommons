#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Merge newly collected monitoring runs into an existing WfFormat corpus.

Unlike `convert_traces`, this never clears WFFORMAT_DIR: it converts only the
directories that have no instance there yet, so the recipe can be re-cooked over
a growing set of traces. Use --force to reconvert a directory anyway.

    python -m wfcommons.wfstream.update_traces monitoring_multi_64
"""

import json
import pathlib

from .config import INPUT_FILES, OUTPUT_FILES, WFFORMAT_DIR, trace_dirs
from .convert_traces import describe
from .dispel_fwd_converter import convert


def read_traces(traces_path: pathlib.Path = None) -> list:
    """Load every converted WfFormat instance under a directory."""
    traces_path = traces_path or WFFORMAT_DIR
    return [json.loads(p.read_text()) for p in sorted(traces_path.glob("**/*.json"))]


def existing(output_dir: pathlib.Path = None) -> dict:
    """Map each monitoring directory name already converted to its instances.

    Instances are written as ``<dir name>-<task count>.json``, so the stem's
    prefix identifies the trace it came from.
    """
    output_dir = output_dir or WFFORMAT_DIR
    found = {}
    for path in sorted(output_dir.glob("*.json")):
        found.setdefault(path.stem.rsplit("-", 1)[0], []).append(path)
    return found


def merge_traces(dirs,
                 output_dir: pathlib.Path = None,
                 force: bool = False,
                 input_files=None,
                 output_files=None) -> tuple:
    """Convert the monitoring directories missing from the corpus.

    :param input_files: real files each workflow reads, as a list or a
        {directory name: [files]} dict; defaults to the config maps.
    :return: (written paths, skipped directory names)
    """
    output_dir = output_dir or WFFORMAT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)
    already = existing(output_dir)

    todo, skipped = [], []
    for monitoring_dir in (pathlib.Path(d) for d in dirs):
        if not force and monitoring_dir.name in already:
            skipped.append(monitoring_dir.name)
        else:
            todo.append(monitoring_dir)

    if not todo:
        return [], skipped
    written = convert(todo,
                      output_dir=output_dir,
                      input_files=INPUT_FILES if input_files is None else input_files,
                      output_files=OUTPUT_FILES if output_files is None else output_files)
    return written, skipped


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("monitoring_dirs", nargs="*", type=pathlib.Path,
                        help="monitoring directories (default: config.TRACE_DIRS)")
    parser.add_argument("-o", "--out", type=pathlib.Path, default=WFFORMAT_DIR)
    parser.add_argument("-f", "--force", action="store_true",
                        help="reconvert directories already in the corpus")
    args = parser.parse_args()

    dirs = args.monitoring_dirs or trace_dirs()
    for name in (d.name for d in dirs if d.name not in INPUT_FILES):
        print(f"  warning: no input/output files configured for {name}")

    written, skipped = merge_traces(dirs, args.out, force=args.force)
    for name in skipped:
        print(f"  skipped {name} (already converted)")
    for path in written:
        print(f"  added {describe(path)}")
    print(f"corpus now holds {len(list(args.out.glob('*.json')))} instance(s)")


if __name__ == "__main__":
    main()

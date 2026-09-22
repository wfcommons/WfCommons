#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 1: dispel4py monitoring directories -> WfFormat instances.

Rebuilds WFFORMAT_DIR from scratch. WfChef cooks over *every* JSON in that
directory, so a leftover instance from an earlier naming scheme would silently
join the recipe; clearing first avoids that.

    python -m wfcommons.wfstream.convert_traces [monitoring_dir ...]
"""

import json
import pathlib

from .config import INPUT_FILES, OUTPUT_FILES, WFFORMAT_DIR, trace_dirs
from .dispel_fwd_converter import convert


def describe(path: pathlib.Path) -> str:
    """One-line summary of a converted instance: tasks, streams, real files."""
    spec = json.loads(path.read_text())["workflow"]["specification"]
    files = spec.get("files", [])
    streams = [f for f in files if ":" in f["id"]]
    reals = [f["id"] for f in files if ":" not in f["id"]]
    return (f"{path.name}: {len(spec['tasks'])} tasks, "
            f"{len(streams)} streams, files={reals}")


def convert_traces(dirs=None,
                   output_dir: pathlib.Path = None,
                   clean: bool = True,
                   input_files=None,
                   output_files=None) -> list:
    """Convert monitoring directories into WfFormat instances.

    :param input_files: real files each workflow reads, as a list applied to
        every directory or a {directory name: [files]} dict. Defaults to the
        config maps; a caller supplying its own traces supplies its own maps.
    """
    dirs = dirs or trace_dirs()
    output_dir = output_dir or WFFORMAT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    if clean:
        for old in sorted(output_dir.glob("*.json")):
            old.unlink()

    return convert(dirs,
                   output_dir=output_dir,
                   input_files=INPUT_FILES if input_files is None else input_files,
                   output_files=OUTPUT_FILES if output_files is None else output_files)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("monitoring_dirs", nargs="*", type=pathlib.Path,
                        help="monitoring directories (default: config.TRACE_DIRS)")
    parser.add_argument("-o", "--out", type=pathlib.Path, default=WFFORMAT_DIR)
    parser.add_argument("--keep", action="store_true",
                        help="do not clear the output directory first")
    args = parser.parse_args()

    dirs = args.monitoring_dirs or trace_dirs()
    print(f"converting {len(dirs)} trace(s) -> {args.out}")
    for path in convert_traces(dirs, args.out, clean=not args.keep):
        print(f"  {describe(path)}")


if __name__ == "__main__":
    main()

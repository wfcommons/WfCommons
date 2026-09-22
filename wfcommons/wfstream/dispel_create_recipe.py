#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Rebuild the climate WfChef recipe from dispel4py monitoring traces.

Runs the three steps in order; each one is also a script in its own right:

  1. convert_traces      monitoring dirs -> WfFormat instances (WFFORMAT_DIR)
  2. build_recipe        cook the recipe, optionally install it
  3. generate_workflows  write synthetic instances (SYNTHETIC_DIR) and check them

Paths and per-trace settings live in config.py. Rendering instances back into a
runnable dispel4py script is not wired up here -- the reverse converter is out of
the tree for now.

    python -m wfcommons.wfstream.dispel_create_recipe --install

It is idempotent: each run reconverts and rebuilds from scratch.
"""

from . import build_recipe, convert_traces, generate_workflows
from .config import GENERATE_SIZES, RECIPE_DIR, RECIPE_NAME, WFFORMAT_DIR


def run(install: bool = False, sizes=None) -> None:
    print(f"[1/3] converting traces -> {WFFORMAT_DIR}")
    for path in convert_traces.convert_traces():
        print(f"  {convert_traces.describe(path)}")

    print(f"[2/3] cooking recipe {RECIPE_NAME!r}")
    cooked = build_recipe.cook()
    build_recipe.summarize(cooked)

    if not install:
        print("[3/3] skipped: pass --install to install the recipe and generate")
        return

    print(f"  installing into {RECIPE_DIR}")
    build_recipe.install(cooked)
    generate_workflows.verify()

    print(f"[3/3] generating instances of sizes {tuple(sizes or GENERATE_SIZES)}")
    generate_workflows.generate(sizes)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--install", action="store_true",
                        help="install the cooked recipe and generate instances")
    parser.add_argument("sizes", nargs="*", type=int,
                        help=f"task counts to generate (default: {GENERATE_SIZES})")
    args = parser.parse_args()
    run(install=args.install, sizes=args.sizes)


if __name__ == "__main__":
    main()

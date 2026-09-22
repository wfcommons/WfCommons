#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 2: cook a WfChef recipe from the converted instances, and install it.

Cooking finds the microstructures shared across the instances in WFFORMAT_DIR and
writes a recipe package under BUILD_DIR. Installing copies that package's data
over the recipe inside the wfcommons package, so imports pick it up.

    python -m wfcommons.wfstream.build_recipe --install
"""

import json
import pathlib
import shutil

from wfcommons.wfchef.chef import create_recipe

from .config import BUILD_DIR, RECIPE_DIR, RECIPE_NAME, WFFORMAT_DIR


def cook(instances_dir: pathlib.Path = None,
         build_dir: pathlib.Path = None,
         name: str = RECIPE_NAME) -> pathlib.Path:
    """Build a recipe package from the instances; returns the package path."""
    instances_dir = instances_dir or WFFORMAT_DIR
    build_dir = build_dir or BUILD_DIR
    if build_dir.exists():
        shutil.rmtree(build_dir)
    create_recipe(instances_dir, build_dir, name)
    return build_dir / f"wfchef_recipe_{name}"


def summarize(cooked: pathlib.Path) -> None:
    """Print the microstructures found per base graph, plus the error table."""
    micro = cooked / "microstructures"
    for base in sorted(d for d in micro.iterdir() if d.is_dir() and d.name != "metric"):
        found = json.loads((base / "microstructures.json").read_text())
        shapes = ", ".join(f"{v['nodes'][0][0].split('_ID')[0]}x{len(v['nodes'])}"
                           for v in found.values())
        print(f"  {base.name}: {len(found)} microstructure(s) {shapes}")
    table = (micro / "metric" / "err.csv").read_text().strip()
    print("  error table:\n    " + table.replace("\n", "\n    "))


def install(cooked: pathlib.Path, recipe_dir: pathlib.Path = None) -> None:
    """Copy the cooked recipe's data over the installed recipe.

    Only the data moves. recipe.py / __init__.py are left alone: the generated
    ones define the same class, and the installed copy may carry local edits.
    """
    recipe_dir = recipe_dir or RECIPE_DIR
    if not recipe_dir.exists():
        raise SystemExit(f"no installed recipe at {recipe_dir}; copy the whole "
                         f"package from {cooked} instead")
    shutil.rmtree(recipe_dir / "microstructures", ignore_errors=True)
    shutil.copytree(cooked / "microstructures", recipe_dir / "microstructures")
    shutil.copy2(cooked / "task_type_stats.json", recipe_dir / "task_type_stats.json")


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("-i", "--instances", type=pathlib.Path, default=WFFORMAT_DIR)
    parser.add_argument("-b", "--build-dir", type=pathlib.Path, default=BUILD_DIR)
    parser.add_argument("-n", "--name", default=RECIPE_NAME)
    parser.add_argument("--install", action="store_true",
                        help=f"copy the result into {RECIPE_DIR}")
    args = parser.parse_args()

    print(f"cooking recipe {args.name!r} -> {args.build_dir}")
    cooked = cook(args.instances, args.build_dir, args.name)
    summarize(cooked)
    if args.install:
        print(f"installing into {RECIPE_DIR}")
        install(cooked)


if __name__ == "__main__":
    main()

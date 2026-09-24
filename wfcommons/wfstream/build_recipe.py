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

from wfcommons.wfchef.chef import create_recipe, get_recipe, install_recipe

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


def register(build_dir: pathlib.Path = None, verbose: bool = False) -> str:
    """Install the cooked package so the recipe becomes a WfChef entry point.

    `create_recipe` writes a pyproject.toml at the build directory's root that
    declares ``<name>_recipe``; until something pip-installs it, the recipe is
    invisible to ``wfchef ls`` and to `get_recipe`, which is how everything else
    in WfCommons finds a recipe. Copying the data into the installed package (see
    `install`) makes the recipe usable but never registers that name.

    Registering is what makes a wfstream recipe behave like a built-in one: it
    survives into a fresh environment, so a notebook that pip-installs WfCommons
    and then this package gets the recipe without knowing where it came from.

    :return: the entry point name now registered.
    """
    build_dir = pathlib.Path(build_dir or BUILD_DIR)
    if not build_dir.joinpath("pyproject.toml").exists():
        raise SystemExit(f"{build_dir} has no pyproject.toml; cook the recipe "
                         f"first so create_recipe writes one")
    install_recipe(build_dir, verbose=verbose)
    return f"{_package_name(build_dir)}_recipe"


def _package_name(build_dir: pathlib.Path) -> str:
    """The recipe name a cooked package declares."""
    for child in sorted(pathlib.Path(build_dir).iterdir()):
        if child.is_dir() and child.name.startswith("wfchef_recipe_"):
            return child.name[len("wfchef_recipe_"):]
    raise SystemExit(f"no wfchef_recipe_* package under {build_dir}")


def registered(name: str = RECIPE_NAME):
    """The recipe class WfChef resolves for ``name``, or None if unregistered."""
    return get_recipe(f"{name}_recipe")


def install(cooked: pathlib.Path,
            recipe_dir: pathlib.Path = None,
            name: str = None) -> pathlib.Path:
    """Copy the cooked recipe's data over the installed recipe.

    The target is wherever the recipe actually resolves from -- site-packages if
    it was registered, the wfcommons tree otherwise. Writing to a fixed path
    instead would leave two copies of one recipe drifting apart, with whichever
    one `load_recipe` finds deciding what gets generated.

    Only the data moves. recipe.py / __init__.py are left alone: the generated
    ones define the same class, and the installed copy may carry local edits.

    :return: the directory written to.
    """
    if recipe_dir is None:
        from .streaming_recipe import recipe_path
        try:
            recipe_dir = recipe_path(name or RECIPE_NAME)
        except SystemExit:
            recipe_dir = RECIPE_DIR          # nothing installed yet
    if not recipe_dir.exists():
        raise SystemExit(f"no installed recipe at {recipe_dir}; register the "
                         f"cooked package instead: build_recipe.register()")
    shutil.rmtree(recipe_dir / "microstructures", ignore_errors=True)
    shutil.copytree(cooked / "microstructures", recipe_dir / "microstructures")
    shutil.copy2(cooked / "task_type_stats.json", recipe_dir / "task_type_stats.json")
    return recipe_dir


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("-i", "--instances", type=pathlib.Path, default=WFFORMAT_DIR)
    parser.add_argument("-b", "--build-dir", type=pathlib.Path, default=BUILD_DIR)
    parser.add_argument("-n", "--name", default=RECIPE_NAME)
    parser.add_argument("--install", action="store_true",
                        help=f"copy the result into {RECIPE_DIR}")
    parser.add_argument("--register", action="store_true",
                        help="pip-install the cooked package so the recipe "
                             "shows up in `wfchef ls`")
    args = parser.parse_args()

    print(f"cooking recipe {args.name!r} -> {args.build_dir}")
    cooked = cook(args.instances, args.build_dir, args.name)
    summarize(cooked)
    if args.install:
        print("installed into", install(cooked, name=args.name))
    if args.register:
        name = register(args.build_dir)
        print(f"registered as {name}; "
              f"{'resolves' if registered(args.name) else 'DOES NOT resolve'} "
              f"through wfchef")


if __name__ == "__main__":
    main()

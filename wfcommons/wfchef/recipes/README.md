# wfchef-recipe-climate

WfChef recipe for climate workflow.

## Installation

```bash
pip install -e .
```

## Usage

```python
from wfchef_recipe_climate import ClimateRecipe

recipe = ClimateRecipe()
# Use the recipe...
```

## Entry Point

This package registers the following workflow recipe entry point:
- `climate_recipe` -> `wfchef_recipe_climate:ClimateRecipe`

You can load it using:
```python
from wfcommons.wfchef import get_recipe

Recipe = get_recipe("climate_recipe")
recipe = Recipe()
```

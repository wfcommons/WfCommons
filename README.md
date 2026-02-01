[![Build][build-badge]][build-link]&nbsp;
[![PyPI version][pypi-badge]][pypi-link]&nbsp;
[![License: LGPL v3][license-badge]](LICENSE)&nbsp;
[![CodeFactor][codefactor-badge]][codefactor-link]&nbsp;
[![Documentation Status][rtd-badge]][rtd-link]&nbsp;
[![Codecov][cov-badge]][cov-link]&nbsp;
[![Downloads](https://static.pepy.tech/personalized-badge/wfcommons?period=total&units=international_system&left_color=grey&right_color=yellowgreen&left_text=Downloads)](https://pepy.tech/project/wfcommons)

<a href="https://wfcommons.org" target="_blank"><img src="https://wfcommons.org/images/wfcommons-logo-white.png" width="350" /></a>
<br/>_An open-source ecosystem of workflow execution instances, synthetic workflow generators, and benchmark specifications_

WfCommons is a Python framework that helps you:

- Analyze real workflow execution instances;
- Build workflow recipes from production traces;
- Generate realistic synthetic workflow instances; and
- Create workflow benchmark specifications.

Quick links: [Documentation](https://wfcommons.readthedocs.io/en/latest/) · [Website](https://wfcommons.org) · [GitHub Issues](https://github.com/wfcommons/wfcommons/issues)

## Quickstart

WfCommons requires Python 3.11+ and has been tested on Linux and macOS.

```
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install wfcommons
```

Generate a synthetic workflow in a few lines:

```python
import pathlib
from wfcommons.wfchef.recipes import SeismologyRecipe
from wfcommons import WorkflowGenerator

generator = WorkflowGenerator(SeismologyRecipe.from_num_tasks(250))
workflow = generator.build_workflow()
workflow.write_json(pathlib.Path("seismology-workflow.json"))
```

Next steps:

- Learn how to build recipes in the [WfChef guide](https://wfcommons.readthedocs.io/en/latest/generating_workflows_recipe.html).
- Generate larger workflow families in the [WfGen guide](https://wfcommons.readthedocs.io/en/latest/generating_workflows.html).
- Produce benchmark specs in the [WfBench guide](https://wfcommons.readthedocs.io/en/latest/generating_workflow_benchmarks.html).

## Installation

WfCommons is available on [PyPI](https://pypi.org/project/wfcommons).

```
python3 -m pip install wfcommons
```

### Installing from source (latest)

```
git clone https://github.com/wfcommons/wfcommons
cd wfcommons
python3 -m pip install .
```

### Optional Requirements

#### Graphviz
WfCommons uses _pygraphviz_ for generating visualizations for the workflow task graph.
If you want to enable this feature, you will have to install the
[graphviz](https://www.graphviz.org/) package (version 2.16 or later).
You can install graphviz easily on Linux with your favorite package manager,
for example for Debian-based distributions:
```
sudo apt-get install graphviz libgraphviz-dev
```
and for RedHat-based distributions:
```
sudo yum install python-devel graphviz-devel
```

On macOS you can use `brew` package manager:
```
brew install graphviz
```

Then you can install pygraphviz by running:
```
python3 -m pip install pygraphviz
```

#### pydot
WfCommons uses _pydot_ for reading and writing DOT files. If you want to enable
this feature, you will have to install the pydot package:
```
python3 -m pip install pydot
```

## Get in Touch

The main channel to reach the WfCommons team is via the support email:
[support@wfcommons.org](mailto:support@wfcommons.org).

**Bug Report / Feature Request:** our preferred channel to report a bug or request a feature is via
WfCommons's [Github Issues Track](https://github.com/wfcommons/wfcommons/issues).

## Citing WfCommons
When citing WfCommons, please use the following paper. You should also actually read
that paper, as it provides a recent and general overview on the framework.

```
@article{wfcommons,
    title = {{WfCommons: A Framework for Enabling Scientific Workflow Research and Development}},
    author = {Coleman, Taina and Casanova, Henri and Pottier, Loic and Kaushik, Manav and Deelman, Ewa and Ferreira da Silva, Rafael},
    journal = {Future Generation Computer Systems},
    volume = {128},
    number = {},
    pages = {16--27},
    doi = {10.1016/j.future.2021.09.043},
    year = {2022},
}
```

[build-badge]:         https://github.com/wfcommons/wfcommons/workflows/Build/badge.svg
[build-link]:          https://github.com/wfcommons/wfcommons/actions
[pypi-badge]:          https://badge.fury.io/py/wfcommons.svg
[pypi-link]:           https://badge.fury.io/py/wfcommons
[license-badge]:       https://img.shields.io/badge/License-LGPL%20v3-blue.svg
[codefactor-badge]:    https://www.codefactor.io/repository/github/wfcommons/wfcommons/badge
[codefactor-link]:     https://www.codefactor.io/repository/github/wfcommons/wfcommons
[rtd-badge]:           https://readthedocs.org/projects/wfcommons/badge/?version=latest
[rtd-link]:            https://wfcommons.readthedocs.io/en/latest/?badge=latest
[cov-badge]:           https://codecov.io/gh/wfcommons/WfCommons/graph/badge.svg?token=PJTXMLCIXD
[cov-link]:            https://codecov.io/gh/wfcommons/WfCommons

[![Build][build-badge]][build-link]
[![PyPI version][pypi-badge]][pypi-link]
[![License: LGPL v3][license-badge]](LICENSE)
[![CodeFactor][codefactor-badge]][codefactor-link]
[![Documentation Status][rtd-badge]][rtd-link]

<a href="https://wfcommons.org" target="_blank"><img src="https://wfcommons.org/images/wfcommons-horizontal.png" width="350" /></a>
<br/>_A Framework for Enabling Scientific Workflow Research and Development_

This Python package provides a collection of tools for:

- Analyzing instances of actual workflow executions;
- Producing recipes structures for creating workflow recipes for workflow generation; and
- Generating synthetic realistic workflow instances.

## Installation

WfCommons is available on [PyPI](https://pypi.org/project/wfcommons).
WfCommons requires Python3.6+ and has been tested on Linux and macOS.

### Requirements

#### Graphviz
WfCommons uses _pygraphviz_ and thus needs the graphviz package installed (version 2.16 or later).
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

### Installation using pip

While `pip` can be used to install WfCommons, we suggest the following
approach for reliable installation when many Python environments are available:

```
$ python3 -m pip install wfcommons
```

### Retrieving the latest unstable version

If you want to use the latest WfCommons unstable version, that will contain
brand new features (but also contain bugs as the stabilization work is still
underway), you may consider retrieving the latest unstable version.

Cloning from [WfCommons's GitHub](https://github.com/wfcommons/wfcommons)
repository:

```
$ git clone https://github.com/wfcommons/wfcommons
$ cd wfcommons
$ pip install .
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
@article{coleman2021wfcommons,
  title={WfCommons: A Framework for Enabling Scientific Workflow Research and Development},
  author={Coleman, Tain{\~a} and Casanova, Henri and Pottier, Lo\"ic and Kaushik, Manav and Deelman, Ewa and Ferreira da Silva, Rafael},
  journal={arXiv preprint arXiv:2105.14352},
  year={2021}
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

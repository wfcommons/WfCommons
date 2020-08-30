[![Travis][travis-badge]][travis-link]
[![PyPI version][pypi-badge]][pypi-link]
[![License: LGPL v3][license-badge]](LICENSE)
[![CodeFactor][codefactor-badge]][codefactor-link]
[![Documentation Status][rtd-badge]][rtd-link]

<a href="https://workflowhub.org" target="_blank"><img src="https://workflowhub.org/assets/images/logo-horizontal.png" width="300" /></a>
<br/>_A Community Framework for Enabling Scientific Workflow Research and Development_

This Python package provides a collection of tools for:

- Analyzing traces of actual workflow executions;
- Producing recipes structures for creating workflow recipes for workflow generation; and
- Generating synthetic realistic workflow traces.

## Installation

WorkflowHub is available on [PyPI](https://pypi.org/project/workflowhub).
WorkflowHub requires Python3.5+ and has been tested on Linux and macOS.

### Requirements

#### Graphviz
WorkflowHub uses _pygraphviz_ and thus needs the graphviz package installed (version 2.16 or later).
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

While `pip` can be used to install WorkflowHub, we suggest the following
approach for reliable installation when many Python environments are available:

```
$ python3 -m pip install workflowhub
```

### Retrieving the latest unstable version

If you want to use the latest WorkflowHub unstable version, that will contain
brand new features (but also contain bugs as the stabilization work is still
underway), you may consider retrieving the latest unstable version.

Cloning from [WorkflowHub's GitHub](https://github.com/workflowhub/workflowhub)
repository:

```
$ git clone https://github.com/workflowhub/workflowhub
$ cd workflowhub
$ pip install .
```

[travis-badge]:        https://travis-ci.org/workflowhub/workflowhub.svg?branch=master
[travis-link]:         https://travis-ci.org/workflowhub/workflowhub
[pypi-badge]:          https://badge.fury.io/py/workflowhub.svg
[pypi-link]:           https://badge.fury.io/py/workflowhub
[license-badge]:       https://img.shields.io/badge/License-LGPL%20v3-blue.svg
[codefactor-badge]:    https://www.codefactor.io/repository/github/workflowhub/workflowhub/badge
[codefactor-link]:     https://www.codefactor.io/repository/github/workflowhub/workflowhub
[rtd-badge]:           https://readthedocs.org/projects/workflowhub/badge/?version=latest
[rtd-link]:            https://workflowhub.readthedocs.io/en/latest/?badge=latest

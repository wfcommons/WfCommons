WorkflowHub
***********

|pypi-badge| |travis-badge| |license-badge|

`WorkflowHub <https://workflowhub.org>`__ is a community framework for
enabling scientific workflow research and education. This Python package
provides a collection of tools for:

- Analyzing traces of actual workflow executions;
- Producing recipes structures for creating workflow recipes for workflow
  generation; and
- Generating synthetic realistic workflow traces.

.. figure::  images/workflowhub.png
   :align:   center

   The WorkflowHub conceptual architecture.

----

Support
=======

The source code for this WorkflowHub's Python package is available on
`GitHub <http://github.com/workflowhub/workflowhub>`_.
Our preferred channel to report a bug or request a feature is via WorkflowHub's
Github `Issues Track <https://github.com/workflowhub/workflowhub/issues>`_.

You can also reach the WorkflowHub team via our support email:
support@workflowhub.org.

----

.. toctree::
    :caption: Quickstart
    :maxdepth: 2

    quickstart_installation.rst

.. toctree::
    :caption: User Guide
    :maxdepth: 2

    introduction.rst
    analyzing_traces.rst
    generating_workflows.rst

.. toctree::
    :caption: API Reference
    :maxdepth: 1

    user_api_reference.rst
    dev_api_reference.rst

.. |travis-badge| image:: https://travis-ci.org/workflowhub/workflowhub.svg?branch=master
    :target: https://travis-ci.org/workflowhub/workflowhub
.. |license-badge| image:: https://img.shields.io/badge/License-LGPL%20v3-blue.svg
    :target: https://github.com/workflowhub/workflowhub/blob/master/LICENSE
.. |pypi-badge| image:: https://badge.fury.io/py/workflowhub.svg
    :target: https://badge.fury.io/py/workflowhub

.. figure::  images/wfcommons-horizontal.png
   :scale:   15 %
   :align:   left

|pypi-badge| |build-badge| |license-badge|

`WfCommons <https://wfcommons.org>`__ is an open-source Python framework for
enabling scientific workflow research and development. Use it to analyze
workflow execution instances, build workflow recipes, generate realistic
synthetic workflows, and create benchmark specifications.

Quick links: `Website <https://wfcommons.org>`__ ·
`GitHub <https://github.com/wfcommons/wfcommons>`__

.. figure::  images/wfcommons.png
   :scale:   70 %
   :align:   center

   The WfCommons conceptual architecture.

What you can do
===============

- Analyze workflow instances from real executions.
- Derive reusable recipes that capture workflow structure and performance.
- Generate synthetic workflows at scale for experimentation.
- Build benchmark specs for multiple workflow systems.

Get started
===========

- :doc:`quickstart_installation` for install and first steps.
- :doc:`introduction` for a project overview and WfFormat context.
- :doc:`generating_workflows_recipe` to build recipes from real instances.
- :doc:`generating_workflows` to generate synthetic workflows.
- :doc:`generating_workflow_benchmarks` to produce benchmark specs.

----

Support
=======

The source code for this WfCommons's Python package is available on
`GitHub <http://github.com/wfcommons/wfcommons>`_.
Our preferred channel to report a bug or request a feature is via WfCommons's
Github `Issues Track <https://github.com/wfcommons/wfcommons/issues>`_.

You can also reach the WfCommons team via our support email:
support@wfcommons.org.

----

Contents
========

.. toctree::
    :caption: Quickstart
    :hidden:
    :maxdepth: 2

    quickstart_installation.rst

.. toctree::
    :caption: User Guide
    :hidden:
    :maxdepth: 2

    introduction.rst
    analyzing_instances.rst
    generating_workflows_recipe.rst
    generating_workflows.rst
    generating_workflow_benchmarks.rst

.. toctree::
    :caption: API Reference
    :hidden:
    :maxdepth: 1

    user_api_reference.rst
    dev_api_reference.rst

.. |build-badge| image:: https://github.com/wfcommons/wfcommons/workflows/Build/badge.svg
    :target: https://github.com/wfcommons/wfcommons/actions
.. |license-badge| image:: https://img.shields.io/badge/License-LGPL%20v3-blue.svg
    :target: https://github.com/wfcommons/wfcommons/blob/master/LICENSE
.. |pypi-badge| image:: https://badge.fury.io/py/wfcommons.svg
    :target: https://badge.fury.io/py/wfcommons

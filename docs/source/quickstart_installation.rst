Quickstart
==========

WfCommons is available on `PyPI <https://pypi.org/project/wfcommons>`_.
It requires Python 3.11+ and has been tested on Linux and macOS.

Installation
------------

We recommend installing into a virtual environment to keep dependencies
isolated:

.. code-block:: bash

    $ python3 -m venv .venv
    $ source .venv/bin/activate
    $ python3 -m pip install wfcommons

Verify the installation
-----------------------

Confirm the :code:`wfchef` CLI entry point is available:

.. code-block:: bash

    $ wfchef --help

Or check the Python import:

.. code-block:: bash

    $ python3 -c "import wfcommons; print(wfcommons.__version__)"

Your first synthetic workflow
-----------------------------

WfCommons ships with ready-to-use recipes for ten scientific applications
(see :ref:`recipes-list`). Generating a realistic workflow takes a few lines:

.. code-block:: python

    import pathlib
    from wfcommons.wfchef.recipes import SeismologyRecipe
    from wfcommons import WorkflowGenerator

    generator = WorkflowGenerator(SeismologyRecipe.from_num_tasks(250))
    workflow = generator.build_workflow()
    workflow.write_json(pathlib.Path("seismology-workflow.json"))

The resulting JSON file follows :ref:`json-format-label` and can be consumed
by any tool or simulator that supports it.

Your first workflow benchmark
-----------------------------

The same recipes can also produce **runnable benchmarks** with controlled CPU
and I/O behavior, translated to the workflow system of your choice:

.. code-block:: python

    import pathlib
    from wfcommons import BlastRecipe
    from wfcommons.wfbench import WorkflowBenchmark, BashTranslator

    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=50)
    benchmark.create_benchmark(pathlib.Path("/tmp/"), cpu_work=100, data=10, percent_cpu=0.6)

    translator = BashTranslator(benchmark.workflow)
    translator.translate(output_folder=pathlib.Path("./bash-wf/"))

See :ref:`generating-workflow-benchmarks-label` for the full list of supported
workflow systems.

Installing from source (latest)
-------------------------------

If you want the latest development version (potentially unstable), clone the
repository and install locally:

.. code-block:: bash

    $ git clone https://github.com/wfcommons/wfcommons
    $ cd wfcommons
    $ python3 -m pip install .

Optional requirements
---------------------

Visualization support (drawing workflow task graphs and reading/writing DOT
files) is available as an extra:

.. code-block:: bash

    $ python3 -m pip install wfcommons[viz]

Graphviz
^^^^^^^^

WfCommons uses `pygraphviz <https://pygraphviz.github.io/documentation/latest/install.html>`_
for generating visualizations of the workflow task graph. Building pygraphviz
requires the `graphviz <https://www.graphviz.org/>`_ package (version 2.16 or
later) and its development headers. You can install graphviz easily on Linux
with your favorite package manager, for example for Debian-based
distributions:

.. code-block:: bash

    $ sudo apt-get install graphviz libgraphviz-dev

and for RedHat-based distributions:

.. code-block:: bash

    $ sudo yum install python-devel graphviz-devel

On macOS you can use the :code:`brew` package manager:

.. code-block:: bash

    $ brew install graphviz

Then install the visualization extra (or pygraphviz directly):

.. code-block:: bash

    $ python3 -m pip install wfcommons[viz]

Next steps
----------

- :doc:`introduction` — understand the WfCommons components and WfFormat.
- :doc:`analyzing_instances` — parse execution logs and analyze real instances.
- :doc:`generating_workflows_recipe` — build a recipe for your own application.
- :doc:`generating_workflow_benchmarks` — generate benchmarks for real systems.

Installation
============

WfCommons is available on `PyPI <https://pypi.org/project/wfcommons>`_.
WfCommons requires Python3.11+ and has been tested on Linux and macOS.

Recommended setup
----------------------

We recommend installing into a virtual environment to keep dependencies isolated:

.. code-block:: bash

    $ python3 -m venv .venv
    $ source .venv/bin/activate
    $ python3 -m pip install wfcommons

Verify installation
-------------------

You can confirm the CLI entry point is available:

.. code-block:: bash

    $ wfchef --help

Or check the Python import:

.. code-block:: bash

    $ python3 -c "import wfcommons; print(wfcommons.__version__)"

Installing from source (latest)
-------------------------------

If you want the latest development version (potentially unstable), clone the
repository and install locally:

.. code-block:: bash

    $ git clone https://github.com/wfcommons/wfcommons
    $ cd wfcommons
    $ python3 -m pip install .

Optional Requirements
---------------------

Graphviz
^^^^^^^^

WfCommons uses `pygraphviz <https://pygraphviz.github.io/documentation/latest/install.html>`_
for generating visualizations for the workflow task graph . If you want to enable this 
feature, you will have to install the `graphviz <https://www.graphviz.org/>`_ package (version 2.16
or later). You can install graphviz easily on Linux with your favorite package manager,
for example for Debian-based distributions:

.. code-block:: bash

    $ sudo apt-get install graphviz libgraphviz-dev

and for RedHat-based distributions:

.. code-block:: bash

    $ sudo yum install python-devel graphviz-devel

On macOS you can use the :code:`brew` package manager:

.. code-block:: bash

    $ brew install graphviz

Then you can install pygraphviz by running:

.. code-block:: bash

    $ python3 -m pip install pygraphviz

Installation
============

WfCommons is available on `PyPI <https://pypi.org/project/wfcommons>`_.
WfCommons requires Python3.6+ and has been tested on Linux and macOS.

Installation using pip
----------------------

While :code:`pip` can be used to install WfCommons, we suggest the following
approach for reliable installation when many Python environments are available:

.. code-block:: bash

    $ python3 -m pip install wfcommons

Retrieving the latest unstable version
--------------------------------------

If you want to use the latest WfCommons unstable version, that will contain
brand new features (but also contain bugs as the stabilization work is still
underway), you may consider retrieving the latest unstable version.

Cloning from `WfCommons <https://github.com/wfcommons/wfcommons>`_'s GitHub
repository: ::

    $ git clone https://github.com/wfcommons/wfcommons
    $ cd wfcommons
    $ pip install .

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

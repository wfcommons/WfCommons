Installation
============

WorkflowHub is available on `PyPI <https://pypi.org/project/workflowhub>`_.
WorkflowHub requires Python3.3+ and has been tested on Linux and macOS.

Installation using pip
----------------------

While :code:`pip` can be used to install WorkflowHub, we suggest the following
approach for reliable installation when many Python environments are available:

.. code-block:: bash

    $ python3 -m pip install workflowhub

Retrieving the latest unstable version
--------------------------------------

If you want to use the latest WorkflowHub unstable version, that will contain
brand new features (but also contain bugs as the stabilization work is still
underway), you may consider retrieving the latest unstable version.

Cloning from `WorkflowHub <https://github.com/workflowhub/workflowhub>`_'s GitHub
repository: ::

    $ git clone https://github.com/workflowhub/workflowhub
    $ cd workflowhub
    $ pip install .

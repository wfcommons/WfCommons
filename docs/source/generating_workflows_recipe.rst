.. _generating-workflows-recipe-label:

WfChef: Workflow Recipes
========================

**WfChef** is the WfCommons component that automates the construction of
synthetic workflow generators — no manual modeling, and no code to write. The
input to this component is a set of real workflow instances described in the
:ref:`json-format-label` (e.g., instances available in
:ref:`wfinstances-label`). WfChef automatically analyzes the set of real
workflow instances for two purposes:

1. It discovers workflow subgraphs that represent fundamental task dependency
   patterns of the application.
2. It derives statistical models of the workflow tasks' performance
   characteristics (see :ref:`instances-label`).

WfChef then outputs a **recipe** that is used by **WfGen** (see
:ref:`generating-workflows-label`) to generate realistic synthetic workflow
instances with any arbitrary number of tasks.

.. _workflow-recipe-label:

Workflow Recipes
----------------

A **workflow recipe** is a data structure that encodes the discovered pattern
occurrences as well as the statistical models of workflow task
characteristics. More precisely, a recipe embodies results from statistical
analysis and distribution fitting performed for each workflow task type so as
to characterize task runtime and input/output data sizes. The recipe also
incorporates information regarding the graph structure of the workflows (task
dependencies and frequency of occurrences), automatically derived from the
analysis of the workflow instances.

This Python package already provides several *workflow recipes* (see
:ref:`recipes-list`) for generating realistic synthetic workflow instances —
if one of those covers your application, you can skip straight to
:ref:`generating-workflows-label`.

.. _workflow-recipe-generator-label:

Creating a Recipe from Real Instances
-------------------------------------

To create a recipe, WfChef analyzes the real workflow graphs in order to
identify subgraphs that represent fundamental task dependency patterns. Based
on the identified subgraphs and on measured task type frequencies in the real
workflows, WfChef outputs a generator that can produce realistic synthetic
workflow instances with an arbitrary number of tasks (see
:ref:`generating-workflows-label`).

Recipes are created with the :code:`wfchef` command-line tool, which is
installed with the package. The example below creates a recipe for the
Epigenomics application from a folder of real instances in
:ref:`json-format-label`:

.. code-block:: bash

    $ wfchef create /path/to/real/instances -o ./epigenomics -v --name Epigenomics

The following flags can be used with this command:

- :code:`-o` or :code:`--out` (required): name of the directory to be created
  that will contain the recipe.
- :code:`-n` or :code:`--name` (required): name of the recipe. Typically, the
  format used is *ApplicationName*.
- :code:`-v` or :code:`--verbose`: if set, activates status messages.
- :code:`--no-install`: if set, does not install the recipe automatically.
- :code:`-r` or :code:`--runs`: number of runs used to compute the mean RMSE
  of the fitted models (default: 1).
- :code:`-c` or :code:`--cutoff`: only consider real instances with at most
  this number of tasks when creating the recipe (default: 4000).

*Example*: with :code:`--cutoff 4000`, all real instances considered for the
creation of the recipe will have 4000 or fewer tasks. This is a useful flag
when there is confidence that all patterns present in the application can
already be found in the smaller instances (recipe creation is then faster).

Installing and Managing Recipes
-------------------------------

Workflow recipes are automatically installed (as small Python packages) and
can then be used throughout the system. WfCommons creates a Python package in
the directory specified by the :code:`--out` flag, in which the
:file:`setup.py` and :file:`recipe.py` files are stored. If the
:code:`--no-install` flag was set when creating the recipe, you will need to
manually install the package before using it:

.. code-block:: bash

    # installing the recipe package
    $ pip install /path/to/the/package

To list the recipes installed in a system (and how to import them):

.. code-block:: bash

    $ wfchef ls

To uninstall a recipe:

.. code-block:: bash

    $ wfchef uninstall -n recipe_name

Using a Recipe
--------------

Once installed, a recipe is imported like any of the bundled ones and handed
to the :class:`~wfcommons.wfgen.generator.WorkflowGenerator` class: ::

    from wfcommons.wfchef.recipes import EpigenomicsRecipe
    from wfcommons import WorkflowGenerator

    generator = WorkflowGenerator(EpigenomicsRecipe.from_num_tasks(1000))
    workflow = generator.build_workflow()

See :ref:`generating-workflows-label` for the full generation guide.

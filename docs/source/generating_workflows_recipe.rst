.. _generating-workflows-recipe-label:

Generating Workflows Recipes
============================

**WfChef** is the WfCommons component that automates the construction of
synthetic workflow generators for any given workflow application. The input
to this component is a set of real workflow instances described in the
:ref:`json-format-label` (e.g., instances available in :ref:`wfinstances-label`).
WfChef automatically analyzes a set of real workflow instances for
two purposes. First, it discovers workflow subgraphs that represent
fundamental task dependency patterns. Second, it derives
statistical models of the workflow tasks' performance characteristics
(see :ref:`instances-label`).
WfChef then outputs a **recipe** that will be used by **WfGen**
(see :ref:`generating-workflows-label`) to generate realistic synthetic  
workflow instances with any arbitrary number of tasks.

.. _workflow-recipe-label:

Workflow Recipes
----------------

A **workflow recipe** is a data structure that encodes the discovered
pattern occurrences as well as the statistical models of workflow task
characteristics. More precisely, a recipe embodies results from statistical
analysis and distribution fitting performed for each workflow task type
so as to characterize task runtime and input/output data sizes. The
recipes also incorporates information regarding the graph structure of
the workflows (tasks dependencies and frequency of occurrences), which are
automatically derived from the analysis of the workflow instances.

This Python package provides several *workflow recipes* (see :ref:`recipes-list`)
for generating realistic synthetic workflow instances.

.. _workflow-recipe-generator-label:

Workflow Recipe Generator
--------------------------

To create a recipe, WfChef analyzes the real workflow graphs in order to
identify subgraphs that represent fundamental task dependency patterns.
Based on the identified subgraphs and on measured task type frequencies in the real
workflows, WfChef outputs a generator that can generate realistic synthetic
workflow instances with an arbitrary numbers of tasks (see :ref:`generating-workflows-label`).

The code snippet below shows an example of how to create a recipe for the
Epigenomics application: ::

    $ wfchef create /path/to/real/instances -o ./epigenomics -v --name Epigenomics

The following flags can be used with this command:

- :code:`-o` or :code:`--out` is a required flag that stands for the name of the directory to be created that is going to contain the recipe.
- :code:`-n` or :code:`--name` is a required flag that stands for the name of the recipe. Tipically, the format used is *ApplicationName*.
- :code:`-v` or :code:`--verbose` if set, activates status messages.
- :code:`--no-install` if set, does not install the recipe automatically.
- :code:`-c` or :code:`--cutoff` takes a number of tasks in the samples the user wants to consider to create the recipe.

*Example*: :code:`--cutoff 4000`, it means that all real world instances
that will be consider for the creation of the recipe will have 4000 or
less tasks. This is a useful flag to use when there is trust that all
possible patterns present in this application can be already found in the
smaller instances.

Workflow recipes are automatically installed and can be used throughout the
system. WfCommons creates a Python package in the directory specified by the
flag :code:`--out` in which the :file:`setup.py` and :file:`recipe.py` files
are stored. If the flag :code:`--no-install` is set when creating a package
for a specific application, the user will need to manually install the package
before using it. The code bellow is an example of how to install/uninstall a
package for an application in WfCommons: ::

    # installing the package
    $ pip install /path/to/the/package

    # uninstalling a package
    $ pip uninstall wfcommons.wfchef.recipes.appication_name_workflow

The snippet below shows an example of how to import the recipes: ::

    from wfcommons.wfchef.recipes import EpigenomicsRecipe

To check which recipes are installed in a system and how to import them use: ::
    
    $ wfchef ls

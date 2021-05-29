.. _generating-workflows-recipe-label:

Generating Workflows Recipe
============================

**WfChef** is the WfCommons component that automates the construction of
synthetic workflow generators for any given workflow application. The input
to this component is a set of real workflow instances described in the
*WfFormat* (e.g., instances available in **WfInstances**).
WfChef automatically analyzes the real workflow instances for
two purposes. First, it discovers workflow subgraphs that represent
fundamental task dependency patterns. Second, it derives
statistical models of the workflow tasks' performance characteristics (more details :ref:`.. _traces-label:`).
WfChef then outputs a **recipe** that will be used by **WfGen** 
(see :ref:`generating-workflows-label`) to generate realistic synthetic  
workflow instances with any arbitrary number of tasks.

.. _workflow-recipe-label:

Workflow Recipes
----------------

A **workfflow recipe** is a data structure that encodes the discovered pattern occurrences 
as well as the statistical models of workflow task characteristics.
The WfCommons package provides a number of *workflow recipes* for generating realistic 
synthetic workflow instances. 

All workflow recipes provide a common method, :code:`from_num_tasks`, that defines the upper 
bound for the total number of tasks in the synthetic workflow.


.. _workflow-recipe-generator-label:

Workflow Recipe Generator
--------------------------

To create a recipe the method :meth:`find_patterns` analyzes the real workflow graphs 
in order to identify subgraphs that represent fundamental task dependency patterns. 
Based on the identified subgraphs and on measured task type frequencies in the real
workflows, WfChef outputs a generator that can generate realistic synthetic
workflow instances with an arbitrary numbers of tasks. The code snippet below shows 
how to create a recipe for the Epigenomics application: ::

    $ wfchef create path/to/real/instances -o ./epigenomics -v --name Epigenomics

The flags that can be used with this command are:

:code:`-o` or :code:`--out` is a required flag that stands for the name of the directory to be created that is going to 
contain the recipe.

:code:`-n` or :code:`--name` is a required flag that stands for the name of the recipe. Tipically, the format used is 
*ApplicationName*. 

:code:`-v` or :code:`--verbose` if set, activates status messages.

:code:`--no-install` if set, does not install the recipe automatically.

:code:`-c` or :code:`--cutoff` takes a number of tasks in the samples the user wants to consider to create the recipe. 
Example: :code:`--cutoff 4000`, it means that all real world instances that will be consider for the creation of the 
recipe have 4000 or less tasks. This is a useful flag to use when there is trust that all possible patterns present
in this application can be already found in the smaller instances. 

The recipes are automatically installed and can be used throughout the system. WfCommons creates a Python package in the directory 
specified by the flag :code:`--out` in which the :file:`setup.py` and :file:`recipe.py` files are stored. If the flag :code:`--no-install` is set 
when creating a package for a specific application, the user will need to install the package before using it. The code 
bellow is an example of how to install/uninstall a package for an application in WfCommons: ::

    # installing the package
    $ pip install path/to/the/package

    # uninstalling a package
    $ pip uninstall wfcommons.wfchef.recipes.appication_name_workflow



The snippet below show an example of how to import the recipes: ::

    # creating an Epigenomics workflow recipe
    from wfcommons.wfchef.recipes import EpigenomicsRecipe


To check which recipes are installed in a system and how to import them use: ::
    
    $ wfchef ls


The current list of available workflow recipes include:

- :class:`~wfcommons.wfchef.recipes.blast_recipe.BlastRecipe`: :code:`from wfcommons.wfchef.recipes import BlastWorkflowRecipe`
- :class:`~wfcommons.wfchef.recipes.bwa_recipe.BwaRecipe`: :code:`from wfcommons.wfchef.recipes import BwaRecipe`
- :class:`~wfcommons.wfchef.recipes.cycles_recipe.CyclesRecipe`: :code:`from wfcommons.wfchef.recipes import CyclesRecipe`
- :class:`~wfcommons.wfchef.recipes.epigenomics_recipe.EpigenomicsRecipe`: :code:`from wfcommons.wfchef.recipes import EpigenomicsRecipe`
- :class:`~wfcommons.wfchef.recipes.genome_recipe.GenomeRecipe`: :code:`from wfcommons.wfchef.recipes import GenomeRecipe`
- :class:`~wfcommons.wfchef.recipes.montage_recipe.MontageRecipe`: :code:`from wfcommons.wfchef.recipes import MontageRecipe`
- :class:`~wfcommons.wfchef.recipes.seismology_recipe.SeismologyRecipe`: :code:`from wfcommons.wfchef.recipes import SeismologyRecipe`
- :class:`~wfcommons.wfchef.recipes.soykb_recipe.SoykbRecipe`: :code:`from wfcommons.wfchef.recipes import SoykbRecipe`
- :class:`~wfcommons.wfchef.recipes.srasearch_recipe.SrasearchRecipe`: :code:`from wfcommons.wfchef.recipes import SrasearchRecipe`



Examples
--------

The following example generates 10 *Epigenomics* synthetic workflow instances
based on the number of tasks entered by the user (1000), builds the synthetic workflow instances, and writes the
synthetic instances to JSON files. ::

    from wfcommons.wfchef.recipes import EpigenomicsRecipe
    from wfcommons.generator import WorkflowGenerator

    generator = WorkflowGenerator(EpigenomicsRecipe.from_num_tasks(1000)) 
    for i, workflow in enumerate(generator.build_workflows(10)):
        workflow.write_json(f'epigenomics-workflow-{i}.json')

The example below generates a *Cycles* (agroecosystem) synthetic workflow instance based on the number 
of tasks entered by the user (250), builds the synthetic workflow instance, and writes the synthetic 
instance to a JSON file. ::
    
    from wfcommons.wfchef.recipes import CyclesRecipe
    from wfcommons.generator import WorkflowGenerator

    generator = WorkflowGenerator(CyclesRecipe.from_num_tasks(250)) 
    workflow = generator.build_workflow()
    workflow.write_json(f'cycles-workflow.json')

..
    maybe we should pout examples only on generator, because we need it

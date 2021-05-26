.. _generating-workflows-recipe-label:

Generating Workflows Recipe
============================

The third axis of the WfCommons project is called **WfChef**, a component 
responsible for automating the process of constructing a synthetic workflow
generator for any given workflow application. WfChef takes as input a set
of real workflow instances from an application, and outputs the 
code of a synthetic workflow generator for that application.

.. _workflow-recipe-generator-label:

Workflow Recipe Generator
--------------------------

To create a recipe the method :meth:`find_patterns` analyzes the real workflow graphs 
in order to identify subgraphs that represent fundamental task dependency patterns. 
Based on the identified subgraphs and on measured task type frequencies in the real
workflows, WfChef outputs a generator that can generate realistic synthetic
workflow instances with an arbitrary numbers of tasks. The code snippet below shows 
how to create a recipe for the Epigenomics application: ::

    $ wfchef create path/to/real/instances -o ./epigenomics -v --name EpigenomicsWorkflow

The flags that can be used with this command are:

:code:`-o` or :code:`--out` is a required flag that stands for the name of the directory to be created that is going to 
contain the recipe.

:code:`-n` or :code:`--name` is a required flag that stands for the name of the recipe. Tipically, the format used is 
*ApplicationNameWorkflow*. 

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
    $ pip uninstall wfchef.recipe.appication_name_workflow

The *wfchef.recipe.appication_name_workflow* can be found in :file:`setup.py` under the tag :code:`name`. The example bellow shows
the :meth:`setup` method for makeflow's BWA application: ::

    # BWA's setup method in setup.py
    setup(
    name='wfchef.recipe.b_w_a_workflow', # information necessary to install and uninstall a package
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'networkx',
        'wfcommons'
    ],
    )


The snippet below show an example of how to import the recipes: ::

    # creating an Epigenomics workflow recipe
    from epigenomics_workflow_recipe import EpigenomicsWorkflowRecipe


To check which recipes are installed in a system and how to import them use: ::
    
    $ wfchef ls


The current list of available workflow recipes include:

- :class:`~wfcommons.generator.workflow.blast_recipe.BLASTRecipe`: :code:`from blast_workflow_recipe import BlastWorkflowRecipe`
- :class:`~wfcommons.generator.workflow.bwa_recipe.BWARecipe`: :code:`from bwa_workflow_recipe import BWAWorkflowRecipe`
- :class:`~wfcommons.generator.workflow.cycles_recipe.CyclesRecipe`: :code:`from cycles_workflow_recipe import CyclesWorkflowRecipe`
- :class:`~wfcommons.generator.workflow.epigenomics_recipe.EpigenomicsRecipe`: :code:`from epigenomics_workflow_recipe import EpigenomicsWorkflowRecipe`
- :class:`~wfcommons.generator.workflow.genome_recipe.GenomeRecipe`: :code:`from genome_workflow_recipe import GenomeWorkflowRecipe`
- :class:`~wfcommons.generator.workflow.montage_recipe.MontageRecipe`: :code:`from montage_workflow_recipe import MontageWorkflowRecipe`
- :class:`~wfcommons.generator.workflow.seismology_recipe.SeismologyRecipe`: :code:`from seismology_workflow_recipe import SeismologyWorkflowRecipe`
- :class:`~wfcommons.generator.workflow.soykb_recipe.SoyKBRecipe`: :code:`from soykb_workflow_recipe import SoyKbWorkflowRecipe`
- :class:`~wfcommons.generator.workflow.srasearch_recipe.SRASearchRecipe`: :code:`from srasearch_workflow_recipe import SRASearchWorkflowRecipe`



Examples
--------

The following example generates 10 *Epigenomics* synthetic workflow instances
based on the number of tasks entered by the user (1000), builds the synthetic workflow instances, and writes the
synthetic instances to JSON files. ::

    from epigenomics_workflow_recipe import EpigenomicsWorkflowRecipe
    from wfcommons.generator import WorkflowGenerator

    generator = WorkflowGenerator(EpigenomicsWorkflowRecipeRecipe.from_num_tasks(1000)) 
    for i, workflow in enumerate(generator.build_workflows(10)):
        workflow.write_json(f'epigenomics-workflow-{i}.json')

The example below generates a *Cycles* (agroecosystem) synthetic workflow instance based on the number 
of tasks entered by the user (250), builds the synthetic workflow instance, and writes the synthetic 
instance to a JSON file. ::
    
    from cycles_workflow_recipe import CyclesWorkflowRecipe
    from wfcommons.generator import WorkflowGenerator

    generator = WorkflowGenerator(CyclesWorkflowRecipeRecipe.from_num_tasks(250)) 
    workflow = generator.build_workflows(1)
    workflow.write_json(f'cycles-workflow.json')

..
    maybe we should pout examples only on generator, because we need it

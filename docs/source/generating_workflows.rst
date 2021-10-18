.. _generating-workflows-label:

Generating Workflows
====================

WfGen is a component of WfCommons project that targets the generation of realistic
synthetic workflow instances with a variety of characteristics. The
:class:`~wfcommons.wfgen.generator.WorkflowGenerator` class uses recipes
of workflows (as described in :ref:`workflow-recipe-generator-label`) 
for creating the realistic synthetic instances. The resulting workflows are represented in the 
:ref:`json-format-label`, which is already supported by simulation frameworks such as
`WRENCH <https://wrench-project.org>`_.

.. _recipes-list:

WfCommons Workflows Recipes
---------------------------

This Python package provides several *workflow recipes* for generating realistic
synthetic workflow instances. The current list of available workflow recipes include:

- :class:`~wfcommons.wfchef.recipes.blast.recipe.BlastRecipe`: :code:`from wfcommons.wfchef.recipes import BlastRecipe`
- :class:`~wfcommons.wfchef.recipes.bwa.recipe.BwaRecipe`: :code:`from wfcommons.wfchef.recipes import BwaRecipe`
- :class:`~wfcommons.wfchef.recipes.cycles.recipe.CyclesRecipe`: :code:`from wfcommons.wfchef.recipes import CyclesRecipe`
- :class:`~wfcommons.wfchef.recipes.epigenomics.recipe.EpigenomicsRecipe`: :code:`from wfcommons.wfchef.recipes import EpigenomicsRecipe`
- :class:`~wfcommons.wfchef.recipes.genome.recipe.GenomeRecipe`: :code:`from wfcommons.wfchef.recipes import GenomeRecipe`
- :class:`~wfcommons.wfchef.recipes.montage.recipe.MontageRecipe`: :code:`from wfcommons.wfchef.recipes import MontageRecipe`
- :class:`~wfcommons.wfchef.recipes.seismology.recipe.SeismologyRecipe`: :code:`from wfcommons.wfchef.recipes import SeismologyRecipe`
- :class:`~wfcommons.wfchef.recipes.soykb.recipe.SoykbRecipe`: :code:`from wfcommons.wfchef.recipes import SoykbRecipe`
- :class:`~wfcommons.wfchef.recipes.srasearch.recipe.SrasearchRecipe`: :code:`from wfcommons.wfchef.recipes import SrasearchRecipe`

The Workflow Instances Generator
--------------------------------

Synthetic workflow instances are generated using the
:class:`~wfcommons.wfgen.generator.WorkflowGenerator` class. This class takes
as input a :class:`~wfcommons.wfgen.abstract_recipe.WorkflowRecipe`
object (see in :ref:`workflow-recipe-generator-label`), and provides two methods
for generating synthetic workflow instances:

- :meth:`~wfcommons.wfgen.generator.WorkflowGenerator.build_workflow`: generates a single synthetic workflow
  instance based on the workflow recipe used to instantiate the generator.
- :meth:`~wfcommons.wfgen.generator.WorkflowGenerator.build_workflows`: generates a number of synthetic workflow
  instances based on the workflow recipe used to instantiate the generator.

The build methods use the workflow recipe for generating realistic synthetic
workflow instances, in which the workflow structure follows workflow composition
rules defined in the workflow recipe, and tasks runtime, and input and output
data sizes are generated according to distributions obtained from actual workflow
execution instances (see :ref:`instances-label`).

Each generated instance is represented as a :class:`~wfcommons.common.workflow.Workflow`
object (which in itself is an extension of the
`NetworkX DiGraph <https://networkx.github.io/documentation/stable/reference/classes/digraph.html>`_
class). The :class:`~wfcommons.common.workflow.Workflow` class provides two
methods for writing the generated workflow instance into files:

- :meth:`~wfcommons.common.workflow.Workflow.write_dot`: write a DOT file of a workflow instance.
- :meth:`~wfcommons.common.workflow.Workflow.write_json`: write a JSON file of a workflow instance.

All workflow recipes provide a common method, :code:`from_num_tasks`, that defines the lower
bound for the total number of tasks in the generated synthetic workflow.

Increasing/Reducing Runtime and File Sizes
******************************************

Workflow recipes also allow the generation of synthetic workflows with increased/reduced
runtimes and/or files sizes determined by a factor provided by the user:

- :code:`runtime_factor`: The factor of which tasks runtime will be increased/decreased.
- :code:`input_file_size_factor`: The factor of which tasks input files size will be increased/decreased.
- :code:`output_file_size_factor`: The factor of which tasks output files size will be increased/decreased.

The following example shows how to create a Seismology workflow recipe in which task
runtime is increased by 10%, input files by 50%, and output files reduced by 20%: ::

    from wfcommons.wfchef.recipes import SeismologyRecipe

    # creating a Seismology workflow recipe with increased/decreased runtime and file sizes
    recipe = SeismologyRecipe.from_num_tasks(num_tasks=100, runtime_factor=1.1, input_file_size_factor=1.5, output_file_size_factor=0.8)

Examples
--------

The following example generates a *Seismology* synthetic workflow instance
os 300 tasks, builds a synthetic workflow instance, and writes the
synthetic instance to a JSON file. ::

    import pathlib
    from wfcommons.wfchef.recipes import SeismologyRecipe
    from wfcommons import WorkflowGenerator

    generator = WorkflowGenerator(SeismologyRecipe.from_num_tasks(250))
    workflow = generator.build_workflow()
    workflow.write_json(pathlib.Path('seismology-workflow.json'))


The example below generates a number of 10 *Blast* synthetic
workflow instances for every size defined in the array :code:`num_tasks`: ::

    import pathlib
    from wfcommons.wfchef.recipes import BlastRecipe
    from wfcommons import WorkflowGenerator

    num_tasks = [100, 250, 370, 800]

    for task in num_tasks:
        generator = WorkflowGenerator(BlastRecipe.from_num_tasks(task))
        workflows = generator.build_workflows(10)

        for i, workflow in enumerate(workflows):
            workflow.write_json(pathlib.Path(f'blast-workflow-{task}-{i}.json'))


The following example generates 10 *Epigenomics* synthetic workflow instances
based on the number of tasks entered by the user (1000), builds the synthetic
workflow instances, and writes the synthetic instances to JSON files. ::

    import pathlib
    from wfcommons.wfchef.recipes import EpigenomicsRecipe
    from wfcommons import WorkflowGenerator

    generator = WorkflowGenerator(EpigenomicsRecipe.from_num_tasks(1000))
    for i, workflow in enumerate(generator.build_workflows(10)):
        workflow.write_json(pathlib.Path(f'epigenomics-workflow-{i}.json'))

The example below generates a *Cycles* (agroecosystem) synthetic workflow instance
based on the number of tasks entered by the user (250), builds the synthetic workflow
instance, and writes the synthetic instance to a JSON file. ::

    import pathlib
    from wfcommons.wfchef.recipes import CyclesRecipe
    from wfcommons import WorkflowGenerator

    generator = WorkflowGenerator(CyclesRecipe.from_num_tasks(250))
    workflow = generator.build_workflow()
    workflow.write_json(pathlib.Path('cycles-workflow.json'))

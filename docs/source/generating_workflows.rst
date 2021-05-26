.. _generating-workflows-label:

Generating Workflows
====================

WfGen is a component of WfCommons project that targets the generation of realistic
synthetic workflow instances with a variety of characteristics. The
:class:`~wfcommons.generator.generator.WorkflowGenerator` class uses recipes
of workflows, many different synthetic workflows based on distributions of workflow 
task runtime,  and input and output file sizes (as described in :ref:`workflow-recipe-generator-label`) 
for creating the realistic synthetic instances. The resulting workflows are represented in the 
WfCommons JSON format (WfFormat), which is already supported by simulation frameworks such as
`WRENCH <https://wrench-project.org>`_.

.. _workflow-recipe-label:

Workflow Recipes
----------------

The WfCommons package provides a number of *workflow recipes* for generating
realistic synthetic workflow instances. Each recipe may provide their own methods
for instantiating a :class:`~wfcommons.generator.workflow.abstract_recipe.WorkflowRecipe`
object depending on the properties that define the structure of the actual
workflow. For instance, the code snippet below shows how to instantiate a recipe
of the Epigenomics and 1000Genome workflows: ::

    from wfcommons.generator import EpigenomicsRecipe, GenomeRecipe

    # creating an Epigenomics workflow recipe
    epigenomics_recipe = EpigenomicsRecipe.from_sequences(num_sequence_files=2, num_lines=100, bin_size=10)

    # creating a 1000Genome workflow recipe
    genome_recipe = GenomeRecipe.from_num_chromosomes(num_chromosomes=3, num_sequences=10000, num_populations=1)


All workflow recipes also provide a common method (:code:`from_num_tasks`) for
instantiating a :class:`~wfcommons.generator.workflow.abstract_recipe.WorkflowRecipe`
object as follows: ::

    from wfcommons.generator import EpigenomicsRecipe, GenomeRecipe

    # creating an Epigenomics workflow recipe
    epigenomics_recipe = EpigenomicsRecipe.from_num_tasks(num_tasks=9)

    # creating a 1000Genome workflow recipe
    genome_recipe = GenomeRecipe.from_num_tasks(num_tasks=5)

Note that :code:`num_tasks` defines the upper bound for the total number of tasks in the
workflow, and that each workflow recipe may define different lower bound values so
that the workflow structure is guaranteed. Please, refer to the :ref:`documentation of
each workflow recipe <wfcommons-generator-label>` for the lower bound values.

The current list of available workflow recipes include:

- :class:`~wfcommons.generator.workflow.blast_recipe.BLASTRecipe`: :code:`from wfcommons.generator import BLASTRecipe`
- :class:`~wfcommons.generator.workflow.bwa_recipe.BWARecipe`: :code:`from wfcommons.generator import BWARecipe`
- :class:`~wfcommons.generator.workflow.cycles_recipe.CyclesRecipe`: :code:`from wfcommons.generator import CyclesRecipe`
- :class:`~wfcommons.generator.workflow.epigenomics_recipe.EpigenomicsRecipe`: :code:`from wfcommons.generator import EpigenomicsRecipe`
- :class:`~wfcommons.generator.workflow.genome_recipe.GenomeRecipe`: :code:`from wfcommons.generator import GenomeRecipe`
- :class:`~wfcommons.generator.workflow.montage_recipe.MontageRecipe`: :code:`from wfcommons.generator import MontageRecipe`
- :class:`~wfcommons.generator.workflow.seismology_recipe.SeismologyRecipe`: :code:`from wfcommons.generator import SeismologyRecipe`
- :class:`~wfcommons.generator.workflow.soykb_recipe.SoyKBRecipe`: :code:`from wfcommons.generator import SoyKBRecipe`
- :class:`~wfcommons.generator.workflow.srasearch_recipe.SRASearchRecipe`: :code:`from wfcommons.generator import SRASearchRecipe`

Increasing/Reducing Runtime and File Sizes
******************************************

Workflow recipes also allow the generation of synthetic workflows with increased/reduced
runtimes and/or files sizes determined by a factor provided by the user:

- :code:`runtime_factor`: The factor of which tasks runtime will be increased/decreased.
- :code:`input_file_size_factor`: The factor of which tasks input files size will be increased/decreased.
- :code:`output_file_size_factor`: The factor of which tasks output files size will be increased/decreased.

The following example shows how to create a Seismology workflow recipe in which task
runtime is increased by 10%, input files by 50%, and output files reduced by 20%: ::

    from wfcommons.generator import SeismologyRecipe

    # creating a Seismology workflow recipe with increased/decreased runtime and file sizes
    recipe = SeismologyRecipe.from_num_tasks(num_tasks=100, runtime_factor=1.1, input_file_size_factor=1.5, output_file_size_factor=0.8)

The Workflow Generator
----------------------

Synthetic workflow instances are generated using the
:class:`~wfcommons.generator.generator.WorkflowGenerator` class. This
class takes as input a :class:`~wfcommons.generator.workflow.abstract_recipe.WorkflowRecipe`
object (see above), and provides two methods for generating synthetic
workflow instances:

- :meth:`~wfcommons.generator.generator.WorkflowGenerator.build_workflow`: generates a single synthetic workflow
  instance based on the workflow recipe used to instantiate the generator.
- :meth:`~wfcommons.generator.generator.WorkflowGenerator.build_workflows`: generates a number of synthetic workflow
  instances based on the workflow recipe used to instantiate the generator.

The build methods use the workflow recipe for generating realistic synthetic
workflow instances, in which the workflow structure follows workflow composition
rules defined in the workflow recipe, and tasks runtime, and input and output
data sizes are generated according to distributions obtained from actual workflow
execution instances (see :ref:`traces-label`).

Each generated instance is a represented as a :class:`~wfcommons.common.workflow.Workflow`
object (which in itself is an extension of the
`NetworkX DiGraph <https://networkx.github.io/documentation/stable/reference/classes/digraph.html>`_
class). The :class:`~wfcommons.common.workflow.Workflow` class provides two
methods for writing the generated workflow instance into files:

- :meth:`~wfcommons.common.workflow.Workflow.write_dot`: write a DOT file of a workflow instance.
- :meth:`~wfcommons.common.workflow.Workflow.write_json`: write a JSON file of a workflow instance.

Examples
--------

The following example generates a *Seismology* synthetic workflow instance
based on the number of pair of signals to estimate earthquake STFs
(:code:`num_pairs`), builds a synthetic workflow instance, and writes the
synthetic instance to a JSON file. ::

    from wfcommons import WorkflowGenerator
    from wfcommons.generator import SeismologyRecipe

    # creating a Seismology workflow recipe based on the number
    # of pair of signals to estimate earthquake STFs
    recipe = SeismologyRecipe.from_num_pairs(num_pairs=10)

    # creating an instance of the workflow generator with the
    # Seismology workflow recipe
    generator = WorkflowGenerator(recipe)

    # generating a synthetic workflow instance of the Seismology workflow
    workflow = generator.build_workflow()

    # writing the synthetic workflow instance into a JSON file
    workflow.write_json('seismology-workflow.json')


The example below generates a number of *Cycles* (agroecosystem) synthetic
workflow instances based on the upper bound number of tasks allowed per workflow. ::

    from wfcommons import WorkflowGenerator
    from wfcommons.generator import CyclesRecipe

    # creating a Cycles workflow recipe based on the number of tasks per workflow
    recipe = CyclesRecipe.from_num_tasks(num_tasks=1000)

    # creating an instance of the workflow generator with the
    # Cycles workflow recipe
    generator = WorkflowGenerator(recipe)

    # generating 10 synthetic workflow instances of the Cycles workflow
    workflows_list = generator.build_workflows(num_workflows=10)

    # writing each synthetic workflow instance into a JSON file
    count = 1
    for workflow in workflows_list:
        workflow.write_json('cycles-workflow-{:02}.json'.format(count))
        count += 1

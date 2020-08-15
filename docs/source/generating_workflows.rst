Generating Workflows
====================

The second axis of the WorkflowHub project targets the generation of realistic
synthetic workflow traces with a variety of characteristics. The
:class:`~workflowhub.generator.generator.WorkflowGenerator` class uses recipes
of workflows (as described in :ref:`traces-label`) for creating many different
synthetic workflows based on distributions of workflow job runtime, and input
and output file sizes.
The resulting workflows are represented in the WorkflowHub JSON format, which
is already supported by simulation frameworks such as
`WRENCH <https://wrench-project.org>`_.

Workflow Recipes
----------------

The WorkflowHub package provides a number of *workflow recipes* for generating
realistic synthetic workflow traces. Each recipe may provide their own methods
for instantiating a :class:`~workflowhub.generator.workflow.abstract_recipe.WorkflowRecipe`
object depending on the properties that define the structure of the actual
workflow. For instance, the code snippet below shows how to instantiate a recipe
of the Epigenomics and 1000Genome workflows: ::

    from workflowhub.generator import EpigenomicsRecipe, GenomeRecipe

    # creating an Epigenomics workflow recipe
    epigenomics_recipe = EpigenomicsRecipe.from_sequences(num_sequence_files=2, num_lines=100, bin_size=10)

    # creating a 1000Genome workflow recipe
    genome_recipe = GenomeRecipe.from_num_chromosomes(num_chromosomes=3, num_sequences=10000, num_populations=1)


All workflow recipes also provide a common method (:code:`from_num_jobs`) for
instantiating a :class:`~workflowhub.generator.workflow.abstract_recipe.WorkflowRecipe`
object as follows: ::

    from workflowhub.generator import EpigenomicsRecipe, GenomeRecipe

    # creating an Epigenomics workflow recipe
    epigenomics_recipe = EpigenomicsRecipe.from_num_jobs(num_jobs=9)

    # creating a 1000Genome workflow recipe
    genome_recipe = GenomeRecipe.from_num_jobs(num_jobs=5)

Note that :code:`num_jobs` defines the upper bound for the total number of jobs in the
workflow, and that each workflow recipe may define different lower bound values so
that the workflow structure is guaranteed. Please, refer to the :ref:`documentation of
each workflow recipe <workflowhub-generator-label>` for the lower bound values.

The current list of available workflow recipes include:

- :class:`~workflowhub.generator.workflow.cycles_recipe.CyclesRecipe`: :code:`from workflowhub.generator import CyclesRecipe`
- :class:`~workflowhub.generator.workflow.epigenomics_recipe.EpigenomicsRecipe`: :code:`from workflowhub.generator import EpigenomicsRecipe`
- :class:`~workflowhub.generator.workflow.genome_recipe.GenomeRecipe`: :code:`from workflowhub.generator import GenomeRecipe`
- :class:`~workflowhub.generator.workflow.montage_recipe.MontageRecipe`: :code:`from workflowhub.generator import MontageRecipe`
- :class:`~workflowhub.generator.workflow.seismology_recipe.SeismologyRecipe`: :code:`from workflowhub.generator import SeismologyRecipe`
- :class:`~workflowhub.generator.workflow.soykb_recipe.SoyKBRecipe`: :code:`from workflowhub.generator import SoyKBRecipe`

Generating Synthetic Workflows
------------------------------

Synthetic workflow traces are generated using the
:class:`~workflowhub.generator.generator.WorkflowGenerator` class. This
class takes as input a :class:`~workflowhub.generator.workflow.abstract_recipe.WorkflowRecipe`
object (see above), and provides two methods for generating synthetic
workflow traces:

- :meth:`~workflowhub.generator.generator.WorkflowGenerator.build_workflow`: generates a single synthetic workflow
  trace based on the workflow recipe used to instantiate the generator.
- :meth:`~workflowhub.generator.generator.WorkflowGenerator.build_workflows`: generates a number of synthetic workflow
  traces based on the workflow recipe used to instantiate the generator.

The build methods use the workflow recipe for generating realistic synthetic
workflow traces, in which the workflow structure follows workflow composition
rules defined in the workflow recipe, and jobs runtime, and input and output
data sizes are generated according to distributions obtained from actual workflow
execution traces (see :ref:`traces-label`).

Each generated trace is a represented as a :class:`~workflowhub.common.workflow.Workflow`
object (which in itself is an extension of the
`NetworkX DiGraph <https://networkx.github.io/documentation/stable/reference/classes/digraph.html>`_
class). The :class:`~workflowhub.common.workflow.Workflow` class provides two
methods for writing the generated workflow trace into files:

- :meth:`~workflowhub.common.workflow.Workflow.write_dot`: write a DOT file of a workflow trace.
- :meth:`~workflowhub.common.workflow.Workflow.write_json`: write a JSON file of a workflow trace.

Examples
--------

The following example generates a *Seismology* synthetic workflow trace
based on the number of pair of signals to estimate earthquake STFs
(:code:`num_pairs`), builds a synthetic workflow trace, and writes the
synthetic trace to a JSON file. ::

    from workflowhub import WorkflowGenerator
    from workflowhub.generator import SeismologyRecipe

    # creating a Seismology workflow recipe based on the number
    # of pair of signals to estimate earthquake STFs
    recipe = SeismologyRecipe.from_num_pairs(num_pairs=10)

    # creating an instance of the workflow generator with the
    # Seismology workflow recipe
    generator = WorkflowGenerator(recipe)

    # generating a synthetic workflow trace of the Seismology workflow
    workflow = generator.build_workflow()

    # writing the synthetic workflow trace into a JSON file
    workflow.write_json('seismology-workflow.json')


The example below generates a number of *Cycles* (agroecosystem) synthetic
workflow traces based on the upper bound number of jobs allowed per workflow. ::

    from workflowhub import WorkflowGenerator
    from workflowhub.generator import CyclesRecipe

    # creating a Cycles workflow recipe based on the number of jobs per workflow
    recipe = CyclesRecipe.from_num_jobs(num_jobs=1000)

    # creating an instance of the workflow generator with the
    # Cycles workflow recipe
    generator = WorkflowGenerator(recipe)

    # generating 10 synthetic workflow traces of the Cycles workflow
    workflows_list = generator.build_workflows(num_workflows=10)

    # writing each synthetic workflow trace into a JSON file
    count = 1
    for workflow in workflows_list:
        workflow.write_json('cycles-workflow-{:02}.json'.format(count))
        count += 1

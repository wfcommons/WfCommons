.. _generating-workflows-label:

WfGen: Generating Workflows
===========================

WfGen is the WfCommons component that targets the generation of realistic
synthetic workflow instances with a variety of characteristics. Because the
generated workflows preserve the structure and performance distributions of
real executions, you can run experiments at scales (and in quantities) that
would be impossible to obtain from production systems alone.

The :class:`~wfcommons.wfgen.generator.WorkflowGenerator` class uses workflow
recipes (as described in :ref:`workflow-recipe-generator-label`) for creating
realistic synthetic instances. The resulting workflows are represented in the
:ref:`json-format-label`, which is already supported by simulation frameworks
such as `WRENCH <https://wrench-project.org>`_.

.. _recipes-list:

Bundled Workflow Recipes
------------------------

This Python package provides *workflow recipes* for ten scientific
applications, covering bioinformatics, astronomy, seismology, and
agroecosystem domains:

.. list-table::
   :header-rows: 1
   :widths: 22 38 40

   * - Application
     - Domain
     - Import
   * - BLAST
     - Bioinformatics (sequence alignment)
     - :code:`from wfcommons.wfchef.recipes import BlastRecipe`
   * - BWA
     - Bioinformatics (read mapping)
     - :code:`from wfcommons.wfchef.recipes import BwaRecipe`
   * - Cycles
     - Agroecosystem simulation
     - :code:`from wfcommons.wfchef.recipes import CyclesRecipe`
   * - Epigenomics
     - Bioinformatics (DNA methylation)
     - :code:`from wfcommons.wfchef.recipes import EpigenomicsRecipe`
   * - 1000Genome
     - Bioinformatics (population genomics)
     - :code:`from wfcommons.wfchef.recipes import GenomeRecipe`
   * - Montage
     - Astronomy (image mosaics)
     - :code:`from wfcommons.wfchef.recipes import MontageRecipe`
   * - RNA-seq
     - Bioinformatics (transcriptomics)
     - :code:`from wfcommons.wfchef.recipes import RnaseqRecipe`
   * - Seismology
     - Seismic cross-correlation
     - :code:`from wfcommons.wfchef.recipes import SeismologyRecipe`
   * - SoyKB
     - Bioinformatics (soybean knowledge base)
     - :code:`from wfcommons.wfchef.recipes import SoykbRecipe`
   * - SRA Search
     - Bioinformatics (sequence read archive)
     - :code:`from wfcommons.wfchef.recipes import SrasearchRecipe`

You can also :ref:`create your own recipe
<workflow-recipe-generator-label>` from real instances of any other
application.

The Workflow Instances Generator
--------------------------------

Synthetic workflow instances are generated using the
:class:`~wfcommons.wfgen.generator.WorkflowGenerator` class. This class takes
as input a :class:`~wfcommons.wfgen.abstract_recipe.WorkflowRecipe` object
(see :ref:`workflow-recipe-generator-label`), and provides two methods for
generating synthetic workflow instances:

- :meth:`~wfcommons.wfgen.generator.WorkflowGenerator.build_workflow`: generates a single synthetic workflow
  instance based on the workflow recipe used to instantiate the generator.
- :meth:`~wfcommons.wfgen.generator.WorkflowGenerator.build_workflows`: generates a number of synthetic workflow
  instances based on the workflow recipe used to instantiate the generator.

The build methods use the workflow recipe for generating realistic synthetic
workflow instances, in which the workflow structure follows workflow
composition rules defined in the recipe, and task runtimes and input and
output data sizes are generated according to distributions obtained from
actual workflow execution instances (see :ref:`instances-label`).

All workflow recipes provide a common constructor, :code:`from_num_tasks`,
that defines the lower bound for the total number of tasks in the generated
synthetic workflow.

Each generated instance is represented as a
:class:`~wfcommons.common.workflow.Workflow` object (which is itself an
extension of the `NetworkX DiGraph
<https://networkx.github.io/documentation/stable/reference/classes/digraph.html>`_
class — all NetworkX graph algorithms work on it directly). The
:class:`~wfcommons.common.workflow.Workflow` class provides two methods for
writing the generated workflow instance to files:

- :meth:`~wfcommons.common.workflow.Workflow.write_json`: write a JSON file of a workflow instance (:ref:`json-format-label`).
- :meth:`~wfcommons.common.workflow.Workflow.write_dot`: write a DOT file of a workflow instance.

Scaling Runtimes and File Sizes
-------------------------------

Workflow recipes also allow the generation of synthetic workflows with
increased/reduced runtimes and/or file sizes, determined by user-provided
factors — useful for what-if experiments (e.g., "what if the input data were
50% larger?"):

- :code:`runtime_factor`: the factor by which task runtimes are increased/decreased.
- :code:`input_file_size_factor`: the factor by which task input file sizes are increased/decreased.
- :code:`output_file_size_factor`: the factor by which task output file sizes are increased/decreased.

The following example creates a Seismology workflow recipe in which task
runtime is increased by 10%, input files by 50%, and output files reduced by
20%: ::

    from wfcommons.wfchef.recipes import SeismologyRecipe

    # creating a Seismology workflow recipe with increased/decreased runtime and file sizes
    recipe = SeismologyRecipe.from_num_tasks(num_tasks=100,
                                             runtime_factor=1.1,
                                             input_file_size_factor=1.5,
                                             output_file_size_factor=0.8)

Examples
--------

The following example generates a *Seismology* synthetic workflow instance
of 250 tasks and writes it to a JSON file: ::

    import pathlib
    from wfcommons.wfchef.recipes import SeismologyRecipe
    from wfcommons import WorkflowGenerator

    generator = WorkflowGenerator(SeismologyRecipe.from_num_tasks(250))
    workflow = generator.build_workflow()
    workflow.write_json(pathlib.Path('seismology-workflow.json'))

The example below generates 10 *Blast* synthetic workflow instances for every
size defined in the array :code:`num_tasks` — 40 workflows in total: ::

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
with (at least) 1000 tasks each, and writes them to JSON files: ::

    import pathlib
    from wfcommons.wfchef.recipes import EpigenomicsRecipe
    from wfcommons import WorkflowGenerator

    generator = WorkflowGenerator(EpigenomicsRecipe.from_num_tasks(1000))
    for i, workflow in enumerate(generator.build_workflows(10)):
        workflow.write_json(pathlib.Path(f'epigenomics-workflow-{i}.json'))

The example below generates a *Cycles* (agroecosystem) synthetic workflow
instance with 250 tasks and writes it to a JSON file: ::

    import pathlib
    from wfcommons.wfchef.recipes import CyclesRecipe
    from wfcommons import WorkflowGenerator

    generator = WorkflowGenerator(CyclesRecipe.from_num_tasks(250))
    workflow = generator.build_workflow()
    workflow.write_json(pathlib.Path('cycles-workflow.json'))

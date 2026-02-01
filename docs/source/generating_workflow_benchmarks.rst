.. _generating-workflow-benchmarks-label:

WfBench: Workflow Benchmarks
============================

**WfBench** generates realistic workflow benchmark specifications that can be
translated into runnable benchmarks for current workflow systems. It produces
tasks with tunable performance characteristics (CPU, memory, and I/O usage)
and realistic dependency structures derived from production workflows.

Benchmark generation is twofold: first, a specification is produced in the
:ref:`json-format-label`; then, that specification is translated into
executable benchmark code for a target workflow system.

Generating Workflow Benchmark Specifications
--------------------------------------------

The :class:`~wfcommons.wfbench.bench.WorkflowBenchmark` class uses recipes
of workflows (as described in :ref:`workflow-recipe-generator-label`) for 
generating workflow benchmarks with an arbitrary number of tasks::

    import pathlib

    from wfcommons import BlastRecipe
    from wfcommons.wfbench import WorkflowBenchmark

    # create a workflow benchmark object to generate specifications based on a recipe
    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=500)

    # generate a specification based on performance characteristics
    path = benchmark.create_benchmark(pathlib.Path("/tmp/"), cpu_work=100, data=10, percent_cpu=0.6)

In the example above, the workflow benchmark generator first invokes the WfChef 
recipe to generate a task graph. Once the task graph has been generated, each task 
is set to be an instance of the workflow task benchmark. For each task, the following 
values for the parameters of the workflow task benchmark can be specified:

- :code:`cpu_work`: CPU work per workflow task. The :code:`cpu-benchmark` executable 
  (compiled C++) calculates an increasingly precise value of π up until the specified 
  total amount of computation (cpu_work) has been performed.
- :code:`data`: Individual data volumes for each task in a way that is coherent 
  with respect to task data dependencies (in the form of a dictionary of input 
  size files per workflow task type). Alternatively, a total data footprint (in MB)
  can be defined, i.e., the sum of the sizes of all data files read/written by 
  workflow tasks, in which case uniform I/O volumes are computed for each workflow 
  task benchmark.
- :code:`percent_cpu`: The fraction of the computation's instructions that
  correspond to non-memory operations. 

Generate from synthetic workflow instances
++++++++++++++++++++++++++++++++++++++++++

WfCommons also allows you to convert synthetic workflow instances into benchmarks directly.
The generated benchmark will have exactly the same structure as the synthetic workflow instance::

    import pathlib

    from wfcommons import BlastRecipe
    from wfcommons.wfbench import WorkflowBenchmark

    # create a synthetic workflow instance with 500 tasks or use one that you already have
    workflow = BlastRecipe.from_num_tasks(500).build_workflow()
    # create a workflow benchmark object to generate specifications based on a recipe
    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=500)
    # generate a specification based on performance characteristics and the structure of the synthetic workflow instance
    path = benchmark.create_benchmark_from_synthetic_workflow(pathlib.Path("/tmp/"), workflow, cpu_work=100, percent_cpu=0.6)

This is useful when you want to generate a benchmark with a specific structure or when you want
benchmarks with the more detailed structure provided by WfChef workflow generation.

Translating Specifications into Benchmark Code
----------------------------------------------

WfCommons provides a collection of translators that turn benchmark specifications
into runnable workflow code. All translators inherit from
:class:`~wfcommons.wfbench.translator.abstract_translator.Translator` and accept
either a :class:`~wfcommons.common.workflow.Workflow` object or a path to a
benchmark specification in :ref:`json-format-label`.

Supported translators (alphabetical)
++++++++++++++++++++++++++++++++++++

- Airflow
- Bash
- CWL
- Dask
- Makeflow
- Nextflow
- Parsl
- Pegasus
- PyCOMPSs
- Swift/T
- TaskVine

.. warning::

    WfBench leverages :code:`stress-ng` (https://github.com/ColinIanKing/stress-ng)
    to execute memory-intensive threads. Ensure :code:`stress-ng` is installed on
    all worker nodes.

Airflow
+++++++

`Apache Airflow <https://airflow.apache.org/>`_ is a platform for authoring,
scheduling, and monitoring workflows as code. Use the Airflow translator to
produce DAGs that can be executed by an Airflow scheduler::

    import pathlib
    from wfcommons import BlastRecipe
    from wfcommons.wfbench import WorkflowBenchmark, AirflowTranslator

    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=200)
    benchmark.create_benchmark(pathlib.Path("/tmp/"), cpu_work=100, data=10, percent_cpu=0.6)

    translator = AirflowTranslator(benchmark.workflow)
    translator.translate(output_folder=pathlib.Path("./airflow-wf/"))

Bash
++++

The Bash translator generates a simple, runnable shell workflow for quick local
validation and debugging::

    import pathlib
    from wfcommons import BlastRecipe
    from wfcommons.wfbench import WorkflowBenchmark, BashTranslator

    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=100)
    benchmark.create_benchmark(pathlib.Path("/tmp/"), cpu_work=50, data=5, percent_cpu=0.7)

    translator = BashTranslator(benchmark.workflow)
    translator.translate(output_folder=pathlib.Path("./bash-wf/"))

CWL
+++

`CWL <https://www.commonwl.org/>`_ is a community standard for describing command-line
tools and workflows. The CWL translator emits portable CWL definitions::

    import pathlib
    from wfcommons import BlastRecipe
    from wfcommons.wfbench import WorkflowBenchmark, CWLTranslator

    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=150)
    benchmark.create_benchmark(pathlib.Path("/tmp/"), cpu_work=80, data=8, percent_cpu=0.6)

    translator = CWLTranslator(benchmark.workflow)
    translator.translate(output_folder=pathlib.Path("./cwl-wf/"))

Dask
++++

`Dask <https://www.dask.org/>`_ is an open-source library for parallel computing
in Python. It supports local execution, HPC schedulers, and cloud environments::

    import pathlib
    from wfcommons import BlastRecipe
    from wfcommons.wfbench import WorkflowBenchmark, DaskTranslator

    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=500)
    benchmark.create_benchmark(pathlib.Path("/tmp/"), cpu_work=100, data=10, percent_cpu=0.6)

    translator = DaskTranslator(benchmark.workflow)
    translator.translate(output_folder=pathlib.Path("./dask-wf/"))

Makeflow
++++++++

`Makeflow <http://ccl.cse.nd.edu/software/makeflow/>`_ targets large, DAG-shaped
workflows on clusters, grids, and clouds. The translator emits Makeflow workflows::

    import pathlib
    from wfcommons import BlastRecipe
    from wfcommons.wfbench import WorkflowBenchmark, MakeflowTranslator

    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=200)
    benchmark.create_benchmark(pathlib.Path("/tmp/"), cpu_work=100, data=10, percent_cpu=0.6)

    translator = MakeflowTranslator(benchmark.workflow)
    translator.translate(output_folder=pathlib.Path("./makeflow-wf/"))

Nextflow
++++++++

`Nextflow <https://www.nextflow.io/>`_ enables portable, reproducible workflows
across local, HPC, and cloud environments::

    import pathlib
    from wfcommons import BlastRecipe
    from wfcommons.wfbench import WorkflowBenchmark, NextflowTranslator

    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=500)
    benchmark.create_benchmark(pathlib.Path("/tmp/"), cpu_work=100, data=10, percent_cpu=0.6)

    translator = NextflowTranslator(
        benchmark.workflow,
        use_subworkflows=False,
        max_tasks_per_subworkflow=1000,
    )
    translator.translate(output_folder=pathlib.Path("./nextflow-wf/"))

If you want to split large workflows across multiple Nextflow module files, enable
subworkflows and set the maximum number of tasks per module. This produces a
``modules/`` directory plus a top-level ``workflow.nf`` that includes and runs
the modules sequentially::

    translator = NextflowTranslator(
        benchmark.workflow,
        use_subworkflows=True,
        max_tasks_per_subworkflow=250,
    )
    translator.translate(output_folder=pathlib.Path("./nextflow-wf/"))

.. warning::

    Nextflow does not support tasks with iterations (tasks that depend on another
    instance of the same abstract task). Translation fails for workflows that
    include iterations.

.. note::

    If you plan to run Nextflow on an HPC system using Slurm, we **strongly
    recommend** using the `HyperQueue <https://github.com/It4innovations/hyperqueue>`_
    executor. HyperQueue efficiently distributes workflow tasks across all allocated
    compute nodes, improving scalability and resource utilization.

    The :class:`~wfcommons.wfbench.translator.nextflow.NextflowTranslator`
    class includes functionality to automatically generate a Slurm script
    template for running the workflow on HPC systems.

Parsl
+++++

`Parsl <https://parsl-project.org/>`_ is a parallel scripting library for Python.
The translator emits a Parsl workflow suitable for local or distributed execution::

    import pathlib
    from wfcommons import BlastRecipe
    from wfcommons.wfbench import WorkflowBenchmark, ParslTranslator

    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=200)
    benchmark.create_benchmark(pathlib.Path("/tmp/"), cpu_work=100, data=10, percent_cpu=0.6)

    translator = ParslTranslator(benchmark.workflow)
    translator.translate(output_folder=pathlib.Path("./parsl-wf/"))

Pegasus
+++++++

`Pegasus <http://pegasus.isi.edu>`_ orchestrates complex scientific workflows on
clusters, grids, and clouds by mapping tasks onto distributed resources::

    import pathlib
    from wfcommons import BlastRecipe
    from wfcommons.wfbench import WorkflowBenchmark, PegasusTranslator

    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=500)
    benchmark.create_benchmark(pathlib.Path("/tmp/"), cpu_work=100, data=10, percent_cpu=0.6)

    translator = PegasusTranslator(benchmark.workflow)
    translator.translate(output_folder=pathlib.Path("./pegasus-wf/"))

.. warning::

    Pegasus uses `HTCondor <https://htcondor.org/>`_ to orchestrate tasks. By
    default, HTCondor does not implement CPU affinity for program threads.
    To enable CPU affinity, specify :code:`lock_files_folder` when using
    :meth:`~wfcommons.wfbench.bench.WorkflowBenchmark.create_benchmark`.

PyCOMPSs
++++++++

`PyCOMPSs <https://compss.bsc.es/>`_ is a programming model and runtime for
parallel Python applications on distributed infrastructures::

    import pathlib
    from wfcommons import CyclesRecipe
    from wfcommons.wfbench import WorkflowBenchmark, PyCompssTranslator

    benchmark = WorkflowBenchmark(recipe=CyclesRecipe, num_tasks=200)
    benchmark.create_benchmark(pathlib.Path("/tmp/"), cpu_work=500, data=1000, percent_cpu=0.8)

    translator = PyCompssTranslator(benchmark.workflow)
    translator.translate(output_folder=pathlib.Path("./pycompss-wf/"))

Swift/T
+++++++

`Swift/T <http://swift-lang.org/Swift-T/>`_ is a workflow system for HPC environments,
designed to scale to large task graphs::

    import pathlib
    from wfcommons import BlastRecipe
    from wfcommons.wfbench import WorkflowBenchmark, SwiftTTranslator

    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=500)
    benchmark.create_benchmark(pathlib.Path("/tmp/"), cpu_work=100, data=10, percent_cpu=1.0)

    translator = SwiftTTranslator(benchmark.workflow)
    translator.translate(output_folder=pathlib.Path("./swift-t-wf/"))

TaskVine
++++++++

`TaskVine <https://ccl.cse.nd.edu/software/taskvine/>`_ is a task scheduler for
data-intensive dynamic workflows across HPC clusters, GPU clusters, and clouds::

    import pathlib
    from wfcommons import BlastRecipe
    from wfcommons.wfbench import WorkflowBenchmark, TaskVineTranslator

    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=500)
    benchmark.create_benchmark(save_dir=pathlib.Path("/tmp/"), cpu_work=100, data=10, percent_cpu=1.0)

    translator = TaskVineTranslator(benchmark.workflow)
    translator.translate(output_folder=pathlib.Path("./taskvine-wf/"))

WfBench will generate a folder containing the TaskVine workflow
:code:`taskvine_workflow.py`, workflow input data (:code:`./taskvine-wf/data/`),
workflow binaries (:code:`./taskvine-wf/bin/`), and the Poncho package specification
(:code:`./taskvine-wf/taskvine_poncho.json`).

.. warning::
    This TaskVine workflow requires :code:`stress-ng` to be installed and accessible
    in the system's :code:`$PATH` where the manager runs.

.. _generating-workflow-benchmarks-label:

WfBench: Workflow Benchmarks
============================

**WfBench** is a generator of realistic workflow benchmark specifications that 
can be translated into benchmark code to be executed with current workflow 
systems. it generates workflow tasks with arbitrary performance characteristics 
(CPU, memory, and I/O usage) and with realistic task dependency structures 
based on those seen in production workflows.

The generation of workflow benchmakrs is twofold. First, a realistic workflow 
benchmark specification is generated in the :ref:`json-format-label`. Then, 
this specification is translated into benchmark code to be executed with a 
workflow system.

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
The generated benchmark will have exactly the same structure as the synthetic workflow instance.

    import pathlib

    from wfcommons import BlastRecipe
    from wfcommons.wfbench import WorkflowBenchmark

    # create a synthetic workflow instance with 500 tasks or use one that you already have
    workflow = BlastRecipe.from_num_tasks(500).build_workflow()
    # create a workflow benchmark object to generate specifications based on a recipe
    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=500)
    # generate a specification based on performance characteristics and the structure of the synthetic workflow instance
    path = benchmark.create_benchmark_from_synthetic_workflow(pathlib.Path("/tmp/"), workflow, cpu_work=100, data=10, percent_cpu=0.6)

This is useful when you want to generate a benchmark with a specific structure or when you want
benchmarks with the more detailed structure provided by WfChef workflow generation.

Translating Specifications into Benchmark Codes
-----------------------------------------------

WfCommons provides a collection of translators for executing the benchmarks as actual
workflow applications. Below, we provide illustrative examples on how to generate 
workflow benchmarks for the currently supported workflow systems.

The :class:`~wfcommons.wfbench.translator.abstract_translator.Translator` class is 
the foundation for each translator class. This class takes as input either a 
:class:`~wfcommons.common.workflow.Workflow` object or a path to a workflow benchmark
description in :ref:`json-format-label`.

.. warning::
    
    WfBench leverages :code:`stress-ng` (https://github.com/ColinIanKing/stress-ng) 
    to execute memory-intensive threads. Therefore, it is crucial to ensure that 
    :code:`stress-ng` is installed on all worker nodes.

Nextflow
++++++++
`Nextflow <https://www.nextflow.io/>`_ is a workflow management system that enables
the development of portable and reproducible workflows. It supports deploying workflows
on a variety of execution platforms including local, HPC schedulers, and cloud-based
and container-based environments. Below, we provide an example on how to generate
workflow benchmark for running with Nextflow:

    import pathlib

    from wfcommons import BlastRecipe
    from wfcommons.wfbench import WorkflowBenchmark, NextflowTranslator

    # create a workflow benchmark object to generate specifications based on a recipe
    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=500)

    # generate a specification based on performance characteristics
    benchmark.create_benchmark(pathlib.Path("/tmp/"), cpu_work=100, data=10, percent_cpu=0.6)

    # generate a Nextflow workflow
    translator = NextflowTranslator(benchmark.workflow)
    translator.translate(output_file_name=pathlib.Path("/tmp/benchmark-workflow.nf"))

.. warning::

    Nextflow's way of defining workflows does not support tasks with iterations i.e. tasks 
    that depend on another instance of the same abstract task. Thus, the translator
    fails when you try to translate a workflow with iterations.

Pegasus
+++++++

`Pegasus <http://pegasus.isi.edu>`_ orchestrates the execution of complex scientific 
workflows by providing a platform to define, organize, and automate computational 
tasks and data dependencies. Pegasus handles the complexity of large-scale workflows 
by automatically mapping tasks onto distributed computing resources, such as clusters, 
grids, or clouds. Below, we provide an example on how to generate workflow benchmark 
for running with Pegasus::

    import pathlib

    from wfcommons import BlastRecipe
    from wfcommons.wfbench import WorkflowBenchmark, PegasusTranslator

    # create a workflow benchmark object to generate specifications based on a recipe
    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=500)

    # generate a specification based on performance characteristics
    benchmark.create_benchmark(pathlib.Path("/tmp/"), cpu_work=100, data=10, percent_cpu=0.6)

    # generate a Pegasus workflow
    translator = PegasusTranslator(benchmark.workflow)
    translator.translate(output_file_name=pathlib.Path("/tmp/benchmark-workflow.py"))

.. warning::

    Pegasus utilizes the `HTCondor <https://htcondor.org/>`_ framework to orchestrate 
    the execution of workflow tasks. By default, HTCondor does not implement CPU affinity 
    for program threads. However, WfBench offers an extra capability to enforce CPU 
    affinity during benchmark execution. To enable this feature, you need to specify 
    the :code:`lock_files_folder` parameter when using 
    :meth:`~wfcommons.wfbench.bench.WorkflowBenchmark.create_benchmark`.

Swift/T
+++++++

`Swift/T <http://swift-lang.org/Swift-T/>`_ is an advanced workflow system designed 
specifically for high-performance computing (HPC) environments. It dynamically manages 
task dependencies and resource allocation, enabling efficient utilization of HPC 
systems. It provides a seamless interface to diverse tools, libraries, and scientific 
applications, making it easy to integrate existing codes into workflows. Below, we 
provide an example on how to generate workflow benchmark for running with Swift/T::

    import pathlib

    from wfcommons import BlastRecipe
    from wfcommons.wfbench import WorkflowBenchmark, SwiftTTranslator

    # create a workflow benchmark object to generate specifications based on a recipe
    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=500)

    # generate a specification based on performance characteristics
    benchmark.create_benchmark(pathlib.Path("/tmp/"), cpu_work=100, data=10, percent_cpu=0.6)

    # generate a Swift/T workflow
    translator = SwiftTTranslator(benchmark.workflow)
    translator.translate(output_file_name=pathlib.Path("/tmp/benchmark-workflow.swift"))

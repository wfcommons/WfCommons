.. _generating-workflow-benchmarks-label:

Generating Workflow Benchmarks
==============================

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
generating workflow benchmarks.::

    import pathlib

    from wfcommons import BlastRecipe
    from wfcommons.wfbench import WorkflowBenchmark

    # create a workflow benchmark object to generate specifications based on a recipe
    benchmark = WorkflowBenchmark(recipe=BlastRecipe, num_tasks=500)
    # generate a specification based on performance characteristics
    path = benchmark.create_benchmark(pathlib.Path("/tmp/"), cpu_work=100, data=10, percent_cpu=0.6)

Translating Specifications into Benchmark Codes
-----------------------------------------------

WfCommons provides a collection of translators for executing the benchmarks as actual
workflow applications. Below, we provide illustrative examples on how to generate 
workflow benchmakrs for the currently supported workflow systems.

The :class:`~wfcommons.wfbench.translator.abstract_translator.Translator` class is 
the foundation for each translator class. This class takes as input either a 
:class:`~wfcommons.common.workflow.Workflow` object or a path to a workflow benchmark
description in :ref:`json-format-label`.

Pegasus
+++++++

`Pegasus <http://pegasus.isi.edu>`_ orchestrates the execution of complex scientific 
workflows by providing a platform to define, organize, and automate computational 
tasks and data dependencies. Pegasus handles the complexity of large-scale workflows 
by automatically mapping tasks onto distributed computing resources, such as clusters, 
grids, or clouds. Below, we provide an example on how to generate workflow benchmark 
for running with Pegasus:::

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

Swift/T
+++++++

`Swift/T <http://swift-lang.org/Swift-T/>`_ is an advanced workflow system designed 
specifically for high-performance computing (HPC) environments. It dynamically manages 
task dependencies and resource allocation, enabling efficient utilization of HPC 
systems. It provides a seamless interface to diverse tools, libraries, and scientific 
applications, making it easy to integrate existing codes into workflows. Below, we 
provide an example on how to generate workflow benchmark for running with Swift/T:::

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

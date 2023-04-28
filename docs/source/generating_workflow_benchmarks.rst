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

Generating a Workflow Benchmark Specification
---------------------------------------------

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
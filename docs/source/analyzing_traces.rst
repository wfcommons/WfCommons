.. _traces-label:

Analyzing Traces
================

Workflow execution traces have been widely used to profile and characterize
workflow executions, and to build distributions of workflow execution behaviors,
which are used to evaluate methods and techniques in simulation or in real
conditions.

The first axis of the WorkflowHub project targets the analysis of actual workflow
execution traces (i.e., the workflow execution profile data and characterizations)
in order to build **recipes** of workflow applications. These recipes contain
the necessary information for generating synthetic, yet realistic, workflow
traces that resemble the structure and distribution of the original workflow
executions.

Workflow Execution Traces
-------------------------

A workflow execution trace represents an actual execution of a scientific
workflow on a distributed platform (e.g., clouds, grids, HPC, etc.). In the
WorkflowHub project, a trace is represented in a JSON file following the
schema described in :ref:`json-format-label` section. This Python package
provides a *trace loader* tool for importing workflow execution traces
for analysis. For instance, the code snippet below shows how a trace can
be loaded using the :class:`~workflowhub.trace.trace.Trace` class: ::

    from workflowhub import Trace
    trace = Trace(input_trace='/path/to/trace/file.json')

The :class:`~workflowhub.trace.trace.Trace` class provides a number of
methods for interacting with the workflow trace, including:

- :meth:`~workflowhub.trace.trace.Trace.draw`: produces an image or a pdf file representing the trace.
- :meth:`~workflowhub.trace.trace.Trace.leaves`: gets the leaves of the workflow (i.e., the tasks without any successors).
- :meth:`~workflowhub.trace.trace.Trace.roots`: gets the roots of the workflow (i.e., the tasks without any predecessors).
- :meth:`~workflowhub.trace.trace.Trace.write_dot`: writes a dot file of the trace.

The Trace Analyzer
------------------

The :class:`~workflowhub.trace.trace_analyzer.TraceAnalyzer` class provides
a number of tools for analyzing collection of workflow execution traces. The
goal of the :class:`~workflowhub.trace.trace_analyzer.TraceAnalyzer` is to
perform analyzes of one or multiple workflow execution traces, and build
summaries of the analyzes per workflow' job type prefix.

.. note::

    Although any workflow execution trace represented as a
    :class:`~workflowhub.trace.trace.Trace` object (i.e., compatible with
    :ref:`json-format-label`) can be appended to the
    :class:`~workflowhub.trace.trace_analyzer.TraceAnalyzer`, we strongly
    recommend that only traces of a single workflow application type be
    appended to an analyzer object. You may though create several analyzer
    objects per workflow application.



Examples
--------

The following example shows the analysis of a set of traces, stored in a local folder,
of a Seismology workflow. In this example, ::

    from workflowhub import Trace, TraceAnalyzer
    from os import listdir
    from os.path import isfile, join

    # obtaining list of trace files in the folder
    TRACES_PATH = "/Path/to/some/trace/folder/"
    trace_files = [f for f in listdir(TRACES_PATH) if isfile(join(TRACES_PATH, f))]

    # creating the trace analyzer object
    analyzer = TraceAnalyzer()

    # appending trace files to the trace analyzer
    for trace_file in trace_files:
        trace = Trace(input_trace=TRACES_PATH + trace_file)
        analyzer.append_trace(trace)

    # list of workflow job name prefixes to be analyzed in each trace
    workflow_jobs = ['sG1IterDecon', 'wrapper_siftSTFByMisfit']

    # building the trace summary
    traces_summary = analyzer.build_summary(workflow_jobs, include_raw_data=True)

    # generating all fit plots (runtime, and input and output files)
    analyzer.generate_all_fit_plots()

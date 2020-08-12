Analyzing Traces
================

Examples
--------

The following example shows the analysis of a set of traces, stored in a local folder,
of a Seismology workflow. In the example, ::

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

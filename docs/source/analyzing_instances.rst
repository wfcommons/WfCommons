.. _instances-label:

Analyzing Instances
===================

Workflow execution instances have been widely used to profile and characterize
workflow executions, and to build distributions of workflow execution behaviors,
which are used to evaluate methods and techniques in simulation or in real
conditions.

The WfCommons project targets the analysis of actual workflow execution instances
(i.e., the workflow execution profile data and characterizations)
in order to build :ref:`workflow-recipe-label` of workflow applications.
These recipes contain the necessary information for generating synthetic, yet
realistic, workflow instances that resemble the structure and distribution of
the original workflow executions.

A `list of workflow execution instances <https://wfcommons.org/instances>`_
that are compatible with :ref:`json-format-label` is kept constantly updated
in our project website.

.. _wfinstances-label:

WfInstances
-----------

A workflow execution instance represents an actual execution of a scientific
workflow on a distributed platform (e.g., clouds, grids, HPC, etc.). In the
WfCommons project, an instance is represented in a JSON file following the
schema described in :ref:`json-format-label`. This Python package
provides an *instance loader* tool for importing workflow execution instances
for analysis. For instance, the code snippet below shows how an instance can
be loaded using the :class:`~wfcommons.trace.trace.Trace` class: ::

    from wfcommons import Trace
    trace = Trace(input_trace='/path/to/instance/file.json')

The :class:`~wfcommons.trace.trace.Trace` class provides a number of
methods for interacting with the workflow instance, including:

- :meth:`~wfcommons.trace.trace.Trace.draw`: produces an image or a pdf file representing the instance.
- :meth:`~wfcommons.trace.trace.Trace.leaves`: gets the leaves of the workflow (i.e., the tasks without any successors).
- :meth:`~wfcommons.trace.trace.Trace.roots`: gets the roots of the workflow (i.e., the tasks without any predecessors).
- :meth:`~wfcommons.trace.trace.Trace.write_dot`: writes a dot file of the instance.

.. note::
    Although the analysis methods are inherently used by WfCommons (specifically
    WfChef) for :ref:`generating-workflows-recipe-label`, they can also be used
    in a standalone manner.

The Instance Analyzer
---------------------

The :class:`~wfcommons.trace.trace_analyzer.TraceAnalyzer` class provides
a number of tools for analyzing collection of workflow execution instances. The
goal of the :class:`~wfcommons.trace.trace_analyzer.TraceAnalyzer` is to
perform analyzes of one or multiple workflow execution instances, and build
summaries of the analyzes per workflow' task type prefix.

.. warning::

    Although any workflow execution instance represented as a
    :class:`~wfcommons.trace.trace.Trace` object (i.e., compatible with
    :ref:`json-format-label`) can be appended to the
    :class:`~wfcommons.trace.trace_analyzer.TraceAnalyzer`, we strongly
    recommend that only instances of a single workflow application type be
    appended to an analyzer object. You may though create several analyzer
    objects per workflow application.

The :meth:`~wfcommons.trace.trace_analyzer.TraceAnalyzer.append_trace` method
allows you to include instances for analysis. The
:meth:`~wfcommons.trace.trace_analyzer.TraceAnalyzer.build_summary` method
processes all appended instances. The method applies probability distributions fitting
to a series of data to find the *best* (i.e., minimizes the mean square error)
probability distribution that represents the analyzed data. The method returns
a summary of the analysis of instances in the form of a Python dictionary object in
which keys are task prefixes (provided when invoking the method) and values
describe the best probability distribution fit for tasks' runtime, and input and
output data file sizes. The code excerpt below shows an example of an analysis
summary showing the best fit probability distribution for runtime of the
:code:`individuals` tasks (1000Genome workflow): ::

    "individuals": {
        "runtime": {
            "min": 48.846,
            "max": 192.232,
            "distribution": {
                "name": "skewnorm",
                "params": [
                    11115267.652937062,
                    -2.9628504044929433e-05,
                    56.03957070238482
                ]
            }
        },
        ...
    }

Workflow analysis summaries are used by WfChef to develop :ref:`workflow-recipe-label`,
in which themselves are used to :ref:`generate realistic synthetic workflow instances
<generating-workflows-label>`.

Probability distribution fits can also be plotted by using the
:meth:`~wfcommons.trace.trace_analyzer.TraceAnalyzer.generate_fit_plots` or
:meth:`~wfcommons.trace.trace_analyzer.TraceAnalyzer.generate_all_fit_plots`
methods -- plots will be saved as :code:`png` files.

Examples
--------

The following example shows the analysis of a set of instances, stored in a local folder,
of a Seismology workflow. In this example, we seek for finding the best probability
distribution fitting for task *prefixes* of the Seismology workflow
(:code:`sG1IterDecon`, and :code:`wrapper_siftSTFByMisfit`), and generate all fit
plots (runtime, and input and output files) into the :code:`fits` folder using
:code:`seismology` as a prefix for each generated plot: ::

    from wfcommons import Trace, TraceAnalyzer
    from os import listdir
    from os.path import isfile, join

    # obtaining list of instance files in the folder
    INSTANCES_PATH = "/path/to/some/instance/folder/"
    instance_files = [f for f in listdir(INSTANCES_PATH) if isfile(join(INSTANCES_PATH, f))]

    # creating the instance analyzer object
    analyzer = TraceAnalyzer()

    # appending instance files to the instance analyzer
    for instance_file in instance_files:
        instance = Trace(input_trace=INSTANCES_PATH + instance_file)
        analyzer.append_trace(instance)

    # list of workflow task name prefixes to be analyzed in each instance
    workflow_tasks = ['sG1IterDecon', 'wrapper_siftSTFByMisfit']

    # building the instance summary
    instances_summary = analyzer.build_summary(workflow_tasks, include_raw_data=True)

    # generating all fit plots (runtime, and input and output files)
    analyzer.generate_all_fit_plots(outfile_prefix='fits/seismology')

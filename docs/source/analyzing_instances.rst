.. _instances-label:

WfInstances: Workflow Instances
===============================

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

A `list of workflow execution instances <https://github.com/wfcommons/wfinstances>`_
that are compatible with :ref:`json-format-label` is kept constantly updated
in our project GitHub.

.. _wfinstances-label:

WfInstances
-----------

A workflow execution instance represents an actual execution of a scientific
workflow on a distributed platform (e.g., clouds, grids, HPC, etc.). In the
WfCommons project, an instance is represented in a JSON file following the
schema described in :ref:`json-format-label`. This Python package
provides an *instance loader* tool for importing workflow execution instances
for analysis. For instance, the code snippet below shows how an instance can
be loaded using the :class:`~wfcommons.wfinstances.instance.Instance` class: ::

    import pathlib
    from wfcommons import Instance
    input_instance = pathlib.Path('/path/to/instance/file.json')
    instance = Instance(input_instance=input_instance)

The :class:`~wfcommons.wfinstances.instance.Instance` class provides a number of
methods for interacting with the workflow instance, including:

- :meth:`~wfcommons.wfinstances.instance.Instance.draw`: produces an image or a pdf file representing the instance.
- :meth:`~wfcommons.wfinstances.instance.Instance.leaves`: gets the leaves of the workflow (i.e., the tasks without any successors).
- :meth:`~wfcommons.wfinstances.instance.Instance.roots`: gets the roots of the workflow (i.e., the tasks without any predecessors).
- :meth:`~wfcommons.wfinstances.instance.Instance.write_dot`: writes a dot file of the instance.

.. note::
    Although the analysis methods are inherently used by WfCommons (specifically
    WfChef) for :ref:`generating-workflows-recipe-label`, they can also be used
    in a standalone manner.

Parsing Workflow Execution Logs
-------------------------------

The most common way for obtaining **workflow instances** from actual workflow
executions is to parse execution logs. As part of the WfCommons project, we
are constantly developing parsers for commonly used workflow management systems.
The parsers provided in this Python package automatically scans execution logs
to produce instances using :ref:`json-format-label`.

Each parser class is derived from the abstract
:class:`~wfcommons.wfinstances.logs.abstract_logs_parser.LogsParser` class. Thus, each
parser provides a
:meth:`~wfcommons.wfinstances.logs.abstract_logs_parser.LogsParser.build_workflow`
method.

Supported log parsers
+++++++++++++++++++++

- :class:`~wfcommons.wfinstances.logs.pegasusrec.HierarchicalPegasusLogsParser`
- :class:`~wfcommons.wfinstances.logs.makeflow.MakeflowLogsParser`
- :class:`~wfcommons.wfinstances.logs.nextflow.NextflowLogsParser`
- :class:`~wfcommons.wfinstances.logs.pegasus.PegasusLogsParser`
- :class:`~wfcommons.wfinstances.logs.taskvine.TaskVineLogsParser`

Examples
++++++++

Hierarchical Pegasus
++++++++++++++++++++

This parser targets Pegasus submit directories that contain hierarchical workflows.
It recursively parses sub-workflows and rebuilds a coherent workflow instance::

    import pathlib
    from wfcommons.wfinstances import HierarchicalPegasusLogsParser

    submit_dir = pathlib.Path('/path/to/pegasus/hierarchical/submit/dir/')
    parser = HierarchicalPegasusLogsParser(submit_dir=submit_dir)
    workflow = parser.build_workflow('pegasus-hierarchical-workflow-test')
    workflow.write_json(pathlib.Path('./pegasus-hierarchical-workflow.json'))

Makeflow
++++++++

`Makeflow <http://ccl.cse.nd.edu/software/makeflow/>`_ is a workflow system for
executing large complex workflows on clusters, clouds, and grids. The Makeflow
language is similar to traditional "Make", so if you can write a Makefile, then you
can write a Makeflow. A workflow can be just a few commands chained together, or
it can be a complex application consisting of thousands of tasks. It can have an
arbitrary DAG structure and is not limited to specific patterns. Makeflow is used
on a daily basis to execute complex scientific applications in fields such as data
mining, high energy physics, image processing, and bioinformatics. It has run on
campus clusters, the Open Science Grid, NSF XSEDE machines, NCSA Blue Waters, and
Amazon Web Services. Makeflow logs provide time-stamped event instances from these
executions. The following example shows the analysis of Makeflow execution logs,
stored in a local folder (:code:`execution_dir`), for a workflow execution using the
:class:`~wfcommons.wfinstances.logs.makeflow.MakeflowLogsParser` class: ::

    import pathlib
    from wfcommons.wfinstances import MakeflowLogsParser

    # creating the parser for the Makeflow workflow
    execution_dir = pathlib.Path('/path/to/makeflow/execution/dir/blast/chameleon-small-001/')
    resource_monitor_logs_dir = pathlib.Path('/path/to/makeflow/resource/monitor/logs/dir')
    parser = MakeflowLogsParser(execution_dir=execution_dir,
                                resource_monitor_logs_dir=resource_monitor_logs_dir)

    # generating the workflow instance object
    workflow = parser.build_workflow('makeflow-workflow-test')

    # writing the workflow instance to a JSON file
    workflow_path = pathlib.Path('./makeflow-workflow.json')
    workflow.write_json(workflow_path)

.. note::
    The :class:`~wfcommons.wfinstances.logs.makeflow.MakeflowLogsParser` class requires
    that Makeflow workflows to run with the
    `Resource Monitor <https://cctools.readthedocs.io/en/latest/resource_monitor/>`_
    tool (e.g., execute the workflow using the :code:`--monitor=logs`).

Nextflow
++++++++

`Nextflow <https://nextflow.io>`_ is a reactive workflow framework and a programming DSL
that eases the writing of data-intensive computational pipelines. It is designed around
the idea that the Linux platform is the lingua franca of data science. Linux provides
many simple but powerful command-line and scripting tools that, when chained together,
facilitate complex data manipulations. Nextflow extends this approach, adding the ability
to define complex program interactions and a high-level parallel computational environment
based on the dataflow programming model. The following example shows the analysis of
Nextflow execution logs, stored in a local folder (:code:`execution_dir`), for a workflow
execution using the :class:`~wfcommons.wfinstances.logs.nextflow.NextflowLogsParser` class: ::

    import pathlib
    from wfcommons.wfinstances import NextflowLogsParser

    # creating the parser for the Nextflow workflow
    execution_dir = pathlib.Path('/path/to/nextflow/execution/dir/')
    parser = NextflowLogsParser(execution_dir=execution_dir)

    # generating the workflow instance object
    workflow = parser.build_workflow('nextflow-workflow-test')

    # writing the workflow instance to a JSON file
    workflow_path = pathlib.Path('./nextflow-workflow.json')
    workflow.write_json(workflow_path)

.. note::
    The :class:`~wfcommons.wfinstances.logs.nextflow.NextflowLogsParser` class assumes
    that workflow executions will produce an :code:`execution_report_*.html` and an
    :code:`execution_timeline_*.html` files.

Pegasus
+++++++

`Pegasus <http://pegasus.isi.edu>`_ is being used in production to execute workflows
for dozens of high-profile applications in a wide range of scientific domains. Pegasus
provides the necessary abstractions for scientists to create workflows and allows for
transparent execution of these workflows on a range of compute platforms including
clusters, clouds, and national cyberinfrastructures. Workflow execution with Pegasus
includes data management, monitoring, and failure handling, and is managed by HTCondor
DAGMan. Individual workflow tasks are managed by a workload management framework,
HTCondor, which supervises task executions on local and remote resources. Pegasus
logs provide time-stamped event instances from these executions. The following example shows
the analysis of Pegasus execution logs, stored in a local folder (:code:`submit_dir`), for a
workflow execution using the :class:`~wfcommons.wfinstances.logs.pegasus.PegasusLogsParser`
class: ::

    import pathlib
    from wfcommons.wfinstances import PegasusLogsParser

    # creating the parser for the Pegasus workflow
    submit_dir = pathlib.Path('/path/to/pegasus/submit/dir/seismology/chameleon-100p-001/')
    parser = PegasusLogsParser(submit_dir=submit_dir)

    # generating the workflow instance object
    workflow = parser.build_workflow('pegasus-workflow-test')

    # writing the workflow instance to a JSON file
    workflow_path = pathlib.Path('./pegasus-workflow.json')
    workflow.write_json(workflow_path)

TaskVine
++++++++

`TaskVine <https://ccl.cse.nd.edu/software/taskvine/>`_ is a task scheduler for
data-intensive dynamic workflows. The TaskVine logs parser translates TaskVine
execution logs into workflow instances compatible with :ref:`json-format-label`::

    import pathlib
    from wfcommons.wfinstances import TaskVineLogsParser

    execution_dir = pathlib.Path('/path/to/taskvine/execution/dir/')
    parser = TaskVineLogsParser(execution_dir=execution_dir)
    workflow = parser.build_workflow('taskvine-workflow-test')
    workflow.write_json(pathlib.Path('./taskvine-workflow.json'))

The Instance Analyzer
---------------------

The :class:`~wfcommons.wfinstances.instance_analyzer.InstanceAnalyzer` class provides
a number of tools for analyzing collection of workflow execution instances. The
goal of the :class:`~wfcommons.wfinstances.instance_analyzer.InstanceAnalyzer` is to
perform analyzes of one or multiple workflow execution instances, and build
summaries of the analyzes per workflow' task type prefix.

.. warning::

    Although any workflow execution instance represented as a
    :class:`~wfcommons.wfinstances.instance.Instance` object (i.e., compatible with
    :ref:`json-format-label`) can be appended to the
    :class:`~wfcommons.wfinstances.instance_analyzer.InstanceAnalyzer`, we strongly
    recommend that only instances of a single workflow application type be
    appended to an analyzer object. You may though create several analyzer
    objects per workflow application.

The :meth:`~wfcommons.wfinstances.instance_analyzer.InstanceAnalyzer.append_instance` method
allows you to include instances for analysis. The
:meth:`~wfcommons.wfinstances.instance_analyzer.InstanceAnalyzer.build_summary` method
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
:meth:`~wfcommons.wfinstances.instance_analyzer.InstanceAnalyzer.generate_fit_plots` or
:meth:`~wfcommons.wfinstances.instance_analyzer.InstanceAnalyzer.generate_all_fit_plots`
methods -- plots will be saved as :code:`png` files.

Examples
--------

The following example shows the analysis of a set of instances, stored in a local folder,
of a Seismology workflow. In this example, we seek for finding the best probability
distribution fitting for task *prefixes* of the Seismology workflow
(:code:`sG1IterDecon`, and :code:`wrapper_siftSTFByMisfit`), and generate all fit
plots (runtime, and input and output files) into the :code:`fits` folder using
:code:`seismology` as a prefix for each generated plot: ::

    import pathlib
    from wfcommons import Instance, InstanceAnalyzer

    # obtaining list of instance files in the folder
    INSTANCES_PATH = pathlib.Path('/path/to/some/instance/folder/')
    instance_files = [f for f in INSTANCES_PATH.glob('*') if INSTANCES_PATH.joinpath(f).is_file()]

    # creating the instance analyzer object
    analyzer = InstanceAnalyzer()

    # appending instance files to the instance analyzer
    for instance_file in instance_files:
        instance = Instance(input_instance=INSTANCES_PATH.joinpath(instance_file))
        analyzer.append_instance(instance)

    # list of workflow task name prefixes to be analyzed in each instance
    workflow_tasks = ['sG1IterDecon', 'wrapper_siftSTFByMisfit']

    # building the instance summary
    instances_summary = analyzer.build_summary(workflow_tasks, include_raw_data=True)

    # generating all fit plots (runtime, and input and output files)
    analyzer.generate_all_fit_plots(outfile_prefix='fits/seismology')

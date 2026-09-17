.. _instances-label:

WfInstances: Workflow Instances
===============================

Workflow execution instances are widely used to profile and characterize
workflow executions, and to build distributions of workflow execution
behaviors, which are used to evaluate methods and techniques in simulation or
in real conditions.

The WfCommons project targets the analysis of actual workflow execution
instances (i.e., the workflow execution profile data and characterizations) in
order to build :ref:`workflow-recipe-label` of workflow applications. These
recipes contain the necessary information for generating synthetic, yet
realistic, workflow instances that resemble the structure and distribution of
the original workflow executions.

With the tools on this page you can:

- **Load** any workflow instance in :ref:`json-format-label` and inspect its
  task graph programmatically.
- **Parse** execution logs from production workflow systems into WfFormat
  instances.
- **Analyze** collections of instances to obtain statistical characterizations
  of task runtimes and data sizes.

A `list of workflow execution instances <https://github.com/wfcommons/wfinstances>`_
that are compatible with :ref:`json-format-label` is kept constantly updated
in our project GitHub — you can start from those instead of running workflows
yourself.

.. _wfinstances-label:

Loading Workflow Instances
--------------------------

A workflow execution instance represents an actual execution of a scientific
workflow on a distributed platform (e.g., clouds, grids, HPC). In the
WfCommons project, an instance is represented in a JSON file following the
schema described in :ref:`json-format-label`. This Python package provides an
*instance loader* tool for importing workflow execution instances for
analysis. The code snippet below shows how an instance can be loaded using the
:class:`~wfcommons.wfinstances.instance.Instance` class: ::

    import pathlib
    from wfcommons import Instance

    input_instance = pathlib.Path('/path/to/instance/file.json')
    instance = Instance(input_instance=input_instance)

The :class:`~wfcommons.wfinstances.instance.Instance` class provides a number
of methods for interacting with the workflow instance, including:

- :meth:`~wfcommons.wfinstances.instance.Instance.draw`: produces an image or a pdf file representing the instance.
- :meth:`~wfcommons.wfinstances.instance.Instance.leaves`: gets the leaves of the workflow (i.e., the tasks without any successors).
- :meth:`~wfcommons.wfinstances.instance.Instance.roots`: gets the roots of the workflow (i.e., the tasks without any predecessors).
- :meth:`~wfcommons.wfinstances.instance.Instance.write_dot`: writes a dot file of the instance.

.. note::
    Although the analysis methods are inherently used by WfCommons
    (specifically WfChef) for :ref:`generating-workflows-recipe-label`, they
    can also be used in a standalone manner.

Parsing Workflow Execution Logs
-------------------------------

The most common way of obtaining **workflow instances** from actual workflow
executions is to parse execution logs. As part of the WfCommons project, we
are constantly developing parsers for commonly used workflow management
systems. The parsers provided in this Python package automatically scan
execution logs to produce instances using :ref:`json-format-label`.

Each parser class is derived from the abstract
:class:`~wfcommons.wfinstances.logs.abstract_logs_parser.LogsParser` class, and
provides a
:meth:`~wfcommons.wfinstances.logs.abstract_logs_parser.LogsParser.build_workflow`
method that returns a :class:`~wfcommons.common.workflow.Workflow` object,
which can then be written to a WfFormat JSON file with
:meth:`~wfcommons.common.workflow.Workflow.write_json`.

Supported log parsers
+++++++++++++++++++++

.. list-table::
   :header-rows: 1
   :widths: 18 42 40

   * - System
     - Parser class
     - Main input
   * - Makeflow
     - :class:`~wfcommons.wfinstances.logs.makeflow.MakeflowLogsParser`
     - Execution directory + Resource Monitor logs
   * - Nextflow
     - :class:`~wfcommons.wfinstances.logs.nextflow.NextflowLogsParser`
     - Execution directory with a trace file
   * - Pegasus
     - :class:`~wfcommons.wfinstances.logs.pegasus.PegasusLogsParser`
     - Submit directory
   * - RO-Crate
     - :class:`~wfcommons.wfinstances.logs.ro_crate.ROCrateLogsParser`
     - RO-Crate directory
   * - Snakemake
     - :class:`~wfcommons.wfinstances.logs.snakemake.SnakemakeLogsParser`
     - Execution directory + snkmt SQLite database
   * - StreamFlow
     - :class:`~wfcommons.wfinstances.logs.streamflow.StreamflowLogsParser`
     - RO-Crate directory
   * - TaskVine
     - :class:`~wfcommons.wfinstances.logs.taskvine.TaskVineLogsParser`
     - :code:`vine-logs` directory

Below we give basic examples for using the log parsers. Review the API
documentation of each log parser for additional parameters that configure its
behavior.

Makeflow
++++++++

`Makeflow <http://ccl.cse.nd.edu/software/makeflow/>`_ is a workflow system
for executing large complex workflows on clusters, clouds, and grids. The
Makeflow language is similar to traditional "Make", and a workflow can be just
a few commands chained together or a complex application consisting of
thousands of tasks. The following example shows the analysis of Makeflow
execution logs, stored in a local folder (:code:`execution_dir`), using the
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
    workflow.write_json(pathlib.Path('./makeflow-workflow.json'))

.. note::
    The :class:`~wfcommons.wfinstances.logs.makeflow.MakeflowLogsParser` class
    requires that Makeflow workflows run with the
    `Resource Monitor <https://cctools.readthedocs.io/en/latest/resource_monitor/>`_
    tool (e.g., execute the workflow using :code:`--monitor=logs`).

Nextflow
++++++++

`Nextflow <https://nextflow.io>`_ is a reactive workflow framework and a
programming DSL that eases the writing of data-intensive computational
pipelines. The following example shows the analysis of Nextflow execution
logs, stored in a local folder (:code:`execution_dir`), using the
:class:`~wfcommons.wfinstances.logs.nextflow.NextflowLogsParser` class: ::

    import pathlib
    from wfcommons.wfinstances import NextflowLogsParser

    # creating the parser for the Nextflow workflow
    execution_dir = pathlib.Path('/path/to/nextflow/execution/dir/')
    parser = NextflowLogsParser(execution_dir=execution_dir,
                                nextflow_version='24.10.0')

    # generating the workflow instance object
    workflow = parser.build_workflow('nextflow-workflow-test')

    # writing the workflow instance to a JSON file
    workflow.write_json(pathlib.Path('./nextflow-workflow.json'))

.. note::
    The :class:`~wfcommons.wfinstances.logs.nextflow.NextflowLogsParser` class
    expects the execution directory to contain a Nextflow trace file (by
    default, any file matching :code:`*trace*.txt`, as produced by running
    Nextflow with the :code:`-with-trace` option or with tracing enabled in
    the configuration).

Pegasus
+++++++

`Pegasus <http://pegasus.isi.edu>`_ is used in production to execute workflows
for dozens of high-profile applications in a wide range of scientific domains.
It provides the necessary abstractions for scientists to create workflows and
allows for transparent execution of these workflows on a range of compute
platforms, with execution managed by HTCondor DAGMan. The following example
shows the analysis of Pegasus execution logs, stored in a local submit
directory, using the
:class:`~wfcommons.wfinstances.logs.pegasus.PegasusLogsParser` class: ::

    import pathlib
    from wfcommons.wfinstances import PegasusLogsParser

    # creating the parser for the Pegasus workflow
    submit_dir = pathlib.Path('/path/to/pegasus/submit/dir/seismology/chameleon-100p-001/')
    parser = PegasusLogsParser(submit_dir=submit_dir)

    # generating the workflow instance object
    workflow = parser.build_workflow('pegasus-workflow-test')

    # writing the workflow instance to a JSON file
    workflow.write_json(pathlib.Path('./pegasus-workflow.json'))

RO-Crate
++++++++

`RO-Crate <https://www.researchobject.org/ro-crate/>`_ is a format for
packaging research data so as to promote open and reproducible science. The
RO-Crate logs parser processes RO-Crate artifacts created by executing
workflows (e.g., with the `StreamFlow <https://streamflow.di.unito.it/>`_
workflow management system). It takes as input the path to the RO-Crate
directory, which contains the :code:`ro-crate-metadata.json` file, and the
name of the workflow system that produced the crate: ::

    import pathlib
    from wfcommons.wfinstances import ROCrateLogsParser

    crate_dir = pathlib.Path('/path/to/ro-crate/dir/')
    parser = ROCrateLogsParser(crate_dir=crate_dir, wms_name='streamflow')
    workflow = parser.build_workflow('ro-crate-workflow-test')
    workflow.write_json(pathlib.Path('./ro-crate-workflow.json'))

Snakemake
+++++++++

`Snakemake <https://snakemake.readthedocs.io>`_ is a popular and easy-to-use
workflow system. The Snakemake logs parser processes execution logs generated
by the Snakemake
`snkmt plugin <https://snakemake.github.io/snakemake-plugin-catalog/plugins/logger/snkmt.html>`_.
It takes as input the path of the directory where all workflow data files
reside (input and output of workflow tasks), and the path of the SQLite
database created by the snkmt plugin: ::

    import pathlib
    from wfcommons.wfinstances import SnakemakeLogsParser

    execution_dir = pathlib.Path('/path/to/snakemake/execution/dir/')
    snkmt_db = pathlib.Path('/path/to/snkmt/sqlite/database/file')
    parser = SnakemakeLogsParser(execution_dir=execution_dir, snkmt_db=snkmt_db)
    workflow = parser.build_workflow('snakemake-workflow-test')
    workflow.write_json(pathlib.Path('./snakemake-workflow.json'))

StreamFlow
++++++++++

`StreamFlow <https://streamflow.di.unito.it/>`_ is a container-native workflow
management system. The StreamFlow logs parser processes the RO-Crate archive
produced by a StreamFlow execution: ::

    import pathlib
    from wfcommons.wfinstances import StreamflowLogsParser

    crate_dir = pathlib.Path('/path/to/streamflow/ro-crate/dir/')
    parser = StreamflowLogsParser(crate_dir=crate_dir,
                                  streamflow_version='0.2.0')
    workflow = parser.build_workflow('streamflow-workflow-test')
    workflow.write_json(pathlib.Path('./streamflow-workflow.json'))

TaskVine
++++++++

`TaskVine <https://ccl.cse.nd.edu/software/taskvine/>`_ is a task scheduler
for data-intensive dynamic workflows. The TaskVine logs parser translates the
logs found in a TaskVine :code:`vine-logs` directory into workflow instances
compatible with :ref:`json-format-label`: ::

    import pathlib
    from wfcommons.wfinstances import TaskVineLogsParser

    vine_logs_dir = pathlib.Path('/path/to/taskvine/vine-logs/dir/')
    parser = TaskVineLogsParser(vine_logs_dir=vine_logs_dir)
    workflow = parser.build_workflow('taskvine-workflow-test')
    workflow.write_json(pathlib.Path('./taskvine-workflow.json'))

The Instance Analyzer
---------------------

The :class:`~wfcommons.wfinstances.instance_analyzer.InstanceAnalyzer` class
provides a number of tools for analyzing collections of workflow execution
instances. Its goal is to analyze one or multiple workflow execution instances
and build summaries of the analysis per workflow task type prefix. These
summaries are what allows WfChef recipes — and thus WfGen's synthetic
workflows — to be *realistic* rather than merely random.

.. warning::

    Although any workflow execution instance represented as an
    :class:`~wfcommons.wfinstances.instance.Instance` object (i.e., compatible
    with :ref:`json-format-label`) can be appended to the
    :class:`~wfcommons.wfinstances.instance_analyzer.InstanceAnalyzer`, we
    strongly recommend that only instances of a single workflow application
    type be appended to an analyzer object. You may, though, create several
    analyzer objects per workflow application.

The :meth:`~wfcommons.wfinstances.instance_analyzer.InstanceAnalyzer.append_instance`
method allows you to include instances for analysis. The
:meth:`~wfcommons.wfinstances.instance_analyzer.InstanceAnalyzer.build_summary`
method processes all appended instances. It applies probability distribution
fitting to a series of data to find the *best* probability distribution
representing the analyzed data (i.e., the one that minimizes the mean square
error). The method returns a summary of the analysis of instances in the form
of a Python dictionary object, in which keys are task prefixes (provided when
invoking the method) and values describe the best probability distribution fit
for tasks' runtime and input and output data file sizes. The code excerpt
below shows an example of an analysis summary showing the best-fit probability
distribution for the runtime of the :code:`individuals` tasks (1000Genome
workflow): ::

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

Workflow analysis summaries are used by WfChef to develop
:ref:`workflow-recipe-label`, which in turn are used to :ref:`generate
realistic synthetic workflow instances <generating-workflows-label>`.

Probability distribution fits can also be plotted by using the
:meth:`~wfcommons.wfinstances.instance_analyzer.InstanceAnalyzer.generate_fit_plots` or
:meth:`~wfcommons.wfinstances.instance_analyzer.InstanceAnalyzer.generate_all_fit_plots`
methods — plots will be saved as :code:`png` files.

Example: analyzing a set of Seismology instances
++++++++++++++++++++++++++++++++++++++++++++++++

The following example shows the analysis of a set of instances, stored in a
local folder, of a Seismology workflow. In this example, we seek the best
probability distribution fitting for task *prefixes* of the Seismology
workflow (:code:`sG1IterDecon` and :code:`wrapper_siftSTFByMisfit`), and
generate all fit plots (runtime, and input and output files) into the
:code:`fits` folder using :code:`seismology` as a prefix for each generated
plot: ::

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

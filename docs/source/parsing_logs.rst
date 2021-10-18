.. _logs-label:

Parsing Workflow Execution Logs
===============================

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

Makeflow
--------

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
--------

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

Pegasus WMS
-----------

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

.. warning::
    By default, the :class:`~wfcommons.wfinstances.logs.pegasus.PegasusLogsParser`
    class assumes that the submit dir is from a Pegasus execution with **version 5.0**
    or later. To enable parsing of Pegasus execution logs from version 4.9 or earlier,
    the option :code:`legacy=True` should be used.

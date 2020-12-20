.. _logs-label:

Parsing Workflow Execution Logs
===============================

The most common way for obtaining traces from actual workflow executions is to parse
execution logs. As part of the WorkflowHub project, we are constantly developing
parsers for commonly used workflow management systems.

Each parser class is derived from the :class:`~workflowhub.trace.logs.abstract_logs_parser.LogsParser`
class. Thus, each parser provides a
:meth:`~workflowhub.trace.logs.abstract_logs_parser.LogsParser.build_workflow`
method.

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
logs provide time-stamped event traces from these executions. The following example shows
the analysis of Pegasus execution logs, stored in a local folder (submit dir), for a
workflow execution using the :class:`~workflowhub.trace.logs.pegasus.PegasusLogsParser`
class: ::

    from workflowhub.trace import PegasusLogsParser

    # creating the parser for the Pegasus workflow
    parser = PegasusLogsParser(submit_dir='/path/to/pegasus/submit/dir/seismology/chameleon-100p-001/')

    # generating the workflow trace object
    workflow = parser.build_workflow('workflow-test')

    # writing the workflow trace to a JSON file
    workflow.write_json('workflow.json')

.. note::

    By default, the :class:`~workflowhub.trace.logs.pegasus.PegasusLogsParser`
    class assumes that the submit dir is from a Pegasus execution with version 5.0 or later.
    To enable parsing of Pegasus execution logs from version 4.9 or earlier, the option
    :code:`legacy=True` should be used.

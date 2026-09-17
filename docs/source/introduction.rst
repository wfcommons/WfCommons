The WfCommons Project
=====================

The `WfCommons project <https://wfcommons.org>`_ is an open-source framework
for enabling scientific workflow research and development. It provides
foundational tools for analyzing workflow execution instances and for
generating synthetic, yet realistic, workflow instances and benchmarks. These
can be used to develop, evaluate, and verify new techniques, algorithms, and
systems that target the efficient and robust execution of ever-larger
workflows on increasingly complex distributed infrastructures.

The figure below shows an overview of the research and development life cycle
that integrates the major components of WfCommons: (i) workflow execution
instances (**WfInstances**), (ii) workflow recipes (**WfChef**), (iii) workflow
generation (**WfGen**), (iv) workflow benchmarks (**WfBench**), and (v)
workflow simulation (**WfSim**).

.. figure::  images/wfcommons.png
   :scale:   70 %
   :align:   center

   The WfCommons conceptual architecture.

The components
--------------

:ref:`WfInstances <instances-label>` — *learn from real executions.*
The WfInstances component provides a collection and curation of open-access
production workflow instances from various scientific applications, all made
available in a common format (:ref:`json-format-label`). A workflow instance
is built from the logs of an actual execution of a scientific workflow on a
distributed platform (e.g., clouds, grids, clusters, HPC) using a workflow
system. This Python package ships log parsers for widely used workflow
systems, and we maintain a `GitHub repository of workflow execution instances
<https://github.com/wfcommons/WfInstances>`_ ready to use.

:ref:`WfChef <generating-workflows-recipe-label>` — *capture what makes a
workflow realistic.* The WfChef component automates the construction of
synthetic workflow generators (**recipes**) for any workflow application. Given
a set of real instances in :ref:`json-format-label`, WfChef discovers the
recurring subgraph patterns of the application and derives statistical models
of task runtimes and data sizes. The result is a self-contained recipe package
— no manual modeling is required.

:ref:`WfGen <generating-workflows-label>` — *generate workflows at any
scale.* The WfGen component generates realistic synthetic workflow instances.
It takes a recipe produced by WfChef and a desired number of tasks, and
produces randomized workflow instances that preserve the structure and
performance distributions of the original application — at scales far beyond
(or below) the available real instances.

:ref:`WfBench <generating-workflow-benchmarks-label>` — *run realistic
experiments on real systems.* The WfBench component generates workflow
benchmark specifications with tunable performance characteristics (CPU,
memory, and I/O usage) and realistic dependency structures. Translators turn
these specifications into runnable benchmarks for production workflow systems
(e.g., Airflow, Dask, Nextflow, Parsl, Pegasus, Snakemake, TaskVine).

**WfSim** — *simulate what you cannot run.* The WfCommons project fosters the
use of simulation for the development, evaluation, and verification of
scheduling and resource provisioning algorithms, and for the evaluation of
current and emerging computing platforms. We do not develop simulators as part
of the WfCommons project; instead, the WfSim component catalogs open-source
workflow management system simulators that support :ref:`json-format-label`,
such as `WRENCH <https://wrench-project.org>`_. Because real and synthetic
instances share the same format, simulators can use them interchangeably.

.. _json-format-label:

WfFormat
--------

All WfCommons components communicate through **WfFormat**, a common JSON
format for representing workflow execution instances and generated synthetic
workflow instances. A WfFormat instance describes the workflow task graph
(tasks, files, and dependencies), per-task performance measurements (runtime,
input/output data sizes, memory, energy, etc.), and information about the
machines on which the workflow was executed.

WfFormat uses a JSON schema available in the
`WfFormat Schema GitHub <https://github.com/wfcommons/WfFormat>`_ repository.
The current version of the WfCommons Python package uses schema version
:code:`1.6`. The schema repository provides a detailed explanation of WfFormat
(including required fields) and a validator script for verifying the
compatibility of instances.

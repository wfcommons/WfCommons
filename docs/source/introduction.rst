The WfCommons Project
=======================

The `WfCommons project <https://wfcommons.org>`_ is an open source framework
for enabling scientific workflow research and development by providing foundational
tools for analyzing workflow execution instances, and generating synthetic, yet
realistic, workflow instances that can be used to develop new techniques, algorithms
and systems that can overcome the challenges of efficient and robust execution of
ever larger workflows on increasingly complex distributed infrastructures.

The figure below shows an overview of the research and development life cycle that
integrates the four major components WfCommons: (i) workflow execution instances
(**WfInstances**), (ii) workflow recipes (**WfChef**), (iii) workflow generator
(**WfGen**), and (iv) workflow simulator (**WfSim**).

.. figure::  images/wfcommons.png
   :scale:   70 %
   :align:   center

   The WfCommons conceptual architecture.

:ref:`WfInstances <instances-label>`.
The WfInstances component provides a collection and curation of open-access
production workflow instances from various scientific applications, all made
available using a common format (i.e., :ref:`json-format-label`).
A workflow instance is built based on logs of an actual execution of a scientific
workflow on a distributed platform (e.g., clouds, grids, clusters) using a
workflow system. We keep a `GitHub repository of workflow execution instances
<https://github.com/wfcommons/WfInstances>`_.

:ref:`WfChef <generating-workflows-recipe-label>`.
The WfChef component automates the construction of synthetic workflow generators
(recipes) for any given workflow application. The input to this component is a set
of real workflow instances described in the :ref:`json-format-label` (e.g.,
instances available in WfInstances).

:ref:`WfGen <generating-workflows-label>`.
The WfGen component targets the generation of realistic synthetic workflow instances.
WfGen takes as input a workflow recipe produced by WfChef for a particular application
and a desired number of tasks. WfGen then automatically generates synthetic, yet
realistic, randomized workflow instances with (approximately) the desired number of
tasks.

:ref:`WfBench <generating-workflow-benchmarks-label>`.
The WfBench component is a generator of realistic workflow benchmark specifications 
that can be translated into benchmark code to be executed with current workflow 
systems. it generates workflow tasks with arbitrary performance characteristics (CPU,
memory, and I/O usage) and with realistic task dependency structures based on those 
seen in production workflows.

**WfSim.**
The WfCommons project fosters the use of simulation for the development, evaluation,
and verification of scheduling and resource provisioning algorithms (e.g.,
multi-objective function optimization, etc.), evaluation of current and emerging
computing platforms (e.g., clouds, IoT, extreme scale, etc.), among others.
We do not develop simulators as part of the WfCommons project. Instead, the WfSim
component catalogs open source WMS simulators that provide support for
:ref:`json-format-label`.

.. _json-format-label:

WfFormat
--------

The WfCommons project uses a common format for representing workflow execution
instances and generated synthetic workflows instances. Workflow simulators and
simulation frameworks that support WfFormat can then use both types of instances
interchangeably. WfFormat uses a JSON specification available in the
`WfFormat Schema GitHub <https://github.com/wfcommons/WfFormat>`_
repository. The current version of the WfCommons Python package uses the schema
version :code:`1.5`. The schema GitHub repository provides detailed explanation
of WfFormat (including required fields), and also a validator script for verifying
the compatibility of instances.

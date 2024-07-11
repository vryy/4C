.. 4C-documentation documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

========
|FOURC|
========

Mission statement
=================

|FOURC| is a parallel multi-physics research code
to analyze and solve a plethora of physical problems
described by ordinary or partial differential equations.
Its development is driven by challenging research questions and real-world problems,
for which existing tools do not suffice,
either due to the lack of capabilities or due to falling short of accuracy or performance.

|FOURC| not only provides ready-to-use simulation capabilities for a variety of physical models,
including single fields such as solids and structures, fluids, or scalar transport,
and multi-physics coupling and interactions between several fields,
but also a modular software environment for research in mathematical modeling and numerical methods.
Pre- and post-processing tools facilitate the use of |FOURC| within streamlined application workflows in science and engineering.
For spatial discretization, |FOURC| mostly relies on finite element methods (FEM, CutFEM).
It leverages the `Trilinos project <https://trilinos.github.io>`_ for sparse linear algebra, nonlinear solvers, and linear solvers and preconditioners
to be executed on MPI-parallel computing clusters.
Through its comprehensive set of physics modules available to all users without coding effort,
|FOURC| facilitates the advancement of research in all areas of science, engineering, and biomedicine.

Content
=======

This guide to |FOURC| is structured as follows:

:ref:`Part I — About 4C<about>`
   Learn about the capabilities and history of |FOURC|.

:ref:`Part II — Tutorials<tutorials>`
   A series of beginner-level tutorials showcases the setup procedure for specific application scenarios.

:ref:`Part III - Analysis guide<analysisguide>`
   Detailed explanations on the whole tool chain from model generation (pre-processing)
   over running a simulation to the evaluation of results (post-processing) offers deep insight into using |FOURC|
   for advanced simulation scenarios.
   This guide includes background information and detailed descriptions
   for the specification of elements, boundary conditions, constitutive laws
   as well as options for linear and nonlinear solvers.

:ref:`Part IV - Developer guide<developerguide>`
   This guide gets you started on actively developing and contributing to |FOURC|.
   It covers the build process, our CI/CD testing infrastructure,
   coding guidelines, and useful tools for the daily development of |FOURC|.

:ref:`Part V - Input Parameter Reference<inputparameterreference>`
   A comprehensive list of all input parameters, elements, materials, and boundary conditions
   with short descriptions for each option

:ref:`Part VI - Tools and Scripts<toolsAndScripts>`
   A collection of useful scripts for working with |FOURC|

:ref:`Appendix<appendix>`
   Information on contributing to this documentation as well as seleceted topics of interest

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents:

   about
   tutorials
   analysisguide
   developmentguide
   parameterreference
   tools
   appendix

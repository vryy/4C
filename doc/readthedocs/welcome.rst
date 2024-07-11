.. _WelcomeTo4C:

===================
Welcome to |FOURC|
===================

Purpose and objectives
---------------------

|FOURC| is a parallel multi-physics research code
to analyze and solve a plethora of physical problems
described by ordinary or partial differential equations.
Its development is driven by challenging research questions and real-world problems,
for which existing tools do not suffice, either due to the lack of capabilities or due to falling short of accuracy or performance.
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

History
-------

The development of |FOURC| arose from the need to tackle challenging research questions
in and with numerical methods for ordinary and partial differential equations
and to advance complex models for real-world applications,
all based on a proper theoretical foundation and with verified and state-of-the-art methods and software implementations.
Since suitable tools are often not available either in commercial or in other (academic) research codes,
we want to close this gap by developing |FOURC|, a comprehensive multi-physics simulation framework.
|FOURC| originated in the 2000s at the `Institute for Computational Mechanics <https://www.epc.ed.tum.de/lnm/home/>`_ of the Technical University of Munich [#f1]_.

.. _items-to-be-added:

List of items to be added
--------------------------

This list shall help the developers to find sections,
where action is desperately needed to make this documentation a useful piece of work.
Unfortunately, also this list is still under development, as all the rest of the site.
Note that some basic information on the format of restructuredText can be found :ref:`here<writingdocumentation>`.

The following list only contains a few parts and sections, which should be filled as soon as possible.

- Add some information to the tutorials, where content is missing

  - :ref:`3dsolidtutorial`
  - :ref:`multiphysicstutorial`

- In the analysis guide, which is the main part of the documentation where information about the procedures and features of the program are contained,
  several sections do not contain all necessary information.

  - In the :ref:`Problem types<problemtypes>` section, all single field and multi field problems are mentioned,
    but many of the subsections are incomplete or even empty. Actually, only very few information should be given here in the beginning,
    but we definitely need the purpose of the problem type, the element types used and the degrees of freedom to be solved.

  - The whole :ref:`Elements<elements>` section is lacking useful information.
    At this point, we only have the materials that may be used with each element type listed.
    Any ideas on useful information is highly appreciated.

  - There is some general information on solid material models in the :ref:`Materials<materials>` section.
    Information on other material models is needed.
    In the more distant future it would also be good to have the underlying equations for some material models in this section.

  - There are two sections that are directly taken from the ``global_report.pdf``:
    :ref:`old_Fluid-structure Interaction<fluid_structure_interaction>` and :ref:`old_Free Surface Flow<free_surface_flow>`.
    Somebody should check, which parts of these sections can still be used, ond which should be removed.

.. rubric:: Footnotes

.. [#f1] A predecessor has been developed under the name BACI, which has been built upon the extensive experience of individual developers of the software ccarat of the Institute for Structural Mechanics of the University of Stuttgart.
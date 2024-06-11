.. _WelcomeTo4C:

===================
Welcome to |FOURC|
===================

History
-------

|FOURC| is a C++ research code that has its roots in the by now abandoned C
code ``ccarat`` and its successor, written in modern C++, BACI (the Bavarian Advanced Computation Initiative).
``ccarat`` was started by Michael Gee in 2000 and has
since grown in several directions. The aim back then was to provide a
framework for parallel multifield finite element simulation, including
moderate parallelization, not massive parallel execution with hundreds of processes.
That is why there is a strong emphasis on data structures in ``ccarat``
but still a lot of data redundancies.
The original ``ccarat`` is a large C project that uses
external linear solver packages (both parallel and serial).

In 2006 it was felt that ``ccarat`` needs a major face lift.
And again it is due to Michael Gee to introduce the trilinos libraries into ``ccarat`` and redesign the whole package on top of that.
This is how the predecessor of |FOURC|, BACI, emerged in spring 2007.
The ``ccarat`` code was replaced step by step (even if some of these steps were XXL sized) until close to no original code remained.
Eventually, after lots and lots of other improvements and extensions, the code was renamed as |FOURC| in spring 2024, and is to become open source soon.


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


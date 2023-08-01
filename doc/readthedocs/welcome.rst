Welcome to BACI
===============

History
-------

BACI (the Bavarian Advanced Computation Initiative) is a C++ research code that has its roots in the by now abandoned C
code ``ccarat``. ``ccarat`` was started by Michael Gee in 2000 and has
since grown in several directions. The aim back then was to provide a
framework for parallel multifield finite element simulation, including 
moderate parallelization, not massive parallel execution with hundreds of processes. 
That is why there is a strong emphasis on data structures in ``ccarat`` 
but still a lot of data redundancies. 
The original ``ccarat`` is a large C project that uses
external linear solver packages (both parallel and serial).

In 2006 it was felt that ``ccarat`` needs a major face lift. And again
it is due to Michael Gee to introduce the trilinos libraries into
``ccarat`` and redesign the whole package on top of that. This is how
BACI emerged in spring 2007. The ``ccarat`` code was replaced step by step (even if
some of these steps were XXL sized) until close to no original code
remained.


Overview
--------

The following Guide to ``BACI`` deals with

Part I — Tutorials
   A number of tutotials are available, which also serve as framework tests during the gitlab CI.

Part II - The BACI Workflow
   Explanations on the whole tool chain from model generation (pre processing) 
   to results evaluation (post processing)

Part III - Analysis guide
   Theoretical background and detailed description on elements, data types, ...

Part IV - Developer guide
   Contains useful information for developers, about adding Doxygen information, C++ basics, etc.

Part V - Input Parameter Reference
   Input parameter (including elements, materials, boundary conditions) with short description

Appendix
   Setup Guide, some additional theory information, 
   and a section about how to add information to the ReadTheDocs documentation.
   Also, you'll find the index and the references therein.
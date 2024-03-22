.. 4C-documentation documentation master file, created by
   sphinx-quickstart on Wed Nov 23 14:27:39 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


========
|FOURC|
========

.. note::

    This is a draft of the |FOURC| documentation, and the content is rather incomplete.
    The faster additional content is added,
    the sooner the documentation will be useful for users to understand the various types of simulations and their application.
    Some documentation parts are highly relevant for many users; those need an update with highest priority.
    A list of items that should be added asap, is given :ref:`here<items-to-be-added>`.

Content
=======

The following Guide to |FOURC| deals with

:ref:`Part I - Welcome <WelcomeTo4C>`
   This section outlines the history of the code.
   Also you'll find a list of documentation sections that needs some update.

:ref:`Part II — Tutorials<tutorials>`
   A number of tutorials are available, which also serve as framework tests during the gitlab CI.

:ref:`Part III - Analysis guide<analysisguide>`
   Explanations on the whole tool chain from model generation (pre processing)
   to results evaluation (post processing)
   Theoretical background and detailed description on elements, data types, ...

:ref:`Part IV - Developer guide<developerguide>`
   Contains useful information for developers, about adding Doxygen information, C++ basics, etc.

:ref:`Part V - Input Parameter Reference<inputparameterreference>`
   Input parameter (including elements, materials, boundary conditions) with short description

:ref:`Part VI - Tools and Scripts<toolsAndScripts>`

:ref:`Appendix<appendix>`
   Some additional theory information,
   and a section about how to add information to the ReadTheDocs documentation.
   Also, you'll find the index and the references therein.


.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents:

   welcome
   tutorials
   analysisguide
   developmentguide
   parameterreference
   tools
   appendix




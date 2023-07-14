
.. _SetupGuidetoBACI:

Setup Guide to BACI
===================

Here you'll find some useful notes on setting up and running BACI, 
the multi purpose multi physics multi-modelling code developed at the Chair of Computational Mechanics, 
the IMCS of the Hochschule der Bunderwehr München, and 
the Institute of Material Systems Modelling of the Helmholtz Zentrum Hereon. 
For information about solving specific types of problems using BACI see the respective sections in the analysis guide chapter. 


Directory structure
~~~~~~~~~~~~~~~~~~~

The ``BACI`` code comes with documentation, example input files and
support scripts. The important subdirectories are the following:


:``src``: contains the real ``BACI`` code in several subdirectories

:``Input``:   contains various valid and running ``*``\ ``.dat`` files, ``BACI``
    input files that are used for (automatic) testing.

:``utilities``:  contains configuration scripts needed to setup a ``BACI`` Makefile,
    
:``tests``:   contains tests that are not directly related to baci, but rather to the input
    and output, balso also comprising the whole tool chain (framework tests).

:``unittests``:  contains the source files for stand alone tests aka unittests.

:``doc``:   contains this documentation and also Doxygen

:``presets, buildconfig``:   contains configuration files (platform specifications) used for cmake 
    to setup a ``BACI`` Makefile

:``docker``: contains the setup file for a docker container for running baci inside of it.



Setup and Run
-------------

``BACI`` is developed and used on Linux. Other Unixes work as well.
Windows versions might be created using cygwin or mingw, but this will
require some small modifications and is not covered here.
People around have already ported the code to windows, but this needs quite a number of changes,
so it is not advised to try on your own, and the topic will not be covered here.




``BACI`` is a non-interactive shell application that reads an input
file and creates a bunch of files in return. To build and run it you
will need a basic understanding of Linux and the Linux shell. You will
also want to choose your favorite text editoror integrated development environment (IDE).

You'll find more information about the ``BACI`` installation in the ``readme.md`` and the
``contribute.md`` files located in the BACI root directory.




Running examples
~~~~~~~~~~~~~~~~

In ``Input`` there are test examples; all necessary “packages” must have
been activated in the defines-file that was used to configure the
``BACI`` at hand. For example,

::

   ./baci-release Input/f2_drivencavity20x20_drt.dat xxx

runs the 2d fluid driven cavity example and writes the output to files
beginning with ``xxx``. 
You can also run the code in parallel with the mpirun
command like this:

::

   mpirun -np 1 ./baci-release Input/f2_drivencavity20x20_drt.dat xxx

Testing
-------

The BACI code is tested by all tests included in the Input directory, whenever modifications enter the code.
However, after changing something in the code locally, it is a good idea to run the test environment yourself.
This is done by

::
    ctest [-L minimal] [-R <regular expression>] [-I Start,End,Stide] [-j <ncpus>]

Without any options, all the more than 2000 tests are run, which will probably take hours.
If ctest is run with the option ``-L minimal``, only some 20 tests are run. 
However, it might be that these tests don't check the part of the code that you modified.
With the options ``-R`` and ``-I`` one can specify granularly which tests should run.



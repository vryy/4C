Directory structure
--------------------
The |FOURC| code comes with documentation, example input files and
support scripts. The important subdirectories are the following:


:``src``: contains the real |FOURC| code in several subdirectories

:``tests/input_files``:   contains various valid and running ``*``\ ``.dat`` files,
    |FOURC| input files that are used for (automatic) testing.

:``utilities``:  contains configuration scripts needed to setup a |FOURC| Makefile,

:``tests``:   contains tests that are not directly related to |FOURC|, but rather to the input
    and output, balso also comprising the whole tool chain (framework tests).

:``unittests``:  contains the source files for stand alone tests aka unittests.

:``doc``:   contains this documentation and also Doxygen

:``presets, buildconfig``:   contains configuration files (platform specifications) used for cmake
    to setup a |FOURC| Makefile

:``docker``: contains the setup file for a docker container for running |FOURC| inside of it.


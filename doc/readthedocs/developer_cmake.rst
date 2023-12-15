.. _cmakepresets:

cmake presets
--------------

CMake presets have been introduced with the Merge Request `!1392 <https://gitlab.lrz.de/baci/baci/-/merge_requests/1392>`_
and can be used to configure and manage different configurations of baci in a very convenient way.
This small article will go through a few of them. The experts should also read the
`official CMake presets documentation <https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html>`_ to get all convenient tricks.
If you find a nice one, feel free to tell your colleagues (via issue, Slack or TGM).

CMake presets are available since cmake 3.19 and have been improved in the following versions.
We are currently using CMake preset version 5, which requires at least cmake 3.24.
Make sure to have a sufficient cmake-version in your path or use the provided onces on the institute's server.

Configuration from a terminal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose you are in the build directory where you want to build baci. You can configure a release build of baci with::

    cmake --preset=lnm_workstation ../path/to/source

**Hint:** The global configurations of cmake are stored in the ``CMakeCache.txt`` within the build folder
and it's sometimes helpful to remove it with ``rm CMakeCache.txt``, before configuring the project new.
Other configurations that may be of interest are those created for the two main host institutions:

**lnm:**

::

    lnm_workstation_debug
    lnm_workstation_relwithdebinfo
    lnm_workstation_coverage
    lnm_bruteforce
    lnm_deep
    lnm_schmarrn

**imcs:**

::

    imcs_workstation
    imcs_workstation_debug
    imcs_workstation_mirco
    imcs_charon

Configuration from the IDE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most ides already support cmake presets (at least VS Code and CLion > 2023.1).

Here is a screenshot taken from VS Code:


.. figure:: figures/vs-code-cmake-preset.png
   :alt: CMake preset selection within VS Code
   :width: 100%

Hints for VS Code: You need to install the extensions "CMake Tools" from Microsoft.
If you don't have a recent cmake-version in your path, you can set the path to the institute cmake-version in the VS Code settings.


Defining your own CMake presets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CMake presets allow you to also create your own configuration.
You need to put a ``CMakeUserPresets.json``-file (**important:** ``User``) in the source-directory of baci.
This file will not be part of the repository (it is listed in ``.gitignore``). There you can define your own configurations.
There you can define the binary directory, so you don't need to go to your binary directory in order to configure baci
and your ide also knows then where to put the build folder.
You can define as many configurations as you need. Note, you can inherit from other configurations by using the keyword ``inherits``.

An example ``CMakeUserPresets.json`` could look like this:

::

    {
      "version": 5,
      "configurePresets": [
        {
          "name": "release",
          "binaryDir": "../build_parallel",
          "generator": "Ninja",
          "inherits": [
            "lnm_workstation"
          ]
        },
        {
          "name": "debug",
          "binaryDir": "../build_debug",
          "generator": "Ninja",
          "inherits": [
            "lnm_workstation_debug"
          ]
        }
      ]
    }


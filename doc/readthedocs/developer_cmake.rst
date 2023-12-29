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

There is a number of available preset files, which can be retrieved by ``cmake <baci_home> --list-presets``.
This is the current output of this command:

.. literalinclude:: baci-cmake-presets.txt

In general, it is highly recommended to create your own preset, see bélow.

Defining your own CMake presets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CMake presets allow you to also create your own configuration.
You need to put a ``CMakeUserPresets.json``-file (**important:** ``User``) in the source-directory of baci.
This file will not be part of the repository (it is listed in ``.gitignore``).
There you can define your own configurations.
Particularly, you may define the binary directory,
so you don't need to go to your binary directory in order to configure BACI.
CMake presets integrate well with recent releases of IDEs.

You can define as many configurations as you need.
Note that you can inherit from other configurations by using the keyword ``inherits``.

Such a local preset could look like this::

    {
      "version": 5,
      "configurePresets": [
        {
          "name": "myworkstation",
          "displayName": "Release build for my workstation",
          "binaryDir": "/home/baci/baci_release",
          "generator": "Ninja",
          "inherits": [
            "lnm_workstation"
          ],
          "cacheVariables": {
            "CMAKE_CXX_COMPILER": "/usr/bin/mpic++",
            "CMAKE_CXX_COMPILER_LAUNCHER": "ccache",
            "BACI_WITH_GOOGLETEST": "OFF",
            "BACI_BUILD_READTHEDOCS": "ON",
            "BACI_SPHINX_THEME": "sphinx_rtd_theme",
            "BACI_BUILD_DOXYGEN": "ON",
          }
        }
      ]
    }


Configuration from the IDE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Our recommended IDEs (VS Code and CLion) already support cmake presets natively.

Here is a screenshot taken from VS Code:


.. figure:: figures/vs-code-cmake-preset.png
   :alt: CMake preset selection within VS Code
   :width: 100%

Hints for VS Code: You need to install the extensions "CMake Tools" from Microsoft.


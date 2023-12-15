.. _ccache:

Speed up recompilation using ``ccache``
----------------------------------------

Working with many different branches in the same local repository is a key feature of Git providing an efficient work flow for code development.

After switching to another branch in the repository, the code must be **recompiled**.
In case of only minimal differences from one to another branch, running ``make`` recompiles only the parts that are actually different and finally relinks them.
However, if there are differences in central header files, `make` will recompile (nearly) everything.

In both ways, ``ccache`` can speed up recompilation time (refer to `ccache <https://ccache.samba.org/manual/latest.html>`__):

.. note::

    ``ccache`` is a compiler cache. It speeds up recompilation by caching the result of previous compilations and detecting when the same compilation is being done again.

Setting up ccache
~~~~~~~~~~~~~~~~~~~

Setting up ``ccache`` depends on your operating system and local environment.
If ``ccache`` is installed and you are using CMake 3.4 or later, ``ccache`` will be set up automatically durinig the BACI configure process
(this can be deactivated with the ``--no-ccache`` configure option).
If your CMake version is outdated, the general procedure is described `here <https://ccache.samba.org/manual/latest.html>`__.
However, it is recommended to talk to your co-workers first to get help.

Useful ccache options
~~~~~~~~~~~~~~~~~~~~~~~~

The following ``ccache`` options are cited from the `official website <https://ccache.samba.org/manual/latest.html>`__:

Print an options summary page.

::

    -h, --help

Print the current statistics summary for the cache.

::

    -s, --show-stats


Clear the entire cache, removing all cached files, but keeping the configuration file.

::

    -C


Zero the cache statistics (but not the configuration options). 

::

    -z, --zero-stats


Enforce a complete recompilation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some cases (e.g., to check for compiler warnings) it is necessary to recompile the complete code:

::

    make clean
    ccache -C -z
    make



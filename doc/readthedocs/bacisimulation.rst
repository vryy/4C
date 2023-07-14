.. _bacisimulation:

BACI Simulation
================

Depending on the compilation (with or without debug information), the resulting executable for BACI will be
either ``baci-release`` or ``baci-debug``. In the following, the name will always be denoted as ``<baci-command>``.

The options when starting BACI from the command line are shown by an output of ``<baci-command> --help``:

.. literalinclude:: baci-help.txt

The output consists at least of one file named ``<output_name>.control``.
If further output has been requested, a number of further files are created.
The file names are all given in the \*.control file.

These results in the output files are by default stored in a proprietary format. 
However, output data can also be stored in a format that is readable by other programs. 
For further information, see :ref:`Postprocessing`.


Restarting an analysis
-----------------------

For restarting an analysis one has to provide restart information in the first simulation by including the parameter
``RESTARTEVRY  <numsteps>`` in the ``--STRUCTURAL DYNAMICS`` section,
so that the information is written every *numsteps* step.

In the second simulation, no additional parameters have be included. The information that it is a restart, is
given on the command line:

::

   <baci-command> <datfile> <outputfile> [restartfrom=<restart_filename>] restart=<step>

Here, one has to provide the step, at which the restart is started from the previous simulation.
If the parameter ``restartfrom`` is given, the initial configuration is read from this file, 
otherwise it is read from ``outfile``. In the latter case the filename of the new output is the same with an appended number, e.g., ``outfile-1``.
Note that the value for ``step`` must be given in the
file ``<outfile>.control`` in one of the step lines: ``step = <step>``.

.. note::

   - The parameters RESTART and RESTARTTIME in the PROBLEMTYP section 
     are not needed anymore, and will probably vanish soon.
   - The parameter MAXTIME indicates the maximum time of all simulations, 
     it is NOT a step time, so be aware that MAXTIME in the subsequent simulation 
     must be larger than in the first one.
   - The parameter TIMEINIT can be 0 also in the subsequent simulation, it is not used, 
     since the initial time is defined by the restart.

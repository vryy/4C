Run a simulation
================

``BACI`` is a non-interactive shell application that reads an input file and creates a bunch of files in return.
To build and run it you will need a basic understanding of Linux and the Linux shell.
You will also want to choose your favorite text editor or integrated development environment (IDE).


Running examples
----------------

In ``tests/input_files`` there are test examples; all necessary “packages” must have
been activated in the defines-file that was used to configure the
``BACI`` at hand. For example,

.. code-block:: bash

   ./baci-release tests/input_files/f2_drivencavity20x20_drt.dat xxx

runs the 2d fluid driven cavity example and writes the output to files
beginning with ``xxx``.
You can also run the code in parallel with the mpirun
command like this:

.. code-block:: bash

   mpirun -np 1 ./baci-release tests/input_files/f2_drivencavity20x20_drt.dat xxx


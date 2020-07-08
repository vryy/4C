# Utility functions for adapting IDEs.


import os


def get_mpi_include_path():
    """
    Get the openmpi include directory.
    The first one of the choices that exists on the system will be used.
    """

    possible_open_mpi_dirs = [
        "/usr/include/openmpi-x86_64",
        "/usr/include/openmpi/1.2.4-gcc",
        "/usr/include/openmpi"
        ]
    for path in possible_open_mpi_dirs:
        if os.path.isdir(path):
            return path
    return None
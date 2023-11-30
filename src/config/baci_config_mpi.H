/*---------------------------------------------------------------------*/
/*! \file
\brief Configuration of MPI
\level 0
*/
/*---------------------------------------------------------------------*/

#ifndef BACI_CONFIG_MPI_H
#define BACI_CONFIG_MPI_H

// Disable the (deprecated) C++ bindings of OpenMPI. We are not using them, and they throw warnings
// that are false positives.
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>  // IWYU pragma: export

#endif
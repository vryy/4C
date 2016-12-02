/*!---------------------------------------------------------------------
\file pss_cpp.h
\brief Proper access to legacy parser module written in C.

\maintainer Martin Kronbichler

\level 1

---------------------------------------------------------------------*/

#ifndef PSS_CPP_H
#define PSS_CPP_H

// this ensures we always obtain the C++ version of mpi and
// not its C version included by pss_table.h
#ifdef PARALLEL
  #include <mpi.h>
#endif

extern "C"
{
  #include "../pss_full/pss_table.h" // access to C methods
}

#endif

/*!
\file
\brief Proper access to legacy parser module written in C.

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bauer
            089 - 289-15252
</pre>

*/

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

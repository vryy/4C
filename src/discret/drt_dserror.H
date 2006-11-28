/*!----------------------------------------------------------------------
\file dserror.H
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifndef DSERROR_H
#define DSERROR_H

#ifdef PARALLEL
#include "mpi.h"
#endif

#ifdef __cplusplus

#include "cstdio"
#include "cstdarg"
void cpp_dsassert_func(char* file, int line, int test, char string[]);
void cpp_dserror_func(char *string, ...);
void cpp_dslatest(char* file, int line);
#define dsassert(test, string) cpp_dsassert_func(__FILE__, __LINE__, test, string)
#define dserror cpp_dslatest(__FILE__, __LINE__), cpp_dserror_func

#endif // #ifdef __cplusplus




#endif  // #ifndef DSERROR_H

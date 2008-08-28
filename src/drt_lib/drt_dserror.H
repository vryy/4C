/*!----------------------------------------------------------------------
\file drt_dserror.H
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

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

#ifdef __cplusplus

#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <string>

void cpp_dsassert_func(const std::string file, const int line, const bool test, const std::string string);
void cpp_dserror_func(const std::string string, ...);
void cpp_dslatest(const std::string file, const int line);
#define dsassert(test, string) cpp_dsassert_func(__FILE__, __LINE__, test, string)
#define dserror cpp_dslatest(__FILE__, __LINE__), cpp_dserror_func

#ifdef THROWELEMENTERRORS

/// notify error subsystem of beginning of element loop
void EnterElementLoop();

/// notify error subsystem of end of element loop
void ExitElementLoop();

/// notify error subsystem of element error
void ElementError(int ele, const std::string& err);

#endif // THROWELEMENTERRORS
#endif // #ifdef __cplusplus
#endif // #ifndef DSERROR_H

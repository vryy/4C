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
#else

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

void cpp_dsassert_func(const char* file, const int line, const int test, const char* string);
void cpp_dserror_func(const char* string, ...);
void cpp_dslatest(const char* file, const int line);
#define dsassert(test, string) cpp_dsassert_func(__FILE__, __LINE__, test, string)
#define dserror cpp_dslatest(__FILE__, __LINE__), cpp_dserror_func

#endif
#endif

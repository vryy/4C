/*!----------------------------------------------------------------------
\file drt_dserror.H
\brief Declaration

<pre>
\brief Declaration
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef DSERROR_H
#define DSERROR_H

#ifdef __cplusplus

#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <string>

extern "C" void cpp_dserror_func(const char* string, ...); // avoid constructing string from char*
void cpp_dserror_func(const std::string string, ...);
void cpp_dslatest(const std::string file, const int line);

inline
void cpp_dsassert_func(const char* file, const int line, const bool test, const char* string)
{
  // inline to avoid function call in release mode
#ifdef DEBUG
  if (!test)
  {
    cpp_dslatest(file, line);
    cpp_dserror_func(string);
  }
#endif
  return;
} /* end of dsassert_func */

// do not emit any code in release mode for dsassert (i.e., 'test' is not evaluated and optimized away)
#ifdef DEBUG
#define dsassert(test, string) cpp_dsassert_func(__FILE__, __LINE__, test, string)
#else
#define dsassert(test, string)
#endif

#define dserror cpp_dslatest(__FILE__, __LINE__), cpp_dserror_func

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

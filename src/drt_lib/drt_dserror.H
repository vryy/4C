/*!----------------------------------------------------------------------
\file drt_dserror.H
\brief

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

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

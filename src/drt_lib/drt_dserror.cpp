/*!----------------------------------------------------------------------
\file drt_dserror.cpp
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

#include "drt_dserror.H"

#include <stdexcept>
#include <string.h>

#include <mpi.h>
#include <execinfo.h>
#include <unistd.h>


static int         latest_line = -1;
static std::string latest_file = "{dserror_func call without prototype}";
/*----------------------------------------------------------------------*
 |  assert function                                          mwgee 11/06|
 | used by macro dsassert in dserror.H                                  |
 *----------------------------------------------------------------------*/
void cpp_dsassert_func(const std::string file, const int line, const bool test, const std::string text)
{
#ifdef DEBUG
  if (!test)
  {
    latest_file = file;
    latest_line = line;
    cpp_dserror_func(text);
  }
#endif
  return;
} /* end of dsassert_func */

/*----------------------------------------------------------------------*
 |  assert function                                          mwgee 11/06|
 | used by macro dsassert in dserror.H                                  |
 *----------------------------------------------------------------------*/
extern "C"
void cpp_dsassert_func(const char* file, const int line, const int test, const char* text)
{
#ifdef DEBUG
  if (!test)
  {
    latest_file = file;
    latest_line = line;
    cpp_dserror_func(text);
  }
#endif
  return;
} /* end of dsassert_func */

/*----------------------------------------------------------------------*
 |  set file and line                                        mwgee 11/06|
 | used by macro dsassert and dserror in dserror.H                      |
 *----------------------------------------------------------------------*/
void cpp_dslatest(const std::string file, const int line)
{
  latest_file = file;
  latest_line = line;
}

/*----------------------------------------------------------------------*
 |  set file and line                                        mwgee 11/06|
 | used by macro dsassert and dserror in dserror.H                      |
 *----------------------------------------------------------------------*/
extern "C"
void cpp_dslatest(const char* file, const int line)
{
  latest_file = file;
  latest_line = line;
}

/*----------------------------------------------------------------------*
 | error function                                            mwgee 11/06|
 | used by macro dsassert in dserror.H                                  |
 *----------------------------------------------------------------------*/
extern "C"
void cpp_dserror_func(const char* text, ...)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  char errbuf[4096];

  va_list ap;
  va_start(ap, text);

  sprintf(errbuf,"PROC %d ERROR in %s, line %i:\n",myrank,latest_file.c_str(),latest_line);
  vsprintf(&errbuf[strlen(errbuf)],text,ap);

#ifdef ENABLE_STACKTR
// print stacktrace
  int nptrs;
  void *buffer[100];
  char **strings;

  nptrs   = backtrace(buffer, 100);
  strings = backtrace_symbols(buffer, nptrs);

  sprintf(&errbuf[strlen(errbuf)], "\n\n--- stacktrace ---");
  for (int j = 0; j < nptrs; ++j)
    sprintf(&errbuf[strlen(errbuf)], "\n%s", strings[j]);

  free(strings);
#endif

  va_end(ap);

  throw std::runtime_error(errbuf);
} /* end of dserror_func */

/*----------------------------------------------------------------------*
 | error function                                            mwgee 11/06|
 | used by macro dsassert in dserror.H                                  |
 *----------------------------------------------------------------------*/
void cpp_dserror_func(const std::string text, ...)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  char errbuf[4096];

  va_list ap;
  va_start(ap, text);

  sprintf(errbuf,"PROC %d ERROR in %s, line %i:\n",myrank,latest_file.c_str(),latest_line);
  vsprintf(&errbuf[strlen(errbuf)],text.c_str(),ap);

#ifdef ENABLE_STACKTR
  // print stacktrace
  int nptrs;
  void *buffer[100];
  char **strings;

  nptrs   = backtrace(buffer, 100);
  strings = backtrace_symbols(buffer, nptrs);

  sprintf(&errbuf[strlen(errbuf)], "\n\n--- stacktrace ---");
  for (int j = 0; j < nptrs; ++j)
    sprintf(&errbuf[strlen(errbuf)], "\n%s", strings[j]);

  free(strings);
#endif

  va_end(ap);

  throw std::runtime_error(errbuf);
} /* end of dserror_func */



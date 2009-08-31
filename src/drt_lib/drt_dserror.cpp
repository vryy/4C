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
#ifdef CCADISCRET

#include "drt_dserror.H"
#include "../drt_lib/standardtypes_cpp.H"

#ifdef THROWELEMENTERRORS

#include "drt_globalproblem.H"
#include "../drt_io/io_control.H"


static bool elementcall = false;
static int elementerrorcount = 0;


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void EnterElementLoop()
{
  if (elementerrorcount!=0)
    dserror("inconsistent error system: %d", elementerrorcount);
  elementcall = true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ExitElementLoop()
{
  elementcall = false;
  if (elementerrorcount!=0)
    dserror("%d element errors occured",elementerrorcount);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ElementError(int ele, const std::string& err)
{
  if (not elementcall)
    dserror("element error outside of element loop");
  elementerrorcount += 1;
  fprintf(DRT::Problem::Instance()->ErrorFile()->Handle(),"element %d error: %s\n",ele,err.c_str());
  fflush(DRT::Problem::Instance()->ErrorFile()->Handle());
}


#endif

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
 |  set file and line                                        mwgee 11/06|
 | used by macro dsassert and dserror in dserror.H                      |
 *----------------------------------------------------------------------*/
void cpp_dslatest(const std::string file, const int line)
{
  latest_file = file;
  latest_line = line;
}

/*----------------------------------------------------------------------*
 | error function                                            mwgee 11/06|
 | used by macro dsassert in dserror.H                                  |
 *----------------------------------------------------------------------*/
void cpp_dserror_func(const std::string text, ...)
{
  int myrank;
#ifdef PARALLEL
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#else
  myrank=0;
#endif

  char errbuf[2048];

  va_list ap;
  va_start(ap, text);

  sprintf(errbuf,"PROC %d ERROR in %s, line %i:\n",myrank,latest_file.c_str(),latest_line);
  vsprintf(&errbuf[strlen(errbuf)],text.c_str(),ap);

  va_end(ap);

#ifdef THROWELEMENTERRORS
  if (elementcall)
  {
    throw std::string(errbuf);
  }
  else
#endif
  {
    char line[] = "=========================================================================\n";
    printf("\n\n");
    printf(line);
    printf(errbuf);
    printf("\n");
    printf(line);
    printf("\n\n");
    fflush(stdout);

#ifdef DSERROR_DUMP
    abort();
#endif

#ifdef PARALLEL
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
    exit(0);
#endif
  }
} /* end of dserror_func */



#endif  // #ifdef CCADISCRET

/*!----------------------------------------------------------------------
\file drt_dserror.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_dserror.H"

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}

#ifdef THROWELEMENTERRORS

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;


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
  fprintf(allfiles.out_err,"element %d error: %s\n",ele,err.c_str());
  fflush(allfiles.out_err);
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
    *((int*)0x0)=0;
#endif

#ifdef PARALLEL
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
    exit(0);
#endif
  }
} /* end of dserror_func */



#endif  // #ifdef CCADISCRET

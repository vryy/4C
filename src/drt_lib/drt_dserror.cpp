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
  va_list ap;
  char line[] = "=========================================================================\n";
  va_start(ap, text);

int myrank;
#ifdef PARALLEL
MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#else
myrank=0;
#endif

  printf("\n");
  printf("\n");
  printf(line);
  printf("PROC %d ERROR in %s, line %i:\n",myrank,latest_file.c_str(),latest_line);
  vprintf(text.c_str(),ap);
  printf("\n");
  printf(line);
  printf("\n");
  printf("\n");

  va_end(ap);

  fflush(stdout);

#ifdef DSERROR_DUMP
  *((int*)0x0)=0;
#endif

#ifdef PARALLEL
  MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
  exit(0);
#endif

  return;
} /* end of dserror_func */



#endif  // #ifdef CCADISCRET

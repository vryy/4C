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

static int   latest_line = -1;
static char* latest_file = "{dserror_func call without prototype}";
/*----------------------------------------------------------------------*
 |  assert function                                          mwgee 11/06|
 | used by macro dsassert in dserror.H                                  |
 *----------------------------------------------------------------------*/
void cpp_dsassert_func(char* file, int line, int test, char string[])
{
#ifdef DEBUG
  if (!test)
  {
    latest_file = file;
    latest_line = line;
    cpp_dserror_func(string);
  }
#endif
  return;
} /* end of dsassert_func */

/*----------------------------------------------------------------------*
 |  set file and line                                        mwgee 11/06|
 | used by macro dsassert and dserror in dserror.H                      |
 *----------------------------------------------------------------------*/
void cpp_dslatest(char* file, int line)
{
  latest_file = file;
  latest_line = line;
}

/*----------------------------------------------------------------------*
 | error function                                            mwgee 11/06|
 | used by macro dsassert in dserror.H                                  |
 *----------------------------------------------------------------------*/
void cpp_dserror_func(char *string, ...)
{
  va_list ap;
  char line[] = "=========================================================================\n";
  va_start(ap, string);

int myrank;
#ifdef PARALLEL
MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#else
myrank=0;
#endif  

  printf("\n");
  printf("\n");
  printf(line);
  printf("PROC %d ERROR in %s, line %i:\n",myrank,latest_file,latest_line);
  vprintf(string,ap);
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

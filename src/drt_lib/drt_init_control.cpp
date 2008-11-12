/*----------------------------------------------------------------------*/
/*!
\file drt_init_control.cpp

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <string>

#include "drt_init_control.H"
#include "standardtypes_cpp.H"


/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | struct _GENPROB       genprob; defined in global_control.c           |
 *----------------------------------------------------------------------*/
extern struct _GENPROB  genprob;


/*----------------------------------------------------------------------*/
/*
  Setup of input and output files. No actual read is performed
  here.
 */
/*----------------------------------------------------------------------*/
void ntaini_ccadiscret(int argc, char** argv)
{
  ds_cputime_init();

  int myrank = 0;

#ifdef PARALLEL
  int nproc  = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#endif

  if (argc <= 1)
  {
    if (myrank==0)
    {
      printf("You forgot to give the input and output file names!\n");
      printf("Try again!\n");
    }
#ifdef PARALLEL
    MPI_Finalize();
#endif
    exit(1);
  }
  else if (argc <= 2)
  {
    if (myrank==0)
    {
      printf("You forgot to give the output file name!\n");
      printf("Try again!\n");
    }
#ifdef PARALLEL
    MPI_Finalize();
#endif
    exit(1);
  }

  allfiles.outlenght = strlen(argv[2]);
  allfiles.outputfile_kenner = argv[2];
  if (allfiles.outlenght>=100)
  {
    if (myrank==0)
    {
      printf("Your outputfile kenner is too long!\n");
      fflush(stdout);
    }
#ifdef PARALLEL
    MPI_Finalize();
#endif
    exit(1);
  }

  allfiles.inputfile_name = argv[1];
  // set error file names
  sprintf(allfiles.outputfile_name, "%s%d.err",
          allfiles.outputfile_kenner, myrank);
  // REMARK:
  // error files are opened by OpenErrorFile()
  // called in ntainp_ccadiscret()

  // inform user
  if (myrank==0)
  {
    printf("input is read from         %s\n", allfiles.inputfile_name);
    printf("errors are reported to     %s\n", allfiles.outputfile_name);
  }

  /*-------------------------------------------------- check for restart */
  genprob.restart = 0;
  if (argc > 3)
  {
    std::string restart(argv[3]);
    if (restart=="restart")
      genprob.restart++;
  }
}

#endif

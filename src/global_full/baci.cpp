/*!
\brief main routine
*/

#include <iostream>
#include <stdexcept>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include <../headers/compile_settings.h>
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/drt_validconditions.H"
#include "../drt_inpar/drt_validmaterials.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_elementdefinition.H"
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_parobjectregister.H"

#ifdef TRAP_FE

#ifdef LINUX_MUENCH
#include <fenv.h>
#endif

#ifdef HPUX_MUENCH
#include <fenv.h>
#endif

#ifdef HPUXITA
#include <fenv.h>
#endif

#ifdef HPUX11
#include <fenv.h>
#endif

#endif /* TRAP_FE */

/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
struct _PAR     par;

void ntam(
    int                 argc,
    char               *argv[]
  );

/*!

\brief main routine

<pre>                                                        m.gee 8/00
main is only printing the ccarat head and the finish
</pre>
\param argc     INT     (i)   number of arguments on command line including exe
\param argv     *char[] (i)   array of arguments from command line

*/
int main(int argc, char *argv[])
{
#ifdef PARALLEL
  char *buff,*dbuff;
  int   buffsize=MPIBUFFSIZE;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &par.myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &par.nprocs);

  /*------------------------------------------------ attach buffer to mpi */
  buff = (char*)malloc(buffsize);
  if (!buff)
  {
    printf("Allocation of memory for mpi buffer failed");
    MPI_Finalize();
    exit(1);
  }
  MPI_Buffer_attach(buff,buffsize);
#else
  par.myrank=0;
  par.nprocs=1;
#endif

  if (par.myrank==0)
  {
    printf("\n"
           "****************************************\n"
           "*                                      *\n"
           "*               B A C I                *\n"
           "*                                      *\n"
           "*                                      *\n"
           "*            revision %5d            *\n"
#ifdef PARALLEL
           "*           parallel version           *\n"
#else
           "*          sequential version          *\n"
#endif
#ifdef DEBUG
           "*            debug version             *\n"
#else
           "*             fast version             *\n"
#endif
           "*                                      *\n"
           "*  Lehrstuhl fuer Numerische Mechanik  *\n"
           "*                 LNM                  *\n"
           "*   Technische Universitaet Muenchen   *\n"
           "*                                      *\n"
           "*    (c) 2010 All Rights Reserved.     *\n"
           "*                                      *\n"
           "****************************************\n\n",
           CHANGEDREVISION+0);
#ifdef PARALLEL
    printf("number of processors: %d\n",par.nprocs);
#endif
  }

  if ((argc == 2) && (strcmp(argv[1], "-v") == 0)) {
    if (par.myrank==0) {
      PrintParObjectList();
      printf("\n\n");
    }
  }
  else if ((argc == 2) &&
           ((strcmp(argv[1], "-p") == 0) ||
            (strcmp(argv[1], "--parameters") == 0)))
  {
    if (par.myrank==0)
    {
      printf("\n\n");
      PrintValidParameters();
      printf("\n\n");
    }
  }
  else if ((argc == 2) &&
           ((strcmp(argv[1], "-d") == 0) ||
            (strcmp(argv[1], "--datfile") == 0)))
  {
    if (par.myrank==0)
    {
      printf("\n\n");
      PrintDefaultDatHeader();
      PrintConditionDatHeader();
      PrintMaterialDatHeader();
      PrintElementDatHeader();
      PrintFunctionDatHeader();
      PrintTimeCurveDatHeader();
      PrintResultDescrDatHeader();
      printf("\n\n");
    }
  }
  else {
    /* Here we turn the NaN and inf numbers of. No need to calculate
     * those. If those appear the calculation needs much (!) more
     * time. Better stop immediately if some illegal operation occurs. */
#ifdef TRAP_FE

    /* Sadly, it seems the functions needed for this are different on
     * different maschines. */
#ifdef LINUX_MUENCH

    /* This is a GNU extention thus it's only available on linux. But
     * it's exactly what we want: SIGFPE just for the given
     * exceptions. We don't care about FE_INEXACT. (It happens all the
     * time.) */
    /* Over- and underflow seem to happen sometimes. Does it worry us?
     * Will that spoil the results? */
    /*feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW);*/
    feenableexcept(FE_INVALID | FE_DIVBYZERO);

    /* The hard GNU way. But it does too much. */
    /*fesetenv((fenv_t*)-2);*/

#endif

#ifdef HPUX_MUENCH
    /*
     * Don't ask me why they want this. The man page said it's needed on
     * itanium maschines. */
#pragma STDC FENV_ACCESS ON
    /*fesettrapenable(FE_INVALID | FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW);*/
    fesettrapenable(FE_INVALID | FE_DIVBYZERO);
#endif

#ifdef HPUXITA
#pragma STDC FENV_ACCESS ON
    fesettrapenable(FE_INVALID | FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW);
#endif

#ifdef HPUX11
#pragma STDC FENV_ACCESS ON
    fesettrapenable(FE_INVALID | FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW);
#endif

#endif /* TRAP_FE */

/*----------------------------------------------- everything is in here */
#ifdef DSERROR_DUMP
      ntam(argc,argv);
#else
    try
    {
      ntam(argc,argv);
    }
    catch ( std::runtime_error & err )
    {
      DRT::Problem::Done();

      char line[] = "=========================================================================\n";
      std::cout << "\n\n"
                << line
                << err.what()
                << "\n"
                << line
                << "\n" << std::endl;

#ifdef PARALLEL
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
      exit(1);
#endif
    }
#endif
/*----------------------------------------------------------------------*/
  }

  DRT::Problem::Done();

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  printf("processor %d finished normally\n",par.myrank);
  MPI_Buffer_detach(&dbuff,&buffsize);
  if (dbuff!=buff || buffsize != MPIBUFFSIZE)
    dserror("Illegal modification of mpi buffer adress or size appeared");
  free(dbuff);
  MPI_Finalize();
#else
  printf("processor %d finished normally\n",par.myrank);
#endif
  return(0);
}

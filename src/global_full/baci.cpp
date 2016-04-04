/*----------------------------------------------------------------------*/
/*!
\file baci.cpp

\maintainer Martin Kronbichler

\brief Main routine
*/
/*----------------------------------------------------------------------*/

#include <iostream>
#include <stdexcept>
#include <Epetra_MpiComm.h>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
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
#include "../drt_lib/drt_utils_createdis.H"

#ifdef TRAP_FE
#include <fenv.h>
#endif /* TRAP_FE */


/*----------------------------------------------------------------------*
 | size of buffer to attach to intra-communicator in byte               |
 *----------------------------------------------------------------------*/
#define MPIBUFFSIZE      (52428800) /* this is 50 MB */


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
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
  char *buff,*dbuff;
  int   buffsize=MPIBUFFSIZE;

  MPI_Init(&argc,&argv);

  COMM_UTILS::CreateComm(argc,argv);

  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<Epetra_Comm> lcomm = Teuchos::rcp(problem->GetNPGroup()->LocalComm().get(), false);
  Teuchos::RCP<Epetra_Comm> gcomm = Teuchos::rcp(problem->GetNPGroup()->GlobalComm().get(), false);
  int ngroups = problem->GetNPGroup()->NumGroups();

  if (strcmp(argv[argc-1], "--interactive") == 0)
  {
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("Global rank %d with PID %d on %s is ready for attach\n", gcomm->MyPID(), getpid(), hostname);
    if (gcomm->MyPID() == 0)
    {
      printf( "\n** Enter a character to continue > \n"); fflush(stdout);
      char go = ' ';
      scanf("%c",&go);
    }
  }

  gcomm->Barrier();

  /*------------------------------------------------ attach buffer to mpi */
  buff = (char*)malloc(buffsize);
  if (!buff)
  {
    printf("Allocation of memory for mpi buffer failed");
    MPI_Finalize();
    exit(1);
  }
  MPI_Buffer_attach(buff,buffsize);

  if (gcomm->MyPID()==0)
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
    printf("Baci SHA1: %s\n", Baci_SHA1);
    printf("Trilinos SHA1: %s\n",Trilinos_SHA1);
    printf("Total number of processors: %d\n",gcomm->NumProc());
  }

  if ((argc == 2) && (strcmp(argv[1], "-v") == 0)) {
    if (lcomm->MyPID()==0) {
      PrintParObjectList();
      printf("\n\n");
    }
  }
  else if ((argc == 2) &&
           ((strcmp(argv[1], "-h") == 0) ||
            (strcmp(argv[1], "--help") == 0)))
  {
    if (lcomm->MyPID()==0)
    {
      printf("\n\n");
      PrintHelpMessage();
      printf("\n\n");
    }
  }
  else if ((argc == 2) &&
           ((strcmp(argv[1], "-p") == 0) ||
            (strcmp(argv[1], "--parameters") == 0)))
  {
    if (lcomm->MyPID()==0)
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
    if (lcomm->MyPID()==0)
    {
      printf("\n\n");
      PrintDefaultDatHeader();
      PrintConditionDatHeader();
      PrintMaterialDatHeader();
      DRT::UTILS::PrintCloningMaterialMapDatHeader();
      PrintElementDatHeader();
      PrintFunctionDatHeader();
      PrintTimeCurveDatHeader();
      PrintResultDescrDatHeader();
      printf("\n\n");
    }
  }
  else {
    /* Here we turn the NaN and INF numbers off. No need to calculate
     * those. If those appear, the calculation needs much (!) more
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
    feclearexcept(FE_ALL_EXCEPT);
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
      char line[] = "=========================================================================\n";
      std::cout << "\n\n"
                << line
                << err.what()
                << "\n"
                << line
                << "\n" << std::endl;

      if(ngroups > 1)
      {
        printf("Global processor %d has thrown an error and is waiting for the remaining procs\n\n",gcomm->MyPID());
        gcomm->Barrier();
      }

      DRT::Problem::Done();

      MPI_Buffer_detach(&dbuff,&buffsize);
      if (dbuff!=buff || buffsize != MPIBUFFSIZE)
        dserror("Illegal modification of mpi buffer adress or size appeared");
      free(dbuff);
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
#endif
/*----------------------------------------------------------------------*/
  }

  lcomm->Barrier();
  if(ngroups > 1)
  {
    printf("Global processor %d with local rank %d finished normally\n",gcomm->MyPID(),lcomm->MyPID());
    gcomm->Barrier();
  }
  else
  {
    printf("processor %d finished normally\n",lcomm->MyPID());
  }

  DRT::Problem::Done();

  MPI_Buffer_detach(&dbuff,&buffsize);
  if (dbuff!=buff || buffsize != MPIBUFFSIZE)
    dserror("Illegal modification of mpi buffer adress or size appeared");
  free(dbuff);
  MPI_Finalize();

  return(0);
}

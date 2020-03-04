/*----------------------------------------------------------------------*/
/*! \file

\maintainer Martin Kronbichler

\brief Main routine

\level 0
*/
/*----------------------------------------------------------------------*/

#include <iostream>
#include <stdexcept>
#include <csignal>
#include <Epetra_MpiComm.h>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include <../headers/compile_settings.h>
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/drt_validconditions.H"
#include "../drt_inpar/drt_validmaterials.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_elementdefinition.H"
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_parobjectregister.H"
#include "../drt_lib/drt_utils_createdis.H"

#include <revision.H>

#ifdef TRAP_FE
#include <fenv.h>
#endif /* TRAP_FE */


/*----------------------------------------------------------------------*
 | size of buffer to attach to intra-communicator in byte               |
 *----------------------------------------------------------------------*/
#define MPIBUFFSIZE (52428800) /* this is 50 MB */

/* Collect and print data on memory high water mark of this run
 *
 * 1. Ask the operating system for memory usage.
 * 2. Compute min/max/average and total memory usage across all MPI ranks.
 * 3. Print a summary to the screen.
 *
 * If status file can't be opened, issue a message to the screen. Do not throw an error, since this
 * is not considered a critical failure during a simulation.
 *
 * @param[in] global_comm Global Epetra_Comm object
 */
void GetMemoryHighWaterMark(const Epetra_Comm &comm)
{
  const std::string status_match = "VmHWM";
  const std::string status_filename = "/proc/" + std::to_string(getpid()) + "/status";
  std::ifstream status_file(status_filename, std::ios_base::in);
  double local_mem = std::nan("0");

  // Retrieve local memory use on each process
  if (status_file.is_open())
  {
    std::string line;
    while (std::getline(status_file, line))
    {
      if (line.find(status_match) != std::string::npos)
      {
        size_t start = line.find_first_of("1234567890");
        size_t stop = line.find_last_of("1234567890");

        std::stringstream(line.substr(start, stop + 1)) >> local_mem;
        break;
      }
    }
    status_file.close();

    // Convert memory from KB to GB
    local_mem /= (1 << 20);

    // Gather values
    const int num_procs = comm.NumProc();
    double *recvbuf = new double[num_procs];
    MPI_Gather(&local_mem, 1, MPI_DOUBLE, recvbuf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Compute and output statistics on proc 0
    if (comm.MyPID() == 0)
    {
      double mem_min = recvbuf[0];
      double mem_max = recvbuf[0];
      double mem_tot = 0.0;
      double mem_avg = 0.0;

      for (int rank = 0; rank < num_procs; ++rank)
      {
        mem_min = std::min(mem_min, recvbuf[rank]);
        mem_max = std::max(mem_max, recvbuf[rank]);
        mem_tot += recvbuf[rank];
      }

      mem_avg = mem_tot / num_procs;

      if (num_procs > 1)
      {
        std::cout << std::scientific << std::setprecision(4) << "\nMemory High Water Mark Summary:"
                  << "\t\tMinOverProcs\tMeanOverProcs\tMaxOverProcs\tSumOverProcs\n"
                  << "(in GB)\t\t\t\t\t" << mem_min << "\t" << mem_avg << "\t" << mem_max << "\t"
                  << mem_tot << "\n"
                  << std::endl;
      }
      else
      {
        std::cout << std::scientific << std::setprecision(4)
                  << "\nMemory High Water Mark Summary:\t\tTotal\n"
                  << "(in GB)\t\t\t\t\t" << mem_tot << "\n"
                  << std::endl;
      }
    }
  }
  else  // Failed to open the status file
  {
    std::cout << "Memory High Water Mark summary can not be generated,\nsince status file "
              << status_filename << " could not be opened.\n"
              << std::endl;
  }
}

/*!
 * \brief FPE signal handle
 *
 * A function to handle floating point exceptions by raising a dserror.
 * So we get a stack-trace also on systems where this is not provided
 * through core-dumps from MPI_Abort() (e.g. OpenMPI does whereas
 * Intel MPI doesn't).
 */
void sigfpe_handler(int sig) { dserror("Baci produced a floating point exception."); }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ntam(int argc, char *argv[]);

/*!

\brief main routine

\author m.gee \date 8/00
main is only printing the ccarat head and the finish
\param argc     INT     (i)   number of arguments on command line including exe
\param argv     *char[] (i)   array of arguments from command line

*/
int main(int argc, char *argv[])
{
  char *buff, *dbuff;
  int buffsize = MPIBUFFSIZE;

  MPI_Init(&argc, &argv);

  COMM_UTILS::CreateComm(argc, argv);

  DRT::Problem *problem = DRT::Problem::Instance();
  Teuchos::RCP<Epetra_Comm> lcomm = Teuchos::rcp(problem->GetNPGroup()->LocalComm().get(), false);
  Teuchos::RCP<Epetra_Comm> gcomm = Teuchos::rcp(problem->GetNPGroup()->GlobalComm().get(), false);
  int ngroups = problem->GetNPGroup()->NumGroups();

  if (strcmp(argv[argc - 1], "--interactive") == 0)
  {
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("Global rank %d with PID %d on %s is ready for attach\n", gcomm->MyPID(), getpid(),
        hostname);
    if (gcomm->MyPID() == 0)
    {
      printf("\n** Enter a character to continue > \n");
      fflush(stdout);
      char go = ' ';
      DRT::UTILS::Checkscanf(scanf("%c", &go));
    }
  }

  gcomm->Barrier();

  /*------------------------------------------------ attach buffer to mpi */
  buff = (char *)malloc(buffsize);
  if (!buff)
  {
    printf("Allocation of memory for mpi buffer failed");
    MPI_Finalize();
    exit(1);
  }
  MPI_Buffer_attach(buff, buffsize);

  if (gcomm->MyPID() == 0)
  {
    printf(
        "\n"
        "**********************************************\n"
        "*                                            *\n"
        "*                  B A C I                   *\n"
        "*                                            *\n"
        "*                                            *\n"
        "*             version (git SHA1):            *\n"
        "*  %s  *\n"
        "*                                            *\n"
#ifdef PARALLEL
        "*              parallel version              *\n"
#else
        "*             sequential version             *\n"
#endif
#ifdef DEBUG
        "*               debug version                *\n"
#else
        "*                fast version                *\n"
#endif
        "*                                            *\n"
        "*     Lehrstuhl fuer Numerische Mechanik     *\n"
        "*                    LNM                     *\n"
        "*      Technische Universitaet Muenchen      *\n"
        "*                                            *\n"
        "*       (c) 2010 All Rights Reserved.        *\n"
        "*                                            *\n"
        "**********************************************\n\n",
        BaciGitHash.c_str());
    printf("Trilinos SHA1: %s\n", Trilinos_SHA1);
    printf("Total number of processors: %d\n", gcomm->NumProc());
  }

  if ((argc == 2) && (strcmp(argv[1], "-v") == 0))
  {
    if (lcomm->MyPID() == 0)
    {
      PrintParObjectList();
      printf("\n\n");
    }
  }
  else if ((argc == 2) && ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0)))
  {
    if (lcomm->MyPID() == 0)
    {
      printf("\n\n");
      PrintHelpMessage();
      printf("\n\n");
    }
  }
  else if ((argc == 2) && ((strcmp(argv[1], "-p") == 0) || (strcmp(argv[1], "--parameters") == 0)))
  {
    if (lcomm->MyPID() == 0)
    {
      printf("\n\n");
      PrintValidParameters();
      printf("\n\n");
    }
  }
  else if ((argc == 2) && ((strcmp(argv[1], "-d") == 0) || (strcmp(argv[1], "--datfile") == 0)))
  {
    if (lcomm->MyPID() == 0)
    {
      printf("\n\n");
      PrintDefaultDatHeader();
      PrintConditionDatHeader();
      PrintMaterialDatHeader();
      DRT::UTILS::PrintCloningMaterialMapDatHeader();
      PrintElementDatHeader();
      PrintFunctionDatHeader();
      PrintResultDescrDatHeader();
      printf("\n\n");
    }
  }
  else
  {
    /* Here we turn the NaN and INF numbers off. No need to calculate
     * those. If those appear, the calculation needs much (!) more
     * time. Better stop immediately if some illegal operation occurs. */
#ifdef TRAP_FE

    /* Sadly, it seems the functions needed for this are different on
     * different machines. */
#ifdef LINUX_MUENCH

    /* This is a GNU extension thus it's only available on linux. But
     * it's exactly what we want: SIGFPE just for the given
     * exceptions. We don't care about FE_INEXACT. (It happens all the
     * time.) */
    /* Over- and underflow seem to happen sometimes. Does it worry us?
     * Will that spoil the results? */
    /*feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW);*/
    feclearexcept(FE_ALL_EXCEPT);
    feenableexcept(FE_INVALID | FE_DIVBYZERO);

    // Initialize a signal handle for SIGFPE
    struct sigaction act;
    act.sa_handler = sigfpe_handler;
    sigemptyset(&act.sa_mask);
    act.sa_flags = 0;
    sigaction(SIGFPE, &act, 0);

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
    ntam(argc, argv);
#else
    try
    {
      ntam(argc, argv);
    }
    catch (std::runtime_error &err)
    {
      char line[] = "=========================================================================\n";
      std::cout << "\n\n" << line << err.what() << "\n" << line << "\n" << std::endl;

      if (ngroups > 1)
      {
        printf("Global processor %d has thrown an error and is waiting for the remaining procs\n\n",
            gcomm->MyPID());
        gcomm->Barrier();
      }

      DRT::Problem::Done();

      MPI_Buffer_detach(&dbuff, &buffsize);
      if (dbuff != buff || buffsize != MPIBUFFSIZE)
        dserror("Illegal modification of mpi buffer adress or size appeared");
      free(dbuff);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
#endif
    /*----------------------------------------------------------------------*/
  }

  GetMemoryHighWaterMark(*gcomm);

  lcomm->Barrier();
  if (ngroups > 1)
  {
    printf("Global processor %d with local rank %d finished normally\n", gcomm->MyPID(),
        lcomm->MyPID());
    gcomm->Barrier();
  }
  else
  {
    printf("processor %d finished normally\n", lcomm->MyPID());
  }

  DRT::Problem::Done();

  MPI_Buffer_detach(&dbuff, &buffsize);
  if (dbuff != buff || buffsize != MPIBUFFSIZE)
    dserror("Illegal modification of mpi buffer adress or size appeared");
  free(dbuff);
  MPI_Finalize();

  return (0);
}

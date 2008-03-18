
#ifdef CCADISCRET

#include <string>
#include <iostream>

#ifdef PARALLEL
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "elch_dyn.H"
#include "elch_algorithm.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*/
// entry point for ELCH in DRT
/*----------------------------------------------------------------------*/
void elch_dyn()
{

// create a communicator
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

if (comm.MyPID() == 0)
cout<<"ELCH problemtype under development..."<<endl;

Teuchos::RCP<ELCH::Algorithm> elch = Teuchos::rcp(new ELCH::Algorithm(comm));
//elch->FluidField()->Update();

// summarize performance measurements
Teuchos::TimeMonitor::summarize();

} // elch_dyn()

#endif  // CCADISCRET

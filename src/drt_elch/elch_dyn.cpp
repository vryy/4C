
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
#include "elch_create_condif.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include "../drt_lib/drt_globalproblem.H"
#include <Epetra_Time.h>
#if 0
#include "../drt_io/io_gmsh.H"
#endif
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

  // print ELCH-Logo to screen
  if (comm.MyPID()==0) printlogo();

  // get discretization ids
  int disnumff = genprob.numff; // typically 0
  int disnumscatra = genprob.numscatra; // typically 1
  
  // access the fluid discretization
  RefCountPtr<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(disnumff,0);
  if (!fluiddis->Filled()) fluiddis->FillComplete();

  if (fluiddis->NumGlobalNodes()==0)
    dserror("No fluid discretization found!");
#if 0
  std::ofstream f_system("mydiscretization.pos");
  f_system<<IO::GMSH::disToString("Fluid",0,fluiddis);
#endif  
  // access the (typically empty) condif discretization
  RefCountPtr<DRT::Discretization> condifdis = DRT::Problem::Instance()->Dis(disnumscatra,0);
  if (!condifdis->Filled()) condifdis->FillComplete();

  // create condif elements if the condif discretization is empty
  if (condifdis->NumGlobalNodes()==0)
  {
    Epetra_Time time(comm);
    ELCH::CreateConDifDiscretization(disnumff,disnumscatra);
    if (comm.MyPID()==0)
    cout<<"Created necessary condif discretization from fluid field in...."
    <<time.ElapsedTime() << " secs\n\n";
  }
  else
    dserror("Fluid AND ConDif discretization present. This is not supported.");

  // access the problem-speific parameter list
  const Teuchos::ParameterList& elchdyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();

  // create an ELCH::Algorithm instance
  Teuchos::RCP<ELCH::Algorithm> elch = Teuchos::rcp(new ELCH::Algorithm(comm,elchdyn));

  if (genprob.restart)
  {
    // read the restart information, set vectors and variables
    //elch->ReadRestart(genprob.restart);
    dserror("restart not yet available");
    exit(1);
  }
  
  // solve the whole electrochemistry problem
  elch->TimeLoop();

  // summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // perform the result test
  DRT::ResultTestManager testmanager(comm);
  testmanager.AddFieldTest(elch->FluidField().CreateFieldTest());
  testmanager.AddFieldTest(elch->CreateScaTraFieldTest());
  testmanager.TestAll();

  return;

} // elch_dyn()


/*----------------------------------------------------------------------*/
// print ELCH-Module logo
/*----------------------------------------------------------------------*/
void printlogo()
{
    // more at http://www.ascii-art.de under entry "moose" (or "elk")
    cout<<"     ___            ___    "<<endl;
    cout<<"    /   \\          /   \\ "<<endl;
    cout<<"    \\_   \\        /  __/ "<<endl;
    cout<<"     _\\   \\      /  /__  "<<"     _____ _     _____  _   _   "<<endl;
    cout<<"     \\___  \\____/   __/  "<<"    |  ___| |   /  __ \\| | | |  "<<endl;
    cout<<"         \\_       _/      "<<"   | |__ | |   | /  \\/| |_| |  "<<endl;
    cout<<"           | @ @  \\_      "<<"   |  __|| |   | |    |  _  |   "<<endl;
    cout<<"           |               "<<"  | |___| |___| \\__/\\| | | | "<<endl;
    cout<<"         _/     /\\        "<<"   \\____/\\_____/\\____/\\_| |_/ "<<endl;
    cout<<"        /o)  (o/\\ \\_     "<<endl;
    cout<<"        \\_____/ /         "<<endl;
    cout<<"          \\____/          "<<endl;
    cout<<"                           "<<endl;
}

#endif  // CCADISCRET

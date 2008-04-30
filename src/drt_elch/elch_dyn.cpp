
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
  int disnumff = genprob.numff;
  int disnumcdf = genprob.numcdf;
  
  // get the fluid discretization
  RefCountPtr<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(disnumff,0);
  if (!fluiddis->Filled()) fluiddis->FillComplete();

  // get the condif discretization
  RefCountPtr<DRT::Discretization> condifdis = DRT::Problem::Instance()->Dis(disnumcdf,0);
  if (!condifdis->Filled()) condifdis->FillComplete();

  if (fluiddis->NumGlobalNodes()==0)
    dserror("No fluid discretization found!");

  // create condif elements if the condif discretization is empty
  if (condifdis->NumGlobalNodes()==0)
  {
    ELCH::CreateConDifDiscretization(disnumff,disnumcdf);
    if (comm.MyPID()==0)
    cout<<"Created necessary condif discretization from fluid field.\n\n";
  }
  else
    dserror("Fluid AND ConDif discretization present. Not supported.");

  // create an ELCH::Algorithm instance
  Teuchos::RCP<ELCH::Algorithm> elch = Teuchos::rcp(new ELCH::Algorithm(comm));

  // solve the whole electrochemistry problem
  elch->TimeLoop();

  // summarize the performance measurements
  //Teuchos::TimeMonitor::summarize();

  // perform the result test (only fluid field up to now!)
  DRT::ResultTestManager testmanager(comm);
  testmanager.AddFieldTest(elch->FluidField().CreateFieldTest());
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
    //cout << "\n" GRAY_LIGHT "Checking results ..." YELLOW "tttt" END_COLOR "\n";
}

#endif  // CCADISCRET

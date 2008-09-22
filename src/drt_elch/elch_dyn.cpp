/*!----------------------------------------------------------------------
\file elch_dyn.cpp
\brief Control routine for Electrochemistry module.


<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/

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
#include "../drt_scatra/scatra_utils.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include "../drt_lib/drt_globalproblem.H"
#include <Epetra_Time.h>
#if 0
#include "../drt_io/io_gmsh.H"
#endif


/*----------------------------------------------------------------------*/
// entry point for ELCH in DRT
/*----------------------------------------------------------------------*/
void elch_dyn(int disnumff,int disnumscatra, int restart)
{
  // create a communicator
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // print ELCH-Logo to screen
  if (comm.MyPID()==0) printlogo();

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
    std::map<string,string> conditions_to_copy;
    conditions_to_copy.insert(pair<string,string>("TransportDirichlet","Dirichlet"));
    //conditions_to_copy.insert("FluidStressCalc","FluxCalculation"); // a hack
    SCATRA::CreateScaTraDiscretization(fluiddis,condifdis,conditions_to_copy,false);
    if (comm.MyPID()==0)
    cout<<"Created scalar transport discretization from fluid field in...."
    <<time.ElapsedTime() << " secs\n\n";
  }
  else
    dserror("Fluid AND ConDif discretization present. This is not supported.");

  // access the problem-specific parameter list
  const Teuchos::ParameterList& elchcontrol = DRT::Problem::Instance()->ELCHControlParams();

  // create an ELCH::Algorithm instance
  Teuchos::RCP<ELCH::Algorithm> elch = Teuchos::rcp(new ELCH::Algorithm(comm,elchcontrol));

  if (restart)
  {
    // read the restart information, set vectors and variables
    elch->ReadRestart(restart);
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

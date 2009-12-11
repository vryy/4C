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


#ifdef PARALLEL
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "elch_dyn.H"
#include "elch_algorithm.H"
#include "elch_moving_boundary_algorithm.H"
#include "../drt_scatra/scatra_utils.H"
#include "../drt_fsi/fsi_utils.H"
#include "../drt_lib/drt_utils_createdis.H"
#include <Teuchos_TimeMonitor.hpp>
#include "../drt_lib/drt_globalproblem.H"
#if 0
#include "../drt_io/io_gmsh.H"
#endif


/*----------------------------------------------------------------------*/
// entry point for ELCH in DRT
/*----------------------------------------------------------------------*/
void elch_dyn(int disnumff,int disnumscatra,int disnumale,int restart)
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
  if (!fluiddis->Filled() or !fluiddis->HaveDofs()) fluiddis->FillComplete();

#if 0
  std::ofstream f_system("mydiscretization.pos");
  f_system<<IO::GMSH::disToString("Fluid",0,fluiddis);
#endif
  // access the scatra discretization
  RefCountPtr<DRT::Discretization> scatradis = DRT::Problem::Instance()->Dis(disnumscatra,0);
  if (!scatradis->Filled()) scatradis->FillComplete();

  // access the problem-specific parameter list
  const Teuchos::ParameterList& elchcontrol = DRT::Problem::Instance()->ELCHControlParams();

  // access the scalar transport parameter list
  const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();
  const INPAR::SCATRA::VelocityField veltype
    = Teuchos::getIntegralValue<INPAR::SCATRA::VelocityField>(scatradyn,"VELOCITYFIELD");

  // choose algorithm depending on velocity field type
  switch (veltype)
  {
  case INPAR::SCATRA::velocity_zero:  // zero  (see case 1)
  case INPAR::SCATRA::velocity_function:  // function
  {
    // we directly use the elements from the scalar transport elements section
    if (scatradis->NumGlobalNodes()==0)
      dserror("No elements in the ---TRANSPORT ELEMENTS section");

    // create instance of scalar transport basis algorithm (empty fluid discretization)
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatraonly = rcp(new ADAPTER::ScaTraBaseAlgorithm(elchcontrol,false));

    // read the restart information, set vectors and variables
    if (restart) scatraonly->ScaTraField().ReadRestart(restart);

    // set velocity field
    //(this is done only once. Time-dependent velocity fields are not supported)
    (scatraonly->ScaTraField()).SetVelocityField();

    // enter time loop to solve problem with given convective velocity
    (scatraonly->ScaTraField()).TimeLoop();

    // perform the result test if required
    DRT::Problem::Instance()->AddFieldTest(scatraonly->CreateScaTraFieldTest());
    DRT::Problem::Instance()->TestAll(comm);

    break;
  }
  case INPAR::SCATRA::velocity_Navier_Stokes:  // Navier_Stokes
  {
    // we use the fluid discretization as layout for the scalar transport discretization
    if (fluiddis->NumGlobalNodes()==0) dserror("Fluid discretization is empty!");

    // create scatra elements if the scatra discretization is empty
    if (scatradis->NumGlobalNodes()==0)
    {
      Epetra_Time time(comm);

      // fetch the desired material id for the transport elements
      const int matid = scatradyn.get<int>("MATID");

      // create the scatra discretization
      {
      Teuchos::RCP<DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy> > clonewizard =
            Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy>() );
      clonewizard->CreateMatchingDiscretization(fluiddis,scatradis,matid);
      }

      if (comm.MyPID()==0)
        cout<<"Created scalar transport discretization from fluid field in...."
        <<time.ElapsedTime() << " secs\n\n";
    }
    else
      dserror("Fluid AND ScaTra discretization present. This is not supported.");

    RefCountPtr<DRT::Discretization> aledis = DRT::Problem::Instance()->Dis(disnumale,0);
    if (!aledis->Filled()) aledis->FillComplete();
    // is ALE needed or not?
    const int withale = Teuchos::getIntegralValue<int>(elchcontrol,"MOVINGBOUNDARY");

    if (withale==1)
    {
      // create ale elements only if the ale discretization is empty
      if (aledis->NumGlobalNodes()==0)
      {
        Epetra_Time time(comm);
        FSI::UTILS::CreateAleDiscretization();

        if(comm.MyPID()==0)
        {
          cout << "\n\nCreated ALE discretization from fluid field in........"
          <<time.ElapsedTime() << " secs\n\n";
        }
      }
      else
        dserror("Providing an ALE mesh is not supported for problemtype Electrochemistry.");

      // create an ELCH::MovingBoundaryAlgorithm instance
      Teuchos::RCP<ELCH::MovingBoundaryAlgorithm> elch
      = Teuchos::rcp(new ELCH::MovingBoundaryAlgorithm(comm,elchcontrol));

      // read the restart information, set vectors and variables
      if (restart) elch->ReadRestart(restart);

      // solve the whole electrochemistry problem
      elch->TimeLoop();

      // summarize the performance measurements
      Teuchos::TimeMonitor::summarize();

      // perform the result test
      DRT::Problem::Instance()->AddFieldTest(elch->FluidField().CreateFieldTest());
      DRT::Problem::Instance()->AddFieldTest(elch->CreateScaTraFieldTest());
      DRT::Problem::Instance()->TestAll(comm);
    }
    else
    {
      // create an ELCH::Algorithm instance
      Teuchos::RCP<ELCH::Algorithm> elch = Teuchos::rcp(new ELCH::Algorithm(comm,elchcontrol));

      // read the restart information, set vectors and variables
      if (restart) elch->ReadRestart(restart);

      // solve the whole electrochemistry problem
      elch->TimeLoop();

      // summarize the performance measurements
      Teuchos::TimeMonitor::summarize();

      // perform the result test
      DRT::Problem::Instance()->AddFieldTest(elch->FluidField().CreateFieldTest());
      DRT::Problem::Instance()->AddFieldTest(elch->CreateScaTraFieldTest());
      DRT::Problem::Instance()->TestAll(comm);
    }

    break;
  } // case 2
  default:
    dserror("Unknown velocity field type for transport of passive scalar: %d",veltype);
  }

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

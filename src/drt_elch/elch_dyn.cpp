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


#include "elch_dyn.H"
#include "elch_algorithm.H"
#include "elch_moving_boundary_algorithm.H"
#include "../drt_inpar/inpar_elch.H"
#include "../drt_scatra/scatra_utils_clonestrategy.H"
#include "../drt_fsi/fsi_utils.H"
#include "../drt_lib/drt_utils_createdis.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_Time.h>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#if 0
#include "../drt_io/io_gmsh.H"
#endif


/*----------------------------------------------------------------------*/
// entry point for ELCH in DRT
/*----------------------------------------------------------------------*/
void elch_dyn(int restart)
{
  // pointer to problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // access the communicator
  const Epetra_Comm& comm = problem->GetDis("fluid")->Comm();

  // print ELCH-Logo to screen
  if (comm.MyPID()==0) printlogo();

  // access the fluid discretization
  RefCountPtr<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  // access the scatra discretization
  RefCountPtr<DRT::Discretization> scatradis = problem->GetDis("scatra");

  // ensure that all dofs are assigned in the right order; this creates dof numbers with
  //       fluid dof < scatra/elch dof
  fluiddis->FillComplete();
  scatradis->FillComplete();

#if 0
  std::ofstream f_system("mydiscretization.pos");
  f_system<<IO::GMSH::disToString("Fluid",0,fluiddis);
#endif

  // access the problem-specific parameter list
  const Teuchos::ParameterList& elchcontrol = problem->ELCHControlParams();

  // print default parameters to screen
  if (comm.MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, elchcontrol);

  // access the scalar transport parameter list
  const Teuchos::ParameterList& scatradyn = problem->ScalarTransportDynamicParams();
  const INPAR::SCATRA::VelocityField veltype
    = DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatradyn,"VELOCITYFIELD");

  // choose algorithm depending on velocity field type
  switch (veltype)
  {
  case INPAR::SCATRA::velocity_zero:  // zero  (see case 1)
  case INPAR::SCATRA::velocity_function:  // spatial function
  case INPAR::SCATRA::velocity_function_and_curve:  // spatial function and time curve
  {
    // we directly use the elements from the scalar transport elements section
    if (scatradis->NumGlobalNodes()==0)
      dserror("No elements in the ---TRANSPORT ELEMENTS section");

    // get linear solver id from SCALAR TRANSPORT DYNAMIC
    const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
    if (linsolvernumber == (-1))
      dserror("no linear solver defined for ELCH problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

    // create instance of scalar transport basis algorithm (empty fluid discretization)
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatraonly = rcp(new ADAPTER::ScaTraBaseAlgorithm(elchcontrol,false,"scatra",DRT::Problem::Instance()->SolverParams(linsolvernumber),false));

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

    // we need a non-const list in order to be able to add sublists below!
    Teuchos::ParameterList prbdyn(elchcontrol);
    // support for turbulent flow statistics
    const Teuchos::ParameterList& fdyn = (problem->FluidDynamicParams());
    prbdyn.sublist("TURBULENCE MODEL")=fdyn.sublist("TURBULENCE MODEL");

    RefCountPtr<DRT::Discretization> aledis = problem->GetDis("ale");
    if (!aledis->Filled()) aledis->FillComplete();
    // is ALE needed or not?
    const INPAR::ELCH::ElchMovingBoundary withale
      = DRT::INPUT::IntegralValue<INPAR::ELCH::ElchMovingBoundary>(elchcontrol,"MOVINGBOUNDARY");

    if (withale!=INPAR::ELCH::elch_mov_bndry_no)
    {
      // create ale elements only if the ale discretization is empty
      if (aledis->NumGlobalNodes()==0)
      {
        Epetra_Time time(comm);
        {
          // get material cloning map
          std::map<std::pair<string,string>,std::map<int,int> > clonefieldmatmap = problem->ClonedMaterialMap();
          if (clonefieldmatmap.size() == 0)
            dserror("No CLONING MATERIAL MAP defined in input file. "
                "This is necessary to assign a material to the ALE elements.");

          std::pair<string,string> key("fluid","ale");
          std::map<int,int> fluidmatmap = clonefieldmatmap[key];
          if (fluidmatmap.size() == 0)
            dserror("Key pair 'fluid/ale' was not found in input file.");

          Teuchos::RCP<DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy> > alecreator =
            Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy>() );

          alecreator->CreateMatchingDiscretization(fluiddis,aledis,fluidmatmap);
        }

        if(comm.MyPID()==0)
        {
          cout << "\n\nCreated ALE discretization from fluid field in........"
          <<time.ElapsedTime() << " secs\n\n";
        }
      }
      else
        dserror("Providing an ALE mesh is not supported for problemtype Electrochemistry.");

      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == (-1))
        dserror("no linear solver defined for ELCH problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

      // create an ELCH::MovingBoundaryAlgorithm instance
      Teuchos::RCP<ELCH::MovingBoundaryAlgorithm> elch
        = Teuchos::rcp(new ELCH::MovingBoundaryAlgorithm(comm,prbdyn,problem->SolverParams(linsolvernumber)));

      // read the restart information, set vectors and variables
      if (restart) elch->ReadRestart(restart);

      // solve the whole electrochemistry problem
      elch->TimeLoop();

      // summarize the performance measurements
      Teuchos::TimeMonitor::summarize();

      // perform the result test
      problem->AddFieldTest(elch->FluidField().CreateFieldTest());
      problem->AddFieldTest(elch->AleField().CreateFieldTest());
      problem->AddFieldTest(elch->CreateScaTraFieldTest());
      problem->TestAll(comm);
    }
    else
    {
      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == (-1))
        dserror("no linear solver defined for ELCH problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

      // create an ELCH::Algorithm instance
      Teuchos::RCP<ELCH::Algorithm> elch = Teuchos::rcp(new ELCH::Algorithm(comm,prbdyn,problem->SolverParams(linsolvernumber)));

      // read the restart information, set vectors and variables
      if (restart) elch->ReadRestart(restart);

      // solve the whole electrochemistry problem
      elch->TimeLoop();

      // summarize the performance measurements
      Teuchos::TimeMonitor::summarize();

      // perform the result test
      problem->AddFieldTest(elch->FluidField().CreateFieldTest());
      problem->AddFieldTest(elch->CreateScaTraFieldTest());
      problem->TestAll(comm);
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


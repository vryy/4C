/*!----------------------------------------------------------------------
\file elch_dyn.cpp
\brief Control routine for Electrochemistry module.


<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>

*----------------------------------------------------------------------*/


#include "elch_dyn.H"
#include "elch_algorithm.H"
#include "elch_moving_boundary_algorithm.H"
#include "../drt_inpar/inpar_elch.H"
#include "../drt_scatra/scatra_utils_clonestrategy.H"
#include "../drt_ale/ale_utils_clonestrategy.H"
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
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  // access the scatra discretization
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");

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
    DRT::INPUT::PrintDefaultParameters(IO::cout, elchcontrol);

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
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatraonly = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(scatradyn,false,"scatra",DRT::Problem::Instance()->SolverParams(linsolvernumber)));

    // read the restart information, set vectors and variables
    if (restart) (scatraonly->ScaTraField())->ReadRestart(restart);

    // set velocity field
    // note: The order ReadRestart() before SetVelocityField() is important here!!
    // for time-dependent velocity fields, SetVelocityField() is additionally called in each PrepareTimeStep()-call
    (scatraonly->ScaTraField())->SetVelocityField();

    // enter time loop to solve problem with given convective velocity
    (scatraonly->ScaTraField())->TimeLoop();

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
      DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(fluiddis,scatradis);
    }
    else
      dserror("Fluid AND ScaTra discretization present. This is not supported.");

    // support for turbulent flow statistics
    const Teuchos::ParameterList& fdyn = (problem->FluidDynamicParams());

    Teuchos::RCP<DRT::Discretization> aledis = problem->GetDis("ale");
    if (!aledis->Filled()) aledis->FillComplete();
    // is ALE needed or not?
    const INPAR::ELCH::ElchMovingBoundary withale
      = DRT::INPUT::IntegralValue<INPAR::ELCH::ElchMovingBoundary>(elchcontrol,"MOVINGBOUNDARY");

    if (withale!=INPAR::ELCH::elch_mov_bndry_no)
    {
      // create ale elements only if the ale discretization is empty
      if (aledis->NumGlobalNodes()==0)
      {
        DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(fluiddis,aledis);
      }
      else
        dserror("Providing an ALE mesh is not supported for problemtype Electrochemistry.");

      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == (-1))
        dserror("no linear solver defined for ELCH problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

      // create an ELCH::MovingBoundaryAlgorithm instance
      Teuchos::RCP<ELCH::MovingBoundaryAlgorithm> elch
        = Teuchos::rcp(new ELCH::MovingBoundaryAlgorithm(comm,elchcontrol,scatradyn,problem->SolverParams(linsolvernumber)));

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
      Teuchos::RCP<ELCH::Algorithm> elch = Teuchos::rcp(new ELCH::Algorithm(comm,elchcontrol,scatradyn,fdyn,problem->SolverParams(linsolvernumber)));

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
  default: dserror("Unknown velocity field type for transport of passive scalar: %d",veltype); break;
  }

  return;

} // elch_dyn()


/*----------------------------------------------------------------------*/
// print ELCH-Module logo
/*----------------------------------------------------------------------*/
void printlogo()
{
    // more at http://www.ascii-art.de under entry "moose" (or "elk")
    std::cout<<"     ___            ___    "<<std::endl;
    std::cout<<"    /   \\          /   \\ "<<std::endl;
    std::cout<<"    \\_   \\        /  __/ "<<std::endl;
    std::cout<<"     _\\   \\      /  /__  "<<"     _____ _     _____  _   _   "<<std::endl;
    std::cout<<"     \\___  \\____/   __/  "<<"    |  ___| |   /  __ \\| | | |  "<<std::endl;
    std::cout<<"         \\_       _/      "<<"   | |__ | |   | /  \\/| |_| |  "<<std::endl;
    std::cout<<"           | @ @  \\_      "<<"   |  __|| |   | |    |  _  |   "<<std::endl;
    std::cout<<"           |               "<<"  | |___| |___| \\__/\\| | | | "<<std::endl;
    std::cout<<"         _/     /\\        "<<"   \\____/\\_____/\\____/\\_| |_/ "<<std::endl;
    std::cout<<"        /o)  (o/\\ \\_     "<<std::endl;
    std::cout<<"        \\_____/ /         "<<std::endl;
    std::cout<<"          \\____/          "<<std::endl;
    std::cout<<"                           "<<std::endl;
}


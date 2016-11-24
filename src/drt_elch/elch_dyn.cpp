/*!----------------------------------------------------------------------
\file elch_dyn.cpp
\brief Control routine for Electrochemistry module.

\level 2

<pre>
\maintainer Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>

*----------------------------------------------------------------------*/
#include "../drt_ale/ale_utils_clonestrategy.H"

#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/inpar_elch.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_dofset_aux_proxy.H"

#include "../drt_scatra/scatra_resulttest_elch.H"
#include "../drt_scatra/scatra_timint_elch.H"
#include "../drt_scatra/scatra_utils_clonestrategy.H"
#include "../drt_scatra_ele/scatra_ele.H"

#include <Epetra_Time.h>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "elch_algorithm.H"
#include "elch_moving_boundary_algorithm.H"
#include "elch_dyn.H"


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

    // add proxy of velocity related degrees of freedom to scatra discretization
    Teuchos::RCP<DRT::DofSetInterface> dofsetaux =
        Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(DRT::Problem::Instance()->NDim()+1, 0, 0, true));
    if ( scatradis->AddDofSet(dofsetaux)!= 1 )
      dserror("Scatra discretization has illegal number of dofsets!");


    // get linear solver id from SCALAR TRANSPORT DYNAMIC
    const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
    if (linsolvernumber == (-1))
      dserror("no linear solver defined for ELCH problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

    // create instance of scalar transport basis algorithm (empty fluid discretization)
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatraonly =
        Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());

    // now we can call Init() on the base algorithm
    // scatra time integrator is constructed and initialized inside
    scatraonly->Init(
        scatradyn,
        scatradyn,
        DRT::Problem::Instance()->SolverParams(linsolvernumber));

    // now me may redistribute or ghost the scatra discretization
    // finalize discretization
    scatradis->FillComplete(true, true, true);

    // only now we must call Setup() on the base algorithm.
    // all objects relying on the parallel distribution are
    // created and pointers are set.
    // calls Setup() on time integrator inside.
    scatraonly->Setup();

    // read the restart information, set vectors and variables
    if (restart) (scatraonly->ScaTraField())->ReadRestart(restart);

    // set velocity field
    // note: The order ReadRestart() before SetVelocityField() is important here!!
    // for time-dependent velocity fields, SetVelocityField() is additionally called in each PrepareTimeStep()-call
    (scatraonly->ScaTraField())->SetVelocityField(1);

    // enter time loop to solve problem with given convective velocity
    (scatraonly->ScaTraField())->TimeLoop();

    // perform the result test if required
    Teuchos::RCP<SCATRA::ScaTraTimIntElch> elchtimint = Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntElch>(scatraonly->ScaTraField());
    if(elchtimint == Teuchos::null)
      dserror("Time integrator is not of electrochemistry type!");
    DRT::Problem::Instance()->AddFieldTest(Teuchos::rcp(new SCATRA::ElchResultTest(elchtimint)));
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
      // fill scatra discretization by cloning fluid discretization
      DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(fluiddis,scatradis);
      scatradis->FillComplete();
      // determine implementation type of cloned scatra elements
      INPAR::SCATRA::ImplType impltype(INPAR::SCATRA::impltype_undefined);
      if(DRT::INPUT::IntegralValue<int>(elchcontrol,"DIFFCOND_FORMULATION"))
        impltype = INPAR::SCATRA::impltype_elch_diffcond;
      else
        impltype = INPAR::SCATRA::impltype_elch_NP;

      // set implementation type
      for(int i=0; i<scatradis->NumMyColElements(); ++i)
      {
        DRT::ELEMENTS::Transport* element = dynamic_cast<DRT::ELEMENTS::Transport*>(scatradis->lColElement(i));
        if(element == NULL)
          dserror("Invalid element type!");
        else
          element->SetImplType(impltype);
      }
    }

    else
      dserror("Fluid AND ScaTra discretization present. This is not supported.");

    // support for turbulent flow statistics
    const Teuchos::ParameterList& fdyn = (problem->FluidDynamicParams());

    Teuchos::RCP<DRT::Discretization> aledis = problem->GetDis("ale");
    if (!aledis->Filled()) aledis->FillComplete(false,false,false);
    // is ALE needed or not?
    const INPAR::ELCH::ElchMovingBoundary withale
      = DRT::INPUT::IntegralValue<INPAR::ELCH::ElchMovingBoundary>(elchcontrol,"MOVINGBOUNDARY");

    if (withale!=INPAR::ELCH::elch_mov_bndry_no)
    {
      // create ale elements only if the ale discretization is empty
      if (aledis->NumGlobalNodes()==0)
      {
        // clone ALE discretization from fluid discretization
        DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(fluiddis,aledis);

        aledis->FillComplete(false,true,false);
        // setup material in every ALE element
        Teuchos::ParameterList params;
        params.set<std::string>("action", "setup_material");
        aledis->Evaluate(params);
      }
      else
        dserror("Providing an ALE mesh is not supported for problemtype Electrochemistry.");

      // add proxy of fluid degrees of freedom to scatra discretization
      if(scatradis->AddDofSet(fluiddis->GetDofSetProxy()) != 1)
        dserror("Scatra discretization has illegal number of dofsets!");

      // add proxy of ALE degrees of freedom to scatra discretization
      if(scatradis->AddDofSet(aledis->GetDofSetProxy()) != 2)
        dserror("Scatra discretization has illegal number of dofsets!");

      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == (-1))
        dserror("no linear solver defined for ELCH problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

      // create an ELCH::MovingBoundaryAlgorithm instance
      Teuchos::RCP<ELCH::MovingBoundaryAlgorithm> elch
        = Teuchos::rcp(new ELCH::MovingBoundaryAlgorithm(comm,elchcontrol,scatradyn,problem->SolverParams(linsolvernumber)));

      // now we must call Init()
      // NOTE : elch reads time parameters from scatra dynamic section !
      elch->Init(
          scatradyn,
          scatradyn,
          problem->SolverParams(linsolvernumber) );

      // NOTE : At this point we may redistribute and/or
      //        ghost our discretizations at will.
      scatradis->FillComplete();
      fluiddis->FillComplete();
      aledis->FillComplete();

      // now we can call Setup() on the scatra time integrator
      elch->Setup();

      // read the restart information, set vectors and variables
      if (restart) elch->ReadRestart(restart);

      // solve the whole electrochemistry problem
      elch->TimeLoop();

      // summarize the performance measurements
      Teuchos::TimeMonitor::summarize();

      // perform the result test
      problem->AddFieldTest(elch->FluidField()->CreateFieldTest());
      problem->AddFieldTest(elch->AleField()->CreateFieldTest());
      Teuchos::RCP<SCATRA::ScaTraTimIntElch> elchtimint = Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntElch>(elch->ScaTraField());
      if(elchtimint == Teuchos::null)
        dserror("Time integrator is not of electrochemistry type!");
      problem->AddFieldTest(Teuchos::rcp(new SCATRA::ElchResultTest(elchtimint)));
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

      // add proxy of fluid degrees of freedom to scatra discretization
      if(scatradis->AddDofSet(fluiddis->GetDofSetProxy()) != 1)
        dserror("Scatra discretization has illegal number of dofsets!");

      // now we must call Init()
      // NOTE : elch reads time parameters from scatra dynamic section !
      elch->Init(
          scatradyn,
          scatradyn,
          problem->SolverParams(linsolvernumber));

      // NOTE : At this point we may redistribute and/or
      //        ghost our discretizations at will.
      scatradis->FillComplete();
      fluiddis->FillComplete();
      aledis->FillComplete();

      // discretizations are done, now we can call Setup() on the algorithm
      elch->Setup();


      // read the restart information, set vectors and variables
      if (restart) elch->ReadRestart(restart);

      // solve the whole electrochemistry problem
      elch->TimeLoop();

      // summarize the performance measurements
      Teuchos::TimeMonitor::summarize();

      // perform the result test
      problem->AddFieldTest(elch->FluidField()->CreateFieldTest());
      Teuchos::RCP<SCATRA::ScaTraTimIntElch> elchtimint = Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntElch>(elch->ScaTraField());
      if(elchtimint == Teuchos::null)
        dserror("Time integrator is not of electrochemistry type!");
      problem->AddFieldTest(Teuchos::rcp(new SCATRA::ElchResultTest(elchtimint)));
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


/*!----------------------------------------------------------------------
\file loma_dyn.cpp

\brief Control routine for low-Mach-number flow module.

\level 2

\maintainer Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915245


*----------------------------------------------------------------------*/


#include <string>
#include <iostream>

#include "loma_dyn.H"
#include "loma_algorithm.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_utils_clonestrategy.H"
#include "../drt_scatra_ele/scatra_ele.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"
#include <Epetra_Time.h>
#include "../drt_lib/drt_dofset_predefineddofnumber.H"


/*----------------------------------------------------------------------*/
// entry point for LOMA in DRT
/*----------------------------------------------------------------------*/
void loma_dyn(int restart)
{
  // create a communicator
#ifdef PARALLEL
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("fluid")->Comm();
#else
  Epetra_SerialComm comm;
#endif

  // print warning to screen
  if (comm.MyPID()==0)
    std::cout << "You are now about to enter the module for low-Mach-number flow!" <<std::endl;

  // define abbreviation
  DRT::Problem* problem = DRT::Problem::Instance();

  // access fluid and (typically empty) scatra discretization
  Teuchos::RCP<DRT::Discretization> fluiddis  = problem->GetDis("fluid");
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");

  // ensure that all dofs are assigned in the right order such that
  // dof numbers are created with fluid dof < scatra/elch dof
  fluiddis->FillComplete();
  scatradis->FillComplete();

  // access problem-specific parameter list for LOMA
  const Teuchos::ParameterList& lomacontrol = problem->LOMAControlParams();

  // access parameter list for scatra
  const Teuchos::ParameterList& scatradyn = problem->ScalarTransportDynamicParams();

  // access parameter list for fluid
  const Teuchos::ParameterList& fdyn = problem->FluidDynamicParams();

  // identify type of velocity field
  const INPAR::SCATRA::VelocityField veltype = DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatradyn,"VELOCITYFIELD");

  // choose algorithm depending on type of velocity field
  switch (veltype)
  {
  case INPAR::SCATRA::velocity_zero:  // zero velocity field (see case 1)
  case INPAR::SCATRA::velocity_function:  // velocity field prescribed by function
  case INPAR::SCATRA::velocity_function_and_curve: // velocity field prescribed by function and time curve
  {
    // directly use elements from input section 'transport elements'
    if (scatradis->NumGlobalNodes()==0)
      dserror("No elements in input section ---TRANSPORT ELEMENTS!");

    // get linear solver id from SCALAR TRANSPORT DYNAMIC
    const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
    if (linsolvernumber == (-1))
      dserror("no linear solver defined for LOMA problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

    // create instance of scalar transport basis algorithm (no fluid discretization)
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatraonly =
        Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());

    // add proxy of velocity related degrees of freedom to scatra discretization
    Teuchos::RCP<DRT::DofSetInterface> dofsetaux =
        Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(DRT::Problem::Instance()->NDim()+1, 0, 0, true));
    if ( scatradis->AddDofSet(dofsetaux)!= 1 )
      dserror("Scatra discretization has illegal number of dofsets!");

    // now we can call Init() on base algo
    scatraonly->Init(
        lomacontrol,
        scatradyn,
        DRT::Problem::Instance()->SolverParams(linsolvernumber));

    // only now we must call Setup() on the scatra time integrator.
    // all objects relying on the parallel distribution are
    // created and pointers are set.
    // calls Setup() on the time integrator inside
    scatraonly->Setup();

    // read restart information
    if (restart) (scatraonly->ScaTraField())->ReadRestart(restart);

    // set initial velocity field
    // note: The order ReadRestart() before SetVelocityField() is important here!!
    // for time-dependent velocity fields, SetVelocityField() is additionally called in each PrepareTimeStep()-call
    (scatraonly->ScaTraField())->SetVelocityField(1);

    // enter time loop to solve problem with given convective velocity field
    (scatraonly->ScaTraField())->TimeLoop();

    // perform result test if required
    problem->AddFieldTest(scatraonly->CreateScaTraFieldTest());
    problem->TestAll(comm);

    break;
  }
  case INPAR::SCATRA::velocity_Navier_Stokes:  // Navier_Stokes
  {
    // use fluid discretization as layout for scatra discretization
    if (fluiddis->NumGlobalNodes()==0) dserror("Fluid discretization is empty!");

    // to generate turbulent flow in the inflow section only, it is not necessary to
    // solve the transport equation for the temperature
    // therefore, use problem type fluid
    if ((DRT::INPUT::IntegralValue<int>(fdyn.sublist("TURBULENT INFLOW"),"TURBULENTINFLOW")==true) and
       (restart<fdyn.sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP")))
      dserror("Choose problem type fluid to generate turbulent flow in the inflow section!");

    // create scatra elements if scatra discretization is empty (typical case)
    if (scatradis->NumGlobalNodes()==0)
    {
      // fill scatra discretization by cloning fluid discretization
      DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(fluiddis,scatradis);

      // set implementation type of cloned scatra elements to loma
      for(int i=0; i<scatradis->NumMyColElements(); ++i)
      {
        DRT::ELEMENTS::Transport* element = dynamic_cast<DRT::ELEMENTS::Transport*>(scatradis->lColElement(i));
        if(element == NULL)
          dserror("Invalid element type!");
        else
          element->SetImplType(INPAR::SCATRA::impltype_loma);
      }
    }
    else dserror("Fluid AND ScaTra discretization present. This is not supported.");

    // get linear solver id from SCALAR TRANSPORT DYNAMIC
    const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
    if (linsolvernumber == (-1))
      dserror("no linear solver defined for LOMA problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

    // create a LOMA::Algorithm instance
    Teuchos::RCP<LOMA::Algorithm> loma = Teuchos::rcp(new LOMA::Algorithm(comm,lomacontrol,DRT::Problem::Instance()->SolverParams(linsolvernumber)));

    // add proxy of fluid transport degrees of freedom to scatra discretization
    if(scatradis->AddDofSet(fluiddis->GetDofSetProxy()) != 1)
      dserror("Scatra discretization has illegal number of dofsets!");

    loma->Init(
        lomacontrol,
        DRT::Problem::Instance()->ScalarTransportDynamicParams(),
        DRT::Problem::Instance()->SolverParams(linsolvernumber) );

    loma->Setup();

    // read restart information
    // in case a inflow generation in the inflow section has been performed, there are not any
    // scatra results available and the initial field is used
    if (restart)
    {
      if ((DRT::INPUT::IntegralValue<int>(fdyn.sublist("TURBULENT INFLOW"),"TURBULENTINFLOW")==true) and
         (restart==fdyn.sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP")))
        loma->ReadInflowRestart(restart);
      else
        loma->ReadRestart(restart);
    }

    // enter LOMA algorithm
    loma->TimeLoop();

    // summarize performance measurements
    Teuchos::TimeMonitor::summarize();

    // perform result test if required
    problem->AddFieldTest(loma->FluidField()->CreateFieldTest());
    problem->AddFieldTest(loma->CreateScaTraFieldTest());
    problem->TestAll(comm);

    break;
  }
  default:
    dserror("Unknown velocity field type for low-Mach-number flow: %d",veltype);
    break;
  }

  return;

} // loma_dyn()



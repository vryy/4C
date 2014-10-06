/*!----------------------------------------------------------------------
\file two_phase_dyn.cpp
\brief Control routine for fluid/xfluid and ScaTra coupled routines.


<pre>
Maintainer: Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915236
</pre>

*----------------------------------------------------------------------*/


#include <string>
#include <iostream>

#include "two_phase_algorithm.H"
#include "../drt_fluid_xfluid/xfluid_levelset_coupling_algorithm.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_scatra/scatra_utils_clonestrategy.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"
#include <Epetra_Time.h>

#include "../drt_lib/drt_dofset_fixed_size.H"
#include "two_phase_dyn.H"

/*----------------------------------------------------------------------*/
// entry point for Two Phase Flow (TPF) in DRT
/*----------------------------------------------------------------------*/
void two_phase_dyn(int restart)
{
  // create a communicator
#ifdef PARALLEL
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("fluid")->Comm();
#else
  Epetra_SerialComm comm;
#endif

  // print warning to screen
  if (comm.MyPID()==0)
  {
    std::cout << "                   T(wo) P(hase) F(low)                  " <<std::endl;
    std::cout << "=========================================================" <<std::endl;
    std::cout << "You are now about to enter the module for two phase flow!" <<std::endl;
    std::cout << "=========================================================" <<std::endl;
  }


  // define abbreviation
  DRT::Problem* problem = DRT::Problem::Instance();

  // access fluid and (typically empty) scatra discretization
  Teuchos::RCP<DRT::Discretization> fluiddis  = problem->GetDis("fluid");
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");

   // ensure that all dofs are assigned in the right order such that
  // dof numbers are created with fluid dof < scatra/elch dof
  fluiddis->FillComplete();
  scatradis->FillComplete();


  //access parameter for two phase flow
  const Teuchos::ParameterList& twophaseflowcontrol = problem->TwoPhaseFlowParams();

  // access parameter for levelset (Not needed as of yet)
  //const Teuchos::ParameterList& levelsetcontrol = problem->LevelSetControl();

  // access parameter list for scatra
  const Teuchos::ParameterList& scatradyn = problem->ScalarTransportDynamicParams();

  // access parameter list for fluid
  const Teuchos::ParameterList& fdyn = problem->FluidDynamicParams();

  // identify type of velocity field
  const INPAR::SCATRA::VelocityField veltype = DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatradyn,"VELOCITYFIELD");

  // choose algorithm depending on type of velocity field
  switch (veltype)
  {
  case INPAR::SCATRA::velocity_Navier_Stokes:  // Navier_Stokes
  {
    // use fluid discretization as layout for scatra discretization
    if (fluiddis->NumGlobalNodes()==0) dserror("Fluid discretization is empty!");

    // to generate turbulent flow in the inflow section only, it is not necessary to
    // solve the levelset equation
    // therefore, use problem type fluid
    if ((DRT::INPUT::IntegralValue<int>(fdyn.sublist("TURBULENT INFLOW"),"TURBULENTINFLOW")==true) and
       (restart<fdyn.sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP")))
      dserror("Choose problem type fluid to generate turbulent flow in the inflow section!");

    // create scatra elements if scatra discretization is empty (typical case)
    if (scatradis->NumGlobalNodes()==0)
    {
      DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(fluiddis,scatradis);
    }
    else dserror("Fluid AND ScaTra discretization present. This is not supported.");

    // get linear solver id from SCALAR TRANSPORT DYNAMIC
    const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
    if (linsolvernumber == (-1))
      dserror("no linear solver defined for two phase flow (TPF) problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

    // create a TWOPHASEFLOW::Algorithm instance
    Teuchos::RCP<TWOPHASEFLOW::Algorithm> twophase = Teuchos::rcp(new TWOPHASEFLOW::Algorithm(comm,twophaseflowcontrol,DRT::Problem::Instance()->SolverParams(linsolvernumber)));

    // read restart information
    // in case an inflow generation in the inflow section has been performed, there are not any
    // scatra results available and the initial field is used
    if (restart)
    {
      // turn on/off read scatra restart from input file
      const bool restartscatrainput = (bool)DRT::INPUT::IntegralValue<int>(twophaseflowcontrol,"RESTART_SCATRA_INPUT");

      if ((DRT::INPUT::IntegralValue<int>(fdyn.sublist("TURBULENT INFLOW"),"TURBULENTINFLOW")==true) and
         (restart==fdyn.sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP")))
      {
        std::cout << "/!\\ warning === Restart with turbulent inflow has not been tested for TPF problems. Proceed with caution." << std::endl;
        //This function is untested! Might need to be modified when turbulent inflow is needed.
        twophase->ReadInflowRestart(restart);
      }
      else
      {
        // read the restart information, set vectors and variables
        twophase->Restart(restart, restartscatrainput);
      }
    }

    INPAR::FLUID::TimeIntegrationScheme timeintscheme = DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fdyn,"TIMEINTEGR");


    if (timeintscheme == INPAR::FLUID::timeint_one_step_theta or
        timeintscheme == INPAR::FLUID::timeint_afgenalpha or
        timeintscheme == INPAR::FLUID::timeint_bdf2)
    {
      // solve the two phase problem utilizing the smoothing function for parameter values.
      twophase->TimeLoop();
    }
    else if (timeintscheme == INPAR::FLUID::timeint_stationary)
    {
      // solve a stationary two phase problem utilizing the smoothing function for parameter values.
      twophase->SolveStationaryProblem();
    }
    else
    {
      // every time integration scheme must be either static or dynamic
      dserror("the two phase module can not handle this type of time integration scheme");
    }

    //------------------------------------------------------------------------------------------------
    // validate the results
    //------------------------------------------------------------------------------------------------
    // summarize the performance measurements
    Teuchos::TimeMonitor::summarize();


    // perform the result test
    twophase->TestResults();

    break;
  }
  default:
    dserror("Unknown velocity field type for two phase flow: %d",veltype);
    break;
  }


  return;

} // two_phase_dyn()


void fluid_xfem_ls_drt()
{
  // create a communicator
#ifdef PARALLEL
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("fluid")->Comm();
#else
  Epetra_SerialComm comm;
#endif

  // print warning to screen
  if (comm.MyPID()==0)
  {
    std::cout << "=========================================================" <<std::endl;
    std::cout << "|                  XFluid with levelset                 |" <<std::endl;
    std::cout << "=========================================================" <<std::endl;
    std::cout << "|   Cut is done with level set. Calculations in Xfluid  | " <<std::endl;
    std::cout << "=========================================================" <<std::endl;
    std::cout << "|          XFEM is utilized for the computations        | " <<std::endl;
    std::cout << "=========================================================" <<std::endl;
  }

  // define abbreviation
   DRT::Problem* problem = DRT::Problem::Instance();

   // access fluid and (typically empty) scatra discretization
   Teuchos::RCP<DRT::Discretization> fluiddis  = problem->GetDis("fluid");
   Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");

   fluiddis->FillComplete();

//   access parameter for two phase flow
   const Teuchos::ParameterList& xdyn = problem->XFEMGeneralParams();

   // Reserve DoF's for fluid
   int numglobalnodes = fluiddis->NumGlobalNodes();
   int maxNumMyReservedDofs = numglobalnodes*(xdyn.get<int>("MAX_NUM_DOFSETS"))*4;
   Teuchos::RCP<DRT::FixedSizeDofSet> maxdofset = Teuchos::rcp(new DRT::FixedSizeDofSet(maxNumMyReservedDofs));
   fluiddis->ReplaceDofSet(maxdofset,true);
   fluiddis->FillComplete();

   // access parameter for levelset (Not needed as of yet)
   //const Teuchos::ParameterList& levelsetcontrol = problem->LevelSetControl();

   // access parameter list for scatra
   const Teuchos::ParameterList& scatradyn = problem->ScalarTransportDynamicParams();

   // access parameter list for fluid
   const Teuchos::ParameterList& fdyn = problem->FluidDynamicParams();

   // use fluid discretization as layout for scatra discretization
   if (fluiddis->NumGlobalNodes()==0) dserror("Fluid discretization is empty!");

   // create scatra elements if scatra discretization is empty (typical case)
   if (scatradis->NumGlobalNodes()==0)
   {
     DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(fluiddis,scatradis);

     // Give ScaTra new dofset (starts after fluid)
     Teuchos::RCP<DRT::DofSet> newdofset = Teuchos::rcp(new DRT::DofSet());
     scatradis->ReplaceDofSet(newdofset,true);
     scatradis->FillComplete();

     // print all dofsets
     //---FLUID---|---SCATRA---|
     fluiddis->GetDofSetProxy()->PrintAllDofsets(fluiddis->Comm());

   }
   else dserror("Fluid AND ScaTra discretization present. This is not supported.");

   // get linear solver id from SCALAR TRANSPORT DYNAMIC
   const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
   if (linsolvernumber == (-1))
     dserror("no linear solver defined for two phase flow (TPF) problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

   Teuchos::RCP<XFLUIDLEVELSET::Algorithm>  xfluid_levelset = Teuchos::rcp(new XFLUIDLEVELSET::Algorithm(comm,fdyn,DRT::Problem::Instance()->SolverParams(linsolvernumber)));

   INPAR::FLUID::TimeIntegrationScheme timeintscheme = DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fdyn,"TIMEINTEGR");

   if (timeintscheme == INPAR::FLUID::timeint_stationary)
   {
     xfluid_levelset->SolveStationaryProblem();
   }
   else
   {
     // every time integration scheme must be either static or dynamic
     dserror("Only stationary time integration is currently implemented for Fluid_XFEM_LevelSet!");
   }

   //------------------------------------------------------------------------------------------------
   // validate the results
   //------------------------------------------------------------------------------------------------
   // summarize the performance measurements (already done in Xfluid SolveStationaryProblem())
   //Teuchos::TimeMonitor::summarize();


   // perform the result test
   xfluid_levelset->TestResults();


   return;

} // fluid_xfem_ls_drt()

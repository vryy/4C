/*!----------------------------------------------------------------------
\file two_phase_dyn.cpp
\brief Control routine for fluid/xfluid and ScaTra coupled routines.
\level 2

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915236
</pre>

*----------------------------------------------------------------------*/
#include "../drt_fluid_xfluid/xfluid_levelset_coupling_algorithm.H"

#include "../drt_inpar/drt_validparameters.H"

#include "../drt_lib/drt_discret_xfem.H"
#include "../drt_lib/drt_dofset_fixed_size.H"
#include "../drt_lib/drt_dofset_proxy.H"
#include "../drt_lib/drt_dofset_interface.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_scatra/scatra_utils_clonestrategy.H"

#include "../drt_scatra_ele/scatra_ele.H"

#include <Epetra_Time.h>
#include <iostream>
#include <string>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "two_phase_algorithm.H"
#include "two_phase_dyn.H"

#include "../drt_lib/drt_dofset_predefineddofnumber.H"

#include "../drt_adapter/ad_fld_fluid.H"
#include "../drt_fluid_xfluid/xfluid.H"



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
  if (comm.MyPID() == 0)
  {
    //    std::cout << "                   T(wo) P(hase) F(low)                  " <<std::endl;
    //    std::cout << "=========================================================" <<std::endl;
    //    std::cout << "You are now about to enter the module for two phase flow!" <<std::endl;
    //    std::cout << "=========================================================" <<std::endl;

    //    std::cout <<"                                      \n"
    //              <<"           T(wo) P(hase) F(low)       \n"
    //              <<"                __ \\ / __            \n"
    //              <<"               /  \\ | /  \\          \n"
    //              <<"                   \\|/               \n"
    //              <<"              _,.---v---._            \n"
    //              <<"     /\\__/\\  /            \\        \n"
    //              <<"     \\_  _/ /              \\        \n"
    //              <<"       \\ \\_|           @ __|        \n"
    //              <<"        \\                \\_         \n"
    //              <<"         \\     ,__/       /          \n"
    //              <<" ~~~~~~~~~`~~~~~~~~~~~~~~/~~~~~~~ " << std::endl;


    std::cout << "         .-``'.     T(wo) P(hase) F(low)   .'''-.          " << std::endl;
    std::cout << "       .`   .`~        SMEARED (CSF)       ~`.   '.        " << std::endl;
    std::cout << "   _.-'     '._                            _.'     '-._    " << std::endl;
    std::cout << "                                                           " << std::endl;
  }


  // define abbreviation
  DRT::Problem* problem = DRT::Problem::Instance();

  // access fluid and (typically empty) scatra discretization
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");

  // ensure that all dofs are assigned in the right order such that
  // dof numbers are created with fluid dof < scatra/elch dof
  fluiddis->FillComplete();
  scatradis->FillComplete();


  // access parameter for two phase flow
  const Teuchos::ParameterList& twophaseflowcontrol = problem->TwoPhaseFlowParams();

  // access parameter for levelset (Not needed as of yet)
  // const Teuchos::ParameterList& levelsetcontrol = problem->LevelSetControl();

  // access parameter list for scatra
  const Teuchos::ParameterList& scatradyn = problem->ScalarTransportDynamicParams();

  // access parameter list for fluid
  const Teuchos::ParameterList& fdyn = problem->FluidDynamicParams();

  // identify type of velocity field
  const INPAR::SCATRA::VelocityField veltype =
      DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatradyn, "VELOCITYFIELD");

  // choose algorithm depending on type of velocity field
  switch (veltype)
  {
    case INPAR::SCATRA::velocity_Navier_Stokes:  // Navier_Stokes
    {
      // use fluid discretization as layout for scatra discretization
      if (fluiddis->NumGlobalNodes() == 0) dserror("Fluid discretization is empty!");

      // to generate turbulent flow in the inflow section only, it is not necessary to
      // solve the levelset equation
      // therefore, use problem type fluid
      if ((DRT::INPUT::IntegralValue<int>(fdyn.sublist("TURBULENT INFLOW"), "TURBULENTINFLOW") ==
              true) and
          (restart < fdyn.sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP")))
        dserror("Choose problem type fluid to generate turbulent flow in the inflow section!");

      // create scatra elements if scatra discretization is empty (typical case)
      if (scatradis->NumGlobalNodes() == 0)
      {
        // fill scatra discretization by cloning fluid discretization
        DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(fluiddis, scatradis);

        // set implementation type of cloned scatra elements to levelset
        for (int i = 0; i < scatradis->NumMyColElements(); ++i)
        {
          DRT::ELEMENTS::Transport* element =
              dynamic_cast<DRT::ELEMENTS::Transport*>(scatradis->lColElement(i));
          if (element == NULL)
            dserror("Invalid element type!");
          else
            element->SetImplType(INPAR::SCATRA::impltype_levelset);
        }
      }
      else
        dserror("Fluid AND ScaTra discretization present. This is not supported.");

      // add proxy of fluid degrees of freedom to scatra discretization
      if (scatradis->AddDofSet(fluiddis->GetDofSetProxy()) != 1)
        dserror("Scatra discretization has illegal number of dofsets!");

      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == (-1))
        dserror(
            "no linear solver defined for two phase flow (TPF) problem. Please set LINEAR_SOLVER "
            "in SCALAR TRANSPORT DYNAMIC to a valid number!");

      // create a TWOPHASEFLOW::Algorithm instance
      Teuchos::RCP<TWOPHASEFLOW::Algorithm> twophase = Teuchos::rcp(new TWOPHASEFLOW::Algorithm(
          comm, twophaseflowcontrol, DRT::Problem::Instance()->SolverParams(linsolvernumber)));
      twophase->Init(twophaseflowcontrol, DRT::Problem::Instance()->ScalarTransportDynamicParams(),
          DRT::Problem::Instance()->SolverParams(linsolvernumber));
      twophase->Setup();

      // read restart information
      // in case an inflow generation in the inflow section has been performed, there are not any
      // scatra results available and the initial field is used
      if (restart)
      {
        // turn on/off read scatra restart from input file
        const bool restartscatrainput =
            (bool)DRT::INPUT::IntegralValue<int>(twophaseflowcontrol, "RESTART_SCATRA_INPUT");

        if ((DRT::INPUT::IntegralValue<int>(fdyn.sublist("TURBULENT INFLOW"), "TURBULENTINFLOW") ==
                true) and
            (restart == fdyn.sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP")))
        {
          std::cout << "/!\\ warning === Restart with turbulent inflow has not been tested for TPF "
                       "problems. Proceed with caution."
                    << std::endl;
          // This function is untested! Might need to be modified when turbulent inflow is needed.
          twophase->ReadInflowRestart(restart);
        }
        else
        {
          // read the restart information, set vectors and variables
          twophase->Restart(restart, restartscatrainput);
        }
      }

      INPAR::FLUID::TimeIntegrationScheme timeintscheme =
          DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fdyn, "TIMEINTEGR");


      if (timeintscheme == INPAR::FLUID::timeint_one_step_theta or
          timeintscheme == INPAR::FLUID::timeint_afgenalpha or
          timeintscheme == INPAR::FLUID::timeint_bdf2)
      {
        // solve the two phase problem utilizing the smoothing function for parameter values.
        twophase->TimeLoop();
      }
      else if (timeintscheme == INPAR::FLUID::timeint_stationary)
      {
        // solve a stationary two phase problem utilizing the smoothing function for parameter
        // values.
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
      dserror("Unknown velocity field type for two phase flow: %d", veltype);
      break;
  }


  return;

}  // two_phase_dyn()


// void fluid_xfem_ls_drt(int restart)
//{
//  // create a communicator
//#ifdef PARALLEL
//  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("fluid")->Comm();
//#else
//  Epetra_SerialComm comm;
//#endif
//
//  // print warning to screen
//  if (comm.MyPID()==0)
//  {
//    //    std::cout << "=========================================================" <<std::endl;
//    //    std::cout << "|                  XFluid with levelset                 |" <<std::endl;
//    //    std::cout << "=========================================================" <<std::endl;
//    //    std::cout << "|   Cut is done with level set. Calculations in Xfluid  | " <<std::endl;
//    //    std::cout << "=========================================================" <<std::endl;
//    //    std::cout << "|          XFEM is utilized for the computations        | " <<std::endl;
//    //    std::cout << "=========================================================" <<std::endl;
//        // A whale when the time is right!
//        std::cout <<"                                        \n"
//                  <<"             T(wo) P(hase) F(low)       \n"
//                  <<"                  __ \\ / __            \n"
//                  <<"                 /  \\ | /  \\          \n"
//                  <<"                     \\|/               \n"
//                  <<"                _,.---v---._            \n"
//                  <<"       /\\__/\\  /            \\        \n"
//                  <<"       \\_  _/ /    X F E M   \\        \n"
//                  <<"         \\ \\_|           @ __|        \n"
//                  <<"          \\                \\_         \n"
//                  <<"           \\     ,__/       /          \n"
//                  <<"   ~~~~~~~~~`~~~~~~~~~~~~~~/~~~~~~~~~~  \n"
//                  <<"                                        "<< std::endl;
//  }
//
//  // define abbreviation
//   DRT::Problem* problem = DRT::Problem::Instance();
//
//   // access fluid and (typically empty) scatra discretization
//   Teuchos::RCP<DRT::DiscretizationXFEM> fluiddis =
//       Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(problem->GetDis("fluid"), true);
//   fluiddis->FillComplete();
//
//   Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");
//
//   //access parameter for two phase flow
//   // these parameters are controling the coupling algorithm.
//   const Teuchos::ParameterList& twophasedyn = problem->TwoPhaseFlowParams();
//
////   access parameter for XFEM
//   const Teuchos::ParameterList& xdyn = problem->XFEMGeneralParams();
//
//   // Reserve DoF's for fluid
//   int numglobalnodes = fluiddis->NumGlobalNodes();
//   int maxNumMyReservedDofsperNode = (xdyn.get<int>("MAX_NUM_DOFSETS"))*4;
//   Teuchos::RCP<DRT::FixedSizeDofSet> maxdofset = Teuchos::rcp(new
//   DRT::FixedSizeDofSet(maxNumMyReservedDofsperNode,numglobalnodes)); fluiddis->ReplaceDofSet(0,
//   maxdofset, true ); //fluid dofset has nds = 0 fluiddis->InitialFillComplete();
//
//   // access parameter for levelset (Not needed as of yet)
//   //const Teuchos::ParameterList& levelsetcontrol = problem->LevelSetControl();
//
//   // access parameter list for scatra
//   const Teuchos::ParameterList& scatradyn = problem->ScalarTransportDynamicParams();
//
//   // access parameter list for fluid
//   //const Teuchos::ParameterList& fdyn = problem->FluidDynamicParams();
//
//   // Needed for altgeo-creation
//   if(!scatradis->Filled())
//     scatradis->FillComplete();
//
//   // use fluid discretization as layout for scatra discretization
//   if (fluiddis->NumGlobalNodes()==0) dserror("Fluid discretization is empty!");
//
//   // create scatra elements if scatra discretization is empty (typical case)
//   if (scatradis->NumGlobalNodes()==0)
//   {
//     // fill scatra discretization by cloning fluid discretization
//     DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(fluiddis,scatradis);
//
//     // Give ScaTra new dofset (starts after fluid)
//     Teuchos::RCP<DRT::DofSet> newdofset = Teuchos::rcp(new DRT::DofSet());
//     scatradis->ReplaceDofSet(newdofset,true);
//     scatradis->FillComplete();
//
//     // set implementation type of cloned scatra elements to levelset
//     for(int i=0; i<scatradis->NumMyColElements(); ++i)
//     {
//       DRT::ELEMENTS::Transport* element =
//       dynamic_cast<DRT::ELEMENTS::Transport*>(scatradis->lColElement(i)); if(element == NULL)
//         dserror("Invalid element type!");
//       else
//         element->SetImplType(INPAR::SCATRA::impltype_levelset);
//     }
//
//     // add proxy of fluid degrees of freedom to scatra discretization
//     if(scatradis->AddDofSet(Teuchos::rcp(new
//     DRT::DofSet(Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(fluiddis)->InitialDofSet())))
//     != 1)
//       dserror("Scatra discretization has illegal number of dofsets!");
//
//     // print all dofsets
//     //---FLUID---|---SCATRA---|
//     fluiddis->GetDofSetProxy()->PrintAllDofsets(fluiddis->Comm());
//   }
//   else dserror("Fluid AND ScaTra discretization present. This is not supported.");
//
//   // ---------------------------------------------
//   //
//   //              SAFETY CHECKS!
//   //   make sure input is not conflicting
//   //
//   // ---------------------------------------------
//
//   // get linear solver id from SCALAR TRANSPORT DYNAMIC
//   const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
//   if (linsolvernumber == (-1))
//     dserror("no linear solver defined for xfem two phase flow (XTPF) problem. Please set
//     LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");
//
//   // Check for coupling condition.
//   {
//     std::string condition_to_check="XFEMLevelsetTwophase";
//
//     std::vector<std::string> names;
//     fluiddis->GetConditionNames( names );
//
//     if(std::find(names.begin(), names.end(), condition_to_check) != names.end())
//     {
//       int levelsetfunctnumberxfem   =
//       fluiddis->GetCondition(condition_to_check)->GetInt("levelsetfieldno"); int
//       levelsetfunctnumberscatra = scatradyn.get<int>("INITFUNCNO"); if(levelsetfunctnumberxfem !=
//       levelsetfunctnumberscatra)
//         dserror("Function number for level-set in SCALAR TRANSPORT DYNAMIC is not same as
//         provided in DESIGN XFEM LEVELSET TWOPHASE VOL CONDITIONS!");
//     }
//     else
//       dserror("No DESIGN XFEM LEVELSET TWOPHASE VOL CONDITIONS provided, please add to input
//       file.");
//   }
//
//   // ---------------------------------------------- //Safety check done.
//
//
//   // Test replacing fdyn in Algorithm with prbdyn
//   Teuchos::RCP<XFLUIDLEVELSET::Algorithm>  xfluid_levelset = Teuchos::rcp(new
//   XFLUIDLEVELSET::Algorithm(comm,twophasedyn,DRT::Problem::Instance()->SolverParams(linsolvernumber)));
//   xfluid_levelset->Init(
//       twophasedyn,
//       problem->ScalarTransportDynamicParams(),
//       DRT::Problem::Instance()->SolverParams(linsolvernumber));
//   xfluid_levelset->Setup();
//
//   // read restart information
//   // in case an inflow generation in the inflow section has been performed, there are not any
//   // scatra results available and the initial field is used
//   // FOR NOW ONLY RESTART WITH FLUID AND SCATRA FROM OUTPUT!!!!
//   if (restart)
//   {
//     xfluid_levelset->Restart(restart);
//   }
//
//   INPAR::FLUID::TimeIntegrationScheme timeintscheme =
//   DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(scatradyn,"TIMEINTEGR");
//
//   if (timeintscheme == INPAR::FLUID::timeint_one_step_theta)
//   {
//     // solve the two phase problem utilizing the smoothing function for parameter values.
//     xfluid_levelset->TimeLoop();
//   }
//   else if (timeintscheme == INPAR::FLUID::timeint_stationary)
//   {
//     xfluid_levelset->SolveStationaryProblem();
//   }
//   else
//   {
//     // every time integration scheme must be either static or dynamic
//     dserror("Time integration schemes not supported. Only OST or Stationary allowed for XFEM.");
//   }
//
//   //------------------------------------------------------------------------------------------------
//   // validate the results
//   //------------------------------------------------------------------------------------------------
//   // summarize the performance measurements (already done in Xfluid SolveStationaryProblem())
//   //Teuchos::TimeMonitor::summarize();
//
//
//   // perform the result test
//   xfluid_levelset->TestResults();
//
//
//   return;
//
//} // fluid_xfem_ls_drt()

// TODO: combine this with scatra_dyn and use the ScatraAlgorithm instead of the
// XFLUIDLEVELSET::Algorithm
// deep cleanup and refactoring required!!!
void fluid_xfem_ls_drt(int restart)
{
  // access the communicator
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("fluid")->Comm();

  // print warning to screen
  if (comm.MyPID() == 0)
  {
    // A whale when the time is right!
    std::cout << "                                        \n"
              << "             T(wo) P(hase) F(low)       \n"
              << "                  __ \\ / __            \n"
              << "                 /  \\ | /  \\          \n"
              << "                     \\|/               \n"
              << "                _,.---v---._            \n"
              << "       /\\__/\\  /            \\        \n"
              << "       \\_  _/ /    X F E M   \\        \n"
              << "         \\ \\_|           @ __|        \n"
              << "          \\                \\_         \n"
              << "           \\     ,__/       /          \n"
              << "   ~~~~~~~~~`~~~~~~~~~~~~~~/~~~~~~~~~~  \n"
              << "                                        " << std::endl;
  }


  // TODO: this might get a container class later (with enums as key and safety checks for access!)
  std::map<std::string, int>
      dofset_coupling_map;  ///< key=meaning of a dofset in a dis, nodal dofset number

  dofset_coupling_map.insert(std::pair<std::string, int>("vel_in_xfluid", 0));
  dofset_coupling_map.insert(std::pair<std::string, int>("phi_in_scatra", 0));

  const int xfluid_dofset = dofset_coupling_map["vel_in_xfluid"];
  const int cutter_nds_phi = dofset_coupling_map["phi_in_scatra"];


  // define abbreviation
  DRT::Problem* problem = DRT::Problem::Instance();

  // access the solid discretization
  Teuchos::RCP<DRT::Discretization> soliddis = problem->GetDis("structure");
  // access the fluid discretization
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  // access the scatra discretization
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");

  //--------------------------------------------------------------------
  // prepare solid
  soliddis->FillComplete();

  //--------------------------------------------------------------------
  // dofset set for fluid and initial fillcomplete

  // access fluid and (typically empty) scatra discretization
  Teuchos::RCP<DRT::DiscretizationXFEM> xfluiddis =
      Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(problem->GetDis("fluid"), true);
  xfluiddis->FillComplete();

  // access parameter for XFEM
  const Teuchos::ParameterList& xdyn = problem->XFEMGeneralParams();

  // Reserve DoF's for fluid
  const Epetra_Map* noderowmap = xfluiddis->NodeRowMap();
  if (noderowmap == NULL) dserror("we expect a fill-complete call before!");
  int nodeindexrange =
      noderowmap->MaxAllGID() - noderowmap->MinAllGID() + 1;  // if id's are not continuous numbered

  int maxNumMyReservedDofsperNode = (xdyn.get<int>("MAX_NUM_DOFSETS")) * 4;
  Teuchos::RCP<DRT::FixedSizeDofSet> maxdofset =
      Teuchos::rcp(new DRT::FixedSizeDofSet(maxNumMyReservedDofsperNode, nodeindexrange));

  fluiddis->ReplaceDofSet(xfluid_dofset, maxdofset, true);

  std::vector<int> nds;
  nds.push_back(xfluid_dofset);  // dofsets to be initialized as initial dofset

  xfluiddis->InitialFillComplete(nds);

  //--------------------------------------------------------------------
  // prepare scatra

  // access parameter list for scatra
  const Teuchos::ParameterList& scatradyn = problem->ScalarTransportDynamicParams();

  // Needed for altgeo-creation
  if (!scatradis->Filled()) scatradis->FillComplete();

  // determine coupling type
  const INPAR::SCATRA::FieldCoupling fieldcoupling =
      DRT::INPUT::IntegralValue<INPAR::SCATRA::FieldCoupling>(
          DRT::Problem::Instance()->ScalarTransportDynamicParams(), "FIELDCOUPLING");

  // determine velocity type
  const INPAR::SCATRA::VelocityField veltype =
      DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatradyn, "VELOCITYFIELD");

  if (scatradis->NumGlobalNodes() == 0)
  {
    if (fieldcoupling != INPAR::SCATRA::coupling_match and
        veltype != INPAR::SCATRA::velocity_Navier_Stokes)
      dserror(
          "If you want matching fluid and scatra meshes, do clone you fluid mesh and use "
          "FIELDCOUPLING match!");
  }
  else
  {
    if (fieldcoupling != INPAR::SCATRA::coupling_volmortar and
        veltype == INPAR::SCATRA::velocity_Navier_Stokes)
      dserror(
          "If you want non-matching fluid and scatra meshes, "
          "you need to use FIELDCOUPLING volmortar!");
  }

  // use fluid discretization as layout for scatra discretization
  if (fluiddis->NumGlobalNodes() == 0) dserror("Fluid discretization is empty!");

  // create scatra elements by cloning from fluid dis in matching case
  if (fieldcoupling == INPAR::SCATRA::coupling_match)
  {
    // fill scatra discretization by cloning fluid discretization
    DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(fluiddis, scatradis);

    // Give ScaTra new dofset (starts after fluid)
    Teuchos::RCP<DRT::DofSet> newdofset = Teuchos::rcp(new DRT::DofSet());
    scatradis->ReplaceDofSet(cutter_nds_phi, newdofset, true);
    scatradis->FillComplete();

    // set implementation type of cloned scatra elements
    for (int i = 0; i < scatradis->NumMyColElements(); ++i)
    {
      DRT::ELEMENTS::Transport* element =
          dynamic_cast<DRT::ELEMENTS::Transport*>(scatradis->lColElement(i));
      if (element == NULL)
        dserror("Invalid element type!");
      else
        element->SetImplType(INPAR::SCATRA::impltype_levelset);
    }
  }
  else
    dserror("not supported so far");

  // add dofset proxies from the respective other fields
  if (fieldcoupling == INPAR::SCATRA::coupling_match)
  {
    // add proxy of fluid transport degrees of freedom to scatra discretization
    const int nds_fluid_proxy_in_scatra =
        scatradis->AddDofSet(xfluiddis->GetInitialDofSetProxy(xfluid_dofset));

    // scatradis: dofset(0) = scalar fields
    // scatradis: dofset(1) = the initial fluid dofset
    if (nds_fluid_proxy_in_scatra != 1)
      dserror("Scatra discretization has illegal number of dofsets!");

    dofset_coupling_map.insert(
        std::pair<std::string, int>("vel_fluid_proxy_in_scatra", nds_fluid_proxy_in_scatra));

    // add proxy of scatra level set degrees of freedom to fluid discretization
    const int nds_scatra_proxy_in_fluid =
        xfluiddis->AddDofSet(scatradis->GetDofSetProxy(cutter_nds_phi));
    if (nds_scatra_proxy_in_fluid != 1)
      dserror("Fluid discretization has illegal number of dofsets!");

    dofset_coupling_map.insert(
        std::pair<std::string, int>("phi_scatra_proxy_in_fluid", nds_scatra_proxy_in_fluid));
  }
  else
    dserror("not supported so far");

  //--------------------------------------------------------------------

  // print all dofsets
  //---FLUID---|---SCATRA---|
  fluiddis->GetDofSetProxy()->PrintAllDofsets(fluiddis->Comm());


  // ---------------------------------------------
  //
  //              SAFETY CHECKS!
  //   make sure input is not conflicting
  //
  // ---------------------------------------------

  // get linear solver id from SCALAR TRANSPORT DYNAMIC
  const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    dserror(
        "no linear solver defined for xfem two phase flow (XTPF) problem. Please set LINEAR_SOLVER "
        "in SCALAR TRANSPORT DYNAMIC to a valid number!");

  DoSafetyCheck();

  // ---------------------------------------------- //Safety check done.

  // TODO: why not SCATRA::ScaTraAlgorithm??? Can we unify both?!

  // access parameter for two phase flow
  // these parameters are controling the coupling algorithm.
  const Teuchos::ParameterList& twophasedyn = problem->TwoPhaseFlowParams();

  // Test replacing fdyn in Algorithm with prbdyn
  Teuchos::RCP<XFLUIDLEVELSET::Algorithm> algo = Teuchos::rcp(new XFLUIDLEVELSET::Algorithm(comm,
      twophasedyn, DRT::Problem::Instance()->SolverParams(linsolvernumber), dofset_coupling_map));

  const Teuchos::RCP<ADAPTER::Fluid>& fluid = algo->FluidField();
  Teuchos::RCP<FLD::XFluid> xfluid = Teuchos::rcp_dynamic_cast<FLD::XFluid>(fluid, true);

  xfluid->SetDofSetCouplingMap(dofset_coupling_map);

  algo->Init(twophasedyn, problem->ScalarTransportDynamicParams(),
      DRT::Problem::Instance()->SolverParams(linsolvernumber));

  algo->Setup();

  // read restart information
  // in case an inflow generation in the inflow section has been performed, there are not any
  // scatra results available and the initial field is used
  // FOR NOW ONLY RESTART WITH FLUID AND SCATRA FROM OUTPUT!!!!
  if (restart)
  {
    algo->Restart(restart);
  }

  INPAR::FLUID::TimeIntegrationScheme timeintscheme =
      DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(scatradyn, "TIMEINTEGR");

  if (timeintscheme == INPAR::FLUID::timeint_one_step_theta)
  {
    // solve the two phase problem utilizing the smoothing function for parameter values.
    algo->TimeLoop();
  }
  else if (timeintscheme == INPAR::FLUID::timeint_stationary)
  {
    algo->SolveStationaryProblem();
  }
  else
  {
    // every time integration scheme must be either static or dynamic
    dserror("Time integration schemes not supported. Only OST or Stationary allowed for XFEM.");
  }

  //------------------------------------------------------------------------------------------------
  // validate the results
  //------------------------------------------------------------------------------------------------

  // perform the result test
  algo->TestResults();

  return;
}  // fluid_xfem_ls_drt()

void DoSafetyCheck()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // access the fluid discretization
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");

  // access parameter list for scatra
  const Teuchos::ParameterList& scatradyn = problem->ScalarTransportDynamicParams();

  // Safety check for xfem coupling condition and parameters set in scatra field.
  {
    std::vector<std::string> conditions_to_check;

    // check for two-phase flow and combustion interfaces - however, we allow just one interface so
    // far
    conditions_to_check.push_back("XFEMLevelsetTwophase");
    conditions_to_check.push_back("XFEMLevelsetCombustion");


    std::vector<std::string> names;
    fluiddis->GetConditionNames(names);

    bool condition_found = false;

    for (size_t c = 0; c < conditions_to_check.size(); c++)
    {
      if (std::find(names.begin(), names.end(), conditions_to_check[c]) !=
          names.end())  // condition found
      {
        std::vector<DRT::Condition*> conditions;
        fluiddis->GetCondition(conditions_to_check[c], conditions);
        if (conditions.size() != 1)
          dserror(
              "number of conditions not unique, are there multiple conditions of the same type? - "
              "not allowed so far!");

        // check for unique interface
        if (condition_found == false)
          condition_found = true;
        else
          dserror(
              "did you specify a TWOPHASE and a COMBUSTION condition at the same time? - not "
              "allowed so far!");

        int levelsetfunctnumberxfem = conditions[0]->GetInt("levelsetfieldno");
        int levelsetfunctnumberscatra = scatradyn.get<int>("INITFUNCNO");
        if (levelsetfunctnumberxfem != levelsetfunctnumberscatra)
          dserror(
              "Function number for level-set in SCALAR TRANSPORT DYNAMIC is not same as provided "
              "in DESIGN XFEM LEVELSET TWOPHASE VOL CONDITIONS!");

        // continue for safety - only a unique interface allowed so far!
      }
    }

    if (condition_found == false)
      dserror(
          "No DESIGN XFEM LEVELSET TWOPHASE or COMBUSTION VOL CONDITIONS provided, please add to "
          "input file.");
  }
}

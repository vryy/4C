/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid Base Algorithm

\level 1


 */
/*----------------------------------------------------------------------*/


#include "baci_adapter_fld_base_algorithm.H"

#include "baci_adapter_fld_fbi_wrapper.H"
#include "baci_adapter_fld_fluid_ac_fsi.H"
#include "baci_adapter_fld_fluid_fluid_fsi.H"
#include "baci_adapter_fld_fluid_fpsi.H"
#include "baci_adapter_fld_fluid_fsi.H"
#include "baci_adapter_fld_fluid_fsi_msht.H"
#include "baci_adapter_fld_fluid_xfsi.H"
#include "baci_adapter_fld_lung.H"
#include "baci_adapter_fld_poro.H"
#include "baci_fluid_implicit_integration.H"
#include "baci_fluid_timint_ac_ost.H"
#include "baci_fluid_timint_hdg.H"
#include "baci_fluid_timint_hdg_weak_comp.H"
#include "baci_fluid_timint_loma_bdf2.H"
#include "baci_fluid_timint_loma_genalpha.H"
#include "baci_fluid_timint_loma_ost.H"
#include "baci_fluid_timint_poro_genalpha.H"
#include "baci_fluid_timint_poro_ost.H"
#include "baci_fluid_timint_poro_stat.H"
#include "baci_fluid_timint_red_bdf2.H"
#include "baci_fluid_timint_red_genalpha.H"
#include "baci_fluid_timint_red_ost.H"
#include "baci_fluid_timint_red_stat.H"
#include "baci_fluid_timint_stat_hdg.H"
#include "baci_fluid_xfluid.H"
#include "baci_fluid_xfluid_fluid.H"
#include "baci_inpar_elch.H"
#include "baci_inpar_fluid.H"
#include "baci_inpar_fsi.H"
#include "baci_inpar_poroelast.H"
#include "baci_inpar_solver.H"
#include "baci_inpar_validparameters.H"
#include "baci_inpar_xfem.H"
#include "baci_io.H"
#include "baci_io_control.H"
#include "baci_io_pstream.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_periodicbc.H"
#include "baci_linear_solver_method_linalg.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidBaseAlgorithm::FluidBaseAlgorithm(const Teuchos::ParameterList& prbdyn,
    const Teuchos::ParameterList& fdyn, const std::string& disname, bool isale, bool init)
{
  SetupFluid(prbdyn, fdyn, disname, isale, init);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidBaseAlgorithm::FluidBaseAlgorithm(
    const Teuchos::ParameterList& prbdyn, const Teuchos::RCP<DRT::Discretization> discret)
{
  SetupInflowFluid(prbdyn, discret);
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidBaseAlgorithm::SetupFluid(const Teuchos::ParameterList& prbdyn,
    const Teuchos::ParameterList& fdyn, const std::string& disname, bool isale, bool init)
{
  Teuchos::RCP<Teuchos::Time> t =
      Teuchos::TimeMonitor::getNewTimer("ADAPTER::FluidBaseAlgorithm::SetupFluid");
  Teuchos::TimeMonitor monitor(*t);

  // -------------------------------------------------------------------
  // what's the current problem type?
  // -------------------------------------------------------------------
  ProblemType probtype = DRT::Problem::Instance()->GetProblemType();

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> actdis = DRT::Problem::Instance()->GetDis(disname);

  // -------------------------------------------------------------------
  // connect degrees of freedom for periodic boundary conditions
  // -------------------------------------------------------------------
  if (probtype != ProblemType::fsi)
  {
    PeriodicBoundaryConditions pbc(actdis);
    pbc.UpdateDofsForPeriodicBoundaryConditions();
  }

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->HaveDofs())
  {
    if (probtype == ProblemType::fsi_xfem or probtype == ProblemType::fluid_xfem or
        (probtype == ProblemType::fpsi_xfem and disname == "fluid") or
        probtype == ProblemType::fluid_xfem_ls)
    {
      actdis->FillComplete(false, false, false);
    }
    else
    {
      actdis->FillComplete();
    }
  }

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<IO::DiscretizationWriter> output = actdis->Writer();
  output->WriteMesh(0, 0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  // const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  // const Teuchos::ParameterList& probsize = DRT::Problem::Instance()->ProblemSizeParams();
  // const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();

  if (actdis->Comm().MyPID() == 0) INPUT::PrintDefaultParameters(IO::cout, fdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  Teuchos::RCP<CORE::LINALG::Solver> solver = Teuchos::null;

  switch (INPUT::IntegralValue<INPAR::FLUID::MeshTying>(fdyn, "MESHTYING"))
  {
    case INPAR::FLUID::condensed_bmat:
    {
      // FIXME: The solver should not be taken from the contact dynamic section here,
      // but must be specified in the fluid dynamic section instead (popp 11/2012)

      const Teuchos::ParameterList& mshparams = DRT::Problem::Instance()->ContactDynamicParams();
      const int mshsolver = mshparams.get<int>(
          "LINEAR_SOLVER");  // meshtying solver (with block preconditioner, e.g. BGS 2x2)
      const int fluidsolver = fdyn.get<int>("LINEAR_SOLVER");           // fluid solver
      const int fluidpressuresolver = fdyn.get<int>("SIMPLER_SOLVER");  // fluid pressure solver
      if (mshsolver == (-1))
        dserror(
            "no linear solver defined for fluid meshtying problem. Please set LINEAR_SOLVER in "
            "CONTACT DYNAMIC to a valid number!");
      if (fluidsolver == (-1))
        dserror(
            "no linear solver defined for fluid meshtying problem. Please set LINEAR_SOLVER in "
            "FLUID DYNAMIC to a valid number! This solver is used within block preconditioner "
            "(e.g. BGS2x2) as \"Inverse 1\".");
      if (fluidpressuresolver == (-1))
        dserror(
            "no linear solver defined for fluid meshtying problem. Please set SIMPLER_SOLVER in "
            "FLUID DYNAMIC to a valid number! This solver is used within block preconditioner "
            "(e.g. BGS2x2) as \"Inverse 2\".");

      // check, if meshtying solver is used with a valid block preconditioner
      const auto azprectype = Teuchos::getIntegralValue<INPAR::SOLVER::PreconditionerType>(
          DRT::Problem::Instance()->SolverParams(mshsolver), "AZPREC");

      // plausibility check
      switch (azprectype)
      {
        case INPAR::SOLVER::PreconditionerType::cheap_simple:
        case INPAR::SOLVER::PreconditionerType::block_gauss_seidel_2x2:  // block preconditioners,
                                                                         // that are implemented in
                                                                         // BACI
          break;
        default:
          dserror(
              "Block Gauss-Seidel BGS2x2 preconditioner expected for fluid meshtying problem. "
              "Please set AZPREC to BGS2x2 in solver block %i",
              mshsolver);
          break;
      }

      // create solver objects
      solver = Teuchos::rcp(new CORE::LINALG::Solver(
          DRT::Problem::Instance()->SolverParams(mshsolver), actdis->Comm()));

      // add sub block solvers/smoothers to block preconditioners
      switch (azprectype)
      {
        case INPAR::SOLVER::PreconditionerType::cheap_simple:
          break;  // CheapSIMPLE adds its own Inverse1 and Inverse2 blocks
        case INPAR::SOLVER::PreconditionerType::block_gauss_seidel_2x2:  // block preconditioners,
                                                                         // that are implemented in
                                                                         // BACI
        {
          // set Inverse blocks for block preconditioner
          // for BGS preconditioner
          // This is only necessary for BGS. CheapSIMPLE has a more modern framework
          solver->PutSolverParamsToSubParams(
              "Inverse1", DRT::Problem::Instance()->SolverParams(fluidsolver));

          solver->PutSolverParamsToSubParams(
              "Inverse2", DRT::Problem::Instance()->SolverParams(fluidpressuresolver));
        }
        break;
        default:
          dserror(
              "Block Gauss-Seidel BGS2x2 preconditioner expected for fluid meshtying problem. "
              "Please set AZPREC to BGS2x2 in solver block %i",
              mshsolver);
          break;
      }

      solver->Params().set<bool>("MESHTYING", true);  // mark it as meshtying problem
    }
    break;
    case INPAR::FLUID::condensed_smat:
    case INPAR::FLUID::condensed_bmat_merged:
    {
      // meshtying (no saddle point problem)
      const Teuchos::ParameterList& mshparams = DRT::Problem::Instance()->ContactDynamicParams();
      const int mshsolver = mshparams.get<int>(
          "LINEAR_SOLVER");  // meshtying solver (with block preconditioner, e.g. BGS 2x2)
      if (mshsolver == (-1))
        dserror(
            "no linear solver defined for fluid meshtying problem. Please set LINEAR_SOLVER in "
            "CONTACT DYNAMIC to a valid number!");

      solver = Teuchos::rcp(new CORE::LINALG::Solver(
          DRT::Problem::Instance()->SolverParams(mshsolver), actdis->Comm()));
    }
    break;
    case INPAR::FLUID::no_meshtying:  // no meshtying -> use FLUID SOLVER
    default:
    {
      // default: create solver using the fluid solver params from FLUID SOLVER block

      // get the solver number used for linear fluid solver
      const int linsolvernumber = fdyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == (-1))
        dserror(
            "no linear solver defined for fluid problem. Please set LINEAR_SOLVER in FLUID DYNAMIC "
            "to a valid number!");
      solver = Teuchos::rcp(new CORE::LINALG::Solver(
          DRT::Problem::Instance()->SolverParams(linsolvernumber), actdis->Comm()));

      break;
    }
  }

  // compute null space information
  if (probtype != ProblemType::fsi_xfem and probtype != ProblemType::fpsi_xfem and
      probtype != ProblemType::fluid_xfem and probtype != ProblemType::fluid_xfem_ls and
      !(probtype == ProblemType::fsi and
          INPUT::IntegralValue<bool>(
              DRT::Problem::Instance()->XFluidDynamicParams().sublist("GENERAL"), "XFLUIDFLUID")))
  {
    switch (INPUT::IntegralValue<int>(fdyn, "MESHTYING"))
    {
      // switch types
      case INPAR::FLUID::condensed_bmat:
      {
        const Teuchos::ParameterList& mshparams = DRT::Problem::Instance()->ContactDynamicParams();
        const int mshsolver = mshparams.get<int>(
            "LINEAR_SOLVER");  // meshtying solver (with block preconditioner, e.g. BGS 2x2)

        // check, if meshtying solver is used with a valid block preconditioner
        const auto azprectype = Teuchos::getIntegralValue<INPAR::SOLVER::PreconditionerType>(
            DRT::Problem::Instance()->SolverParams(mshsolver), "AZPREC");

        switch (azprectype)
        {
          // block preconditioners, that are implemented in BACI
          case INPAR::SOLVER::PreconditionerType::cheap_simple:
          {
            actdis->ComputeNullSpaceIfNecessary(
                solver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1"), true);
            actdis->ComputeNullSpaceIfNecessary(
                solver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse2"), true);
          }
          break;
          case INPAR::SOLVER::PreconditionerType::block_gauss_seidel_2x2:
          {
            actdis->ComputeNullSpaceIfNecessary(solver->Params().sublist("Inverse1"), true);
            actdis->ComputeNullSpaceIfNecessary(solver->Params().sublist("Inverse2"), true);
          }
          break;
          default:
          {
          }
        }  // end switch azprectype
      }
      break;
      default:
        // no block matrix
        actdis->ComputeNullSpaceIfNecessary(solver->Params(), true);
        break;
    }
  }

  // create a second solver for SIMPLER preconditioner if chosen from input
  CreateSecondSolver(solver, fdyn);

  // -------------------------------------------------------------------
  // set parameters in list
  // -------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> fluidtimeparams = Teuchos::rcp(new Teuchos::ParameterList());

  // physical type of fluid flow (incompressible, Boussinesq Approximation, varying density, loma,
  // temperature-dependent water, poro)
  fluidtimeparams->set<int>(
      "Physical Type", INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn, "PHYSICAL_TYPE"));
  // and  check correct setting
  if (probtype == ProblemType::loma and (INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(
                                             fdyn, "PHYSICAL_TYPE") != INPAR::FLUID::loma and
                                            INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn,
                                                "PHYSICAL_TYPE") != INPAR::FLUID::tempdepwater))
    dserror(
        "Input parameter PHYSICAL_TYPE in section FLUID DYNAMIC needs to be 'Loma' or "
        "'Temp_dep_water' for low-Mach-number flow!");
  if ((probtype == ProblemType::thermo_fsi) and
      (INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn, "PHYSICAL_TYPE") !=
              INPAR::FLUID::loma and
          INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn, "PHYSICAL_TYPE") !=
              INPAR::FLUID::tempdepwater))
    dserror(
        "Input parameter PHYSICAL_TYPE in section FLUID DYNAMIC needs to be 'Loma' or "
        "'Temp_dep_water' for Thermo-fluid-structure interaction!");
  if ((probtype == ProblemType::poroelast or probtype == ProblemType::poroscatra or
          probtype == ProblemType::fpsi or probtype == ProblemType::fps3i or
          probtype == ProblemType::fpsi_xfem) and
      disname == "porofluid")
  {
    const Teuchos::ParameterList& pedyn = DRT::Problem::Instance()->PoroelastDynamicParams();
    fluidtimeparams->set<int>(
        "Physical Type", INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(pedyn, "PHYSICAL_TYPE"));
    if (fluidtimeparams->get<int>("Physical Type") != INPAR::FLUID::poro and
        fluidtimeparams->get<int>("Physical Type") != INPAR::FLUID::poro_p1)
      dserror(
          "Input parameter PHYSICAL_TYPE in section POROELASTICITY DYNAMIC needs to be 'Poro' or "
          "'Poro_P1' for poro-elasticity!");

    fluidtimeparams->set<int>("Transient Terms Poro Fluid",
        INPUT::IntegralValue<INPAR::POROELAST::TransientEquationsOfPoroFluid>(
            pedyn, "TRANSIENT_TERMS"));
  }

  // now, set general parameters required for all problems
  SetGeneralParameters(fluidtimeparams, prbdyn, fdyn);

  // and, finally, add problem specific parameters

  // for poro problems, use POROUS-FLOW STABILIZATION
  if ((probtype == ProblemType::poroelast or probtype == ProblemType::poroscatra or
          probtype == ProblemType::fpsi or probtype == ProblemType::fps3i or
          probtype == ProblemType::fpsi_xfem) and
      disname == "porofluid")
  {
    fluidtimeparams->sublist("RESIDUAL-BASED STABILIZATION") =
        fdyn.sublist("POROUS-FLOW STABILIZATION");
    fluidtimeparams->sublist("RESIDUAL-BASED STABILIZATION")
        .set<bool>("POROUS-FLOW STABILIZATION", true);
  }

  // add some loma specific parameters
  // get also scatra stabilization sublist
  const Teuchos::ParameterList& lomadyn = DRT::Problem::Instance()->LOMAControlParams();
  fluidtimeparams->sublist("LOMA").set<bool>(
      "update material", INPUT::IntegralValue<int>(lomadyn, "SGS_MATERIAL_UPDATE"));

  // ----------------------------- sublist for general xfem-specific parameters
  if (probtype == ProblemType::fluid_xfem or probtype == ProblemType::fsi_xfem or
      (probtype == ProblemType::fpsi_xfem and disname == "fluid") or
      (probtype == ProblemType::fluid_ale and
          INPUT::IntegralValue<bool>(
              DRT::Problem::Instance()->XFluidDynamicParams().sublist("GENERAL"), "XFLUIDFLUID")) or
      (probtype == ProblemType::fsi and
          INPUT::IntegralValue<bool>(
              DRT::Problem::Instance()->XFluidDynamicParams().sublist("GENERAL"), "XFLUIDFLUID")) or
      probtype == ProblemType::fluid_xfem_ls)
  {
    // get also scatra stabilization sublist
    const Teuchos::ParameterList& xdyn = DRT::Problem::Instance()->XFEMGeneralParams();

    fluidtimeparams->sublist("XFEM") = xdyn;
    // ----------------------------- sublist for xfem-specific fluid parameters
    const Teuchos::ParameterList& xfdyn = DRT::Problem::Instance()->XFluidDynamicParams();

    fluidtimeparams->sublist("XFLUID DYNAMIC/GENERAL") = xfdyn.sublist("GENERAL");
    fluidtimeparams->sublist("XFLUID DYNAMIC/STABILIZATION") = xfdyn.sublist("STABILIZATION");
    fluidtimeparams->sublist("XFLUID DYNAMIC/XFPSI MONOLITHIC") = xfdyn.sublist("XFPSI MONOLITHIC");

    fluidtimeparams->sublist("XFLUID DYNAMIC/GENERAL")
        .set<std::string>("MONOLITHIC_XFFSI_APPROACH",
            xfdyn.sublist("GENERAL").get<std::string>("MONOLITHIC_XFFSI_APPROACH"));
    fluidtimeparams->sublist("XFLUID DYNAMIC/GENERAL")
        .set<double>("XFLUIDFLUID_SEARCHRADIUS",
            xfdyn.sublist("GENERAL").get<double>("XFLUIDFLUID_SEARCHRADIUS"));
  }

  // -------------------------------------------------------------------
  // additional parameters and algorithm call depending on respective
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  INPAR::FLUID::TimeIntegrationScheme timeint =
      INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fdyn, "TIMEINTEGR");

  // sanity checks and default flags
  if (probtype == ProblemType::fsi or probtype == ProblemType::fsi_lung or
      probtype == ProblemType::gas_fsi or probtype == ProblemType::ac_fsi or
      probtype == ProblemType::biofilm_fsi or probtype == ProblemType::thermo_fsi or
      probtype == ProblemType::fsi_xfem or
      (probtype == ProblemType::fpsi_xfem and disname == "fluid") or
      probtype == ProblemType::fsi_redmodels)
  {
    // in case of FSI calculations we do not want a stationary fluid solver
    if (timeint == INPAR::FLUID::timeint_stationary)
      dserror("Stationary fluid solver not allowed for FSI.");

    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

    fluidtimeparams->set<bool>(
        "interface second order", INPUT::IntegralValue<int>(fsidyn, "SECONDORDER"));
    fluidtimeparams->set<bool>(
        "shape derivatives", INPUT::IntegralValue<int>(fsimono, "SHAPEDERIVATIVES"));

    const int coupling = INPUT::IntegralValue<int>(fsidyn, "COUPALGO");

    if (coupling == fsi_iter_lung_monolithicstructuresplit or
        coupling == fsi_iter_lung_monolithicfluidsplit or
        coupling == fsi_iter_constr_monolithicstructuresplit or
        coupling == fsi_iter_constr_monolithicfluidsplit)
    {
      // No explicit predictor for these monolithic FSI schemes, yet.
      // Check, whether fluid predictor is 'steady_state'. Otherwise, throw
      // an error.
      if (fluidtimeparams->get<std::string>("predictor") != "steady_state")
        dserror(
            "No fluid predictor allowed for current monolithic FSI scheme, yet. Use "
            "'steady_state', instead!");
    }
  }
  if (probtype == ProblemType::freesurf)
  {
    // in case of FSI calculations we do not want a stationary fluid solver
    if (timeint == INPAR::FLUID::timeint_stationary)
      dserror("Stationary fluid solver not allowed for free surface problem.");

    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

    fluidtimeparams->set<bool>(
        "interface second order", INPUT::IntegralValue<int>(fsidyn, "SECONDORDER"));
    fluidtimeparams->set<bool>(
        "shape derivatives", INPUT::IntegralValue<int>(fsimono, "SHAPEDERIVATIVES"));

    const int coupling = INPUT::IntegralValue<int>(fsidyn, "COUPALGO");
    if (coupling == fsi_iter_monolithicfluidsplit or coupling == fsi_iter_monolithicstructuresplit)
    {
      // No explicit predictor for monolithic free surface flow schemes, yet.
      // Check, whether fluid predictor is 'steady_state'. Otherwise, throw
      // an error.
      if (fluidtimeparams->get<std::string>("predictor") != "steady_state")
        dserror(
            "No fluid predictor allowed for current monolithic free surface scheme, yet. Use "
            "'steady_state', instead!");
    }
  }

  // sanity checks and default flags
  if (probtype == ProblemType::fluid_xfem)
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    fluidtimeparams->set<bool>(
        "interface second order", INPUT::IntegralValue<int>(fsidyn, "SECONDORDER"));
  }

  // sanity checks and default flags
  if (probtype == ProblemType::fsi_xfem or
      (probtype == ProblemType::fpsi_xfem and disname == "fluid"))
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

    const int coupling = INPUT::IntegralValue<int>(fsidyn, "COUPALGO");

    if (coupling == fsi_iter_monolithicfluidsplit or coupling == fsi_iter_monolithicstructuresplit)
    {
      // there are a couple of restrictions in monolithic FSI
      dserror(
          "for XFSI there is no monolithicfluidsplit or monolithicstructuresplit, use "
          "monolithicxfem or any partitioned algorithm instead");
    }
  }

  // sanity checks and default flags
  if (probtype == ProblemType::fluid_xfem or probtype == ProblemType::fsi_xfem or
      (probtype == ProblemType::fpsi_xfem and disname == "fluid"))
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    const int coupling = INPUT::IntegralValue<int>(fsidyn, "COUPALGO");
    fluidtimeparams->set<int>("COUPALGO", coupling);
  }

  if (probtype == ProblemType::elch)
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    fluidtimeparams->set<bool>(
        "interface second order", INPUT::IntegralValue<int>(fsidyn, "SECONDORDER"));
  }

  if (probtype == ProblemType::poroelast or probtype == ProblemType::poroscatra or
      (probtype == ProblemType::fpsi and disname == "porofluid") or
      (probtype == ProblemType::fps3i and disname == "porofluid") or
      (probtype == ProblemType::fpsi_xfem and disname == "porofluid"))
  {
    const Teuchos::ParameterList& porodyn = DRT::Problem::Instance()->PoroelastDynamicParams();
    fluidtimeparams->set<bool>("poroelast", true);
    fluidtimeparams->set<bool>(
        "interface second order", INPUT::IntegralValue<int>(porodyn, "SECONDORDER"));
    fluidtimeparams->set<bool>("shape derivatives", false);
    fluidtimeparams->set<bool>(
        "conti partial integration", INPUT::IntegralValue<int>(porodyn, "CONTIPARTINT"));
    fluidtimeparams->set<bool>(
        "convective term", INPUT::IntegralValue<bool>(porodyn, "CONVECTIVE_TERM"));
  }
  else if ((probtype == ProblemType::fpsi and disname == "fluid") or
           (probtype == ProblemType::fps3i and disname == "fluid"))
  {
    if (timeint == INPAR::FLUID::timeint_stationary)
      dserror("Stationary fluid solver not allowed for FPSI.");

    fluidtimeparams->set<bool>(
        "interface second order", INPUT::IntegralValue<int>(prbdyn, "SECONDORDER"));
    fluidtimeparams->set<bool>(
        "shape derivatives", INPUT::IntegralValue<int>(prbdyn, "SHAPEDERIVATIVES"));
  }

  // =================================================================================
  // Safety Check for usage of DESIGN SURF VOLUMETRIC FLOW CONDITIONS       AN 06/2014
  // =================================================================================
  if (nullptr != actdis->GetCondition("VolumetricSurfaceFlowCond"))
  {
    if (not(ProblemType::fluid_redmodels == probtype or ProblemType::fsi_redmodels == probtype))
    {
      dserror(
          "ERROR: Given Volumetric Womersly infow condition only works with Problemtyp "
          "Fluid_RedModels or Fluid_Structure_Interaction_RedModels. \n"
          " --> If you want to use this conditions change Problemtyp to Fluid_RedModels or "
          "Fluid_Structure_Interaction_RedModels. \n"
          " --> If you don't want to use this condition comment the respective bcFluid section.");
    }
  }

  // -------------------------------------------------------------------
  // additional parameters and algorithm call depending on respective
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  if (timeint == INPAR::FLUID::timeint_stationary or
      timeint == INPAR::FLUID::timeint_one_step_theta or timeint == INPAR::FLUID::timeint_bdf2 or
      timeint == INPAR::FLUID::timeint_afgenalpha or timeint == INPAR::FLUID::timeint_npgenalpha)
  {
    // -----------------------------------------------------------------
    // set additional parameters in list for
    // one-step-theta/BDF2/af-generalized-alpha/stationary scheme
    // -----------------------------------------------------------------
    // type of time-integration (or stationary) scheme
    fluidtimeparams->set<int>("time int algo", timeint);
    // parameter theta for time-integration schemes
    fluidtimeparams->set<double>("theta", fdyn.get<double>("THETA"));
    // number of steps for potential start algorithm
    fluidtimeparams->set<int>("number of start steps", fdyn.get<int>("NUMSTASTEPS"));
    // parameter theta for potential start algorithm
    fluidtimeparams->set<double>("start theta", fdyn.get<double>("START_THETA"));
    // parameter for grid velocity interpolation
    fluidtimeparams->set<int>("order gridvel", INPUT::IntegralValue<int>(fdyn, "GRIDVEL"));
    // handling of pressure and continuity discretization in new one step theta framework
    fluidtimeparams->set<int>("ost cont and press",
        INPUT::IntegralValue<INPAR::FLUID::OST_Cont_and_Press>(fdyn, "OST_CONT_PRESS"));
    // flag to switch on the new One Step Theta implementation
    bool ostnew = INPUT::IntegralValue<bool>(fdyn, "NEW_OST");
    // if the time integration strategy is not even a one step theta strategy, it cannot be the
    // new one step theta strategy either. As it seems, so far there is no sanity check of the
    // input file
    if (timeint != INPAR::FLUID::timeint_one_step_theta and ostnew)
    {
#ifdef DEBUG
      dserror(
          "You are not using the One Step Theta Integration Strategy in the Fluid solver,\n"
          "but you set the flag NEW_OST to use the new implementation of the One Step Theta "
          "Strategy. \n"
          "This is impossible. \n"
          "Please change your input file!\n");
#endif
      printf(
          "You are not using the One Step Theta Integration Strategy in the Fluid solver,\n"
          "but you set the flag NEW_OST to use the new implementation of the One Step Theta "
          "Strategy. \n"
          "This is impossible. \n"
          "Please change your input file! In this run, NEW_OST is set to false!\n");
      ostnew = false;
    }
    fluidtimeparams->set<bool>("ost new", ostnew);

    bool dirichletcond = true;
    if (probtype == ProblemType::fsi or probtype == ProblemType::fsi_lung or
        probtype == ProblemType::gas_fsi or probtype == ProblemType::ac_fsi or
        probtype == ProblemType::biofilm_fsi or probtype == ProblemType::thermo_fsi or
        probtype == ProblemType::fsi_redmodels)
    {
      // FSI input parameters
      const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
      const int coupling = INPUT::IntegralValue<int>(fsidyn, "COUPALGO");
      if (coupling == fsi_iter_monolithicfluidsplit or
          coupling == fsi_iter_monolithicstructuresplit or
          coupling == fsi_iter_lung_monolithicstructuresplit or
          coupling == fsi_iter_lung_monolithicfluidsplit or
          coupling == fsi_iter_constr_monolithicstructuresplit or
          coupling == fsi_iter_constr_monolithicfluidsplit or
          coupling == fsi_iter_mortar_monolithicstructuresplit or
          coupling == fsi_iter_mortar_monolithicfluidsplit or
          coupling == fsi_iter_mortar_monolithicfluidsplit_saddlepoint or
          coupling == fsi_iter_fluidfluid_monolithicstructuresplit or
          coupling == fsi_iter_fluidfluid_monolithicfluidsplit or
          coupling == fsi_iter_fluidfluid_monolithicstructuresplit_nonox or
          coupling == fsi_iter_fluidfluid_monolithicfluidsplit_nonox or
          coupling == fsi_iter_sliding_monolithicfluidsplit or
          coupling == fsi_iter_sliding_monolithicstructuresplit)
      {
        dirichletcond = false;
      }
    }

    if (probtype == ProblemType::poroelast or probtype == ProblemType::poroscatra or
        probtype == ProblemType::fpsi or probtype == ProblemType::fps3i or
        (probtype == ProblemType::fpsi_xfem and disname == "porofluid"))
      dirichletcond = false;

    //------------------------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor);
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom

    switch (probtype)
    {
      case ProblemType::fluid:
      case ProblemType::scatra:
      {
        // HDG implements all time stepping schemes within gen-alpha
        if (DRT::Problem::Instance()->SpatialApproximationType() ==
                CORE::FE::ShapeFunctionType::hdg &&
            timeint != INPAR::FLUID::timeint_stationary &&
            INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn, "PHYSICAL_TYPE") !=
                INPAR::FLUID::weakly_compressible_dens_mom &&
            INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn, "PHYSICAL_TYPE") !=
                INPAR::FLUID::weakly_compressible_stokes_dens_mom)
          fluid_ = Teuchos::rcp(new FLD::TimIntHDG(actdis, solver, fluidtimeparams, output, isale));
        else if (DRT::Problem::Instance()->SpatialApproximationType() ==
                     CORE::FE::ShapeFunctionType::hdg &&
                 (INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn, "PHYSICAL_TYPE") ==
                         INPAR::FLUID::weakly_compressible_dens_mom ||
                     INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn, "PHYSICAL_TYPE") ==
                         INPAR::FLUID::weakly_compressible_stokes_dens_mom))
          fluid_ = Teuchos::rcp(
              new FLD::TimIntHDGWeakComp(actdis, solver, fluidtimeparams, output, isale));
        else if (DRT::Problem::Instance()->SpatialApproximationType() ==
                     CORE::FE::ShapeFunctionType::hdg &&
                 timeint == INPAR::FLUID::timeint_stationary)
          fluid_ = Teuchos::rcp(
              new FLD::TimIntStationaryHDG(actdis, solver, fluidtimeparams, output, isale));
        else if (timeint == INPAR::FLUID::timeint_stationary)
          fluid_ = Teuchos::rcp(
              new FLD::TimIntStationary(actdis, solver, fluidtimeparams, output, isale));
        else if (timeint == INPAR::FLUID::timeint_one_step_theta)
          fluid_ = Teuchos::rcp(
              new FLD::TimIntOneStepTheta(actdis, solver, fluidtimeparams, output, isale));
        else if (timeint == INPAR::FLUID::timeint_bdf2)
          fluid_ =
              Teuchos::rcp(new FLD::TimIntBDF2(actdis, solver, fluidtimeparams, output, isale));
        else if (timeint == INPAR::FLUID::timeint_afgenalpha or
                 timeint == INPAR::FLUID::timeint_npgenalpha)
          fluid_ =
              Teuchos::rcp(new FLD::TimIntGenAlpha(actdis, solver, fluidtimeparams, output, isale));
        else
          dserror("Unknown time integration for this fluid problem type\n");
      }
      break;
      case ProblemType::fluid_redmodels:
      {
        if (timeint == INPAR::FLUID::timeint_stationary)
          fluid_ = Teuchos::rcp(
              new FLD::TimIntRedModelsStat(actdis, solver, fluidtimeparams, output, isale));
        else if (timeint == INPAR::FLUID::timeint_one_step_theta)
          fluid_ = Teuchos::rcp(
              new FLD::TimIntRedModelsOst(actdis, solver, fluidtimeparams, output, isale));
        else if (timeint == INPAR::FLUID::timeint_afgenalpha or
                 timeint == INPAR::FLUID::timeint_npgenalpha)
          fluid_ = Teuchos::rcp(
              new FLD::TimIntRedModelsGenAlpha(actdis, solver, fluidtimeparams, output, isale));
        else if (timeint == INPAR::FLUID::timeint_bdf2)
          fluid_ = Teuchos::rcp(
              new FLD::TimIntRedModelsBDF2(actdis, solver, fluidtimeparams, output, isale));
        else
          dserror("Unknown time integration for this fluid problem type\n");
      }
      break;
      case ProblemType::loma:
      {
        if (INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn, "PHYSICAL_TYPE") ==
            INPAR::FLUID::tempdepwater)
        {
          if (timeint == INPAR::FLUID::timeint_afgenalpha or
              timeint == INPAR::FLUID::timeint_npgenalpha)
            fluid_ = Teuchos::rcp(
                new FLD::TimIntGenAlpha(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_one_step_theta)
            fluid_ = Teuchos::rcp(
                new FLD::TimIntOneStepTheta(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_bdf2)
            fluid_ =
                Teuchos::rcp(new FLD::TimIntBDF2(actdis, solver, fluidtimeparams, output, isale));
          else
            dserror("Unknown time integration for this fluid problem type\n");
        }
        else
        {
          if (timeint == INPAR::FLUID::timeint_afgenalpha or
              timeint == INPAR::FLUID::timeint_npgenalpha)
            fluid_ = Teuchos::rcp(
                new FLD::TimIntLomaGenAlpha(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_one_step_theta)
            fluid_ = Teuchos::rcp(
                new FLD::TimIntLomaOst(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_bdf2)
            fluid_ = Teuchos::rcp(
                new FLD::TimIntLomaBDF2(actdis, solver, fluidtimeparams, output, isale));
          else
            dserror("Unknown time integration for this fluid problem type\n");
        }
      }
      break;
      case ProblemType::fluid_xfem:
      {
        if (INPUT::IntegralValue<bool>(
                DRT::Problem::Instance()->XFluidDynamicParams().sublist("GENERAL"), "XFLUIDFLUID"))
        {
          // actdis is the embedded fluid discretization
          Teuchos::RCP<DRT::Discretization> xfluiddis = DRT::Problem::Instance()->GetDis("xfluid");

          Teuchos::RCP<FLD::FluidImplicitTimeInt> tmpfluid;
          if (timeint == INPAR::FLUID::timeint_stationary)
            tmpfluid = Teuchos::rcp(
                new FLD::TimIntStationary(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_one_step_theta)
            tmpfluid = Teuchos::rcp(
                new FLD::TimIntOneStepTheta(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_bdf2)
            tmpfluid =
                Teuchos::rcp(new FLD::TimIntBDF2(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_afgenalpha or
                   timeint == INPAR::FLUID::timeint_npgenalpha)
            tmpfluid = Teuchos::rcp(
                new FLD::TimIntGenAlpha(actdis, solver, fluidtimeparams, output, isale));
          else
            dserror("Unknown time integration for this fluid problem type\n");

          fluid_ = Teuchos::rcp(
              new FLD::XFluidFluid(tmpfluid, xfluiddis, solver, fluidtimeparams, isale));
          break;
        }

        Teuchos::RCP<DRT::Discretization> soliddis = DRT::Problem::Instance()->GetDis("structure");
        Teuchos::RCP<DRT::Discretization> scatradis = Teuchos::null;

        if (DRT::Problem::Instance()->DoesExistDis("scatra"))
          scatradis = DRT::Problem::Instance()->GetDis("scatra");

        Teuchos::RCP<FLD::XFluid> tmpfluid = Teuchos::rcp(
            new FLD::XFluid(actdis, soliddis, scatradis, solver, fluidtimeparams, output, isale));

        std::string condition_name = "";

        // TODO: actually in case of ale fluid with e.g. only level-set we do not want to use the
        // XFluidFSI class since not always
        // a boundary discretization is necessary.
        // however, the xfluid-class itself does not support the full ALE-functionality without the
        // FSI itself ALE-fluid with level-set/without mesh discretization not supported yet
        if (isale)  // in ale case
          fluid_ = Teuchos::rcp(
              new XFluidFSI(tmpfluid, condition_name, solver, fluidtimeparams, output));
        else
          fluid_ = tmpfluid;
      }
      break;
      case ProblemType::fsi_xfem:
      {
        std::string condition_name;

        // FSI input parameters
        const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
        const int coupling = INPUT::IntegralValue<int>(fsidyn, "COUPALGO");
        if (coupling == fsi_iter_xfem_monolithic)
        {
          condition_name = "XFEMSurfFSIMono";  // not used anymore!
        }
        else if (coupling == fsi_iter_stagg_fixed_rel_param or
                 coupling == fsi_iter_stagg_AITKEN_rel_param or
                 coupling == fsi_iter_stagg_steep_desc or
                 coupling == fsi_iter_stagg_CHEB_rel_param or
                 coupling == fsi_iter_stagg_AITKEN_rel_force or
                 coupling == fsi_iter_stagg_steep_desc_force or
                 coupling == fsi_iter_stagg_steep_desc_force or
                 coupling == fsi_iter_stagg_steep_desc_force)
        {
          condition_name = "XFEMSurfFSIPart";
        }
        else
          dserror("non supported COUPALGO for FSI");

        Teuchos::RCP<DRT::Discretization> soliddis = DRT::Problem::Instance()->GetDis("structure");
        Teuchos::RCP<FLD::XFluid> tmpfluid;
        if (INPUT::IntegralValue<bool>(
                DRT::Problem::Instance()->XFluidDynamicParams().sublist("GENERAL"), "XFLUIDFLUID"))
        {
          dserror(
              "XFLUIDFLUID with XFSI framework not supported via FLD::XFluidFluid but via "
              "FLD::XFluid");

          // actdis is the embedded fluid discretization
          Teuchos::RCP<DRT::Discretization> xfluiddis = DRT::Problem::Instance()->GetDis("xfluid");

          Teuchos::RCP<FLD::FluidImplicitTimeInt> tmpfluid_emb;
          if (timeint == INPAR::FLUID::timeint_stationary)
            tmpfluid_emb = Teuchos::rcp(
                new FLD::TimIntStationary(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_one_step_theta)
            tmpfluid_emb = Teuchos::rcp(
                new FLD::TimIntOneStepTheta(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_bdf2)
            tmpfluid_emb =
                Teuchos::rcp(new FLD::TimIntBDF2(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_afgenalpha or
                   timeint == INPAR::FLUID::timeint_npgenalpha)
            tmpfluid_emb = Teuchos::rcp(
                new FLD::TimIntGenAlpha(actdis, solver, fluidtimeparams, output, isale));
          else
            dserror("Unknown time integration for this fluid problem type\n");

          tmpfluid = Teuchos::rcp(new FLD::XFluidFluid(
              tmpfluid_emb, xfluiddis, soliddis, solver, fluidtimeparams, isale));
        }
        else
        {
          Teuchos::RCP<DRT::Discretization> scatradis = Teuchos::null;

          if (DRT::Problem::Instance()->DoesExistDis("scatra"))
            scatradis = DRT::Problem::Instance()->GetDis("scatra");

          tmpfluid = Teuchos::rcp(
              new FLD::XFluid(actdis, soliddis, scatradis, solver, fluidtimeparams, output, isale));
        }

        if (coupling == fsi_iter_xfem_monolithic)
          fluid_ = tmpfluid;
        else
          fluid_ = Teuchos::rcp(
              new XFluidFSI(tmpfluid, condition_name, solver, fluidtimeparams, output));
      }
      break;
      case ProblemType::fluid_xfem_ls:
      {
        Teuchos::RCP<DRT::Discretization> soliddis = DRT::Problem::Instance()->GetDis("structure");
        Teuchos::RCP<DRT::Discretization> scatradis = Teuchos::null;

        if (DRT::Problem::Instance()->DoesExistDis("scatra"))
          scatradis = DRT::Problem::Instance()->GetDis("scatra");

        fluid_ = Teuchos::rcp(
            new FLD::XFluid(actdis, soliddis, scatradis, solver, fluidtimeparams, output));
      }
      break;
      case ProblemType::fsi:
      case ProblemType::immersed_fsi:
      case ProblemType::gas_fsi:
      case ProblemType::biofilm_fsi:
      case ProblemType::fbi:
      case ProblemType::fluid_ale:
      case ProblemType::freesurf:
      {  //
        Teuchos::RCP<FLD::FluidImplicitTimeInt> tmpfluid;
        if (DRT::Problem::Instance()->SpatialApproximationType() ==
                CORE::FE::ShapeFunctionType::hdg &&
            (INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn, "PHYSICAL_TYPE") ==
                    INPAR::FLUID::weakly_compressible_dens_mom ||
                INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn, "PHYSICAL_TYPE") ==
                    INPAR::FLUID::weakly_compressible_stokes_dens_mom))
          tmpfluid = Teuchos::rcp(
              new FLD::TimIntHDGWeakComp(actdis, solver, fluidtimeparams, output, isale));
        else if (timeint == INPAR::FLUID::timeint_stationary)
          tmpfluid = Teuchos::rcp(
              new FLD::TimIntStationary(actdis, solver, fluidtimeparams, output, isale));
        else if (timeint == INPAR::FLUID::timeint_one_step_theta)
          tmpfluid = Teuchos::rcp(
              new FLD::TimIntOneStepTheta(actdis, solver, fluidtimeparams, output, isale));
        else if (timeint == INPAR::FLUID::timeint_bdf2)
          tmpfluid =
              Teuchos::rcp(new FLD::TimIntBDF2(actdis, solver, fluidtimeparams, output, isale));
        else if (timeint == INPAR::FLUID::timeint_afgenalpha or
                 timeint == INPAR::FLUID::timeint_npgenalpha)
          tmpfluid =
              Teuchos::rcp(new FLD::TimIntGenAlpha(actdis, solver, fluidtimeparams, output, isale));
        else
          dserror("Unknown time integration for this fluid problem type\n");

        const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
        int coupling = INPUT::IntegralValue<int>(fsidyn, "COUPALGO");

        if (INPUT::IntegralValue<bool>(
                DRT::Problem::Instance()->XFluidDynamicParams().sublist("GENERAL"), "XFLUIDFLUID"))
        {
          fluidtimeparams->set<bool>("shape derivatives", false);
          // actdis is the embedded fluid discretization
          Teuchos::RCP<DRT::Discretization> xfluiddis = DRT::Problem::Instance()->GetDis("xfluid");
          Teuchos::RCP<FLD::XFluidFluid> xffluid = Teuchos::rcp(
              new FLD::XFluidFluid(tmpfluid, xfluiddis, solver, fluidtimeparams, false, isale));
          fluid_ = Teuchos::rcp(
              new FluidFluidFSI(xffluid, tmpfluid, solver, fluidtimeparams, isale, dirichletcond));
        }
        else if (coupling == fsi_iter_sliding_monolithicfluidsplit or
                 coupling == fsi_iter_sliding_monolithicstructuresplit)
          fluid_ = Teuchos::rcp(new FluidFSIMsht(
              tmpfluid, actdis, solver, fluidtimeparams, output, isale, dirichletcond));
        else if (probtype == ProblemType::fbi)
          fluid_ = Teuchos::rcp(new FluidFBI(
              tmpfluid, actdis, solver, fluidtimeparams, output, isale, dirichletcond));
        else
          fluid_ = Teuchos::rcp(new FluidFSI(
              tmpfluid, actdis, solver, fluidtimeparams, output, isale, dirichletcond));
      }
      break;
      case ProblemType::thermo_fsi:
      {
        Teuchos::RCP<FLD::FluidImplicitTimeInt> tmpfluid;
        if (INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn, "PHYSICAL_TYPE") ==
            INPAR::FLUID::tempdepwater)
        {
          if (timeint == INPAR::FLUID::timeint_afgenalpha or
              timeint == INPAR::FLUID::timeint_npgenalpha)
            tmpfluid = Teuchos::rcp(
                new FLD::TimIntGenAlpha(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_one_step_theta)
            tmpfluid = Teuchos::rcp(
                new FLD::TimIntOneStepTheta(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_bdf2)
            tmpfluid =
                Teuchos::rcp(new FLD::TimIntBDF2(actdis, solver, fluidtimeparams, output, isale));
          else
            dserror("Unknown time integration for this fluid problem type\n");
        }
        else
        {
          if (timeint == INPAR::FLUID::timeint_afgenalpha or
              timeint == INPAR::FLUID::timeint_npgenalpha)
            tmpfluid = Teuchos::rcp(
                new FLD::TimIntLomaGenAlpha(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_one_step_theta)
            tmpfluid = Teuchos::rcp(
                new FLD::TimIntLomaOst(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_bdf2)
            tmpfluid = Teuchos::rcp(
                new FLD::TimIntLomaBDF2(actdis, solver, fluidtimeparams, output, isale));
          else
            dserror("Unknown time integration for this fluid problem type\n");
        }

        const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
        int coupling = INPUT::IntegralValue<int>(fsidyn, "COUPALGO");

        if (coupling == fsi_iter_sliding_monolithicfluidsplit or
            coupling == fsi_iter_sliding_monolithicstructuresplit)
          fluid_ = Teuchos::rcp(new FluidFSIMsht(
              tmpfluid, actdis, solver, fluidtimeparams, output, isale, dirichletcond));
        else
          fluid_ = Teuchos::rcp(new FluidFSI(
              tmpfluid, actdis, solver, fluidtimeparams, output, isale, dirichletcond));
      }
      break;
      case ProblemType::ac_fsi:
      {  //
        Teuchos::RCP<FLD::FluidImplicitTimeInt> tmpfluid;
        if (timeint == INPAR::FLUID::timeint_one_step_theta)
          tmpfluid =
              Teuchos::rcp(new FLD::TimIntACOst(actdis, solver, fluidtimeparams, output, isale));
        else
          dserror("Unknown time integration for this fluid problem type\n");

        fluid_ = Teuchos::rcp(new FluidACFSI(
            tmpfluid, actdis, solver, fluidtimeparams, output, isale, dirichletcond));
      }
      break;
      case ProblemType::fsi_redmodels:
      {  // give a warning
        if (actdis->Comm().MyPID() == 0)
          std::cout << "\n Warning: FSI_RedModels is little tested. Keep testing! \n" << std::endl;

        // create the fluid time integration object
        Teuchos::RCP<FLD::FluidImplicitTimeInt> tmpfluid;
        if (timeint == INPAR::FLUID::timeint_stationary)
          tmpfluid = Teuchos::rcp(
              new FLD::TimIntRedModelsStat(actdis, solver, fluidtimeparams, output, isale));
        else if (timeint == INPAR::FLUID::timeint_one_step_theta)
          tmpfluid = Teuchos::rcp(
              new FLD::TimIntRedModelsOst(actdis, solver, fluidtimeparams, output, isale));
        else if (timeint == INPAR::FLUID::timeint_afgenalpha or
                 timeint == INPAR::FLUID::timeint_npgenalpha)
          tmpfluid = Teuchos::rcp(
              new FLD::TimIntRedModelsGenAlpha(actdis, solver, fluidtimeparams, output, isale));
        else if (timeint == INPAR::FLUID::timeint_bdf2)
          tmpfluid = Teuchos::rcp(
              new FLD::TimIntRedModelsBDF2(actdis, solver, fluidtimeparams, output, isale));
        else
          dserror("Unknown time integration for this fluid problem type\n");
        fluid_ = Teuchos::rcp(
            new FluidFSI(tmpfluid, actdis, solver, fluidtimeparams, output, isale, dirichletcond));
      }
      break;
      case ProblemType::fsi_lung:
      {
        Teuchos::RCP<FLD::FluidImplicitTimeInt> tmpfluid;
        if (timeint == INPAR::FLUID::timeint_stationary)
          tmpfluid = Teuchos::rcp(
              new FLD::TimIntRedModelsStat(actdis, solver, fluidtimeparams, output, isale));
        else if (timeint == INPAR::FLUID::timeint_one_step_theta)
          tmpfluid = Teuchos::rcp(
              new FLD::TimIntRedModelsOst(actdis, solver, fluidtimeparams, output, isale));
        else if (timeint == INPAR::FLUID::timeint_afgenalpha or
                 timeint == INPAR::FLUID::timeint_npgenalpha)
          tmpfluid = Teuchos::rcp(
              new FLD::TimIntRedModelsGenAlpha(actdis, solver, fluidtimeparams, output, isale));
        else if (timeint == INPAR::FLUID::timeint_bdf2)
          tmpfluid = Teuchos::rcp(
              new FLD::TimIntRedModelsBDF2(actdis, solver, fluidtimeparams, output, isale));
        else
          dserror("Unknown time integration for this fluid problem type\n");
        fluid_ = Teuchos::rcp(
            new FluidLung(tmpfluid, actdis, solver, fluidtimeparams, output, isale, dirichletcond));
      }
      break;
      case ProblemType::poroelast:
      case ProblemType::poroscatra:
      case ProblemType::fpsi:
      case ProblemType::fps3i:
      case ProblemType::fpsi_xfem:
      {
        Teuchos::RCP<FLD::FluidImplicitTimeInt> tmpfluid;
        if (disname == "porofluid")
        {
          if (timeint == INPAR::FLUID::timeint_stationary)
            tmpfluid = Teuchos::rcp(
                new FLD::TimIntPoroStat(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_one_step_theta)
            tmpfluid = Teuchos::rcp(
                new FLD::TimIntPoroOst(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_afgenalpha or
                   timeint == INPAR::FLUID::timeint_npgenalpha)
            tmpfluid = Teuchos::rcp(
                new FLD::TimIntPoroGenAlpha(actdis, solver, fluidtimeparams, output, isale));
          else
            dserror("Unknown time integration for this fluid problem type\n");
          fluid_ = Teuchos::rcp(new FluidPoro(
              tmpfluid, actdis, solver, fluidtimeparams, output, isale, dirichletcond));
        }
        else if (disname == "fluid")
        {
          if (probtype == ProblemType::fpsi or probtype == ProblemType::fps3i)
          {
            if (timeint == INPAR::FLUID::timeint_stationary)
              tmpfluid = Teuchos::rcp(
                  new FLD::TimIntStationary(actdis, solver, fluidtimeparams, output, isale));
            else if (timeint == INPAR::FLUID::timeint_one_step_theta)
              tmpfluid = Teuchos::rcp(
                  new FLD::TimIntOneStepTheta(actdis, solver, fluidtimeparams, output, isale));
            else
              dserror("Unknown time integration for this fluid problem type\n");
            fluid_ = Teuchos::rcp(new FluidFPSI(
                tmpfluid, actdis, solver, fluidtimeparams, output, isale, dirichletcond));
          }
          else if (probtype == ProblemType::fpsi_xfem)
          {
            Teuchos::RCP<DRT::Discretization> soliddis =
                DRT::Problem::Instance()->GetDis("structure");
            Teuchos::RCP<DRT::Discretization> scatradis = Teuchos::null;

            if (DRT::Problem::Instance()->DoesExistDis("scatra"))
              scatradis = DRT::Problem::Instance()->GetDis("scatra");

            fluid_ = Teuchos::rcp(new FLD::XFluid(
                actdis, soliddis, scatradis, solver, fluidtimeparams, output, isale));
          }
        }
      }
      break;
      case ProblemType::elch:
      {
        // access the problem-specific parameter list
        const Teuchos::ParameterList& elchcontrol = DRT::Problem::Instance()->ELCHControlParams();
        // is ALE needed or not?
        const INPAR::ELCH::ElchMovingBoundary withale =
            INPUT::IntegralValue<INPAR::ELCH::ElchMovingBoundary>(elchcontrol, "MOVINGBOUNDARY");
        if (withale != INPAR::ELCH::elch_mov_bndry_no)
        {
          Teuchos::RCP<FLD::FluidImplicitTimeInt> tmpfluid;
          if (timeint == INPAR::FLUID::timeint_stationary)
            tmpfluid = Teuchos::rcp(
                new FLD::TimIntStationary(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_one_step_theta)
            tmpfluid = Teuchos::rcp(
                new FLD::TimIntOneStepTheta(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_bdf2)
            tmpfluid =
                Teuchos::rcp(new FLD::TimIntBDF2(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_afgenalpha or
                   timeint == INPAR::FLUID::timeint_npgenalpha)
            tmpfluid = Teuchos::rcp(
                new FLD::TimIntGenAlpha(actdis, solver, fluidtimeparams, output, isale));
          else
            dserror("Unknown time integration for this fluid problem type\n");
          fluid_ = Teuchos::rcp(new FluidFSI(
              tmpfluid, actdis, solver, fluidtimeparams, output, isale, dirichletcond));
        }
        else
        {
          if (timeint == INPAR::FLUID::timeint_stationary)
            fluid_ = Teuchos::rcp(
                new FLD::TimIntStationary(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_one_step_theta)
            fluid_ = Teuchos::rcp(
                new FLD::TimIntOneStepTheta(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_bdf2)
            fluid_ =
                Teuchos::rcp(new FLD::TimIntBDF2(actdis, solver, fluidtimeparams, output, isale));
          else if (timeint == INPAR::FLUID::timeint_afgenalpha or
                   timeint == INPAR::FLUID::timeint_npgenalpha)
            fluid_ = Teuchos::rcp(
                new FLD::TimIntGenAlpha(actdis, solver, fluidtimeparams, output, isale));
          else
            dserror("Unknown time integration for this fluid problem type\n");
        }
      }
      break;
      default:
      {
        dserror("Undefined problem type.");
      }
      break;
    }  // end switch (probtype)
  }
  else
  {
    dserror("Unknown time integration for fluid\n");
  }

  // initialize algorithm for specific time-integration scheme
  if (init)
  {
    fluid_->Init();

    SetInitialFlowField(fdyn);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidBaseAlgorithm::SetInitialFlowField(const Teuchos::ParameterList& fdyn)
{
  // set initial field by given function
  // we do this here, since we have direct access to all necessary parameters
  INPAR::FLUID::InitialField initfield =
      INPUT::IntegralValue<INPAR::FLUID::InitialField>(fdyn, "INITIALFIELD");
  if (initfield != INPAR::FLUID::initfield_zero_field)
  {
    int startfuncno = fdyn.get<int>("STARTFUNCNO");
    if (initfield != INPAR::FLUID::initfield_field_by_function and
        initfield != INPAR::FLUID::initfield_disturbed_field_from_function)
    {
      startfuncno = -1;
    }
    fluid_->SetInitialFlowField(initfield, startfuncno);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidBaseAlgorithm::SetInitialInflowField(const Teuchos::ParameterList& fdyn)
{
  // set initial field for inflow section by given function
  // we do this here, since we have direct access to all necessary parameters
  INPAR::FLUID::InitialField initfield = INPUT::IntegralValue<INPAR::FLUID::InitialField>(
      fdyn.sublist("TURBULENT INFLOW"), "INITIALINFLOWFIELD");
  if (initfield != INPAR::FLUID::initfield_zero_field)
  {
    int startfuncno = fdyn.sublist("TURBULENT INFLOW").get<int>("INFLOWFUNC");
    if (initfield != INPAR::FLUID::initfield_field_by_function and
        initfield != INPAR::FLUID::initfield_disturbed_field_from_function)
    {
      startfuncno = -1;
    }
    fluid_->SetInitialFlowField(initfield, startfuncno);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidBaseAlgorithm::SetupInflowFluid(
    const Teuchos::ParameterList& prbdyn, const Teuchos::RCP<DRT::Discretization> discret)
{
  Teuchos::RCP<Teuchos::Time> t =
      Teuchos::TimeMonitor::getNewTimer("ADAPTER::FluidBaseAlgorithm::SetupFluid");
  Teuchos::TimeMonitor monitor(*t);

  // -------------------------------------------------------------------
  // what's the current problem type?
  // -------------------------------------------------------------------
  ProblemType probtype = DRT::Problem::Instance()->GetProblemType();

  // the inflow computation can only deal with standard fluid problems so far
  // extensions for xfluid, fsi problems have to be added if necessary
  // they should not pose any additional problem
  // meshtying or xfem related parameters are not supported, yet
  if (probtype != ProblemType::fluid)
    dserror("Only fluid problems supported! Read comment and add your problem type!");

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!discret->HaveDofs())
  {
    dserror("FillComplete shouldn't be necessary!");
  }

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<IO::DiscretizationWriter> output = discret->Writer();

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();

  if (discret->Comm().MyPID() == 0) INPUT::PrintDefaultParameters(IO::cout, fdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // get the solver number used for linear fluid solver
  const int linsolvernumber = fdyn.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    dserror(
        "no linear solver defined for fluid problem. Please set LINEAR_SOLVER in FLUID DYNAMIC to "
        "a valid number!");
  Teuchos::RCP<CORE::LINALG::Solver> solver = Teuchos::rcp(new CORE::LINALG::Solver(
      DRT::Problem::Instance()->SolverParams(linsolvernumber), discret->Comm()));

  discret->ComputeNullSpaceIfNecessary(solver->Params(), true);

  // create a second solver for SIMPLER preconditioner if chosen from input
  CreateSecondSolver(solver, fdyn);

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> fluidtimeparams = Teuchos::rcp(new Teuchos::ParameterList());

  // physical type of fluid flow (incompressible, Boussinesq Approximation, varying density, loma,
  // temperature-dependent water)
  fluidtimeparams->set<int>(
      "Physical Type", INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn, "PHYSICAL_TYPE"));

  // now, set general parameters required for all problems
  SetGeneralParameters(fluidtimeparams, prbdyn, fdyn);

  // overwrite canonical flow parameters by inflow type
  fluidtimeparams->sublist("TURBULENCE MODEL")
      .set<std::string>(
          "CANONICAL_FLOW", fdyn.sublist("TURBULENT INFLOW").get<std::string>("CANONICAL_INFLOW"));
  fluidtimeparams->sublist("TURBULENCE MODEL")
      .set<std::string>(
          "HOMDIR", fdyn.sublist("TURBULENT INFLOW").get<std::string>("INFLOW_HOMDIR"));
  fluidtimeparams->sublist("TURBULENCE MODEL")
      .set<int>(
          "DUMPING_PERIOD", fdyn.sublist("TURBULENT INFLOW").get<int>("INFLOW_DUMPING_PERIOD"));
  fluidtimeparams->sublist("TURBULENCE MODEL")
      .set<int>(
          "SAMPLING_START", fdyn.sublist("TURBULENT INFLOW").get<int>("INFLOW_SAMPLING_START"));
  fluidtimeparams->sublist("TURBULENCE MODEL")
      .set<int>("SAMPLING_STOP", fdyn.sublist("TURBULENT INFLOW").get<int>("INFLOW_SAMPLING_STOP"));
  fluidtimeparams->sublist("TURBULENCE MODEL")
      .set<double>(
          "CHAN_AMPL_INIT_DIST", fdyn.sublist("TURBULENT INFLOW").get<double>("INFLOW_INIT_DIST"));

  // -------------------------------------------------------------------
  // additional parameters and algorithm call depending on respective
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  INPAR::FLUID::TimeIntegrationScheme timeint =
      INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fdyn, "TIMEINTEGR");

  // -------------------------------------------------------------------
  // additional parameters and algorithm call depending on respective
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  if (timeint == INPAR::FLUID::timeint_stationary or
      timeint == INPAR::FLUID::timeint_one_step_theta or timeint == INPAR::FLUID::timeint_bdf2 or
      timeint == INPAR::FLUID::timeint_afgenalpha or timeint == INPAR::FLUID::timeint_npgenalpha)
  {
    // -----------------------------------------------------------------
    // set additional parameters in list for
    // one-step-theta/BDF2/af-generalized-alpha/stationary scheme
    // -----------------------------------------------------------------
    // type of time-integration (or stationary) scheme
    fluidtimeparams->set<int>("time int algo", timeint);
    // parameter theta for time-integration schemes
    fluidtimeparams->set<double>("theta", fdyn.get<double>("THETA"));
    // number of steps for potential start algorithm
    fluidtimeparams->set<int>("number of start steps", fdyn.get<int>("NUMSTASTEPS"));
    // parameter theta for potential start algorithm
    fluidtimeparams->set<double>("start theta", fdyn.get<double>("START_THETA"));
    // parameter for grid velocity interpolation
    fluidtimeparams->set<int>("order gridvel", INPUT::IntegralValue<int>(fdyn, "GRIDVEL"));
    // handling of pressure and continuity discretization in new one step theta framework
    fluidtimeparams->set<int>("ost cont and press",
        INPUT::IntegralValue<INPAR::FLUID::OST_Cont_and_Press>(fdyn, "OST_CONT_PRESS"));
    // flag to switch on the new One Step Theta implementation
    bool ostnew = INPUT::IntegralValue<bool>(fdyn, "NEW_OST");
    // if the time integration strategy is not even a one step theta strategy, it cannot be the
    // new one step theta strategy either. As it seems, so far there is no sanity check of the
    // input file
    if (timeint != INPAR::FLUID::timeint_one_step_theta and ostnew)
    {
#ifdef DEBUG
      dserror(
          "You are not using the One Step Theta Integration Strategy in the Fluid solver,\n"
          "but you set the flag NEW_OST to use the new implementation of the One Step Theta "
          "Strategy. \n"
          "This is impossible. \n"
          "Please change your input file!\n");
#endif
      printf(
          "You are not using the One Step Theta Integration Strategy in the Fluid solver,\n"
          "but you set the flag NEW_OST to use the new implementation of the One Step Theta "
          "Strategy. \n"
          "This is impossible. \n"
          "Please change your input file! In this run, NEW_OST is set to false!\n");
      ostnew = false;
    }
    fluidtimeparams->set<bool>("ost new", ostnew);

    //------------------------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor);
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom
    //    fluid_ = Teuchos::rcp(new FLD::FluidImplicitTimeInt(discret, solver, fluidtimeparams,
    //    output, false));
    if (timeint == INPAR::FLUID::timeint_stationary)
      fluid_ =
          Teuchos::rcp(new FLD::TimIntStationary(discret, solver, fluidtimeparams, output, false));
    else if (timeint == INPAR::FLUID::timeint_one_step_theta)
      fluid_ = Teuchos::rcp(
          new FLD::TimIntOneStepTheta(discret, solver, fluidtimeparams, output, false));
    else if (timeint == INPAR::FLUID::timeint_bdf2)
      fluid_ = Teuchos::rcp(new FLD::TimIntBDF2(discret, solver, fluidtimeparams, output, false));
    else if (timeint == INPAR::FLUID::timeint_afgenalpha or
             timeint == INPAR::FLUID::timeint_npgenalpha)
      fluid_ =
          Teuchos::rcp(new FLD::TimIntGenAlpha(discret, solver, fluidtimeparams, output, false));
  }
  else
  {
    dserror("Unknown time integration for fluid\n");
  }

  // initialize algorithm for specific time-integration scheme
  fluid_->Init();

  SetInitialInflowField(fdyn);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidBaseAlgorithm::SetGeneralParameters(
    const Teuchos::RCP<Teuchos::ParameterList> fluidtimeparams,
    const Teuchos::ParameterList& prbdyn, const Teuchos::ParameterList& fdyn)
{
  fluidtimeparams->set<int>("Simple Preconditioner", INPUT::IntegralValue<int>(fdyn, "SIMPLER"));

  // -------------------------------------- number of degrees of freedom
  // number of degrees of freedom
  const int ndim = DRT::Problem::Instance()->NDim();
  fluidtimeparams->set<int>("number of velocity degrees of freedom", ndim);

  // -------------------------------------------------- time integration
  // note: here, the values are taken out of the problem-dependent ParameterList prbdyn
  // (which also can be fluiddyn itself!)

  // the default time step size
  fluidtimeparams->set<double>("time step size", prbdyn.get<double>("TIMESTEP"));
  // maximum simulation time
  fluidtimeparams->set<double>("total time", prbdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  fluidtimeparams->set<int>("max number timesteps", prbdyn.get<int>("NUMSTEP"));
  // sublist for adaptive time stepping
  fluidtimeparams->sublist("TIMEADAPTIVITY") = fdyn.sublist("TIMEADAPTIVITY");

  // -------- additional parameters in list for generalized-alpha scheme
  // parameter alpha_M
  fluidtimeparams->set<double>("alpha_M", fdyn.get<double>("ALPHA_M"));
  // parameter alpha_F
  fluidtimeparams->set<double>("alpha_F", fdyn.get<double>("ALPHA_F"));
  // parameter gamma
  fluidtimeparams->set<double>("gamma", fdyn.get<double>("GAMMA"));

  // ---------------------------------------------- nonlinear iteration
  // type of predictor
  fluidtimeparams->set<std::string>("predictor", fdyn.get<std::string>("PREDICTOR"));
  // set linearisation scheme
  fluidtimeparams->set<int>(
      "Linearisation", INPUT::IntegralValue<INPAR::FLUID::LinearisationAction>(fdyn, "NONLINITER"));
  // maximum number of nonlinear iteration steps
  fluidtimeparams->set<int>("max nonlin iter steps", fdyn.get<int>("ITEMAX"));
  // maximum number of nonlinear iteration steps for initial stationary solution
  fluidtimeparams->set<int>("max nonlin iter steps init stat sol", fdyn.get<int>("INITSTATITEMAX"));

  // parameter list containing the nonlinear solver tolerances
  const Teuchos::ParameterList& nonlinsolvertolerances =
      fdyn.sublist("NONLINEAR SOLVER TOLERANCES");

  // stop nonlinear iteration when the velocity residual is below this tolerance
  fluidtimeparams->set<double>(
      "velocity residual tolerance", nonlinsolvertolerances.get<double>("TOL_VEL_RES"));
  // stop nonlinear iteration when the pressure residual is below this tolerance
  fluidtimeparams->set<double>(
      "pressure residual tolerance", nonlinsolvertolerances.get<double>("TOL_PRES_RES"));
  // stop nonlinear iteration when the relative velocity increment is below this tolerance
  fluidtimeparams->set<double>(
      "velocity increment tolerance", nonlinsolvertolerances.get<double>("TOL_VEL_INC"));
  // stop nonlinear iteration when the relative pressure increment is below this tolerance
  fluidtimeparams->set<double>(
      "pressure increment tolerance", nonlinsolvertolerances.get<double>("TOL_PRES_INC"));

  // set convergence check
  fluidtimeparams->set<std::string>("CONVCHECK", fdyn.get<std::string>("CONVCHECK"));
  // set recomputation of residual after solution has convergenced
  fluidtimeparams->set<bool>(
      "INCONSISTENT_RESIDUAL", INPUT::IntegralValue<int>(fdyn, "INCONSISTENT_RESIDUAL") == 1);
  // set solver for L2 projection of gradients for reconstruction of consistent residual
  fluidtimeparams->set<int>("VELGRAD_PROJ_SOLVER", fdyn.get<int>("VELGRAD_PROJ_SOLVER"));
  // set adaptive linear solver tolerance
  fluidtimeparams->set<bool>("ADAPTCONV", INPUT::IntegralValue<int>(fdyn, "ADAPTCONV") == 1);
  fluidtimeparams->set<double>("ADAPTCONV_BETTER", fdyn.get<double>("ADAPTCONV_BETTER"));
  fluidtimeparams->set<bool>(
      "INFNORMSCALING", (INPUT::IntegralValue<int>(fdyn, "INFNORMSCALING") == 1));

  // ----------------------------------------------- restart and output
  const Teuchos::ParameterList& ioflags = DRT::Problem::Instance()->IOParams();
  // restart
  fluidtimeparams->set<int>("write restart every", prbdyn.get<int>("RESTARTEVRY"));
  // solution output
  fluidtimeparams->set<int>("write solution every", prbdyn.get<int>("RESULTSEVRY"));
  // flag for writing stresses
  fluidtimeparams->set<int>("write stresses", INPUT::IntegralValue<int>(ioflags, "FLUID_STRESS"));
  // flag for writing wall shear stress
  fluidtimeparams->set<int>(
      "write wall shear stresses", INPUT::IntegralValue<int>(ioflags, "FLUID_WALL_SHEAR_STRESS"));
  // flag for writing element data in every step and not only once (i.e. at step == 0 or step ==
  // upres)
  fluidtimeparams->set<int>("write element data in every step",
      INPUT::IntegralValue<int>(ioflags, "FLUID_ELEDATA_EVRY_STEP"));
  // flag for writing node data in the first time step
  fluidtimeparams->set<int>("write node data in first step",
      INPUT::IntegralValue<int>(ioflags, "FLUID_NODEDATA_FIRST_STEP"));
  // flag for writing fluid field to gmsh
  if (INPUT::IntegralValue<bool>(DRT::Problem::Instance()->IOParams(), "OUTPUT_GMSH") == false)
  {
    fluidtimeparams->set<bool>("GMSH_OUTPUT", false);
    if (INPUT::IntegralValue<bool>(fdyn, "GMSH_OUTPUT") == true)
      std::cout << "WARNING! Conflicting GMSH parameter in IO and fluid sections. No GMSH output "
                   "is written!"
                << std::endl;
  }
  else
    fluidtimeparams->set<bool>("GMSH_OUTPUT", INPUT::IntegralValue<bool>(fdyn, "GMSH_OUTPUT"));
  // flag for computing divergence
  fluidtimeparams->set<bool>("COMPUTE_DIVU", INPUT::IntegralValue<bool>(fdyn, "COMPUTE_DIVU"));
  // flag for computing kinetix energy
  fluidtimeparams->set<bool>("COMPUTE_EKIN", INPUT::IntegralValue<bool>(fdyn, "COMPUTE_EKIN"));
  // flag for computing lift and drag values
  fluidtimeparams->set<bool>("LIFTDRAG", INPUT::IntegralValue<bool>(fdyn, "LIFTDRAG"));

  // -------------------------------------------------- Oseen advection
  // set function number of given Oseen advective field
  fluidtimeparams->set<int>("OSEENFIELDFUNCNO", fdyn.get<int>("OSEENFIELDFUNCNO"));

  // ---------------------------------------------------- lift and drag
  fluidtimeparams->set<int>("liftdrag", INPUT::IntegralValue<int>(fdyn, "LIFTDRAG"));

  // -----------evaluate error for test flows with analytical solutions
  INPAR::FLUID::InitialField initfield =
      INPUT::IntegralValue<INPAR::FLUID::InitialField>(fdyn, "INITIALFIELD");
  fluidtimeparams->set<int>("eval err for analyt sol", initfield);

  // ------------------------------------------ form of convective term
  fluidtimeparams->set<std::string>("form of convective term", fdyn.get<std::string>("CONVFORM"));

  // -------------------------- potential nonlinear boundary conditions
  fluidtimeparams->set<std::string>(
      "Nonlinear boundary conditions", fdyn.get<std::string>("NONLINEARBC"));

  // ------------------------------------ potential reduced_D 3D coupling method
  fluidtimeparams->set<std::string>(
      "Strong 3D_redD coupling", fdyn.get<std::string>("STRONG_REDD_3D_COUPLING_TYPE"));

  //--------------------------------------  mesh tying for fluid
  fluidtimeparams->set<int>(
      "MESHTYING", INPUT::IntegralValue<INPAR::FLUID::MeshTying>(fdyn, "MESHTYING"));

  fluidtimeparams->set<bool>("ALLDOFCOUPLED", INPUT::IntegralValue<bool>(fdyn, "ALLDOFCOUPLED"));

  //--------------------------------------analytical error evaluation
  fluidtimeparams->set<int>("calculate error", Teuchos::getIntegralValue<int>(fdyn, "CALCERROR"));
  fluidtimeparams->set<int>("error function number", fdyn.get<int>("CALCERRORFUNCNO"));

  // -----------------------sublist containing stabilization parameters
  fluidtimeparams->sublist("RESIDUAL-BASED STABILIZATION") =
      fdyn.sublist("RESIDUAL-BASED STABILIZATION");
  fluidtimeparams->sublist("EDGE-BASED STABILIZATION") = fdyn.sublist("EDGE-BASED STABILIZATION");

  // -----------------------------get also scatra stabilization sublist
  const Teuchos::ParameterList& scatradyn =
      DRT::Problem::Instance()->ScalarTransportDynamicParams();
  fluidtimeparams->sublist("SCATRA STABILIZATION") = scatradyn.sublist("STABILIZATION");

  // --------------------------sublist containing turbulence parameters
  {
    fluidtimeparams->sublist("TURBULENCE MODEL") = fdyn.sublist("TURBULENCE MODEL");
    fluidtimeparams->sublist("SUBGRID VISCOSITY") = fdyn.sublist("SUBGRID VISCOSITY");
    fluidtimeparams->sublist("MULTIFRACTAL SUBGRID SCALES") =
        fdyn.sublist("MULTIFRACTAL SUBGRID SCALES");
    fluidtimeparams->sublist("TURBULENT INFLOW") = fdyn.sublist("TURBULENT INFLOW");
    fluidtimeparams->sublist("WALL MODEL") = fdyn.sublist("WALL MODEL");

    fluidtimeparams->sublist("TURBULENCE MODEL")
        .set<std::string>(
            "statistics outfile", DRT::Problem::Instance()->OutputControlFile()->FileName());
  }

  // ---------------------------parallel evaluation
  fluidtimeparams->set<bool>(
      "OFF_PROC_ASSEMBLY", INPUT::IntegralValue<int>(fdyn, "OFF_PROC_ASSEMBLY") == 1);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidBaseAlgorithm::CreateSecondSolver(
    const Teuchos::RCP<CORE::LINALG::Solver> solver, const Teuchos::ParameterList& fdyn)
{
  // The SIMPLER (yes,no) parameter only controls whether the fluid matrix is
  // assembled into a 2x2 blocked operator or a plain 1x1 block matrix
  // A "second solver" for the preconditioner is only needed if SIMPLER == yes
  if (INPUT::IntegralValue<int>(fdyn, "SIMPLER"))
  {
    const int linsolvernumber = fdyn.get<int>("LINEAR_SOLVER");
    const auto prec = Teuchos::getIntegralValue<INPAR::SOLVER::PreconditionerType>(
        DRT::Problem::Instance()->SolverParams(linsolvernumber), "AZPREC");
    switch (prec)
    {
      case INPAR::SOLVER::PreconditionerType::cheap_simple:
      {
        // add Inverse1 block for velocity dofs
        // tell Inverse1 block about NodalBlockInformation
        // In contrary to contact/meshtying problems this is necessary here, since we originally
        // have built the null space for the whole problem (velocity and pressure dofs). However, if
        // we split the matrix into velocity and pressure block, we have to adapt the null space
        // information for the subblocks. Therefore we need the nodal block information in the first
        // subblock for the velocities. The pressure null space is trivial to be built using a
        // constant vector
        Teuchos::ParameterList& inv1 =
            solver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1");
        inv1.sublist("NodalBlockInformation") = solver->Params().sublist("NodalBlockInformation");

        // CheapSIMPLE is somewhat hardwired here
        solver->Params().sublist("CheapSIMPLE Parameters").set("Prec Type", "CheapSIMPLE");
        solver->Params().set("FLUID", true);
      }
      break;
      case INPAR::SOLVER::PreconditionerType::multigrid_muelu_fluid:
      {
        // add Inverse1 block for velocity dofs
        // tell Inverse1 block about NodalBlockInformation
        // In contrary to contact/meshtying problems this is necessary here, since we originally
        // have built the null space for the whole problem (velocity and pressure dofs). However, if
        // we split the matrix into velocity and pressure block, we have to adapt the null space
        // information for the subblocks. Therefore we need the nodal block information in the first
        // subblock for the velocities. The pressure null space is trivial to be built using a
        // constant vector
        solver->Params().sublist("MueLu (Fluid) Parameters").sublist("NodalBlockInformation") =
            solver->Params().sublist("NodalBlockInformation");
      }
      break;
      default:
        dserror(
            "If SIMPLER flag is set to YES you can only use CheapSIMPLE as preconditioners in your "
            "fluid solver. Choose CheapSIMPLE in the SOLVER %i block in your dat file. "
            "Alternatively you can also try a multigrid block preconditioner. Use then "
            "\"MueLu_fluid\" as preconditioner and provide a parameter xml file.",
            linsolvernumber);
        break;
    }
  }

  return;
}

BACI_NAMESPACE_CLOSE

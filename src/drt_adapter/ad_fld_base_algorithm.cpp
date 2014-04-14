/*----------------------------------------------------------------------*/
/*!
\file ad_fld_base_algorithm.cpp

\brief Fluid Base Algorithm

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
 */
/*----------------------------------------------------------------------*/


#include "ad_fld_base_algorithm.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_solver.H"
#include "../drt_inpar/inpar_elch.H"
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_inpar/inpar_combust.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_poroelast.H"
#include "../drt_inpar/inpar_topopt.H"
#include "../drt_fluid/drt_periodicbc.H"
#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid/fluid_timint_loma_genalpha.H"
#include "../drt_fluid/fluid_timint_loma_ost.H"
#include "../drt_fluid/fluid_timint_loma_bdf2.H"
#include "../drt_fluid/fluid_timint_poro_ost.H"
#include "../drt_fluid/fluid_timint_poro_stat.H"
#include "../drt_fluid/fluid_timint_topopt_genalpha.H"
#include "../drt_fluid/fluid_timint_topopt_bdf2.H"
#include "../drt_fluid/fluid_timint_topopt_ost.H"
#include "../drt_fluid/fluid_timint_topopt_stat.H"
#include "../drt_fluid/fluid_timint_red_genalpha.H"
#include "../drt_fluid/fluid_timint_red_bdf2.H"
#include "../drt_fluid/fluid_timint_red_ost.H"
#include "../drt_fluid/fluid_timint_red_stat.H"
#include "../drt_fluid_xfluid/xfluid.H"
#include "../drt_fluid_xfluid/xfluidfluid.H"
#include "../drt_combust/combust_fluidimplicitintegration.H"
#include "../linalg/linalg_solver.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "ad_fld_fluid_fluid_fsi.H"
#include "ad_fld_fluid_fpsi.H"
#include "ad_fld_fluid_fsi.H"
#include "ad_fld_lung.H"
#include "ad_fld_poro.H"
#include "ad_fld_xfluid_fsi.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidBaseAlgorithm::FluidBaseAlgorithm(
  const Teuchos::ParameterList& prbdyn,
  const Teuchos::ParameterList& fdyn,
  const std::string& disname,
  bool isale,
  bool init)
{
  SetupFluid(prbdyn, fdyn, disname, isale, init);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidBaseAlgorithm::FluidBaseAlgorithm(
  const Teuchos::ParameterList& prbdyn,
  const Teuchos::RCP<DRT::Discretization> discret)
{
  SetupInflowFluid(prbdyn, discret);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidBaseAlgorithm::~FluidBaseAlgorithm()
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidBaseAlgorithm::SetupFluid(
  const Teuchos::ParameterList& prbdyn,
  const Teuchos::ParameterList& fdyn,
  const std::string& disname,
  bool isale,
  bool init)
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("ADAPTER::FluidBaseAlgorithm::SetupFluid");
  Teuchos::TimeMonitor monitor(*t);

  // -------------------------------------------------------------------
  // what's the current problem type?
  // -------------------------------------------------------------------
  PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> actdis
    = DRT::Problem::Instance()->GetDis(disname);

  // -------------------------------------------------------------------
  // connect degrees of freedom for periodic boundary conditions
  // -------------------------------------------------------------------
  if((probtype != prb_fsi) and
     (probtype != prb_combust))
  {
    PeriodicBoundaryConditions pbc(actdis);
    pbc.UpdateDofsForPeriodicBoundaryConditions();
  }

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->HaveDofs())
  {
    if (probtype == prb_fsi_xfem or
        probtype == prb_fluid_xfem or
        probtype == prb_combust or
        probtype == prb_fsi_crack)
    {
      actdis->FillComplete(false,false,false);
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
  if (probtype != prb_combust and
      probtype != prb_fluid_fluid and
      probtype != prb_fluid_fluid_ale and
      probtype != prb_fluid_fluid_fsi)
  {
    output->WriteMesh(0,0.0);
  }
  // in case Init() is not called in this, we have to be carefull with
  // the order fields are written into the control file
  // this order has to be respected by WriteMesh
  // if problems with post_drt_ensight occur, a possible solution might be
  // to move output->WriteMesh(0,0.0) from the Init()-function to the
  // contructor of the xfluidfluid-class (see also combust algorithm)
  if (init == false and (probtype == prb_fluid_fluid or
      probtype == prb_fluid_fluid_ale or
      probtype == prb_fluid_fluid_fsi))
    dserror("Read remark!");

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  //const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  //const Teuchos::ParameterList& probsize = DRT::Problem::Instance()->ProblemSizeParams();
  //const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();

  if (actdis->Comm().MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(IO::cout, fdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::null;

  switch(DRT::INPUT::IntegralValue<INPAR::FLUID::MeshTying>(fdyn,"MESHTYING"))
  {
    case INPAR::FLUID::condensed_bmat:
    {
      // FIXME: The solver should not be taken from the contact dynamic section here,
      // but must be specified in the fluid dynamic section instead (popp 11/2012)

      const Teuchos::ParameterList& mshparams = DRT::Problem::Instance()->ContactDynamicParams();
      const int mshsolver = mshparams.get<int>("LINEAR_SOLVER");        // meshtying solver (with block preconditioner, e.g. BGS 2x2)
      const int fluidsolver = fdyn.get<int>("LINEAR_SOLVER");           // fluid solver
      const int fluidpressuresolver = fdyn.get<int>("SIMPLER_SOLVER");  // fluid pressure solver
      if (mshsolver == (-1))
        dserror("no linear solver defined for fluid meshtying problem. Please set LINEAR_SOLVER in CONTACT DYNAMIC to a valid number!");
      if (fluidsolver == (-1))
        dserror("no linear solver defined for fluid meshtying problem. Please set LINEAR_SOLVER in FLUID DYNAMIC to a valid number! This solver is used within block preconditioner (e.g. BGS2x2) as \"Inverse 1\".");
      if (fluidpressuresolver == (-1))
        dserror("no linear solver defined for fluid meshtying problem. Please set SIMPLER_SOLVER in FLUID DYNAMIC to a valid number! This solver is used within block preconditioner (e.g. BGS2x2) as \"Inverse 2\".");

      // check, if meshtying solver is used with a valid block preconditioner
      const int azprectype
        = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(
            DRT::Problem::Instance()->SolverParams(mshsolver),
            "AZPREC"
            );

      // plausibility check
      switch (azprectype)
      {
        case INPAR::SOLVER::azprec_CheapSIMPLE:
        case INPAR::SOLVER::azprec_BGS2x2:      // block preconditioners, that are implemented in BACI
          break;
        case INPAR::SOLVER::azprec_BGSnxn:      // block preconditioners from Teko
        case INPAR::SOLVER::azprec_TekoSIMPLE:
        {
    #ifdef HAVE_TEKO
          // check if structural solver and thermal solver are Stratimikos based (Teko expects stratimikos)
          int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(DRT::Problem::Instance()->SolverParams(fluidsolver), "SOLVER");
          if (solvertype != INPAR::SOLVER::stratimikos_amesos &&
              solvertype != INPAR::SOLVER::stratimikos_aztec  &&
              solvertype != INPAR::SOLVER::stratimikos_belos)
          dserror("Teko expects a STRATIMIKOS solver object in SOLVER %i", fluidsolver);

          solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(DRT::Problem::Instance()->SolverParams(fluidpressuresolver), "SOLVER");
          if (solvertype != INPAR::SOLVER::stratimikos_amesos &&
              solvertype != INPAR::SOLVER::stratimikos_aztec  &&
              solvertype != INPAR::SOLVER::stratimikos_belos)
            dserror("Teko expects a STRATIMIKOS solver object in SOLVER %i",fluidpressuresolver);
    #else
          dserror("Teko preconditioners only available with HAVE_TEKO flag for TRILINOS_DEV (>Q1/2011)");
    #endif
        }
        break;
        default:
              dserror("Block Gauss-Seidel BGS2x2 preconditioner expected for fluid meshtying problem. Please set AZPREC to BGS2x2 in solver block %i",mshsolver);
              break;
      }

      // create solver objects
      solver =
        Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(mshsolver),
                               actdis->Comm(),
                               DRT::Problem::Instance()->ErrorFile()->Handle()));

      // add sub block solvers/smoothers to block preconditioners
      switch (azprectype)
      {
        case INPAR::SOLVER::azprec_CheapSIMPLE:
          break; // CheapSIMPLE adds its own Inverse1 and Inverse2 blocks
        case INPAR::SOLVER::azprec_BGS2x2:      // block preconditioners, that are implemented in BACI
        case INPAR::SOLVER::azprec_BGSnxn:      // block preconditioners from Teko
        case INPAR::SOLVER::azprec_TekoSIMPLE:
        {
          // set Inverse blocks for block preconditioner
          // for BGS preconditioner
          // This is only necessary for BGS. CheapSIMPLE has a more modern framework
          solver->PutSolverParamsToSubParams("Inverse1",
              DRT::Problem::Instance()->SolverParams(fluidsolver));

          solver->PutSolverParamsToSubParams("Inverse2",
              DRT::Problem::Instance()->SolverParams(fluidpressuresolver));
        }
        break;
        default:
              dserror("Block Gauss-Seidel BGS2x2 preconditioner expected for fluid meshtying problem. Please set AZPREC to BGS2x2 in solver block %i",mshsolver);
              break;
      }

      solver->Params().set<bool>("MESHTYING",true);   // mark it as meshtying problem
    }
    break;
    case INPAR::FLUID::condensed_smat:
    case INPAR::FLUID::condensed_bmat_merged:
    {
      // meshtying (no saddle point problem)
      const Teuchos::ParameterList& mshparams = DRT::Problem::Instance()->ContactDynamicParams();
      const int mshsolver = mshparams.get<int>("LINEAR_SOLVER");             // meshtying solver (with block preconditioner, e.g. BGS 2x2)
      if (mshsolver == (-1))
        dserror("no linear solver defined for fluid meshtying problem. Please set LINEAR_SOLVER in CONTACT DYNAMIC to a valid number!");

      solver =
        Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(mshsolver),
                               actdis->Comm(),
                               DRT::Problem::Instance()->ErrorFile()->Handle()));
    }
    break;
    case INPAR::FLUID::no_meshtying: // no meshtying -> use FLUID SOLVER
    default:
    {
      // default: create solver using the fluid solver params from FLUID SOLVER block

      // get the solver number used for linear fluid solver
      const int linsolvernumber = fdyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == (-1))
        dserror("no linear solver defined for fluid problem. Please set LINEAR_SOLVER in FLUID DYNAMIC to a valid number!");
      solver =
        Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                               actdis->Comm(),
                               DRT::Problem::Instance()->ErrorFile()->Handle()));

      break;
    }
  }

  // compute null space information
  if (probtype != prb_fsi_xfem and
      probtype != prb_fluid_xfem and
      probtype != prb_combust and
      probtype != prb_fluid_fluid and
      probtype != prb_fluid_fluid_ale and
      probtype != prb_fluid_fluid_fsi and
      probtype != prb_fsi_crack )
  {
    switch(DRT::INPUT::IntegralValue<int>(fdyn,"MESHTYING"))
    {
      // switch types
      case INPAR::FLUID::condensed_bmat:
      {
        const Teuchos::ParameterList& mshparams = DRT::Problem::Instance()->ContactDynamicParams();
        const int mshsolver = mshparams.get<int>("LINEAR_SOLVER");        // meshtying solver (with block preconditioner, e.g. BGS 2x2)

        // check, if meshtying solver is used with a valid block preconditioner
        const int azprectype
          = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(
              DRT::Problem::Instance()->SolverParams(mshsolver),
              "AZPREC"
              );

        switch (azprectype)
        {
          // block preconditioners, that are implemented in BACI
          case INPAR::SOLVER::azprec_CheapSIMPLE:
          {
            actdis->ComputeNullSpaceIfNecessary(solver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1"),true);
            actdis->ComputeNullSpaceIfNecessary(solver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse2"),true);
          }
          break;
          case INPAR::SOLVER::azprec_BGS2x2:
          case INPAR::SOLVER::azprec_BGSnxn:      // block preconditioners from Teko
          case INPAR::SOLVER::azprec_TekoSIMPLE:
          {
            actdis->ComputeNullSpaceIfNecessary(solver->Params().sublist("Inverse1"),true);
            actdis->ComputeNullSpaceIfNecessary(solver->Params().sublist("Inverse2"),true);
          }
          break;
        } // end switch azprectype
      }
      break;
      default:
        // no block matrix
        actdis->ComputeNullSpaceIfNecessary(solver->Params(),true);
        break;
    }
  }

  // create a second solver for SIMPLER preconditioner if chosen from input
  CreateSecondSolver(solver,fdyn);

  // -------------------------------------------------------------------
  // set parameters in list
  // -------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> fluidtimeparams = Teuchos::rcp(new Teuchos::ParameterList());

  // physical type of fluid flow (incompressible, Boussinesq Approximation, varying density, loma, poro)
  fluidtimeparams->set<int>("Physical Type",
      DRT::INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn,"PHYSICAL_TYPE"));
  // and  check correct setting
  if ((probtype == prb_thermo_fsi or probtype == prb_loma) and
      DRT::INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn,"PHYSICAL_TYPE")
      != INPAR::FLUID::loma)
    dserror("Input parameter PHYSICAL_TYPE in section FLUID DYNAMIC needs to be 'Loma' for low-Mach-number flow and Thermo-fluid-structure interaction!");
  if ( (probtype == prb_poroelast or probtype == prb_poroscatra or (probtype == prb_fpsi and disname == "porofluid") ) and
      DRT::INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn,"PHYSICAL_TYPE")!= INPAR::FLUID::poro and
      DRT::INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn,"PHYSICAL_TYPE")!= INPAR::FLUID::poro_p1 and
      DRT::INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn,"PHYSICAL_TYPE")!= INPAR::FLUID::poro_p2 )
    dserror("Input parameter PHYSICAL_TYPE in section FLUID DYNAMIC needs to be 'Poro' or 'Poro_P1' for poro-elasticity!");

  // now, set general parameters required for all problems
  SetGeneralParameters(fluidtimeparams,prbdyn,fdyn);

  // and, finally, add problem specific parameters

  // add some loma specific parameters
  // get also scatra stabilization sublist
  const Teuchos::ParameterList& lomadyn =
    DRT::Problem::Instance()->LOMAControlParams();
  fluidtimeparams->sublist("LOMA").set<bool>("update material",DRT::INPUT::IntegralValue<int>(lomadyn,"SGS_MATERIAL_UPDATE"));

  // ----------------------------- sublist for general xfem-specific parameters
  if (   probtype == prb_fluid_xfem
      or probtype == prb_fsi_xfem
      or probtype == prb_fluid_fluid
      or probtype == prb_fluid_fluid_ale
      or probtype == prb_fluid_fluid_fsi
      or probtype == prb_fsi_crack )
  {
    // get also scatra stabilization sublist
    const Teuchos::ParameterList& xdyn =
      DRT::Problem::Instance()->XFEMGeneralParams();

    fluidtimeparams->sublist("XFEM").set<int>("MAX_NUM_DOFSETS", xdyn.get<int>("MAX_NUM_DOFSETS"));
    fluidtimeparams->sublist("XFEM").set<string>("VOLUME_GAUSS_POINTS_BY", xdyn.get<std::string>("VOLUME_GAUSS_POINTS_BY"));
    fluidtimeparams->sublist("XFEM").set<string>("BOUNDARY_GAUSS_POINTS_BY", xdyn.get<std::string>("BOUNDARY_GAUSS_POINTS_BY"));

    // GMSH solution output
    fluidtimeparams->sublist("XFEM").set<int>("GMSH_DEBUG_OUT",        DRT::INPUT::IntegralValue<int>(xdyn, "GMSH_DEBUG_OUT"));
    fluidtimeparams->sublist("XFEM").set<int>("GMSH_DEBUG_OUT_SCREEN", DRT::INPUT::IntegralValue<int>(xdyn, "GMSH_DEBUG_OUT_SCREEN"));
    fluidtimeparams->sublist("XFEM").set<int>("GMSH_EOS_OUT",          DRT::INPUT::IntegralValue<int>(xdyn, "GMSH_EOS_OUT"));
    fluidtimeparams->sublist("XFEM").set<int>("GMSH_SOL_OUT",          DRT::INPUT::IntegralValue<int>(xdyn, "GMSH_SOL_OUT"));
    fluidtimeparams->sublist("XFEM").set<int>("GMSH_DISCRET_OUT",      DRT::INPUT::IntegralValue<int>(xdyn, "GMSH_DISCRET_OUT"));
    fluidtimeparams->sublist("XFEM").set<int>("GMSH_CUT_OUT",          DRT::INPUT::IntegralValue<int>(xdyn, "GMSH_CUT_OUT"));
  }

  // ----------------------------- sublist for xfem-specific fluid parameters
  if (   probtype == prb_fluid_xfem
      or probtype == prb_fsi_xfem
      or probtype == prb_fluid_fluid
      or probtype == prb_fluid_fluid_ale
      or probtype == prb_fluid_fluid_fsi
      or probtype == prb_fsi_crack )
  {

    const Teuchos::ParameterList& xfdyn     = DRT::Problem::Instance()->XFluidDynamicParams();

    fluidtimeparams->sublist("XFLUID DYNAMIC/GENERAL")       = xfdyn.sublist("GENERAL");
    fluidtimeparams->sublist("XFLUID DYNAMIC/STABILIZATION") = xfdyn.sublist("STABILIZATION");

    fluidtimeparams->sublist("XFLUID DYNAMIC/GENERAL").set<string>("MONOLITHIC_XFFSI_APPROACH",xfdyn.sublist("GENERAL").get<std::string>("MONOLITHIC_XFFSI_APPROACH"));
    fluidtimeparams->sublist("XFLUID DYNAMIC/GENERAL").set<string>("XFLUIDFLUID_TIMEINT",xfdyn.sublist("GENERAL").get<std::string>("XFLUIDFLUID_TIMEINT"));
    fluidtimeparams->sublist("XFLUID DYNAMIC/GENERAL").set<double>("XFLUIDFLUID_SEARCHRADIUS",  xfdyn.sublist("GENERAL").get<double>("XFLUIDFLUID_SEARCHRADIUS"));
  }


  // sublist for combustion-specific fluid parameters
  /* This sublist COMBUSTION FLUID contains parameters for the fluid field
   * which are only relevant for a combustion problem.                 07/08 henke */
  if (probtype == prb_combust)
  {
    fluidtimeparams->sublist("COMBUSTION FLUID")=prbdyn.sublist("COMBUSTION FLUID");
    // parameter COMBUSTTYPE from sublist COMBUSTION FLUID is also added to sublist XFEM
    fluidtimeparams->sublist("XFEM").set<int>("combusttype", DRT::INPUT::IntegralValue<INPAR::COMBUST::CombustionType>(prbdyn.sublist("COMBUSTION FLUID"),"COMBUSTTYPE"));
  }

  if (probtype == prb_fluid_topopt)
    fluidtimeparams->set<int>("opti testcase",DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(prbdyn.sublist("TOPOLOGY OPTIMIZER"),"TESTCASE"));

  // -------------------------------------------------------------------
  // additional parameters and algorithm call depending on respective
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  INPAR::FLUID::TimeIntegrationScheme timeint = DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fdyn,"TIMEINTEGR");

  // sanity checks and default flags
  if (probtype == prb_fsi or
      probtype == prb_fsi_lung or
      probtype == prb_gas_fsi or
      probtype == prb_biofilm_fsi or
      probtype == prb_thermo_fsi or
      probtype == prb_fluid_fluid_fsi or
      probtype == prb_fsi_xfem or
      probtype == prb_fsi_crack or
      probtype == prb_fsi_redmodels)
  {
    // in case of FSI calculations we do not want a stationary fluid solver
    if (timeint == INPAR::FLUID::timeint_stationary)
      dserror("Stationary fluid solver not allowed for FSI.");

    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

    fluidtimeparams->set<bool>("interface second order", DRT::INPUT::IntegralValue<int>(fsidyn,"SECONDORDER"));
    fluidtimeparams->set<bool>("shape derivatives",
                               DRT::INPUT::IntegralValue<int>(fsimono,"SHAPEDERIVATIVES"));

    const int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");

    if (coupling == fsi_iter_lung_monolithicstructuresplit or
        coupling == fsi_iter_lung_monolithicfluidsplit or
        coupling == fsi_iter_constr_monolithicstructuresplit or
        coupling == fsi_iter_constr_monolithicfluidsplit or
        coupling == fsi_iter_fluidfluid_monolithicstructuresplit or
        coupling == fsi_iter_fluidfluid_monolithicfluidsplit or
        coupling == fsi_iter_fluidfluid_monolithicstructuresplit_nox or
        coupling == fsi_iter_xfem_monolithic)
    {
      // No explicit predictor for these monolithic FSI schemes, yet.
      // Check, whether fluid predictor is 'steady_state'. Otherwise, throw
      // an error.
      if (fluidtimeparams->get<std::string>("predictor") != "steady_state")
        dserror("No fluid predictor allowed for current monolithic FSI scheme, yet. Use 'steady_state', instead!");
    }
  }
  if (probtype == prb_freesurf)
  {
    // in case of FSI calculations we do not want a stationary fluid solver
    if (timeint == INPAR::FLUID::timeint_stationary)
      dserror("Stationary fluid solver not allowed for free surface problem.");

    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

    fluidtimeparams->set<bool>("interface second order", DRT::INPUT::IntegralValue<int>(fsidyn,"SECONDORDER"));
    fluidtimeparams->set<bool>("shape derivatives",
                               DRT::INPUT::IntegralValue<int>(fsimono,"SHAPEDERIVATIVES"));

    const int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
    if (coupling == fsi_iter_monolithicfluidsplit or
        coupling == fsi_iter_monolithicstructuresplit)
    {
      // No explicit predictor for monolithic free surface flow schemes, yet.
      // Check, whether fluid predictor is 'steady_state'. Otherwise, throw
      // an error.
      if (fluidtimeparams->get<std::string>("predictor") != "steady_state")
        dserror("No fluid predictor allowed for current monolithic free surface scheme, yet. Use 'steady_state', instead!");
    }
  }

  // sanity checks and default flags
  if ( probtype == prb_fluid_xfem )
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    fluidtimeparams->set<bool>("interface second order", DRT::INPUT::IntegralValue<int>(fsidyn,"SECONDORDER"));
  }

  // sanity checks and default flags
  if ( probtype == prb_fsi_xfem or probtype == prb_fsi_crack )
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

    const int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");

    if (coupling == fsi_iter_monolithicfluidsplit or
        coupling == fsi_iter_monolithicstructuresplit)
    {
      // there are a couple of restrictions in monolithic FSI
      dserror("for XFSI there is no monolithicfluidsplit or monolithicstructuresplit, use monolithicxfem or any partitioned algorithm instead");
    }
  }

  // sanity checks and default flags
  if ( probtype == prb_fluid_xfem or probtype == prb_fsi_xfem or probtype == prb_fsi_crack )
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    const int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
    fluidtimeparams->set<int>("COUPALGO", coupling);
  }

  if (probtype == prb_elch)
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    fluidtimeparams->set<bool>("interface second order", DRT::INPUT::IntegralValue<int>(fsidyn,"SECONDORDER"));
  }

  if ( probtype == prb_poroelast or
       probtype == prb_poroscatra or
      (probtype == prb_fpsi and disname == "porofluid")
     )
  {
    const Teuchos::ParameterList& porodyn = DRT::Problem::Instance()->PoroelastDynamicParams();
    fluidtimeparams->set<bool>("poroelast",true);
    fluidtimeparams->set<bool>("interface second order", DRT::INPUT::IntegralValue<int>(porodyn,"SECONDORDER"));
    fluidtimeparams->set<bool>("shape derivatives",false);
    fluidtimeparams->set<bool>("conti partial integration", DRT::INPUT::IntegralValue<int>(porodyn,"CONTIPARTINT"));
  }
  else if (probtype == prb_fpsi and disname == "fluid")
  {
    if (timeint == INPAR::FLUID::timeint_stationary)
      dserror("Stationary fluid solver not allowed for FPSI.");

    fluidtimeparams->set<bool>("interface second order", DRT::INPUT::IntegralValue<int>(prbdyn,"SECONDORDER"));
    fluidtimeparams->set<bool>("shape derivatives", DRT::INPUT::IntegralValue<int>(prbdyn,"SHAPEDERIVATIVES"));
  }

  // -------------------------------------------------------------------
  // additional parameters and algorithm call depending on respective
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  if(timeint == INPAR::FLUID::timeint_stationary or
     timeint == INPAR::FLUID::timeint_one_step_theta or
     timeint == INPAR::FLUID::timeint_bdf2 or
     timeint == INPAR::FLUID::timeint_afgenalpha or
     timeint == INPAR::FLUID::timeint_npgenalpha
    )
  {
    // -----------------------------------------------------------------
    // set additional parameters in list for
    // one-step-theta/BDF2/af-generalized-alpha/stationary scheme
    // -----------------------------------------------------------------
    // type of time-integration (or stationary) scheme
    fluidtimeparams->set<int>("time int algo",timeint);
    // parameter theta for time-integration schemes
    fluidtimeparams->set<double>           ("theta"                    ,fdyn.get<double>("THETA"));
    // number of steps for potential start algorithm
    fluidtimeparams->set<int>              ("number of start steps"    ,fdyn.get<int>("NUMSTASTEPS"));
    // parameter theta for potential start algorithm
    fluidtimeparams->set<double>           ("start theta"              ,fdyn.get<double>("START_THETA"));
    // parameter for grid velocity interpolation
    fluidtimeparams->set<int>              ("order gridvel"            ,DRT::INPUT::IntegralValue<int>(fdyn,"GRIDVEL"));

    fluidtimeparams->set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());
    bool dirichletcond = true;
    if (probtype == prb_fsi or
        probtype == prb_fsi_lung or
        probtype == prb_gas_fsi or
        probtype == prb_biofilm_fsi or
        probtype == prb_thermo_fsi or
        probtype == prb_fluid_fluid_fsi or
        probtype == prb_fsi_redmodels)
    {
      // FSI input parameters
      const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
      const int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
      if (coupling == fsi_iter_monolithicfluidsplit or
          coupling == fsi_iter_monolithicstructuresplit or
          coupling == fsi_iter_lung_monolithicstructuresplit or
          coupling == fsi_iter_lung_monolithicfluidsplit or
          coupling == fsi_iter_constr_monolithicstructuresplit or
          coupling == fsi_iter_constr_monolithicfluidsplit or
          coupling == fsi_iter_mortar_monolithicstructuresplit or
          coupling == fsi_iter_mortar_monolithicfluidsplit or
          coupling == fsi_iter_fluidfluid_monolithicstructuresplit 	or
          coupling == fsi_iter_fluidfluid_monolithicfluidsplit or
          coupling == fsi_iter_fluidfluid_monolithicstructuresplit_nox)
      {
        dirichletcond = false;
      }
    }

    if (probtype == prb_poroelast or probtype == prb_poroscatra or probtype == prb_fpsi)
      dirichletcond = false;

    //------------------------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor);
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom

    switch (probtype)
    {
    case prb_fluid:
    case prb_scatra:
    case prb_cavitation:
    {
      if(timeint == INPAR::FLUID::timeint_stationary)
        fluid_ = Teuchos::rcp(new FLD::TimIntStationary(actdis, solver, fluidtimeparams, output, isale));
      else if(timeint == INPAR::FLUID::timeint_one_step_theta)
        fluid_ = Teuchos::rcp(new FLD::TimIntOneStepTheta(actdis, solver, fluidtimeparams, output, isale));
      else if(timeint == INPAR::FLUID::timeint_bdf2)
        fluid_ = Teuchos::rcp(new FLD::TimIntBDF2(actdis, solver, fluidtimeparams, output, isale));
      else if(timeint == INPAR::FLUID::timeint_afgenalpha or
          timeint == INPAR::FLUID::timeint_npgenalpha)
        fluid_ = Teuchos::rcp(new FLD::TimIntGenAlpha(actdis, solver, fluidtimeparams, output, isale));
      else
        dserror("Unknown time integration for this fluid problem type\n");

    }
    break;
    case prb_fluid_redmodels:
    {
        if(timeint == INPAR::FLUID::timeint_stationary)
          fluid_ = Teuchos::rcp(new FLD::TimIntRedModelsStat(actdis, solver, fluidtimeparams, output, isale));
        else if(timeint == INPAR::FLUID::timeint_one_step_theta)
          fluid_ = Teuchos::rcp(new FLD::TimIntRedModelsOst(actdis, solver, fluidtimeparams, output, isale));
        else if(timeint == INPAR::FLUID::timeint_afgenalpha or
            timeint == INPAR::FLUID::timeint_npgenalpha)
          fluid_ = Teuchos::rcp(new FLD::TimIntRedModelsGenAlpha(actdis, solver, fluidtimeparams, output, isale));
        else if(timeint == INPAR::FLUID::timeint_bdf2)
          fluid_ = Teuchos::rcp(new FLD::TimIntRedModelsBDF2(actdis, solver, fluidtimeparams, output, isale));
        else
          dserror("Unknown time integration for this fluid problem type\n");

    }
    break;
    case prb_loma:
    {
      if(timeint == INPAR::FLUID::timeint_afgenalpha
         or timeint == INPAR::FLUID::timeint_npgenalpha)
        fluid_ = Teuchos::rcp(new FLD::TimIntLomaGenAlpha(actdis, solver, fluidtimeparams, output, isale));
      else if(timeint == INPAR::FLUID::timeint_one_step_theta)
        fluid_ = Teuchos::rcp(new FLD::TimIntLomaOst(actdis, solver, fluidtimeparams, output, isale));
      else if(timeint == INPAR::FLUID::timeint_bdf2)
        fluid_ = Teuchos::rcp(new FLD::TimIntLomaBDF2(actdis, solver, fluidtimeparams, output, isale));
      else
        dserror("Unknown time integration for this fluid problem type\n");
    }
      break;
    case prb_fluid_topopt:
    {
      if(timeint == INPAR::FLUID::timeint_stationary)
        fluid_ = Teuchos::rcp(new FLD::TimIntTopOptStat(actdis, solver, fluidtimeparams, output, isale));
      else if(timeint == INPAR::FLUID::timeint_one_step_theta)
        fluid_ = Teuchos::rcp(new FLD::TimIntTopOptOst(actdis, solver, fluidtimeparams, output, isale));
      else if(timeint == INPAR::FLUID::timeint_bdf2)
        fluid_ = Teuchos::rcp(new FLD::TimIntTopOptBDF2(actdis, solver, fluidtimeparams, output, isale));
      else if(timeint == INPAR::FLUID::timeint_afgenalpha
              or timeint == INPAR::FLUID::timeint_npgenalpha)
        fluid_ = Teuchos::rcp(new FLD::TimIntTopOptGenAlpha(actdis, solver, fluidtimeparams, output, isale));
      else
        dserror("Unknown time integration for this fluid problem type\n");
    }
      break;
    case prb_fsi_xfem:
    case prb_fluid_xfem:
    case prb_fsi_crack:
    {
      Teuchos::RCP<DRT::Discretization> soliddis = DRT::Problem::Instance()->GetDis("structure");
      Teuchos::RCP<FLD::XFluid> tmpfluid = Teuchos::rcp( new FLD::XFluid( actdis, soliddis, solver, fluidtimeparams, output));
      fluid_ = Teuchos::rcp(new XFluidFSI(tmpfluid, actdis, soliddis, solver, fluidtimeparams, output));
    }
    break;
    case prb_combust:
    {
      fluid_ = Teuchos::rcp(new FLD::CombustFluidImplicitTimeInt(actdis, solver, fluidtimeparams, output));
    }
    break;
    case prb_fluid_fluid_ale:
    case prb_fluid_fluid:
    {
      Teuchos::RCP<DRT::Discretization> embfluiddis  =  DRT::Problem::Instance()->GetDis("xfluid");
      bool monolithicfluidfluidfsi = false;
      Teuchos::RCP<FLD::XFluidFluid> tmpfluid = Teuchos::rcp(new FLD::XFluidFluid(actdis,embfluiddis,solver,fluidtimeparams,output,isale,monolithicfluidfluidfsi));
      fluid_ = Teuchos::rcp(new FluidFluidFSI(tmpfluid,embfluiddis,actdis,solver,fluidtimeparams,isale,dirichletcond,monolithicfluidfluidfsi));
    }
    break;
    case prb_fluid_fluid_fsi:
    {
      Teuchos::RCP<DRT::Discretization> embfluiddis  =  DRT::Problem::Instance()->GetDis("xfluid");
      const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
      const int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
      bool monolithicfluidfluidfsi;
      if(coupling == fsi_iter_fluidfluid_monolithicstructuresplit
         or coupling == fsi_iter_fluidfluid_monolithicfluidsplit
         or coupling == fsi_iter_fluidfluid_monolithicstructuresplit_nox)
        monolithicfluidfluidfsi = true;
      else
        monolithicfluidfluidfsi = false;

      Teuchos::RCP<FLD::XFluidFluid> tmpfluid = Teuchos::rcp(new FLD::XFluidFluid(actdis,embfluiddis,solver,fluidtimeparams,output,isale,monolithicfluidfluidfsi));
      fluid_ = Teuchos::rcp(new FluidFluidFSI(tmpfluid,embfluiddis,actdis,solver,fluidtimeparams,isale,dirichletcond,monolithicfluidfluidfsi));
    }
    break;
    case prb_fsi:
    case prb_immersed_fsi:
    case prb_gas_fsi:
    case prb_biofilm_fsi:
    case prb_thermo_fsi:
    case prb_fluid_ale:
    case prb_freesurf:
    { //
      Teuchos::RCP<FLD::FluidImplicitTimeInt> tmpfluid;
      if(timeint == INPAR::FLUID::timeint_stationary)
        tmpfluid = Teuchos::rcp(new FLD::TimIntStationary(actdis, solver, fluidtimeparams, output, isale));
      else if(timeint == INPAR::FLUID::timeint_one_step_theta)
        tmpfluid = Teuchos::rcp(new FLD::TimIntOneStepTheta(actdis, solver, fluidtimeparams, output, isale));
      else if(timeint == INPAR::FLUID::timeint_bdf2)
        tmpfluid = Teuchos::rcp(new FLD::TimIntBDF2(actdis, solver, fluidtimeparams, output, isale));
      else if(timeint == INPAR::FLUID::timeint_afgenalpha or
          timeint == INPAR::FLUID::timeint_npgenalpha)
        tmpfluid = Teuchos::rcp(new FLD::TimIntGenAlpha(actdis, solver, fluidtimeparams, output, isale));
      else
        dserror("Unknown time integration for this fluid problem type\n");

      fluid_ = Teuchos::rcp(new FluidFSI(tmpfluid,actdis,solver,fluidtimeparams,output,isale,dirichletcond));
    }
    break;
    case prb_fsi_redmodels:
    { //
      Teuchos::RCP<FLD::FluidImplicitTimeInt> tmpfluid;
      std::cout << "\n Warning: FSI_RedModels has never been tested. Test it! \n" << std::endl;
      if(timeint == INPAR::FLUID::timeint_stationary)
        tmpfluid = Teuchos::rcp(new FLD::TimIntRedModelsStat(actdis, solver, fluidtimeparams, output, isale));
      else if(timeint == INPAR::FLUID::timeint_one_step_theta)
        tmpfluid = Teuchos::rcp(new FLD::TimIntRedModelsOst(actdis, solver, fluidtimeparams, output, isale));
      else if(timeint == INPAR::FLUID::timeint_afgenalpha or
          timeint == INPAR::FLUID::timeint_npgenalpha)
        tmpfluid = Teuchos::rcp(new FLD::TimIntRedModelsGenAlpha(actdis, solver, fluidtimeparams, output, isale));
      else if(timeint == INPAR::FLUID::timeint_bdf2)
        tmpfluid = Teuchos::rcp(new FLD::TimIntRedModelsBDF2(actdis, solver, fluidtimeparams, output, isale));
      else
        dserror("Unknown time integration for this fluid problem type\n");
      fluid_ = Teuchos::rcp(new FluidFSI(tmpfluid,actdis,solver,fluidtimeparams,output,isale,dirichletcond));
    }
    break;
    case prb_fsi_lung:
    {
      Teuchos::RCP<FLD::FluidImplicitTimeInt> tmpfluid;
      if(timeint == INPAR::FLUID::timeint_stationary)
        tmpfluid = Teuchos::rcp(new FLD::TimIntRedModelsStat(actdis, solver, fluidtimeparams, output, isale));
      else if(timeint == INPAR::FLUID::timeint_one_step_theta)
        tmpfluid = Teuchos::rcp(new FLD::TimIntRedModelsOst(actdis, solver, fluidtimeparams, output, isale));
      else if(timeint == INPAR::FLUID::timeint_afgenalpha or
          timeint == INPAR::FLUID::timeint_npgenalpha)
        tmpfluid = Teuchos::rcp(new FLD::TimIntRedModelsGenAlpha(actdis, solver, fluidtimeparams, output, isale));
      else if(timeint == INPAR::FLUID::timeint_bdf2)
        tmpfluid = Teuchos::rcp(new FLD::TimIntRedModelsBDF2(actdis, solver, fluidtimeparams, output, isale));
      else
        dserror("Unknown time integration for this fluid problem type\n");
      fluid_ = Teuchos::rcp(new FluidLung(tmpfluid,actdis,solver,fluidtimeparams,output,isale,dirichletcond));
    }
    break;
    case prb_poroelast:
    case prb_poroscatra:
    case prb_fpsi:
    {
      Teuchos::RCP<FLD::FluidImplicitTimeInt> tmpfluid;
      if(disname == "porofluid")
      {
        if(timeint == INPAR::FLUID::timeint_stationary)
          tmpfluid = Teuchos::rcp(new FLD::TimIntPoroStat(actdis, solver, fluidtimeparams, output, isale));
        else if(timeint == INPAR::FLUID::timeint_one_step_theta)
          tmpfluid = Teuchos::rcp(new FLD::TimIntPoroOst(actdis, solver, fluidtimeparams, output, isale));
        else
          dserror("Unknown time integration for this fluid problem type\n");
        fluid_ = Teuchos::rcp(new FluidPoro(tmpfluid, actdis, solver, fluidtimeparams, output, isale, dirichletcond));
      }
      else if(disname == "fluid")
      {
        fluidtimeparams->set<int>("Physical Type",INPAR::FLUID::incompressible);
        if(timeint == INPAR::FLUID::timeint_stationary)
          tmpfluid = Teuchos::rcp(new FLD::TimIntPoroStat(actdis, solver, fluidtimeparams, output, isale));
        else if(timeint == INPAR::FLUID::timeint_one_step_theta)
          tmpfluid = Teuchos::rcp(new FLD::TimIntPoroOst(actdis, solver, fluidtimeparams, output, isale));
        else
          dserror("Unknown time integration for this fluid problem type\n");
        fluid_ = Teuchos::rcp(new FluidFPSI(tmpfluid,actdis,solver,fluidtimeparams,output,isale,dirichletcond));
      }
    }
    break;
    case prb_elch:
    {
      // access the problem-specific parameter list
      const Teuchos::ParameterList& elchcontrol = DRT::Problem::Instance()->ELCHControlParams();
      // is ALE needed or not?
      const INPAR::ELCH::ElchMovingBoundary withale = DRT::INPUT::IntegralValue<INPAR::ELCH::ElchMovingBoundary>(elchcontrol,"MOVINGBOUNDARY");
      if (withale!= INPAR::ELCH::elch_mov_bndry_no)
      {
        Teuchos::RCP<FLD::FluidImplicitTimeInt> tmpfluid;
        if(timeint == INPAR::FLUID::timeint_stationary)
          tmpfluid = Teuchos::rcp(new FLD::TimIntStationary(actdis, solver, fluidtimeparams, output, isale));
        else if(timeint == INPAR::FLUID::timeint_one_step_theta)
          tmpfluid = Teuchos::rcp(new FLD::TimIntOneStepTheta(actdis, solver, fluidtimeparams, output, isale));
        else if(timeint == INPAR::FLUID::timeint_bdf2)
          tmpfluid = Teuchos::rcp(new FLD::TimIntBDF2(actdis, solver, fluidtimeparams, output, isale));
        else if(timeint == INPAR::FLUID::timeint_afgenalpha or
            timeint == INPAR::FLUID::timeint_npgenalpha)
          tmpfluid = Teuchos::rcp(new FLD::TimIntGenAlpha(actdis, solver, fluidtimeparams, output, isale));
        else
          dserror("Unknown time integration for this fluid problem type\n");
        fluid_ = Teuchos::rcp(new FluidFSI(tmpfluid, actdis, solver, fluidtimeparams, output, isale, dirichletcond));
      }
      else
      {
        if(timeint == INPAR::FLUID::timeint_stationary)
          fluid_ = Teuchos::rcp(new FLD::TimIntStationary(actdis, solver, fluidtimeparams, output, isale));
        else if(timeint == INPAR::FLUID::timeint_one_step_theta)
          fluid_ = Teuchos::rcp(new FLD::TimIntOneStepTheta(actdis, solver, fluidtimeparams, output, isale));
        else if(timeint == INPAR::FLUID::timeint_bdf2)
          fluid_ = Teuchos::rcp(new FLD::TimIntBDF2(actdis, solver, fluidtimeparams, output, isale));
        else if(timeint == INPAR::FLUID::timeint_afgenalpha or
            timeint == INPAR::FLUID::timeint_npgenalpha)
          fluid_ = Teuchos::rcp(new FLD::TimIntGenAlpha(actdis, solver, fluidtimeparams, output, isale));
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
    } // end switch (probtype)
  }
  else
  {
    dserror("Unknown time integration for fluid\n");
  }

  // initialize algorithm for specific time-integration scheme
  if (init)
  {
    fluid_->Init();

    // set initial field by given function
    // we do this here, since we have direct access to all necessary parameters
    INPAR::FLUID::InitialField initfield = DRT::INPUT::IntegralValue<INPAR::FLUID::InitialField>(fdyn,"INITIALFIELD");
    if(initfield != INPAR::FLUID::initfield_zero_field)
    {
      int startfuncno = fdyn.get<int>("STARTFUNCNO");
      if (initfield != INPAR::FLUID::initfield_field_by_function and
          initfield != INPAR::FLUID::initfield_disturbed_field_from_function)
      {
        startfuncno=-1;
      }
      fluid_->SetInitialFlowField(initfield,startfuncno);
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidBaseAlgorithm::SetupInflowFluid(
  const Teuchos::ParameterList& prbdyn,
  const Teuchos::RCP<DRT::Discretization> discret)
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("ADAPTER::FluidBaseAlgorithm::SetupFluid");
  Teuchos::TimeMonitor monitor(*t);

  // -------------------------------------------------------------------
  // what's the current problem type?
  // -------------------------------------------------------------------
  PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();

  // the inflow computation can only deal with standard fluid problems so far
  // extensions for xfluid, fsi or combust problems have to be added if necessary
  // they should not pose any additional problem
  // meshtying or xfem related parameters are not supported, yet
  if (probtype != prb_fluid)
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
  const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();

  if (discret->Comm().MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(IO::cout, fdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // get the solver number used for linear fluid solver
  const int linsolvernumber = fdyn.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for fluid problem. Please set LINEAR_SOLVER in FLUID DYNAMIC to a valid number!");
  Teuchos::RCP<LINALG::Solver> solver =
    Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                           discret->Comm(),
                           DRT::Problem::Instance()->ErrorFile()->Handle()));

  discret->ComputeNullSpaceIfNecessary(solver->Params(),true);

  // create a second solver for SIMPLER preconditioner if chosen from input
  CreateSecondSolver(solver,fdyn);

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> fluidtimeparams = Teuchos::rcp(new Teuchos::ParameterList());

  // physical type of fluid flow (incompressible, Boussinesq Approximation, varying density, loma)
  fluidtimeparams->set<int>("Physical Type",
      DRT::INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn,"PHYSICAL_TYPE"));

  // now, set general parameters required for all problems
  SetGeneralParameters(fluidtimeparams,prbdyn,fdyn);

  // overwrite canonical flow parameters by inflow type
  fluidtimeparams->sublist("TURBULENCE MODEL").set<string>("CANONICAL_FLOW",fdyn.sublist("TURBULENT INFLOW").get<std::string>("CANONICAL_INFLOW"));
  fluidtimeparams->sublist("TURBULENCE MODEL").set<string>("HOMDIR",fdyn.sublist("TURBULENT INFLOW").get<std::string>("INFLOW_HOMDIR"));
  fluidtimeparams->sublist("TURBULENCE MODEL").set<int>("DUMPING_PERIOD",fdyn.sublist("TURBULENT INFLOW").get<int>("INFLOW_DUMPING_PERIOD"));
  fluidtimeparams->sublist("TURBULENCE MODEL").set<int>("SAMPLING_START",fdyn.sublist("TURBULENT INFLOW").get<int>("INFLOW_SAMPLING_START"));
  fluidtimeparams->sublist("TURBULENCE MODEL").set<int>("SAMPLING_STOP",fdyn.sublist("TURBULENT INFLOW").get<int>("INFLOW_SAMPLING_STOP"));
  fluidtimeparams->sublist("TURBULENCE MODEL").set<double>("CHAN_AMPL_INIT_DIST",fdyn.sublist("TURBULENT INFLOW").get<double>("INFLOW_INIT_DIST"));

  // -------------------------------------------------------------------
  // additional parameters and algorithm call depending on respective
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  INPAR::FLUID::TimeIntegrationScheme timeint = DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fdyn,"TIMEINTEGR");

  // -------------------------------------------------------------------
  // additional parameters and algorithm call depending on respective
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  if(timeint == INPAR::FLUID::timeint_stationary or
     timeint == INPAR::FLUID::timeint_one_step_theta or
     timeint == INPAR::FLUID::timeint_bdf2 or
     timeint == INPAR::FLUID::timeint_afgenalpha or
     timeint == INPAR::FLUID::timeint_npgenalpha
    )
  {
    // -----------------------------------------------------------------
    // set additional parameters in list for
    // one-step-theta/BDF2/af-generalized-alpha/stationary scheme
    // -----------------------------------------------------------------
    // type of time-integration (or stationary) scheme
    fluidtimeparams->set<int>("time int algo",timeint);
    // parameter theta for time-integration schemes
    fluidtimeparams->set<double>           ("theta"                    ,fdyn.get<double>("THETA"));
    // number of steps for potential start algorithm
    fluidtimeparams->set<int>              ("number of start steps"    ,fdyn.get<int>("NUMSTASTEPS"));
    // parameter theta for potential start algorithm
    fluidtimeparams->set<double>           ("start theta"              ,fdyn.get<double>("START_THETA"));
    // parameter for grid velocity interpolation
    fluidtimeparams->set<int>              ("order gridvel"            ,DRT::INPUT::IntegralValue<int>(fdyn,"GRIDVEL"));

    fluidtimeparams->set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

    //------------------------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor);
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom
//    fluid_ = Teuchos::rcp(new FLD::FluidImplicitTimeInt(discret, solver, fluidtimeparams, output, false));
    if(timeint == INPAR::FLUID::timeint_stationary)
      fluid_ = Teuchos::rcp(new FLD::TimIntStationary(discret, solver, fluidtimeparams, output, false));
    else if(timeint == INPAR::FLUID::timeint_one_step_theta)
      fluid_ = Teuchos::rcp(new FLD::TimIntOneStepTheta(discret, solver, fluidtimeparams, output, false));
    else if(timeint == INPAR::FLUID::timeint_bdf2)
      fluid_ = Teuchos::rcp(new FLD::TimIntBDF2(discret, solver, fluidtimeparams, output, false));
    else if(timeint == INPAR::FLUID::timeint_afgenalpha or
        timeint == INPAR::FLUID::timeint_npgenalpha)
      fluid_ = Teuchos::rcp(new FLD::TimIntGenAlpha(discret, solver, fluidtimeparams, output, false));

  }
  else
  {
    dserror("Unknown time integration for fluid\n");
  }

  // initialize algorithm for specific time-integration scheme
  fluid_->Init();

  // set initial field for inflow section by given function
  // we do this here, since we have direct access to all necessary parameters
  INPAR::FLUID::InitialField initfield = DRT::INPUT::IntegralValue<INPAR::FLUID::InitialField>(fdyn,"INITIALFIELD");
  initfield = DRT::INPUT::IntegralValue<INPAR::FLUID::InitialField>(fdyn.sublist("TURBULENT INFLOW"),"INITIALINFLOWFIELD");
  if(initfield != INPAR::FLUID::initfield_zero_field)
  {
    int startfuncno = fdyn.sublist("TURBULENT INFLOW").get<int>("INFLOWFUNC");
    if (initfield != INPAR::FLUID::initfield_field_by_function and
        initfield != INPAR::FLUID::initfield_disturbed_field_from_function)
    {
      startfuncno=-1;
    }
    fluid_->SetInitialFlowField(initfield,startfuncno);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidBaseAlgorithm::SetGeneralParameters(
 const Teuchos::RCP<Teuchos::ParameterList> fluidtimeparams,
 const Teuchos::ParameterList& prbdyn,
 const Teuchos::ParameterList& fdyn)
{
  fluidtimeparams->set<int>("Simple Preconditioner",DRT::INPUT::IntegralValue<int>(fdyn,"SIMPLER"));
  // fluidtimeparams->set<int>("AMG(BS) Preconditioner",DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(DRT::Problem::Instance()->FluidSolverParams(),"AZPREC")); // probably can be removed

  // -------------------------------------- number of degrees of freedom
  // number of degrees of freedom
  const int ndim = DRT::Problem::Instance()->NDim();
  fluidtimeparams->set<int>("number of velocity degrees of freedom" ,ndim);

  // -------------------------------------------------- time integration
  // note: here, the values are taken out of the problem-dependent ParameterList prbdyn
  // (which also can be fluiddyn itself!)

  // the default time step size
  fluidtimeparams->set<double> ("time step size"      ,prbdyn.get<double>("TIMESTEP"));
  // maximum simulation time
  fluidtimeparams->set<double> ("total time"          ,prbdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  fluidtimeparams->set<int>    ("max number timesteps",prbdyn.get<int>("NUMSTEP"));

  // -------- additional parameters in list for generalized-alpha scheme
#if 1
  // parameter alpha_M
  fluidtimeparams->set<double> ("alpha_M", fdyn.get<double>("ALPHA_M"));
  // parameter alpha_F
  fluidtimeparams->set<double> ("alpha_F", fdyn.get<double>("ALPHA_F"));
  // parameter gamma
  fluidtimeparams->set<double> ("gamma",   fdyn.get<double>("GAMMA"));
#else
  // parameter alpha_M
  fluidtimeparams->set<double> ("alpha_M", 1.-prbdyn.get<double>("ALPHA_M"));
  // parameter alpha_F
  fluidtimeparams->set<double> ("alpha_F", 1.-prbdyn.get<double>("ALPHA_F"));
  // parameter gamma
  fluidtimeparams->set<double> ("gamma",   prbdyn.get<double>("GAMMA"));
#endif

  // ---------------------------------------------- nonlinear iteration
  // type of predictor
  fluidtimeparams->set<string>          ("predictor"                 ,fdyn.get<std::string>("PREDICTOR"));
  // set linearisation scheme
  fluidtimeparams->set<int>("Linearisation", DRT::INPUT::IntegralValue<INPAR::FLUID::LinearisationAction>(fdyn,"NONLINITER"));
  // set bool flag "Newton true or false" for combustion formulation and XFEM
  //fluidtimeparams->set<bool>("Use reaction terms for linearisation",
  //                           DRT::INPUT::IntegralValue<INPAR::FLUID::LinearisationAction>(fdyn,"NONLINITER")== INPAR::FLUID::Newton);
  // maximum number of nonlinear iteration steps
  fluidtimeparams->set<int>             ("max nonlin iter steps"     ,fdyn.get<int>("ITEMAX"));
  // maximum number of nonlinear iteration steps for initial stationary solution
  fluidtimeparams->set<int>             ("max nonlin iter steps init stat sol",fdyn.get<int>("INITSTATITEMAX"));
  // stop nonlinear iteration when both incr-norms are below this bound
  fluidtimeparams->set<double>          ("tolerance for nonlin iter" ,fdyn.get<double>("CONVTOL"));
  // set convergence check
  fluidtimeparams->set<string>          ("CONVCHECK"  ,fdyn.get<std::string>("CONVCHECK"));
  // set recomputation of residual after solution has convergenced
  fluidtimeparams->set<bool>            ("INCONSISTENT_RESIDUAL",DRT::INPUT::IntegralValue<int>(fdyn,"INCONSISTENT_RESIDUAL")==1);
  // set adaptive linear solver tolerance
  fluidtimeparams->set<bool>            ("ADAPTCONV",DRT::INPUT::IntegralValue<int>(fdyn,"ADAPTCONV")==1);
  fluidtimeparams->set<double>          ("ADAPTCONV_BETTER",fdyn.get<double>("ADAPTCONV_BETTER"));
  fluidtimeparams->set<bool>            ("INFNORMSCALING", (DRT::INPUT::IntegralValue<int>(fdyn,"INFNORMSCALING")==1));

  // ----------------------------------------------- restart and output
  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  // restart
  fluidtimeparams->set<int>("write restart every", prbdyn.get<int>("RESTARTEVRY"));
  // solution output
  fluidtimeparams->set<int>("write solution every", prbdyn.get<int>("UPRES"));
  // flag for writing stresses
  fluidtimeparams->set<int>("write stresses"  ,DRT::INPUT::IntegralValue<int>(ioflags,"FLUID_STRESS"));
  // flag for writing wall shear stress
  fluidtimeparams->set<int>("write wall shear stresses"  ,DRT::INPUT::IntegralValue<int>(ioflags,"FLUID_WALL_SHEAR_STRESS"));
  // flag for writing fluid field to gmsh
  fluidtimeparams->set<bool>("GMSH_OUTPUT", DRT::INPUT::IntegralValue<bool>(fdyn,"GMSH_OUTPUT"));
  // flag for computing divergence
  fluidtimeparams->set<bool>("COMPUTE_DIVU", DRT::INPUT::IntegralValue<bool>(fdyn,"COMPUTE_DIVU"));
  // flag for computing lift and drag values
  fluidtimeparams->set<bool>("LIFTDRAG", DRT::INPUT::IntegralValue<bool>(fdyn,"LIFTDRAG"));

  // -------------------------------------------------- Oseen advection
  // set function number of given Oseen advective field
  fluidtimeparams->set<int>("OSEENFIELDFUNCNO", fdyn.get<int>("OSEENFIELDFUNCNO"));

  // ---------------------------------------------------- lift and drag
  fluidtimeparams->set<int>("liftdrag",DRT::INPUT::IntegralValue<int>(fdyn,"LIFTDRAG"));

  // -----------evaluate error for test flows with analytical solutions
  INPAR::FLUID::InitialField initfield = DRT::INPUT::IntegralValue<INPAR::FLUID::InitialField>(fdyn,"INITIALFIELD");
  fluidtimeparams->set<int>("eval err for analyt sol", initfield);

  // ------------------------------------------ form of convective term
  fluidtimeparams->set<string> ("form of convective term", fdyn.get<std::string>("CONVFORM"));

  // -------------------------- potential nonlinear boundary conditions
  fluidtimeparams->set<string> ("Nonlinear boundary conditions",fdyn.get<std::string>("NONLINEARBC"));

  // ------------------------------------ potential reduced_D 3D coupling method
  fluidtimeparams->set<string> ("Strong 3D_redD coupling",fdyn.get<std::string>("STRONG_REDD_3D_COUPLING_TYPE"));

  //--------------------------------------  mesh tying for fluid
  fluidtimeparams->set<int>("MESHTYING",
      DRT::INPUT::IntegralValue<INPAR::FLUID::MeshTying>(fdyn,"MESHTYING"));

  fluidtimeparams->set<bool>("ALLDOFCOUPLED",DRT::INPUT::IntegralValue<bool>(fdyn,"ALLDOFCOUPLED"));

  //--------------------------------------analytical error evaluation
  fluidtimeparams->set<int>("calculate error",
      Teuchos::getIntegralValue<int>(fdyn,"CALCERROR"));

  // -----------------------sublist containing stabilization parameters
  fluidtimeparams->sublist("RESIDUAL-BASED STABILIZATION") =fdyn.sublist("RESIDUAL-BASED STABILIZATION");
  fluidtimeparams->sublist("EDGE-BASED STABILIZATION")     =fdyn.sublist("EDGE-BASED STABILIZATION");
  fluidtimeparams->sublist("POROUS-FLOW STABILIZATION")    =fdyn.sublist("POROUS-FLOW STABILIZATION");

  // -----------------------------get also scatra stabilization sublist
  const Teuchos::ParameterList& scatradyn =
    DRT::Problem::Instance()->ScalarTransportDynamicParams();
  fluidtimeparams->sublist("SCATRA STABILIZATION")=scatradyn.sublist("STABILIZATION");

  // --------------------------sublist containing turbulence parameters
  {
    fluidtimeparams->sublist("TURBULENCE MODEL")=fdyn.sublist("TURBULENCE MODEL");
    fluidtimeparams->sublist("SUBGRID VISCOSITY")=fdyn.sublist("SUBGRID VISCOSITY");
    fluidtimeparams->sublist("MULTIFRACTAL SUBGRID SCALES")=fdyn.sublist("MULTIFRACTAL SUBGRID SCALES");
    fluidtimeparams->sublist("TURBULENT INFLOW")=fdyn.sublist("TURBULENT INFLOW");

    fluidtimeparams->sublist("TURBULENCE MODEL").set<string>("statistics outfile",DRT::Problem::Instance()->OutputControlFile()->FileName());
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidBaseAlgorithm::CreateSecondSolver(
  const Teuchos::RCP<LINALG::Solver> solver,
  const Teuchos::ParameterList& fdyn)
{
  // The SIMPLER (yes,no) parameter only controls whether the fluid matrix is
  // assembled into a 2x2 blocked operator or a plain 1x1 block matrix
  // A "second solver" for the preconditioner is only needed if SIMPLER == yes
  if (DRT::INPUT::IntegralValue<int>(fdyn,"SIMPLER"))
  {
    const int linsolvernumber = fdyn.get<int>("LINEAR_SOLVER");
    INPAR::SOLVER::AzPrecType prec = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(DRT::Problem::Instance()->SolverParams(linsolvernumber),"AZPREC");
    switch (prec) {
    case INPAR::SOLVER::azprec_CheapSIMPLE:
    case INPAR::SOLVER::azprec_TekoSIMPLE:
    {
      // add Inverse1 block for velocity dofs
      // tell Inverse1 block about NodalBlockInformation
      // In contrary to contact/meshtying problems this is necessary here, since we originally have built the
      // null space for the whole problem (velocity and pressure dofs). However, if we split the matrix into
      // velocity and pressure block, we have to adapt the null space information for the subblocks. Therefore
      // we need the nodal block information in the first subblock for the velocities. The pressure null space
      // is trivial to be built using a constant vector
      Teuchos::ParameterList& inv1 = solver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1");
      inv1.sublist("NodalBlockInformation") = solver->Params().sublist("NodalBlockInformation");

      // CheapSIMPLE is somewhat hardwired here
      solver->Params().sublist("CheapSIMPLE Parameters").set("Prec Type","CheapSIMPLE");
      solver->Params().set("FLUID",true);
    }
    break;
    case INPAR::SOLVER::azprec_MueLuAMG_sym:
    {
      // add Inverse1 block for velocity dofs
      // tell Inverse1 block about NodalBlockInformation
      // In contrary to contact/meshtying problems this is necessary here, since we originally have built the
      // null space for the whole problem (velocity and pressure dofs). However, if we split the matrix into
      // velocity and pressure block, we have to adapt the null space information for the subblocks. Therefore
      // we need the nodal block information in the first subblock for the velocities. The pressure null space
      // is trivial to be built using a constant vector
      Teuchos::ParameterList& inv1 = solver->Params().sublist("MueLu Parameters").sublist("SubSmoother1");
      inv1.sublist("NodalBlockInformation") = solver->Params().sublist("NodalBlockInformation");

      solver->Params().sublist("MueLu Parameters").set("FLUID",true);
    }
    break;
    default:
      dserror("If SIMPLER flag is set to YES you can only use CheapSIMPLE or TekoSIMPLE as preconditioners in your fluid solver. Choose CheapSIMPLE or TekoSIMPLE in the SOLVER %i block in your dat file. Alternatively you can also try a multigrid block preconditioner. Use then \"MueLu_sym\" as preconditioner and provide a parameter xml file.",linsolvernumber);
      break;
    }
  }

  return;
}

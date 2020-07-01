/*----------------------------------------------------------------------------*/
/*! \file

 \brief ALE field adapter

 \level 1

 */
/*----------------------------------------------------------------------------*/

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "ad_ale.H"
#include "ad_ale_fluid.H"
#include "ad_ale_fpsi.H"
#include "ad_ale_fsi_msht.H"
#include "ad_ale_fsi.H"
#include "ad_ale_wear.H"
#include "ad_ale_xffsi.H"
#include "../drt_ale/ale.H"

#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"

#include "../drt_inpar/inpar_ale.H"
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_inpar/inpar_fpsi.H"
#include "../drt_lib/drt_periodicbc.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::AleBaseAlgorithm::AleBaseAlgorithm(
    const Teuchos::ParameterList& prbdyn, Teuchos::RCP<DRT::Discretization> actdis)
{
  SetupAle(prbdyn, actdis);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::AleBaseAlgorithm::~AleBaseAlgorithm() {}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::AleBaseAlgorithm::SetupAle(
    const Teuchos::ParameterList& prbdyn, Teuchos::RCP<DRT::Discretization> actdis)
{
  Teuchos::RCP<Teuchos::Time> t =
      Teuchos::TimeMonitor::getNewTimer("ALE::AleBaseAlgorithm::SetupAle");
  Teuchos::TimeMonitor monitor(*t);

  // what's the current problem type?
  const ProblemType probtype = DRT::Problem::Instance()->GetProblemType();

  // ---------------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // ---------------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // ---------------------------------------------------------------------------
  // connect degrees of freedom for coupled nodes
  // ---------------------------------------------------------------------------
  PeriodicBoundaryConditions pbc(actdis);
  pbc.UpdateDofsForPeriodicBoundaryConditions();

  // ---------------------------------------------------------------------------
  // context for output and restart
  // ---------------------------------------------------------------------------
  Teuchos::RCP<IO::DiscretizationWriter> output = actdis->Writer();

  // Output for these problems are not necessary because we write
  // restart data at each time step for visualization
  bool write_output = true;
  if (probtype == prb_uq) write_output = false;


  if (write_output)
  {
    output->WriteMesh(0, 0.0);
  }
  // ---------------------------------------------------------------------------
  // set some pointers and variables
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> adyn =
      Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->AleDynamicParams()));

  // ---------------------------------------------------------------------------
  // create a linear solver
  // ---------------------------------------------------------------------------
  // get the linear solver number
  const int linsolvernumber = adyn->get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    dserror(
        "No linear solver defined for ALE problems. Please set "
        "LINEAR_SOLVER in ALE DYNAMIC to a valid number!");

  Teuchos::RCP<LINALG::Solver> solver =
      Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
          actdis->Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // ---------------------------------------------------------------------------
  // overwrite certain parameters when ALE is part of a multi-field problem
  // ---------------------------------------------------------------------------
  adyn->set<int>("NUMSTEP", prbdyn.get<int>("NUMSTEP"));
  adyn->set<double>("MAXTIME", prbdyn.get<double>("MAXTIME"));
  adyn->set<double>("TIMESTEP", prbdyn.get<double>("TIMESTEP"));
  adyn->set<int>("RESTARTEVRY", prbdyn.get<int>("RESTARTEVRY"));
  adyn->set<int>("RESULTSEVRY", prbdyn.get<int>("RESULTSEVRY"));

  if (probtype == prb_fpsi)
  {
    // FPSI input parameters
    const Teuchos::ParameterList& fpsidyn = DRT::Problem::Instance()->FPSIDynamicParams();
    int coupling = DRT::INPUT::IntegralValue<int>(fpsidyn, "COUPALGO");
    if (coupling == partitioned)
    {
      dserror("partitioned fpsi solution scheme has not been implemented yet.");
    }
  }


  // create the ALE time integrator
  INPAR::ALE::AleDynamic aletype =
      DRT::INPUT::IntegralValue<INPAR::ALE::AleDynamic>(*adyn, "ALE_TYPE");
  Teuchos::RCP<ALE::Ale> ale = Teuchos::null;
  switch (aletype)
  {
    // catch all nonlinear cases
    case INPAR::ALE::solid:
    case INPAR::ALE::laplace_spatial:
    case INPAR::ALE::springs_spatial:
    {
      ale = Teuchos::rcp(new ALE::Ale(actdis, solver, adyn, output));

      break;
    }
    // catch and linear cases
    case INPAR::ALE::solid_linear:
    case INPAR::ALE::laplace_material:
    case INPAR::ALE::springs_material:
    {
      ale = Teuchos::rcp(new ALE::AleLinear(actdis, solver, adyn, output));

      break;
    }
    default:
    {
      dserror("Decide, whether ALE_TYPE = '%s' is linear or nonlinear.",
          adyn->get<std::string>("ALE_TYPE").c_str());
      break;
    }
  }

  /* determine problem type and then wrap the ALE time integrator into a
   * problem-specific wrapper */
  switch (probtype)
  {
    case prb_ale:
    {
      ale_ = ale;
      break;
    }
    case prb_fsi:
    case prb_gas_fsi:
    case prb_thermo_fsi:
    case prb_ac_fsi:
    case prb_biofilm_fsi:
    {
      const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
      int coupling = DRT::INPUT::IntegralValue<int>(fsidyn, "COUPALGO");
      if (coupling == fsi_iter_monolithicfluidsplit or
          coupling == fsi_iter_monolithicstructuresplit or
          coupling == fsi_iter_constr_monolithicfluidsplit or
          coupling == fsi_iter_constr_monolithicstructuresplit or
          coupling == fsi_iter_lung_monolithicfluidsplit or
          coupling == fsi_iter_lung_monolithicstructuresplit or
          coupling == fsi_iter_mortar_monolithicstructuresplit or
          coupling == fsi_iter_mortar_monolithicfluidsplit)
      {
        ale_ = Teuchos::rcp(new ADAPTER::AleFsiWrapper(ale));
      }
      else if (coupling == fsi_iter_sliding_monolithicfluidsplit or
               coupling == fsi_iter_sliding_monolithicstructuresplit)
      {
        ale_ = Teuchos::rcp(new ADAPTER::AleFsiMshtWrapper(ale));
      }
      else if (coupling == fsi_iter_fluidfluid_monolithicstructuresplit or
               coupling == fsi_iter_fluidfluid_monolithicfluidsplit or
               coupling == fsi_iter_fluidfluid_monolithicstructuresplit_nonox or
               coupling == fsi_iter_fluidfluid_monolithicfluidsplit_nonox)
      {
        ale_ = Teuchos::rcp(new ADAPTER::AleXFFsiWrapper(ale));
      }
      else if (coupling == fsi_iter_stagg_AITKEN_rel_force or
               coupling == fsi_iter_stagg_AITKEN_rel_param or
               coupling == fsi_iter_stagg_CHEB_rel_param or coupling == fsi_iter_stagg_MFNK_FD or
               coupling == fsi_iter_stagg_MFNK_FSI or coupling == fsi_iter_stagg_MPE or
               coupling == fsi_iter_stagg_NLCG or coupling == fsi_iter_stagg_Newton_FD or
               coupling == fsi_iter_stagg_Newton_I or coupling == fsi_iter_stagg_RRE or
               coupling == fsi_iter_stagg_fixed_rel_param or
               coupling == fsi_iter_stagg_steep_desc or coupling == fsi_iter_stagg_steep_desc_force)
      {
        ale_ = Teuchos::rcp(new ADAPTER::AleFluidWrapper(ale));
      }
      else
      {
        dserror(
            "No ALE adapter available yet for your chosen FSI coupling "
            "algorithm!");
      }
      break;
    }
    case prb_fsi_lung:
    {
      const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
      int coupling = DRT::INPUT::IntegralValue<int>(fsidyn, "COUPALGO");
      if (coupling == fsi_iter_lung_monolithicfluidsplit or
          coupling == fsi_iter_lung_monolithicstructuresplit)
      {
        ale_ = Teuchos::rcp(new ADAPTER::AleFsiWrapper(ale));
      }
      else
      {
        dserror(
            "No ALE adapter available yet for your chosen FSI coupling "
            "algorithm!");
      }
      break;
    }
    case prb_fsi_redmodels:
    {
      const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
      int coupling = DRT::INPUT::IntegralValue<int>(fsidyn, "COUPALGO");
      if (coupling == fsi_iter_monolithicfluidsplit or
          coupling == fsi_iter_monolithicstructuresplit or
          coupling == fsi_iter_constr_monolithicfluidsplit or
          coupling == fsi_iter_constr_monolithicstructuresplit or
          coupling == fsi_iter_mortar_monolithicstructuresplit or
          coupling == fsi_iter_mortar_monolithicfluidsplit or
          coupling == fsi_iter_sliding_monolithicfluidsplit or
          coupling == fsi_iter_sliding_monolithicstructuresplit)
      {
        ale_ = Teuchos::rcp(new ADAPTER::AleFsiWrapper(ale));
      }
      else if (coupling == fsi_iter_stagg_AITKEN_rel_force or
               coupling == fsi_iter_stagg_AITKEN_rel_param or
               coupling == fsi_iter_stagg_CHEB_rel_param or coupling == fsi_iter_stagg_MFNK_FD or
               coupling == fsi_iter_stagg_MFNK_FSI or coupling == fsi_iter_stagg_MPE or
               coupling == fsi_iter_stagg_NLCG or coupling == fsi_iter_stagg_Newton_FD or
               coupling == fsi_iter_stagg_Newton_I or coupling == fsi_iter_stagg_RRE or
               coupling == fsi_iter_stagg_fixed_rel_param or
               coupling == fsi_iter_stagg_steep_desc or coupling == fsi_iter_stagg_steep_desc_force)
      {
        ale_ = Teuchos::rcp(new ADAPTER::AleFluidWrapper(ale));
      }
      else
      {
        dserror(
            "No ALE adapter available yet for your chosen FSI coupling "
            "algorithm!");
      }
      break;
    }
    case prb_fpsi:
    case prb_fps3i:
    case prb_fsi_xfem:
    case prb_fpsi_xfem:
    {
      ale_ = Teuchos::rcp(new ADAPTER::AleFpsiWrapper(ale));
      break;
    }
    case prb_struct_ale:
    {
      ale_ = Teuchos::rcp(new ADAPTER::AleWearWrapper(ale));
      break;
    }
    case prb_uq:
    {
      ale_ = ale;
      break;
    }
    case prb_freesurf:
    case prb_fluid_ale:
    case prb_elch:
    case prb_fluid_xfem:
    {
      ale_ = Teuchos::rcp(new ADAPTER::AleFluidWrapper(ale));
      break;
    }
    default:
      dserror("ALE type not implemented yet!!");
      break;
  }

  return;
}

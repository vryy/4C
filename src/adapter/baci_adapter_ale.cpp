/*----------------------------------------------------------------------------*/
/*! \file

 \brief ALE field adapter

 \level 1

 */
/*----------------------------------------------------------------------------*/

#include "baci_adapter_ale.hpp"

#include "baci_adapter_ale_fluid.hpp"
#include "baci_adapter_ale_fpsi.hpp"
#include "baci_adapter_ale_fsi.hpp"
#include "baci_adapter_ale_fsi_msht.hpp"
#include "baci_adapter_ale_wear.hpp"
#include "baci_adapter_ale_xffsi.hpp"
#include "baci_ale.hpp"
#include "baci_global_data.hpp"
#include "baci_inpar_ale.hpp"
#include "baci_inpar_fpsi.hpp"
#include "baci_inpar_fsi.hpp"
#include "baci_inpar_validparameters.hpp"
#include "baci_io.hpp"
#include "baci_io_control.hpp"
#include "baci_io_pstream.hpp"
#include "baci_lib_discret.hpp"
#include "baci_lib_periodicbc.hpp"
#include "baci_linalg_sparsematrix.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"
#include "baci_linear_solver_method_linalg.hpp"
#include "baci_utils_parameter_list.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::AleBaseAlgorithm::AleBaseAlgorithm(
    const Teuchos::ParameterList& prbdyn, Teuchos::RCP<DRT::Discretization> actdis)
{
  SetupAle(prbdyn, actdis);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::AleBaseAlgorithm::SetupAle(
    const Teuchos::ParameterList& prbdyn, Teuchos::RCP<DRT::Discretization> actdis)
{
  Teuchos::RCP<Teuchos::Time> t =
      Teuchos::TimeMonitor::getNewTimer("ALE::AleBaseAlgorithm::SetupAle");
  Teuchos::TimeMonitor monitor(*t);

  // what's the current problem type?
  const GLOBAL::ProblemType probtype = GLOBAL::Problem::Instance()->GetProblemType();

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
  output->WriteMesh(0, 0.0);

  // ---------------------------------------------------------------------------
  // set some pointers and variables
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> adyn =
      Teuchos::rcp(new Teuchos::ParameterList(GLOBAL::Problem::Instance()->AleDynamicParams()));

  // ---------------------------------------------------------------------------
  // create a linear solver
  // ---------------------------------------------------------------------------
  // get the linear solver number
  const int linsolvernumber = adyn->get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    dserror(
        "No linear solver defined for ALE problems. Please set "
        "LINEAR_SOLVER in ALE DYNAMIC to a valid number!");

  Teuchos::RCP<CORE::LINALG::Solver> solver = Teuchos::rcp(new CORE::LINALG::Solver(
      GLOBAL::Problem::Instance()->SolverParams(linsolvernumber), actdis->Comm()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // ---------------------------------------------------------------------------
  // overwrite certain parameters when ALE is part of a multi-field problem
  // ---------------------------------------------------------------------------
  adyn->set<int>("NUMSTEP", prbdyn.get<int>("NUMSTEP"));
  adyn->set<double>("MAXTIME", prbdyn.get<double>("MAXTIME"));
  adyn->set<double>("TIMESTEP", prbdyn.get<double>("TIMESTEP"));
  adyn->set<int>("RESTARTEVRY", prbdyn.get<int>("RESTARTEVRY"));
  adyn->set<int>("RESULTSEVRY", prbdyn.get<int>("RESULTSEVRY"));

  if (probtype == GLOBAL::ProblemType::fpsi)
  {
    // FPSI input parameters
    const Teuchos::ParameterList& fpsidyn = GLOBAL::Problem::Instance()->FPSIDynamicParams();
    int coupling = CORE::UTILS::IntegralValue<int>(fpsidyn, "COUPALGO");
    if (coupling == partitioned)
    {
      dserror("partitioned fpsi solution scheme has not been implemented yet.");
    }
  }


  // create the ALE time integrator
  INPAR::ALE::AleDynamic aletype =
      CORE::UTILS::IntegralValue<INPAR::ALE::AleDynamic>(*adyn, "ALE_TYPE");
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
    case GLOBAL::ProblemType::ale:
    {
      ale_ = ale;
      break;
    }
    case GLOBAL::ProblemType::fsi:
    case GLOBAL::ProblemType::gas_fsi:
    case GLOBAL::ProblemType::thermo_fsi:
    case GLOBAL::ProblemType::ac_fsi:
    case GLOBAL::ProblemType::biofilm_fsi:
    {
      const Teuchos::ParameterList& fsidyn = GLOBAL::Problem::Instance()->FSIDynamicParams();
      int coupling = CORE::UTILS::IntegralValue<int>(fsidyn, "COUPALGO");
      if (coupling == fsi_iter_monolithicfluidsplit or
          coupling == fsi_iter_monolithicstructuresplit or
          coupling == fsi_iter_constr_monolithicfluidsplit or
          coupling == fsi_iter_constr_monolithicstructuresplit or
          coupling == fsi_iter_lung_monolithicfluidsplit or
          coupling == fsi_iter_lung_monolithicstructuresplit or
          coupling == fsi_iter_mortar_monolithicstructuresplit or
          coupling == fsi_iter_mortar_monolithicfluidsplit or
          coupling == fsi_iter_mortar_monolithicfluidsplit_saddlepoint)
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
    case GLOBAL::ProblemType::fsi_lung:
    {
      const Teuchos::ParameterList& fsidyn = GLOBAL::Problem::Instance()->FSIDynamicParams();
      int coupling = CORE::UTILS::IntegralValue<int>(fsidyn, "COUPALGO");
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
    case GLOBAL::ProblemType::fsi_redmodels:
    {
      const Teuchos::ParameterList& fsidyn = GLOBAL::Problem::Instance()->FSIDynamicParams();
      int coupling = CORE::UTILS::IntegralValue<int>(fsidyn, "COUPALGO");
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
    case GLOBAL::ProblemType::fpsi:
    case GLOBAL::ProblemType::fps3i:
    case GLOBAL::ProblemType::fsi_xfem:
    case GLOBAL::ProblemType::fpsi_xfem:
    {
      ale_ = Teuchos::rcp(new ADAPTER::AleFpsiWrapper(ale));
      break;
    }
    case GLOBAL::ProblemType::struct_ale:
    {
      ale_ = Teuchos::rcp(new ADAPTER::AleWearWrapper(ale));
      break;
    }
    case GLOBAL::ProblemType::freesurf:
    case GLOBAL::ProblemType::fluid_ale:
    case GLOBAL::ProblemType::elch:
    case GLOBAL::ProblemType::fluid_xfem:
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

FOUR_C_NAMESPACE_CLOSE

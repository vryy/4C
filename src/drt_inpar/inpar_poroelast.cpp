/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for poro elasticity

\level 2


*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_poroelast.H"
#include "inpar_fluid.H"

#include "../linalg/linalg_equilibrate.H"


void INPAR::POROELAST::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::Array<std::string> yesnotuple =
      tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

  Teuchos::ParameterList& poroelastdyn =
      list->sublist("POROELASTICITY DYNAMIC", false, "Poroelasticity");

  // Coupling strategy for (monolithic) porous media solvers
  setStringToIntegralParameter<int>("COUPALGO", "poro_monolithic",
      "Coupling strategies for poroelasticity solvers",
      tuple<std::string>("poro_partitioned", "poro_monolithic", "poro_monolithicstructuresplit",
          "poro_monolithicfluidsplit", "poro_monolithicnopenetrationsplit",
          "poro_monolithicmeshtying"),
      tuple<int>(Partitioned, Monolithic, Monolithic_structuresplit, Monolithic_fluidsplit,
          Monolithic_nopenetrationsplit, Monolithic_meshtying),
      &poroelastdyn);

  // physical type of poro fluid flow (incompressible, varying density, loma, Boussinesq
  // approximation)
  setStringToIntegralParameter<int>("PHYSICAL_TYPE", "Poro", "Physical Type of Porofluid",
      tuple<std::string>("Poro", "Poro_P1"), tuple<int>(INPAR::FLUID::poro, INPAR::FLUID::poro_p1),
      &poroelastdyn);

  // physical type of poro fluid flow (incompressible, varying density, loma, Boussinesq
  // approximation)
  setStringToIntegralParameter<int>("TRANSIENT_TERMS", "all",
      "which equation includes transient terms",
      tuple<std::string>("none", "momentum", "continuity", "all"),
      tuple<int>(transient_none, transient_momentum_only, transient_continuity_only, transient_all),
      &poroelastdyn);

  // Output type
  IntParameter(
      "RESTARTEVRY", 1, "write restart possibility every RESTARTEVRY steps", &poroelastdyn);
  // Time loop control
  IntParameter("NUMSTEP", 200, "maximum number of Timesteps", &poroelastdyn);
  DoubleParameter("MAXTIME", 1000.0, "total simulation time", &poroelastdyn);
  DoubleParameter("TIMESTEP", 0.05, "time step size dt", &poroelastdyn);
  IntParameter("ITEMAX", 10, "maximum number of iterations over fields", &poroelastdyn);
  IntParameter("ITEMIN", 1, "minimal number of iterations over fields", &poroelastdyn);
  IntParameter("RESULTSEVRY", 1, "increment for writing solution", &poroelastdyn);

  // Iterationparameters
  DoubleParameter("TOLRES_GLOBAL", 1e-8, "tolerance in the residual norm for the Newton iteration",
      &poroelastdyn);
  DoubleParameter("TOLINC_GLOBAL", 1e-8, "tolerance in the increment norm for the Newton iteration",
      &poroelastdyn);
  DoubleParameter("TOLRES_DISP", 1e-8, "tolerance in the residual norm for the Newton iteration",
      &poroelastdyn);
  DoubleParameter("TOLINC_DISP", 1e-8, "tolerance in the increment norm for the Newton iteration",
      &poroelastdyn);
  DoubleParameter("TOLRES_PORO", 1e-8, "tolerance in the residual norm for the Newton iteration",
      &poroelastdyn);
  DoubleParameter("TOLINC_PORO", 1e-8, "tolerance in the increment norm for the Newton iteration",
      &poroelastdyn);
  DoubleParameter(
      "TOLRES_VEL", 1e-8, "tolerance in the residual norm for the Newton iteration", &poroelastdyn);
  DoubleParameter("TOLINC_VEL", 1e-8, "tolerance in the increment norm for the Newton iteration",
      &poroelastdyn);
  DoubleParameter("TOLRES_PRES", 1e-8, "tolerance in the residual norm for the Newton iteration",
      &poroelastdyn);
  DoubleParameter("TOLINC_PRES", 1e-8, "tolerance in the increment norm for the Newton iteration",
      &poroelastdyn);
  DoubleParameter("TOLRES_NCOUP", 1e-8, "tolerance in the residual norm for the Newton iteration",
      &poroelastdyn);

  setStringToIntegralParameter<int>("NORM_INC", "AbsSingleFields",
      "type of norm for primary variables convergence check",
      tuple<std::string>("AbsGlobal", "AbsSingleFields"),
      tuple<int>(convnorm_abs_global, convnorm_abs_singlefields), &poroelastdyn);

  setStringToIntegralParameter<int>("NORM_RESF", "AbsSingleFields",
      "type of norm for residual convergence check",
      tuple<std::string>("AbsGlobal", "AbsSingleFields"),
      tuple<int>(convnorm_abs_global, convnorm_abs_singlefields), &poroelastdyn);

  setStringToIntegralParameter<int>("NORMCOMBI_RESFINC", "And",
      "binary operator to combine primary variables and residual force values",
      tuple<std::string>("And", "Or"), tuple<int>(bop_and, bop_or), &poroelastdyn);

  setStringToIntegralParameter<int>("VECTORNORM_RESF", "L2",
      "type of norm to be applied to residuals",
      tuple<std::string>("L1", "L1_Scaled", "L2", "Rms", "Inf"),
      tuple<int>(norm_l1, norm_l1_scaled, norm_l2, norm_rms, norm_inf), &poroelastdyn);

  setStringToIntegralParameter<int>("VECTORNORM_INC", "L2",
      "type of norm to be applied to residuals",
      tuple<std::string>("L1", "L1_Scaled", "L2", "Rms", "Inf"),
      tuple<int>(norm_l1, norm_l1_scaled, norm_l2, norm_rms, norm_inf), &poroelastdyn);

  setStringToIntegralParameter<int>("SECONDORDER", "Yes", "Second order coupling at the interface.",
      yesnotuple, yesnovalue, &poroelastdyn);

  setStringToIntegralParameter<int>("CONTIPARTINT", "No",
      "Partial integration of porosity gradient in continuity equation", yesnotuple, yesnovalue,
      &poroelastdyn);

  setStringToIntegralParameter<int>("CONTACTNOPEN", "No",
      "No-Penetration Condition on active contact surface in case of poro contact problem!",
      yesnotuple, yesnovalue, &poroelastdyn);

  BoolParameter("MATCHINGGRID", "Yes", "is matching grid", &poroelastdyn);

  // number of linear solver used for poroelasticity
  IntParameter("LINEAR_SOLVER", -1, "number of linear solver used for poroelasticity problems",
      &poroelastdyn);

  // flag for equilibration of global system of equations
  setStringToIntegralParameter<LINALG::EquilibrationMethod>("EQUILIBRATION", "none",
      "flag for equilibration of global system of equations",
      tuple<std::string>("none", "rows_full", "rows_maindiag", "columns_full", "columns_maindiag",
          "rowsandcolumns_full", "rowsandcolumns_maindiag"),
      tuple<LINALG::EquilibrationMethod>(LINALG::EquilibrationMethod::none,
          LINALG::EquilibrationMethod::rows_full, LINALG::EquilibrationMethod::rows_maindiag,
          LINALG::EquilibrationMethod::columns_full, LINALG::EquilibrationMethod::columns_maindiag,
          LINALG::EquilibrationMethod::rowsandcolumns_full,
          LINALG::EquilibrationMethod::rowsandcolumns_maindiag),
      &poroelastdyn);
}

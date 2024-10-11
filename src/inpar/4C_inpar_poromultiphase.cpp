/*----------------------------------------------------------------------*/
/*! \file
 \brief input parameters for porous multiphase problem

   \level 3

 *----------------------------------------------------------------------*/



#include "4C_inpar_poromultiphase.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_linalg_equilibrate.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


void Inpar::POROMULTIPHASE::set_valid_parameters(Teuchos::ParameterList& list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  // ----------------------------------------------------------------------
  // (1) general control parameters
  Teuchos::ParameterList& poromultiphasedyn = list.sublist(
      "POROMULTIPHASE DYNAMIC", false, "Control paramters for multiphase porous medium");

  // Output type
  Core::UTILS::int_parameter(
      "RESTARTEVRY", 1, "write restart possibility every RESTARTEVRY steps", &poromultiphasedyn);
  // Time loop control
  Core::UTILS::int_parameter("NUMSTEP", 200, "maximum number of Timesteps", &poromultiphasedyn);
  Core::UTILS::double_parameter("MAXTIME", 1000.0, "total simulation time", &poromultiphasedyn);
  Core::UTILS::double_parameter("TIMESTEP", -1, "time step size dt", &poromultiphasedyn);
  Core::UTILS::int_parameter(
      "RESULTSEVRY", 1, "increment for writing solution", &poromultiphasedyn);
  Core::UTILS::int_parameter(
      "ITEMAX", 10, "maximum number of iterations over fields", &poromultiphasedyn);

  // here the computation of the structure can be skipped, this is helpful if only fluid-scatra
  // coupling should be calculated
  Core::UTILS::bool_parameter(
      "SOLVE_STRUCTURE", "yes", "Flag to skip computation of structural field", &poromultiphasedyn);


  // Coupling strategy for solvers
  setStringToIntegralParameter<SolutionSchemeOverFields>("COUPALGO", "twoway_partitioned",
      "Coupling strategies for poro multiphase solvers",
      tuple<std::string>("twoway_partitioned", "twoway_monolithic"),
      tuple<SolutionSchemeOverFields>(solscheme_twoway_partitioned, solscheme_twoway_monolithic),
      &poromultiphasedyn);

  // coupling with 1D artery network active
  Core::UTILS::bool_parameter(
      "ARTERY_COUPLING", "No", "Coupling with 1D blood vessels.", &poromultiphasedyn);

  // ----------------------------------------------------------------------
  // (2) monolithic parameters
  Teuchos::ParameterList& poromultiphasedynmono = poromultiphasedyn.sublist(
      "MONOLITHIC", false, "Parameters for monolithic Poro-Multiphase-Scatra Interaction");

  // convergence tolerances for monolithic coupling
  Core::UTILS::double_parameter("TOLRES_GLOBAL", 1e-8,
      "tolerance in the residual norm for the Newton iteration", &poromultiphasedynmono);
  Core::UTILS::double_parameter("TOLINC_GLOBAL", 1e-8,
      "tolerance in the increment norm for the Newton iteration", &poromultiphasedynmono);

  // number of linear solver used for poroelasticity
  Core::UTILS::int_parameter("LINEAR_SOLVER", -1,
      "number of linear solver used for poroelasticity problems", &poromultiphasedynmono);

  // parameters for finite difference check
  setStringToIntegralParameter<FdCheck>("FDCHECK", "none",
      "flag for finite difference check: none or global",
      tuple<std::string>("none",
          "global"),  // perform finite difference check on time integrator level
      tuple<FdCheck>(fdcheck_none, fdcheck_global), &poromultiphasedynmono);

  setStringToIntegralParameter<VectorNorm>("VECTORNORM_RESF", "L2",
      "type of norm to be applied to residuals",
      tuple<std::string>("L1", "L1_Scaled", "L2", "Rms", "Inf"),
      tuple<VectorNorm>(Inpar::POROMULTIPHASE::norm_l1, Inpar::POROMULTIPHASE::norm_l1_scaled,
          Inpar::POROMULTIPHASE::norm_l2, Inpar::POROMULTIPHASE::norm_rms,
          Inpar::POROMULTIPHASE::norm_inf),
      &poromultiphasedynmono);

  setStringToIntegralParameter<VectorNorm>("VECTORNORM_INC", "L2",
      "type of norm to be applied to residuals",
      tuple<std::string>("L1", "L1_Scaled", "L2", "Rms", "Inf"),
      tuple<VectorNorm>(Inpar::POROMULTIPHASE::norm_l1, Inpar::POROMULTIPHASE::norm_l1_scaled,
          Inpar::POROMULTIPHASE::norm_l2, Inpar::POROMULTIPHASE::norm_rms,
          Inpar::POROMULTIPHASE::norm_inf),
      &poromultiphasedynmono);

  // flag for equilibration of global system of equations
  setStringToIntegralParameter<Core::LinAlg::EquilibrationMethod>("EQUILIBRATION", "none",
      "flag for equilibration of global system of equations",
      tuple<std::string>("none", "rows_full", "rows_maindiag", "columns_full", "columns_maindiag",
          "rowsandcolumns_full", "rowsandcolumns_maindiag"),
      tuple<Core::LinAlg::EquilibrationMethod>(Core::LinAlg::EquilibrationMethod::none,
          Core::LinAlg::EquilibrationMethod::rows_full,
          Core::LinAlg::EquilibrationMethod::rows_maindiag,
          Core::LinAlg::EquilibrationMethod::columns_full,
          Core::LinAlg::EquilibrationMethod::columns_maindiag,
          Core::LinAlg::EquilibrationMethod::rowsandcolumns_full,
          Core::LinAlg::EquilibrationMethod::rowsandcolumns_maindiag),
      &poromultiphasedynmono);

  // convergence criteria adaptivity --> note ADAPTCONV_BETTER set pretty small
  Core::UTILS::bool_parameter("ADAPTCONV", "No",
      "Switch on adaptive control of linear solver tolerance for nonlinear solution",
      &poromultiphasedynmono);
  Core::UTILS::double_parameter("ADAPTCONV_BETTER", 0.001,
      "The linear solver shall be this much better "
      "than the current nonlinear residual in the nonlinear convergence limit",
      &poromultiphasedynmono);

  // ----------------------------------------------------------------------
  // (3) partitioned parameters
  Teuchos::ParameterList& poromultiphasedynpart = poromultiphasedyn.sublist(
      "PARTITIONED", false, "Parameters for partitioned Poro-Multiphase-Scatra Interaction");

  // convergence tolerance of outer iteration loop
  Core::UTILS::double_parameter("CONVTOL", 1e-6,
      "tolerance for convergence check of outer iteration", &poromultiphasedynpart);

  // flag for relaxation of partitioned scheme
  setStringToIntegralParameter<RelaxationMethods>("RELAXATION", "none",
      "flag for relaxation of partitioned scheme", tuple<std::string>("none", "Constant", "Aitken"),
      tuple<RelaxationMethods>(relaxation_none, relaxation_constant, relaxation_aitken),
      &poromultiphasedynpart);

  // parameters for relaxation of partitioned coupling
  Core::UTILS::double_parameter(
      "STARTOMEGA", 1.0, "fixed relaxation parameter", &poromultiphasedynpart);
  Core::UTILS::double_parameter(
      "MINOMEGA", 0.1, "smallest omega allowed for Aitken relaxation", &poromultiphasedynpart);
  Core::UTILS::double_parameter(
      "MAXOMEGA", 10.0, "largest omega allowed for Aitken relaxation", &poromultiphasedynpart);
}

FOUR_C_NAMESPACE_CLOSE

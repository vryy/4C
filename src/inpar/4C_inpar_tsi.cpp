/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameters for tsi

\level 1

*/

/*----------------------------------------------------------------------*/



#include "4C_inpar_tsi.hpp"

#include "4C_inpar_contact.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::TSI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& tsidyn = list->sublist("TSI DYNAMIC", false,
      "Thermo Structure Interaction\n"
      "Dynamic section for TSI solver with various coupling methods");

  // coupling strategy for (partitioned and monolithic) TSI solvers
  setStringToIntegralParameter<int>("COUPALGO", "tsi_monolithic",
      "Coupling strategies for TSI solvers",
      tuple<std::string>("tsi_oneway", "tsi_sequstagg", "tsi_iterstagg", "tsi_iterstagg_aitken",
          "tsi_iterstagg_aitkenirons", "tsi_iterstagg_fixedrelax", "tsi_monolithic"),
      tuple<int>(OneWay, SequStagg, IterStagg, IterStaggAitken, IterStaggAitkenIrons,
          IterStaggFixedRel, Monolithic),
      &tsidyn);

  Core::UTILS::BoolParameter("MATCHINGGRID", "Yes", "is matching grid", &tsidyn);

  // output type
  Core::UTILS::IntParameter(
      "RESTARTEVRY", 1, "write restart possibility every RESTARTEVRY steps", &tsidyn);

  // time loop control
  Core::UTILS::IntParameter("NUMSTEP", 200, "maximum number of Timesteps", &tsidyn);
  Core::UTILS::DoubleParameter("MAXTIME", 1000.0, "total simulation time", &tsidyn);
  Core::UTILS::DoubleParameter("TIMESTEP", 0.05, "time step size dt", &tsidyn);
  Core::UTILS::IntParameter("ITEMAX", 10, "maximum number of iterations over fields", &tsidyn);
  Core::UTILS::IntParameter("ITEMIN", 1, "minimal number of iterations over fields", &tsidyn);
  Core::UTILS::IntParameter("RESULTSEVRY", 1, "increment for writing solution", &tsidyn);

  setStringToIntegralParameter<int>("NORM_INC", "Abs",
      "type of norm for convergence check of primary variables in TSI",
      tuple<std::string>("Abs", "Rel", "Mix"), tuple<int>(convnorm_abs, convnorm_rel, convnorm_mix),
      &tsidyn);

  /*----------------------------------------------------------------------*/
  /* parameters for monolithic TSI */
  Teuchos::ParameterList& tsidynmono = tsidyn.sublist("MONOLITHIC", false,
      "Monolithic Thermo Structure Interaction\n"
      "Dynamic section for monolithic TSI");

  // convergence tolerance of tsi residual
  Core::UTILS::DoubleParameter(
      "CONVTOL", 1e-6, "tolerance for convergence check of TSI", &tsidynmono);
  // Iterationparameters
  Core::UTILS::DoubleParameter("TOLINC", 1.0e-6,
      "tolerance for convergence check of TSI-increment in monolithic TSI", &tsidynmono);

  setStringToIntegralParameter<int>("NORM_RESF", "Abs",
      "type of norm for residual convergence check", tuple<std::string>("Abs", "Rel", "Mix"),
      tuple<int>(convnorm_abs, convnorm_rel, convnorm_mix), &tsidynmono);

  setStringToIntegralParameter<int>("NORMCOMBI_RESFINC", "Coupl_And_Singl",
      "binary operator to combine primary variables and residual force values",
      tuple<std::string>("And", "Or", "Coupl_Or_Singl", "Coupl_And_Singl", "And_Singl", "Or_Singl"),
      tuple<int>(
          bop_and, bop_or, bop_coupl_or_singl, bop_coupl_and_singl, bop_and_singl, bop_or_singl),
      &tsidynmono);

  setStringToIntegralParameter<int>("ITERNORM", "Rms", "type of norm to be applied to residuals",
      tuple<std::string>("L1", "L1_Scaled", "L2", "Rms", "Inf"),
      tuple<int>(norm_l1, norm_l1_scaled, norm_l2, norm_rms, norm_inf), &tsidynmono);

  setStringToIntegralParameter<int>("NLNSOL", "fullnewton", "Nonlinear solution technique",
      tuple<std::string>("fullnewton", "ptc"), tuple<int>(soltech_newtonfull, soltech_ptc),
      &tsidynmono);

  Core::UTILS::DoubleParameter("PTCDT", 0.1,
      "pseudo time step for pseudo-transient continuation (PTC) stabilised Newton procedure",
      &tsidynmono);

  // number of linear solver used for monolithic TSI
  Core::UTILS::IntParameter(
      "LINEAR_SOLVER", -1, "number of linear solver used for monolithic TSI problems", &tsidynmono);

  // convergence criteria adaptivity of monolithic TSI solver
  Core::UTILS::BoolParameter("ADAPTCONV", "No",
      "Switch on adaptive control of linear solver tolerance for nonlinear solution", &tsidynmono);
  Core::UTILS::DoubleParameter("ADAPTCONV_BETTER", 0.1,
      "The linear solver shall be this much better than the current nonlinear residual in the "
      "nonlinear convergence limit",
      &tsidynmono);

  Core::UTILS::BoolParameter(
      "INFNORMSCALING", "yes", "Scale blocks of matrix with row infnorm?", &tsidynmono);

  // merge TSI block matrix to enable use of direct solver in monolithic TSI
  // default: "No", i.e. use block matrix
  Core::UTILS::BoolParameter("MERGE_TSI_BLOCK_MATRIX", "No", "Merge TSI block matrix", &tsidynmono);

  // in case of monolithic TSI nodal values (displacements, temperatures and
  // reaction forces) at fix points of the body can be calculated
  // default: "No", i.e. nothing is calculated
  Core::UTILS::BoolParameter("CALC_NECKING_TSI_VALUES", "No",
      "Calculate nodal values for evaluation and validation of necking", &tsidynmono);

  /*----------------------------------------------------------------------*/
  /* parameters for partitioned TSI */
  Teuchos::ParameterList& tsidynpart = tsidyn.sublist("PARTITIONED", false,
      "Partitioned Thermo Structure Interaction\n"
      "Dynamic section for partitioned TSI");

  // decide in partitioned TSI which one-way coupling or predictor should be used
  setStringToIntegralParameter<int>("COUPVARIABLE", "Displacement", "Coupling variable",
      tuple<std::string>("Displacement", "Temperature"), tuple<int>(0, 1), &tsidynpart);

  // Solver parameter for relaxation of iterative staggered partitioned TSI
  Core::UTILS::DoubleParameter("MAXOMEGA", 0.0,
      "largest omega allowed for Aitken relaxation (0.0 means no constraint)", &tsidynpart);
  Core::UTILS::DoubleParameter("FIXEDOMEGA", 1.0, "fixed relaxation parameter", &tsidynpart);

  // convergence tolerance of outer iteration loop
  Core::UTILS::DoubleParameter("CONVTOL", 1e-6,
      "tolerance for convergence check of outer iteraiton within partitioned TSI", &tsidynpart);

  /*----------------------------------------------------------------------*/
  /* parameters for tsi contact */
  Teuchos::ParameterList& tsic = list->sublist("TSI CONTACT", false, "");

  Core::UTILS::DoubleParameter(
      "HEATTRANSSLAVE", 0.0, "Heat transfer parameter for slave side in thermal contact", &tsic);
  Core::UTILS::DoubleParameter(
      "HEATTRANSMASTER", 0.0, "Heat transfer parameter for master side in thermal contact", &tsic);
  Core::UTILS::DoubleParameter("TEMP_DAMAGE", 1.0e12,
      "damage temperatue at contact interface: friction coefficient zero there", &tsic);
  Core::UTILS::DoubleParameter("TEMP_REF", 0.0,
      "reference temperatue at contact interface: friction coefficient equals the given value",
      &tsic);

  Core::UTILS::DoubleParameter(
      "NITSCHE_THETA_TSI", 0.0, "+1: symmetric, 0: non-symmetric, -1: skew-symmetric", &tsic);

  setStringToIntegralParameter<int>("NITSCHE_WEIGHTING_TSI", "harmonic",
      "how to weight consistency terms in Nitsche contact formulation",
      tuple<std::string>("slave", "master", "harmonic", "physical"),
      tuple<int>(Inpar::CONTACT::NitWgt_slave, Inpar::CONTACT::NitWgt_master,
          Inpar::CONTACT::NitWgt_harmonic, Inpar::CONTACT::NitWgt_phyiscal),
      &tsic);

  Core::UTILS::BoolParameter("NITSCHE_PENALTY_ADAPTIVE_TSI", "yes",
      "adapt penalty parameter after each converged time step", &tsic);

  Core::UTILS::DoubleParameter(
      "PENALTYPARAM_THERMO", 0.0, "Penalty parameter for Nitsche solution strategy", &tsic);

  setStringToIntegralParameter<int>("NITSCHE_METHOD_TSI", "nitsche",
      "how to treat thermal interface problem: strong substitution or Nitsche for general "
      "interface conditions",
      tuple<std::string>("nitsche", "substitution"),
      tuple<int>(Inpar::CONTACT::NitThr_nitsche, Inpar::CONTACT::NitThr_substitution), &tsic);

  setStringToIntegralParameter<int>("TSI_LINE_SEARCH", "none", "line-search strategy",
      tuple<std::string>("none", "structure", "thermo", "and", "or"),
      tuple<int>(LS_none, LS_structure, LS_thermo, LS_and, LS_or), &tsidynmono);
}

FOUR_C_NAMESPACE_CLOSE

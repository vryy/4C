/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameters for contact

\level 2


*/
/*----------------------------------------------------------------------*/
#include "4C_inpar_contact.hpp"

#include "4C_inpar_structure.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::CONTACT::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  /* parameters for structural meshtying and contact */
  Teuchos::ParameterList& scontact = list->sublist("CONTACT DYNAMIC", false, "");

  Core::UTILS::IntParameter(
      "LINEAR_SOLVER", -1, "number of linear solver used for meshtying and contact", &scontact);

  Core::UTILS::BoolParameter("RESTART_WITH_CONTACT", "No",
      "Must be chosen if a non-contact simulation is to be restarted with contact", &scontact);

  setStringToIntegralParameter<int>("ADHESION", "None", "Type of adhesion law",
      tuple<std::string>("None", "none", "bounded", "b"),
      tuple<int>(adhesion_none, adhesion_none, adhesion_bound, adhesion_bound), &scontact);

  setStringToIntegralParameter<int>("FRICTION", "None", "Type of friction law",
      tuple<std::string>("None", "Stick", "Tresca", "Coulomb"),
      tuple<int>(friction_none, friction_stick, friction_tresca, friction_coulomb), &scontact);

  Core::UTILS::BoolParameter("FRLESS_FIRST", "No",
      "If chosen the first time step of a newly in contact slave node is regarded as frictionless",
      &scontact);

  Core::UTILS::BoolParameter("GP_SLIP_INCR", "No",
      "If chosen the slip increment is computed gp-wise which results to a non-objective quantity, "
      "but this would be consistent to wear and tsi calculations.",
      &scontact);

  setStringToIntegralParameter<int>("STRATEGY", "LagrangianMultipliers",
      "Type of employed solving strategy",
      tuple<std::string>("LagrangianMultipliers", "lagrange", "Lagrange", "penalty", "Penalty",
          "Uzawa", "Nitsche", "Ehl", "MultiScale"),
      tuple<int>(solution_lagmult, solution_lagmult, solution_lagmult, solution_penalty,
          solution_penalty, solution_uzawa, solution_nitsche, solution_ehl, solution_multiscale),
      &scontact);

  setStringToIntegralParameter<int>("SYSTEM", "Condensed", "Type of linear system setup / solution",
      tuple<std::string>("Condensed", "condensed", "cond", "Condensedlagmult", "condensedlagmult",
          "condlm", "SaddlePoint", "Saddlepoint", "saddlepoint", "sp", "none"),
      tuple<int>(system_condensed, system_condensed, system_condensed, system_condensed_lagmult,
          system_condensed_lagmult, system_condensed_lagmult, system_saddlepoint,
          system_saddlepoint, system_saddlepoint, system_saddlepoint, system_none),
      &scontact);

  Core::UTILS::DoubleParameter("PENALTYPARAM", 0.0,
      "Penalty parameter for penalty / Uzawa augmented solution strategy", &scontact);
  Core::UTILS::DoubleParameter("PENALTYPARAMTAN", 0.0,
      "Tangential penalty parameter for penalty / Uzawa augmented solution strategy", &scontact);
  Core::UTILS::IntParameter(
      "UZAWAMAXSTEPS", 10, "Maximum no. of Uzawa steps for Uzawa solution strategy", &scontact);
  Core::UTILS::DoubleParameter("UZAWACONSTRTOL", 1.0e-8,
      "Tolerance of constraint norm for Uzawa solution strategy", &scontact);

  Core::UTILS::BoolParameter(
      "SEMI_SMOOTH_NEWTON", "Yes", "If chosen semi-smooth Newton concept is applied", &scontact);

  Core::UTILS::DoubleParameter(
      "SEMI_SMOOTH_CN", 1.0, "Weighting factor cn for semi-smooth PDASS", &scontact);
  Core::UTILS::DoubleParameter(
      "SEMI_SMOOTH_CT", 1.0, "Weighting factor ct for semi-smooth PDASS", &scontact);

  Core::UTILS::BoolParameter("CONTACTFORCE_ENDTIME", "No",
      "If chosen, the contact force is not evaluated at the generalized midpoint, but at the end "
      "of the time step",
      &scontact);

  Core::UTILS::BoolParameter(
      "VELOCITY_UPDATE", "No", "If chosen, velocity update method is applied", &scontact);

  setStringToIntegralParameter<int>("EMOUTPUT", "None", "Type of energy and momentum output",
      tuple<std::string>(
          "None", "none", "No", "no", "Screen", "screen", "File", "file", "Both", "both"),
      tuple<int>(output_none, output_none, output_none, output_none, output_screen, output_screen,
          output_file, output_file, output_both, output_both),
      &scontact);

  Core::UTILS::BoolParameter(
      "INITCONTACTBYGAP", "No", "Initialize init contact by weighted gap vector", &scontact);

  Core::UTILS::DoubleParameter("INITCONTACTGAPVALUE", 0.0,
      "Value for initialization of init contact set with gap vector", &scontact);

  // solver convergence test parameters for contact/meshtying in saddlepoint formulation
  setStringToIntegralParameter<int>("NORMCOMBI_RESFCONTCONSTR", "And",
      "binary operator to combine contact constraints and residual force values",
      tuple<std::string>("And", "Or"), tuple<int>(Inpar::Solid::bop_and, Inpar::Solid::bop_or),
      &scontact);

  setStringToIntegralParameter<int>("NORMCOMBI_DISPLAGR", "And",
      "binary operator to combine displacement increments and Lagrange multiplier increment values",
      tuple<std::string>("And", "Or"), tuple<int>(Inpar::Solid::bop_and, Inpar::Solid::bop_or),
      &scontact);

  Core::UTILS::DoubleParameter("TOLCONTCONSTR", 1.0E-6,
      "tolerance in the contact constraint norm for the newton iteration (saddlepoint formulation "
      "only)",
      &scontact);
  Core::UTILS::DoubleParameter("TOLLAGR", 1.0E-6,
      "tolerance in the LM norm for the newton iteration (saddlepoint formulation only)",
      &scontact);

  setStringToIntegralParameter<int>("CONSTRAINT_DIRECTIONS", "ntt",
      "formulation of constraints in normal/tangential or xyz-direction",
      tuple<std::string>("ntt", "xyz"), tuple<int>(constr_ntt, constr_xyz), &scontact);

  setStringToIntegralParameter<int>("CONTACT_REGULARIZATION", "no", "use regularized contact",
      tuple<std::string>("no", "tanh"), tuple<int>(reg_none, reg_tanh), &scontact);

  Core::UTILS::BoolParameter("NONSMOOTH_GEOMETRIES", "No",
      "If chosen the contact algorithm combines mortar and nts formulations. This is needed if "
      "contact between entities of different geometric dimension (such as contact between surfaces "
      "and lines, or lines and nodes) can occur",
      &scontact);

  Core::UTILS::BoolParameter("NONSMOOTH_CONTACT_SURFACE", "No",
      "This flag is used to alter the criterion for the evaluation of the so-called qualified "
      "vectors in the case of a self contact scenario. This is needed as the standard criterion is "
      "only valid for smooth surfaces and thus has to be altered, if the surface that is defined "
      "to be a self contact surface is non-smooth!",
      &scontact);

  Core::UTILS::DoubleParameter("HYBRID_ANGLE_MIN", -1.0,
      "Non-smooth contact: angle between cpp normal and element normal: begin transition (Mortar)",
      &scontact);
  Core::UTILS::DoubleParameter("HYBRID_ANGLE_MAX", -1.0,
      "Non-smooth contact: angle between cpp normal and element normal: end transition (NTS)",
      &scontact);

  Core::UTILS::BoolParameter("CPP_NORMALS", "No",
      "If chosen the nodal normal field is created as averaged CPP normal field.", &scontact);

  Core::UTILS::BoolParameter(
      "TIMING_DETAILS", "No", "Enable and print detailed contact timings to screen.", &scontact);

  // --------------------------------------------------------------------------
  Core::UTILS::DoubleParameter(
      "NITSCHE_THETA", 0.0, "+1: symmetric, 0: non-symmetric, -1: skew-symmetric", &scontact);
  Core::UTILS::DoubleParameter("NITSCHE_THETA_2", 1.0,
      "+1: Chouly-type, 0: Burman penalty-free (only with theta=-1)", &scontact);

  setStringToIntegralParameter<int>("NITSCHE_WEIGHTING", "harmonic",
      "how to weight consistency terms in Nitsche contact formulation",
      tuple<std::string>("slave", "master", "harmonic"),
      tuple<int>(NitWgt_slave, NitWgt_master, NitWgt_harmonic), &scontact);

  Core::UTILS::BoolParameter("NITSCHE_PENALTY_ADAPTIVE", "yes",
      "adapt penalty parameter after each converged time step", &scontact);

  Core::UTILS::BoolParameter("REGULARIZED_NORMAL_CONTACT", "No",
      "add a regularized normal contact formulation", &scontact);
  Core::UTILS::DoubleParameter(
      "REGULARIZATION_THICKNESS", -1., "maximum contact penetration", &scontact);
  Core::UTILS::DoubleParameter("REGULARIZATION_STIFFNESS", -1.,
      "initial contact stiffness (i.e. initial \"penalty parameter\")", &scontact);
}

FOUR_C_NAMESPACE_CLOSE

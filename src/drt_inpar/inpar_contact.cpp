/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameters for contact

\level 2


*/
/*----------------------------------------------------------------------*/
#include "drt_validparameters.H"
#include "inpar_contact.H"
#include "inpar_structure.H"



void INPAR::CONTACT::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::Array<std::string> yesnotuple =
      tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

  /* parameters for structural meshtying and contact */
  Teuchos::ParameterList& scontact = list->sublist("CONTACT DYNAMIC", false, "");

  IntParameter(
      "LINEAR_SOLVER", -1, "number of linear solver used for meshtying and contact", &scontact);

  setStringToIntegralParameter<int>("RESTART_WITH_CONTACT", "No",
      "Must be chosen if a non-contact simulation is to be restarted with contact", yesnotuple,
      yesnovalue, &scontact);

  setStringToIntegralParameter<int>("ADHESION", "None", "Type of adhesion law",
      tuple<std::string>("None", "none", "bounded", "b"),
      tuple<int>(adhesion_none, adhesion_none, adhesion_bound, adhesion_bound), &scontact);

  setStringToIntegralParameter<int>("FRICTION", "None", "Type of friction law",
      tuple<std::string>("None", "Stick", "Tresca", "Coulomb"),
      tuple<int>(friction_none, friction_stick, friction_tresca, friction_coulomb), &scontact);

  setStringToIntegralParameter<int>("FRLESS_FIRST", "No",
      "If chosen the first time step of a newly in contact slave node is regarded as frictionless",
      yesnotuple, yesnovalue, &scontact);

  setStringToIntegralParameter<int>("GP_SLIP_INCR", "No",
      "If chosen the slip increment is computed gp-wise which results to a non-objective quantity, "
      "but this would be consistent to wear and tsi calculations.",
      yesnotuple, yesnovalue, &scontact);

  setStringToIntegralParameter<int>("STRATEGY", "LagrangianMultipliers",
      "Type of employed solving strategy",
      tuple<std::string>("LagrangianMultipliers", "lagrange", "Lagrange", "penalty", "Penalty",
          "Uzawa", "Augmented", "SteepestAscent", "SteepestAscent_SP", "StdLagrange", "Combo",
          "XContact", "Nitsche", "Ehl", "MultiScale"),
      tuple<int>(solution_lagmult, solution_lagmult, solution_lagmult, solution_penalty,
          solution_penalty, solution_uzawa, solution_augmented, solution_steepest_ascent,
          solution_steepest_ascent_sp, solution_std_lagrange, solution_combo, solution_xcontact,
          solution_nitsche, solution_ehl, solution_multiscale),
      &scontact);

  setStringToIntegralParameter<int>("SYSTEM", "Condensed", "Type of linear system setup / solution",
      tuple<std::string>("Condensed", "condensed", "cond", "Condensedlagmult", "condensedlagmult",
          "condlm", "SaddlePoint", "Saddlepoint", "saddlepoint", "sp", "none"),
      tuple<int>(system_condensed, system_condensed, system_condensed, system_condensed_lagmult,
          system_condensed_lagmult, system_condensed_lagmult, system_saddlepoint,
          system_saddlepoint, system_saddlepoint, system_saddlepoint, system_none),
      &scontact);

  DoubleParameter("PENALTYPARAM", 0.0,
      "Penalty parameter for penalty / Uzawa augmented solution strategy", &scontact);
  DoubleParameter("PENALTYPARAMTAN", 0.0,
      "Tangential penalty parameter for penalty / Uzawa augmented solution strategy", &scontact);
  IntParameter(
      "UZAWAMAXSTEPS", 10, "Maximum no. of Uzawa steps for Uzawa solution strategy", &scontact);
  DoubleParameter("UZAWACONSTRTOL", 1.0e-8,
      "Tolerance of constraint norm for Uzawa solution strategy", &scontact);

  setStringToIntegralParameter<int>("SEMI_SMOOTH_NEWTON", "Yes",
      "If chosen semi-smooth Newton concept is applied", yesnotuple, yesnovalue, &scontact);

  DoubleParameter("SEMI_SMOOTH_CN", 1.0, "Weighting factor cn for semi-smooth PDASS", &scontact);
  DoubleParameter("SEMI_SMOOTH_CT", 1.0, "Weighting factor ct for semi-smooth PDASS", &scontact);

  setStringToIntegralParameter<int>("CONTACTFORCE_ENDTIME", "No",
      "If chosen, the contact force is not evaluated at the generalized midpoint, but at the end "
      "of the time step",
      yesnotuple, yesnovalue, &scontact);

  setStringToIntegralParameter<int>("VELOCITY_UPDATE", "No",
      "If chosen, velocity update method is applied", yesnotuple, yesnovalue, &scontact);

  setStringToIntegralParameter<int>("EMOUTPUT", "None", "Type of energy and momentum output",
      tuple<std::string>(
          "None", "none", "No", "no", "Screen", "screen", "File", "file", "Both", "both"),
      tuple<int>(output_none, output_none, output_none, output_none, output_screen, output_screen,
          output_file, output_file, output_both, output_both),
      &scontact);

  setStringToIntegralParameter<int>("ERROR_NORMS", "None",
      "Choice of analytical solution for error norm computation",
      tuple<std::string>("None", "none", "No", "no", "Zero", "zero", "Bending", "bending", "Sphere",
          "sphere", "Thick", "thick", "Plate", "plate"),
      tuple<int>(errornorms_none, errornorms_none, errornorms_none, errornorms_none,
          errornorms_zero, errornorms_zero, errornorms_bending, errornorms_bending,
          errornorms_sphere, errornorms_sphere, errornorms_thicksphere, errornorms_thicksphere,
          errornorms_infiniteplate, errornorms_infiniteplate),
      &scontact);

  setStringToIntegralParameter<int>("INITCONTACTBYGAP", "No",
      "Initialize init contact by weighted gap vector", yesnotuple, yesnovalue, &scontact);

  DoubleParameter("INITCONTACTGAPVALUE", 0.0,
      "Value for initialization of init contact set with gap vector", &scontact);

  // solver convergence test parameters for contact/meshtying in saddlepoint formulation
  setStringToIntegralParameter<int>("NORMCOMBI_RESFCONTCONSTR", "And",
      "binary operator to combine contact constraints and residual force values",
      tuple<std::string>("And", "Or"), tuple<int>(INPAR::STR::bop_and, INPAR::STR::bop_or),
      &scontact);

  setStringToIntegralParameter<int>("NORMCOMBI_DISPLAGR", "And",
      "binary operator to combine displacement increments and Lagrange multiplier increment values",
      tuple<std::string>("And", "Or"), tuple<int>(INPAR::STR::bop_and, INPAR::STR::bop_or),
      &scontact);

  DoubleParameter("TOLCONTCONSTR", 1.0E-6,
      "tolerance in the contact constraint norm for the newton iteration (saddlepoint formulation "
      "only)",
      &scontact);
  DoubleParameter("TOLLAGR", 1.0E-6,
      "tolerance in the LM norm for the newton iteration (saddlepoint formulation only)",
      &scontact);

  setStringToIntegralParameter<int>("CONSTRAINT_DIRECTIONS", "ntt",
      "formulation of constraints in normal/tangential or xyz-direction",
      tuple<std::string>("ntt", "xyz"), tuple<int>(constr_ntt, constr_xyz), &scontact);

  setStringToIntegralParameter<int>("CONTACT_REGULARIZATION", "no", "use regularized contact",
      tuple<std::string>("no", "tanh"), tuple<int>(reg_none, reg_tanh), &scontact);

  setStringToIntegralParameter<int>("NONSMOOTH_GEOMETRIES", "No",
      "If chosen the contact algorithm combines mortar and nts formulations. This is needed if "
      "contact between entities of different geometric dimension (such as contact between surfaces "
      "and lines, or lines and nodes) can occur",
      yesnotuple, yesnovalue, &scontact);

  setStringToIntegralParameter<int>("NONSMOOTH_CONTACT_SURFACE", "No",
      "This flag is used to alter the criterion for the evaluation of the so-called qualified "
      "vectors in the case of a self contact scenario. This is needed as the standard criterion is "
      "only valid for smooth surfaces and thus has to be altered, if the surface that is defined "
      "to be a self contact surface is non-smooth!",
      yesnotuple, yesnovalue, &scontact);

  DoubleParameter("HYBRID_ANGLE_MIN", -1.0,
      "Non-smooth contact: angle between cpp normal and element normal: begin transition (Mortar)",
      &scontact);
  DoubleParameter("HYBRID_ANGLE_MAX", -1.0,
      "Non-smooth contact: angle between cpp normal and element normal: end transition (NTS)",
      &scontact);

  setStringToIntegralParameter<int>("CPP_NORMALS", "No",
      "If chosen the nodal normal field is created as averaged CPP normal field.", yesnotuple,
      yesnovalue, &scontact);

  BoolParameter(
      "TIMING_DETAILS", "No", "Enable and print detailed contact timings to screen.", &scontact);

  // --------------------------------------------------------------------------
  // sub-list "Augmented"
  Teuchos::ParameterList& augcontact = scontact.sublist("AUGMENTED");

  setStringToIntegralParameter<int>("PRINT_LINEAR_CONSERVATION", "No",
      "Do and print the linear momentum conservation check.", yesnotuple, yesnovalue, &augcontact);

  setStringToIntegralParameter<int>("PRINT_ANGULAR_CONSERVATION", "No",
      "Do and print the angular momentum conservation check.", yesnotuple, yesnovalue, &augcontact);

  setStringToIntegralParameter<int>("VARIATIONAL_APPROACH", "incomplete",
      "Type of employed variational approach", tuple<std::string>("complete", "incomplete"),
      tuple<int>(var_complete, var_incomplete), &augcontact);

  setStringToIntegralParameter<int>("ADD_INACTIVE_FORCE_CONTRIBUTIONS", "No",
      "Add the contribution from the inactive Lagrange multipliers to the"
      "force balance.",
      yesnotuple, yesnovalue, &augcontact);

  setStringToIntegralParameter<int>("ASSEMBLE_STRATEGY", "node_based",
      "Type of employed assemble strategy", tuple<std::string>("node_based"),
      tuple<int>(assemble_node_based), &augcontact);

  setStringToIntegralParameter<FDCheck>("FD_CHECK", "off",
      "Switch finite difference check on and off.",
      tuple<std::string>("off", "global", "gauss_point"),
      tuple<FDCheck>(FDCheck::off, FDCheck::global, FDCheck::gauss_point), &augcontact);

  IntParameter("PARALLEL_REDIST_INTERVAL", -1,
      "Specifies the Newton iteration interval in which the parallel "
      "redistribution is controlled. An interval value equal to or smaller than "
      "zero disables the dynamic redistribution control mechanism during the Newton"
      " iteration.",
      &augcontact);

  // --------------------------------------------------------------------------
  // sub-sub-list "Augmented/SteepestAscent"
  Teuchos::ParameterList& sacontact = augcontact.sublist("STEEPESTASCENT");

  DoubleParameter("CORRECTION_PARAMETER", 0.0,
      "Some penalty update methods use a user-specified correction parameter, "
      "which accelerates the convergence. This parameter can be specified here.",
      &sacontact);

  DoubleParameter("DECREASE_CORRECTION_PARAMETER", 1.0,
      "Some penalty update methods use a user-specified decrease correction "
      "parameter, which can be used to initiate a decrease of the "
      "regularization parameter in cumbersome situations.",
      &sacontact);

  setStringToIntegralParameter<PenaltyUpdate>("PENALTY_UPDATE", "Vague",
      "Which kind of "
      "penalty update should be used?",
      tuple<std::string>("Vague", "SufficientLinearReduction", "SufficientAngle", "None"),
      tuple<PenaltyUpdate>(PenaltyUpdate::vague, PenaltyUpdate::sufficient_lin_reduction,
          PenaltyUpdate::sufficient_angle, PenaltyUpdate::none),
      &sacontact);

  // --------------------------------------------------------------------------
  // sub-sub-list "Augmented/Combo"
  Teuchos::ParameterList& combo_contact = augcontact.sublist("COMBO");

  setStringToIntegralParameter<int>("STRATEGY_0", "Vague",
      "Type of"
      " first solving strategy",
      tuple<std::string>(
          "Vague", "Augmented", "StdLagrange", "SteepestAscent", "SteepestAscent_SP"),
      tuple<int>(solution_vague, solution_augmented, solution_std_lagrange,
          solution_steepest_ascent, solution_steepest_ascent_sp),
      &combo_contact);

  setStringToIntegralParameter<int>("STRATEGY_1", "Vague",
      "Type of"
      " second solving strategy",
      tuple<std::string>(
          "Vague", "Augmented", "StdLagrange", "SteepestAscent", "SteepestAscent_SP"),
      tuple<int>(solution_vague, solution_augmented, solution_std_lagrange,
          solution_steepest_ascent, solution_steepest_ascent_sp),
      &combo_contact);

  IntParameter("LINEAR_SOLVER_STRATEGY_0", -1,
      "Linear solver for STRATEGY_0 of the COMBO strategy.", &combo_contact);

  IntParameter("LINEAR_SOLVER_STRATEGY_1", -1,
      "Linear solver for STRATEGY_1 of the COMBO strategy.", &combo_contact);

  setStringToIntegralParameter<int>("SWITCHING_STRATEGY", "PreAsymptotic",
      "Type of"
      " switching strategy to switch between the different solving strategies",
      tuple<std::string>("PreAsymptotic"), tuple<int>(switch_preasymptotic), &combo_contact);

  setStringToIntegralParameter<int>("PRINT2SCREEN", "Yes",
      "Activate the screen output of the COMBO strategy and the different "
      "switching strategies.",
      yesnotuple, yesnovalue, &combo_contact);

  // --------------------------------------------------------------------------
  // sub-sub-list "Augmented/Lagrange_Multiplier_Function"
  Teuchos::ParameterList& lm_funct = augcontact.sublist("LAGRANGE_MULTIPLIER_FUNCTION");

  IntParameter("LINEAR_SOLVER", -1,
      "Linear Solver number for the "
      "least squares problem.",
      &lm_funct);

  // --------------------------------------------------------------------------
  // sub-sub-list "Augmented/Plot"
  Teuchos::ParameterList& plot_contact = augcontact.sublist("PLOT");

  IntParameter("RESOLUTION_X", 10, "Plot resolution in x/displacement direction", &plot_contact);
  IntParameter(
      "RESOLUTION_Y", 10, "Plot resolution in y/Lagrange multiplier direction", &plot_contact);

  IntParameter("STEP", -1, "Plot this step.", &plot_contact);
  IntParameter("ITER", -1, "Plot this iteration of the specified step.", &plot_contact);


  IntParameter("OUTPUT_PRECISION", 16, "Precision for scientific numbers.", &plot_contact);

  setStringToIntegralParameter<PlotSupportType>("X_TYPE", "vague",
      "Support type for the x-direction.",
      tuple<std::string>("vague", "step_length", "characteristic_element_length", "position_angle",
          "position_distance"),
      tuple<PlotSupportType>(PlotSupportType::vague, PlotSupportType::step_length,
          PlotSupportType::characteristic_element_length, PlotSupportType::position_angle,
          PlotSupportType::position_distance),
      &plot_contact);

  DoubleParameter("MIN_X", 1.0, "", &plot_contact);
  DoubleParameter("MAX_X", 1.0, "", &plot_contact);

  StringParameter(
      "FIRST_REF_POINT", "0.0 0.0 0.0", "coordinates of the first reference point", &plot_contact);

  StringParameter("SECOND_REF_POINT", "0.0 0.0 0.0", "coordinates of the second reference point",
      &plot_contact);

  setStringToIntegralParameter<PlotSupportType>("Y_TYPE", "vague",
      "Support type"
      " for the y-direction. Only relevant for multi-dimensional plots.",
      tuple<std::string>("vague", "step_length", "characteristic_element_length", "position_angle"),
      tuple<PlotSupportType>(PlotSupportType::vague, PlotSupportType::step_length,
          PlotSupportType::characteristic_element_length, PlotSupportType::position_angle),
      &plot_contact);

  DoubleParameter("MIN_Y", 1.0, "", &plot_contact);
  DoubleParameter("MAX_Y", 1.0, "", &plot_contact);

  setStringToIntegralParameter<PlotFuncName>("FUNC_NAME", "vague", "Plot type.",
      tuple<std::string>("vague", "Lagrangian", "Infeasibility", "Energy", "Energy_Gradient",
          "Weighted_Gap", "Weighted_Gap_Gradient", "Weighted_Gap_Modified_Gradient",
          "Weighted_Gap_Gradient_Error", "Weighted_Gap_Gradient_Nodal_Jacobian_Error",
          "Weighted_Gap_Gradient_Nodal_Ma_Proj_Error"),
      tuple<PlotFuncName>(PlotFuncName::vague, PlotFuncName::lagrangian,
          PlotFuncName::infeasibility, PlotFuncName::energy, PlotFuncName::energy_gradient,
          PlotFuncName::weighted_gap, PlotFuncName::weighted_gap_gradient,
          PlotFuncName::weighted_gap_mod_gradient, PlotFuncName::weighted_gap_gradient_error,
          PlotFuncName::weighted_gap_gradient_nodal_jacobian_error,
          PlotFuncName::weighted_gap_gradient_nodal_ma_proj_error),
      &plot_contact);

  setStringToIntegralParameter<PlotFileFormat>("FILE_FORMAT", "matlab",
      "Format of the written text file.", tuple<std::string>("matlab", "pgfplot"),
      tuple<PlotFileFormat>(PlotFileFormat::matlab, PlotFileFormat::pgfplot), &plot_contact);

  setStringToIntegralParameter<std::ios_base::openmode>("FILE_OPEN_MODE", "truncate",
      "File opening mode.", tuple<std::string>("truncate", "append"),
      tuple<std::ios_base::openmode>(
          std::ios_base::out | std::ios_base::trunc, std::ios_base::out | std::ios_base::app),
      &plot_contact);

  IntParameter("WGAP_NODE_GID", -1,
      "Weighted gap of the slave node with "
      "this global ID will be considered, if FUNC_NAME == \"Weighted_GAP\".",
      &plot_contact);

  setStringToIntegralParameter<PlotMode>("MODE", "off", "Plot mode.",
      tuple<std::string>("off", "write_single_iteration_of_step", "write_each_iteration_of_step",
          "write_last_iteration_of_step"),
      tuple<PlotMode>(PlotMode::off, PlotMode::write_single_iteration_of_step,
          PlotMode::write_each_iteration_of_step, PlotMode::write_last_iteration_of_step),
      &plot_contact);

  setStringToIntegralParameter<PlotType>("TYPE", "line", "Plot type.",
      tuple<std::string>("vague", "scalar", "line", "surface", "vector_field_2d"),
      tuple<PlotType>(PlotType::vague, PlotType::scalar, PlotType::line, PlotType::surface,
          PlotType::vector_field_2d),
      &plot_contact);

  setStringToIntegralParameter<PlotDirection>("DIRECTION", "vague",
      "Choose the direction for the plot.",
      tuple<std::string>("vague", "current", "read_from_file", "zero"),
      tuple<PlotDirection>(PlotDirection::vague, PlotDirection::current_search_direction,
          PlotDirection::read_from_file, PlotDirection::zero),
      &plot_contact);

  DRT::INPUT::StringParameter("DIRECTION_FILE", "none",
      "Filename of the "
      "text file containing the plot direction.",
      &plot_contact);

  setStringToIntegralParameter<PlotDirectionSplit>("DIRECTION_SPLIT", "vague",
      "Choose how to split the search direction for a multi-dimensional plot,"
      " e.g. a surface plot.",
      tuple<std::string>("vague", "displ_lm", "slave_master_displ"),
      tuple<PlotDirectionSplit>(PlotDirectionSplit::vague,
          PlotDirectionSplit::displacement_lagrange_multiplier,
          PlotDirectionSplit::slave_master_displacements),
      &plot_contact);

  setStringToIntegralParameter<PlotReferenceType>("REFERENCE_TYPE", "vague",
      "Reference state for the plot.",
      tuple<std::string>("vague", "current_solution", "previous_solution"),
      tuple<PlotReferenceType>(PlotReferenceType::vague, PlotReferenceType::current_solution,
          PlotReferenceType::previous_solution),
      &plot_contact);

  // --------------------------------------------------------------------------
  // sub-list "eXtended contact formulation"
  Teuchos::ParameterList& xcontact = scontact.sublist("XCONTACT");
  // TODO
  setStringToIntegralParameter<int>("CONST_CPP_NORMAL", "No",
      "If chosen, closest point normal on master is assumed to be constant during"
      " variation and linearization.",
      yesnotuple, yesnovalue, &xcontact);

  // TODO
  setStringToIntegralParameter<int>("H1_DUALITY_PAIRING", "Yes",
      "If chosen, H1 duality pairing for contact potential is used.", yesnotuple, yesnovalue,
      &xcontact);

  // --------------------------------------------------------------------------
  DoubleParameter(
      "NITSCHE_THETA", 0.0, "+1: symmetric, 0: non-symmetric, -1: skew-symmetric", &scontact);
  DoubleParameter("NITSCHE_THETA_2", 1.0,
      "+1: Chouly-type, 0: Burman penalty-free (only with theta=-1)", &scontact);

  setStringToIntegralParameter<int>("NITSCHE_WEIGHTING", "harmonic",
      "how to weight consistency terms in Nitsche contact formulation",
      tuple<std::string>("slave", "master", "harmonic"),
      tuple<int>(NitWgt_slave, NitWgt_master, NitWgt_harmonic), &scontact);

  setStringToIntegralParameter<int>("NITSCHE_PENALTY_ADAPTIVE", "yes",
      "adapt penalty parameter after each converged time step", yesnotuple, yesnovalue, &scontact);

  setStringToIntegralParameter<int>("REGULARIZED_NORMAL_CONTACT", "No",
      "add a regularized normal contact formulation", yesnotuple, yesnovalue, &scontact);
  DoubleParameter("REGULARIZATION_THICKNESS", -1., "maximum contact penetration", &scontact);
  DoubleParameter("REGULARIZATION_STIFFNESS", -1.,
      "initial contact stiffness (i.e. initial \"penalty parameter\")", &scontact);
}

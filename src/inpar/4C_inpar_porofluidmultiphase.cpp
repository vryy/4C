#include "4C_inpar_porofluidmultiphase.hpp"

#include "4C_inpar_bio.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

void Inpar::POROFLUIDMULTIPHASE::set_valid_parameters(Teuchos::ParameterList& list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& porofluidmultiphasedyn = list.sublist("POROFLUIDMULTIPHASE DYNAMIC",
      false, "control parameters for porofluidmultiphase problems\n");

  Core::Utils::double_parameter(
      "MAXTIME", 1000.0, "Total simulation time", &porofluidmultiphasedyn);
  Core::Utils::int_parameter("NUMSTEP", 20, "Total number of time steps", &porofluidmultiphasedyn);
  Core::Utils::double_parameter("TIMESTEP", 0.1, "Time increment dt", &porofluidmultiphasedyn);
  Core::Utils::int_parameter(
      "RESULTSEVRY", 1, "Increment for writing solution", &porofluidmultiphasedyn);
  Core::Utils::int_parameter(
      "RESTARTEVRY", 1, "Increment for writing restart", &porofluidmultiphasedyn);

  Core::Utils::double_parameter(
      "THETA", 0.5, "One-step-theta time integration factor", &porofluidmultiphasedyn);

  setStringToIntegralParameter<TimeIntegrationScheme>("TIMEINTEGR", "One_Step_Theta",
      "Time Integration Scheme", tuple<std::string>("One_Step_Theta"),
      tuple<TimeIntegrationScheme>(timeint_one_step_theta), &porofluidmultiphasedyn);

  setStringToIntegralParameter<CalcError>("CALCERROR", "No",
      "compute error compared to analytical solution",
      tuple<std::string>("No", "error_by_function"),
      tuple<CalcError>(calcerror_no, calcerror_byfunction), &porofluidmultiphasedyn);

  Core::Utils::int_parameter("CALCERRORNO", -1,
      "function number for porofluidmultiphase error computation", &porofluidmultiphasedyn);

  // linear solver id used for porofluidmultiphase problems
  Core::Utils::int_parameter("LINEAR_SOLVER", -1,
      "number of linear solver used for the porofluidmultiphase problem", &porofluidmultiphasedyn);

  Core::Utils::int_parameter(
      "ITEMAX", 10, "max. number of nonlin. iterations", &porofluidmultiphasedyn);
  Core::Utils::double_parameter("ABSTOLRES", 1e-14,
      "Absolute tolerance for deciding if residual of nonlinear problem is already zero",
      &porofluidmultiphasedyn);

  // convergence criteria adaptivity
  Core::Utils::bool_parameter("ADAPTCONV", "No",
      "Switch on adaptive control of linear solver tolerance for nonlinear solution",
      &porofluidmultiphasedyn);
  Core::Utils::double_parameter("ADAPTCONV_BETTER", 0.1,
      "The linear solver shall be this much better than the current nonlinear residual in the "
      "nonlinear convergence limit",
      &porofluidmultiphasedyn);

  // parameters for finite difference check
  setStringToIntegralParameter<FdCheck>("FDCHECK", "none",
      "flag for finite difference check: none, local, or global",
      tuple<std::string>("none",
          "global"),  // perform finite difference check on time integrator level
      tuple<FdCheck>(fdcheck_none, fdcheck_global), &porofluidmultiphasedyn);
  Core::Utils::double_parameter("FDCHECKEPS", 1.e-6,
      "dof perturbation magnitude for finite difference check (1.e-6 seems to work very well, "
      "whereas smaller values don't)",
      &porofluidmultiphasedyn);
  Core::Utils::double_parameter("FDCHECKTOL", 1.e-6,
      "relative tolerance for finite difference check", &porofluidmultiphasedyn);
  Core::Utils::bool_parameter("SKIPINITDER", "yes",
      "Flag to skip computation of initial time derivative", &porofluidmultiphasedyn);
  Core::Utils::bool_parameter("OUTPUT_SATANDPRESS", "yes",
      "Flag if output of saturations and pressures should be calculated", &porofluidmultiphasedyn);
  Core::Utils::bool_parameter("OUTPUT_SOLIDPRESS", "yes",
      "Flag if output of solid pressure should be calculated", &porofluidmultiphasedyn);
  Core::Utils::bool_parameter("OUTPUT_POROSITY", "yes",
      "Flag if output of porosity should be calculated", &porofluidmultiphasedyn);
  Core::Utils::bool_parameter("OUTPUT_PHASE_VELOCITIES", "yes",
      "Flag if output of phase velocities should be calculated", &porofluidmultiphasedyn);

  // Biot stabilization
  Core::Utils::bool_parameter(
      "STAB_BIOT", "No", "Flag to (de)activate BIOT stabilization.", &porofluidmultiphasedyn);
  Core::Utils::double_parameter("STAB_BIOT_SCALING", 1.0,
      "Scaling factor for stabilization parameter for biot stabilization of porous flow.",
      &porofluidmultiphasedyn);

  setStringToIntegralParameter<VectorNorm>("VECTORNORM_RESF", "L2",
      "type of norm to be applied to residuals",
      tuple<std::string>("L1", "L1_Scaled", "L2", "Rms", "Inf"),
      tuple<VectorNorm>(Inpar::POROFLUIDMULTIPHASE::norm_l1,
          Inpar::POROFLUIDMULTIPHASE::norm_l1_scaled, Inpar::POROFLUIDMULTIPHASE::norm_l2,
          Inpar::POROFLUIDMULTIPHASE::norm_rms, Inpar::POROFLUIDMULTIPHASE::norm_inf),
      &porofluidmultiphasedyn);

  setStringToIntegralParameter<VectorNorm>("VECTORNORM_INC", "L2",
      "type of norm to be applied to residuals",
      tuple<std::string>("L1", "L1_Scaled", "L2", "Rms", "Inf"),
      tuple<VectorNorm>(Inpar::POROFLUIDMULTIPHASE::norm_l1,
          Inpar::POROFLUIDMULTIPHASE::norm_l1_scaled, Inpar::POROFLUIDMULTIPHASE::norm_l2,
          Inpar::POROFLUIDMULTIPHASE::norm_rms, Inpar::POROFLUIDMULTIPHASE::norm_inf),
      &porofluidmultiphasedyn);

  // Iterationparameters
  Core::Utils::double_parameter("TOLRES", 1e-6,
      "tolerance in the residual norm for the Newton iteration", &porofluidmultiphasedyn);
  Core::Utils::double_parameter("TOLINC", 1e-6,
      "tolerance in the increment norm for the Newton iteration", &porofluidmultiphasedyn);

  setStringToIntegralParameter<InitialField>("INITIALFIELD", "zero_field",
      "Initial Field for transport problem",
      tuple<std::string>("zero_field", "field_by_function", "field_by_condition"),
      tuple<InitialField>(
          initfield_zero_field, initfield_field_by_function, initfield_field_by_condition),
      &porofluidmultiphasedyn);

  Core::Utils::int_parameter("INITFUNCNO", -1, "function number for scalar transport initial field",
      &porofluidmultiphasedyn);

  setStringToIntegralParameter<DivContAct>("DIVERCONT", "stop",
      "What to do with time integration when Newton-Raphson iteration failed",
      tuple<std::string>("stop", "continue"), tuple<DivContAct>(divcont_stop, divcont_continue),
      &porofluidmultiphasedyn);

  Core::Utils::int_parameter("FLUX_PROJ_SOLVER", -1,
      "Number of linear solver used for L2 projection", &porofluidmultiphasedyn);

  setStringToIntegralParameter<FluxReconstructionMethod>("FLUX_PROJ_METHOD", "none",
      "Flag to (de)activate flux reconstruction.", tuple<std::string>("none", "L2_projection"),
      tuple<std::string>("no gradient reconstruction", "gracient reconstruction via l2-projection"),
      tuple<FluxReconstructionMethod>(
          gradreco_none,  // no convective streamline edge-based stabilization
          gradreco_l2     // pressure edge-based stabilization as ghost penalty around cut elements
          ),
      &porofluidmultiphasedyn);

  // functions used for domain integrals
  setNumericStringParameter(
      "DOMAININT_FUNCT", "-1.0", "functions used for domain integrals", &porofluidmultiphasedyn);

  // coupling with 1D artery network active
  Core::Utils::bool_parameter(
      "ARTERY_COUPLING", "No", "Coupling with 1D blood vessels.", &porofluidmultiphasedyn);

  Core::Utils::double_parameter("STARTING_DBC_TIME_END", -1.0,
      "End time for the starting Dirichlet BC.", &porofluidmultiphasedyn);

  setNumericStringParameter("STARTING_DBC_ONOFF", "0",
      "Switching the starting Dirichlet BC on or off.", &porofluidmultiphasedyn);

  setNumericStringParameter("STARTING_DBC_FUNCT", "0",
      "Function prescribing the starting Dirichlet BC.", &porofluidmultiphasedyn);

  // ----------------------------------------------------------------------
  // artery mesh tying
  Teuchos::ParameterList& porofluidmultiphasemshtdyn =
      porofluidmultiphasedyn.sublist("ARTERY COUPLING", false, "Parameters for artery mesh tying");

  // maximum number of segments per artery element for 1D-3D artery coupling
  Core::Utils::int_parameter("MAXNUMSEGPERARTELE", 5,
      "maximum number of segments per artery element for 1D-3D artery coupling",
      &porofluidmultiphasemshtdyn);

  // penalty parameter
  Core::Utils::double_parameter(
      "PENALTY", 1000.0, "Penalty parameter for line-based coupling", &porofluidmultiphasemshtdyn);

  setStringToIntegralParameter<Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod>(
      "ARTERY_COUPLING_METHOD", "None", "Coupling method for artery coupling.",
      tuple<std::string>("None", "Nodal", "GPTS", "MP", "NTP"),
      tuple<std::string>("none", "Nodal Coupling", "Gauss-Point-To-Segment Approach",
          "Mortar Penalty Approach", "1D node-to-point in 2D/3D Approach"),
      tuple<Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod>(
          Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::none,   // none
          Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::nodal,  // Nodal Coupling
          Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::
              gpts,  // Gauss-Point-To-Segment
          Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::mp,  // Mortar Penalty
          Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::ntp  // 1Dnode-to-point in
                                                                               // 2D/3D
          ),
      &porofluidmultiphasemshtdyn);

  // coupled artery dofs for mesh tying
  setNumericStringParameter(
      "COUPLEDDOFS_ART", "-1.0", "coupled artery dofs for mesh tying", &porofluidmultiphasemshtdyn);

  // coupled porofluid dofs for mesh tying
  setNumericStringParameter("COUPLEDDOFS_PORO", "-1.0", "coupled porofluid dofs for mesh tying",
      &porofluidmultiphasemshtdyn);

  // functions for coupling (artery part)
  setNumericStringParameter(
      "REACFUNCT_ART", "-1", "functions for coupling (artery part)", &porofluidmultiphasemshtdyn);

  // scale for coupling (artery part)
  setNumericStringParameter(
      "SCALEREAC_ART", "0", "scale for coupling (artery part)", &porofluidmultiphasemshtdyn);

  // functions for coupling (porofluid part)
  setNumericStringParameter("REACFUNCT_CONT", "-1", "functions for coupling (porofluid part)",
      &porofluidmultiphasemshtdyn);

  // scale for coupling (porofluid part)
  setNumericStringParameter(
      "SCALEREAC_CONT", "0", "scale for coupling (porofluid part)", &porofluidmultiphasemshtdyn);

  // Flag if artery elements are evaluated in reference or current configuration
  Core::Utils::bool_parameter("EVALUATE_IN_REF_CONFIG", "yes",
      "Flag if artery elements are evaluated in reference or current configuration",
      &porofluidmultiphasemshtdyn);

  // Flag if 1D-3D coupling should be evaluated on lateral (cylinder) surface of embedded artery
  // elements
  Core::Utils::bool_parameter("LATERAL_SURFACE_COUPLING", "no",
      "Flag if 1D-3D coupling should be evaluated on lateral (cylinder) surface of embedded artery "
      "elements",
      &porofluidmultiphasemshtdyn);

  // Number of integration patches per 1D element in axial direction for lateral surface coupling
  Core::Utils::int_parameter("NUMPATCH_AXI", 1,
      "Number of integration patches per 1D element in axial direction for lateral surface "
      "coupling",
      &porofluidmultiphasemshtdyn);

  // Number of integration patches per 1D element in radial direction for lateral surface coupling
  Core::Utils::int_parameter("NUMPATCH_RAD", 1,
      "Number of integration patches per 1D element in radial direction for lateral surface "
      "coupling",
      &porofluidmultiphasemshtdyn);

  // Flag if blood vessel volume fraction should be output
  Core::Utils::bool_parameter("OUTPUT_BLOODVESSELVOLFRAC", "no",
      "Flag if output of blood vessel volume fraction should be calculated",
      &porofluidmultiphasemshtdyn);

  // Flag if summary of coupling-pairs should be printed
  Core::Utils::bool_parameter("PRINT_OUT_SUMMARY_PAIRS", "no",
      "Flag if summary of coupling-pairs should be printed", &porofluidmultiphasemshtdyn);

  // Flag if free-hanging elements (after blood vessel collapse) should be deleted
  Core::Utils::bool_parameter("DELETE_FREE_HANGING_ELES", "no",
      "Flag if free-hanging elements (after blood vessel collapse) should be deleted",
      &porofluidmultiphasemshtdyn);

  // components whose size is smaller than this fraction of the total network size are also deleted
  Core::Utils::double_parameter("DELETE_SMALL_FREE_HANGING_COMPS", -1.0,
      "Small connected components whose size is smaller than this fraction of the overall network "
      "size are additionally deleted (a valid choice of this parameter should lie between 0 and 1)",
      &porofluidmultiphasemshtdyn);
}

FOUR_C_NAMESPACE_CLOSE

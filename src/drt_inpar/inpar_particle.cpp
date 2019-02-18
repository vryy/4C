/*---------------------------------------------------------------------------*/
/*!
\file inpar_particle.cpp

\brief input parameters for particle problems

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include "inpar_particle.H"
#include "drt_validparameters.H"
#include "inpar_parameterlist_utils.H"

/*---------------------------------------------------------------------------*
 | set the particle parameters                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void INPAR::PARTICLE::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  /*-------------------------------------------------------------------------*
   | general control parameters for particle simulations                     |
   *-------------------------------------------------------------------------*/
  Teuchos::ParameterList& particledyn =
      list->sublist("PARTICLE DYNAMIC", false, "control parameters for particle simulations\n");

  // type of particle time integration
  setStringToIntegralParameter<int>("DYNAMICTYP", "VelocityVerlet",
      "type of particle time integration",
      tuple<std::string>("SemiImplicitEuler", "VelocityVerlet"),
      tuple<int>(INPAR::PARTICLE::dyna_semiimpliciteuler, INPAR::PARTICLE::dyna_velocityverlet),
      &particledyn);

  // type of particle interaction
  setStringToIntegralParameter<int>("INTERACTION", "None", "type of particle interaction",
      tuple<std::string>("None", "SPH", "DEM"),
      tuple<int>(INPAR::PARTICLE::interaction_none, INPAR::PARTICLE::interaction_sph,
          INPAR::PARTICLE::interaction_dem),
      &particledyn);

  // output type
  IntParameter(
      "RESULTSEVRY", 1, "write particle runtime output every RESULTSEVRY steps", &particledyn);
  IntParameter("RESTARTEVRY", 1, "write restart possibility every RESTARTEVRY steps", &particledyn);

  // data format for written numeric data via vtp
  setStringToIntegralParameter<int>("OUTPUT_DATA_FORMAT", "Binary",
      "data format for written numeric data", tuple<std::string>("Binary", "ASCII"),
      tuple<int>(INPAR::PARTICLE::binary, INPAR::PARTICLE::ascii), &particledyn);

  // write ghosted particles
  BoolParameter(
      "WRITE_GHOSTED_PARTICLES", "no", "write ghosted particles (debug feature)", &particledyn);

  // time loop control
  DoubleParameter("TIMESTEP", 0.01, "time step size", &particledyn);
  IntParameter("NUMSTEP", 100, "maximum number of steps", &particledyn);
  DoubleParameter("MAXTIME", 1.0, "maximum time", &particledyn);

  // gravity acceleration control
  setNumericStringParameter(
      "GRAVITY_ACCELERATION", "0.0 0.0 0.0", "acceleration due to gravity", &particledyn);
  IntParameter("GRAVITY_RAMP_FUNCT", -1, "number of function governing gravity ramp", &particledyn);

  // viscous damping factor
  DoubleParameter("VISCOUS_DAMPING", -1.0,
      "apply viscous damping force to determine static equilibrium solutions", &particledyn);

  // transfer particles to new bins every time step
  BoolParameter(
      "TRANSFER_EVERY", "no", "transfer particles to new bins every time step", &particledyn);

  // relate particle phase to material id
  StringParameter("PHASE_TO_MATERIAL_ID", "", "relate particle phase to material id", &particledyn);

  // amplitude of noise added to initial position
  DoubleParameter("INITIAL_POSITION_AMPLITUDE", 0.0, "amplitude of noise added to initial position",
      &particledyn);

  // type of particle wall source
  setStringToIntegralParameter<int>("PARTICLE_WALL_SOURCE", "NoParticleWall",
      "type of particle wall source",
      tuple<std::string>("NoParticleWall", "DiscretCondition", "BoundingBox"),
      tuple<int>(INPAR::PARTICLE::NoParticleWall, INPAR::PARTICLE::DiscretCondition,
          INPAR::PARTICLE::BoundingBox),
      &particledyn);

  /*-------------------------------------------------------------------------*
   | control parameters for initial/boundary conditions                      |
   *-------------------------------------------------------------------------*/
  Teuchos::ParameterList& particledynconditions =
      particledyn.sublist("INITIAL AND BOUNDARY CONDITIONS", false,
          "control parameters for initial/boundary conditions in particle simulations\n");

  // initial velocity field of particle phase given by function
  StringParameter("INITIAL_VELOCITY_FIELD", "",
      "initial velocity field of particle phase given by function", &particledynconditions);

  // initial acceleration field of particle phase given by function
  StringParameter("INITIAL_ACCELERATION_FIELD", "",
      "initial acceleration field of particle phase given by function", &particledynconditions);

  // dirichlet boundary condition of particle phase given by function
  StringParameter("DIRICHLET_BOUNDARY_CONDITION", "",
      "dirichlet boundary condition of particle phase given by function", &particledynconditions);

  // temperature boundary condition of particle phase given by function
  StringParameter("TEMPERATURE_BOUNDARY_CONDITION", "",
      "temperature boundary condition of particle phase given by function", &particledynconditions);

  /*-------------------------------------------------------------------------*
   | smoothed particle hydrodynamics (SPH) specific control parameters       |
   *-------------------------------------------------------------------------*/
  Teuchos::ParameterList& particledynsph = particledyn.sublist(
      "SPH", false, "control parameters for smoothed particle hydrodynamics (SPH) simulations\n");

  // type of smoothed particle hydrodynamics kernel
  setStringToIntegralParameter<int>("KERNEL", "CubicSpline",
      "type of smoothed particle hydrodynamics kernel",
      tuple<std::string>("CubicSpline", "QuinticSpline"),
      tuple<int>(INPAR::PARTICLE::CubicSpline, INPAR::PARTICLE::QuinticSpline), &particledynsph);

  // kernel space dimension number
  setStringToIntegralParameter<int>("KERNEL_SPACE_DIM", "Kernel3D", "kernel space dimension number",
      tuple<std::string>("Kernel1D", "Kernel2D", "Kernel3D"),
      tuple<int>(INPAR::PARTICLE::Kernel1D, INPAR::PARTICLE::Kernel2D, INPAR::PARTICLE::Kernel3D),
      &particledynsph);

  // type of smoothed particle hydrodynamics equation of state
  setStringToIntegralParameter<int>("EQUATIONOFSTATE", "GenTait",
      "type of smoothed particle hydrodynamics equation of state",
      tuple<std::string>("GenTait", "IdealGas"),
      tuple<int>(INPAR::PARTICLE::GenTait, INPAR::PARTICLE::IdealGas), &particledynsph);

  // type of smoothed particle hydrodynamics momentum formulation
  setStringToIntegralParameter<int>("MOMENTUMFORMULATION", "AdamiMomentumFormulation",
      "type of smoothed particle hydrodynamics momentum formulation",
      tuple<std::string>("AdamiMomentumFormulation", "MonaghanMomentumFormulation"),
      tuple<int>(
          INPAR::PARTICLE::AdamiMomentumFormulation, INPAR::PARTICLE::MonaghanMomentumFormulation),
      &particledynsph);

  // type of density evaluation scheme
  setStringToIntegralParameter<int>("DENSITYEVALUATION", "DensitySummation",
      "type of density evaluation scheme",
      tuple<std::string>("DensitySummation", "DensityIntegration", "DensityPredictCorrect"),
      tuple<int>(INPAR::PARTICLE::DensitySummation, INPAR::PARTICLE::DensityIntegration,
          INPAR::PARTICLE::DensityPredictCorrect),
      &particledynsph);

  // type of density correction scheme
  setStringToIntegralParameter<int>("DENSITYCORRECTION", "NoCorrection",
      "type of density correction scheme",
      tuple<std::string>(
          "NoCorrection", "InteriorCorrection", "NormalizedCorrection", "RandlesCorrection"),
      tuple<int>(INPAR::PARTICLE::NoCorrection, INPAR::PARTICLE::InteriorCorrection,
          INPAR::PARTICLE::NormalizedCorrection, INPAR::PARTICLE::RandlesCorrection),
      &particledynsph);


  //! type of boundary particle formulation
  setStringToIntegralParameter<int>("BOUNDARYPARTICLEFORMULATION", "NoBoundaryFormulation",
      "type of boundary particle formulation",
      tuple<std::string>("NoBoundaryFormulation", "AdamiBoundaryFormulation"),
      tuple<int>(INPAR::PARTICLE::NoBoundaryFormulation, INPAR::PARTICLE::AdamiBoundaryFormulation),
      &particledynsph);

  // type of boundary particle interaction
  setStringToIntegralParameter<int>("BOUNDARYPARTICLEINTERACTION", "NoSlipBoundaryParticle",
      "type of boundary particle interaction",
      tuple<std::string>("NoSlipBoundaryParticle", "FreeSlipBoundaryParticle"),
      tuple<int>(
          INPAR::PARTICLE::NoSlipBoundaryParticle, INPAR::PARTICLE::FreeSlipBoundaryParticle),
      &particledynsph);

  DoubleParameter("CONSISTENTPROBLEMVOLUME", 0.0,
      "prescribe problem volume filled by (non-boundary) particles to consistently initialize "
      "particle masses",
      &particledynsph);

  // type of transport velocity formulation
  setStringToIntegralParameter<int>("TRANSPORTVELOCITYFORMULATION", "NoTransportVelocity",
      "type of transport velocity formulation",
      tuple<std::string>(
          "NoTransportVelocity", "StandardTransportVelocity", "GeneralizedTransportVelocity"),
      tuple<int>(INPAR::PARTICLE::NoTransportVelocity, INPAR::PARTICLE::StandardTransportVelocity,
          INPAR::PARTICLE::GeneralizedTransportVelocity),
      &particledynsph);

  BoolParameter("NO_RELVEL_TERM", "no",
      "do not apply convection of momentum with relative velocity in case of transport velocity "
      "formulation",
      &particledynsph);

  //! type of temperature evaluation scheme
  setStringToIntegralParameter<int>("TEMPERATUREEVALUATION", "NoTemperatureEvaluation",
      "type of temperature evaluation scheme",
      tuple<std::string>("NoTemperatureEvaluation", "TemperatureIntegration"),
      tuple<int>(INPAR::PARTICLE::NoTemperatureEvaluation, INPAR::PARTICLE::TemperatureIntegration),
      &particledynsph);

  // type of surface tension formulation
  setStringToIntegralParameter<int>("SURFACETENSIONFORMULATION", "NoSurfaceTension",
      "type of surface tension formulation",
      tuple<std::string>("NoSurfaceTension", "ContinuumSurfaceForce"),
      tuple<int>(INPAR::PARTICLE::NoSurfaceTension, INPAR::PARTICLE::ContinuumSurfaceForce),
      &particledynsph);

  IntParameter("SURFACETENSION_RAMP_FUNCT", -1, "number of function governing surface tension ramp",
      &particledynsph);

  DoubleParameter("SURFACETENSIONCOEFFICIENT", -1.0,
      "surface tension coefficient in continuum surface force formulation", &particledynsph);
  DoubleParameter("STATICCONTACTANGLE", 0.0,
      "static contact angle in degree in continuum surface force formulation with wetting effects",
      &particledynsph);

  // type of phase change
  setStringToIntegralParameter<int>("PHASECHANGETYPE", "NoPhaseChange", "type of phase change",
      tuple<std::string>("NoPhaseChange", "TwoWayScalarPhaseChange"),
      tuple<int>(INPAR::PARTICLE::NoPhaseChange, INPAR::PARTICLE::TwoWayScalarPhaseChange),
      &particledynsph);

  // definition of phase change
  StringParameter("PHASECHANGEDEFINITION", "", "phase change definition", &particledynsph);

  /*-------------------------------------------------------------------------*
   | discrete element method (DEM) specific control parameters               |
   *-------------------------------------------------------------------------*/
  Teuchos::ParameterList& particledyndem = particledyn.sublist(
      "DEM", false, "control parameters for discrete element method (DEM) simulations\n");

  // type of normal contact law
  setStringToIntegralParameter<int>("NORMALCONTACTLAW", "LinearSpring",
      "normal contact law for particles",
      tuple<std::string>(
          "LinearSpring", "LinearSpringDamp", "Hertz", "LeeHerrmann", "KuwabaraKono", "Tsuji"),
      tuple<int>(INPAR::PARTICLE::LinSpring, INPAR::PARTICLE::LinSpringDamp, INPAR::PARTICLE::Hertz,
          INPAR::PARTICLE::LeeHerrmann, INPAR::PARTICLE::KuwabaraKono, INPAR::PARTICLE::Tsuji),
      &particledyndem);

  DoubleParameter("MAX_RADIUS", -1.0, "maximum expected particle radius", &particledyndem);
  DoubleParameter("MAX_VELOCITY", -1.0, "maximum expected particle velocity", &particledyndem);
  DoubleParameter("REL_PENETRATION", -1.0,
      "maximum allowed relative penetration (particle-particle)", &particledyndem);
  DoubleParameter(
      "NORMAL_STIFF", -1.0, "normal contact stiffness (particle-particle)", &particledyndem);
  DoubleParameter(
      "NORMAL_DAMP", -1.0, "normal contact damping parameter (particle-particle)", &particledyndem);
  DoubleParameter(
      "COEFF_RESTITUTION", -1.0, "coefficient of restitution (particle-particle)", &particledyndem);
  DoubleParameter("DAMP_REG_FAC", -1.0,
      "linearly regularized damping normal force in the interval |g| < (DAMP_REG_FAC * r_min)",
      &particledyndem);
  BoolParameter("TENSION_CUTOFF", "no", "switch on/off tension cutoff", &particledyndem);
}

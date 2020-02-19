/*---------------------------------------------------------------------------*/
/*! \file
\brief input parameters for particle problems

\level 3

\maintainer  Sebastian Fuchs

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include "inpar_particle.H"
#include "drt_validparameters.H"
#include "inpar_parameterlist_utils.H"

#include "../drt_lib/drt_conditiondefinition.H"

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
  StringParameter(
      "GRAVITY_ACCELERATION", "0.0 0.0 0.0", "acceleration due to gravity", &particledyn);
  IntParameter("GRAVITY_RAMP_FUNCT", -1, "number of function governing gravity ramp", &particledyn);

  // viscous damping factor
  DoubleParameter("VISCOUS_DAMPING", -1.0,
      "apply viscous damping force to determine static equilibrium solutions", &particledyn);

  // transfer particles to new bins every time step
  BoolParameter(
      "TRANSFER_EVERY", "no", "transfer particles to new bins every time step", &particledyn);

  // considered particle phases with dynamic load balance weighting factor
  StringParameter("PHASE_TO_DYNLOADBALFAC", "none",
      "considered particle phases with dynamic load balance weighting factor", &particledyn);

  // relate particle phase to material id
  StringParameter(
      "PHASE_TO_MATERIAL_ID", "none", "relate particle phase to material id", &particledyn);

  // amplitude of noise added to initial position for each spatial direction
  StringParameter("INITIAL_POSITION_AMPLITUDE", "0.0 0.0 0.0",
      "amplitude of noise added to initial position for each spatial direction", &particledyn);

  // type of particle wall source
  setStringToIntegralParameter<int>("PARTICLE_WALL_SOURCE", "NoParticleWall",
      "type of particle wall source",
      tuple<std::string>("NoParticleWall", "DiscretCondition", "BoundingBox"),
      tuple<int>(INPAR::PARTICLE::NoParticleWall, INPAR::PARTICLE::DiscretCondition,
          INPAR::PARTICLE::BoundingBox),
      &particledyn);

  // material id for particle wall from bounding box source
  IntParameter("PARTICLE_WALL_MAT", -1, "material id for particle wall from bounding box source",
      &particledyn);

  // flags defining considered states of particle wall
  BoolParameter("PARTICLE_WALL_MOVING", "no", "consider a moving particle wall", &particledyn);
  BoolParameter("PARTICLE_WALL_LOADED", "no", "consider loading on particle wall", &particledyn);

  /*-------------------------------------------------------------------------*
   | control parameters for initial/boundary conditions                      |
   *-------------------------------------------------------------------------*/
  Teuchos::ParameterList& particledynconditions =
      particledyn.sublist("INITIAL AND BOUNDARY CONDITIONS", false,
          "control parameters for initial/boundary conditions in particle simulations\n");

  // initial temperature field of particle phase given by function
  StringParameter("INITIAL_TEMP_FIELD", "none",
      "initial temperature field of particle phase given by function", &particledynconditions);

  // initial velocity field of particle phase given by function
  StringParameter("INITIAL_VELOCITY_FIELD", "none",
      "initial velocity field of particle phase given by function", &particledynconditions);

  // initial acceleration field of particle phase given by function
  StringParameter("INITIAL_ACCELERATION_FIELD", "none",
      "initial acceleration field of particle phase given by function", &particledynconditions);

  // dirichlet boundary condition of particle phase given by function
  StringParameter("DIRICHLET_BOUNDARY_CONDITION", "none",
      "dirichlet boundary condition of particle phase given by function", &particledynconditions);

  // temperature boundary condition of particle phase given by function
  StringParameter("TEMPERATURE_BOUNDARY_CONDITION", "none",
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

  DoubleParameter("INITIALPARTICLESPACING", 0.0, "initial spacing of particles", &particledynsph);

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

  // type of boundary particle formulation
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

  // type of temperature evaluation scheme
  setStringToIntegralParameter<int>("TEMPERATUREEVALUATION", "NoTemperatureEvaluation",
      "type of temperature evaluation scheme",
      tuple<std::string>("NoTemperatureEvaluation", "TemperatureIntegration"),
      tuple<int>(INPAR::PARTICLE::NoTemperatureEvaluation, INPAR::PARTICLE::TemperatureIntegration),
      &particledynsph);

  BoolParameter("TEMPERATUREGRADIENT", "no", "evaluate temperature gradient", &particledynsph);

  // type of heat source
  setStringToIntegralParameter<int>("HEATSOURCETYPE", "NoHeatSource", "type of heat source",
      tuple<std::string>("NoHeatSource", "VolumeHeatSource", "SurfaceHeatSource"),
      tuple<int>(INPAR::PARTICLE::NoHeatSource, INPAR::PARTICLE::VolumeHeatSource,
          INPAR::PARTICLE::SurfaceHeatSource),
      &particledynsph);

  IntParameter("HEATSOURCE_FUNCT", -1, "number of function governing heat source", &particledynsph);

  StringParameter(
      "HEATSOURCE_DIRECTION", "0.0 0.0 0.0", "direction of surface heat source", &particledynsph);

  // type of surface tension formulation
  setStringToIntegralParameter<int>("SURFACETENSIONFORMULATION", "NoSurfaceTension",
      "type of surface tension formulation",
      tuple<std::string>("NoSurfaceTension", "ContinuumSurfaceForce"),
      tuple<int>(INPAR::PARTICLE::NoSurfaceTension, INPAR::PARTICLE::ContinuumSurfaceForce),
      &particledynsph);

  IntParameter("SURFACETENSION_RAMP_FUNCT", -1, "number of function governing surface tension ramp",
      &particledynsph);

  DoubleParameter("SURFACETENSIONCOEFFICIENT", -1.0,
      "constant factor of surface tension coefficient", &particledynsph);
  DoubleParameter("STATICCONTACTANGLE", 0.0,
      "static contact angle in degree in continuum surface force formulation with wetting effects",
      &particledynsph);
  DoubleParameter("SURFACETENSIONTEMPFAC", 0.0,
      "factor of dependence of surface tension coefficient on temperature", &particledynsph);
  DoubleParameter("SURFACETENSIONREFTEMP", 0.0,
      "reference temperature for surface tension coefficient", &particledynsph);

  // type of dirichlet open boundary
  setStringToIntegralParameter<int>("DIRICHLETBOUNDARYTYPE", "NoDirichletOpenBoundary",
      "type of dirichlet open boundary",
      tuple<std::string>("NoDirichletOpenBoundary", "DirichletNormalToPlane"),
      tuple<int>(INPAR::PARTICLE::NoDirichletOpenBoundary, INPAR::PARTICLE::DirichletNormalToPlane),
      &particledynsph);

  IntParameter("DIRICHLET_FUNCT", -1,
      "number of function governing velocity condition on dirichlet open boundary",
      &particledynsph);

  StringParameter("DIRICHLET_OUTWARD_NORMAL", "0.0 0.0 0.0",
      "direction of outward normal on dirichlet open boundary", &particledynsph);
  StringParameter("DIRICHLET_PLANE_POINT", "0.0 0.0 0.0", "point on dirichlet open boundary plane",
      &particledynsph);

  // type of neumann open boundary
  setStringToIntegralParameter<int>("NEUMANNBOUNDARYTYPE", "NoNeumannOpenBoundary",
      "type of neumann open boundary",
      tuple<std::string>("NoNeumannOpenBoundary", "NeumannNormalToPlane"),
      tuple<int>(INPAR::PARTICLE::NoNeumannOpenBoundary, INPAR::PARTICLE::NeumannNormalToPlane),
      &particledynsph);

  IntParameter("NEUMANN_FUNCT", -1,
      "number of function governing pressure condition on neumann open boundary", &particledynsph);

  StringParameter("NEUMANN_OUTWARD_NORMAL", "0.0 0.0 0.0",
      "direction of outward normal on neumann open boundary", &particledynsph);
  StringParameter("NEUMANN_PLANE_POINT", "0.0 0.0 0.0", "point on neumann open boundary plane",
      &particledynsph);

  // type of phase change
  setStringToIntegralParameter<int>("PHASECHANGETYPE", "NoPhaseChange", "type of phase change",
      tuple<std::string>("NoPhaseChange", "TwoWayScalarPhaseChange"),
      tuple<int>(INPAR::PARTICLE::NoPhaseChange, INPAR::PARTICLE::TwoWayScalarPhaseChange),
      &particledynsph);

  // definition of phase change
  StringParameter("PHASECHANGEDEFINITION", "none", "phase change definition", &particledynsph);

  /*-------------------------------------------------------------------------*
   | discrete element method (DEM) specific control parameters               |
   *-------------------------------------------------------------------------*/
  Teuchos::ParameterList& particledyndem = particledyn.sublist(
      "DEM", false, "control parameters for discrete element method (DEM) simulations\n");

  // write particle-wall interaction output
  BoolParameter("WRITE_PARTICLE_WALL_INTERACTION", "no", "write particle-wall interaction output",
      &particledyndem);

  // type of normal contact law
  setStringToIntegralParameter<int>("NORMALCONTACTLAW", "NormalLinearSpring",
      "normal contact law for particles",
      tuple<std::string>("NormalLinearSpring", "NormalLinearSpringDamp", "NormalHertz",
          "NormalLeeHerrmann", "NormalKuwabaraKono", "NormalTsuji"),
      tuple<int>(INPAR::PARTICLE::NormalLinSpring, INPAR::PARTICLE::NormalLinSpringDamp,
          INPAR::PARTICLE::NormalHertz, INPAR::PARTICLE::NormalLeeHerrmann,
          INPAR::PARTICLE::NormalKuwabaraKono, INPAR::PARTICLE::NormalTsuji),
      &particledyndem);

  // type of tangential contact law
  setStringToIntegralParameter<int>("TANGENTIALCONTACTLAW", "NoTangentialContact",
      "tangential contact law for particles",
      tuple<std::string>("NoTangentialContact", "TangentialLinSpringDamp"),
      tuple<int>(INPAR::PARTICLE::NoTangentialContact, INPAR::PARTICLE::TangentialLinSpringDamp),
      &particledyndem);

  // type of rolling contact law
  setStringToIntegralParameter<int>("ROLLINGCONTACTLAW", "NoRollingContact",
      "rolling contact law for particles",
      tuple<std::string>("NoRollingContact", "RollingViscous", "RollingCoulomb"),
      tuple<int>(INPAR::PARTICLE::NoRollingContact, INPAR::PARTICLE::RollingViscous,
          INPAR::PARTICLE::RollingCoulomb),
      &particledyndem);

  // type of normal adhesion law
  setStringToIntegralParameter<int>("ADHESIONLAW", "NoAdhesion",
      "type of adhesion law for particles",
      tuple<std::string>("NoAdhesion", "AdhesionVdWDMT", "AdhesionRegDMT"),
      tuple<int>(INPAR::PARTICLE::NoAdhesion, INPAR::PARTICLE::AdhesionVdWDMT,
          INPAR::PARTICLE::AdhesionRegDMT),
      &particledyndem);

  // type of (random) surface energy distribution
  setStringToIntegralParameter<int>("ADHESION_SURFACE_ENERGY_DISTRIBUTION", "ConstantSurfaceEnergy",
      "type of (random) surface energy distribution",
      tuple<std::string>("ConstantSurfaceEnergy", "NormalSurfaceEnergyDistribution",
          "LogNormalSurfaceEnergyDistribution"),
      tuple<int>(INPAR::PARTICLE::ConstantSurfaceEnergy,
          INPAR::PARTICLE::NormalSurfaceEnergyDistribution,
          INPAR::PARTICLE::LogNormalSurfaceEnergyDistribution),
      &particledyndem);

  DoubleParameter("MIN_RADIUS", 0.0, "minimum allowed particle radius", &particledyndem);
  DoubleParameter("MAX_RADIUS", 0.0, "maximum allowed particle radius", &particledyndem);
  DoubleParameter("MAX_VELOCITY", -1.0, "maximum expected particle velocity", &particledyndem);

  // type of initial particle radius assignment
  setStringToIntegralParameter<int>("INITIAL_RADIUS", "RadiusFromParticleMaterial",
      "type of initial particle radius assignment",
      tuple<std::string>("RadiusFromParticleMaterial", "RadiusFromParticleInput",
          "NormalRadiusDistribution", "LogNormalRadiusDistribution"),
      tuple<int>(INPAR::PARTICLE::RadiusFromParticleMaterial,
          INPAR::PARTICLE::RadiusFromParticleInput, INPAR::PARTICLE::NormalRadiusDistribution,
          INPAR::PARTICLE::LogNormalRadiusDistribution),
      &particledyndem);

  DoubleParameter("RADIUSDISTRIBUTION_SIGMA", -1.0, "sigma of random particle radius distribution",
      &particledyndem);

  DoubleParameter("REL_PENETRATION", -1.0, "maximum allowed relative penetration", &particledyndem);
  DoubleParameter("NORMAL_STIFF", -1.0, "normal contact stiffness", &particledyndem);
  DoubleParameter("NORMAL_DAMP", -1.0, "normal contact damping parameter", &particledyndem);
  DoubleParameter("COEFF_RESTITUTION", -1.0, "coefficient of restitution", &particledyndem);
  DoubleParameter("DAMP_REG_FAC", -1.0,
      "linearly regularized damping normal force in the interval |g| < (DAMP_REG_FAC * r_min)",
      &particledyndem);
  BoolParameter("TENSION_CUTOFF", "no", "switch on/off tension cutoff", &particledyndem);

  DoubleParameter("POISSON_RATIO", -1.0, "poisson ratio", &particledyndem);
  DoubleParameter("YOUNG_MODULUS", -1.0, "young's modulus", &particledyndem);

  DoubleParameter(
      "FRICT_COEFF_TANG", -1.0, "friction coefficient for tangential contact", &particledyndem);
  DoubleParameter(
      "FRICT_COEFF_ROLL", -1.0, "friction coefficient for rolling contact", &particledyndem);

  DoubleParameter(
      "ADHESION_DISTANCE", -1.0, "adhesion distance between interacting surfaces", &particledyndem);

  DoubleParameter(
      "ADHESION_MAX_CONTACT_PRESSURE", 0.0, "adhesion maximum contact pressure", &particledyndem);
  DoubleParameter(
      "ADHESION_MAX_CONTACT_FORCE", 0.0, "adhesion maximum contact force", &particledyndem);
  BoolParameter("ADHESION_USE_MAX_CONTACT_FORCE", "no",
      "use maximum contact force instead of maximum contact pressure", &particledyndem);

  BoolParameter(
      "ADHESION_VDW_CURVE_SHIFT", "no", "shifts van-der-Waals-curve to g = 0", &particledyndem);

  DoubleParameter(
      "ADHESION_HAMAKER", -1.0, "hamaker constant of van-der-Waals interaction", &particledyndem);

  DoubleParameter("ADHESION_SURFACE_ENERGY", -1.0,
      "adhesion surface energy for the calculation of the pull-out force", &particledyndem);
  DoubleParameter("ADHESION_SURFACE_ENERGY_DISTRIBUTION_VAR", -1.0,
      "variance of adhesion surface energy distribution", &particledyndem);
  DoubleParameter("ADHESION_SURFACE_ENERGY_DISTRIBUTION_CUTOFF_FACTOR", -1.0,
      "adhesion surface energy distribution limited by multiple of variance", &particledyndem);
  DoubleParameter("ADHESION_SURFACE_ENERGY_FACTOR", 1.0,
      "factor to calculate minimum adhesion surface energy", &particledyndem);
}

/*---------------------------------------------------------------------------*
 | set the particle conditions                                sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
void INPAR::PARTICLE::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  /*-------------------------------------------------------------------------*
   | particle wall condition                                                 |
   *-------------------------------------------------------------------------*/
  std::vector<Teuchos::RCP<ConditionComponent>> particlewallcomponents;
  particlewallcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("MAT")));
  particlewallcomponents.push_back(Teuchos::rcp(new IntConditionComponent("MAT")));

  Teuchos::RCP<ConditionDefinition> surfpartwall =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE PARTICLE WALL", "ParticleWall",
          "Wall for particle interaction with (optional) material definition",
          DRT::Condition::ParticleWall, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < particlewallcomponents.size(); ++i)
    surfpartwall->AddComponent(particlewallcomponents[i]);

  condlist.push_back(surfpartwall);
}

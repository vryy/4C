/*---------------------------------------------------------------------------*/
/*! \file
\brief input parameters for particle problems

\level 3


*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include "baci_inpar_particle.hpp"

#include "baci_lib_conditiondefinition.hpp"
#include "baci_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | set the particle parameters                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void INPAR::PARTICLE::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace INPUT;
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
  CORE::UTILS::IntParameter(
      "RESULTSEVRY", 1, "write particle runtime output every RESULTSEVRY steps", &particledyn);
  CORE::UTILS::IntParameter(
      "RESTARTEVRY", 1, "write restart possibility every RESTARTEVRY steps", &particledyn);

  // write ghosted particles
  CORE::UTILS::BoolParameter(
      "WRITE_GHOSTED_PARTICLES", "no", "write ghosted particles (debug feature)", &particledyn);

  // time loop control
  CORE::UTILS::DoubleParameter("TIMESTEP", 0.01, "time step size", &particledyn);
  CORE::UTILS::IntParameter("NUMSTEP", 100, "maximum number of steps", &particledyn);
  CORE::UTILS::DoubleParameter("MAXTIME", 1.0, "maximum time", &particledyn);

  // gravity acceleration control
  CORE::UTILS::StringParameter(
      "GRAVITY_ACCELERATION", "0.0 0.0 0.0", "acceleration due to gravity", &particledyn);
  CORE::UTILS::IntParameter(
      "GRAVITY_RAMP_FUNCT", -1, "number of function governing gravity ramp", &particledyn);

  // viscous damping factor
  CORE::UTILS::DoubleParameter("VISCOUS_DAMPING", -1.0,
      "apply viscous damping force to determine static equilibrium solutions", &particledyn);

  // transfer particles to new bins every time step
  CORE::UTILS::BoolParameter(
      "TRANSFER_EVERY", "no", "transfer particles to new bins every time step", &particledyn);

  // considered particle phases with dynamic load balance weighting factor
  CORE::UTILS::StringParameter("PHASE_TO_DYNLOADBALFAC", "none",
      "considered particle phases with dynamic load balance weighting factor", &particledyn);

  // relate particle phase to material id
  CORE::UTILS::StringParameter(
      "PHASE_TO_MATERIAL_ID", "none", "relate particle phase to material id", &particledyn);

  // amplitude of noise added to initial position for each spatial direction
  CORE::UTILS::StringParameter("INITIAL_POSITION_AMPLITUDE", "0.0 0.0 0.0",
      "amplitude of noise added to initial position for each spatial direction", &particledyn);

  // type of particle wall source
  setStringToIntegralParameter<int>("PARTICLE_WALL_SOURCE", "NoParticleWall",
      "type of particle wall source",
      tuple<std::string>("NoParticleWall", "DiscretCondition", "BoundingBox"),
      tuple<int>(INPAR::PARTICLE::NoParticleWall, INPAR::PARTICLE::DiscretCondition,
          INPAR::PARTICLE::BoundingBox),
      &particledyn);

  // material id for particle wall from bounding box source
  CORE::UTILS::IntParameter("PARTICLE_WALL_MAT", -1,
      "material id for particle wall from bounding box source", &particledyn);

  // flags defining considered states of particle wall
  CORE::UTILS::BoolParameter(
      "PARTICLE_WALL_MOVING", "no", "consider a moving particle wall", &particledyn);
  CORE::UTILS::BoolParameter(
      "PARTICLE_WALL_LOADED", "no", "consider loading on particle wall", &particledyn);

  // consider rigid body motion
  CORE::UTILS::BoolParameter("RIGID_BODY_MOTION", "no", "consider rigid body motion", &particledyn);

  CORE::UTILS::DoubleParameter("RIGID_BODY_PHASECHANGE_RADIUS", -1.0,
      "search radius for neighboring rigid bodies in case of phase change", &particledyn);

  /*-------------------------------------------------------------------------*
   | control parameters for initial/boundary conditions                      |
   *-------------------------------------------------------------------------*/
  Teuchos::ParameterList& particledynconditions =
      particledyn.sublist("INITIAL AND BOUNDARY CONDITIONS", false,
          "control parameters for initial/boundary conditions in particle simulations\n");

  // initial temperature field of particle phase given by function
  CORE::UTILS::StringParameter("INITIAL_TEMP_FIELD", "none",
      "initial temperature field of particle phase given by function", &particledynconditions);

  // initial velocity field of particle phase given by function
  CORE::UTILS::StringParameter("INITIAL_VELOCITY_FIELD", "none",
      "initial velocity field of particle phase given by function", &particledynconditions);

  // initial acceleration field of particle phase given by function
  CORE::UTILS::StringParameter("INITIAL_ACCELERATION_FIELD", "none",
      "initial acceleration field of particle phase given by function", &particledynconditions);

  // dirichlet boundary condition of particle phase given by function
  CORE::UTILS::StringParameter("DIRICHLET_BOUNDARY_CONDITION", "none",
      "dirichlet boundary condition of particle phase given by function", &particledynconditions);

  // temperature boundary condition of particle phase given by function
  CORE::UTILS::StringParameter("TEMPERATURE_BOUNDARY_CONDITION", "none",
      "temperature boundary condition of particle phase given by function", &particledynconditions);

  /*-------------------------------------------------------------------------*
   | smoothed particle hydrodynamics (SPH) specific control parameters       |
   *-------------------------------------------------------------------------*/
  Teuchos::ParameterList& particledynsph = particledyn.sublist(
      "SPH", false, "control parameters for smoothed particle hydrodynamics (SPH) simulations\n");

  // write particle-wall interaction output
  CORE::UTILS::BoolParameter("WRITE_PARTICLE_WALL_INTERACTION", "no",
      "write particle-wall interaction output", &particledynsph);

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

  CORE::UTILS::DoubleParameter(
      "INITIALPARTICLESPACING", 0.0, "initial spacing of particles", &particledynsph);

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

  // type of wall formulation
  setStringToIntegralParameter<int>("WALLFORMULATION", "NoWallFormulation",
      "type of wall formulation",
      tuple<std::string>("NoWallFormulation", "VirtualParticleWallFormulation"),
      tuple<int>(
          INPAR::PARTICLE::NoWallFormulation, INPAR::PARTICLE::VirtualParticleWallFormulation),
      &particledynsph);

  // type of transport velocity formulation
  setStringToIntegralParameter<int>("TRANSPORTVELOCITYFORMULATION", "NoTransportVelocity",
      "type of transport velocity formulation",
      tuple<std::string>(
          "NoTransportVelocity", "StandardTransportVelocity", "GeneralizedTransportVelocity"),
      tuple<int>(INPAR::PARTICLE::NoTransportVelocity, INPAR::PARTICLE::StandardTransportVelocity,
          INPAR::PARTICLE::GeneralizedTransportVelocity),
      &particledynsph);

  // type of temperature evaluation scheme
  setStringToIntegralParameter<int>("TEMPERATUREEVALUATION", "NoTemperatureEvaluation",
      "type of temperature evaluation scheme",
      tuple<std::string>("NoTemperatureEvaluation", "TemperatureIntegration"),
      tuple<int>(INPAR::PARTICLE::NoTemperatureEvaluation, INPAR::PARTICLE::TemperatureIntegration),
      &particledynsph);

  CORE::UTILS::BoolParameter(
      "TEMPERATUREGRADIENT", "no", "evaluate temperature gradient", &particledynsph);

  // type of heat source
  setStringToIntegralParameter<int>("HEATSOURCETYPE", "NoHeatSource", "type of heat source",
      tuple<std::string>("NoHeatSource", "VolumeHeatSource", "SurfaceHeatSource"),
      tuple<int>(INPAR::PARTICLE::NoHeatSource, INPAR::PARTICLE::VolumeHeatSource,
          INPAR::PARTICLE::SurfaceHeatSource),
      &particledynsph);

  CORE::UTILS::IntParameter(
      "HEATSOURCE_FUNCT", -1, "number of function governing heat source", &particledynsph);

  CORE::UTILS::StringParameter(
      "HEATSOURCE_DIRECTION", "0.0 0.0 0.0", "direction of surface heat source", &particledynsph);

  // evaporation induced heat loss
  CORE::UTILS::BoolParameter(
      "VAPOR_HEATLOSS", "no", "evaluate evaporation induced heat loss", &particledynsph);
  CORE::UTILS::DoubleParameter(
      "VAPOR_HEATLOSS_LATENTHEAT", 0.0, "latent heat in heat loss formula", &particledynsph);
  CORE::UTILS::DoubleParameter("VAPOR_HEATLOSS_ENTHALPY_REFTEMP", 0.0,
      "enthalpy reference temperature in heat loss formula", &particledynsph);
  CORE::UTILS::DoubleParameter(
      "VAPOR_HEATLOSS_PFAC", 0.0, "pressure factor in heat loss formula", &particledynsph);
  CORE::UTILS::DoubleParameter(
      "VAPOR_HEATLOSS_TFAC", 0.0, "temperature factor in heat loss formula", &particledynsph);

  // evaporation induced recoil pressure
  CORE::UTILS::BoolParameter(
      "VAPOR_RECOIL", "no", "evaluate evaporation induced recoil pressure", &particledynsph);
  CORE::UTILS::DoubleParameter("VAPOR_RECOIL_BOILINGTEMPERATURE", 0.0,
      "boiling temperature in recoil pressure formula", &particledynsph);
  CORE::UTILS::DoubleParameter(
      "VAPOR_RECOIL_PFAC", 0.0, "pressure factor in recoil pressure formula", &particledynsph);
  CORE::UTILS::DoubleParameter(
      "VAPOR_RECOIL_TFAC", 0.0, "temperature factor in recoil pressure formula", &particledynsph);

  // type of surface tension formulation
  setStringToIntegralParameter<int>("SURFACETENSIONFORMULATION", "NoSurfaceTension",
      "type of surface tension formulation",
      tuple<std::string>("NoSurfaceTension", "ContinuumSurfaceForce"),
      tuple<int>(INPAR::PARTICLE::NoSurfaceTension, INPAR::PARTICLE::ContinuumSurfaceForce),
      &particledynsph);

  CORE::UTILS::IntParameter("SURFACETENSION_RAMP_FUNCT", -1,
      "number of function governing surface tension ramp", &particledynsph);

  CORE::UTILS::DoubleParameter("SURFACETENSIONCOEFFICIENT", -1.0,
      "constant part of surface tension coefficient", &particledynsph);
  CORE::UTILS::DoubleParameter("SURFACETENSIONMINIMUM", 0.0,
      "minimum surface tension coefficient in case of temperature dependence", &particledynsph);
  CORE::UTILS::DoubleParameter("SURFACETENSIONTEMPFAC", 0.0,
      "factor of dependence of surface tension coefficient on temperature", &particledynsph);
  CORE::UTILS::DoubleParameter("SURFACETENSIONREFTEMP", 0.0,
      "reference temperature for surface tension coefficient", &particledynsph);

  // wetting
  CORE::UTILS::DoubleParameter("STATICCONTACTANGLE", 0.0,
      "static contact angle in degree with wetting effects", &particledynsph);
  CORE::UTILS::DoubleParameter("TRIPLEPOINTNORMAL_CORR_CF_LOW", 0.0,
      "triple point normal correction wall color field low", &particledynsph);
  CORE::UTILS::DoubleParameter("TRIPLEPOINTNORMAL_CORR_CF_UP", 0.0,
      "triple point normal correction wall color field up", &particledynsph);

  // interface viscosity
  CORE::UTILS::BoolParameter(
      "INTERFACE_VISCOSITY", "no", "evaluate interface viscosity", &particledynsph);
  CORE::UTILS::DoubleParameter("INTERFACE_VISCOSITY_LIQUIDGAS", 0.0,
      "artificial viscosity on liquid-gas interface", &particledynsph);
  CORE::UTILS::DoubleParameter("INTERFACE_VISCOSITY_SOLIDLIQUID", 0.0,
      "artificial viscosity on solid-liquid interface", &particledynsph);

  // barrier force
  CORE::UTILS::BoolParameter("BARRIER_FORCE", "no", "evaluate barrier force", &particledynsph);
  CORE::UTILS::DoubleParameter(
      "BARRIER_FORCE_DISTANCE", 0.0, "barrier force distance", &particledynsph);
  CORE::UTILS::DoubleParameter(
      "BARRIER_FORCE_TEMPSCALE", 0.0, "barrier force temperature scaling", &particledynsph);
  CORE::UTILS::DoubleParameter(
      "BARRIER_FORCE_STIFF_HEAVY", -1.0, "barrier force stiffness of heavy phase", &particledynsph);
  CORE::UTILS::DoubleParameter("BARRIER_FORCE_DAMP_HEAVY", 0.0,
      "barrier force damping parameter of heavy phase", &particledynsph);
  CORE::UTILS::DoubleParameter(
      "BARRIER_FORCE_STIFF_GAS", -1.0, "barrier force stiffness of gas phase", &particledynsph);
  CORE::UTILS::DoubleParameter("BARRIER_FORCE_DAMP_GAS", 0.0,
      "barrier force damping parameter of gas phase", &particledynsph);

  // linear transition in surface tension evaluation
  CORE::UTILS::DoubleParameter(
      "TRANS_REF_TEMPERATURE", 0.0, "transition reference temperature", &particledynsph);
  CORE::UTILS::DoubleParameter("TRANS_DT_SURFACETENSION", 0.0,
      "transition temperature difference for surface tension evaluation", &particledynsph);
  CORE::UTILS::DoubleParameter("TRANS_DT_MARANGONI", 0.0,
      "transition temperature difference for marangoni evaluation", &particledynsph);
  CORE::UTILS::DoubleParameter("TRANS_DT_CURVATURE", 0.0,
      "transition temperature difference for curvature evaluation", &particledynsph);
  CORE::UTILS::DoubleParameter("TRANS_DT_WETTING", 0.0,
      "transition temperature difference for wetting evaluation", &particledynsph);
  CORE::UTILS::DoubleParameter("TRANS_DT_INTVISC", 0.0,
      "transition temperature difference for interface viscosity evaluation", &particledynsph);
  CORE::UTILS::DoubleParameter("TRANS_DT_BARRIER", 0.0,
      "transition temperature difference for barrier force evaluation", &particledynsph);

  // type of dirichlet open boundary
  setStringToIntegralParameter<int>("DIRICHLETBOUNDARYTYPE", "NoDirichletOpenBoundary",
      "type of dirichlet open boundary",
      tuple<std::string>("NoDirichletOpenBoundary", "DirichletNormalToPlane"),
      tuple<int>(INPAR::PARTICLE::NoDirichletOpenBoundary, INPAR::PARTICLE::DirichletNormalToPlane),
      &particledynsph);

  CORE::UTILS::IntParameter("DIRICHLET_FUNCT", -1,
      "number of function governing velocity condition on dirichlet open boundary",
      &particledynsph);

  CORE::UTILS::StringParameter("DIRICHLET_OUTWARD_NORMAL", "0.0 0.0 0.0",
      "direction of outward normal on dirichlet open boundary", &particledynsph);
  CORE::UTILS::StringParameter("DIRICHLET_PLANE_POINT", "0.0 0.0 0.0",
      "point on dirichlet open boundary plane", &particledynsph);

  // type of neumann open boundary
  setStringToIntegralParameter<int>("NEUMANNBOUNDARYTYPE", "NoNeumannOpenBoundary",
      "type of neumann open boundary",
      tuple<std::string>("NoNeumannOpenBoundary", "NeumannNormalToPlane"),
      tuple<int>(INPAR::PARTICLE::NoNeumannOpenBoundary, INPAR::PARTICLE::NeumannNormalToPlane),
      &particledynsph);

  CORE::UTILS::IntParameter("NEUMANN_FUNCT", -1,
      "number of function governing pressure condition on neumann open boundary", &particledynsph);

  CORE::UTILS::StringParameter("NEUMANN_OUTWARD_NORMAL", "0.0 0.0 0.0",
      "direction of outward normal on neumann open boundary", &particledynsph);
  CORE::UTILS::StringParameter("NEUMANN_PLANE_POINT", "0.0 0.0 0.0",
      "point on neumann open boundary plane", &particledynsph);

  // type of phase change
  setStringToIntegralParameter<int>("PHASECHANGETYPE", "NoPhaseChange", "type of phase change",
      tuple<std::string>("NoPhaseChange", "OneWayScalarBelowToAbovePhaseChange",
          "OneWayScalarAboveToBelowPhaseChange", "TwoWayScalarPhaseChange"),
      tuple<int>(INPAR::PARTICLE::NoPhaseChange,
          INPAR::PARTICLE::OneWayScalarBelowToAbovePhaseChange,
          INPAR::PARTICLE::OneWayScalarAboveToBelowPhaseChange,
          INPAR::PARTICLE::TwoWayScalarPhaseChange),
      &particledynsph);

  // definition of phase change
  CORE::UTILS::StringParameter(
      "PHASECHANGEDEFINITION", "none", "phase change definition", &particledynsph);

  // type of rigid particle contact
  setStringToIntegralParameter<int>("RIGIDPARTICLECONTACTTYPE", "NoRigidParticleContact",
      "type of rigid particle contact",
      tuple<std::string>("NoRigidParticleContact", "ElasticRigidParticleContact"),
      tuple<int>(
          INPAR::PARTICLE::NoRigidParticleContact, INPAR::PARTICLE::ElasticRigidParticleContact),
      &particledynsph);

  CORE::UTILS::DoubleParameter(
      "RIGIDPARTICLECONTACTSTIFF", -1.0, "rigid particle contact stiffness", &particledynsph);
  CORE::UTILS::DoubleParameter(
      "RIGIDPARTICLECONTACTDAMP", 0.0, "rigid particle contact damping parameter", &particledynsph);

  /*-------------------------------------------------------------------------*
   | discrete element method (DEM) specific control parameters               |
   *-------------------------------------------------------------------------*/
  Teuchos::ParameterList& particledyndem = particledyn.sublist(
      "DEM", false, "control parameters for discrete element method (DEM) simulations\n");

  // write particle energy output
  CORE::UTILS::BoolParameter(
      "WRITE_PARTICLE_ENERGY", "no", "write particle energy output", &particledyndem);

  // write particle-wall interaction output
  CORE::UTILS::BoolParameter("WRITE_PARTICLE_WALL_INTERACTION", "no",
      "write particle-wall interaction output", &particledyndem);

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

  CORE::UTILS::DoubleParameter(
      "MIN_RADIUS", 0.0, "minimum allowed particle radius", &particledyndem);
  CORE::UTILS::DoubleParameter(
      "MAX_RADIUS", 0.0, "maximum allowed particle radius", &particledyndem);
  CORE::UTILS::DoubleParameter(
      "MAX_VELOCITY", -1.0, "maximum expected particle velocity", &particledyndem);

  // type of initial particle radius assignment
  setStringToIntegralParameter<int>("INITIAL_RADIUS", "RadiusFromParticleMaterial",
      "type of initial particle radius assignment",
      tuple<std::string>("RadiusFromParticleMaterial", "RadiusFromParticleInput",
          "NormalRadiusDistribution", "LogNormalRadiusDistribution"),
      tuple<int>(INPAR::PARTICLE::RadiusFromParticleMaterial,
          INPAR::PARTICLE::RadiusFromParticleInput, INPAR::PARTICLE::NormalRadiusDistribution,
          INPAR::PARTICLE::LogNormalRadiusDistribution),
      &particledyndem);

  CORE::UTILS::DoubleParameter("RADIUSDISTRIBUTION_SIGMA", -1.0,
      "sigma of random particle radius distribution", &particledyndem);

  CORE::UTILS::DoubleParameter(
      "REL_PENETRATION", -1.0, "maximum allowed relative penetration", &particledyndem);
  CORE::UTILS::DoubleParameter("NORMAL_STIFF", -1.0, "normal contact stiffness", &particledyndem);
  CORE::UTILS::DoubleParameter(
      "NORMAL_DAMP", -1.0, "normal contact damping parameter", &particledyndem);
  CORE::UTILS::DoubleParameter(
      "COEFF_RESTITUTION", -1.0, "coefficient of restitution", &particledyndem);
  CORE::UTILS::DoubleParameter("DAMP_REG_FAC", -1.0,
      "linearly regularized damping normal force in the interval "
      "\f$|g| < (\\text{DAMP_REG_FAC} \\cdot r_{\\min})\f$",
      &particledyndem);
  CORE::UTILS::BoolParameter(
      "TENSION_CUTOFF", "yes", "evaluate tension cutoff of normal contact force", &particledyndem);

  CORE::UTILS::DoubleParameter("POISSON_RATIO", -1.0, "poisson ratio", &particledyndem);
  CORE::UTILS::DoubleParameter("YOUNG_MODULUS", -1.0, "young's modulus", &particledyndem);

  CORE::UTILS::DoubleParameter(
      "FRICT_COEFF_TANG", -1.0, "friction coefficient for tangential contact", &particledyndem);
  CORE::UTILS::DoubleParameter(
      "FRICT_COEFF_ROLL", -1.0, "friction coefficient for rolling contact", &particledyndem);

  CORE::UTILS::DoubleParameter(
      "ADHESION_DISTANCE", -1.0, "adhesion distance between interacting surfaces", &particledyndem);

  CORE::UTILS::DoubleParameter(
      "ADHESION_MAX_CONTACT_PRESSURE", 0.0, "adhesion maximum contact pressure", &particledyndem);
  CORE::UTILS::DoubleParameter(
      "ADHESION_MAX_CONTACT_FORCE", 0.0, "adhesion maximum contact force", &particledyndem);
  CORE::UTILS::BoolParameter("ADHESION_USE_MAX_CONTACT_FORCE", "no",
      "use maximum contact force instead of maximum contact pressure", &particledyndem);

  CORE::UTILS::BoolParameter(
      "ADHESION_VDW_CURVE_SHIFT", "no", "shifts van-der-Waals-curve to g = 0", &particledyndem);

  CORE::UTILS::DoubleParameter(
      "ADHESION_HAMAKER", -1.0, "hamaker constant of van-der-Waals interaction", &particledyndem);

  CORE::UTILS::DoubleParameter("ADHESION_SURFACE_ENERGY", -1.0,
      "adhesion surface energy for the calculation of the pull-out force", &particledyndem);
  CORE::UTILS::DoubleParameter("ADHESION_SURFACE_ENERGY_DISTRIBUTION_VAR", -1.0,
      "variance of adhesion surface energy distribution", &particledyndem);
  CORE::UTILS::DoubleParameter("ADHESION_SURFACE_ENERGY_DISTRIBUTION_CUTOFF_FACTOR", -1.0,
      "adhesion surface energy distribution limited by multiple of variance", &particledyndem);
  CORE::UTILS::DoubleParameter("ADHESION_SURFACE_ENERGY_FACTOR", 1.0,
      "factor to calculate minimum adhesion surface energy", &particledyndem);
}

/*---------------------------------------------------------------------------*
 | set the particle conditions                                sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
void INPAR::PARTICLE::SetValidConditions(
    std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist)
{
  using namespace INPUT;

  /*-------------------------------------------------------------------------*
   | particle wall condition                                                 |
   *-------------------------------------------------------------------------*/
  std::vector<Teuchos::RCP<INPUT::LineComponent>> particlewallcomponents;
  particlewallcomponents.push_back(Teuchos::rcp(new INPUT::SeparatorComponent("MAT")));
  particlewallcomponents.push_back(Teuchos::rcp(new INPUT::IntComponent("MAT")));

  Teuchos::RCP<ConditionDefinition> surfpartwall =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE PARTICLE WALL", "ParticleWall",
          "Wall for particle interaction with (optional) material definition",
          DRT::Condition::ParticleWall, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < particlewallcomponents.size(); ++i)
    surfpartwall->AddComponent(particlewallcomponents[i]);

  condlist.push_back(surfpartwall);
}

FOUR_C_NAMESPACE_CLOSE

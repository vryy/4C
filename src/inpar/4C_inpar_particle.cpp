#include "4C_inpar_particle.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | set the particle parameters                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
void Inpar::PARTICLE::set_valid_parameters(Teuchos::ParameterList& list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  /*-------------------------------------------------------------------------*
   | general control parameters for particle simulations                     |
   *-------------------------------------------------------------------------*/
  Teuchos::ParameterList& particledyn =
      list.sublist("PARTICLE DYNAMIC", false, "control parameters for particle simulations\n");

  // type of particle time integration
  setStringToIntegralParameter<DynamicType>("DYNAMICTYP", "VelocityVerlet",
      "type of particle time integration",
      tuple<std::string>("SemiImplicitEuler", "VelocityVerlet"),
      tuple<DynamicType>(
          Inpar::PARTICLE::dyna_semiimpliciteuler, Inpar::PARTICLE::dyna_velocityverlet),
      &particledyn);

  // type of particle interaction
  setStringToIntegralParameter<InteractionType>("INTERACTION", "None",
      "type of particle interaction", tuple<std::string>("None", "SPH", "DEM"),
      tuple<InteractionType>(Inpar::PARTICLE::interaction_none, Inpar::PARTICLE::interaction_sph,
          Inpar::PARTICLE::interaction_dem),
      &particledyn);

  // output type
  Core::Utils::int_parameter(
      "RESULTSEVRY", 1, "write particle runtime output every RESULTSEVRY steps", &particledyn);
  Core::Utils::int_parameter(
      "RESTARTEVRY", 1, "write restart possibility every RESTARTEVRY steps", &particledyn);

  // write ghosted particles
  Core::Utils::bool_parameter(
      "WRITE_GHOSTED_PARTICLES", "no", "write ghosted particles (debug feature)", &particledyn);

  // time loop control
  Core::Utils::double_parameter("TIMESTEP", 0.01, "time step size", &particledyn);
  Core::Utils::int_parameter("NUMSTEP", 100, "maximum number of steps", &particledyn);
  Core::Utils::double_parameter("MAXTIME", 1.0, "maximum time", &particledyn);

  // gravity acceleration control
  Core::Utils::string_parameter(
      "GRAVITY_ACCELERATION", "0.0 0.0 0.0", "acceleration due to gravity", &particledyn);
  Core::Utils::int_parameter(
      "GRAVITY_RAMP_FUNCT", -1, "number of function governing gravity ramp", &particledyn);

  // viscous damping factor
  Core::Utils::double_parameter("VISCOUS_DAMPING", -1.0,
      "apply viscous damping force to determine static equilibrium solutions", &particledyn);

  // transfer particles to new bins every time step
  Core::Utils::bool_parameter(
      "TRANSFER_EVERY", "no", "transfer particles to new bins every time step", &particledyn);

  // considered particle phases with dynamic load balance weighting factor
  Core::Utils::string_parameter("PHASE_TO_DYNLOADBALFAC", "none",
      "considered particle phases with dynamic load balance weighting factor", &particledyn);

  // relate particle phase to material id
  Core::Utils::string_parameter(
      "PHASE_TO_MATERIAL_ID", "none", "relate particle phase to material id", &particledyn);

  // amplitude of noise added to initial position for each spatial direction
  Core::Utils::string_parameter("INITIAL_POSITION_AMPLITUDE", "0.0 0.0 0.0",
      "amplitude of noise added to initial position for each spatial direction", &particledyn);

  // type of particle wall source
  setStringToIntegralParameter<ParticleWallSource>("PARTICLE_WALL_SOURCE", "NoParticleWall",
      "type of particle wall source",
      tuple<std::string>("NoParticleWall", "DiscretCondition", "BoundingBox"),
      tuple<ParticleWallSource>(Inpar::PARTICLE::NoParticleWall, Inpar::PARTICLE::DiscretCondition,
          Inpar::PARTICLE::BoundingBox),
      &particledyn);

  // material id for particle wall from bounding box source
  Core::Utils::int_parameter("PARTICLE_WALL_MAT", -1,
      "material id for particle wall from bounding box source", &particledyn);

  // flags defining considered states of particle wall
  Core::Utils::bool_parameter(
      "PARTICLE_WALL_MOVING", "no", "consider a moving particle wall", &particledyn);
  Core::Utils::bool_parameter(
      "PARTICLE_WALL_LOADED", "no", "consider loading on particle wall", &particledyn);

  // consider rigid body motion
  Core::Utils::bool_parameter(
      "RIGID_BODY_MOTION", "no", "consider rigid body motion", &particledyn);

  Core::Utils::double_parameter("RIGID_BODY_PHASECHANGE_RADIUS", -1.0,
      "search radius for neighboring rigid bodies in case of phase change", &particledyn);

  /*-------------------------------------------------------------------------*
   | control parameters for initial/boundary conditions                      |
   *-------------------------------------------------------------------------*/
  Teuchos::ParameterList& particledynconditions =
      particledyn.sublist("INITIAL AND BOUNDARY CONDITIONS", false,
          "control parameters for initial/boundary conditions in particle simulations\n");

  // initial temperature field of particle phase given by function
  Core::Utils::string_parameter("INITIAL_TEMP_FIELD", "none",
      "Refer to the function ID describing the initial temperature field of particle phase",
      &particledynconditions);

  // initial velocity field of particle phase given by function
  Core::Utils::string_parameter("INITIAL_VELOCITY_FIELD", "none",
      "Refer to the function ID describing the initial velocity field of particle phase",
      &particledynconditions);

  // initial angular velocity field of particle phase given by function
  Core::Utils::string_parameter("INITIAL_ANGULAR_VELOCITY_FIELD", "none",
      "Refer to the function ID describing the initial angular velocity field of rigid body "
      "phase/DEM particle",
      &particledynconditions);

  // initial acceleration field of particle phase given by function
  Core::Utils::string_parameter("INITIAL_ACCELERATION_FIELD", "none",
      "Refer to the function ID describing the initial acceleration field of particle phase",
      &particledynconditions);

  // initial angular acceleration field of particle phase given by function
  Core::Utils::string_parameter("INITIAL_ANGULAR_ACCELERATION_FIELD", "none",
      "Refer to the function ID describing the initial angular acceleration field of rigid body "
      "phase/DEM particle",
      &particledynconditions);

  // dirichlet boundary condition of particle phase given by function
  Core::Utils::string_parameter("DIRICHLET_BOUNDARY_CONDITION", "none",
      "Refer to the function ID describing the dirichlet boundary condition of particle phase",
      &particledynconditions);

  // temperature boundary condition of particle phase given by function
  Core::Utils::string_parameter("TEMPERATURE_BOUNDARY_CONDITION", "none",
      "Refer to the function ID describing the temperature boundary condition of particle phase",
      &particledynconditions);

  /*-------------------------------------------------------------------------*
   | smoothed particle hydrodynamics (SPH) specific control parameters       |
   *-------------------------------------------------------------------------*/
  Teuchos::ParameterList& particledynsph = particledyn.sublist(
      "SPH", false, "control parameters for smoothed particle hydrodynamics (SPH) simulations\n");

  // write particle-wall interaction output
  Core::Utils::bool_parameter("WRITE_PARTICLE_WALL_INTERACTION", "no",
      "write particle-wall interaction output", &particledynsph);

  // type of smoothed particle hydrodynamics kernel
  setStringToIntegralParameter<KernelType>("KERNEL", "CubicSpline",
      "type of smoothed particle hydrodynamics kernel",
      tuple<std::string>("CubicSpline", "QuinticSpline"),
      tuple<KernelType>(Inpar::PARTICLE::CubicSpline, Inpar::PARTICLE::QuinticSpline),
      &particledynsph);

  // kernel space dimension number
  setStringToIntegralParameter<KernelSpaceDimension>("KERNEL_SPACE_DIM", "Kernel3D",
      "kernel space dimension number", tuple<std::string>("Kernel1D", "Kernel2D", "Kernel3D"),
      tuple<KernelSpaceDimension>(
          Inpar::PARTICLE::Kernel1D, Inpar::PARTICLE::Kernel2D, Inpar::PARTICLE::Kernel3D),
      &particledynsph);

  Core::Utils::double_parameter(
      "INITIALPARTICLESPACING", 0.0, "initial spacing of particles", &particledynsph);

  // type of smoothed particle hydrodynamics equation of state
  setStringToIntegralParameter<EquationOfStateType>("EQUATIONOFSTATE", "GenTait",
      "type of smoothed particle hydrodynamics equation of state",
      tuple<std::string>("GenTait", "IdealGas"),
      tuple<EquationOfStateType>(Inpar::PARTICLE::GenTait, Inpar::PARTICLE::IdealGas),
      &particledynsph);

  // type of smoothed particle hydrodynamics momentum formulation
  setStringToIntegralParameter<MomentumFormulationType>("MOMENTUMFORMULATION",
      "AdamiMomentumFormulation", "type of smoothed particle hydrodynamics momentum formulation",
      tuple<std::string>("AdamiMomentumFormulation", "MonaghanMomentumFormulation"),
      tuple<MomentumFormulationType>(
          Inpar::PARTICLE::AdamiMomentumFormulation, Inpar::PARTICLE::MonaghanMomentumFormulation),
      &particledynsph);

  // type of density evaluation scheme
  setStringToIntegralParameter<DensityEvaluationScheme>("DENSITYEVALUATION", "DensitySummation",
      "type of density evaluation scheme",
      tuple<std::string>("DensitySummation", "DensityIntegration", "DensityPredictCorrect"),
      tuple<DensityEvaluationScheme>(Inpar::PARTICLE::DensitySummation,
          Inpar::PARTICLE::DensityIntegration, Inpar::PARTICLE::DensityPredictCorrect),
      &particledynsph);

  // type of density correction scheme
  setStringToIntegralParameter<DensityCorrectionScheme>("DENSITYCORRECTION", "NoCorrection",
      "type of density correction scheme",
      tuple<std::string>(
          "NoCorrection", "InteriorCorrection", "NormalizedCorrection", "RandlesCorrection"),
      tuple<DensityCorrectionScheme>(Inpar::PARTICLE::NoCorrection,
          Inpar::PARTICLE::InteriorCorrection, Inpar::PARTICLE::NormalizedCorrection,
          Inpar::PARTICLE::RandlesCorrection),
      &particledynsph);

  // type of boundary particle formulation
  setStringToIntegralParameter<BoundaryParticleFormulationType>("BOUNDARYPARTICLEFORMULATION",
      "NoBoundaryFormulation", "type of boundary particle formulation",
      tuple<std::string>("NoBoundaryFormulation", "AdamiBoundaryFormulation"),
      tuple<BoundaryParticleFormulationType>(
          Inpar::PARTICLE::NoBoundaryFormulation, Inpar::PARTICLE::AdamiBoundaryFormulation),
      &particledynsph);

  // type of boundary particle interaction
  setStringToIntegralParameter<BoundaryParticleInteraction>("BOUNDARYPARTICLEINTERACTION",
      "NoSlipBoundaryParticle", "type of boundary particle interaction",
      tuple<std::string>("NoSlipBoundaryParticle", "FreeSlipBoundaryParticle"),
      tuple<BoundaryParticleInteraction>(
          Inpar::PARTICLE::NoSlipBoundaryParticle, Inpar::PARTICLE::FreeSlipBoundaryParticle),
      &particledynsph);

  // type of wall formulation
  setStringToIntegralParameter<WallFormulationType>("WALLFORMULATION", "NoWallFormulation",
      "type of wall formulation",
      tuple<std::string>("NoWallFormulation", "VirtualParticleWallFormulation"),
      tuple<WallFormulationType>(
          Inpar::PARTICLE::NoWallFormulation, Inpar::PARTICLE::VirtualParticleWallFormulation),
      &particledynsph);

  // type of transport velocity formulation
  setStringToIntegralParameter<TransportVelocityFormulation>("TRANSPORTVELOCITYFORMULATION",
      "NoTransportVelocity", "type of transport velocity formulation",
      tuple<std::string>(
          "NoTransportVelocity", "StandardTransportVelocity", "GeneralizedTransportVelocity"),
      tuple<TransportVelocityFormulation>(Inpar::PARTICLE::NoTransportVelocity,
          Inpar::PARTICLE::StandardTransportVelocity,
          Inpar::PARTICLE::GeneralizedTransportVelocity),
      &particledynsph);

  // type of temperature evaluation scheme
  setStringToIntegralParameter<TemperatureEvaluationScheme>("TEMPERATUREEVALUATION",
      "NoTemperatureEvaluation", "type of temperature evaluation scheme",
      tuple<std::string>("NoTemperatureEvaluation", "TemperatureIntegration"),
      tuple<TemperatureEvaluationScheme>(
          Inpar::PARTICLE::NoTemperatureEvaluation, Inpar::PARTICLE::TemperatureIntegration),
      &particledynsph);

  Core::Utils::bool_parameter(
      "TEMPERATUREGRADIENT", "no", "evaluate temperature gradient", &particledynsph);

  // type of heat source
  setStringToIntegralParameter<HeatSourceType>("HEATSOURCETYPE", "NoHeatSource",
      "type of heat source",
      tuple<std::string>("NoHeatSource", "VolumeHeatSource", "SurfaceHeatSource"),
      tuple<HeatSourceType>(Inpar::PARTICLE::NoHeatSource, Inpar::PARTICLE::VolumeHeatSource,
          Inpar::PARTICLE::SurfaceHeatSource),
      &particledynsph);

  Core::Utils::int_parameter(
      "HEATSOURCE_FUNCT", -1, "number of function governing heat source", &particledynsph);

  Core::Utils::string_parameter(
      "HEATSOURCE_DIRECTION", "0.0 0.0 0.0", "direction of surface heat source", &particledynsph);

  // evaporation induced heat loss
  Core::Utils::bool_parameter(
      "VAPOR_HEATLOSS", "no", "evaluate evaporation induced heat loss", &particledynsph);
  Core::Utils::double_parameter(
      "VAPOR_HEATLOSS_LATENTHEAT", 0.0, "latent heat in heat loss formula", &particledynsph);
  Core::Utils::double_parameter("VAPOR_HEATLOSS_ENTHALPY_REFTEMP", 0.0,
      "enthalpy reference temperature in heat loss formula", &particledynsph);
  Core::Utils::double_parameter(
      "VAPOR_HEATLOSS_PFAC", 0.0, "pressure factor in heat loss formula", &particledynsph);
  Core::Utils::double_parameter(
      "VAPOR_HEATLOSS_TFAC", 0.0, "temperature factor in heat loss formula", &particledynsph);

  // evaporation induced recoil pressure
  Core::Utils::bool_parameter(
      "VAPOR_RECOIL", "no", "evaluate evaporation induced recoil pressure", &particledynsph);
  Core::Utils::double_parameter("VAPOR_RECOIL_BOILINGTEMPERATURE", 0.0,
      "boiling temperature in recoil pressure formula", &particledynsph);
  Core::Utils::double_parameter(
      "VAPOR_RECOIL_PFAC", 0.0, "pressure factor in recoil pressure formula", &particledynsph);
  Core::Utils::double_parameter(
      "VAPOR_RECOIL_TFAC", 0.0, "temperature factor in recoil pressure formula", &particledynsph);

  // type of surface tension formulation
  setStringToIntegralParameter<SurfaceTensionFormulation>("SURFACETENSIONFORMULATION",
      "NoSurfaceTension", "type of surface tension formulation",
      tuple<std::string>("NoSurfaceTension", "ContinuumSurfaceForce"),
      tuple<SurfaceTensionFormulation>(
          Inpar::PARTICLE::NoSurfaceTension, Inpar::PARTICLE::ContinuumSurfaceForce),
      &particledynsph);

  Core::Utils::int_parameter("SURFACETENSION_RAMP_FUNCT", -1,
      "number of function governing surface tension ramp", &particledynsph);

  Core::Utils::double_parameter("SURFACETENSIONCOEFFICIENT", -1.0,
      "constant part of surface tension coefficient", &particledynsph);
  Core::Utils::double_parameter("SURFACETENSIONMINIMUM", 0.0,
      "minimum surface tension coefficient in case of temperature dependence", &particledynsph);
  Core::Utils::double_parameter("SURFACETENSIONTEMPFAC", 0.0,
      "factor of dependence of surface tension coefficient on temperature", &particledynsph);
  Core::Utils::double_parameter("SURFACETENSIONREFTEMP", 0.0,
      "reference temperature for surface tension coefficient", &particledynsph);

  // wetting
  Core::Utils::double_parameter("STATICCONTACTANGLE", 0.0,
      "static contact angle in degree with wetting effects", &particledynsph);
  Core::Utils::double_parameter("TRIPLEPOINTNORMAL_CORR_CF_LOW", 0.0,
      "triple point normal correction wall color field low", &particledynsph);
  Core::Utils::double_parameter("TRIPLEPOINTNORMAL_CORR_CF_UP", 0.0,
      "triple point normal correction wall color field up", &particledynsph);

  // interface viscosity
  Core::Utils::bool_parameter(
      "INTERFACE_VISCOSITY", "no", "evaluate interface viscosity", &particledynsph);
  Core::Utils::double_parameter("INTERFACE_VISCOSITY_LIQUIDGAS", 0.0,
      "artificial viscosity on liquid-gas interface", &particledynsph);
  Core::Utils::double_parameter("INTERFACE_VISCOSITY_SOLIDLIQUID", 0.0,
      "artificial viscosity on solid-liquid interface", &particledynsph);

  // barrier force
  Core::Utils::bool_parameter("BARRIER_FORCE", "no", "evaluate barrier force", &particledynsph);
  Core::Utils::double_parameter(
      "BARRIER_FORCE_DISTANCE", 0.0, "barrier force distance", &particledynsph);
  Core::Utils::double_parameter(
      "BARRIER_FORCE_TEMPSCALE", 0.0, "barrier force temperature scaling", &particledynsph);
  Core::Utils::double_parameter(
      "BARRIER_FORCE_STIFF_HEAVY", -1.0, "barrier force stiffness of heavy phase", &particledynsph);
  Core::Utils::double_parameter("BARRIER_FORCE_DAMP_HEAVY", 0.0,
      "barrier force damping parameter of heavy phase", &particledynsph);
  Core::Utils::double_parameter(
      "BARRIER_FORCE_STIFF_GAS", -1.0, "barrier force stiffness of gas phase", &particledynsph);
  Core::Utils::double_parameter("BARRIER_FORCE_DAMP_GAS", 0.0,
      "barrier force damping parameter of gas phase", &particledynsph);

  // linear transition in surface tension evaluation
  Core::Utils::double_parameter(
      "TRANS_REF_TEMPERATURE", 0.0, "transition reference temperature", &particledynsph);
  Core::Utils::double_parameter("TRANS_DT_SURFACETENSION", 0.0,
      "transition temperature difference for surface tension evaluation", &particledynsph);
  Core::Utils::double_parameter("TRANS_DT_MARANGONI", 0.0,
      "transition temperature difference for marangoni evaluation", &particledynsph);
  Core::Utils::double_parameter("TRANS_DT_CURVATURE", 0.0,
      "transition temperature difference for curvature evaluation", &particledynsph);
  Core::Utils::double_parameter("TRANS_DT_WETTING", 0.0,
      "transition temperature difference for wetting evaluation", &particledynsph);
  Core::Utils::double_parameter("TRANS_DT_INTVISC", 0.0,
      "transition temperature difference for interface viscosity evaluation", &particledynsph);
  Core::Utils::double_parameter("TRANS_DT_BARRIER", 0.0,
      "transition temperature difference for barrier force evaluation", &particledynsph);

  // type of dirichlet open boundary
  setStringToIntegralParameter<DirichletOpenBoundaryType>("DIRICHLETBOUNDARYTYPE",
      "NoDirichletOpenBoundary", "type of dirichlet open boundary",
      tuple<std::string>("NoDirichletOpenBoundary", "DirichletNormalToPlane"),
      tuple<DirichletOpenBoundaryType>(
          Inpar::PARTICLE::NoDirichletOpenBoundary, Inpar::PARTICLE::DirichletNormalToPlane),
      &particledynsph);

  Core::Utils::int_parameter("DIRICHLET_FUNCT", -1,
      "number of function governing velocity condition on dirichlet open boundary",
      &particledynsph);

  Core::Utils::string_parameter("DIRICHLET_OUTWARD_NORMAL", "0.0 0.0 0.0",
      "direction of outward normal on dirichlet open boundary", &particledynsph);
  Core::Utils::string_parameter("DIRICHLET_PLANE_POINT", "0.0 0.0 0.0",
      "point on dirichlet open boundary plane", &particledynsph);

  // type of neumann open boundary
  setStringToIntegralParameter<NeumannOpenBoundaryType>("NEUMANNBOUNDARYTYPE",
      "NoNeumannOpenBoundary", "type of neumann open boundary",
      tuple<std::string>("NoNeumannOpenBoundary", "NeumannNormalToPlane"),
      tuple<NeumannOpenBoundaryType>(
          Inpar::PARTICLE::NoNeumannOpenBoundary, Inpar::PARTICLE::NeumannNormalToPlane),
      &particledynsph);

  Core::Utils::int_parameter("NEUMANN_FUNCT", -1,
      "number of function governing pressure condition on neumann open boundary", &particledynsph);

  Core::Utils::string_parameter("NEUMANN_OUTWARD_NORMAL", "0.0 0.0 0.0",
      "direction of outward normal on neumann open boundary", &particledynsph);
  Core::Utils::string_parameter("NEUMANN_PLANE_POINT", "0.0 0.0 0.0",
      "point on neumann open boundary plane", &particledynsph);

  // type of phase change
  setStringToIntegralParameter<PhaseChangeType>("PHASECHANGETYPE", "NoPhaseChange",
      "type of phase change",
      tuple<std::string>("NoPhaseChange", "OneWayScalarBelowToAbovePhaseChange",
          "OneWayScalarAboveToBelowPhaseChange", "TwoWayScalarPhaseChange"),
      tuple<PhaseChangeType>(Inpar::PARTICLE::NoPhaseChange,
          Inpar::PARTICLE::OneWayScalarBelowToAbovePhaseChange,
          Inpar::PARTICLE::OneWayScalarAboveToBelowPhaseChange,
          Inpar::PARTICLE::TwoWayScalarPhaseChange),
      &particledynsph);

  // definition of phase change
  Core::Utils::string_parameter(
      "PHASECHANGEDEFINITION", "none", "phase change definition", &particledynsph);

  // type of rigid particle contact
  setStringToIntegralParameter<RigidParticleContactType>("RIGIDPARTICLECONTACTTYPE",
      "NoRigidParticleContact", "type of rigid particle contact",
      tuple<std::string>("NoRigidParticleContact", "ElasticRigidParticleContact"),
      tuple<RigidParticleContactType>(
          Inpar::PARTICLE::NoRigidParticleContact, Inpar::PARTICLE::ElasticRigidParticleContact),
      &particledynsph);

  Core::Utils::double_parameter(
      "RIGIDPARTICLECONTACTSTIFF", -1.0, "rigid particle contact stiffness", &particledynsph);
  Core::Utils::double_parameter(
      "RIGIDPARTICLECONTACTDAMP", 0.0, "rigid particle contact damping parameter", &particledynsph);

  /*-------------------------------------------------------------------------*
   | discrete element method (DEM) specific control parameters               |
   *-------------------------------------------------------------------------*/
  Teuchos::ParameterList& particledyndem = particledyn.sublist(
      "DEM", false, "control parameters for discrete element method (DEM) simulations\n");

  // write particle energy output
  Core::Utils::bool_parameter(
      "WRITE_PARTICLE_ENERGY", "no", "write particle energy output", &particledyndem);

  // write particle-wall interaction output
  Core::Utils::bool_parameter("WRITE_PARTICLE_WALL_INTERACTION", "no",
      "write particle-wall interaction output", &particledyndem);

  // type of normal contact law
  setStringToIntegralParameter<NormalContact>("NORMALCONTACTLAW", "NormalLinearSpring",
      "normal contact law for particles",
      tuple<std::string>("NormalLinearSpring", "NormalLinearSpringDamp", "NormalHertz",
          "NormalLeeHerrmann", "NormalKuwabaraKono", "NormalTsuji"),
      tuple<NormalContact>(Inpar::PARTICLE::NormalLinSpring, Inpar::PARTICLE::NormalLinSpringDamp,
          Inpar::PARTICLE::NormalHertz, Inpar::PARTICLE::NormalLeeHerrmann,
          Inpar::PARTICLE::NormalKuwabaraKono, Inpar::PARTICLE::NormalTsuji),
      &particledyndem);

  // type of tangential contact law
  setStringToIntegralParameter<TangentialContact>("TANGENTIALCONTACTLAW", "NoTangentialContact",
      "tangential contact law for particles",
      tuple<std::string>("NoTangentialContact", "TangentialLinSpringDamp"),
      tuple<TangentialContact>(
          Inpar::PARTICLE::NoTangentialContact, Inpar::PARTICLE::TangentialLinSpringDamp),
      &particledyndem);

  // type of rolling contact law
  setStringToIntegralParameter<RollingContact>("ROLLINGCONTACTLAW", "NoRollingContact",
      "rolling contact law for particles",
      tuple<std::string>("NoRollingContact", "RollingViscous", "RollingCoulomb"),
      tuple<RollingContact>(Inpar::PARTICLE::NoRollingContact, Inpar::PARTICLE::RollingViscous,
          Inpar::PARTICLE::RollingCoulomb),
      &particledyndem);

  // type of normal adhesion law
  setStringToIntegralParameter<AdhesionLaw>("ADHESIONLAW", "NoAdhesion",
      "type of adhesion law for particles",
      tuple<std::string>("NoAdhesion", "AdhesionVdWDMT", "AdhesionRegDMT"),
      tuple<AdhesionLaw>(Inpar::PARTICLE::NoAdhesion, Inpar::PARTICLE::AdhesionVdWDMT,
          Inpar::PARTICLE::AdhesionRegDMT),
      &particledyndem);

  // type of (random) surface energy distribution
  setStringToIntegralParameter<SurfaceEnergyDistribution>("ADHESION_SURFACE_ENERGY_DISTRIBUTION",
      "ConstantSurfaceEnergy", "type of (random) surface energy distribution",
      tuple<std::string>("ConstantSurfaceEnergy", "NormalSurfaceEnergyDistribution",
          "LogNormalSurfaceEnergyDistribution"),
      tuple<SurfaceEnergyDistribution>(Inpar::PARTICLE::ConstantSurfaceEnergy,
          Inpar::PARTICLE::NormalSurfaceEnergyDistribution,
          Inpar::PARTICLE::LogNormalSurfaceEnergyDistribution),
      &particledyndem);

  Core::Utils::double_parameter(
      "MIN_RADIUS", 0.0, "minimum allowed particle radius", &particledyndem);
  Core::Utils::double_parameter(
      "MAX_RADIUS", 0.0, "maximum allowed particle radius", &particledyndem);
  Core::Utils::double_parameter(
      "MAX_VELOCITY", -1.0, "maximum expected particle velocity", &particledyndem);

  // type of initial particle radius assignment
  setStringToIntegralParameter<InitialRadiusAssignment>("INITIAL_RADIUS",
      "RadiusFromParticleMaterial", "type of initial particle radius assignment",
      tuple<std::string>("RadiusFromParticleMaterial", "RadiusFromParticleInput",
          "NormalRadiusDistribution", "LogNormalRadiusDistribution"),
      tuple<InitialRadiusAssignment>(Inpar::PARTICLE::RadiusFromParticleMaterial,
          Inpar::PARTICLE::RadiusFromParticleInput, Inpar::PARTICLE::NormalRadiusDistribution,
          Inpar::PARTICLE::LogNormalRadiusDistribution),
      &particledyndem);

  Core::Utils::double_parameter("RADIUSDISTRIBUTION_SIGMA", -1.0,
      "sigma of random particle radius distribution", &particledyndem);

  Core::Utils::double_parameter(
      "REL_PENETRATION", -1.0, "maximum allowed relative penetration", &particledyndem);
  Core::Utils::double_parameter("NORMAL_STIFF", -1.0, "normal contact stiffness", &particledyndem);
  Core::Utils::double_parameter(
      "NORMAL_DAMP", -1.0, "normal contact damping parameter", &particledyndem);
  Core::Utils::double_parameter(
      "COEFF_RESTITUTION", -1.0, "coefficient of restitution", &particledyndem);
  Core::Utils::double_parameter("DAMP_REG_FAC", -1.0,
      "linearly regularized damping normal force in the interval "
      "\f$|g| < (\\text{DAMP_REG_FAC} \\cdot r_{\\min})\f$",
      &particledyndem);
  Core::Utils::bool_parameter(
      "TENSION_CUTOFF", "yes", "evaluate tension cutoff of normal contact force", &particledyndem);

  Core::Utils::double_parameter("POISSON_RATIO", -1.0, "poisson ratio", &particledyndem);
  Core::Utils::double_parameter("YOUNG_MODULUS", -1.0, "young's modulus", &particledyndem);

  Core::Utils::double_parameter(
      "FRICT_COEFF_TANG", -1.0, "friction coefficient for tangential contact", &particledyndem);
  Core::Utils::double_parameter(
      "FRICT_COEFF_ROLL", -1.0, "friction coefficient for rolling contact", &particledyndem);

  Core::Utils::double_parameter(
      "ADHESION_DISTANCE", -1.0, "adhesion distance between interacting surfaces", &particledyndem);

  Core::Utils::double_parameter(
      "ADHESION_MAX_CONTACT_PRESSURE", 0.0, "adhesion maximum contact pressure", &particledyndem);
  Core::Utils::double_parameter(
      "ADHESION_MAX_CONTACT_FORCE", 0.0, "adhesion maximum contact force", &particledyndem);
  Core::Utils::bool_parameter("ADHESION_USE_MAX_CONTACT_FORCE", "no",
      "use maximum contact force instead of maximum contact pressure", &particledyndem);

  Core::Utils::bool_parameter(
      "ADHESION_VDW_CURVE_SHIFT", "no", "shifts van-der-Waals-curve to g = 0", &particledyndem);

  Core::Utils::double_parameter(
      "ADHESION_HAMAKER", -1.0, "hamaker constant of van-der-Waals interaction", &particledyndem);

  Core::Utils::double_parameter("ADHESION_SURFACE_ENERGY", -1.0,
      "adhesion surface energy for the calculation of the pull-out force", &particledyndem);
  Core::Utils::double_parameter("ADHESION_SURFACE_ENERGY_DISTRIBUTION_VAR", -1.0,
      "variance of adhesion surface energy distribution", &particledyndem);
  Core::Utils::double_parameter("ADHESION_SURFACE_ENERGY_DISTRIBUTION_CUTOFF_FACTOR", -1.0,
      "adhesion surface energy distribution limited by multiple of variance", &particledyndem);
  Core::Utils::double_parameter("ADHESION_SURFACE_ENERGY_FACTOR", 1.0,
      "factor to calculate minimum adhesion surface energy", &particledyndem);
}

/*---------------------------------------------------------------------------*
 | set the particle conditions                                sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
void Inpar::PARTICLE::set_valid_conditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  /*-------------------------------------------------------------------------*
   | particle wall condition                                                 |
   *-------------------------------------------------------------------------*/
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfpartwall =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>("DESIGN SURFACE PARTICLE WALL",
          "ParticleWall", "Wall for particle interaction with (optional) material definition",
          Core::Conditions::ParticleWall, true, Core::Conditions::geometry_type_surface);

  add_named_int(surfpartwall, "MAT");

  condlist.push_back(surfpartwall);
}

FOUR_C_NAMESPACE_CLOSE

/*----------------------------------------------------------------------*/
/*!
\file inpar_particle.cpp

\brief Input parameters for particle problems

\level 1

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*-----------------------------------------------------------------------*/

#include "drt_validparameters.H"
#include "inpar_particle.H"
#include "inpar_parameterlist_utils.H"
#include "../drt_lib/drt_conditiondefinition.H"


void INPAR::PARTICLE::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& particledyn =
      list->sublist("PARTICLE DYNAMIC", false, "control parameters for particle problems\n");

  setStringToIntegralParameter<int>("DYNAMICTYP", "CentrDiff", "type of time integration control",
      tuple<std::string>("ExplicitEuler", "CentrDiff", "KickDrift", "RK2", "RK4"),
      tuple<int>(INPAR::PARTICLE::dyna_expleuler, INPAR::PARTICLE::dyna_centrdiff,
          INPAR::PARTICLE::dyna_kickdrift, INPAR::PARTICLE::dyna_rk2, INPAR::PARTICLE::dyna_rk4),
      &particledyn);

  setStringToIntegralParameter<int>("REPARTITIONSTRATEGY", "Adaptive",
      "type of employed repartitioning strategy", tuple<std::string>("Everydt", "Adaptive"),
      tuple<int>(repstr_everydt, repstr_adaptive), &particledyn);

  // Output type
  IntParameter("RESULTSEVRY", 1, "save displacements and contact forces every RESULTSEVRY steps",
      &particledyn);
  IntParameter("RESEVRYERGY", 0, "write system energies every requested step", &particledyn);
  IntParameter("RESEVRYREND", 0, "write SPH rendering every requested step", &particledyn);
  IntParameter("RESTARTEVRY", 1, "write restart possibility every RESTARTEVRY steps", &particledyn);
  IntParameter("CONTACTFORCESEVRY", 0,
      "output particle-particle and particle-wall contact forces to *.csv file every "
      "CONTACTFORCESEVRY steps",
      &particledyn);
  IntParameter("PARTICLESTATSEVRY", 0,
      "output particle statistics to *.csv file every PARTICLESTATSEVRY steps", &particledyn);

  // Time loop control
  DoubleParameter("TIMESTEP", 0.05, "time step size", &particledyn);
  IntParameter("NUMSTEP", 200, "maximum number of steps", &particledyn);
  DoubleParameter("MAXTIME", 5.0, "maximum time", &particledyn);

  setStringToIntegralParameter<int>("PARTICLE_INTERACTION", "None",
      "Interaction types for particle problems",
      tuple<std::string>(
          "None", "NormalContact_DEM", "NormalContact_MD", "NormalAndTangentialContact_DEM", "SPH"),
      tuple<int>(INPAR::PARTICLE::None, INPAR::PARTICLE::Normal_DEM, INPAR::PARTICLE::Normal_MD,
          INPAR::PARTICLE::NormalAndTang_DEM, INPAR::PARTICLE::SPH),
      &particledyn);

  setStringToIntegralParameter<int>("NORMAL_CONTACT_LAW", "LinearSpringDamp",
      "contact law for normal contact of particles",
      tuple<std::string>(
          "LinearSpring", "Hertz", "LinearSpringDamp", "LeeHerrmann", "KuwabaraKono", "Tsuji"),
      tuple<int>(INPAR::PARTICLE::LinSpring, INPAR::PARTICLE::Hertz, INPAR::PARTICLE::LinSpringDamp,
          INPAR::PARTICLE::LeeHerrmann, INPAR::PARTICLE::KuwabaraKono, INPAR::PARTICLE::Tsuji),
      &particledyn);

  setStringToIntegralParameter<int>("NORMAL_ADHESION_LAW", "None",
      "adhesion law governing normal contact of particles",
      tuple<std::string>(
          "None", "LinearSpring", "LinearSpringDamp", "vanderWaalsDMT", "regularDMT"),
      tuple<int>(INPAR::PARTICLE::adhesion_none, INPAR::PARTICLE::adhesion_linspring,
          INPAR::PARTICLE::adhesion_linspringdamp, INPAR::PARTICLE::adhesion_vdWDMT,
          INPAR::PARTICLE::adhesion_regDMT),
      &particledyn);

  setStringToIntegralParameter<int>("ROLLING_CONTACT_LAW", "None",
      "rolling contact law governing rolling contact of particles",
      tuple<std::string>("None", "Viscous", "Constant"),
      tuple<int>(INPAR::PARTICLE::rolling_none, INPAR::PARTICLE::rolling_viscous,
          INPAR::PARTICLE::rolling_constant),
      &particledyn);

  setStringToIntegralParameter<int>("WEIGHT_FUNCTION", "CubicBspline",
      "weight function for SPH interaction dynamics",
      tuple<std::string>("CubicBspline", "QuinticBspline", "SqrtHyperbola", "HyperbolaNoRsz"),
      tuple<int>(INPAR::PARTICLE::CubicBspline, INPAR::PARTICLE::QuinticBspline,
          INPAR::PARTICLE::SqrtHyperbola, INPAR::PARTICLE::HyperbolaNoRsz),
      &particledyn);

  setStringToIntegralParameter<int>("WEIGHT_FUNCTION_DIM", "WF_3D",
      "number of weight function space dimensions for SPH interaction dynamics",
      tuple<std::string>("WF_3D", "WF_2D", "WF_1D"),
      tuple<int>(INPAR::PARTICLE::WF_3D, INPAR::PARTICLE::WF_2D, INPAR::PARTICLE::WF_1D),
      &particledyn);

  setStringToIntegralParameter<int>("EQUATION_OF_STATE", "GenTait",
      "equation of state for SPH interaction dynamics", tuple<std::string>("GenTait", "IdealGas"),
      tuple<int>(INPAR::PARTICLE::GenTait, INPAR::PARTICLE::IdealGas), &particledyn);

  setStringToIntegralParameter<int>("EXTENDED_GHOSTING", "StandardGhosting", "Extend the ghosting",
      tuple<std::string>(
          "StandardGhosting", "BdryParticleGhosting", "WallElementGhosting", "AddLayerGhosting"),
      tuple<int>(INPAR::PARTICLE::StandardGhosting, INPAR::PARTICLE::BdryParticleGhosting,
          INPAR::PARTICLE::WallElementGhosting, INPAR::PARTICLE::AddLayerGhosting),
      &particledyn);

  setStringToIntegralParameter<int>("WALL_INTERACTION_TYPE", "BoundarParticle_NoSlip",
      "wall interaction type for SPH interactions",
      tuple<std::string>("BoundarParticle_NoSlip", "BoundarParticle_FreeSlip", "NoWallInteraction"),
      tuple<int>(INPAR::PARTICLE::BoundarParticle_NoSlip, INPAR::PARTICLE::BoundarParticle_FreeSlip,
          INPAR::PARTICLE::NoWallInteraction),
      &particledyn);

  setStringToIntegralParameter<int>("FREE_SURFACE_TYPE", "None",
      "type of free-surface treatment for SPH interactions",
      tuple<std::string>("None", "DensityIntegration", "InteriorReinitialization",
          "NormalizedReinitialization", "RandlesReinitialization", "TwoPhase"),
      tuple<int>(INPAR::PARTICLE::FreeSurface_None, INPAR::PARTICLE::DensityIntegration,
          INPAR::PARTICLE::InteriorReinitialization, INPAR::PARTICLE::NormalizedReinitialization,
          INPAR::PARTICLE::RandlesReinitialization, INPAR::PARTICLE::TwoPhase),
      &particledyn);

  setStringToIntegralParameter<int>("SURFACE_TENSION_TYPE", "ST_NONE",
      "type of surface tension forces for SPH interactions",
      tuple<std::string>("ST_NONE", "ST_VDW_INDIRECT", "ST_CONTI_ADAMI"),
      tuple<int>(INPAR::PARTICLE::ST_NONE, INPAR::PARTICLE::ST_VDW_INDIRECT,
          INPAR::PARTICLE::ST_CONTI_ADAMI),
      &particledyn);

  DoubleParameter("MIN_RADIUS", -1.0, "smallest particle radius", &particledyn);
  DoubleParameter("MAX_RADIUS", -1.0, "largest particle radius", &particledyn);
  DoubleParameter("REL_PENETRATION", -1.0, "relative particle-particle penetration", &particledyn);
  DoubleParameter("REL_PENETRATION_WALL", -1.0, "relative particle-wall penetration", &particledyn);
  DoubleParameter("MAX_VELOCITY", -1.0, "highest particle velocity", &particledyn);
  DoubleParameter("MAX_ANGULAR_VELOCITY", -1.0, "highest particle angular velocity", &particledyn);
  DoubleParameter("COEFF_RESTITUTION", -1.0, "coefficient of restitution", &particledyn);
  DoubleParameter(
      "COEFF_RESTITUTION_WALL", -1.0, "coefficient of restitution (wall)", &particledyn);
  DoubleParameter(
      "FRICT_COEFF_WALL", -1.0, "friction coefficient for contact particle-wall", &particledyn);
  DoubleParameter("FRICT_COEFF", -1.0, "dynamic friction coefficient for contact particle-particle",
      &particledyn);
  DoubleParameter("ROLL_FRICT_COEFF", -1.0,
      "dynamic rolling friction coefficient for contact particle-particle", &particledyn);
  DoubleParameter("ROLL_FRICT_COEFF_WALL", -1.0,
      "dynamic rolling friction coefficient for contact particle-wall", &particledyn);
  DoubleParameter("NORMAL_STIFF", -1.0, "stiffness for normal contact force", &particledyn);
  DoubleParameter(
      "NORMAL_STIFF_WALL", -1.0, "stiffness for normal contact force (wall)", &particledyn);
  DoubleParameter("TANG_STIFF", -1.0, "stiffness for tangential contact force", &particledyn);
  DoubleParameter(
      "NORMAL_DAMP", -1.0, "damping coefficient for normal contact force", &particledyn);
  DoubleParameter("NORMAL_DAMP_WALL", -1.0, "damping coefficient for normal contact force (wall)",
      &particledyn);
  DoubleParameter(
      "TANG_DAMP", -1.0, "damping coefficient for tangential contact force", &particledyn);
  DoubleParameter("TANG_DAMP_WALL", -1.0, "damping coefficient for tangential contact force (wall)",
      &particledyn);
  DoubleParameter("DAMP_REG_FAC", -1.0,
      "Apply linearly regularized damping normal force in the interval |g|< "
      "PARTICLE_REGDAMDBOUND*r_min",
      &particledyn);
  DoubleParameter("YOUNG_WALL", -1.0, "Young's modulus of wall", &particledyn);
  DoubleParameter("NUE_WALL", -1.0, "Possion's ratio of wall", &particledyn);
  BoolParameter("TENSION_CUTOFF", "no", "switch on/off tension cutoff", &particledyn);
  BoolParameter("MOVING_WALLS", "no", "switch on/off moving walls", &particledyn);
  DoubleParameter("RANDOM_AMPLITUDE", 0.0, "random value for initial position", &particledyn);

  setStringToIntegralParameter<int>("RADIUS_DISTRIBUTION", "None",
      "random distribution of particle radii", tuple<std::string>("None", "LogNormal", "Normal"),
      tuple<int>(INPAR::PARTICLE::radiusdistribution_none,
          INPAR::PARTICLE::radiusdistribution_lognormal,
          INPAR::PARTICLE::radiusdistribution_normal),
      &particledyn);

  setStringToIntegralParameter<int>("ADHESION_SURFACE_ENERGY_DISTRIBUTION", "None",
      "random distribution of adhesion surface energy",
      tuple<std::string>("None", "LogNormal", "Normal"),
      tuple<int>(INPAR::PARTICLE::adhesionsurfaceenergydistribution_none,
          INPAR::PARTICLE::adhesionsurfaceenergydistribution_lognormal,
          INPAR::PARTICLE::adhesionsurfaceenergydistribution_normal),
      &particledyn);

  DoubleParameter("RADIUS_DISTRIBUTION_SIGMA", -1.0,
      "standard deviation of particle radii distribution", &particledyn);
  IntParameter("RADIUS_CHANGE_FUNCT", -1, "number of curve governing radius change", &particledyn);
  BoolParameter("DENSITY_SUMMATION", "no",
      "determine density via summation formula instead of time integration", &particledyn);
  DoubleParameter("CONSISTENT_PROBLEM_VOLUME", -1.0,
      "prescribe problem volume and determine particle masses consistently based on this volume "
      "and the initial density",
      &particledyn);
  DoubleParameter("VISCOUS_DAMPING", -1.0,
      "apply artificial viscous damping force to particles in order to determine static "
      "equilibrium solutions",
      &particledyn);
  BoolParameter("NO_VELDIFF_TERM", "no",
      "Do not apply velocity difference tensor in case of transport velocity formulation",
      &particledyn);
  BoolParameter("CALC_ACC_VAR2", "no",
      "I apply alternative variant for discretization of pressure gradient and viscous forces "
      "according to Adami et al. 2013",
      &particledyn);
  setNumericStringParameter("GRAVITY_ACCELERATION", "0.0 0.0 0.0",
      "Acceleration due to gravity in particle simulations.", &particledyn);
  DoubleParameter("GRAVITY_RAMP_TIME", -1.0,
      "increase gravity force smoothly within a time interval GRAVITY_RAMP_TIME in the beginning "
      "of the simulation",
      &particledyn);
  DoubleParameter("SURFTENSION_RAMP_TIME", -1.0,
      "increase surface tension force smoothly within a time interval SURFTENSION_RAMP_TIME in the "
      "beginning of the simulation",
      &particledyn);
  DoubleParameter("BACKGROUND_PRESSURE", -1.0,
      "constant background pressure employed for modified particle convection "
      "velocity/acceleration in KickDrift time integrator (SPH only)",
      &particledyn);
  DoubleParameter("XSPH_DAMPFAC", -1.0, "damping factor of XSPH scheme", &particledyn);
  DoubleParameter("XSPH_STIFFAC", -1.0, "damping factor of XSPH scheme", &particledyn);
  BoolParameter("TRANSPORT_VELOCITY", "no",
      "apply particle convection velocity that differs from momentum velocity", &particledyn);
  DoubleParameter("SURFTENSFAC_FF", 0.0,
      "surface tension value/factor between fluid and void phase/second fluid", &particledyn);
  DoubleParameter("SURFTENSFAC_FS", 0.0,
      "surface tension value/factor between fluid and solid phase", &particledyn);
  DoubleParameter("SURFTENSPOT_A_FF", 2.0,
      "surface tension scale factor of repulsive part of pairwise interaction force fluid-fluid",
      &particledyn);
  DoubleParameter("SURFTENSPOT_A_FS", 2.0,
      "surface tension scale factor of repulsive part of pairwise interaction force fluid-solid",
      &particledyn);
  DoubleParameter("SURFTENSPOT_RATIO", 0.5,
      "range ratio of repulsive and attractive part of pairwise interaction force", &particledyn);
  IntParameter("MIN_VOIDPARTICLE_ID", -1,
      "all particles with an ID >= MIN_VOIDPARTICLE_ID are considered as void particles with "
      "colorfield = 0",
      &particledyn);
  BoolParameter("CONTACT_ANGLE_VAR2", "no",
      "Apply alternative model to enforce the static contact angle", &particledyn);

  setStringToIntegralParameter<int>("DIMENSION", "3D",
      "number of space dimensions for handling of quasi-2D problems with 3D particles",
      tuple<std::string>("3D", "2Dx", "2Dy", "2Dz"),
      tuple<int>(INPAR::PARTICLE::particle_3D, INPAR::PARTICLE::particle_2Dx,
          INPAR::PARTICLE::particle_2Dy, INPAR::PARTICLE::particle_2Dz),
      &particledyn);

  setStringToIntegralParameter<int>("RENDERING", "None",
      "rendering type for smoothed particle hydrodynamics",
      tuple<std::string>("None", "Standard", "Normalized"),
      tuple<int>(INPAR::PARTICLE::NoRendering, INPAR::PARTICLE::StandardRendering,
          INPAR::PARTICLE::NormalizedRendering),
      &particledyn);

  setStringToIntegralParameter<int>("RENDERING_OUTPUT", "DiscretAndMatlab",
      "rendering output for smoothed particle hydrodynamics",
      tuple<std::string>("DiscretAndMatlab", "Discret", "Matlab"),
      tuple<int>(
          INPAR::PARTICLE::DiscretAndMatlab, INPAR::PARTICLE::Discret, INPAR::PARTICLE::Matlab),
      &particledyn);

  setStringToIntegralParameter<int>("RENDERING_BDRYPARTICLE", "WithBdryParticle",
      "consider boundary particles in rendering for smoothed particle hydrodynamics",
      tuple<std::string>("WithBdryParticle", "NoBdryParticle"),
      tuple<int>(INPAR::PARTICLE::WithBdryParticle, INPAR::PARTICLE::NoBdryParticle), &particledyn);

  IntParameter("AVRG_REND_STEPS", 1,
      "Average rendering vectors over AVRG_REND_STEPS time steps within one output inverval",
      &particledyn);

  DoubleParameter("ADHESION_EQ_GAP", 0.0,
      "gap between two particles or between particle and wall in adhesion equilibrium",
      &particledyn);
  DoubleParameter(
      "ADHESION_NORMAL_STIFF", -1.0, "stiffness for normal adhesion force", &particledyn);
  DoubleParameter(
      "ADHESION_NORMAL_DAMP", -1.0, "damping coefficient for normal adhesion force", &particledyn);
  DoubleParameter(
      "ADHESION_NORMAL_EPS", -1.0, "depth of Lennard-Jones adhesion potential well", &particledyn);
  DoubleParameter("ADHESION_MAX_FORCE", -1.0, "maximum adhesion force", &particledyn);
  DoubleParameter(
      "ADHESION_MAX_DISP", -1.0, "maximum displacement from adhesion equilibrium", &particledyn);
  DoubleParameter("ADHESION_SURFACE_ENERGY", -1.0,
      "surface energy density for the calculation of the pull-out force", &particledyn);
  DoubleParameter("ADHESION_SURFACE_ENERGY_FACTOR", 1.0,
      "factor to calculate minimum surface energy out of surface energy", &particledyn);
  DoubleParameter("ADHESION_SURFACE_ENERGY_WALL", -1.0,
      "surface energy density between particle and wall", &particledyn);
  DoubleParameter("ADHESION_SURFACE_ENERGY_WALL_FACTOR", -1.0,
      "factor to calculate minimum surface energy out of surface energy of wall", &particledyn);
  DoubleParameter("ADHESION_SURFACE_ENERGY_CUTOFF_FACTOR", -1.0,
      "surface energy distribution limited by multiple of standard deviation", &particledyn);
  DoubleParameter("ADHESION_SURFACE_ENERGY_DISTRIBUTION_SIGMA", -1.0,
      "standard deviation of adhesion surface energy distribution", &particledyn);
  DoubleParameter("ADHESION_NORMAL_FORCE_PULLOUT_MAX", -1.0, "maximum pullout force", &particledyn);
  DoubleParameter("ADHESION_NORMAL_FORCE_PULLOUT_MIN", -1.0, "minimum pullout force", &particledyn);
  DoubleParameter("ADHESION_MAX_CONTACT_PRESSURE", 0.0,
      "adhesion maximum contact pressure particle particle contact", &particledyn);
  DoubleParameter("ADHESION_MAX_CONTACT_PRESSURE_WALL", 0.0,
      "adhesion maximum contact pressure particle wall contact", &particledyn);
  BoolParameter("ADHESION_VARIANT_MAX_CONTACT_FORCE", "no",
      "take maximum contact force instead of pressure for adhesion ramp function", &particledyn);
  IntParameter("ADHESION_VDW_CURVE_SHIFT", 0, "shifts van-der-Waals-curve to g = 0", &particledyn);
  IntParameter("ADHESION_MAXWALLELEID", -1,
      "Maximal wall element ID to which adhesion forces are applied (-1: apply to all wall "
      "elements)",
      &particledyn);
  IntParameter("ADHESION_MINWALLELEID", -1,
      "Minimal wall element ID to which adhesion forces are applied (-1: apply to all wall "
      "elements)",
      &particledyn);
}


void INPAR::PARTICLE::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // particle inflow condition

  std::vector<Teuchos::RCP<ConditionComponent>> particleinflowcomponents;
  // two vertices describing the bounding box for the inflow
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("vertex1")));
  particleinflowcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("vertex1", 3)));
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("vertex2")));
  particleinflowcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("vertex2", 3)));
  // number of particles to inflow in each direction in bounding box
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("num_per_dir")));
  particleinflowcomponents.push_back(
      Teuchos::rcp(new IntVectorConditionComponent("num_per_dir", 3)));
  // particle inflow velocity
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("inflow_vel")));
  particleinflowcomponents.push_back(
      Teuchos::rcp(new RealVectorConditionComponent("inflow_vel", 3)));
  // particle inflow velocity can be superposed with a time curve
  particleinflowcomponents.push_back(
      Teuchos::rcp(new SeparatorConditionComponent("inflow_vel_funct")));
  particleinflowcomponents.push_back(Teuchos::rcp(new IntConditionComponent("inflow_vel_funct")));
  // inflow frequency of particles
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("inflow_freq")));
  particleinflowcomponents.push_back(Teuchos::rcp(new RealConditionComponent("inflow_freq")));
  // time delay of inflow particles
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("timedelay")));
  particleinflowcomponents.push_back(Teuchos::rcp(new RealConditionComponent("timedelay")));
  // inflow stopping time
  particleinflowcomponents.push_back(
      Teuchos::rcp(new SeparatorConditionComponent("stopinflowtime")));
  particleinflowcomponents.push_back(Teuchos::rcp(new RealConditionComponent("stopinflowtime")));


  Teuchos::RCP<ConditionDefinition> particlecond = Teuchos::rcp(new ConditionDefinition(
      "DESIGN PARTICLE INFLOW CONDITION", "ParticleInflow", "Particle Inflow Condition",
      DRT::Condition::ParticleInflow, false, DRT::Condition::Particle));


  for (unsigned i = 0; i < particleinflowcomponents.size(); ++i)
  {
    particlecond->AddComponent(particleinflowcomponents[i]);
  }

  condlist.push_back(particlecond);

  /*--------------------------------------------------------------------*/
  // particle init radius condition

  std::vector<Teuchos::RCP<ConditionComponent>> particleinitradiuscomponents;

  particleinitradiuscomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("SCALAR")));
  particleinitradiuscomponents.push_back(
      Teuchos::rcp(new RealVectorConditionComponent("SCALAR", 1)));

  particleinitradiuscomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("FUNCT")));
  particleinitradiuscomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("FUNCT", 1)));

  Teuchos::RCP<ConditionDefinition> particleradiuscond =
      Teuchos::rcp(new ConditionDefinition("DESIGN PARTICLE INIT RADIUS CONDITIONS",
          "InitialParticleRadius", "Particle Initial Radius Condition",
          DRT::Condition::ParticleInitRadius, false, DRT::Condition::Particle));


  for (unsigned i = 0; i < particleinitradiuscomponents.size(); ++i)
  {
    particleradiuscond->AddComponent(particleinitradiuscomponents[i]);
  }

  condlist.push_back(particleradiuscond);

  /*--------------------------------------------------------------------*/
  // particle wall condition

  std::vector<Teuchos::RCP<ConditionComponent>> particlewallcomponents;
  particlewallcomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

  Teuchos::RCP<ConditionDefinition> surfpartwall = Teuchos::rcp(new ConditionDefinition(
      "DESIGN SURFACE PARTICLE WALL", "ParticleWall", "Wall for particle collisions",
      DRT::Condition::ParticleWall, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < particlewallcomponents.size(); ++i)
  {
    surfpartwall->AddComponent(particlewallcomponents[i]);
  }

  condlist.push_back(surfpartwall);

  return;
}

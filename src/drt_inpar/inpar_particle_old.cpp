/*----------------------------------------------------------------------*/
/*!
\brief Input parameters for particle problems (old implementation)

\level 1

\maintainer  Sebastian Fuchs

*-----------------------------------------------------------------------*/

#include "drt_validparameters.H"
#include "inpar_parameterlist_utils.H"
#include "../drt_lib/drt_conditiondefinition.H"
#include "inpar_particle_old.H"


void INPAR::PARTICLEOLD::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& particledyn = list->sublist("PARTICLE DYNAMIC OLD", false,
      "control parameters for particle problems (old implementation)\n");

  setStringToIntegralParameter<int>("DYNAMICTYP", "CentrDiff", "type of time integration control",
      tuple<std::string>("ExplicitEuler", "CentrDiff", "RK2", "RK4"),
      tuple<int>(INPAR::PARTICLEOLD::dyna_expleuler, INPAR::PARTICLEOLD::dyna_centrdiff,
          INPAR::PARTICLEOLD::dyna_rk2, INPAR::PARTICLEOLD::dyna_rk4),
      &particledyn);

  setStringToIntegralParameter<int>("REPARTITIONSTRATEGY", "Adaptive",
      "type of employed repartitioning strategy", tuple<std::string>("Everydt", "Adaptive"),
      tuple<int>(repstr_everydt, repstr_adaptive), &particledyn);

  // Output type
  IntParameter("RESULTSEVRY", 1, "save displacements and contact forces every RESULTSEVRY steps",
      &particledyn);
  IntParameter("RESEVRYERGY", 0, "write system energies every requested step", &particledyn);
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
          "None", "NormalContact_DEM", "NormalContact_MD", "NormalAndTangentialContact_DEM"),
      tuple<int>(INPAR::PARTICLEOLD::None, INPAR::PARTICLEOLD::Normal_DEM,
          INPAR::PARTICLEOLD::Normal_MD, INPAR::PARTICLEOLD::NormalAndTang_DEM),
      &particledyn);

  setStringToIntegralParameter<int>("NORMAL_CONTACT_LAW", "LinearSpringDamp",
      "contact law for normal contact of particles",
      tuple<std::string>(
          "LinearSpring", "Hertz", "LinearSpringDamp", "LeeHerrmann", "KuwabaraKono", "Tsuji"),
      tuple<int>(INPAR::PARTICLEOLD::LinSpring, INPAR::PARTICLEOLD::Hertz,
          INPAR::PARTICLEOLD::LinSpringDamp, INPAR::PARTICLEOLD::LeeHerrmann,
          INPAR::PARTICLEOLD::KuwabaraKono, INPAR::PARTICLEOLD::Tsuji),
      &particledyn);

  setStringToIntegralParameter<int>("NORMAL_ADHESION_LAW", "None",
      "adhesion law governing normal contact of particles",
      tuple<std::string>(
          "None", "LinearSpring", "LinearSpringDamp", "vanderWaalsDMT", "regularDMT"),
      tuple<int>(INPAR::PARTICLEOLD::adhesion_none, INPAR::PARTICLEOLD::adhesion_linspring,
          INPAR::PARTICLEOLD::adhesion_linspringdamp, INPAR::PARTICLEOLD::adhesion_vdWDMT,
          INPAR::PARTICLEOLD::adhesion_regDMT),
      &particledyn);

  setStringToIntegralParameter<int>("ROLLING_CONTACT_LAW", "None",
      "rolling contact law governing rolling contact of particles",
      tuple<std::string>("None", "Viscous", "Constant"),
      tuple<int>(INPAR::PARTICLEOLD::rolling_none, INPAR::PARTICLEOLD::rolling_viscous,
          INPAR::PARTICLEOLD::rolling_constant),
      &particledyn);

  setStringToIntegralParameter<int>("EXTENDED_GHOSTING", "StandardGhosting", "Extend the ghosting",
      tuple<std::string>("StandardGhosting", "WallElementGhosting", "AddLayerGhosting"),
      tuple<int>(INPAR::PARTICLEOLD::StandardGhosting, INPAR::PARTICLEOLD::WallElementGhosting,
          INPAR::PARTICLEOLD::AddLayerGhosting),
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
      tuple<int>(INPAR::PARTICLEOLD::radiusdistribution_none,
          INPAR::PARTICLEOLD::radiusdistribution_lognormal,
          INPAR::PARTICLEOLD::radiusdistribution_normal),
      &particledyn);

  setStringToIntegralParameter<int>("ADHESION_SURFACE_ENERGY_DISTRIBUTION", "None",
      "random distribution of adhesion surface energy",
      tuple<std::string>("None", "LogNormal", "Normal"),
      tuple<int>(INPAR::PARTICLEOLD::adhesionsurfaceenergydistribution_none,
          INPAR::PARTICLEOLD::adhesionsurfaceenergydistribution_lognormal,
          INPAR::PARTICLEOLD::adhesionsurfaceenergydistribution_normal),
      &particledyn);

  DoubleParameter("RADIUS_DISTRIBUTION_SIGMA", -1.0,
      "standard deviation of particle radii distribution", &particledyn);
  IntParameter("RADIUS_CHANGE_FUNCT", -1, "number of curve governing radius change", &particledyn);
  setNumericStringParameter("GRAVITY_ACCELERATION", "0.0 0.0 0.0",
      "Acceleration due to gravity in particle simulations.", &particledyn);

  setStringToIntegralParameter<int>("DIMENSION", "3D",
      "number of space dimensions for handling of quasi-2D problems with 3D particles",
      tuple<std::string>("3D", "2Dx", "2Dy", "2Dz"),
      tuple<int>(INPAR::PARTICLEOLD::particle_3D, INPAR::PARTICLEOLD::particle_2Dx,
          INPAR::PARTICLEOLD::particle_2Dy, INPAR::PARTICLEOLD::particle_2Dz),
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


void INPAR::PARTICLEOLD::SetValidConditions(
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

  return;
}

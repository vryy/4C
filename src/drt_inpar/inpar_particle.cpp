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
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& particledyn = list->sublist(
       "PARTICLE DYNAMIC",
       false,
       "control parameters for particle problems\n");

   setStringToIntegralParameter<int>("DYNAMICTYP","CentrDiff",
                                "type of time integration control",
                                tuple<std::string>(
                                  "ExplicitEuler",
                                  "CentrDiff",
                                  "KickDrift",
                                  "RK2",
                                  "RK4",
                                  "GenAlpha"),
                                tuple<int>(
                                  INPAR::PARTICLE::dyna_expleuler,
                                  INPAR::PARTICLE::dyna_centrdiff,
                                  INPAR::PARTICLE::dyna_kickdrift,
                                  INPAR::PARTICLE::dyna_rk2,
                                  INPAR::PARTICLE::dyna_rk4,
                                  INPAR::PARTICLE::dyna_genAlpha
                                ),
                                &particledyn);

   // Output type
   IntParameter("RESULTSEVRY",1,"save displacements and contact forces every RESULTSEVRY steps",&particledyn);
   IntParameter("RESEVRYERGY",0,"write system energies every requested step",&particledyn);
   IntParameter("RESEVRYREND",0,"write meshfree rendering every requested step",&particledyn);
   IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&particledyn);
   IntParameter("CONTACTFORCESEVRY",0,"output particle-particle and particle-wall contact forces to *.csv file every CONTACTFORCESEVRY steps",&particledyn);
   IntParameter("PARTICLESTATSEVRY",0,"output particle statistics to *.csv file every PARTICLESTATSEVRY steps",&particledyn);

   // Time loop control
   DoubleParameter("TIMESTEP",0.05,"time step size",&particledyn);
   IntParameter("NUMSTEP",200,"maximum number of steps",&particledyn);
   DoubleParameter("MAXTIME",5.0,"maximum time",&particledyn);

   setStringToIntegralParameter<int>(
                               "PARTICLE_INTERACTION","None",
                               "Interaction types for particle or meshfree problems",
                               tuple<std::string>(
                                 "None",
                                 "NormalContact_DEM",
                                 "NormalContact_MD",
                                 "NormalAndTangentialContact_DEM",
                                 "NormalContact_DEM_thermo",
                                 "MeshFree"
                                 ),
                               tuple<int>(
                                 INPAR::PARTICLE::None,
                                 INPAR::PARTICLE::Normal_DEM,
                                 INPAR::PARTICLE::Normal_MD,
                                 INPAR::PARTICLE::NormalAndTang_DEM,
                                 INPAR::PARTICLE::Normal_DEM_thermo,
                                 INPAR::PARTICLE::MeshFree
                                 ),
                               &particledyn);

   setStringToIntegralParameter<int>(
                               "NORMAL_CONTACT_LAW","LinearSpringDamp",
                               "contact law for normal contact of particles",
                               tuple<std::string>(
                                 "LinearSpring",
                                 "Hertz",
                                 "LinearSpringDamp",
                                 "LeeHerrmann",
                                 "KuwabaraKono",
                                 "Tsuji"
                                 ),
                               tuple<int>(
                                 INPAR::PARTICLE::LinSpring,
                                 INPAR::PARTICLE::Hertz,
                                 INPAR::PARTICLE::LinSpringDamp,
                                 INPAR::PARTICLE::LeeHerrmann,
                                 INPAR::PARTICLE::KuwabaraKono,
                                 INPAR::PARTICLE::Tsuji
                                 ),
                               &particledyn);

   setStringToIntegralParameter<int>(
                               "NORMAL_ADHESION_LAW","None",
                               "adhesion law governing normal contact of particles",
                               tuple<std::string>(
                                 "None",
                                 "LennardJones",
                                 "LinearSpring",
                                 "LinearSpringDamp"
                                 ),
                               tuple<int>(
                                 INPAR::PARTICLE::adhesion_none,
                                 INPAR::PARTICLE::adhesion_lennardjones,
                                 INPAR::PARTICLE::adhesion_linspring,
                                 INPAR::PARTICLE::adhesion_linspringdamp
                                 ),
                               &particledyn);

   setStringToIntegralParameter<int>(
                               "WEIGHT_FUNCTION","CubicBspline",
                               "weight function for meshFree interaction dynamics",
                               tuple<std::string>(
                                 "CubicBspline",
                                 "QuinticBspline",
                                 "SqrtHyperbola",
                                 "HyperbolaNoRsz"
                                 ),
                               tuple<int>(
                                 INPAR::PARTICLE::CubicBspline,
                                 INPAR::PARTICLE::QuinticBspline,
                                 INPAR::PARTICLE::SqrtHyperbola,
                                 INPAR::PARTICLE::HyperbolaNoRsz
                                 ),
                               &particledyn);

   setStringToIntegralParameter<int>(
                                  "WEIGHT_FUNCTION_DIM","WF_3D",
                                  "number of weight function space dimensions for meshFree interaction dynamics",
                                  tuple<std::string>(
                                    "WF_3D",
                                    "WF_2D",
                                    "WF_1D"
                                    ),
                                  tuple<int>(
                                    INPAR::PARTICLE::WF_3D,
                                    INPAR::PARTICLE::WF_2D,
                                    INPAR::PARTICLE::WF_1D
                                    ),
                                  &particledyn);

   setStringToIntegralParameter<int>(
                                  "EXTENDED_GHOSTING","StandardGhosting",
                                  "Extend the ghosting",
                                  tuple<std::string>(
                                    "StandardGhosting",
                                    "BdryParticleGhosting",
                                    "AddLayerGhosting"
                                    ),
                                  tuple<int>(
                                    INPAR::PARTICLE::StandardGhosting,
                                    INPAR::PARTICLE::BdryParticleGhosting,
                                    INPAR::PARTICLE::AddLayerGhosting
                                    ),
                                  &particledyn);

   setStringToIntegralParameter<int>(
                               "WALL_INTERACTION_TYPE","InitParticle",
                               "wall interaction type for MeshFree interactions",
                               tuple<std::string>(
                                 "InitParticle",
                                 "Mirror",
                                 "Custom",
                                 "BoundarParticle_NoSlip",
                                 "BoundarParticle_FreeSlip",
                                 "NoWallInteraction"
                                 ),
                               tuple<int>(
                                 INPAR::PARTICLE::InitParticle,
                                 INPAR::PARTICLE::Mirror,
                                 INPAR::PARTICLE::Custom,
                                 INPAR::PARTICLE::BoundarParticle_NoSlip,
                                 INPAR::PARTICLE::BoundarParticle_FreeSlip,
                                 INPAR::PARTICLE::NoWallInteraction
                                 ),
                               &particledyn);

   IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for structural problems",&particledyn);
   DoubleParameter("ERROR_TOLL",1e-6,"tolerance of the error for implicit schemes",&particledyn);
   IntParameter("ITER_MAX",10,"maximum iteration per time step for implicit schemes",&particledyn);
   DoubleParameter("ALPHA_MIN",1e-6,"cap of the alpha parameter for the divergence free scheme for meshfree",&particledyn);
   DoubleParameter("WALL_FAKE_DENSITY",-1.0,"fake density for meshfree dynamics, in case of -1 the coefficients are extracted from the initial values of the particle material parameters",&particledyn);
   DoubleParameter("WALL_FAKE_MASS",-1.0,"fake mass of the wall element for meshfree dynamics, in case of -1 the coefficients are extracted from the initial values of the particle material parameters",&particledyn);
   DoubleParameter("WALL_FAKE_PRESSURE",-1.0,"fake pressure of the wall element for meshfree dynamics, in case of -1 the coefficients are extracted from the initial values of the particle material parameters",&particledyn);
   DoubleParameter("MIN_RADIUS",-1.0,"smallest particle radius",&particledyn);
   DoubleParameter("MAX_RADIUS",-1.0,"largest particle radius",&particledyn);
   DoubleParameter("REL_PENETRATION",-1.0,"relative particle-particle penetration",&particledyn);
   DoubleParameter("REL_PENETRATION_WALL",-1.0,"relative particle-wall penetration",&particledyn);
   DoubleParameter("MAX_VELOCITY",-1.0,"highest particle velocity",&particledyn);
   DoubleParameter("COEFF_RESTITUTION",-1.0,"coefficient of restitution",&particledyn);
   DoubleParameter("COEFF_RESTITUTION_WALL",-1.0,"coefficient of restitution (wall)",&particledyn);
   DoubleParameter("FRICT_COEFF_WALL",-1.0,"friction coefficient for contact particle-wall",&particledyn);
   DoubleParameter("FRICT_COEFF",-1.0,"dynamic friction coefficient for contact particle-particle",&particledyn);
   DoubleParameter("NORMAL_STIFF",-1.0,"stiffness for normal contact force",&particledyn);
   DoubleParameter("NORMAL_STIFF_WALL",-1.0,"stiffness for normal contact force (wall)",&particledyn);
   DoubleParameter("TANG_STIFF",-1.0,"stiffness for tangential contact force",&particledyn);
   DoubleParameter("NORMAL_DAMP",-1.0,"damping coefficient for normal contact force",&particledyn);
   DoubleParameter("NORMAL_DAMP_WALL",-1.0,"damping coefficient for normal contact force (wall)",&particledyn);
   DoubleParameter("TANG_DAMP",-1.0,"damping coefficient for tangential contact force",&particledyn);
   DoubleParameter("TANG_DAMP_WALL",-1.0,"damping coefficient for tangential contact force (wall)",&particledyn);
   DoubleParameter("YOUNG_WALL",-1.0,"Young's modulus of wall",&particledyn);
   DoubleParameter("NUE_WALL",-1.0,"Possion's ratio of wall",&particledyn);
   BoolParameter("TENSION_CUTOFF","no","switch on/off tension cutoff",&particledyn);
   BoolParameter("MOVING_WALLS","no","switch on/off moving walls",&particledyn);
   IntParameter("TRANSFER_EVERY",1,"transfer particles every TRANSFER_EVERY steps",&particledyn);
   DoubleParameter("RANDOM_AMPLITUDE",0.0,"random value for initial position",&particledyn);

   setStringToIntegralParameter<int>(
                               "RADIUS_DISTRIBUTION","None",
                               "random distribution of particle radii",
                               tuple<std::string>(
                                 "None",
                                 "LogNormal",
                                 "Normal"
                                 ),
                               tuple<int>(
                                 INPAR::PARTICLE::radiusdistribution_none,
                                 INPAR::PARTICLE::radiusdistribution_lognormal,
                                 INPAR::PARTICLE::radiusdistribution_normal
                                 ),
                               &particledyn);

   DoubleParameter("RADIUS_DISTRIBUTION_SIGMA",-1.0,"standard deviation of particle radii distribution",&particledyn);
   IntParameter("RADIUS_CHANGE_FUNCT",-1,"number of curve governing radius change",&particledyn);
   BoolParameter("DENSITY_SUMMATION","no","determine density via summation formula instead of time integration",&particledyn);
   DoubleParameter("CONSISTENT_PROBLEM_VOLUME",-1.0,"prescribe problem volume and determine particle masses consistently based on this volume and the initial density",&particledyn);
   DoubleParameter("VISCOUS_DAMPING",-1.0,"apply artificial viscous damping force to particles in order to determine static equilibrium solutions",&particledyn);
   BoolParameter("SOLVE_THERMAL_PROBLEM","yes","solve also the thermal problem?",&particledyn);
   BoolParameter("NO_VELDIFF_TERM","no","Do not apply velocity difference tensor in case of transport velocity formulation",&particledyn);
   BoolParameter("CALC_ACC_VAR2","no","I apply alternative variant for discretization of pressure gradient and viscous forces according to Adami et al. 2013",&particledyn);
   setNumericStringParameter("GRAVITY_ACCELERATION","0.0 0.0 0.0",
                             "Acceleration due to gravity in particle simulations.",
                             &particledyn);
   DoubleParameter("GRAVITY_RAMP_TIME",-1.0,"increase gravity force smoothly within a time interval GRAVITY_RAMP_TIME in the beginning of the simulation",&particledyn);
   DoubleParameter("BACKGROUND_PRESSURE",-1.0,"constant background pressure employed for modified particle convection velocity/acceleration in KickDrift time integrator (SPH only)",&particledyn);

   setStringToIntegralParameter<int>("DIMENSION","3D",
                                "number of space dimensions for handling of quasi-2D problems with 3D particles",
                                tuple<std::string>(
                                  "3D",
                                  "2Dx",
                                  "2Dy",
                                  "2Dz"),
                                tuple<int>(
                                  INPAR::PARTICLE::particle_3D,
                                  INPAR::PARTICLE::particle_2Dx,
                                  INPAR::PARTICLE::particle_2Dy,
                                  INPAR::PARTICLE::particle_2Dz),
                                &particledyn);

   setStringToIntegralParameter<int>(
                               "RENDERING","None",
                               "rendering type for smoothed particle hydrodynamics",
                               tuple<std::string>(
                                 "None",
                                 "Standard",
                                 "Normalized"
                                 ),
                               tuple<int>(
                                 INPAR::PARTICLE::NoRendering,
                                 INPAR::PARTICLE::StandardRendering,
                                 INPAR::PARTICLE::NormalizedRendering
                                 ),
                               &particledyn);

   setStringToIntegralParameter<int>(
                               "RENDERING_OUTPUT","DiscretAndMatlab",
                               "rendering output for smoothed particle hydrodynamics",
                               tuple<std::string>(
                                 "DiscretAndMatlab",
                                 "Discret",
                                 "Matlab"
                                 ),
                               tuple<int>(
                                 INPAR::PARTICLE::DiscretAndMatlab,
                                 INPAR::PARTICLE::Discret,
                                 INPAR::PARTICLE::Matlab
                                 ),
                               &particledyn);

   setStringToIntegralParameter<int>(
                               "RENDERING_BDRYPARTICLE","WithBdryParticle",
                               "consider boundary particles in rendering for smoothed particle hydrodynamics",
                               tuple<std::string>(
                                 "WithBdryParticle",
                                 "NoBdryParticle"
                                 ),
                               tuple<int>(
                                 INPAR::PARTICLE::WithBdryParticle,
                                 INPAR::PARTICLE::NoBdryParticle
                                 ),
                               &particledyn);

   IntParameter("AVRG_REND_STEPS",1,"Average rendering vectors over AVRG_REND_STEPS time steps within one output inverval",&particledyn);

   DoubleParameter("ADHESION_EQ_GAP",0.0,"gap between two particles or between particle and wall in adhesion equilibrium",&particledyn);
   DoubleParameter("ADHESION_NORMAL_STIFF",-1.0,"stiffness for normal adhesion force",&particledyn);
   DoubleParameter("ADHESION_NORMAL_DAMP",-1.0,"damping coefficient for normal adhesion force",&particledyn);
   DoubleParameter("ADHESION_NORMAL_EPS",-1.0,"depth of Lennard-Jones adhesion potential well",&particledyn);
   DoubleParameter("ADHESION_MAX_FORCE",-1.0,"maximum adhesion force",&particledyn);
   DoubleParameter("ADHESION_MAX_DISP",-1.0,"maximum displacement from adhesion equilibrium",&particledyn);
   /*----------------------------------------------------------------------*/
   /* parameters for generalised-alpha integrator */
   Teuchos::ParameterList& genalpha = particledyn.sublist("GENALPHA",false,"");

   setStringToIntegralParameter<int>("GENAVG","TrLike",
                                "mid-average type of internal forces",
                                tuple<std::string>(
                                  "Vague",
                                  "ImrLike",
                                  "TrLike"),
                                tuple<int>(
                                  midavg_vague,
                                  midavg_imrlike,
                                  midavg_trlike),
                                &genalpha);
   setStringToIntegralParameter<int>("GRADRES_APPROX","OnlyHess",
                                "type of approximation of the Generalised-alpha's residual of the gradient",
                                tuple<std::string>(
                                  "Full",
                                  "OnlyHess"),
                                tuple<int>(
                                  gaapprox_full,
                                  gaapprox_onlyhess),
                                &genalpha);


   DoubleParameter("BETA",0.25,"Generalised-alpha factor in (0,1/2]",&genalpha);
   DoubleParameter("GAMMA",0.5,"Generalised-alpha factor in (0,1]",&genalpha);
   DoubleParameter("ALPHA_M",0.5,"Generalised-alpha factor in [0,1)",&genalpha);
   DoubleParameter("ALPHA_F",0.5,"Generalised-alpha factor in [0,1)",&genalpha);
   DoubleParameter("RHO_INF",-1.0,"Generalised-alpha factor in [0,1]",&genalpha);
   DoubleParameter("TOL",1e-5,"Generalised-alpha tollerance",&genalpha);
   IntParameter("MAXIT",10,"Generalised-alpha max number of iterations",&genalpha);
   IntParameter("MAXIT",10,"Generalised-alpha max number of iterations",&genalpha);
}


void INPAR::PARTICLE::SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // particle inflow condition

  std::vector<Teuchos::RCP<ConditionComponent> > particleinflowcomponents;
  // two vertices describing the bounding box for the inflow
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("vertex1")));
  particleinflowcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("vertex1",3)));
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("vertex2")));
  particleinflowcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("vertex2",3)));
  // number of particles to inflow in each direction in bounding box
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("num_per_dir")));
  particleinflowcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("num_per_dir",3)));
  // particle inflow velocity
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("inflow_vel")));
  particleinflowcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("inflow_vel",3)));
  // particle inflow velocity can be superposed with a time curve
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("inflow_vel_funct")));
  particleinflowcomponents.push_back(Teuchos::rcp(new IntConditionComponent("inflow_vel_funct")));
  // inflow frequency of particles
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("inflow_freq")));
  particleinflowcomponents.push_back(Teuchos::rcp(new RealConditionComponent("inflow_freq")));
  // time delay of inflow particles
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("timedelay")));
  particleinflowcomponents.push_back(Teuchos::rcp(new RealConditionComponent("timedelay")));
  // inflow stopping time
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("stopinflowtime")));
  particleinflowcomponents.push_back(Teuchos::rcp(new RealConditionComponent("stopinflowtime")));


  Teuchos::RCP<ConditionDefinition> particlecond =
    Teuchos::rcp(new ConditionDefinition("DESIGN PARTICLE INFLOW CONDITION",
                                         "ParticleInflow",
                                         "Particle Inflow Condition",
                                         DRT::Condition::ParticleInflow,
                                         false,
                                         DRT::Condition::Particle));


  for (unsigned i=0; i<particleinflowcomponents.size(); ++i)
  {
    particlecond->AddComponent(particleinflowcomponents[i]);
  }

  condlist.push_back(particlecond);
  /*--------------------------------------------------------------------*/
    // particle heat source condition

    std::vector<Teuchos::RCP<ConditionComponent> > particleHSComponents;
    // two vertices describing the bounding box for the Heat Source
    particleHSComponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("vertex0")));
    particleHSComponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("vertex0",3)));
    particleHSComponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("vertex1")));
    particleHSComponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("vertex1",3)));
    // intake Q/s
    particleHSComponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("HSQDot")));
    particleHSComponents.push_back(Teuchos::rcp(new RealConditionComponent("HSQDot")));
    // t start
    particleHSComponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("HSTstart")));
    particleHSComponents.push_back(Teuchos::rcp(new RealConditionComponent("HSTstart")));
    // t end
    particleHSComponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("HSTend")));
    particleHSComponents.push_back(Teuchos::rcp(new RealConditionComponent("HSTend")));


    Teuchos::RCP<ConditionDefinition> particleHScond =
      Teuchos::rcp(new ConditionDefinition("DESIGN HEAT SOURCE CONDITION",
                                           "ParticleHeatSource",
                                           "Particle Heat Source Condition",
                                           DRT::Condition::ParticleHeatSource,
                                           false,
                                           DRT::Condition::Particle));


    for (unsigned i=0; i<particleHSComponents.size(); ++i)
    {
      particleHScond->AddComponent(particleHSComponents[i]);
    }

    condlist.push_back(particleHScond);
  /*--------------------------------------------------------------------*/
  // particle periodic boundary condition

  std::vector<Teuchos::RCP<ConditionComponent> > particlepbccomponents;
  // two vertices describing the bounding box for the pbc
  particlepbccomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ONOFF")));
  particlepbccomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("ONOFF",3)));

  Teuchos::RCP<ConditionDefinition> particlepbccond =
    Teuchos::rcp(new ConditionDefinition("DESIGN PARTICLE PERIODIC BOUNDARY CONDITION",
                                         "ParticlePeriodic",
                                         "Particle Periodic Boundary Condition",
                                         DRT::Condition::ParticlePeriodic,
                                         false,
                                         DRT::Condition::Particle));


  for (unsigned i=0; i<particlepbccomponents.size(); ++i)
  {
    particlepbccond->AddComponent(particlepbccomponents[i]);
  }

  condlist.push_back(particlepbccond);

  /*--------------------------------------------------------------------*/
  // particle init radius condition

  std::vector<Teuchos::RCP<ConditionComponent> > particleinitradiuscomponents;

  particleinitradiuscomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("SCALAR")));
  particleinitradiuscomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("SCALAR",1)));

  particleinitradiuscomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("FUNCT")));
  particleinitradiuscomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("FUNCT",1)));

  Teuchos::RCP<ConditionDefinition> particleradiuscond =
    Teuchos::rcp(new ConditionDefinition("DESIGN PARTICLE INIT RADIUS CONDITIONS",
                                         "InitialParticleRadius",
                                         "Particle Initial Radius Condition",
                                         DRT::Condition::ParticleInitRadius,
                                         false,
                                         DRT::Condition::Particle));


  for (unsigned i=0; i<particleinitradiuscomponents.size(); ++i)
  {
    particleradiuscond->AddComponent(particleinitradiuscomponents[i]);
  }

  condlist.push_back(particleradiuscond);

  /*--------------------------------------------------------------------*/
  // particle wall condition

  std::vector<Teuchos::RCP<ConditionComponent> > particlewallcomponents;
  particlewallcomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

  Teuchos::RCP<ConditionDefinition> surfpartwall =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE PARTICLE WALL",
                                         "ParticleWall",
                                         "Wall for particle collisions",
                                         DRT::Condition::ParticleWall,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<particlewallcomponents.size(); ++i)
  {
    surfpartwall->AddComponent(particlewallcomponents[i]);
  }

  condlist.push_back(surfpartwall);

}


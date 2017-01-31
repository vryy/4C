/*----------------------------------------------------------------------*/
/*!
\file inpar_particle.cpp

\brief Input parameters for particle problems

\level 1

\maintainer Georg Hammerl
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
                                  "RK2",
                                  "RK4",
                                  "HybridMeshFreeDivFree",
                                  "GenAlpha"),
                                tuple<int>(
                                  INPAR::PARTICLE::dyna_expleuler,
                                  INPAR::PARTICLE::dyna_centrdiff,
                                  INPAR::PARTICLE::dyna_rk2,
                                  INPAR::PARTICLE::dyna_rk4,
                                  INPAR::PARTICLE::dyna_hybridMeshFreeDivFree,
                                  INPAR::PARTICLE::dyna_genAlpha
                                ),
                                &particledyn);

   setStringToIntegralParameter<int>("TIMESTEPTYPE","Auto_CFL",
                                "type of time step choice",
                                tuple<std::string>(
                                  "Manual",
                                  "Auto_CFL"),
                                tuple<int>(
                                  INPAR::PARTICLE::Manual,
                                  INPAR::PARTICLE::Auto_CFL
                                ),
                                &particledyn);

   // Output type
   IntParameter("RESULTSEVRY",1,"save displacements and contact forces every RESULTSEVRY steps",&particledyn);
   IntParameter("RESEVRYERGY",0,"write system energies every requested step",&particledyn);
   IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&particledyn);
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
                               "WEIGHT_FUNCTION","CubicBspline",
                               "weight function for meshFree interaction dynamics",
                               tuple<std::string>(
                                 "CubicBspline",
                                 "SqrtHyperbola",
                                 "HyperbolaNoRsz"
                                 ),
                               tuple<int>(
                                 INPAR::PARTICLE::CubicBspline,
                                 INPAR::PARTICLE::SqrtHyperbola,
                                 INPAR::PARTICLE::HyperbolaNoRsz
                                 ),
                               &particledyn);

   setStringToIntegralParameter<int>(
                               "WALL_INTERACTION_TYPE","InitParticle",
                               "wall interaction type for MeshFree interactions",
                               tuple<std::string>(
                                 "InitParticle",
                                 "Mirror",
                                 "Custom"
                                 ),
                               tuple<int>(
                                 INPAR::PARTICLE::InitParticle,
                                 INPAR::PARTICLE::Mirror,
                                 INPAR::PARTICLE::Custom
                                 ),
                               &particledyn);

   DoubleParameter("ERROR_TOLL",1e-6,"tolerance of the error for implicit schemes",&particledyn);
   IntParameter("ITER_MAX",10,"maximum iteration per time step for implicit schemes",&particledyn);
   DoubleParameter("ALPHA_MIN",1e-6,"cap of the alpha parameter for the divergence free scheme for meshfree",&particledyn);
   DoubleParameter("CORRECT_DIVERGENCE_TOLL",1e-6,"tolerance for the iterative divergence corrector, divFree integration scheme, meshFree interaction",&particledyn);
   IntParameter("CORRECT_DIVERGENCE_ITER",10,"iterations for the iterative divergence corrector, divFree integration scheme, meshFree interaction",&particledyn);
   DoubleParameter("CORRECT_DENSITY_TOLL",1e-6,"tolerance for the iterative density corrector, divFree integration scheme, meshFree interaction",&particledyn);
   IntParameter("CORRECT_DENSITY_ITER",10,"iterations for the iterative density corrector, divFree integration scheme, meshFree interaction",&particledyn);
   DoubleParameter("WALL_FAKE_DENSITY",-1.0,"fake density for meshfree dynamics, in case of -1 the coefficients are extracted from the initial values of the particle material parameters",&particledyn);
   DoubleParameter("WALL_FAKE_MASS",-1.0,"fake mass of the wall element for meshfree dynamics, in case of -1 the coefficients are extracted from the initial values of the particle material parameters",&particledyn);
   DoubleParameter("WALL_FAKE_PRESSURE",-1.0,"fake pressure of the wall element for meshfree dynamics, in case of -1 the coefficients are extracted from the initial values of the particle material parameters",&particledyn);
   DoubleParameter("MIN_RADIUS",-1.0,"smallest particle radius",&particledyn);
   DoubleParameter("MIN_RADIUS",-1.0,"smallest particle radius",&particledyn);
   DoubleParameter("MAX_RADIUS",-1.0,"largest particle radius",&particledyn);
   DoubleParameter("REL_PENETRATION",-1.0,"relative particle penetration",&particledyn);
   DoubleParameter("MAX_VELOCITY",-1.0,"highest particle velocity",&particledyn);
   DoubleParameter("COEFF_RESTITUTION",-1.0,"coefficient of restitution",&particledyn);
   DoubleParameter("COEFF_RESTITUTION_WALL",-1.0,"coefficient of restitution (wall)",&particledyn);
   DoubleParameter("FRICT_COEFF_WALL",-1.0,"friction coefficient for contact particle-wall",&particledyn);
   DoubleParameter("FRICT_COEFF",-1.0,"dynamic friction coefficient for contact particle-particle",&particledyn);
   DoubleParameter("NORMAL_STIFF",-1.0,"stiffness for normal contact force",&particledyn);
   DoubleParameter("TANG_STIFF",-1.0,"stiffness for tangential contact force",&particledyn);
   DoubleParameter("NORMAL_DAMP",-1.0,"damping coefficient for normal contact force",&particledyn);
   DoubleParameter("TANG_DAMP",-1.0,"damping coefficient for tangential contact force",&particledyn);
   DoubleParameter("NORMAL_DAMP_WALL",-1.0,"damping coefficient for normal contact force (wall)",&particledyn);
   DoubleParameter("TANG_DAMP_WALL",-1.0,"damping coefficient for tangential contact force (wall)",&particledyn);
   BoolParameter("TENSION_CUTOFF","no","switch on/off tension cutoff",&particledyn);
   BoolParameter("MOVING_WALLS","no","switch on/off moving walls",&particledyn);
   IntParameter("TRANSFER_EVERY",1,"transfer particles every TRANSFER_EVERY steps",&particledyn);
   DoubleParameter("RANDOM_AMPLITUDE",0.0,"random value for initial position",&particledyn);
   BoolParameter("RADIUS_DISTRIBUTION","no","switch on/off random normal distribution of particle radii",&particledyn);
   DoubleParameter("RADIUS_DISTRIBUTION_SIGMA",-1.0,"standard deviation of normal distribution of particle radii",&particledyn);
   BoolParameter("RENDERING","no","switch on/off the rendering domain. If it is yes... you better have a RENDERING DOMAIN available",&particledyn);
   setNumericStringParameter("GRAVITY_ACCELERATION","0.0 0.0 0.0",
                             "Acceleration due to gravity in particle simulations.",
                             &particledyn);

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


   DoubleParameter("BETA",0.25,"Generalised-alpha factor in (0,1/2]",&genalpha);
   DoubleParameter("GAMMA",0.5,"Generalised-alpha factor in (0,1]",&genalpha);
   DoubleParameter("ALPHA_M",0.5,"Generalised-alpha factor in [0,1)",&genalpha);
   DoubleParameter("ALPHA_F",0.5,"Generalised-alpha factor in [0,1)",&genalpha);
   DoubleParameter("RHO_INF",-1.0,"Generalised-alpha factor in [0,1]",&genalpha);
   DoubleParameter("TOL",1e-6,"Generalised-alpha tollerance",&genalpha);
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
  particleinflowcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("inflow_vel_curve")));
  particleinflowcomponents.push_back(Teuchos::rcp(new IntConditionComponent("inflow_vel_curve")));
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


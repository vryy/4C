/*----------------------------------------------------------------------*/
/*!
\file particle_timint.cpp

\brief Time integration for particle dynamics

\level 1

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint.H"
#include "particle_algorithm.H"
#include "particle_contact.H"
#include "particle_resulttest.H"
#include "particle_utils.H"
#include "particle_timint_strategy.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/extparticle_mat.H"
#include "../drt_mat/matpar_bundle.H"
#include "particle_sph_interaction.H"
#include "particle_sph_rendering.H"

/*----------------------------------------------------------------------*/
/* print particle time logo */
void PARTICLE::TimInt::Logo()
{
 IO::cout << "Welcome to Particle Time Integration " <<IO::endl;
 IO::cout << "    ---                      ---     " <<IO::endl;
 IO::cout << "  /     \\                  /     \\   " <<IO::endl;
 IO::cout << "  |     |   ---->  <----   |     |   " <<IO::endl;
 IO::cout << "  \\     /                  \\     /   " <<IO::endl;
 IO::cout << "    ---                      ---     " <<IO::endl;
 IO::cout <<IO::endl;
}

/*----------------------------------------------------------------------*/
/* constructor */
PARTICLE::TimInt::TimInt
(
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& particledynparams,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization> actdis,
  Teuchos::RCP<IO::DiscretizationWriter> output
)
: discret_(actdis),
  myrank_(actdis->Comm().MyPID()),
  dbcmaps_(Teuchos::null),
  dbcdofs_(Teuchos::null),
  strategy_(Teuchos::null),
  output_(output),
  printlogo_(true),
  printscreen_(ioparams.get<int>("STDOUTEVRY")),
  errfile_(xparams.get<FILE*>("err file")),
  printerrfile_(errfile_),
  writerestartevery_(particledynparams.get<int>("RESTARTEVRY")),
  writestate_((bool) DRT::INPUT::IntegralValue<int>(ioparams,"STRUCT_DISP")),
  writevelacc_((bool) DRT::INPUT::IntegralValue<int>(ioparams,"STRUCT_VEL_ACC")),
  writeresultsevery_(particledynparams.get<int>("RESULTSEVRY")),
  writeenergyevery_(particledynparams.get<int>("RESEVRYERGY")),
  writerenderingevery_(particledynparams.get<int>("RESEVRYREND")),
  avrgrenderingsteps_(particledynparams.get<int>("AVRG_REND_STEPS")),
  writeparticlestatsevery_(particledynparams.get<int>("PARTICLESTATSEVRY")),
  energyfile_(Teuchos::null),
  particlestatsfile_(Teuchos::null),
  writeorientation_(false),
  kinergy_(0),
  intergy_(0),
  extergy_(0),
  linmomentum_(LINALG::Matrix<3,1>(true)),
  time_(Teuchos::null),
  timen_(0.0),
  dt_(Teuchos::null),
  timemax_(particledynparams.get<double>("MAXTIME")),
  stepmax_(particledynparams.get<int>("NUMSTEP")),
  step_(0),
  stepn_(0),
  restart_(0),
  dis_(Teuchos::null),
  vel_(Teuchos::null),
  acc_(Teuchos::null),
  velmod_(Teuchos::null),
  accmod_(Teuchos::null),
  angVel_(Teuchos::null),
  angAcc_(Teuchos::null),
  radius_(Teuchos::null),
  density_(Teuchos::null),
  densityDot_(Teuchos::null),

  disn_(Teuchos::null),
  veln_(Teuchos::null),
  accn_(Teuchos::null),
  velmodn_(Teuchos::null),
  accmodn_(Teuchos::null),
  angVeln_(Teuchos::null),
  angAccn_(Teuchos::null),
  radiusn_(Teuchos::null),
  densityn_(Teuchos::null),
  densityDotn_(Teuchos::null),

  fifc_(Teuchos::null),
  orient_(Teuchos::null),

  radius0_(Teuchos::null),
  radiusDot_(Teuchos::null),
  mass_(Teuchos::null),
  inertia_(Teuchos::null),
  pressure_(Teuchos::null),
  f_structure_(Teuchos::null),

  colorField_(Teuchos::null),
  colorFieldGrad_(Teuchos::null),
  smoothedColorFieldGrad_(Teuchos::null),
  accSF_(Teuchos::null),
  fspType_(Teuchos::null),
  curvature_(Teuchos::null),
  phaseColor_(Teuchos::null),

  global_num_boundaryparticles_(0),

  dofmapexporter_(Teuchos::null),
  nodemapexporter_(Teuchos::null),

  radiusdistribution_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::RadiusDistribution>(particledynparams,"RADIUS_DISTRIBUTION")),
  variableradius_((bool)DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->CavitationParams(),"COMPUTE_RADIUS_RP_BASED")),
  radiuschangefunct_(particledynparams.get<int>("RADIUS_CHANGE_FUNCT")),
  particle_algorithm_(Teuchos::null),
  collhandler_(Teuchos::null),
  interHandler_(Teuchos::null)
{
  // welcome user
  if ( (printlogo_) and (myrank_ == 0) )
  {
    Logo();
  }

  // check whether discretisation has been completed
  if (not discret_->Filled() || not actdis->HaveDofs())
  {
    dserror("Discretisation is not complete or has no dofs!");
  }

  // time state
  time_ = Teuchos::rcp(new TIMINT::TimIntMStep<double>(0, 0, 0.0));  // HERE SHOULD BE SOMETHING LIKE (particledynparams.get<double>("TIMEINIT"))
  dt_ = Teuchos::rcp(new TIMINT::TimIntMStep<double>(0, 0, particledynparams.get<double>("TIMESTEP")));
  step_ = 0;
  timen_ = (*time_)[0] + (*dt_)[0];  // set target time to initial time plus step size
  stepn_ = step_ + 1;

  return;
}

/*----------------------------------------------------------------------*/
/* initialization of time integration */
void PARTICLE::TimInt::Init()
{
  // initialize time integration strategy
  if(DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat_ellipsoids) >= 0)
    strategy_ = Teuchos::rcp(new TimIntStrategyEllipsoids(this));
  else
    strategy_ = Teuchos::rcp(new TimIntStrategySpheres(this));

  // initialize the vectors
  dis_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
  vel_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
  acc_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
  fifc_ = LINALG::CreateVector(*DofRowMapView(), true);
  mass_ = LINALG::CreateVector(*NodeRowMapView(), true);
  if(writeorientation_)
    orient_ = LINALG::CreateVector(*DofRowMapView());

  switch (particle_algorithm_->ParticleInteractionType())
  {
  case INPAR::PARTICLE::SPH :
  {
    density_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, NodeRowMapView(), true));
    densityDot_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, NodeRowMapView(), true));
    pressure_ = LINALG::CreateVector(*NodeRowMapView(), true);

#ifdef PARTICLE_WRITECOLORFIELD
    colorField_ = LINALG::CreateVector(*NodeRowMapView(), true);
    colorFieldGrad_ = LINALG::CreateVector(*DofRowMapView(), true);
    smoothedColorFieldGrad_ = LINALG::CreateVector(*DofRowMapView(), true);
    accSF_ = LINALG::CreateVector(*DofRowMapView(), true);
    fspType_ = LINALG::CreateVector(*NodeRowMapView(), true);
    curvature_ = LINALG::CreateVector(*NodeRowMapView(), true);
    phaseColor_ = LINALG::CreateVector(*NodeRowMapView(), true);
#endif

    break;
  }
  default : //do nothing
    break;
  }

  if(variableradius_ or radiuschangefunct_ > 0)
  {
    // initial radius of each particle for time dependent radius
    radius0_  = LINALG::CreateVector(*discret_->NodeRowMap(), true);
    // time derivative of radius of each particle for time dependent radius
    radiusDot_  = LINALG::CreateVector(*NodeRowMapView(), true);
  }

  // Apply Dirichlet BC and create dbc map extractor
  {
    dbcdofs_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
    dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
    Teuchos::ParameterList p;
    p.set("total time", (*time_)[0]);
    discret_->EvaluateDirichlet(p, (*dis_)(0), (*vel_)(0), (*acc_)(0), dbcdofs_, dbcmaps_);
  }

  // determine boundary particles
  DetermineBdryParticles();

  // set initial fields
  SetInitialFields();

  // copy everything into the n+1 state vectors
  disn_ = Teuchos::rcp(new Epetra_Vector(*(*dis_)(0)));
  veln_ = Teuchos::rcp(new Epetra_Vector(*(*vel_)(0)));
  accn_ = Teuchos::rcp(new Epetra_Vector(*(*acc_)(0)));

  if (DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"TRANSPORT_VELOCITY")==true)
  {
    velmod_ =  Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
    velmodn_ = Teuchos::rcp(new Epetra_Vector(*(*velmod_)(0)));
    accmod_ =  Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
    accmodn_ = Teuchos::rcp(new Epetra_Vector(*(*accmod_)(0)));
  }

  switch (particle_algorithm_->ParticleInteractionType())
  {
  case INPAR::PARTICLE::SPH :
  {
    radiusn_  = Teuchos::rcp(new Epetra_Vector(*(*radius_)(0)));
    densityn_ = Teuchos::rcp(new Epetra_Vector(*(*density_)(0)));
    densityDotn_ = Teuchos::rcp(new Epetra_Vector(*(*densityDot_)(0)));
    break;
  }
  default : //do nothing
    break;
  }

  // decide whether there is particle contact
  if(particle_algorithm_->ParticleInteractionType() != INPAR::PARTICLE::None)
  {

    // allocate vectors
    angVel_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
    angAcc_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));

    // copy the vectors to the (n+1) state vectors
    angVeln_ = LINALG::CreateVector(*DofRowMapView(),true);
    angAccn_ = LINALG::CreateVector(*DofRowMapView(),true);

    // create and fill inertia
    strategy_->ComputeInertia();
  }

  // output file for energy
  if(writeenergyevery_ != 0 and myrank_ == 0)
    AttachEnergyFile();

  // output file for particle statistics
  if(writeparticlestatsevery_ > 0 and myrank_ == 0)
    AttachParticleStatisticsFile();

  return;
}

/*----------------------------------------------------------------------*/
/* Set initial fields */
void PARTICLE::TimInt::SetInitialFields()
{
  // -----------------------------------------//
  // set material properties
  // -----------------------------------------//
  // set initial particle radii
  strategy_->SetInitialRadii();

  // extract particle density (for now, all particles have identical density)
  const MAT::PAR::ParticleMat* const particlemat = particle_algorithm_->ParticleMat();

  if(particle_algorithm_->ParticleInteractionType()!=INPAR::PARTICLE::SPH)
  {
    const double initDensity = particlemat != NULL ? particlemat->initDensity_ : 0.;
    double consistent_problem_volume=DRT::Problem::Instance()->ParticleParams().get<double>("CONSISTENT_PROBLEM_VOLUME");
    if(consistent_problem_volume<0.0)
      strategy_->ComputeMass();
    else
      dserror("The definition of a CONSISTENT_PROBLEM_VOLUME is only possible for SPH applications!");

    // -----------------------------------------//
    // set initial radius condition if existing
    // -----------------------------------------//

    std::vector<DRT::Condition*> condition;
    discret_->GetCondition("InitialParticleRadius", condition);

    if(consistent_problem_volume>0.0 and condition.size()>0)
      dserror("The combination of InitialParticleRadius and CONSISTENT_PROBLEM_VOLUME not possible so far!");

    // loop over conditions
    for (size_t i=0; i<condition.size(); ++i)
    {
      double scalar  = condition[i]->GetDouble("SCALAR");
      int funct_num  = condition[i]->GetInt("FUNCT");

      const std::vector<int>* nodeids = condition[i]->Nodes();
      //loop over particles in current condition
      for(size_t counter=0; counter<(*nodeids).size(); ++counter)
      {
        int lid = discret_->NodeRowMap()->LID((*nodeids)[counter]);
        if(lid != -1)
        {
          DRT::Node *currparticle = discret_->gNode((*nodeids)[counter]);
          double function_value =  DRT::Problem::Instance()->Funct(funct_num-1).Evaluate(0, currparticle->X(),0.0);
          double r_p = (*(*radius_)(0))[lid];
          r_p *= function_value * scalar;
          (*(*radius_)(0))[lid] = r_p;
          if(r_p <= 0.0)
            dserror("negative initial radius");

          // mass-vector: m = rho * 4/3 * PI * r^3
          (*mass_)[lid] = initDensity * PARTICLE::Utils::Radius2Volume(r_p);
        }
      }
    }

    // -----------------------------------------//
    // evaluate random normal distribution for particle radii if applicable
    // -----------------------------------------//

    switch(radiusdistribution_)
    {
      case INPAR::PARTICLE::radiusdistribution_none:
      {
        // do nothing
        break;
      }
      case INPAR::PARTICLE::radiusdistribution_lognormal:
      case INPAR::PARTICLE::radiusdistribution_normal:
      {
        if(consistent_problem_volume>0.0)
          dserror("The combination of RADIUS_DISTRIBUTION and CONSISTENT_PROBLEM_VOLUME not possible so far!");

        // get minimum and maximum radius for particles
        const double min_radius = DRT::Problem::Instance()->ParticleParams().get<double>("MIN_RADIUS");
        const double max_radius = DRT::Problem::Instance()->ParticleParams().get<double>("MAX_RADIUS");

        // loop over all particles
        for(int n=0; n<discret_->NumMyRowNodes(); ++n)
        {
          // get local ID of current particle
          const int lid = discret_->NodeRowMap()->LID(discret_->lRowNode(n)->Id());

          // initialize random number generator with current particle radius or its natural logarithm as mean and input parameter value as standard deviation
          DRT::Problem::Instance()->Random()->SetMeanVariance(radiusdistribution_ == INPAR::PARTICLE::radiusdistribution_lognormal ? log((*(*radius_)(0))[lid]) : (*(*radius_)(0))[lid],DRT::Problem::Instance()->ParticleParams().get<double>("RADIUS_DISTRIBUTION_SIGMA"));

          // generate normally or log-normally distributed random value for particle radius
          double random_radius = radiusdistribution_ == INPAR::PARTICLE::radiusdistribution_lognormal ? exp(DRT::Problem::Instance()->Random()->Normal()) : DRT::Problem::Instance()->Random()->Normal();

          // check whether random value lies within allowed bounds, and adjust otherwise
          if(random_radius > max_radius)
            random_radius = max_radius;
          else if(random_radius < min_radius)
            random_radius = min_radius;

          // set particle radius to random value
          (*(*radius_)(0))[lid] = random_radius;

          // recompute particle mass
          (*mass_)[lid] = initDensity * PARTICLE::Utils::Radius2Volume(random_radius);
        }

        break;
      }

      default:
      {
        dserror("Invalid random distribution of particle radii!");
        break;
      }
    }
  }

  // -----------------------------------------//
  // initialize displacement field
  // -----------------------------------------//

  const double amplitude = DRT::Problem::Instance()->ParticleParams().get<double>("RANDOM_AMPLITUDE");
  const double initRadius = particlemat != NULL ? particlemat->initRadius_ : 0.;

  for(int n=0; n<discret_->NumMyRowNodes(); ++n)
  {
    DRT::Node* actnode = discret_->lRowNode(n);
    // get the first gid of a node and convert it into a LID
    int gid = discret_->Dof(actnode, 0);
    int lid = discret_->DofRowMap()->LID(gid);
    for (int dim=0; dim<3; ++dim)
    {
      if(amplitude)
      {
        double randomValue = DRT::Problem::Instance()->Random()->Uni();
        (*(*dis_)(0))[lid+dim] = actnode->X()[dim] + randomValue * amplitude * initRadius;
      }
      else
      {
        (*(*dis_)(0))[lid+dim] = actnode->X()[dim];
      }
    }
  }

  // -----------------------------------------//
  // initialize orientation field
  // -----------------------------------------//
  if(writeorientation_)
    strategy_->SetInitialOrientation();

  // -----------------------------------------//
  // set initial velocity field if existing
  // -----------------------------------------//

  const std::string field = "Velocity";
  std::vector<int> localdofs;
  localdofs.push_back(0);
  localdofs.push_back(1);
  localdofs.push_back(2);
  discret_->EvaluateInitialField(field,(*vel_)(0),localdofs);

  // set vector of initial particle radii if necessary
  if(radiuschangefunct_ > 0)
    radius0_->Update(1.,*(*radius_)(0),0.);

  // In case of SPH, the material might be different from particle to particle. Thus, initial mass, density etc. has to be set individually
  if(particle_algorithm_->ParticleInteractionType()==INPAR::PARTICLE::SPH)
  {
    Teuchos::RCP<PARTICLE::ParticleSPHInteractionHandler> interHandler= Teuchos::rcp(new PARTICLE::ParticleSPHInteractionHandler(discret_, particle_algorithm_, DRT::Problem::Instance()->ParticleParams(),true));
    interHandler->InitColParticles();

    double consistent_problem_volume=DRT::Problem::Instance()->ParticleParams().get<double>("CONSISTENT_PROBLEM_VOLUME");
    if(consistent_problem_volume>0.0)
    {
      //boundary particles are excluded for determination of particle mass since they are also not part of the consistent_problem_volume
      const int num_particles = discret_->NumGlobalNodes() - global_num_boundaryparticles_;
      const double particle_volume = consistent_problem_volume / ((double)(num_particles));

      interHandler->InitDensityAndMass(particle_volume,(*density_)(0),mass_);
    }
    else
      dserror("CONSISTENT_PROBLEM_VOLUME problem volume required for SPH simulations!");
  }

  return;
}

/*----------------------------------------------------------------------*/
/* Determine boundary particles */
void PARTICLE::TimInt::DetermineBdryParticles()
{
  const INPAR::PARTICLE::WallInteractionType wallInteractionType = DRT::INPUT::IntegralValue<INPAR::PARTICLE::WallInteractionType>(DRT::Problem::Instance()->ParticleParams(),"WALL_INTERACTION_TYPE");
  if( wallInteractionType == INPAR::PARTICLE::BoundarParticle_FreeSlip or wallInteractionType == INPAR::PARTICLE::BoundarParticle_NoSlip )
  {
    // local number of boundary particles on this processor
    int local_num_boundaryparticles = 0;

    // export dbcdofs (which is in row map format) into overlapping column map format
    Epetra_Vector dbcdofs_col(*(discret_->DofColMap()),true);
    LINALG::Export(*dbcdofs_,dbcdofs_col);

    // loop over all row nodes
    for (int lidNodeCol=0; lidNodeCol < discret_->NodeColMap()->NumMyElements(); ++lidNodeCol)
    {
      DRT::Node* particle_i = discret_->lColNode(lidNodeCol);

      std::vector<int> lm;
      lm.reserve(3);
      discret_->Dof(particle_i, lm);

      if ( dbcdofs_col[discret_->DofColMap()->LID(lm[0])] == 1 )
      {
        if ( dbcdofs_col[discret_->DofColMap()->LID(lm[1])] != 1 or dbcdofs_col[discret_->DofColMap()->LID(lm[2])] != 1 )
          dserror("For boundary particles all three DoFs have to be prescribed by a Dirichlet Condition!");

        // found a boundary particle
        PARTICLE::ParticleNode* particleNode_i = dynamic_cast<PARTICLE::ParticleNode*>(particle_i);
        if (particle_i == NULL)
          dserror("Dynamic cast to ParticleNode failed");
        particleNode_i->Set_bdry_particle(true);

        // boundary particle is in row map of this processor
        if ( not (discret_->NodeRowMap()->LID(particle_i->Id()) < 0) )
          local_num_boundaryparticles += 1;
      }
    }

    // global number of boundary particles on all processors
    discret_->Comm().SumAll(&local_num_boundaryparticles, &global_num_boundaryparticles_, 1);
  }

  return;
}

/*----------------------------------------------------------------------*/
/* prepare time step and apply Dirichlet boundary conditions */
void PARTICLE::TimInt::PrepareTimeStep()
{
  // Update map containing Dirichlet DOFs if existing
  if(dbcmaps_ != Teuchos::null && dbcmaps_->CondMap()->NumGlobalElements() != 0)
  {
    // apply Dirichlet BC and rebuild map extractor
    ApplyDirichletBC(timen_, disn_, veln_, accn_, true);
  }

  // update particle radii if necessary
  if(radiuschangefunct_ > 0)
    (*radius_)(0)->Update(DRT::Problem::Instance()->Funct(radiuschangefunct_-1).EvaluateTime(timen_),*radius0_,0.);

  return;
}

/*----------------------------------------------------------------------*/
/* equilibrate system at initial state and identify consistent accelerations */
void PARTICLE::TimInt::DetermineMassDampConsistAccel()
{
  if(particle_algorithm_->ParticleInteractionType()!=INPAR::PARTICLE::SPH)
    ComputeAcc(Teuchos::null, Teuchos::null, (*acc_)(0), Teuchos::null);
  else
  {
    INPAR::PARTICLE::DynamicType timinttype = DRT::INPUT::IntegralValue<INPAR::PARTICLE::DynamicType>(DRT::Problem::Instance()->ParticleParams(),"DYNAMICTYP");
    if(timinttype==INPAR::PARTICLE::dyna_kickdrift)
    {
      if (DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"TRANSPORT_VELOCITY")==true)
      {
        //Initially, the convection velocity velmodn_ is set equal to the initial velocity veln_
        velmodn_->Update(1.0,*veln_,0.0);
        //Initialize veln_ and velmodn_ since required for (*acc_)(0) and (*accmod_)(0)
        ApplyDirichletBC(0.0, Teuchos::null, veln_, Teuchos::null, false);
        ApplyDirichletBC(0.0, Teuchos::null, velmodn_, Teuchos::null, false);

        DetermineSPHDensAndAcc((*acc_)(0),(*accmod_)(0),velmodn_,Teuchos::null,0.0,0.0);
      }
      else
        DetermineSPHDensAndAcc((*acc_)(0),Teuchos::null,Teuchos::null,Teuchos::null,0.0,0.0);
    }
    else
      dserror("Currently, on the kick-drift time integrator is suitable for SPH!");
  }

  return;
}

/*----------------------------------------------------------------------*/
/* acceleration is applied from given forces */
void PARTICLE::TimInt::ComputeAcc(
  Teuchos::RCP<Epetra_Vector> f_contact,
  Teuchos::RCP<Epetra_Vector> m_contact,
  Teuchos::RCP<Epetra_Vector> global_acc,
  Teuchos::RCP<Epetra_Vector> global_angAcc)
{
  int numrownodes = discret_->NodeRowMap()->NumMyElements();

  // in case of contact, consider corresponding forces and moments
  if(f_contact != Teuchos::null)
  {
    // sum all forces (contact and external)
    fifc_->Update(1.0, *f_contact, 1.0);

    // zero out non-planar entries in case of 2D
    if(particle_algorithm_->BinStrategy()->ParticleDim() == INPAR::PARTICLE::particle_2Dz)
    {
      for(int i=0; i<numrownodes; ++i)
      {
        (*m_contact)[i*3+0] = 0.0;
        (*m_contact)[i*3+1] = 0.0;
      }
    }

    // compute angular acceleration
    if(global_angAcc != Teuchos::null)
      strategy_->ComputeAngularAcceleration(*global_angAcc,*m_contact);
  }

  // zero out non-planar entries in case of 2D
  if(particle_algorithm_->BinStrategy()->ParticleDim() == INPAR::PARTICLE::particle_2Dz)
  {
    for(int i=0; i<numrownodes; ++i)
      (*fifc_)[i*3+2] = 0.0;
  }

  // update of translational acceleration
  for(int i=0; i<numrownodes; ++i)
  {
    const double invmass = 1.0/(*mass_)[i];
    for(int dim=0; dim<3; ++dim)
      (*global_acc)[i*3+dim] = invmass * (*fifc_)[i*3+dim];
  }

  return;
}

/*----------------------------------------------------------------------*/
/* Determine acceleration */
void PARTICLE::TimInt::DetermineSPHDensAndAcc(Teuchos::RCP<Epetra_Vector> acc,
                                            Teuchos::RCP<Epetra_Vector> accmod,
                                            Teuchos::RCP<Epetra_Vector> velConv,
                                            Teuchos::RCP<Epetra_Vector> acc_A,
                                            const double time,
                                            const double dt)
{

  //Initialize all columns and boundary particles, set sate vectors, search for neighbor particles.
  interHandler_->Init(disn_, veln_, radiusn_, mass_);
  //Set also state vector velConv
  if(velConv!=Teuchos::null)
    interHandler_->SetStateVector(velConv, PARTICLE::VelConv);

  const INPAR::PARTICLE::FreeSurfaceType freeSurfaceType=DRT::INPUT::IntegralValue<INPAR::PARTICLE::FreeSurfaceType>(DRT::Problem::Instance()->ParticleParams(),"FREE_SURFACE_TYPE");
  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"DENSITY_SUMMATION")==false)
  {
    //Density update via continuity equation: This is only required as initialization for density summation below in case of free-surface flows
    //Attention: Here, r_{n+1} and v_{n+1/2} are applied for calculation of densityDotn_ instead of r_{n+1/2} and v_{n+1/2} as in Adami 2012 or
    //r_{n} and v_{n+1/2} as in Zhang 2017. This procedure is simpler here since otherwise determination of particle positions and neighbors would have
    //to be done twice. On the other hand, the error should be negligible since densityDotn_ is only required for density initialization before applying the summation formula.
    //TODO: Enable density evaluation at flexible points in time (e.g. at t_{n}, t_{n+1/2}, t_{n+1})

    interHandler_->Inter_pvp_densityDot(densityDotn_);
    densityn_->Update(dt, *densityDotn_, 1.0);
    interHandler_->SetStateVector(densityn_, PARTICLE::Density);

    #ifdef PARTICLE_BOUNDARYDENSITY
    interHandler_->Density2Pressure(densityn_,pressure_);
    interHandler_->SetStateVector(pressure_, PARTICLE::Pressure);
    //Asign accelerations, modified pressures and modified velocities for boundary particles and calculate their mechanical energy contribution
    const INPAR::PARTICLE::WallInteractionType wallInteractionType=DRT::INPUT::IntegralValue<INPAR::PARTICLE::WallInteractionType>(DRT::Problem::Instance()->ParticleParams(),"WALL_INTERACTION_TYPE");
    if(wallInteractionType==INPAR::PARTICLE::BoundarParticle_NoSlip or wallInteractionType==INPAR::PARTICLE::BoundarParticle_FreeSlip)
    {
      interHandler_->InitBoundaryData(acc,particle_algorithm_->GetGravityAcc(time));
    }
    #endif
  }
  else
  {
    if(freeSurfaceType!=INPAR::PARTICLE::FreeSurface_None and freeSurfaceType!=INPAR::PARTICLE::TwoPhase)
      dserror("Density summation not suitable for free-surface problems. Choose input parameter DENSITY_SUMMATION=No!");

    interHandler_->MF_mW(densityn_);
    interHandler_->SetStateVector(densityn_, PARTICLE::Density);
  }

  if(freeSurfaceType!=INPAR::PARTICLE::FreeSurface_None)
  {
    //In case of density summation, the new density and new pressure have been determined as very first step since they are required for all the following calculations
    Teuchos::RCP<Epetra_Vector> density_sum = Teuchos::rcp(new Epetra_Vector(*densityn_));
    density_sum->PutScalar(0.0);
    Teuchos::RCP<Epetra_Vector> colorField = Teuchos::rcp(new Epetra_Vector(*densityn_));
    colorField->PutScalar(0.0);
    Teuchos::RCP<Epetra_Vector> colorFieldGrad = Teuchos::rcp(new Epetra_Vector(*veln_));
    colorFieldGrad->PutScalar(0.0);

    interHandler_->MF_mW(density_sum,colorField,colorFieldGrad,phaseColor_);
    interHandler_->SetStateVector(colorField, PARTICLE::ColorField);
    interHandler_->SetStateVector(density_sum, PARTICLE::DensitySum);
    interHandler_->SetStateVector(colorFieldGrad, PARTICLE::ColorFieldGrad);

    const INPAR::PARTICLE::SurfaceTensionType surfaceTensionType=DRT::INPUT::IntegralValue<INPAR::PARTICLE::SurfaceTensionType>(DRT::Problem::Instance()->ParticleParams(),"SURFACE_TENSION_TYPE");
    if(surfaceTensionType==INPAR::PARTICLE::ST_CONTI_ADAMI or surfaceTensionType==INPAR::PARTICLE::ST_CONTI_HU)
    {
      Teuchos::RCP<Epetra_Vector> smoothedColorFieldGrad = Teuchos::rcp(new Epetra_Vector(*veln_));
      smoothedColorFieldGrad->PutScalar(0.0);
      interHandler_->MF_SmoothedCFG(smoothedColorFieldGrad);
      interHandler_->SetStateVector(smoothedColorFieldGrad, PARTICLE::SmoothedColorFieldGrad);

      if(smoothedColorFieldGrad_!=Teuchos::null)
        smoothedColorFieldGrad_->Update(1.0,*smoothedColorFieldGrad,0.0);
    }

    //Determine free-surface particles
    interHandler_->InitFreeSurfaceParticles(fspType_);

    //Re-initialize density if required in case of free-surface flow
    if(freeSurfaceType!=INPAR::PARTICLE::DensityIntegration and freeSurfaceType!=INPAR::PARTICLE::TwoPhase)
      interHandler_->MF_ReInitDensity(densityn_,freeSurfaceType);

    if(colorField_!=Teuchos::null)
      colorField_->Update(1.0,*colorField,0.0);

    if(colorFieldGrad_!=Teuchos::null)
      colorFieldGrad_->Update(1.0,*colorFieldGrad,0.0);
  }

  // determine also the new pressure and set state vector
  interHandler_->Density2Pressure(densityn_,pressure_);
  interHandler_->SetStateVector(pressure_, PARTICLE::Pressure);

  // assign accelerations, modified pressures and modified velocities for boundary particles and calculate their mechanical energy contribution
  const INPAR::PARTICLE::WallInteractionType wallInteractionType=DRT::INPUT::IntegralValue<INPAR::PARTICLE::WallInteractionType>(DRT::Problem::Instance()->ParticleParams(),"WALL_INTERACTION_TYPE");
  if(wallInteractionType==INPAR::PARTICLE::BoundarParticle_NoSlip or wallInteractionType==INPAR::PARTICLE::BoundarParticle_FreeSlip)
    interHandler_->InitBoundaryData(acc,particle_algorithm_->GetGravityAcc(time));

  // clear acceleration states
  acc->PutScalar(0.0);
  if(accmod!=Teuchos::null)
    accmod->PutScalar(0.0);

  // acceleration contributions due to gravity forces
  GravityAcc(acc,time);

  // acceleration contributions due to internal forces (pressure, viscosity, etc.)
  interHandler_->Inter_pvp_acc(acc,accmod,acc_A,time);

  const INPAR::PARTICLE::SurfaceTensionType surfaceTensionType=DRT::INPUT::IntegralValue<INPAR::PARTICLE::SurfaceTensionType>(DRT::Problem::Instance()->ParticleParams(),"SURFACE_TENSION_TYPE");
  if(surfaceTensionType==INPAR::PARTICLE::ST_CONTI_ADAMI or surfaceTensionType==INPAR::PARTICLE::ST_CONTI_HU)
  {
    Teuchos::RCP<Epetra_Vector> kappa = Teuchos::rcp(new Epetra_Vector(*densityn_));
    kappa->PutScalar(0.0);
    Teuchos::RCP<Epetra_Vector> accSF = Teuchos::rcp(new Epetra_Vector(*acc));
    accSF->PutScalar(0.0);

    if(freeSurfaceType==INPAR::PARTICLE::TwoPhase)
    {
      if(surfaceTensionType==INPAR::PARTICLE::ST_CONTI_ADAMI)
        interHandler_->Inter_fspvp_Adami_1(accSF,kappa,time);

      if(surfaceTensionType==INPAR::PARTICLE::ST_CONTI_HU)
        dserror("This variant is not implemented for two-phase flow yet!");
    }
    else
    {
      if(surfaceTensionType==INPAR::PARTICLE::ST_CONTI_ADAMI)
        interHandler_->Inter_fspvp_Adami_2(accSF,kappa,time);

      if(surfaceTensionType==INPAR::PARTICLE::ST_CONTI_HU)
        interHandler_->Inter_fspvp_Hu(accSF,time);
    }


    acc->Update(1.0,*accSF,1.0);

    if(accSF_!=Teuchos::null)
      accSF_->Update(1.0,*accSF,0.0);

    if(curvature_!=Teuchos::null)
      curvature_->Update(1.0,*kappa,0.0);
  }

  intergy_ = 0.0;
  interHandler_->DetermineIntEnergy(intergy_);

  // clear vectors, keep memory
  interHandler_->Clear();
}

/*---------------------------------------------------------------*/
/* Apply Dirichlet boundary conditions on provided state vectors */
void PARTICLE::TimInt::ApplyDirichletBC
(
  const double time,
  Teuchos::RCP<Epetra_Vector> dis,
  Teuchos::RCP<Epetra_Vector> vel,
  Teuchos::RCP<Epetra_Vector> acc,
  bool recreatemap
)
{
  // needed parameters
  Teuchos::ParameterList p;
  p.set("total time", time);  // target time

  // predicted Dirichlet values
  // \c dis then also holds prescribed new Dirichlet displacements
  discret_->ClearState();
  if (recreatemap)
    discret_->EvaluateDirichlet(p, dis, vel, acc, Teuchos::null, dbcmaps_);
  else
    discret_->EvaluateDirichlet(p, dis, vel, acc, Teuchos::null, Teuchos::null);
  discret_->ClearState();

  return;
}

/*----------------------------------------------------------------------*/
/* Update time and step counter */
void PARTICLE::TimInt::UpdateStepTime()
{
  // update time and step
  time_->UpdateSteps(timen_);  // t_{n} := t_{n+1}, etc
  step_ = stepn_;  // n := n+1
  //
  timen_ += (*dt_)[0];
  stepn_ += 1;

  // new deal
  return;
}

/*----------------------------------------------------------------------*/
/* State vectors are updated according to the new distribution of particles */

void PARTICLE::TimInt::UpdateStatesAfterParticleTransfer()
{
  UpdateExportersIfNecessary(mass_->Map(),(*dis_)(0)->Map());

  UpdateStateVectorMap(dis_);
  UpdateStateVectorMap(vel_);
  UpdateStateVectorMap(acc_);
  UpdateStateVectorMap(velmod_);
  UpdateStateVectorMap(accmod_);
  UpdateStateVectorMap(angVel_);
  UpdateStateVectorMap(angAcc_);
  strategy_->UpdateRadiusVectorMap();
  UpdateStateVectorMap(density_,true);
  UpdateStateVectorMap(densityDot_,true);
  UpdateStateVectorMap(disn_);
  UpdateStateVectorMap(veln_);
  UpdateStateVectorMap(accn_);
  UpdateStateVectorMap(velmodn_);
  UpdateStateVectorMap(accmodn_);
  UpdateStateVectorMap(angVeln_);
  UpdateStateVectorMap(angAccn_);
  UpdateStateVectorMap(radiusn_,true);
  UpdateStateVectorMap(densityn_,true);
  UpdateStateVectorMap(densityDotn_,true);

  UpdateStateVectorMap(fifc_);
  UpdateStateVectorMap(orient_);

  UpdateStateVectorMap(radius0_,true);
  UpdateStateVectorMap(radiusDot_,true);
  UpdateStateVectorMap(mass_,true);
  strategy_->UpdateInertiaVectorMap();
  UpdateStateVectorMap(pressure_,true);

  if(colorField_!=Teuchos::null)
    UpdateStateVectorMap(colorField_,true);

  if(colorFieldGrad_!=Teuchos::null)
    UpdateStateVectorMap(colorFieldGrad_);

  if(smoothedColorFieldGrad_!=Teuchos::null)
    UpdateStateVectorMap(smoothedColorFieldGrad_);

  if(accSF_!=Teuchos::null)
    UpdateStateVectorMap(accSF_);

  if(fspType_!=Teuchos::null)
    UpdateStateVectorMap(fspType_,true);

  if(curvature_!=Teuchos::null)
    UpdateStateVectorMap(curvature_,true);

  if(phaseColor_!=Teuchos::null)
    UpdateStateVectorMap(phaseColor_,true);
}

/*----------------------------------------------------------------------*/
/* Read and set restart values */
void PARTICLE::TimInt::ReadRestart
(
  const int step
)
{
  IO::DiscretizationReader reader(discret_, step);
  if (step != reader.ReadInt("step"))
    dserror("Time step on file not equal to given step");

  restart_ = step;
  step_ = step;
  stepn_ = step_ + 1;
  time_ = Teuchos::rcp(new TIMINT::TimIntMStep<double>(0, 0, reader.ReadDouble("time")));
  timen_ = (*time_)[0] + (*dt_)[0];

  ReadRestartState();

  // short screen output
  if (myrank_ == 0)
    IO::cout << "====== Restart of the particle simulation from step "
        << step_ << IO::endl;

  return;
}

/*----------------------------------------------------------------------*/
/* Read and set restart state */
void PARTICLE::TimInt::ReadRestartState()
{
  IO::DiscretizationReader reader(discret_, step_);
  // maps need to be adapted to restarted discretization

  UpdateStatesAfterParticleTransfer();

  // start with reading displacement in order to find out whether particles exist
  reader.ReadVector(disn_, "displacement");

  // check, in case there is nothing, do not read the file
  if(disn_->GlobalLength() == 0)
    return;

  // now finish the displacement and read the remaining state vectors
  dis_->UpdateSteps(*disn_);
  reader.ReadVector(veln_, "velocity");
  vel_->UpdateSteps(*veln_);

#ifndef PARTICLE_NORESTARTACC
  reader.ReadVector(accn_, "acceleration");
  acc_->UpdateSteps(*accn_);
#endif

  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"TRANSPORT_VELOCITY")==true)
  {
    reader.ReadVector(velmodn_, "modified_velocity");
    velmod_->UpdateSteps(*velmodn_);
    reader.ReadVector(accmodn_, "modified_acceleration");
    accmod_->UpdateSteps(*accmodn_);
  }

  reader.ReadVector(mass_, "mass");


  if ( particle_algorithm_->ParticleInteractionType() != INPAR::PARTICLE::SPH &&
       DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat_ellipsoids) < 0)
  {
    // create a dummy vector to extract the radius vector (radiusn_ does not exist)
    Teuchos::RCP<Epetra_Vector> radius = LINALG::CreateVector(*discret_->NodeRowMap(), true);
    reader.ReadVector(radius, "radius");
    radius_->UpdateSteps(*radius);
  }

  if (particle_algorithm_->ParticleInteractionType() == INPAR::PARTICLE::SPH)
    interHandler_->Density2Pressure(densityn_,pressure_);

  // read in particle collision relevant data
  if(collhandler_ != Teuchos::null)
  {
    // initialize inertia
    strategy_->ComputeInertia();

    reader.ReadVector(angVeln_, "ang_velocity");
    angVel_->UpdateSteps(*angVeln_);
    reader.ReadVector(angAccn_, "ang_acceleration");
    angAcc_->UpdateSteps(*angAccn_);
    if(writeorientation_)
      reader.ReadVector(orient_, "orientation");
  }

  // read in variable radius relevant data
  if(variableradius_ == true)
  {
    reader.ReadVector(radius0_, "radius0");
    // time derivative of radius of each particle for time dependent radius
    reader.ReadVector(radiusDot_, "radiusDot");
  }
}

/*----------------------------------------------------------------------*/
/* Calculate all output quantities that depend on a potential material history */
void PARTICLE::TimInt::PrepareOutput()
{
  DetermineEnergy();
  return;
}

/*----------------------------------------------------------------------*/
/* output displacement */
void PARTICLE::TimInt::OutputDisplacement() const
{
  WriteVector("displacement", dis_);

  return;
}

/*----------------------------------------------------------------------*/
/* output to file */
void PARTICLE::TimInt::OutputStep(bool forced_writerestart)
{
  // this flag is passed along subroutines and prevents
  // repeated initialising of output writer, printing of
  // state vectors, or similar
  bool datawritten = false;

  // output restart (try this first)
  // write restart step
  if ( (writerestartevery_ and ((step_-restart_)%writerestartevery_ == 0)) or forced_writerestart )
  {
    OutputRestart(datawritten);
  }

  // output results (not necessary if restart in same step)
  if ( writestate_
       and writeresultsevery_ and ((step_-restart_)%writeresultsevery_ == 0)
       and (not datawritten) )
  {
    OutputState(datawritten);
  }

  // output energy
  if ( writeenergyevery_ and ((step_-restart_)%writeenergyevery_ == 0) )
  {
    OutputEnergy();
  }

  // output SPH rendering
  if (writerenderingevery_)
  {
    int renderOutputMod = (step_-restart_)%writerenderingevery_;
    int averageRenderingEvery=writerenderingevery_/avrgrenderingsteps_;
    int renderAverageMod = (step_-restart_)%averageRenderingEvery;
    if (renderAverageMod==0)
    {
      // If the last output has been writen at t_n, the next averaging steps are:
      // t_n+1*averageRenderingEvery, t_n+2*averageRenderingEvery, ..., t_n+avrgrenderingsteps_*averageRenderingEvery=t_n+writerenderingevery_.
      // At the first of these averaging steps, the rendering states have to be cleared initially.
      // If avrgrenderingsteps_==1, the vectors have to be cleared in every step!
      bool clearstate = ((renderOutputMod == averageRenderingEvery) or (avrgrenderingsteps_==1));

      // write rendering states at output step
      bool writeoutput = (renderOutputMod==0);

      PerformSPHRendering(clearstate, writeoutput);
    }
  }

  // output particle statistics
  if(writeparticlestatsevery_ and (step_-restart_) % writeparticlestatsevery_ == 0)
    OutputParticleStatistics();

  return;
}

/*----------------------------------------------------------------------*/
/* write restart */
void PARTICLE::TimInt::OutputRestart
(
  bool& datawritten
)
{
  // Yes, we are going to write...
  datawritten = true;

  // mesh is written to disc
  output_->ParticleOutput(step_, (*time_)[0], true);
  output_->NewStep(step_, (*time_)[0]);

  OutputDisplacement();
  WriteVector("velocity", vel_);
  WriteVector("acceleration", acc_);

  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"TRANSPORT_VELOCITY")==true)
  {
    WriteVector("modified_velocity", velmod_);
    WriteVector("modified_acceleration", accmod_);
  }

  WriteVector("radius", radius_, false);
  WriteVector("mass", mass_, false);

  WriteVector("densityDot", densityDot_, false);
  WriteVector("pressure", pressure_, false);

  WriteVector("density", density_, false);

  if(colorField_!=Teuchos::null)
    WriteVector("colorField", colorField_, false);

  if(colorFieldGrad_!=Teuchos::null)
    WriteVector("colorFieldGrad", colorFieldGrad_);

  if(smoothedColorFieldGrad_!=Teuchos::null)
    WriteVector("smoothedColorFieldGrad", smoothedColorFieldGrad_);

  if(accSF_!=Teuchos::null)
    WriteVector("accSF", accSF_);

  if(fspType_!=Teuchos::null)
    WriteVector("fspType", fspType_, false);

  if(curvature_!=Teuchos::null)
    WriteVector("curvature", curvature_, false);

  if(phaseColor_!=Teuchos::null)
    WriteVector("phaseColor", phaseColor_, false);

  if(variableradius_)
  {
    WriteVector("radius0", radius0_, false);
    WriteVector("radiusDot", radiusDot_, false);
  }

  if(collhandler_ != Teuchos::null)
  {
    if(angVeln_ != Teuchos::null)
    {
      WriteVector("ang_velocity", (*angVel_)(0));
      WriteVector("ang_acceleration", (*angAcc_)(0));
    }

    strategy_->OutputOrientation();
  }

  // maps are rebuild in every step so that reuse is not possible
  // keeps memory usage bounded
  output_->ClearMapCache();

  // info dedicated to user's eyes staring at standard out
  if ( (myrank_ == 0) and printscreen_ and ((step_-restart_)%printscreen_==0))
  {
    printf("====== Restart for field 'Particle' written in step %d\n", step_);
    fflush(stdout);
  }

  // info dedicated to processor error file
  if (printerrfile_)
  {
    fprintf(errfile_, "====== Restart for field 'Particle' written in step %d\n", step_);
    fflush(errfile_);
  }

  return;
}

/*----------------------------------------------------------------------*/
/* output states */
void PARTICLE::TimInt::OutputState
(
  bool& datawritten
)
{
  // Yes, we are going to write...
  datawritten = true;

  // mesh is not written to disc, only maximum node id is important for output
  output_->ParticleOutput(step_, (*time_)[0], false);
  output_->NewStep(step_, (*time_)[0]);

  OutputDisplacement();
  WriteVector("velocity", vel_);
  WriteVector("radius", radius_, false);

  if(writevelacc_)
  {
    WriteVector("acceleration", acc_);
    WriteVector("densityDot", densityDot_, false);
  }
  WriteVector("pressure", pressure_, false);
  //WriteVector("densityapprox", densityapproxn_, false);
  WriteVector("density", density_, false);

  if(colorField_!=Teuchos::null)
    WriteVector("colorField", colorField_, false);

  if(colorFieldGrad_!=Teuchos::null)
    WriteVector("colorFieldGrad", colorFieldGrad_);

  if(smoothedColorFieldGrad_!=Teuchos::null)
    WriteVector("smoothedColorFieldGrad", smoothedColorFieldGrad_);

  if(accSF_!=Teuchos::null)
    WriteVector("accSF", accSF_);

  if(fspType_!=Teuchos::null)
    WriteVector("fspType", fspType_, false);

  if(curvature_!=Teuchos::null)
    WriteVector("curvature", curvature_, false);

  if(phaseColor_!=Teuchos::null)
    WriteVector("phaseColor", phaseColor_, false);

  WriteVector("modified_velocity", velmod_);

  if(collhandler_ != Teuchos::null)
    strategy_->OutputOrientation();

  // maps are rebuild in every step so that reuse is not possible
  // keeps memory usage bounded
  output_->ClearMapCache();

  return;
}

/*----------------------------------------------------------------------*/
/* Calculation of internal, external and kinetic energy */
void PARTICLE::TimInt::DetermineEnergy()
{

  if ( writeenergyevery_ and (stepn_%writeenergyevery_ == 0))
  {
    LINALG::Matrix<3,1> gravity_acc = particle_algorithm_->GetGravityAcc(timen_);

    int numrownodes = discret_->NodeRowMap()->NumMyElements();

    //energy for all cases besides SPH
    if (collhandler_ != Teuchos::null and particle_algorithm_->ParticleInteractionType()!=INPAR::PARTICLE::SPH )
    {
      // total kinetic energy
      kinergy_ = strategy_->ComputeKineticEnergy();

      for(int i=0; i<numrownodes; ++i)
      {
        double specific_energy = 0.0;

        for(int dim=0; dim<3; ++dim)
          specific_energy -=  gravity_acc(dim) * (*disn_)[i*3+dim];

        intergy_ += (*mass_)[i] * specific_energy;
      }
      // total external energy not available
      extergy_ = 0.0;
    }//energy for SPH case
    else if (collhandler_ == Teuchos::null and particle_algorithm_->ParticleInteractionType()==INPAR::PARTICLE::SPH )
    {
      // set energies to zero (intergy_ has already been calculated earlier)
      double specific_energy = 0.0;
      kinergy_ = 0.0;
      extergy_ = 0.0;
      linmomentum_.Clear();

      for(int i=0; i<numrownodes; ++i)
      {
        double kinetic_energy = 0.0;

        for(int dim=0; dim<3; ++dim)
        {
          // gravitation
          specific_energy +=  gravity_acc(dim) * (*disn_)[i*3+dim];

          // kinetic energy
          kinetic_energy += pow((*veln_)[i*3+dim], 2.0 );
          linmomentum_(dim) +=(*mass_)[i]*(*veln_)[i*3+dim];
        }

        extergy_ += (*mass_)[i]*specific_energy;
        kinergy_ += 0.5 * (*mass_)[i] * kinetic_energy;
      }
    }
    else
      dserror("Energy output only possible for DEM (with collhandler_ == true) or for SPH interaction (with collhandler_ == false)");

    double global_energy[3] = {0.0, 0.0, 0.0};
    double energies[3] = {intergy_, kinergy_, extergy_};
    discret_->Comm().SumAll(&energies[0], &global_energy[0], 3);

    intergy_ = global_energy[0];
    kinergy_ = global_energy[1];
    extergy_ = global_energy[2];

    double global_linmomentum[3] = {0.0, 0.0, 0.0};
    double linmomentum[3] = {linmomentum_(0), linmomentum_(1), linmomentum_(2)};
    discret_->Comm().SumAll(&linmomentum[0], &global_linmomentum[0], 3);

    linmomentum_(0) = global_linmomentum[0];
    linmomentum_(1) = global_linmomentum[1];
    linmomentum_(2) = global_linmomentum[2];

  }

  return;
}

/*----------------------------------------------------------------------*/
/* output system energies */
void PARTICLE::TimInt::OutputEnergy()
{
  // total energy
  double totergy = kinergy_ + intergy_ - extergy_;

  // the output
  if (myrank_ == 0)
  {
    //energy for all cases besides SPH
    if (collhandler_ != Teuchos::null and particle_algorithm_->ParticleInteractionType()!=INPAR::PARTICLE::SPH )
    {
      *energyfile_ << " " << std::setw(9) << step_
                   << std::scientific  << std::setprecision(16)
                   << " " << (*time_)[0]
                   << " " << totergy
                   << " " << kinergy_
                   << " " << intergy_
                   << " " << extergy_
                   << " " << collhandler_->GetMaxPenetration()
                   << std::endl;
    }//energy for SPH case
    else if (collhandler_ == Teuchos::null and particle_algorithm_->ParticleInteractionType()==INPAR::PARTICLE::SPH )
    {
      *energyfile_ << step_
                   << std::scientific  << std::setprecision(16)
                   << " " << (*time_)[0]
                   << " " << totergy
                   << " " << kinergy_
                   << " " << intergy_
                   << " " << extergy_
                   << " " << linmomentum_(0)
                   << " " << linmomentum_(1)
                   << " " << linmomentum_(2)
                   << std::endl;
    }
    else
      dserror("Energy output only possible for DEM (with collhandler_ == true) or for SPH interaction (with collhandler_ == false)");
  }
  return;
}

/*----------------------------------------------------------------------*
 | SPH rendering output                                    sfuchs 06/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimInt::PerformSPHRendering(bool clearstate, bool writeoutput)
{
  if(interHandler_ != Teuchos::null  and particle_algorithm_->ParticleInteractionType()==INPAR::PARTICLE::SPH)
  {
    Teuchos::RCP<Rendering> rendering = particle_algorithm_->GetRendering();
    if (rendering != Teuchos::null)
    {
      // clear the (averaged) rendering vectors
      if (clearstate)
      {
        std::cout << "Clear rendering state vectors before first averaging step!" << std::endl;
        rendering->ClearState();
      }

      // determine the (averaged) rendering vectors
      std::cout << "Update (averaged) rendering vectors!" << std::endl;
      if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"TRANSPORT_VELOCITY")==true)
        rendering->UpdateRenderingVectors(discret_, (*dis_)(0), (*vel_)(0), (*acc_)(0), (*velmod_)(0), (*density_)(0), (*radius_)(0), pressure_, mass_);
      else
        rendering->UpdateRenderingVectors(discret_, (*dis_)(0), (*vel_)(0), (*acc_)(0), Teuchos::null, (*density_)(0), (*radius_)(0), pressure_, mass_);

      // write the (averaged) rendering vectors
      if (writeoutput)
      {
        std::cout << "Write rendering output!" << std::endl;
        rendering->OutputState();
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | output particle statistics                                fang 04/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimInt::OutputParticleStatistics()
{
  // compute total particle volume
  double myvolume(0.), globalvolume(0.);
  for(int i=0; i<(*radius_)(0)->MyLength(); ++i)
    myvolume += 4./3.*M_PI*pow((*(*radius_)(0))[i],3);
  discret_->Comm().SumAll(&myvolume,&globalvolume,1);

  // determine minimum and maximum particle radius
  double r_min(0.), r_max(0.);
  (*radius_)(0)->MinValue(&r_min);
  (*radius_)(0)->MaxValue(&r_max);

  // output only performed by first processor
  if(myrank_ == 0)
    *particlestatsfile_ << std::setw(10) << step_
                        << std::scientific  << std::setprecision(10) << " " << (*time_)[0]
                        << " " << discret_->NodeRowMap()->NumGlobalElements()
                        << " " << r_min
                        << " " << r_max
                        << " " << globalvolume
                        << std::endl;

  return;
}


/*----------------------------------------------------------------------*/
/* Attach file handle for energy file #energyfile_                      */
void PARTICLE::TimInt::AttachEnergyFile()
{
  if (energyfile_.is_null())
  {
    // energy for all cases besides SPH
    if(particle_algorithm_->ParticleInteractionType() != INPAR::PARTICLE::SPH)
    {
      std::string energyname
        = DRT::Problem::Instance()->OutputControlFile()->FileName()
        + "_particle.energy";
      energyfile_ = Teuchos::rcp(new std::ofstream(energyname.c_str()));
      (*energyfile_) << "# timestep time total_energy"
                     << " kinetic_energy internal_energy external_energy max_particle_penetration"
                     << std::endl;
    }

    // energy for SPH case
    else
    {
      std::string energyname
        = DRT::Problem::Instance()->OutputControlFile()->FileName()
        + "_particle.energy";
      energyfile_ = Teuchos::rcp(new std::ofstream(energyname.c_str()));
      (*energyfile_) << "timestep time total_energy"
                     << " kinetic_energy internal_energy external_energy x_lin_momentum y_lin_momentum z_lin_momentum"
                     << std::endl;
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | attach file handle for particle statistics                fang 04/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimInt::AttachParticleStatisticsFile()
{
  // create file if not yet existent
  if(particlestatsfile_ == Teuchos::null and myrank_ == 0)
  {
    // set file name
    const std::string filename(DRT::Problem::Instance()->OutputControlFile()->FileName()+"_particle.statistics.csv");

    // initialize file
    particlestatsfile_ = Teuchos::rcp(new std::ofstream(filename.c_str()));

    // write header
    *particlestatsfile_ << "# timestep time number r_min r_max volume" << std::endl;
  }

  return;
}


/*----------------------------------------------------------------------*/
/* Creates the field test                                               */
Teuchos::RCP<DRT::ResultTest> PARTICLE::TimInt::CreateFieldTest()
{
  return Teuchos::rcp(new PartResultTest(*this));
}

/*----------------------------------------------------------------------*/
/* dof map of vector of unknowns                                        */
Teuchos::RCP<const Epetra_Map> PARTICLE::TimInt::DofRowMap()
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  return Teuchos::rcp(new Epetra_Map(*dofrowmap));
}

/*----------------------------------------------------------------------*/
/* view of dof map of vector of unknowns                                */
const Epetra_Map* PARTICLE::TimInt::DofRowMapView()
{
  return discret_->DofRowMap();
}

/*----------------------------------------------------------------------*/
/* node map of particles                                                */
Teuchos::RCP<const Epetra_Map> PARTICLE::TimInt::NodeRowMap()
{
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  return Teuchos::rcp(new Epetra_Map(*noderowmap));
}

/*----------------------------------------------------------------------*/
/* view of node map of particles                                        */
const Epetra_Map* PARTICLE::TimInt::NodeRowMapView()
{
  return discret_->NodeRowMap();
}


/*-----------------------------------------------------------------------------*/
/* Update exporter objects if layout has changed                               */
void PARTICLE::TimInt::UpdateExportersIfNecessary(const Epetra_BlockMap& oldnodemap, const Epetra_BlockMap& olddofmap)
{
  if(nodemapexporter_ == Teuchos::null
        or not nodemapexporter_->TargetMap().SameAs(*NodeRowMapView())
        or not nodemapexporter_->SourceMap().SameAs(oldnodemap)
     or dofmapexporter_ == Teuchos::null
        or not dofmapexporter_->TargetMap().SameAs(*DofRowMapView())
        or not dofmapexporter_->SourceMap().SameAs(olddofmap)
     )
  {
    nodemapexporter_ = Teuchos::rcp(new Epetra_Export(oldnodemap, *NodeRowMapView()));
    dofmapexporter_ = Teuchos::rcp(new Epetra_Export(olddofmap, *DofRowMapView()));
  }

  return;
}

/*-----------------------------------------------------------------------------*/
/* Update TimIntMStep state vector with the new (appropriate) map from discret_*/
void PARTICLE::TimInt::UpdateStateVectorMap(Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > &stateVector, bool trg_nodeVectorType)
{
  if (stateVector != Teuchos::null and (*stateVector)(0) != Teuchos::null)
  {
    const Teuchos::RCP<Epetra_Vector> old = Teuchos::rcp(new Epetra_Vector(*(*stateVector)(0)));

    if (trg_nodeVectorType)
    {
      // create new vector
      stateVector->ReplaceMaps(NodeRowMapView());

      // transfer data
      int err = (*stateVector)(0)->Export(*old, *nodemapexporter_, Insert);
      if (err)
        dserror("Export using exporter returned err=%d", err);
    }
    else
    {
      // create new vector
      stateVector->ReplaceMaps(DofRowMapView());

      // transfer data
      int err = (*stateVector)(0)->Export(*old, *dofmapexporter_, Insert);
      if (err)
        dserror("Export using exporter returned err=%d", err);
    }
  }
}

/*-----------------------------------------------------------------------------*/
/* Update state vector with the new (appropriate) map from discret_*/
void PARTICLE::TimInt::UpdateStateVectorMap(Teuchos::RCP<Epetra_Vector> &stateVector, bool trg_nodeVectorType)
{
  if (stateVector != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> old = stateVector;

    if (trg_nodeVectorType)
    {
      // create new vector
      stateVector = LINALG::CreateVector(*discret_->NodeRowMap(),true);

      // transfer data
      int err = stateVector->Export(*old, *nodemapexporter_, Insert);
      if (err)
        dserror("Export using exporter returned err=%d", err);
    }
    else
    {
      // create new vector
      stateVector = LINALG::CreateVector(*discret_->DofRowMap(),true);

      // transfer data
      int err = stateVector->Export(*old, *dofmapexporter_, Insert);
      if (err)
        dserror("Export using exporter returned err=%d", err);
    }
  }
}

/*-----------------------------------------------------------------------------*/
/* Update TimIntMStep state vector with the new dof map from discret_*/
void PARTICLE::TimInt::UpdateStateVectorMap(Teuchos::RCP<Epetra_MultiVector> &stateVector, bool trg_nodeVectorType)
{
  if (stateVector != Teuchos::null)
  {
    Teuchos::RCP<Epetra_MultiVector> old = stateVector;
    stateVector = Teuchos::rcp(new Epetra_MultiVector(*DofRowMapView(), old->NumVectors(), true));
    int err = stateVector->Export(*old, *dofmapexporter_, Insert);
    if (err)
      dserror("Export using exporter returned err=%d", err);

  }
}

/*----------------------------------------------------------------------*/
//! Check exixstence of the state vectors
void PARTICLE::TimInt::CheckStateVector(std::string vecName, const Teuchos::RCP<const Epetra_Vector> vec, bool trg_showVec)
{
  if (vec == Teuchos::null)
  {
    std::cout << "The pointer to " << vecName << " is null\n";
    std::cin.get();
    return;
  }

  std::cout << "Processor " << myrank_ << " owns " << vec->MyLength() << "/" << vec->GlobalLength() << " elements of " << vecName << std::endl;
  if (trg_showVec)
    std::cout << *vec << std::endl;
  std::cin.get();
}


/*----------------------------------------------------------------------------*
 | return maximum particle-particle or particle-wall penetration   fang 02/17 |
 *----------------------------------------------------------------------------*/
double PARTICLE::TimInt::MaximumPenetration() const
{
  return collhandler_ == Teuchos::null ? 0. : collhandler_->GetMaxPenetration();
}


/*----------------------------------------------------------------------*/
/* update step */
void PARTICLE::TimInt::UpdateStepState()
{
  // new state vectors at t_{n+1} -> t_n
  //    SV_{n} := SV_{n+1}, SV_{n-1} := SV_{n}
  UpdateStateVector(dis_, disn_);
  UpdateStateVector(vel_, veln_);
  UpdateStateVector(acc_, accn_);

  if(interHandler_ != Teuchos::null)
  {
    UpdateStateVector(radius_, radiusn_);
    UpdateStateVector(density_, densityn_);
    UpdateStateVector(densityDot_, densityDotn_);

    interHandler_->Density2Pressure(densityn_, pressure_);

    if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"TRANSPORT_VELOCITY")==true)
    {
      UpdateStateVector(velmod_, velmodn_);
      UpdateStateVector(accmod_, accmodn_);
    }
  }

  if(collhandler_ != Teuchos::null)
  {
    // new angular-velocities at t_{n+1} -> t_n
    //    ang_V_{n} := ang_V_{n+1}, ang_V_{n-1} := ang_V_{n}
    UpdateStateVector(angVel_, angVeln_);
    // new angular-accelerations at t_{n+1} -> t_n
    //    ang_A_{n} := ang_A_{n+1}, ang_A_{n-1} := ang_A_{n}
    UpdateStateVector(angAcc_, angAccn_);
  }

  return;
}


/*----------------------------------------------------------------------*/
// wrapper. On top of the output_->WriteVector() it checks that the pointer is not null. In case, it does not write
void PARTICLE::TimInt::WriteVector(const std::string name,
                                           Teuchos::RCP<Epetra_Vector> vec,
                                           const bool isdof) const
{
  if (vec != Teuchos::null)
  {
    if (isdof)
    {
      output_->WriteVector(name, vec, IO::dofvector);
    }
    else
    {
      output_->WriteVector(name, vec, IO::nodevector);
    }
  }
}

/*----------------------------------------------------------------------*/
// wrapper. On top of the output_->WriteVector() it checks that the pointer is not null. In case, it does not write
void PARTICLE::TimInt::WriteVector(const std::string name,
                                           Teuchos::RCP<Epetra_MultiVector> vec,
                                           const bool isdof) const
{
  if (vec != Teuchos::null)
  {
    if (isdof)
    {
      output_->WriteVector(name, vec, IO::dofvector);
    }
    else
    {
      output_->WriteVector(name, vec, IO::nodevector);
    }
  }
}

/*----------------------------------------------------------------------*/
// wrapper. On top of the output_->WriteVector() it checks that the pointer is not null. In case, it does not write
void PARTICLE::TimInt::WriteVector(const std::string name,
                                   Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > vec,
                                   const bool isdof) const
{
  if (vec != Teuchos::null && (*vec)(0) != Teuchos::null)
  {
    if (isdof)
    {
      output_->WriteVector(name, (*vec)(0), IO::dofvector);
    }
    else
    {
      output_->WriteVector(name, (*vec)(0), IO::nodevector);
    }
  }
}


/*----------------------------------------------------------------------*
 | calculate forces on particle and apply it               ghamm 02/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::TimInt::UpdateExtActions(bool init)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::UpdateExtActions");

  LINALG::Matrix<3,1> gravity_acc = particle_algorithm_->GetGravityAcc(-1.0);

  // forces/accelerations
  if (particle_algorithm_->ParticleInteractionType() == INPAR::PARTICLE::SPH)
  {
    //In SPH simulations, gravity and all other interaction forces are considered at once in the methdo IntegrateStep()
  }
  else if (fifc_ != Teuchos::null)
  {
    fifc_->PutScalar(0.0);
    GravityForces(fifc_);
  }

  // densityDot, in case it exists
  if (densityDotn_ != Teuchos::null)
  {
    densityDotn_->PutScalar(0);
  }

  particle_algorithm_->CalculateAndApplyForcesToParticles(init);

  return;
}


/*----------------------------------------------------------------------*
 | calculate and ADD gravity forces (no reset)             katta 01/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::TimInt::GravityForces(Teuchos::RCP<Epetra_Vector> force)
{
  if (force != Teuchos::null)
  {
    LINALG::Matrix<3,1> gravity_acc = particle_algorithm_->GetGravityAcc(-1.0);
    for (int i=0; i<discret_->NodeRowMap()->NumMyElements(); ++i)
    {
      /*------------------------------------------------------------------*/
      //// gravity acc = mass_p * g
      for(int dim=0; dim<3; ++dim)
      {
        (*force)[i*3+dim] = (*mass_)[i] * gravity_acc(dim);
      }
      /*------------------------------------------------------------------*/
    }
  }
  else
  {
    dserror("You are trying to apply gravity forces to a null pointer");
  }
}

/*----------------------------------------------------------------------*
 | calculate and ADD gravity forces (no reset)             katta 01/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::TimInt::GravityAcc(Teuchos::RCP<Epetra_Vector> acc, const double time)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::TimInt::GravityAcc");

  if (acc != Teuchos::null)
  {
    LINALG::Matrix<3,1> gravity_acc = particle_algorithm_->GetGravityAcc(time);

    for (int i=0; i<discret_->NodeRowMap()->NumMyElements(); ++i)
    {
      /*------------------------------------------------------------------*/
      //// gravity acc = g
      for(int dim=0; dim<3; ++dim)
      {
        (*acc)[i*3+dim] = gravity_acc(dim);
      }
      /*------------------------------------------------------------------*/
    }
  }
  else
  {
    dserror("You are trying to apply gravity accelerations to a null pointer");
  }
}

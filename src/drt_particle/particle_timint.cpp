/*----------------------------------------------------------------------*/
/*!
\file particle_timint.cpp

\brief Time integration for particle dynamics

\level 1

\maintainer Georg Hammerl
*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint.H"
#include "particle_algorithm.H"
#include "particle_contact.H"
#include "particle_resulttest.H"
#include "particleMeshFree_interaction.H"
#include "particleMeshFree_weightFunction.H"
#include "particleMeshFree_rendering.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/extparticle_mat.H"

#include <Teuchos_TimeMonitor.hpp>

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
  bpDoFs_(Teuchos::null),
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
  energyfile_(Teuchos::null),
  writeorientation_(false),
  kinergy_(0),
  intergy_(0),
  bpintergy_(0),
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
  angVel_(Teuchos::null),
  angAcc_(Teuchos::null),
  radius_(Teuchos::null),
  density_(Teuchos::null),
  densityDot_(Teuchos::null),
  specEnthalpy_(Teuchos::null),
  specEnthalpyDot_(Teuchos::null),

  disn_(Teuchos::null),
  veln_(Teuchos::null),
  accn_(Teuchos::null),
  angVeln_(Teuchos::null),
  angAccn_(Teuchos::null),
  radiusn_(Teuchos::null),
  densityn_(Teuchos::null),
  densityDotn_(Teuchos::null),
  specEnthalpyn_(Teuchos::null),
  specEnthalpyDotn_(Teuchos::null),

  fifc_(Teuchos::null),
  orient_(Teuchos::null),

  radius0_(Teuchos::null),
  radiusDot_(Teuchos::null),
  mass_(Teuchos::null),
  inertia_(Teuchos::null),
  temperature_(Teuchos::null),
  pressure_(Teuchos::null),
  f_structure_(Teuchos::null),

  initDensity_(-1.0),
  restDensity_(-1.0),
  refdensfac_(-1.0),

  dofmapexporter_(Teuchos::null),
  nodemapexporter_(Teuchos::null),

  radiusdistribution_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::RadiusDistribution>(particledynparams,"RADIUS_DISTRIBUTION")),
  variableradius_((bool)DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->CavitationParams(),"COMPUTE_RADIUS_RP_BASED")),
  radiuschangecurve_(particledynparams.get<int>("RADIUS_CHANGE_CURVE")),
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
  // initialize the vectors
  dis_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
  vel_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
  acc_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
  radius_  = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, NodeRowMapView(), true));

  fifc_ = LINALG::CreateVector(*DofRowMapView(), true);
  mass_ = LINALG::CreateVector(*NodeRowMapView(), true);

  switch (particle_algorithm_->ParticleInteractionType())
  {
  case INPAR::PARTICLE::MeshFree :
  {
    pressure_ = LINALG::CreateVector(*NodeRowMapView(), true);
    densityDot_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, NodeRowMapView(), true));
  }// no break
  case INPAR::PARTICLE::Normal_DEM_thermo :
  {
    density_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, NodeRowMapView(), true));
    specEnthalpy_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, NodeRowMapView(), true));
    specEnthalpyDot_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, NodeRowMapView(), true));
    temperature_ = LINALG::CreateVector(*NodeRowMapView(), true);
    break;
  }
  default : //do nothing
    break;
  }

  if(variableradius_ or radiuschangecurve_ > 0)
  {
    // initial radius of each particle for time dependent radius
    radius0_  = LINALG::CreateVector(*discret_->NodeRowMap(), true);
    // time derivative of radius of each particle for time dependent radius
    radiusDot_  = LINALG::CreateVector(*NodeRowMapView(), true);
  }

  //Vector to identify possible boundary particle DoFs
  const INPAR::PARTICLE::WallInteractionType wallInteractionType=DRT::INPUT::IntegralValue<INPAR::PARTICLE::WallInteractionType>(DRT::Problem::Instance()->ParticleParams(),"WALL_INTERACTION_TYPE");
  if(wallInteractionType == INPAR::PARTICLE::BoundarParticle_FreeSlip or wallInteractionType == INPAR::PARTICLE::BoundarParticle_NoSlip)
    bpDoFs_ = LINALG::CreateVector(*discret_->DofRowMap(), true);

  // Apply Dirichlet BC and create dbc map extractor
  {
    dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
    Teuchos::ParameterList p;
    p.set("total time", (*time_)[0]);
    discret_->EvaluateDirichlet(p, (*dis_)(0), (*vel_)(0), (*acc_)(0), bpDoFs_, dbcmaps_);
  }

  // set initial fields
  SetInitialFields();

  // copy everything into the n+1 state vectors
  disn_ = Teuchos::rcp(new Epetra_Vector(*(*dis_)(0)));
  veln_ = Teuchos::rcp(new Epetra_Vector(*(*vel_)(0)));
  accn_ = Teuchos::rcp(new Epetra_Vector(*(*acc_)(0)));

  switch (particle_algorithm_->ParticleInteractionType())
  {
  case INPAR::PARTICLE::MeshFree :
  {
    densityDotn_ = Teuchos::rcp(new Epetra_Vector(*(*densityDot_)(0)));
  }// no break
  case INPAR::PARTICLE::Normal_DEM_thermo :
  {
    radiusn_  = Teuchos::rcp(new Epetra_Vector(*(*radius_)(0)));
    densityn_ = Teuchos::rcp(new Epetra_Vector(*(*density_)(0)));
    specEnthalpyn_ = Teuchos::rcp(new Epetra_Vector(*(*specEnthalpy_)(0)));
    specEnthalpyDotn_ = Teuchos::rcp(new Epetra_Vector(*(*specEnthalpyDot_)(0)));
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

    if(writeorientation_)
    {
      // initialize orientation-vector for visualization
      orient_ = LINALG::CreateVector(*DofRowMapView(),true);
      InitializeOrientVector();
    }

    // create and fill inertia
    PARTICLE::Utils::ComputeInertia((*radius_)(0), mass_, inertia_, true);
  }

  // output file for energy
  if(writeenergyevery_ != 0 and myrank_ == 0)
    AttachEnergyFile();

  return;
}

/*----------------------------------------------------------------------*/
/* Set initial fields in structure (e.g. initial velocities) */
void PARTICLE::TimInt::SetInitialFields()
{
  // -----------------------------------------//
  // set material properties
  // -----------------------------------------//

  // all particles have identical density and radius (for now)
  const double initRadius = particle_algorithm_->ParticleMat()->initRadius_;
  const double initDensity = particle_algorithm_->ParticleMat()->initDensity_;

  double amplitude = DRT::Problem::Instance()->ParticleParams().get<double>("RANDOM_AMPLITUDE");

  (*radius_)(0)->PutScalar(initRadius);

  double consistent_problem_volume=DRT::Problem::Instance()->ParticleParams().get<double>("CONSISTENT_PROBLEM_VOLUME");
  if(consistent_problem_volume<0.0) // particle mass via m = rho * 4/3 * PI *r^3
    mass_->PutScalar(initDensity * PARTICLE::Utils::Radius2Volume(initRadius));
  else // particle mass via problem volume and density
  {
    //boundary particles are excluded for determination of particle mass since they are also not part of the consistent_problem_volume
    int num_boundaryparticles=0;

    const INPAR::PARTICLE::WallInteractionType wallInteractionType=DRT::INPUT::IntegralValue<INPAR::PARTICLE::WallInteractionType>(DRT::Problem::Instance()->ParticleParams(),"WALL_INTERACTION_TYPE");
    if(wallInteractionType == INPAR::PARTICLE::BoundarParticle_FreeSlip or wallInteractionType == INPAR::PARTICLE::BoundarParticle_NoSlip)
    {
      double num_boundaryDoFs = 0.0;
      bpDoFs_->Norm1(&num_boundaryDoFs);
      num_boundaryparticles=(int)(num_boundaryDoFs/3);
      int tot_num_boundaryparticles=0;
      if(wallInteractionType == INPAR::PARTICLE::BoundarParticle_FreeSlip or wallInteractionType == INPAR::PARTICLE::BoundarParticle_NoSlip)
      discret_->Comm().SumAll(&num_boundaryparticles,&tot_num_boundaryparticles,1);
    }

    const int num_particles=discret_->NumGlobalNodes()-num_boundaryparticles;
    mass_->PutScalar(initDensity * consistent_problem_volume / (double)(num_particles));
  }

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
        double function_value =  DRT::Problem::Instance()->Funct(funct_num-1).Evaluate(0, currparticle->X(),0.0,discret_.get());
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

  // -----------------------------------------//
  // initialize displacement field
  // -----------------------------------------//

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
  // set initial velocity field if existing
  // -----------------------------------------//

  const std::string field = "Velocity";
  std::vector<int> localdofs;
  localdofs.push_back(0);
  localdofs.push_back(1);
  localdofs.push_back(2);
  discret_->EvaluateInitialField(field,(*vel_)(0),localdofs);

  // -----------------------------------------//
  // set the other parameters. In case of meshfree set also pressure and density
  // -----------------------------------------//

  const MAT::PAR::ExtParticleMat* extParticleMat = particle_algorithm_->ExtParticleMat();
  switch (particle_algorithm_->ParticleInteractionType())
  {
  case INPAR::PARTICLE::MeshFree :
  {
    // set the rest density used for pressure-related dynamics
    restDensity_ = extParticleMat->restDensity_;
    refdensfac_ = extParticleMat->refdensfac_;
    initDensity_ = initDensity;

    // set density in the density vector (useful only for thermodynamics)
    (*density_)(0)->PutScalar(initDensity);

    if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"SOLVE_THERMAL_PROBLEM"))
    {
      // initialize temperature of particles
      const double initTemperature = extParticleMat->initTemperature_;
      const double transitionTemperature = extParticleMat->transitionTemperature_;
      const double tempDiff = initTemperature - transitionTemperature;

      if (tempDiff > 0)
        (*specEnthalpy_)(0)->PutScalar(extParticleMat->SpecEnthalpyTL() + tempDiff * extParticleMat->CPL_);
      else if (initTemperature < transitionTemperature)
        (*specEnthalpy_)(0)->PutScalar(initTemperature * extParticleMat->CPS_);
      else
        dserror("TODO: start in the transition point - solid <-> liquid - still not implemented");

      UpdateTemperaturen();
    }

    break;
  }
  case INPAR::PARTICLE::Normal_DEM_thermo :
  {
    // set density in the density vector (useful only for thermodynamics)
    (*density_)(0)->PutScalar(initDensity);

    // initialize temperature of particles
    const double initTemperature = extParticleMat->initTemperature_;
    const double transitionTemperature = extParticleMat->transitionTemperature_;
    const double tempDiff = initTemperature - transitionTemperature;

    if (tempDiff > 0)
      (*specEnthalpy_)(0)->PutScalar(extParticleMat->SpecEnthalpyTL() + tempDiff * extParticleMat->CPL_);
    else if (initTemperature < transitionTemperature)
      (*specEnthalpy_)(0)->PutScalar(initTemperature * extParticleMat->CPS_);
    else
      dserror("TODO: start in the transition point - solid <-> liquid - still not implemented");

    UpdateTemperaturen();
    break;
  }
  default : //do nothing
    break;
  }

  // It is here because it is a slave of the density
  if (particle_algorithm_->ParticleInteractionType() == INPAR::PARTICLE::MeshFree)
  {
    Teuchos::RCP<Epetra_Vector> deltaDensity = Teuchos::rcp(new Epetra_Vector(*(discret_->NodeRowMap()), true));
    deltaDensity->PutScalar(refdensfac_*restDensity_);
    deltaDensity->Update(1.0, *(*density_)(0), -1.0);
    bool solve_thermal_problem=DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"SOLVE_THERMAL_PROBLEM");
    PARTICLE::Utils::Density2Pressure(deltaDensity, (*specEnthalpy_)(0), pressure_, particle_algorithm_->ExtParticleMat(), true, solve_thermal_problem);
  }

  // set vector of initial particle radii if necessary
  if(radiuschangecurve_ > 0)
    radius0_->Update(1.,*(*radius_)(0),0.);

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

    // do particle business
    particle_algorithm_->TransferParticles(true);
  }

  // update particle radii if necessary
  if(radiuschangecurve_ > 0)
    (*radius_)(0)->Update(DRT::Problem::Instance()->Curve(radiuschangecurve_-1).f(timen_),*radius0_,0.);

  return;
}

/*----------------------------------------------------------------------*/
/* equilibrate system at initial state and identify consistent accelerations */
void PARTICLE::TimInt::DetermineMassDampConsistAccel()
{
  ComputeAcc(Teuchos::null, Teuchos::null, (*acc_)(0), Teuchos::null);

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
    for(int i=0; i<numrownodes; ++i)
    {
      const double invinertia = 1.0/(*inertia_)[i];
      for(int dim=0; dim<3; ++dim)
        (*global_angAcc)[i*3+dim] = invinertia * (*m_contact)[i*3+dim];
    }
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
  UpdateExportersIfNecessary((*radius_)(0)->Map(),(*dis_)(0)->Map());

  UpdateStateVectorMap(dis_);
  UpdateStateVectorMap(vel_);
  UpdateStateVectorMap(acc_);
  UpdateStateVectorMap(angVel_);
  UpdateStateVectorMap(angAcc_);
  UpdateStateVectorMap(radius_,true);
  UpdateStateVectorMap(density_,true);
  UpdateStateVectorMap(densityDot_,true);
  UpdateStateVectorMap(specEnthalpy_,true);
  UpdateStateVectorMap(specEnthalpyDot_,true);
  UpdateStateVectorMap(disn_);
  UpdateStateVectorMap(veln_);
  UpdateStateVectorMap(accn_);
  UpdateStateVectorMap(angVeln_);
  UpdateStateVectorMap(angAccn_);
  UpdateStateVectorMap(radiusn_,true);
  UpdateStateVectorMap(densityn_,true);
  UpdateStateVectorMap(densityDotn_,true);
  UpdateStateVectorMap(specEnthalpyn_,true);
  UpdateStateVectorMap(specEnthalpyDotn_,true);

  UpdateStateVectorMap(fifc_);
  UpdateStateVectorMap(orient_);

  UpdateStateVectorMap(radius0_,true);
  UpdateStateVectorMap(radiusDot_,true);
  UpdateStateVectorMap(mass_,true);
  UpdateStateVectorMap(inertia_,true);
  UpdateStateVectorMap(temperature_,true);
  UpdateStateVectorMap(pressure_,true);
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
  reader.ReadVector(accn_, "acceleration");
  acc_->UpdateSteps(*accn_);

  reader.ReadVector(mass_, "mass");


  if (particle_algorithm_->ParticleInteractionType() != INPAR::PARTICLE::MeshFree ||
      particle_algorithm_->ParticleInteractionType() != INPAR::PARTICLE::Normal_DEM_thermo)
  {
    // create a dummy vector to extract the radius vector (radiusn_ does not exist)
    Teuchos::RCP<Epetra_Vector> radius = LINALG::CreateVector(*discret_->NodeRowMap(), true);
    reader.ReadVector(radius, "radius");
    radius_->UpdateSteps(*radius);
  }

  if (particle_algorithm_->ParticleInteractionType() == INPAR::PARTICLE::MeshFree)
  {
    Teuchos::RCP<Epetra_Vector> deltaDensity = Teuchos::rcp(new Epetra_Vector(*(discret_->NodeRowMap()), true));
    deltaDensity->PutScalar(refdensfac_*restDensity_);
    deltaDensity->Update(1.0,*densityn_,-1.0);
    bool solve_thermal_problem=DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"SOLVE_THERMAL_PROBLEM");
    PARTICLE::Utils::Density2Pressure(deltaDensity,specEnthalpyn_,pressure_,particle_algorithm_->ExtParticleMat(),true,solve_thermal_problem);
  }

  // read in particle collision relevant data
  if(collhandler_ != Teuchos::null)
  {
    // initialize inertia
    PARTICLE::Utils::ComputeInertia((*radius_)(0), mass_, inertia_, true);

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

  WriteVector("displacement", dis_);
  WriteVector("velocity", vel_);
  WriteVector("acceleration", acc_);

  WriteVector("radius", radius_, false);
  WriteVector("mass", mass_, false);

  WriteVector("densityDot", densityDot_, false);
  WriteVector("pressure", pressure_, false);

  WriteVector("density", density_, false);
  WriteVector("specEnthalpy", specEnthalpy_, false);
  WriteVector("specEnthalpyDot", specEnthalpyDot_, false);
  WriteVector("temperature", temperature_, false);

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

    if(writeorientation_)
      WriteVector("orientation", orient_);
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

  Teuchos::RCP<Rendering> rendering = particle_algorithm_->GetRendering();
  if (rendering != Teuchos::null)
  {
    rendering->UpdateStateVectors(discret_, (*dis_)(0), (*vel_)(0),
        (*acc_)(0), (*density_)(0), (*specEnthalpy_)(0), temperature_, (*radius_)(0), pressure_, mass_);
    rendering->OutputState();
  }

  return;
}

/*----------------------------------------------------------------------*/
/* output displacements, velocities, accelerations, temperatures, and pressure */
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

  WriteVector("displacement", dis_);
  WriteVector("velocity", vel_);
  WriteVector("radius", radius_, false);

  if(writevelacc_)
  {
    WriteVector("acceleration", acc_);
    WriteVector("densityDot", densityDot_, false);
    WriteVector("specEnthalpyDot", specEnthalpyDot_, false);
  }
  WriteVector("pressure", pressure_, false);
  //WriteVector("densityapprox", densityapproxn_, false);
  WriteVector("density", density_, false);
  WriteVector("specEnthalpy", specEnthalpy_, false);
  WriteVector("temperature", temperature_, false);

  if(collhandler_ != Teuchos::null and writeorientation_)
    WriteVector("orientation", orient_);
  // maps are rebuild in every step so that reuse is not possible
  // keeps memory usage bounded
  output_->ClearMapCache();

  Teuchos::RCP<Rendering> rendering = particle_algorithm_->GetRendering();
  if (rendering != Teuchos::null)
  {
    rendering->UpdateStateVectors(discret_, (*dis_)(0), (*vel_)(0),
        (*acc_)(0), (*density_)(0), (*specEnthalpy_)(0), temperature_, (*radius_)(0), pressure_, mass_);
    rendering->OutputState();
  }
}

/*----------------------------------------------------------------------*/
/* Calculation of internal, external and kinetic energy */
void PARTICLE::TimInt::DetermineEnergy()
{

  if ( writeenergyevery_ and (stepn_%writeenergyevery_ == 0))
  {
    LINALG::Matrix<3,1> gravity_acc = particle_algorithm_->GetGravityAcc();

    // total kinetic energy
    kinergy_ = 0.0;

    int numrownodes = discret_->NodeRowMap()->NumMyElements();

    //energy for all cases besides SPH
    if (collhandler_ != Teuchos::null and particle_algorithm_->ParticleInteractionType()!=INPAR::PARTICLE::MeshFree )
    {
      for(int i=0; i<numrownodes; ++i)
      {
        double specific_energy = 0.0;
        double kinetic_energy = 0.0;
        double rot_energy = 0.0;

        for(int dim=0; dim<3; ++dim)
        {
          // gravitation
          specific_energy -=  gravity_acc(dim) * (*disn_)[i*3+dim];

          // kinetic energy
          kinetic_energy += pow((*veln_)[i*3+dim], 2.0 );

          // rotation
          rot_energy += pow((*angVeln_)[i*3+dim], 2.0);
        }

        intergy_ += (*mass_)[i] * specific_energy;
        kinergy_ += 0.5 * ((*mass_)[i] * kinetic_energy + (*inertia_)[i] * rot_energy);
      }
      // total external energy not available
      extergy_ = 0.0;
    }//energy for SPH case
    else if (collhandler_ == Teuchos::null and particle_algorithm_->ParticleInteractionType()==INPAR::PARTICLE::MeshFree )
    {
      // set internal and kinetic energy to zero
      double specific_energy = 0.0;
      kinergy_ = 0.0;
      intergy_ = 0.0;
      extergy_ = 0.0;
      linmomentum_.Clear();

      //First add energy contributions of boundary particles (bpintergy_ has already been determined in IntegrateStep())
      //intergy_+=bpintergy_;

      std::cout << "bpintergy_: " << bpintergy_ << std::endl;

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

        // thermodynamic energy E with p=-dE/dV, T=dE/dS (see Espanol2003, Eq.(5))
        // Attention: currently, only the first, pressure-dependent contribution of the thermodynamic energy is implemented!
        // Thus, it is only valid for isentrop problems, i.e.dE/dS=0! Furthermore, it is only considered for the fluid phase so far (since SpeedOfSoundL is used)!
        // In the following, for simplicity, all particles including also boundary particles are considered. The boundary particles are represented by their initial
        // density (usually initDensity_=restDensity) in the vector densityn_ and do not yield energy changes in the following lines. However, the actual energy
        // contributions of the boundary particles is contained in bpintergy_!

        //TODO: So far, energy output is only considered in the context of fluid problems (use of SpeedOfSoundL())!
        if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"SOLVE_THERMAL_PROBLEM"))
          dserror("Energy output is only considered in the context of pure fluid problems so far!");

        double c0 = particle_algorithm_->ExtParticleMat()->SpeedOfSoundL();
        intergy_+=PARTICLE::Utils::Density2Energy(c0, (*densityn_)[i], restDensity_, refdensfac_, (*mass_)[i]);
      }
    }
    else
      dserror("Energy output only possible for DEM (with collhandler_ == true) or for meshfree interaction (with collhandler_ == false)");

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
    if (collhandler_ != Teuchos::null and particle_algorithm_->ParticleInteractionType()!=INPAR::PARTICLE::MeshFree )
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
    else if (collhandler_ == Teuchos::null and particle_algorithm_->ParticleInteractionType()==INPAR::PARTICLE::MeshFree )
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
      dserror("Energy output only possible for DEM (with collhandler_ == true) or for meshfree interaction (with collhandler_ == false)");
  }
  return;
}




/*----------------------------------------------------------------------*/
/* Attach file handle for energy file #energyfile_                      */
void PARTICLE::TimInt::AttachEnergyFile()
{
  if (energyfile_.is_null())
  {
    // energy for all cases besides SPH
    if(particle_algorithm_->ParticleInteractionType() != INPAR::PARTICLE::MeshFree)
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
/* initialization of vector for visualization of the particle orientation */
void PARTICLE::TimInt::InitializeOrientVector()
{
  int numrownodes = discret_->NodeRowMap()->NumMyElements();
  for(int i=0; i<numrownodes; ++i)
  {
    (*orient_)[i*3] = 0.0;
    (*orient_)[i*3+1] = 0.0;
    (*orient_)[i*3+2] = 1.0;
  }
}

/*----------------------------------------------------------------------*/
//! update temperatures \f$T_{n}\f$
void PARTICLE::TimInt::UpdateTemperaturen()
{
  PARTICLE::Utils::SpecEnthalpy2Temperature(temperature_,(*specEnthalpy_)(0), particle_algorithm_->ExtParticleMat());
}

/*----------------------------------------------------------------------*/
//! update temperatures \f$T_{n+1}\f$
void PARTICLE::TimInt::UpdateTemperaturenp()
{
  if (temperature_ != Teuchos::null && specEnthalpyn_ != Teuchos::null)
  {
    PARTICLE::Utils::SpecEnthalpy2Temperature(temperature_,specEnthalpyn_, particle_algorithm_->ExtParticleMat());
  }
}

/*----------------------------------------------------------------------*/
//! update temperatures \f$T_{n+1}\f$
void PARTICLE::TimInt::UpdatePressure()
{
  if (pressure_ != Teuchos::null && specEnthalpyn_ != Teuchos::null && densityn_ != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> deltaDensity = Teuchos::rcp(new Epetra_Vector(*(discret_->NodeRowMap()), true));
    deltaDensity->PutScalar(refdensfac_*restDensity_);
    deltaDensity->Update(1.0,*densityn_,-1.0);
    bool solve_thermal_problem=DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"SOLVE_THERMAL_PROBLEM");
    PARTICLE::Utils::Density2Pressure(deltaDensity, specEnthalpyn_, pressure_, particle_algorithm_->ExtParticleMat(),false,solve_thermal_problem);
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
  UpdateStateVector(radius_, radiusn_);
  UpdateStateVector(density_, densityn_);
  UpdateStateVector(densityDot_, densityDotn_);
  UpdateStateVector(specEnthalpy_, specEnthalpyn_);
  UpdateStateVector(specEnthalpyDot_, specEnthalpyDotn_);

  UpdateTemperaturenp();
  UpdatePressure();

  // legacy... we can probably erase this if (not what is inside tho!) but... whatever
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
/* update of vector for visualization of the particle orientation */
void PARTICLE::TimInt::RotateOrientVector(double dt)
{
  int numrownodes = discret_->NodeRowMap()->NumMyElements();
  for(int i=0; i<numrownodes; ++i)
  {
    double angVel[3];
    double r[3];

    for(int dim=0; dim<3; ++dim)
    {
      angVel[dim] = (*angVeln_)[i*3+dim];
      r[dim] = (*orient_)[i*3+dim];
    }

    // visualization only valid for 2D when rotating around y-axis

    // simplified/linearized orient vector - just for visualization
    // delta_r = \Delta t * (ang_vel x r)
    (*orient_)[i*3]   += dt * (angVel[1] * r[2] - angVel[2] * r[1]);
    (*orient_)[i*3+1] += dt * (angVel[2] * r[0] - angVel[0] * r[2]);
    (*orient_)[i*3+2] += dt * (angVel[0] * r[1] - angVel[1] * r[0]);
    //--------------------------------------------------------------

    //more exactly------------------------------------------------
//    double d[3];
//    double norm_ang_vel=0.0;
//    //norm ang_vel
//    for(int dim=0; dim<3; ++dim)
//    {
//      norm_ang_vel += ang_vel[dim]*ang_vel[dim];
//    }
//    norm_ang_vel = sqrt(norm_ang_vel);
//
//    if(norm_ang_vel > 1E-10)
//    {
//      double scalar = 0.0;
//      double invnorm_ang_vel = 1.0 / norm_ang_vel;
//      for(int dim=0; dim<3; ++dim)
//      {
//        d[dim] = invnorm_ang_vel * ang_vel[dim];
//        scalar += d[dim] * r[dim];
//      }
//
//      for(int dim=0; dim<3; ++dim)
//      {
//        (*orient_)[i*3+dim] += (1-cos(norm_ang_vel*dt)) * scalar * d[dim] - (1-cos(norm_ang_vel*dt)) * r[dim]
//              + sin(norm_ang_vel*dt) * ( d[(dim+1)%3] * r[(dim+2)%3] - d[(dim+2)%3] * r[(dim+1)%3] );
//      }
//    }
    //--------------------------------------------------------------

  }

  return;
}

/*----------------------------------------------------------------------*/
// wrapper. On top of the output_->WriteVector() it checks that the pointer is not null. In case, it does not write
void PARTICLE::TimInt::WriteVector(const std::string name,
                                           Teuchos::RCP<Epetra_Vector> vec,
                                           const bool isdof)
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
                                           const bool isdof)
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
                                   const bool isdof)
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

  LINALG::Matrix<3,1> gravity_acc = particle_algorithm_->GetGravityAcc();

  // forces/accelerations
  if (particle_algorithm_->ParticleInteractionType() == INPAR::PARTICLE::MeshFree)
  {
    accn_->PutScalar(0.0);
    GravityAcc(accn_);
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

  // specEnthalpyDot, in case it exists
  if (specEnthalpyDotn_ != Teuchos::null)
  {
    specEnthalpyDotn_->PutScalar(0);
  }

  particle_algorithm_->CalculateAndApplyForcesToParticles(init);

  return;
}


/*----------------------------------------------------------------------*
 | calculate and ADD gravity forces (no reset)             katta 01/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::TimInt::GravityForces(Teuchos::RCP<Epetra_Vector> force, const double extMulti)
{
  if (force != Teuchos::null)
  {
    LINALG::Matrix<3,1> gravity_acc = particle_algorithm_->GetGravityAcc();
    for (int i=0; i<discret_->NodeRowMap()->NumMyElements(); ++i)
    {
      /*------------------------------------------------------------------*/
      //// gravity acc = mass_p * g
      for(int dim=0; dim<3; ++dim)
      {
        (*force)[i*3+dim] = extMulti * (*mass_)[i] * gravity_acc(dim);
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
void PARTICLE::TimInt::GravityAcc(Teuchos::RCP<Epetra_Vector> acc, const double extMulti)
{
  if (acc != Teuchos::null)
  {
    LINALG::Matrix<3,1> gravity_acc = particle_algorithm_->GetGravityAcc();
    for (int i=0; i<discret_->NodeRowMap()->NumMyElements(); ++i)
    {
      /*------------------------------------------------------------------*/
      //// gravity acc = g
      for(int dim=0; dim<3; ++dim)
      {
        (*acc)[i*3+dim] = extMulti * gravity_acc(dim);
      }
      /*------------------------------------------------------------------*/
    }
  }
  else
  {
    dserror("You are trying to apply gravity accelerations to a null pointer");
  }
}

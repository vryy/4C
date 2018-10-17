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
#include "../drt_lib/drt_discret.H"
#include "../drt_mat/particle_mat.H"
#include "../drt_mat/matpar_bundle.H"
#include "../linalg/linalg_utils.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*/
/* print particle time logo */
void PARTICLE::TimInt::Logo()
{
  IO::cout << "Welcome to Particle Time Integration " << IO::endl;
  IO::cout << "    ---                      ---     " << IO::endl;
  IO::cout << "  /     \\                  /     \\   " << IO::endl;
  IO::cout << "  |     |   ---->  <----   |     |   " << IO::endl;
  IO::cout << "  \\     /                  \\     /   " << IO::endl;
  IO::cout << "    ---                      ---     " << IO::endl;
  IO::cout << IO::endl;
}

/*----------------------------------------------------------------------*/
/* constructor */
PARTICLE::TimInt::TimInt(const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& particledynparams, const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<IO::DiscretizationWriter> output)
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
      writestate_((bool)DRT::INPUT::IntegralValue<int>(ioparams, "STRUCT_DISP")),
      writevelacc_((bool)DRT::INPUT::IntegralValue<int>(ioparams, "STRUCT_VEL_ACC")),
      writeresultsevery_(particledynparams.get<int>("RESULTSEVRY")),
      writeenergyevery_(particledynparams.get<int>("RESEVRYERGY")),
      writeparticlestatsevery_(particledynparams.get<int>("PARTICLESTATSEVRY")),
      energyfile_(Teuchos::null),
      particlestatsfile_(Teuchos::null),
      writeorientation_(false),
      kinergy_(0),
      intergy_(0),
      extergy_(0),
      linmomentum_(LINALG::Matrix<3, 1>(true)),
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

      disn_(Teuchos::null),
      veln_(Teuchos::null),
      accn_(Teuchos::null),
      angVeln_(Teuchos::null),
      angAccn_(Teuchos::null),
      radiusn_(Teuchos::null),

      fifc_(Teuchos::null),
      orient_(Teuchos::null),

      radius0_(Teuchos::null),
      radiusDot_(Teuchos::null),
      mass_(Teuchos::null),
      inertia_(Teuchos::null),
      f_structure_(Teuchos::null),

      dofmapexporter_(Teuchos::null),
      nodemapexporter_(Teuchos::null),

      radiusdistribution_(DRT::INPUT::IntegralValue<INPAR::PARTICLEOLD::RadiusDistribution>(
          particledynparams, "RADIUS_DISTRIBUTION")),
      variableradius_((bool)DRT::INPUT::IntegralValue<int>(
          DRT::Problem::Instance()->CavitationParams(), "COMPUTE_RADIUS_RP_BASED")),
      radiuschangefunct_(particledynparams.get<int>("RADIUS_CHANGE_FUNCT")),
      particle_algorithm_(Teuchos::null),
      collhandler_(Teuchos::null)
{
  // welcome user
  if ((printlogo_) and (myrank_ == 0))
  {
    Logo();
  }

  // check whether discretisation has been completed
  if (not discret_->Filled() || not actdis->HaveDofs())
  {
    dserror("Discretisation is not complete or has no dofs!");
  }

  // time state
  time_ = Teuchos::rcp(new TIMINT::TimIntMStep<double>(
      0, 0, 0.0));  // HERE SHOULD BE SOMETHING LIKE (particledynparams.get<double>("TIMEINIT"))
  dt_ = Teuchos::rcp(
      new TIMINT::TimIntMStep<double>(0, 0, particledynparams.get<double>("TIMESTEP")));
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
  if (DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat_ellipsoids) >=
      0)
    strategy_ = Teuchos::rcp(new TimIntStrategyEllipsoids(this));
  else
    strategy_ = Teuchos::rcp(new TimIntStrategySpheres(this));

  // initialize the vectors
  dis_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
  vel_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
  acc_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
  fifc_ = LINALG::CreateVector(*DofRowMapView(), true);
  mass_ = LINALG::CreateVector(*NodeRowMapView(), true);
  if (writeorientation_) orient_ = LINALG::CreateVector(*DofRowMapView());

  if (variableradius_ or radiuschangefunct_ > 0)
  {
    // initial radius of each particle for time dependent radius
    radius0_ = LINALG::CreateVector(*discret_->NodeRowMap(), true);
    // time derivative of radius of each particle for time dependent radius
    radiusDot_ = LINALG::CreateVector(*NodeRowMapView(), true);
  }

  // Apply Dirichlet BC and create dbc map extractor
  {
    dbcdofs_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
    dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
    Teuchos::ParameterList p;
    p.set("total time", (*time_)[0]);
    discret_->EvaluateDirichlet(p, (*dis_)(0), (*vel_)(0), (*acc_)(0), dbcdofs_, dbcmaps_);
  }

  // set initial fields
  SetInitialFields();

  // copy everything into the n+1 state vectors
  disn_ = Teuchos::rcp(new Epetra_Vector(*(*dis_)(0)));
  veln_ = Teuchos::rcp(new Epetra_Vector(*(*vel_)(0)));
  accn_ = Teuchos::rcp(new Epetra_Vector(*(*acc_)(0)));

  // decide whether there is particle contact
  if (particle_algorithm_->ParticleInteractionType() != INPAR::PARTICLEOLD::None)
  {
    // allocate vectors
    angVel_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
    angAcc_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));

    // copy the vectors to the (n+1) state vectors
    angVeln_ = LINALG::CreateVector(*DofRowMapView(), true);
    angAccn_ = LINALG::CreateVector(*DofRowMapView(), true);

    // create and fill inertia
    strategy_->ComputeInertia();
  }

  // output file for energy
  if (writeenergyevery_ != 0 and myrank_ == 0) AttachEnergyFile();

  // output file for particle statistics
  if (writeparticlestatsevery_ > 0 and myrank_ == 0) AttachParticleStatisticsFile();

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
  std::vector<const MAT::PAR::ParticleMat*> particlemat = particle_algorithm_->ParticleMat();


  // safety check
  if (particlemat.size() > 1) dserror("Only one particle material allowed!");

  const double initDensity = particlemat.size() == 1 ? particlemat[0]->initDensity_ : 0.0;

  strategy_->ComputeMass();

  // -----------------------------------------//
  // set initial radius condition if existing
  // -----------------------------------------//

  std::vector<DRT::Condition*> condition;
  discret_->GetCondition("InitialParticleRadius", condition);

  // loop over conditions
  for (size_t i = 0; i < condition.size(); ++i)
  {
    double scalar = condition[i]->GetDouble("SCALAR");
    int funct_num = condition[i]->GetInt("FUNCT");

    const std::vector<int>* nodeids = condition[i]->Nodes();
    // loop over particles in current condition
    for (size_t counter = 0; counter < (*nodeids).size(); ++counter)
    {
      int lid = discret_->NodeRowMap()->LID((*nodeids)[counter]);
      if (lid != -1)
      {
        DRT::Node* currparticle = discret_->gNode((*nodeids)[counter]);
        double function_value =
            DRT::Problem::Instance()->Funct(funct_num - 1).Evaluate(0, currparticle->X(), 0.0);
        double r_p = (*(*radius_)(0))[lid];
        r_p *= function_value * scalar;
        (*(*radius_)(0))[lid] = r_p;
        if (r_p <= 0.0) dserror("negative initial radius");

        // mass-vector: m = rho * 4/3 * PI * r^3
        (*mass_)[lid] = initDensity * PARTICLE::Utils::Radius2Volume(r_p);
      }
    }
  }

  // -----------------------------------------//
  // evaluate random normal distribution for particle radii if applicable
  // -----------------------------------------//

  switch (radiusdistribution_)
  {
    case INPAR::PARTICLEOLD::radiusdistribution_none:
    {
      // do nothing
      break;
    }
    case INPAR::PARTICLEOLD::radiusdistribution_lognormal:
    case INPAR::PARTICLEOLD::radiusdistribution_normal:
    {
      // get minimum and maximum radius for particles
      const double min_radius =
          DRT::Problem::Instance()->ParticleParamsOld().get<double>("MIN_RADIUS");
      const double max_radius =
          DRT::Problem::Instance()->ParticleParamsOld().get<double>("MAX_RADIUS");

      // loop over all particles
      for (int n = 0; n < discret_->NumMyRowNodes(); ++n)
      {
        // get local ID of current particle
        const int lid = discret_->NodeRowMap()->LID(discret_->lRowNode(n)->Id());

        // initialize random number generator with current particle radius or its natural
        // logarithm as mean and input parameter value as standard deviation
        DRT::Problem::Instance()->Random()->SetMeanVariance(
            radiusdistribution_ == INPAR::PARTICLEOLD::radiusdistribution_lognormal
                ? log((*(*radius_)(0))[lid])
                : (*(*radius_)(0))[lid],
            DRT::Problem::Instance()->ParticleParamsOld().get<double>("RADIUS_DISTRIBUTION_SIGMA"));

        // generate normally or log-normally distributed random value for particle radius
        double random_radius =
            radiusdistribution_ == INPAR::PARTICLEOLD::radiusdistribution_lognormal
                ? exp(DRT::Problem::Instance()->Random()->Normal())
                : DRT::Problem::Instance()->Random()->Normal();

        // check whether random value lies within allowed bounds, and adjust otherwise
        if (random_radius > max_radius)
          random_radius = max_radius;
        else if (random_radius < min_radius)
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

  // safety check
  if (particlemat.size() == 2 and particlemat[0]->initRadius_ != particlemat[1]->initRadius_)
    dserror("Both particle materials need to be defined with equal initRadius_!");

  const double amplitude =
      DRT::Problem::Instance()->ParticleParamsOld().get<double>("RANDOM_AMPLITUDE");
  const double initRadius = particlemat.size() > 0 ? particlemat[0]->initRadius_ : 0.0;

  for (int n = 0; n < discret_->NumMyRowNodes(); ++n)
  {
    DRT::Node* actnode = discret_->lRowNode(n);
    // get the first gid of a node and convert it into a LID
    int gid = discret_->Dof(actnode, 0);
    int lid = discret_->DofRowMap()->LID(gid);
    for (int dim = 0; dim < 3; ++dim)
    {
      if (amplitude)
      {
        double randomValue = DRT::Problem::Instance()->Random()->Uni();
        (*(*dis_)(0))[lid + dim] = actnode->X()[dim] + randomValue * amplitude * initRadius;
      }
      else
      {
        (*(*dis_)(0))[lid + dim] = actnode->X()[dim];
      }
    }
  }

  // -----------------------------------------//
  // initialize orientation field
  // -----------------------------------------//
  if (writeorientation_) strategy_->SetInitialOrientation();

  // -----------------------------------------//
  // set initial velocity field if existing
  // -----------------------------------------//

  const std::string field = "Velocity";
  std::vector<int> localdofs;
  localdofs.push_back(0);
  localdofs.push_back(1);
  localdofs.push_back(2);
  discret_->EvaluateInitialField(field, (*vel_)(0), localdofs);

  // set vector of initial particle radii if necessary
  if (radiuschangefunct_ > 0) radius0_->Update(1., *(*radius_)(0), 0.);

  return;
}

/*----------------------------------------------------------------------*/
/* prepare time step and apply Dirichlet boundary conditions */
void PARTICLE::TimInt::PrepareTimeStep()
{
  // Update map containing Dirichlet DOFs if existing
  if (dbcmaps_ != Teuchos::null && dbcmaps_->CondMap()->NumGlobalElements() != 0)
  {
    // apply Dirichlet BC and rebuild map extractor
    ApplyDirichletBC(timen_, disn_, veln_, accn_, true);
  }

  // update particle radii if necessary
  if (radiuschangefunct_ > 0)
    (*radius_)(0)->Update(
        DRT::Problem::Instance()->Funct(radiuschangefunct_ - 1).EvaluateTime(timen_), *radius0_,
        0.);

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
void PARTICLE::TimInt::ComputeAcc(Teuchos::RCP<Epetra_Vector> f_contact,
    Teuchos::RCP<Epetra_Vector> m_contact, Teuchos::RCP<Epetra_Vector> global_acc,
    Teuchos::RCP<Epetra_Vector> global_angAcc)
{
  int numrownodes = discret_->NodeRowMap()->NumMyElements();

  // in case of contact, consider corresponding forces and moments
  if (f_contact != Teuchos::null)
  {
    // sum all forces (contact and external)
    fifc_->Update(1.0, *f_contact, 1.0);

    // zero out non-planar entries in case of 2D
    if (particle_algorithm_->BinStrategy()->ParticleDim() == INPAR::PARTICLEOLD::particle_2Dz)
    {
      for (int i = 0; i < numrownodes; ++i)
      {
        (*m_contact)[i * 3 + 0] = 0.0;
        (*m_contact)[i * 3 + 1] = 0.0;
      }
    }

    // compute angular acceleration
    if (global_angAcc != Teuchos::null)
      strategy_->ComputeAngularAcceleration(*global_angAcc, *m_contact);
  }

  // zero out non-planar entries in case of 2D
  if (particle_algorithm_->BinStrategy()->ParticleDim() == INPAR::PARTICLEOLD::particle_2Dz)
  {
    for (int i = 0; i < numrownodes; ++i) (*fifc_)[i * 3 + 2] = 0.0;
  }

  // update of translational acceleration
  for (int i = 0; i < numrownodes; ++i)
  {
    const double invmass = 1.0 / (*mass_)[i];
    for (int dim = 0; dim < 3; ++dim)
    {
      (*global_acc)[i * 3 + dim] = invmass * (*fifc_)[i * 3 + dim];
    }
  }

  return;
}

/*---------------------------------------------------------------*/
/* Apply Dirichlet boundary conditions on provided state vectors */
void PARTICLE::TimInt::ApplyDirichletBC(const double time, Teuchos::RCP<Epetra_Vector> dis,
    Teuchos::RCP<Epetra_Vector> vel, Teuchos::RCP<Epetra_Vector> acc, bool recreatemap)
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
  step_ = stepn_;              // n := n+1
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
  UpdateExportersIfNecessary(mass_->Map(), (*dis_)(0)->Map());

  UpdateStateVectorMap(dis_);
  UpdateStateVectorMap(vel_);
  UpdateStateVectorMap(acc_);
  UpdateStateVectorMap(angVel_);
  UpdateStateVectorMap(angAcc_);
  strategy_->UpdateRadiusVectorMap();
  UpdateStateVectorMap(disn_);
  UpdateStateVectorMap(veln_);
  UpdateStateVectorMap(accn_);
  UpdateStateVectorMap(angVeln_);
  UpdateStateVectorMap(angAccn_);
  UpdateStateVectorMap(radiusn_, true);

  UpdateStateVectorMap(fifc_);
  UpdateStateVectorMap(orient_);

  UpdateStateVectorMap(radius0_, true);
  UpdateStateVectorMap(radiusDot_, true);
  UpdateStateVectorMap(mass_, true);
  strategy_->UpdateInertiaVectorMap();
}

/*----------------------------------------------------------------------*/
/* Read and set restart values */
void PARTICLE::TimInt::ReadRestart(const int step)
{
  IO::DiscretizationReader reader(discret_, step);
  if (step != reader.ReadInt("step")) dserror("Time step on file not equal to given step");

  restart_ = step;
  step_ = step;
  stepn_ = step_ + 1;
  time_ = Teuchos::rcp(new TIMINT::TimIntMStep<double>(0, 0, reader.ReadDouble("time")));
  timen_ = (*time_)[0] + (*dt_)[0];

  ReadRestartState();

  // short screen output
  if (myrank_ == 0)
    IO::cout << "====== Restart of the particle simulation from step " << step_ << IO::endl;

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
  if (disn_->GlobalLength() == 0) return;

  // now finish the displacement and read the remaining state vectors
  dis_->UpdateSteps(*disn_);
  reader.ReadVector(veln_, "velocity");
  vel_->UpdateSteps(*veln_);

#ifndef PARTICLE_NORESTARTACC
  reader.ReadVector(accn_, "acceleration");
  acc_->UpdateSteps(*accn_);
#endif

  reader.ReadVector(mass_, "mass");


  if (DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat_ellipsoids) <
      0)
  {
    // create a dummy vector to extract the radius vector (radiusn_ does not exist)
    Teuchos::RCP<Epetra_Vector> radius = LINALG::CreateVector(*discret_->NodeRowMap(), true);
    reader.ReadVector(radius, "radius");
    radius_->UpdateSteps(*radius);
  }

  // read in particle collision relevant data
  if (collhandler_ != Teuchos::null)
  {
    // initialize inertia
    strategy_->ComputeInertia();

    reader.ReadVector(angVeln_, "ang_velocity");
    angVel_->UpdateSteps(*angVeln_);
    reader.ReadVector(angAccn_, "ang_acceleration");
    angAcc_->UpdateSteps(*angAccn_);
    if (writeorientation_) reader.ReadVector(orient_, "orientation");
  }

  // read in variable radius relevant data
  if (variableradius_ == true)
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
  if ((writerestartevery_ and ((step_ - restart_) % writerestartevery_ == 0)) or
      forced_writerestart)
  {
    OutputRestart(datawritten);
  }

  // output results (not necessary if restart in same step)
  if (writestate_ and writeresultsevery_ and ((step_ - restart_) % writeresultsevery_ == 0) and
      (not datawritten))
  {
    OutputState(datawritten);
  }

  // output energy
  if (writeenergyevery_ and ((step_ - restart_) % writeenergyevery_ == 0))
  {
    OutputEnergy();
  }

  // output particle statistics
  if (writeparticlestatsevery_ and (step_ - restart_) % writeparticlestatsevery_ == 0)
    OutputParticleStatistics();

  return;
}

/*----------------------------------------------------------------------*/
/* write restart */
void PARTICLE::TimInt::OutputRestart(bool& datawritten)
{
  // Yes, we are going to write...
  datawritten = true;

  // mesh is written to disc
  output_->ParticleOutput(step_, (*time_)[0], true);
  output_->NewStep(step_, (*time_)[0]);

  OutputDisplacement();
  WriteVector("velocity", vel_);
  WriteVector("acceleration", acc_);

  WriteVector("radius", radius_, false);
  WriteVector("mass", mass_, false);

  if (variableradius_)
  {
    WriteVector("radius0", radius0_, false);
    WriteVector("radiusDot", radiusDot_, false);
  }

  if (collhandler_ != Teuchos::null)
  {
    if (angVeln_ != Teuchos::null)
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
  if ((myrank_ == 0) and printscreen_ and ((step_ - restart_) % printscreen_ == 0))
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
void PARTICLE::TimInt::OutputState(bool& datawritten)
{
  // Yes, we are going to write...
  datawritten = true;

  // mesh is not written to disc, only maximum node id is important for output
  output_->ParticleOutput(step_, (*time_)[0], false);
  output_->NewStep(step_, (*time_)[0]);

  OutputDisplacement();
  WriteVector("velocity", vel_);
  WriteVector("radius", radius_, false);

  if (collhandler_ != Teuchos::null)
  {
    if (angVeln_ != Teuchos::null)
    {
      WriteVector("ang_velocity", (*angVel_)(0));
      WriteVector("ang_acceleration", (*angAcc_)(0));
    }
  }

  if (writevelacc_)
  {
    WriteVector("acceleration", acc_);
  }

  if (collhandler_ != Teuchos::null) strategy_->OutputOrientation();

  // maps are rebuild in every step so that reuse is not possible
  // keeps memory usage bounded
  output_->ClearMapCache();

  return;
}

/*----------------------------------------------------------------------*/
/* Calculation of internal, external and kinetic energy */
void PARTICLE::TimInt::DetermineEnergy()
{
  if (writeenergyevery_ and (stepn_ % writeenergyevery_ == 0))
  {
    LINALG::Matrix<3, 1> gravity_acc = particle_algorithm_->GetGravityAcc(timen_);

    int numrownodes = discret_->NodeRowMap()->NumMyElements();

    // energy for all cases
    if (collhandler_ != Teuchos::null)
    {
      // total kinetic energy
      kinergy_ = strategy_->ComputeKineticEnergy();

      for (int i = 0; i < numrownodes; ++i)
      {
        double specific_energy = 0.0;

        for (int dim = 0; dim < 3; ++dim)
          specific_energy -= gravity_acc(dim) * (*disn_)[i * 3 + dim];

        intergy_ += (*mass_)[i] * specific_energy;
      }
      // total external energy not available
      extergy_ = 0.0;
    }
    else
      dserror("Energy output only possible for DEM (with collhandler_ == true)");

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
    // energy for all cases
    if (collhandler_ != Teuchos::null)
    {
      *energyfile_ << " " << std::setw(9) << step_ << std::scientific << std::setprecision(16)
                   << " " << (*time_)[0] << " " << totergy << " " << kinergy_ << " " << intergy_
                   << " " << extergy_ << " " << collhandler_->GetMaxPenetration() << std::endl;
    }
    else
      dserror("Energy output only possible for DEM (with collhandler_ == true)");
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
  for (int i = 0; i < (*radius_)(0)->MyLength(); ++i)
    myvolume += 4. / 3. * M_PI * pow((*(*radius_)(0))[i], 3);
  discret_->Comm().SumAll(&myvolume, &globalvolume, 1);

  // determine minimum and maximum particle radius
  double r_min(0.), r_max(0.);
  (*radius_)(0)->MinValue(&r_min);
  (*radius_)(0)->MaxValue(&r_max);

  // output only performed by first processor
  if (myrank_ == 0)
    *particlestatsfile_ << std::setw(10) << step_ << std::scientific << std::setprecision(10) << " "
                        << (*time_)[0] << " " << discret_->NodeRowMap()->NumGlobalElements() << " "
                        << r_min << " " << r_max << " " << globalvolume << std::endl;

  return;
}


/*----------------------------------------------------------------------*/
/* Attach file handle for energy file #energyfile_                      */
void PARTICLE::TimInt::AttachEnergyFile()
{
  if (energyfile_.is_null())
  {
    // energy for all cases
    std::string energyname =
        DRT::Problem::Instance()->OutputControlFile()->FileName() + "_particle.energy";
    energyfile_ = Teuchos::rcp(new std::ofstream(energyname.c_str()));
    (*energyfile_) << "# timestep time total_energy"
                   << " kinetic_energy internal_energy external_energy max_particle_penetration"
                   << std::endl;
  }

  return;
}


/*----------------------------------------------------------------------*
 | attach file handle for particle statistics                fang 04/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimInt::AttachParticleStatisticsFile()
{
  // create file if not yet existent
  if (particlestatsfile_ == Teuchos::null and myrank_ == 0)
  {
    // set file name
    const std::string filename(
        DRT::Problem::Instance()->OutputControlFile()->FileName() + "_particle.statistics.csv");

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
const Epetra_Map* PARTICLE::TimInt::DofRowMapView() { return discret_->DofRowMap(); }

/*----------------------------------------------------------------------*/
/* node map of particles                                                */
Teuchos::RCP<const Epetra_Map> PARTICLE::TimInt::NodeRowMap()
{
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  return Teuchos::rcp(new Epetra_Map(*noderowmap));
}

/*----------------------------------------------------------------------*/
/* view of node map of particles                                        */
const Epetra_Map* PARTICLE::TimInt::NodeRowMapView() { return discret_->NodeRowMap(); }


/*-----------------------------------------------------------------------------*/
/* Update exporter objects if layout has changed                               */
void PARTICLE::TimInt::UpdateExportersIfNecessary(
    const Epetra_BlockMap& oldnodemap, const Epetra_BlockMap& olddofmap)
{
  if (nodemapexporter_ == Teuchos::null or
      not nodemapexporter_->TargetMap().SameAs(*NodeRowMapView()) or
      not nodemapexporter_->SourceMap().SameAs(oldnodemap) or dofmapexporter_ == Teuchos::null or
      not dofmapexporter_->TargetMap().SameAs(*DofRowMapView()) or
      not dofmapexporter_->SourceMap().SameAs(olddofmap))
  {
    nodemapexporter_ = Teuchos::rcp(new Epetra_Export(oldnodemap, *NodeRowMapView()));
    dofmapexporter_ = Teuchos::rcp(new Epetra_Export(olddofmap, *DofRowMapView()));
  }

  return;
}

/*-----------------------------------------------------------------------------*/
/* Update TimIntMStep state vector with the new (appropriate) map from discret_*/
void PARTICLE::TimInt::UpdateStateVectorMap(
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector>>& stateVector, bool trg_nodeVectorType)
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
      if (err) dserror("Export using exporter returned err=%d", err);
    }
    else
    {
      // create new vector
      stateVector->ReplaceMaps(DofRowMapView());

      // transfer data
      int err = (*stateVector)(0)->Export(*old, *dofmapexporter_, Insert);
      if (err) dserror("Export using exporter returned err=%d", err);
    }
  }
}

/*-----------------------------------------------------------------------------*/
/* Update state vector with the new (appropriate) map from discret_*/
void PARTICLE::TimInt::UpdateStateVectorMap(
    Teuchos::RCP<Epetra_Vector>& stateVector, bool trg_nodeVectorType)
{
  if (stateVector != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> old = stateVector;

    if (trg_nodeVectorType)
    {
      // create new vector
      stateVector = LINALG::CreateVector(*discret_->NodeRowMap(), true);

      // transfer data
      int err = stateVector->Export(*old, *nodemapexporter_, Insert);
      if (err) dserror("Export using exporter returned err=%d", err);
    }
    else
    {
      // create new vector
      stateVector = LINALG::CreateVector(*discret_->DofRowMap(), true);

      // transfer data
      int err = stateVector->Export(*old, *dofmapexporter_, Insert);
      if (err) dserror("Export using exporter returned err=%d", err);
    }
  }
}

/*----------------------------------------------------------------------*/
//! Check exixstence of the state vectors
void PARTICLE::TimInt::CheckStateVector(
    std::string vecName, const Teuchos::RCP<const Epetra_Vector> vec, bool trg_showVec)
{
  if (vec == Teuchos::null)
  {
    std::cout << "The pointer to " << vecName << " is null\n";
    std::cin.get();
    return;
  }

  std::cout << "Processor " << myrank_ << " owns " << vec->MyLength() << "/" << vec->GlobalLength()
            << " elements of " << vecName << std::endl;
  if (trg_showVec) std::cout << *vec << std::endl;
  std::cin.get();
}


/*----------------------------------------------------------------------------*
 | return maximum particle-particle or particle-wall penetration   fang 02/17 |
 *----------------------------------------------------------------------------*/
double PARTICLE::TimInt::MaximumPenetration() const
{
  return collhandler_ == Teuchos::null ? 0. : collhandler_->GetMaxPenetration();
}

/*----------------------------------------------------------------*
 | return maximum particle-particle penetration       meier 12/17 |
 *----------------------------------------------------------------*/
double PARTICLE::TimInt::MaximumPenetrationParticle() const
{
  return collhandler_ == Teuchos::null ? 0. : collhandler_->GetMaxPenetrationParticle();
}

/*------------------------------------------------------------*
 | return maximum particle-wall penetration       meier 12/17 |
 *------------------------------------------------------------*/
double PARTICLE::TimInt::MaximumPenetrationWall() const
{
  return collhandler_ == Teuchos::null ? 0. : collhandler_->GetMaxPenetrationWall();
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

  if (collhandler_ != Teuchos::null)
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
// wrapper. On top of the output_->WriteVector() it checks that the pointer is not null. In case, it
// does not write
void PARTICLE::TimInt::WriteVector(
    const std::string name, Teuchos::RCP<Epetra_Vector> vec, const bool isdof) const
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
// wrapper. On top of the output_->WriteVector() it checks that the pointer is not null. In case, it
// does not write
void PARTICLE::TimInt::WriteVector(const std::string name,
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector>> vec, const bool isdof) const
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

  if (fifc_ != Teuchos::null)
  {
    fifc_->PutScalar(0.0);
    GravityForces(fifc_);
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
    LINALG::Matrix<3, 1> gravity_acc = particle_algorithm_->GetGravityAcc(-1.0);
    for (int i = 0; i < discret_->NodeRowMap()->NumMyElements(); ++i)
    {
      /*------------------------------------------------------------------*/
      //// gravity acc = mass_p * g
      for (int dim = 0; dim < 3; ++dim)
      {
        (*force)[i * 3 + dim] = (*mass_)[i] * gravity_acc(dim);
      }
      /*------------------------------------------------------------------*/
    }
  }
  else
  {
    dserror("You are trying to apply gravity forces to a null pointer");
  }
}

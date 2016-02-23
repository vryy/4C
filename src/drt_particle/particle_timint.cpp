/*----------------------------------------------------------------------*/
/*!
\file particle_timint.cpp
\brief Time integration for particle dynamics

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint.H"
#include "particle_algorithm.H"
#include "particle_contact.H"
#include "particle_resulttest.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/particle_mat.H"
#include "../drt_mat/matpar_bundle.H"

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
  writeorientation_(true),
  time_(Teuchos::null),
  timen_(0.0),
  dt_(Teuchos::null),
  timemax_(particledynparams.get<double>("MAXTIME")),
  stepmax_(particledynparams.get<int>("NUMSTEP")),
  step_(0),
  stepn_(0),
  dis_(Teuchos::null),
  vel_(Teuchos::null),
  acc_(Teuchos::null),
  disn_(Teuchos::null),
  veln_(Teuchos::null),
  accn_(Teuchos::null),
  radius_(Teuchos::null),
  radius0_(Teuchos::null),
  radiusdot_(Teuchos::null),
  mass_(Teuchos::null),
  inertia_(Teuchos::null),
  ang_vel_(Teuchos::null),
  ang_acc_(Teuchos::null),
  ang_veln_(Teuchos::null),
  ang_accn_(Teuchos::null),
  orient_(Teuchos::null),
  density_(-1.0),
  fifc_(Teuchos::null),
  variableradius_((bool)DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->CavitationParams(),"COMPUTE_RADIUS_RP_BASED")),
  collhandler_(Teuchos::null)
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

  // output file for energy
  if ( (writeenergyevery_ != 0) and (myrank_ == 0) )
    AttachEnergyFile();

  return;
}

/*----------------------------------------------------------------------*/
/* initialization of time integration */
void PARTICLE::TimInt::Init()
{
  Teuchos::RCP<Epetra_Vector> zeros = LINALG::CreateVector(*DofRowMapView(), true);

   // displacements D_{n}
   dis_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
   // velocities V_{n}
   vel_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
   // accelerations A_{n}
   acc_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));

   // create empty interface force vector
   fifc_ = LINALG::CreateVector(*DofRowMapView(), true);
   // radius of each particle
   radius_  = LINALG::CreateVector(*discret_->NodeRowMap(), true);
   if(variableradius_)
   {
     // initial radius of each particle for time dependent radius
     radius0_  = LINALG::CreateVector(*discret_->NodeRowMap(), true);
     // time derivative of radius of each particle for time dependent radius
     radiusdot_  = LINALG::CreateVector(*discret_->NodeRowMap(), true);
   }
   // mass of each particle
   mass_  = LINALG::CreateVector(*discret_->NodeRowMap(), true);

   // set initial fields
   SetInitialFields();

   // Apply Dirichlet BC and create dbc map extractor
   {
     dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
     Teuchos::ParameterList p;
     p.set("total time", (*time_)[0]);
     discret_->EvaluateDirichlet(p, (*dis_)(0), (*vel_)(0), (*acc_)(0), Teuchos::null, dbcmaps_);
   }

   // displacements D_{n+1} at t_{n+1}
   disn_ = Teuchos::rcp(new Epetra_Vector(*(*dis_)(0)));
   // velocities V_{n+1} at t_{n+1}
   veln_ = Teuchos::rcp(new Epetra_Vector(*(*vel_)(0)));
   // accelerations A_{n+1} at t_{n+1}
   accn_ = Teuchos::rcp(new Epetra_Vector(*(*acc_)(0)));

  return;
}

/*----------------------------------------------------------------------*/
/* Set intitial fields in structure (e.g. initial velocities */
void PARTICLE::TimInt::SetInitialFields()
{
  // make sure that a particle material is defined in the dat-file
  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat);
  if (id==-1)
    dserror("Could not find particle material");

  const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
  const MAT::PAR::ParticleMat* actmat = static_cast<const MAT::PAR::ParticleMat*>(mat);
  // all particles have identical density and initially the same radius
  density_ = actmat->density_;
  double initial_radius = actmat->initialradius_;

  double amplitude = DRT::Problem::Instance()->ParticleParams().get<double>("RANDOM_AMPLITUDE");

  // initialize displacement field, radius and mass vector
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
        double randomwert = DRT::Problem::Instance()->Random()->Uni();
        (*(*dis_)(0))[lid+dim] = actnode->X()[dim] + randomwert * amplitude * initial_radius;
      }
      else
      {
        (*(*dis_)(0))[lid+dim] = actnode->X()[dim];
      }
    }

    // initialize radii of particles
    (*radius_)[n] = initial_radius;

    // mass-vector: m = rho * 4/3 * PI *r^3
    (*mass_)[n] = density_ * 4.0/3.0 * M_PI * initial_radius * initial_radius * initial_radius;
  }

  // set initial radius condition if existing
  std::vector<DRT::Condition*> condition;
  discret_->GetCondition("InitialParticleRadius", condition);

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
        double r_p = (*radius_)[lid];
        r_p *= function_value * scalar;
        (*radius_)[lid] = r_p;
        if(r_p <= 0.0)
          dserror("negative initial radius");

        // mass-vector: m = rho * 4/3 * PI * r^3
        (*mass_)[lid] = density_ * 4.0/3.0 * M_PI * r_p * r_p * r_p;
      }
    }
  }

  // evaluate random normal distribution for particle radii if applicable
  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"RADIUS_DISTRIBUTION"))
  {
    // get minimum and maximum radius for particles
    const double min_radius = DRT::Problem::Instance()->ParticleParams().get<double>("MIN_RADIUS");
    const double max_radius = DRT::Problem::Instance()->ParticleParams().get<double>("MAX_RADIUS");

    // loop over all particles
    for(int n=0; n<discret_->NumMyRowNodes(); ++n)
    {
      // get local ID of current particle
      int lid = discret_->NodeRowMap()->LID(discret_->lRowNode(n)->Id());

      if(lid != -1)
      {
        // provide random number generator with local ID of current particle as deterministic seed to ensure reproducibility of simulation
        DRT::Problem::Instance()->Random()->SetRandSeed(lid);

        // initialize random number generator with current particle radius as mean and input parameter value as standard deviation
        DRT::Problem::Instance()->Random()->SetMeanVariance((*radius_)[lid],DRT::Problem::Instance()->ParticleParams().get<double>("RADIUS_DISTRIBUTION_SIGMA"));

        // generate normally distributed random value for particle radius
        double random_radius = DRT::Problem::Instance()->Random()->Normal();

        // check whether random value lies within allowed bounds, and adjust otherwise
        if(random_radius > max_radius)
          random_radius = max_radius;
        else if(random_radius < min_radius)
          random_radius = min_radius;

        // set particle radius to random value
        (*radius_)[lid] = random_radius;

        // recompute particle mass
        (*mass_)[lid] = density_*4./3.*M_PI*random_radius*random_radius*random_radius;
      }
    }
  }

  // set initial velocity field if existing
  const std::string field = "Velocity";
  std::vector<int> localdofs;
  localdofs.push_back(0);
  localdofs.push_back(1);
  localdofs.push_back(2);
  discret_->EvaluateInitialField(field,(*vel_)(0),localdofs);

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
  Teuchos::RCP<Epetra_Vector> global_ang_acc)
{
  int numrownodes = discret_->NodeRowMap()->NumMyElements();
  // in case of contact, consider corresponding forces and moments
  if(f_contact != Teuchos::null)
  {
    // sum all forces (contact and external)
    fifc_->Update(1.0, *f_contact, 1.0);

    for(int i=0; i<numrownodes; ++i)
    {
      double invinertia = 1.0/(*inertia_)[i];
      for(int dim=0; dim<3; ++dim)
        (*global_ang_acc)[i*3+dim] = invinertia * (*m_contact)[i*3+dim];
    }
  }

  // update of translational acceleration
  for(int i=0; i<numrownodes; ++i)
  {
    double invmass = 1.0/(*mass_)[i];
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
  {
    discret_->EvaluateDirichlet(p, dis, vel, acc,
                                Teuchos::null, dbcmaps_);
  }
  else
  {
    discret_->EvaluateDirichlet(p, dis, vel, acc,
                               Teuchos::null, Teuchos::null);
  }
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
  Teuchos::RCP<Epetra_Vector> old;

  if (disn_ != Teuchos::null)
  {
    old = disn_;
    disn_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
    LINALG::Export(*old, *disn_);
  }

  if (veln_ != Teuchos::null)
  {
    old = veln_;
    veln_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
    LINALG::Export(*old, *veln_);
  }

  if (ang_veln_ != Teuchos::null)
  {
    old = ang_veln_;
    ang_veln_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
    LINALG::Export(*old, *ang_veln_);
  }

  if (accn_ != Teuchos::null)
  {
    old = accn_;
    accn_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
    LINALG::Export(*old, *accn_);
  }

  if (ang_accn_ != Teuchos::null)
  {
    old = ang_accn_;
    ang_accn_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
    LINALG::Export(*old, *ang_accn_);
  }

  if ((*dis_)(0) != Teuchos::null)
  {
    const Teuchos::RCP<Epetra_Vector> oldvec = Teuchos::rcp(new Epetra_Vector(*(*dis_)(0)));
    dis_->ReplaceMaps(discret_->DofRowMap());
    LINALG::Export(*oldvec, *(*dis_)(0));
  }

  if ((*vel_)(0) != Teuchos::null)
  {
    const Teuchos::RCP<Epetra_Vector> oldvec = Teuchos::rcp(new Epetra_Vector(*(*vel_)(0)));
    vel_->ReplaceMaps(discret_->DofRowMap());
    LINALG::Export(*oldvec, *(*vel_)(0));
  }

  if (ang_vel_ != Teuchos::null and (*ang_vel_)(0) != Teuchos::null)
  {
    const Teuchos::RCP<Epetra_Vector> oldvec = Teuchos::rcp(new Epetra_Vector(*(*ang_vel_)(0)));
    ang_vel_->ReplaceMaps(discret_->DofRowMap());
    LINALG::Export(*oldvec, *(*ang_vel_)(0));
  }

  if ((*acc_)(0) != Teuchos::null)
  {
    const Teuchos::RCP<Epetra_Vector> oldvec = Teuchos::rcp(new Epetra_Vector(*(*acc_)(0)));
    acc_->ReplaceMaps(discret_->DofRowMap());
    LINALG::Export(*oldvec, *(*acc_)(0));
  }

  if (ang_acc_ != Teuchos::null and (*ang_acc_)(0) != Teuchos::null)
  {
    const Teuchos::RCP<Epetra_Vector> oldvec = Teuchos::rcp(new Epetra_Vector(*(*ang_acc_)(0)));
    ang_acc_->ReplaceMaps(discret_->DofRowMap());
    LINALG::Export(*oldvec, *(*ang_acc_)(0));
  }
  if (orient_ != Teuchos::null)
  {
    old = orient_;
    orient_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
    LINALG::Export(*old, *orient_);
  }
  if (radius_ != Teuchos::null)
  {
    old = radius_;
    radius_ = LINALG::CreateVector(*discret_->NodeRowMap(),true);
    LINALG::Export(*old, *radius_);
  }

  if (radius0_ != Teuchos::null)
  {
    old = radius0_;
    radius0_ = LINALG::CreateVector(*discret_->NodeRowMap(),true);
    LINALG::Export(*old, *radius0_);
  }

  if (radiusdot_ != Teuchos::null)
  {
    old = radiusdot_;
    radiusdot_ = LINALG::CreateVector(*discret_->NodeRowMap(),true);
    LINALG::Export(*old, *radiusdot_);
  }

  if (mass_ != Teuchos::null)
  {
    old = mass_;
    mass_ = LINALG::CreateVector(*discret_->NodeRowMap(),true);
    LINALG::Export(*old, *mass_);
  }

  if (inertia_ != Teuchos::null)
  {
    old = inertia_;
    inertia_ = LINALG::CreateVector(*discret_->NodeRowMap(),true);
    LINALG::Export(*old, *inertia_);
  }

  if (fifc_ != Teuchos::null)
  {
    fifc_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
  }

  return;
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

  step_ = step;
  stepn_ = step_ + 1;
  time_ = Teuchos::rcp(new TIMINT::TimIntMStep<double>(0, 0, reader.ReadDouble("time")));
  timen_ = (*time_)[0] + (*dt_)[0];

  ReadRestartState();

  return;
}

/*----------------------------------------------------------------------*/
/* Read and set restart state */
void PARTICLE::TimInt::ReadRestartState()
{
  IO::DiscretizationReader reader(discret_, step_);
  // maps need to be adapted to restarted discretization
  UpdateStatesAfterParticleTransfer();

  // start with reading radius in order to find out whether particles exist
  reader.ReadVector(radius_, "radius");

  if(radius_->GlobalLength() != 0)
  {
    // now, remaining state vectors an be read in
    reader.ReadVector(disn_, "displacement");
    dis_->UpdateSteps(*disn_);
    reader.ReadVector(veln_, "velocity");
    vel_->UpdateSteps(*veln_);
    reader.ReadVector(accn_, "acceleration");
    acc_->UpdateSteps(*accn_);
    reader.ReadVector(mass_, "mass");

    // read in particle collision relevant data
    if(collhandler_ != Teuchos::null)
    {
      // initialize inertia
      for(int lid=0; lid<discret_->NumMyRowNodes(); ++lid)
      {
        const double rad = (*radius_)[lid];
        // inertia-vector: sphere: I = 2/5 * m * r^2
        (*inertia_)[lid] = 0.4 * (*mass_)[lid] * rad * rad;
      }

      reader.ReadVector(ang_veln_, "ang_velocity");
      ang_vel_->UpdateSteps(*ang_veln_);
      reader.ReadVector(ang_accn_, "ang_acceleration");
      ang_acc_->UpdateSteps(*ang_accn_);
      if(writeorientation_)
        reader.ReadVector(orient_, "orientation");
    }

    // read in variable radius relevant data
    if(variableradius_ == true)
    {
      reader.ReadVector(radius0_, "radius0");
      reader.ReadVector(radiusdot_, "radiusdot");
    }
  }

  return;
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
  if ( (writerestartevery_ and (step_%writerestartevery_ == 0)) or forced_writerestart )
  {
    OutputRestart(datawritten);
  }

  // output results (not necessary if restart in same step)
  if ( writestate_
       and writeresultsevery_ and (step_%writeresultsevery_ == 0)
       and (not datawritten) )
  {
    OutputState(datawritten);
  }

  // output energy
  if ( writeenergyevery_ and (step_%writeenergyevery_ == 0) )
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
  output_->WriteVector("displacement", (*dis_)(0));
  output_->WriteVector("velocity", (*vel_)(0));
  output_->WriteVector("acceleration", (*acc_)(0));
  output_->WriteVector("radius", radius_, output_->nodevector);
  output_->WriteVector("mass", mass_, output_->nodevector);
  if(variableradius_)
  {
    output_->WriteVector("radius0", radius0_, output_->nodevector);
    output_->WriteVector("radiusdot", radiusdot_, output_->nodevector);
  }

  if(collhandler_ != Teuchos::null)
  {
    if(ang_veln_ != Teuchos::null)
    {
      output_->WriteVector("ang_velocity", (*ang_vel_)(0));
      output_->WriteVector("ang_acceleration", (*ang_acc_)(0));
    }

    if(writeorientation_)
      output_->WriteVector("orientation", orient_);
  }

  // maps are rebuild in every step so that reuse is not possible
  // keeps memory usage bounded
  output_->ClearMapCache();

  // info dedicated to user's eyes staring at standard out
  if ( (myrank_ == 0) and printscreen_ and (step_%printscreen_==0))
  {
    printf("====== Restart written in step %d\n", step_);
    fflush(stdout);
  }

  // info dedicated to processor error file
  if (printerrfile_)
  {
    fprintf(errfile_, "====== Restart written in step %d\n", step_);
    fflush(errfile_);
  }

  return;
}

/*----------------------------------------------------------------------*/
/* output displacements, velocities and accelerations */
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
  output_->WriteVector("displacement", (*dis_)(0));
  output_->WriteVector("velocity", (*vel_)(0));
  if(writevelacc_)
  {
    output_->WriteVector("acceleration", (*acc_)(0));
  }

  output_->WriteVector("radius", radius_, output_->nodevector);
  if(collhandler_ != Teuchos::null and writeorientation_)
  {
    output_->WriteVector("orientation", orient_);
  }

  // maps are rebuild in every step so that reuse is not possible
  // keeps memory usage bounded
  output_->ClearMapCache();

  return;
}

/*----------------------------------------------------------------------*/
/* Calculation of internal, external and kinetic energy */
void PARTICLE::TimInt::DetermineEnergy()
{
  if ( writeenergyevery_ and (stepn_%writeenergyevery_ == 0) and collhandler_ != Teuchos::null)
  {
    LINALG::Matrix<3,1> gravity_acc = particle_algorithm_->GetGravityAcc();

    // total kinetic energy
    kinergy_ = 0.0;

    int numrownodes = discret_->NodeRowMap()->NumMyElements();
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
        rot_energy += pow((*ang_veln_)[i*3+dim], 2.0);
      }

      intergy_ += (*mass_)[i] * specific_energy;
      kinergy_ += 0.5 * ((*mass_)[i] * kinetic_energy + (*inertia_)[i] * rot_energy);
    }

    double global_energy[2] = {0.0, 0.0};
    double energies[2] = {intergy_, kinergy_};
    discret_->Comm().SumAll(&energies[0], &global_energy[0], 2);

    intergy_ = global_energy[0];
    kinergy_ = global_energy[1];
    // total external energy not available
    extergy_ = 0.0;
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
    *energyfile_ << " " << std::setw(9) << step_
                 << std::scientific  << std::setprecision(16)
                 << " " << (*time_)[0]
                 << " " << totergy
                 << " " << kinergy_
                 << " " << intergy_
                 << " " << extergy_
                 << " " << collhandler_->GetMaxPenetration()
                 << std::endl;
  }

  return;
}

/*----------------------------------------------------------------------*/
/* Set forces due to interface loads, the force is expected external-force-like */
void PARTICLE::TimInt::SetForceInterface
(
  Teuchos::RCP<Epetra_MultiVector> iforce  ///< the force on interface
)
{
  fifc_->Update(1.0, *iforce, 0.0);
  return;
}

/*----------------------------------------------------------------------*/
/* Attach file handle for energy file #energyfile_                      */
void PARTICLE::TimInt::AttachEnergyFile()
{
  if (energyfile_.is_null())
  {
    std::string energyname
      = DRT::Problem::Instance()->OutputControlFile()->FileName()
      + "_particle.energy";
    energyfile_ = Teuchos::rcp(new std::ofstream(energyname.c_str()));
    (*energyfile_) << "# timestep time total_energy"
                   << " kinetic_energy internal_energy external_energy max_particle_penetration"
                   << std::endl;
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

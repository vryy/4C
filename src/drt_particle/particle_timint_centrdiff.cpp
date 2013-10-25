/*----------------------------------------------------------------------*/
/*!
\file particle_timint_centrdiff.cpp
\brief Particle time integration with central difference scheme 2nd order (explicit),
       also known as Velocity-Verlet algorithm

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_centrdiff.H"
#include "particle_algorithm.H"
#include "particle_contact.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_mat/particle_mat.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"

#define orientationvector

/*----------------------------------------------------------------------*/
/* Constructor */
PARTICLE::TimIntCentrDiff::TimIntCentrDiff(
    const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& sdynparams,
    const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<LINALG::Solver> contactsolver,
    Teuchos::RCP<IO::DiscretizationWriter> output
  ) : STR::TimIntCentrDiff
  (
    ioparams,
    sdynparams,
    xparams,
    actdis,
    solver,
    contactsolver,
    output,
    false // DetermineMassDampConsistAccel() must not be called
  ),
  mass_(Teuchos::null),
  inertia_(Teuchos::null),
  ang_vel_(Teuchos::null),
  ang_acc_(Teuchos::null),
  ang_veln_(Teuchos::null),
  ang_accn_(Teuchos::null),
  orient_(Teuchos::null),
  collhandler_(Teuchos::null)
{
  // allocate vectors
  radius_  = LINALG::CreateVector(*discret_->NodeRowMap(), true);
  mass_  = LINALG::CreateVector(*discret_->NodeRowMap(), true);

  // make sure that a particle material is defined in the dat-file
  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat);
  if (id==-1)
    dserror("Could not find particle material");

  const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
  const MAT::PAR::ParticleMat* actmat = static_cast<const MAT::PAR::ParticleMat*>(mat);
  // currently all particles have identical density and radius
  density_ = actmat->density_;
  double initial_radius = actmat->initialradius_;

  const Teuchos::ParameterList& params = DRT::Problem::Instance()->ParticleParams();
  double amplitude = params.get<double>("RANDOM_AMPLITUDE");

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
    (*mass_)[n] = density_ * 4.0/3.0 * M_PI * pow(initial_radius, 3.0);
  }

  // DetermineMassDampConsistAccel() is called at the end of Algorithm::Init() after proper setup of the problem

  return;
}


/*----------------------------------------------------------------------*/
/* mostly init of collision handling  */
void PARTICLE::TimIntCentrDiff::Init()
{
  // decide whether there is particle contact
  const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
  INPAR::PARTICLE::ContactStrategy contact_strategy = DRT::INPUT::IntegralValue<INPAR::PARTICLE::ContactStrategy>(particleparams,"CONTACT_STRATEGY");

  if(contact_strategy != INPAR::PARTICLE::None)
  {
    collhandler_ = Teuchos::rcp(new PARTICLE::ParticleCollisionHandler(discret_, particle_algorithm_));

    // allocate vectors
    inertia_  = LINALG::CreateVector(*discret_->NodeRowMap(), true);

    ang_vel_ = Teuchos::rcp(new STR::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
    ang_acc_ = Teuchos::rcp(new STR::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));

    ang_veln_ = LINALG::CreateVector(*DofRowMapView(),true);
    ang_accn_ = LINALG::CreateVector(*DofRowMapView(),true);

#ifdef orientationvector
    //initialize orientation-vector for visualization
    orient_ = LINALG::CreateVector(*DofRowMapView(),true);
    InitializeOrientVector();
#endif

    // initialize inertia
    for(int n=0; n<discret_->NumMyRowNodes(); ++n)
    {
      double r_p = (*radius_)[n];

      //inertia-vector: sphere: I = 2/5 * m * r^2
      (*inertia_)[n] = 0.4 * (*mass_)[n] * pow(r_p, 2.0);
    }

    // energy file with correct name must be added here (to overwrite base class constructor call)
    AttachEnergyFile();
  }

  //Initial radius condition
  ApplyInitialRadiusCondition();

  return;
}


/*---------------------------------------------------------------------------*/
/* equilibrate system at initial state and identify consistent accelerations */
void PARTICLE::TimIntCentrDiff::DetermineMassDampConsistAccel()
{
  ComputeAcc(Teuchos::null, Teuchos::null, (*acc_)(0), Teuchos::null);

  return;
}


/*----------------------------------------------------------------------*/
/* Integrate step */
int PARTICLE::TimIntCentrDiff::IntegrateStep()
{
  // time this step
  timer_->ResetStartTime();

  const double dt = (*dt_)[0];   // \f$\Delta t_{n}\f$
  const double dthalf = dt/2.0;  // \f$\Delta t_{n+1/2}\f$

  // new velocities \f$V_{n+1/2}\f$
  veln_->Update(1.0, *(*vel_)(0), 0.0);
  veln_->Update(dthalf, *(*acc_)(0), 1.0);

  // new displacements \f$D_{n+1}\f$
  disn_->Update(1.0, *(*dis_)(0), 0.0);
  disn_->Update(dt, *veln_, 1.0);

  // apply Dirichlet BCs
  ApplyDirichletBC(timen_, disn_, Teuchos::null, Teuchos::null, false);
  ApplyDirichletBC(timen_-dthalf, Teuchos::null, veln_, Teuchos::null, false);

  // define vector for contact force and moment
  Teuchos::RCP<Epetra_Vector> f_contact = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> m_contact = Teuchos::null;

  // total internal energy (elastic spring and potential energy)
  intergy_ = 0.0;
  //---------------------Compute Collisions-----------------------
  if(collhandler_ != Teuchos::null)
  {
    // new angular-velocities \f$ang_V_{n+1/2}\f$
    ang_veln_->Update(1.0, *(*ang_vel_)(0), 0.0);
    ang_veln_->Update(dthalf, *(*ang_acc_)(0), 1.0);

    // initialize vectors for contact force and moment
    f_contact = LINALG::CreateVector(*(discret_->DofRowMap()),true);
    m_contact = LINALG::CreateVector(*(discret_->DofRowMap()),true);

    SetStatesForCollision();

    intergy_ = collhandler_->ComputeCollisions(dt, f_contact, m_contact);
  }
  //--------------------------------------------------------------

  ComputeAcc(f_contact, m_contact, accn_, ang_accn_);

  // update of end-velocities \f$V_{n+1}\f$
  veln_->Update(dthalf, *accn_, 1.0);
  if(collhandler_ != Teuchos::null)
  {
    ang_veln_->Update(dthalf,*ang_accn_,1.0);

#ifdef orientationvector
  // for visualization of orientation vector
  RotateOrientVector(dt);
#endif
  }

  // apply Dirichlet BCs
  ApplyDirichletBC(timen_, Teuchos::null, veln_, accn_, false);

  return 0;
}


/*----------------------------------------------------------------------*/
/* State vectors are updated according to the new distribution of particles */
void PARTICLE::TimIntCentrDiff::UpdateStatesAfterParticleTransfer()
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
#ifdef orientationvector
  if (orient_ != Teuchos::null)
  {
    old = orient_;
    orient_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
    LINALG::Export(*old, *orient_);
  }
#endif
  if (radius_ != Teuchos::null)
  {
    old = radius_;
    radius_ = LINALG::CreateVector(*discret_->NodeRowMap(),true);
    LINALG::Export(*old, *radius_);
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
/* Write restart output for particles                                   */
void PARTICLE::TimIntCentrDiff::OutputRestart
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
#ifdef energyoutput
  if(collhandler_ != Teuchos::null)
  {
    output_->WriteVector("ang_velocity", (*ang_vel_)(0));
    output_->WriteVector("ang_acceleration", (*ang_acc_)(0));
#ifdef orientationvector
    output_->WriteVector("orientation", orient_);
#endif
    collhandler_->PrintMaxPenetration(stepn_, timen_);
  }
#endif

   // maps are rebuild in every step so that reuse is not possible
   // keeps memory usage bounded
   output_->ClearMapCache();

   // info dedicated to user's eyes staring at standard out
   if ( (myrank_ == 0) and printscreen_ and (GetStep()%printscreen_==0))
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

   // we will say what we did
   return;

}


/*----------------------------------------------------------------------*/
/* output displacements, velocities, accelerations and radius           */
void PARTICLE::TimIntCentrDiff::OutputState
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
  output_->WriteVector("acceleration", (*acc_)(0));
  output_->WriteVector("radius", radius_, output_->nodevector);
  if(collhandler_ != Teuchos::null)
  {
#ifdef orientationvector
    output_->WriteVector("orientation", orient_);
#endif
    collhandler_->PrintMaxPenetration(step_, (*time_)[0]);
  }

  // maps are rebuild in every step so that reuse is not possible
  // keeps memory usage bounded
  output_->ClearMapCache();

  // leave for good
  return;
}


/*----------------------------------------------------------------------*/
/* Read and set restart state */
void PARTICLE::TimIntCentrDiff::ReadRestartState()
{
  IO::DiscretizationReader reader(discret_, step_);
  // maps need to be adapted to restarted discretization
  UpdateStatesAfterParticleTransfer();

  // now, state vectors an be read in
  reader.ReadVector(disn_, "displacement");
  dis_->UpdateSteps(*disn_);
  reader.ReadVector(veln_, "velocity");
  vel_->UpdateSteps(*veln_);
  reader.ReadVector(accn_, "acceleration");
  acc_->UpdateSteps(*accn_);
  reader.ReadVector(radius_, "radius");

  // initialize radius
  for(int lid=0; lid<discret_->NumMyRowNodes(); ++lid)
  {
    // mass-vector: m = rho * 4/3 * PI *r^3
    (*mass_)[lid] = density_ * 4.0/3.0 * M_PI * pow((*radius_)[lid], 3.0);
  }

  if(collhandler_ != Teuchos::null)
  {
    // initialize inertia
    for(int lid=0; lid<discret_->NumMyRowNodes(); ++lid)
    {
      // inertia-vector: sphere: I = 2/5 * m * r^2
      (*inertia_)[lid] = 0.4 * (*mass_)[lid] * pow((*radius_)[lid], 2.0);
    }

    reader.ReadVector(ang_veln_, "ang_velocity");
    ang_vel_->UpdateSteps(*ang_veln_);
    reader.ReadVector(ang_accn_, "ang_acceleration");
    ang_acc_->UpdateSteps(*ang_accn_);
#ifdef orientationvector
    reader.ReadVector(orient_, "orientation");
#endif
  }

  return;
}


/*----------------------------------------------------------------------*/
/* apply initial condition for particle radius */
void PARTICLE::TimIntCentrDiff::ApplyInitialRadiusCondition()
{
  std::vector<DRT::Condition*> condition;
  discret_->GetCondition("InitialParticleRadius", condition);

  //loop over conditions
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
        (*mass_)[lid] = density_ * 4.0/3.0 * M_PI * pow(r_p, 3.0);

        if(collhandler_ != Teuchos::null)
        {
          //inertia-vector: sphere: I = 2/5 * m * r^2
          (*inertia_)[lid] = 0.4 * (*mass_)[lid] * pow(r_p, 2.0);
        }
      }
    }
  }

  if(collhandler_ != Teuchos::null)
  {
    // check for validity of input data
    double maxradius = 1.0e12;
    radius_->MaxValue(&maxradius);
    if(maxradius > collhandler_->GetMaxRadius())
      dserror("Input parameter MAX_RADIUS invalid");

    double minradius = -1.0;
    radius_->MinValue(&minradius);
    if(minradius < collhandler_->GetMinRadius())
      dserror("Input parameter MIN_RADIUS invalid");
  }

  return;
}


/*----------------------------------------------------------------------*/
/* acceleration is applied from given forces */
void PARTICLE::TimIntCentrDiff::ComputeAcc(
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


/*----------------------------------------------------------------------*/
/* update step */
void PARTICLE::TimIntCentrDiff::UpdateStepState()
{
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1}, D_{n-1} := D_{n}
  dis_->UpdateSteps(*disn_);
  // new velocities at t_{n+1} -> t_n
  //    V_{n} := V_{n+1}, V_{n-1} := V_{n}
  vel_->UpdateSteps(*veln_);
  // new accelerations at t_{n+1} -> t_n
  //    A_{n} := A_{n+1}, A_{n-1} := A_{n}
  acc_->UpdateSteps(*accn_);

  if(collhandler_ != Teuchos::null)
  {
    // new angular-velocities at t_{n+1} -> t_n
    //    ang_V_{n} := ang_V_{n+1}, ang_V_{n-1} := ang_V_{n}
    ang_vel_->UpdateSteps(*ang_veln_);
    // new angular-accelerations at t_{n+1} -> t_n
    //    ang_A_{n} := ang_A_{n+1}, ang_A_{n-1} := ang_A_{n}
    ang_acc_->UpdateSteps(*ang_accn_);
  }

  // update contact and meshtying
  UpdateStepContactMeshtying();

  return;
}


/*----------------------------------------------------------------------*/
/* states are given to the collision handler */
void PARTICLE::TimIntCentrDiff::SetStatesForCollision()
{
  collhandler_->SetState(radius_, mass_);

  // dof based vectors
  discret_->SetState("bubblepos", disn_);
  discret_->SetState("bubblevel", veln_);
  discret_->SetState("bubbleangvel", ang_veln_);

  return;
}


/*----------------------------------------------------------------------*/
/* initialization of vector for visualization of the particle orientation */
void PARTICLE::TimIntCentrDiff::InitializeOrientVector()
{
  int numrownodes = discret_->NodeRowMap()->NumMyElements();
  for(int i=0; i<numrownodes; ++i)
  {
    (*orient_)[i*3] = 0.0;
    (*orient_)[i*3+1] = 0.0;
    (*orient_)[i*3+2] = 1.0;
  }

  return;
}


/*----------------------------------------------------------------------*/
/* update of vector for visualization of the particle orientation */
void PARTICLE::TimIntCentrDiff::RotateOrientVector(double dt)
{
  int numrownodes = discret_->NodeRowMap()->NumMyElements();
  for(int i=0; i<numrownodes; ++i)
  {
    double ang_vel[3];
    double r[3];

    for(int dim=0; dim<3; ++dim)
    {
      ang_vel[dim] = (*ang_veln_)[i*3+dim];
      r[dim] = (*orient_)[i*3+dim];
    }

    // visualization only valid for 2D when rotating around y-axis

    // simplified/linearized orient vector - just for visualization
    // delta_r = \Delta t * (ang_vel x r)
    (*orient_)[i*3]   += dt * (ang_vel[1] * r[2] - ang_vel[2] * r[1]);
    (*orient_)[i*3+1] += dt * (ang_vel[2] * r[0] - ang_vel[0] * r[2]);
    (*orient_)[i*3+2] += dt * (ang_vel[0] * r[1] - ang_vel[1] * r[0]);
    //--------------------------------------------------------------

    //more exactly------------------------------------------------
//  	double d[3];
//  	double norm_ang_vel=0.0;
//  	//norm ang_vel
//  	for(int dim=0; dim<3; ++dim)
//  	{
//  		norm_ang_vel += ang_vel[dim]*ang_vel[dim];
//  	}
//  	norm_ang_vel = sqrt(norm_ang_vel);
//
//  	if(norm_ang_vel > 1E-10)
//  	{
//      double scalar = 0.0;
//  	  double invnorm_ang_vel = 1.0 / norm_ang_vel;
//  		for(int dim=0; dim<3; ++dim)
//  		{
//  			d[dim] = invnorm_ang_vel * ang_vel[dim];
//  			scalar += d[dim] * r[dim];
//  		}
//
//  		for(int dim=0; dim<3; ++dim)
//  		{
//  		  (*orient_)[i*3+dim] += (1-cos(norm_ang_vel*dt)) * scalar * d[dim] - (1-cos(norm_ang_vel*dt)) * r[dim]
//  			      + sin(norm_ang_vel*dt) * ( d[(dim+1)%3] * r[(dim+2)%3] - d[(dim+2)%3] * r[(dim+1)%3] );
//  		}
//  	}
    //--------------------------------------------------------------

  }

  return;
}


/*----------------------------------------------------------------------*/
/* energy of the particle system is calculated */
void PARTICLE::TimIntCentrDiff::DetermineEnergy()
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
    // total external energy not yet computed
    extergy_ = 0.0;
  }

  return;
}

/*----------------------------------------------------------------------*/
/* Attach file handle for particle energy file #energyfile_ */
void PARTICLE::TimIntCentrDiff::AttachEnergyFile(std::string name)
{
  STR::TimInt::AttachEnergyFile("_particle");

  return;
}

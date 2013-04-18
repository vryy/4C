/*----------------------------------------------------------------------*/
/*!
\file particle_timint_centrdiff.cpp
\brief Particle time integration with central difference scheme 2nd order (explicit),
       also known as Velocity-Verlet algorithm

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-152537
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_centrdiff.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_mat/particle_mat.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"

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
  )
{
  // allocate vectors
  radius_  = LINALG::CreateVector(*discret_->NodeRowMap(), true);
  radiusn_  = LINALG::CreateVector(*discret_->NodeRowMap(), true);

  // make sure that a particle material is defined in the dat-file
  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat);
  if (id==-1)
    dserror("Could not find particle material");

  const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
  const MAT::PAR::ParticleMat* actmat = static_cast<const MAT::PAR::ParticleMat*>(mat);
  // currently all particles have identical density and radius
  density_ = actmat->density_;
  double initial_radius = actmat->initialradius_;
  // initialize displacement field with nodal positions
  for(int n=0; n<discret_->NumMyRowNodes(); n++)
  {
    DRT::Node* actnode = discret_->lRowNode(n);
    // get the first gid of a node and convert it into a LID
    int gid = discret_->Dof(actnode, 0);
    int lid = discret_->DofRowMap()->LID(gid);
    for (int dim=0; dim<3; dim++)
    {
      (*(*dis_)(0))[lid+dim] = actnode->X()[dim];
    }

    // initialize radii of particles
    (*radiusn_)[n] = initial_radius;
  }

  // DetermineMassDampConsistAccel() is called at the end of Algorithm::Init() after proper setup of the problem

  return;
}


/*----------------------------------------------------------------------*/
/* Integrate step */
void PARTICLE::TimIntCentrDiff::IntegrateStep()
{
  // time this step
  timer_->ResetStartTime();

  // get state before Dirichlet conditions are applied because ClearState is called there
  Teuchos::RCP<const Epetra_Vector> bubbleforces = discret_->GetState("particleforces");

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

  // build new external forces (This is too much work here because an assembled row vector
  // for the forces would be enough. But in later applications with particle contact
  // it is necessary to evaluate also ghost bubbles.
  // in short: correct col layout is available but only row layout is needed
  // --> LINALG::Export(col->row) from bubbleforces to fextn_
  fextn_->PutScalar(0.0);
  LINALG::Export(*bubbleforces, *fextn_);

  //check whether all bubbles are of the same size
  double maxradius=-1.0;
  double minradius=-0.5;
  radiusn_->MaxValue(&maxradius);
  radiusn_->MinValue(&minradius);
  if(abs(maxradius-minradius) > EPS10)
    dserror("Assumption of equally sized bubbles is not valid. Implementation is missing!");
  // normal case is a proc with at least a few bubbles
  if(radiusn_->MyLength() != 0)
  {
    // as long as it is valid, mass can be calculated just once
    double radius = (*radiusn_)[0];
    double mass = density_ * 4.0/3.0 * M_PI * radius * radius * radius;
    accn_->Update(1/mass, *fextn_, 0.0);
  }
  else
    accn_->PutScalar(0.0);

  // update of end-velocities \f$V_{n+1}\f$
  veln_->Update(dthalf, *accn_, 1.0);

  // apply Dirichlet BCs on accelerations
  ApplyDirichletBC(timen_, Teuchos::null, veln_, accn_, false);

  return;
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

  if (accn_ != Teuchos::null)
  {
    old = accn_;
    accn_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
    LINALG::Export(*old, *accn_);
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

  if ((*acc_)(0) != Teuchos::null)
  {
    const Teuchos::RCP<Epetra_Vector> oldvec = Teuchos::rcp(new Epetra_Vector(*(*acc_)(0)));
    acc_->ReplaceMaps(discret_->DofRowMap());
    LINALG::Export(*oldvec, *(*acc_)(0));
  }

  if (radiusn_ != Teuchos::null)
  {
    old = radiusn_;
    radiusn_ = LINALG::CreateVector(*discret_->NodeRowMap(),true);
    LINALG::Export(*old, *radiusn_);
  }

  if (fextn_ != Teuchos::null)
  {
    fextn_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
  }

  return;
}

/*---------------------------------------------------------------------------*/
/* equilibrate system at initial state and identify consistent accelerations */
void PARTICLE::TimIntCentrDiff::DetermineMassDampConsistAccel()
{
  // build new external forces (This is too much work here because an assembled row vector
  // for the forces would be enough. But in later applications with particle contact
  // it is necessary to evaluate also ghost bubbles.
  // in short: correct col layout is available but only row layout is needed
  // --> LINALG::Export(col->row) from bubbleforces to fextn_
  Teuchos::RCP<const Epetra_Vector> bubbleforces = discret_->GetState("particleforces");
  fextn_->PutScalar(0.0);
  LINALG::Export(*bubbleforces, *fextn_);

  if(radiusn_->MyLength() != 0)
  {
    // as long as it is valid, mass can be calculated just once
    double radius = (*radiusn_)[0];
    double mass = density_ * 4.0/3.0 * M_PI * radius * radius * radius;
    (*acc_)(0)->Update(1/mass, *fextn_, 0.0);
  }
  else
    (*acc_)(0)->PutScalar(0.0);

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
//   output_->WriteVector("fexternal", Fext());
   output_->WriteVector("radius", radiusn_, output_->nodevector);

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
  output_->WriteVector("radius", radiusn_, output_->nodevector);

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
  reader.ReadVector(radiusn_, "radius");

  return;
}

/*----------------------------------------------------------------------*/
/* read restart forces */
void PARTICLE::TimIntCentrDiff::ReadRestartForce()
{
  // currently do nothing
  return;
}
/*----------------------------------------------------------------------*/

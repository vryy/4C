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

  const double dt = (*dt_)[0];   // \f$\Delta t_{n}\f$
  const double dthalf = dt/2.0;  // \f$\Delta t_{n+1/2}\f$

  // new velocities \f$V_{n+1/2}\f$
  veln_->Update(1.0, *(*vel_)(0), 0.0);
  veln_->Update(dthalf, *(*acc_)(0), 1.0);

  // new displacements \f$D_{n+1}\f$
  disn_->Update(1.0, *(*dis_)(0), 0.0);
  disn_->Update(dt, *veln_, 1.0);

  // build new external forces (This is too much work here because an assembled row vector
  // for the forces would be enough. But in later applications with particle contact
  // it is necessary to evaluate also ghost bubbles.
  // in short: correct col layout is available but only row layout is needed
  // --> LINALG::Export(col->row) from bubbleforces to fextn_
  Teuchos::RCP<const Epetra_Vector> bubbleforces = discret_->GetState("particleforces");
  fextn_->PutScalar(0.0);
  LINALG::Export(*bubbleforces, *fextn_);

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

/*----------------------------------------------------------------------*/
/* State vectors are updated according to the new distribution of particles */
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

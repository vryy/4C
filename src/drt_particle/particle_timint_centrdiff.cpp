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
/* headers */
#include "particle_timint_centrdiff.H"
#include "particle_contact.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*/
/* Constructor */
PARTICLE::TimIntCentrDiff::TimIntCentrDiff(
    const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& particledynparams,
    const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<IO::DiscretizationWriter> output
  ) : PARTICLE::TimIntExpl
  (
    ioparams,
    particledynparams,
    xparams,
    actdis,
    output
  )
{
  // DetermineMassDampConsistAccel() is called at the end of Algorithm::Init() after proper setup of the problem
  return;
}


/*----------------------------------------------------------------------*/
/* mostly init of collision handling  */
void PARTICLE::TimIntCentrDiff::Init()
{
  // decide whether there is particle contact
  const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
  INPAR::PARTICLE::ContactStrategy contact_strategy =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::ContactStrategy>(particleparams,"CONTACT_STRATEGY");

  switch(contact_strategy)
  {
  case INPAR::PARTICLE::Normal_DEM:
  case INPAR::PARTICLE::NormalAndTang_DEM:
    collhandler_ = Teuchos::rcp(new PARTICLE::ParticleCollisionHandlerDEM(discret_, particle_algorithm_, particleparams));
  break;
  case INPAR::PARTICLE::Normal_MD:
    dserror("central difference time integrator cannot be combined with molecular dynamics collision mechanism");
  break;
  default:
  {
    if(myrank_ == 0)
      std::cout << "central difference time integrator is not combined with a collision model" << std::endl;
  }
  break;
  }

  // call base class init
  PARTICLE::TimIntExpl::Init();

  return;
}


/*----------------------------------------------------------------------*/
/* Integrate step */
int PARTICLE::TimIntCentrDiff::IntegrateStep()
{
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

    intergy_ = collhandler_->EvaluateParticleContact(dt, f_contact, m_contact);
  }
  //--------------------------------------------------------------

  ComputeAcc(f_contact, m_contact, accn_, ang_accn_);

  // update of end-velocities \f$V_{n+1}\f$
  veln_->Update(dthalf, *accn_, 1.0);
  if(collhandler_ != Teuchos::null)
  {
    ang_veln_->Update(dthalf,*ang_accn_,1.0);

    if(writeorientation_)
    {
      // for visualization of orientation vector
      RotateOrientVector(dt);
    }
  }

  // apply Dirichlet BCs
  ApplyDirichletBC(timen_, Teuchos::null, veln_, accn_, false);

  return 0;
}

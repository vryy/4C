/*----------------------------------------------------------------------*/
/*!
\file particle_timint_expleuler.cpp

\brief Particle time integration with Forward Euler time integration scheme of 1st order (explicit)

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_expleuler.H"
#include "particle_algorithm.H"
#include "particle_contact.H"

/*----------------------------------------------------------------------*/
/* Constructor */
PARTICLE::TimIntExplEuler::TimIntExplEuler(
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
void PARTICLE::TimIntExplEuler::Init()
{
  // safety check
  if(particle_algorithm_->ParticleInteractionType() != INPAR::PARTICLE::Normal_MD)
    dserror("explicit euler time integrator currently just combined with molecular dynamics collision mechanism");

  // call base class init
  PARTICLE::TimIntExpl::Init();

  return;
}


/*----------------------------------------------------------------------*/
/* Integrate step */
int PARTICLE::TimIntExplEuler::IntegrateStep()
{
  // new displacements \f$A_{n+1}\f$ based on external forces
  ComputeAcc(Teuchos::null, Teuchos::null, accn_, Teuchos::null);

  // \f$\Delta t_{n}\f$
  const double dt = (*dt_)[0];

  // new velocities \f$V_{n+1}\f$
  veln_->Update(1.0, *(*vel_)(0), dt, *accn_, 0.0);

  // starting point for new displacements \f$D_{n+1}\f$
  disn_->Update(1.0, *(*dis_)(0), 0.0);

  //Transfer particles into their correct bins
  particle_algorithm_->UpdateConnectivity();

  // total internal energy (elastic spring and potential energy)
  intergy_ = 0.0;
  //---------------------Compute Collisions-----------------------
  if(collhandler_ != Teuchos::null)
  {
//    // angular-velocities \f$ang_V_{n}\f$
//    ang_veln_->Update(1.0, *(*ang_vel_)(0), 0.0);

    collhandler_->Init(disn_, veln_, angVeln_, (*radius_)(0), orient_, mass_);

    intergy_ = collhandler_->EvaluateParticleContact(dt, disn_, veln_);
  }
  else
  {
    // new displacements \f$D_{n+1}\f$
    disn_->Update(1.0, *(*dis_)(0), dt, *veln_, 0.0);
  }

//  if(writeorientation_)
//  {
//    // for visualization of orientation vector
//    strategy_->RotateOrientVector(dt);
//  }

  return 0;
}

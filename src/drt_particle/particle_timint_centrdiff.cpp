/*----------------------------------------------------------------------*/
/*!
\file particle_timint_centrdiff.cpp

\brief Particle time integration with central difference scheme 2nd order (explicit),
       also known as Velocity-Verlet algorithm

\level 2

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/

/* headers */
#include "particle_timint_centrdiff.H"
#include "particle_timint_strategy.H"
#include "particle_algorithm.H"
#include "particle_contact.H"
#include "../linalg/linalg_utils.H"
#include <Teuchos_TimeMonitor.hpp>

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
  // safety checks
  if(particle_algorithm_->ParticleInteractionType() != INPAR::PARTICLE::Normal_DEM
      and particle_algorithm_->ParticleInteractionType() != INPAR::PARTICLE::NormalAndTang_DEM
      and particle_algorithm_->ParticleInteractionType() != INPAR::PARTICLE::None)
    dserror("central difference time integrator currently just combined with cavitation or discrete element collision mechanism");

  // call base class init
  PARTICLE::TimIntExpl::Init();

  return;
}

/*----------------------------------------------------------------------*/
/* Integrate step */
int PARTICLE::TimIntCentrDiff::IntegrateStep()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::TimIntCentrDiff::IntegrateStep");

  //Velocity-Verlet scheme,
  //see e.g. Adami2012, Eqs. (14)-(18)
  /*
  v_{n+1/2}=v_{n}+dt/2*a_{n}                                with              a_{n}=a(r_{n},v_{n-1/2})
  r_{n+1}=r_{n}+dt*v_{n+1/2}
  v_{n+1}=v_{n+1/2}+dt/2*a_{n+1}                            with              a_{n+1}=a(r_{n+1},v_{n+1/2})
  */

  const double dt = (*dt_)[0];   // \f$\Delta t_{n}\f$
  const double dthalf = dt/2.0;  // \f$\Delta t_{n+1/2}\f$

  // new velocities v_{n+1/2}=v_{n}+dt/2*a_{n}
  // In DEM, this procedure is known to reduce the temporal convergence order if velocity-proportional damping terms exist.
  veln_->Update(dthalf, *(*acc_)(0), 1.0);

  // new displacements r_{n+1}=r_{n}+dt*v_{n+1/2}
  disn_->Update(dt, *veln_, 1.0);

  // apply Dirichlet BCs
  ApplyDirichletBC(timen_, disn_, Teuchos::null, accn_, false);
  ApplyDirichletBC(timen_-dthalf, Teuchos::null, veln_, Teuchos::null, false);

  //Transfer particles and heat sources into their correct bins
  particle_algorithm_->UpdateConnectivity();

  // define vector for contact force and moment
  Teuchos::RCP<Epetra_Vector> f_contact = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> m_contact = Teuchos::null;

  // total internal energy (elastic spring and potential energy)
  intergy_ = 0.0;
  //---------------------Compute Collisions-----------------------
  if(collhandler_ != Teuchos::null)
  {
    // new angular-velocities \f$ang_V_{n+1/2}\f$
    angVeln_->Update(dthalf, *(*angAcc_)(0), 1.0);

    // update orientation vector
    if(writeorientation_)
      strategy_->RotateOrientVector(dt);

    // initialize vectors for contact force and moment
    f_contact = LINALG::CreateVector(*(discret_->DofRowMap()),true);
    m_contact = LINALG::CreateVector(*(discret_->DofRowMap()),true);

    collhandler_->Init(disn_, veln_, angVeln_, Radiusn(), orient_, mass_);

    intergy_ = collhandler_->EvaluateParticleContact(dt, f_contact, m_contact, f_structure_);
  }
  //--------------------------------------------------------------

  ComputeAcc(f_contact, m_contact, accn_, angAccn_);

  //--- update with the new accelerations ---//

  // update of end-velocities v_{n+1}=v_{n+1/2}+dt/2*a_{n+1}
  veln_->Update(dthalf, *accn_, 1.0);

  if(collhandler_ != Teuchos::null)
    angVeln_->Update(dthalf,*angAccn_,1.0);

  // apply Dirichlet BCs
  ApplyDirichletBC(timen_, Teuchos::null, veln_, accn_, false);

  return 0;
}

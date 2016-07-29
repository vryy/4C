/*----------------------------------------------------------------------*/
/*!
\file particle_timint_expleuler.cpp
\brief Particle time integration with Forward Euler time integration
       scheme of 1st order (explicit),

\level 3

<pre>
\maintainer Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/


/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_expleuler.H"
#include "particle_contact.H"
#include "../drt_lib/drt_globalproblem.H"


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
  // decide whether there is particle contact
  const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
  INPAR::PARTICLE::ParticleInteractions contact_strategy = DRT::INPUT::IntegralValue<INPAR::PARTICLE::ParticleInteractions>(particleparams,"PARTICLE_INTERACTION");

  switch(contact_strategy)
  {
  case INPAR::PARTICLE::Normal_DEM:
  case INPAR::PARTICLE::NormalAndTang_DEM:
    dserror("explicit euler time integrator is not yet combined with discrete element collision mechanism");
  break;
  case INPAR::PARTICLE::Normal_MD:
    collhandler_ = Teuchos::rcp(new PARTICLE::ParticleCollisionHandlerMD(discret_, particle_algorithm_, particleparams));
  break;
  default:
  {
    if(myrank_ == 0)
      std::cout << "explicit euler time integrator is not combined with a collision model" << std::endl;
  }
  break;
  }

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

  // total internal energy (elastic spring and potential energy)
  intergy_ = 0.0;
  //---------------------Compute Collisions-----------------------
  if(collhandler_ != Teuchos::null)
  {
//    // angular-velocities \f$ang_V_{n}\f$
//    ang_veln_->Update(1.0, *(*ang_vel_)(0), 0.0);

    SetStatesForCollision();

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
//    RotateOrientVector(dt);
//  }

  return 0;
}

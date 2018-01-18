/*----------------------------------------------------------------------*/
/*!
\file particle_timint_rk.cpp
\brief Particle time integration with Runge-Kutta time integration
       scheme of 2nd/4th order (explicit),

\level 3
<pre>
\maintainer Sebastian Fuchs
            fuchs@lnm.mw.tum.de
            http://www.lnm.mw.tum.de

</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_rk.H"
#include "particle_timint_strategy.H"
#include "particle_contact.H"
#include "scatra_particle_coupling.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"

/*----------------------------------------------------------------------*
 | constructor                                          rasthofer 01/14 |
 *----------------------------------------------------------------------*/
PARTICLE::TimIntRK::TimIntRK(
    const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& particledynparams,
    const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<IO::DiscretizationWriter> output
  ) : PARTICLE::TimIntExpl(ioparams,particledynparams,xparams,actdis,output),
  rk_scheme_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::DynamicType>(particledynparams,"DYNAMICTYP")),
  sign_(Teuchos::null)
{
  return;
}


/*----------------------------------------------------------------------*
 | initialization of time integration                   rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntRK::Init()
{
  // call base class routine
  PARTICLE::TimIntExpl::Init();

  if(Teuchos::rcp_dynamic_cast<PARTICLE::ScatraParticleCoupling>(particle_algorithm_) != Teuchos::null)
    sign_  =  Teuchos::rcp(new Epetra_Vector(*discret_->NodeRowMap(),true));

  return;
}


/*----------------------------------------------------------------------*
 | time integration                                     rasthofer 01/14 |
 *----------------------------------------------------------------------*/
int PARTICLE::TimIntRK::IntegrateStep()
{
  // get time integration scheme
  if (rk_scheme_ == INPAR::PARTICLE::dyna_rk2)
    // use Runge-Kutta scheme 2nd order
    Integrate_RK_Second();
  else if (rk_scheme_ == INPAR::PARTICLE::dyna_rk4)
    // use Runge-Kutta scheme 4th order
    Integrate_RK_Fourth();
  else
    dserror("Unknown time-integration order for Runge-Kutta scheme!");

  return 0;
}


/*----------------------------------------------------------------------*
 | Runge-Kutta scheme 2nd order                         rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntRK::Integrate_RK_Second()
{
  const double dt = (*dt_)[0];
  const double dthalf = dt/2.0;

  // perform time integration for scatra-particle coupling
  if(Teuchos::rcp_dynamic_cast<PARTICLE::ScatraParticleCoupling>(particle_algorithm_) != Teuchos::null)
  {
    // -------------------------------------------------------------------
    //                predictor
    // -------------------------------------------------------------------

    // set velocity of particles: using vel
    Teuchos::RCP<Epetra_Vector> vel = Teuchos::rcp_dynamic_cast<PARTICLE::ScatraParticleCoupling>(particle_algorithm_)->GetVelocity(0.0);

    // displacements at time_{n}
    Teuchos::RCP<Epetra_Vector> dis = Teuchos::rcp(new Epetra_Vector(*disn_));

    // compute position at time n+1/2
    disn_->Update(dthalf, *vel, 1.0);

    // transfer particles to their correct bin, if required
    particle_algorithm_->TransferParticles(false, false);

    // transfer state vectors to potentially new maps
    Teuchos::RCP<Epetra_Vector> old;

    old = disn_;
    disn_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
    LINALG::Export(*old, *disn_);

    old = dis;
    dis = LINALG::CreateVector(*discret_->DofRowMap(),true);
    LINALG::Export(*old, *dis);

    // -------------------------------------------------------------------
    //                corrector
    // -------------------------------------------------------------------

    // get velocity at intermediate position computed in predictor step:
    vel = Teuchos::rcp_dynamic_cast<PARTICLE::ScatraParticleCoupling>(particle_algorithm_)->GetVelocity(0.5);

    disn_->Update(1.0, *dis, dt, *vel, 0.0);

    // transfer particles to their correct bin, if required
    particle_algorithm_->TransferParticles(true, false);

    // testing output for convergence study (one particle assumed)
    //std::cout << "position " << std::setprecision(12) << (*disn_)[0] << std::endl;
  }

  // perform time integration for pure particle simulations
  else
  {
    // update particle connectivity
    particle_algorithm_->UpdateConnectivity();

    // initialize vectors for contact forces and moments
    const Teuchos::RCP<Epetra_Vector> f_contact = LINALG::CreateVector(*(discret_->DofRowMap()));
    const Teuchos::RCP<Epetra_Vector> m_contact = LINALG::CreateVector(*(discret_->DofRowMap()));

    // initialize variable for total internal energy, i.e., the sum of elastic and potential energies
    intergy_ = 0.;

    // evaluate particle-particle and particle-wall collisions at time t_n
    if(collhandler_ != Teuchos::null)
    {
      collhandler_->Init(disn_, veln_, angVeln_, Radiusn(), orient_, mass_);
      intergy_ = collhandler_->EvaluateParticleContact(dt, f_contact, m_contact);
    }

    // compute accelerations and angular accelerations at time t_n
    ComputeAcc(f_contact, m_contact, accn_, Teuchos::null);

    //
    // predictor step
    //

    // compute velocities at time t_{n+1/2}
    veln_->Update(dthalf, *accn_, 1.);

    // store orientations at time t_n
    Epetra_Vector orient(*orient_);

    // compute orientations at time t_{n+1/2}
    strategy_->RotateOrientVector(dthalf);

    // compute angular velocities at time t_{n+1/2}
    strategy_->PredictOrCorrectAngularVelocity(dthalf,*m_contact,orient);

    //
    // corrector step
    //

    // compute displacements at time t_{n+1}
    disn_->Update(dt, *veln_, 1.);

    // compute velocities at time t_{n+1}
    veln_->Update(1., *(*vel_)(0), dt, *accn_, 0.);

    // restore orientations at time t_n
    orient_->Update(1.,orient,0.);

    // compute orientations at time t_{n+1}
    strategy_->RotateOrientVector(dt);

    // restore angular velocities at time t_n
    angVeln_->Update(1.,*(*angVel_)(0),0.);

    // compute angular velocities at time t_{n+1}
    strategy_->PredictOrCorrectAngularVelocity(dt,*m_contact,orient);
  }

  return;
}


/*----------------------------------------------------------------------*
 | Runge-Kutta scheme 4th order                         rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntRK::Integrate_RK_Fourth()
{
  dserror("Runge-Kutta scheme of 4th order not yet implemented!");
  return;
}


/*----------------------------------------------------------------------*
 | update state vectors according to new particle distribution          |
 |                                                      rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntRK::UpdateStatesAfterParticleTransfer()
{
  if(Teuchos::rcp_dynamic_cast<PARTICLE::ScatraParticleCoupling>(particle_algorithm_) != Teuchos::null)
  {
    // transfer only available state vectors
    Teuchos::RCP<Epetra_Vector> old;

    if (disn_ != Teuchos::null)
    {
      old = disn_;
      disn_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
      LINALG::Export(*old, *disn_);
    }

    if (radius_ != Teuchos::null)
    {
      old = Teuchos::rcp(new Epetra_Vector(*(*radius_)(0)));
      radius_->ReplaceMaps(NodeRowMapView());

      LINALG::Export(*old, *(*radius_)(0));
    }

    if (sign_ != Teuchos::null)
    {
      old = sign_;
      sign_ = LINALG::CreateVector(*discret_->NodeRowMap(),true);
      LINALG::Export(*old, *sign_);
    }
  }

  else
    // call base class routine
    TimInt::UpdateStatesAfterParticleTransfer();

  return;
}


/*----------------------------------------------------------------------*
 | output displacement                                       fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntRK::OutputDisplacement() const
{
  output_->WriteVector("displacement", disn_);

  return;
}

/*----------------------------------------------------------------------*
 | write restart                                        rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntRK::OutputRestart
(
  bool& datawritten
)
{
  // call base class routine
  TimInt::OutputRestart(datawritten);

  if(sign_ != Teuchos::null)
    output_->WriteVector("sign", sign_, IO::nodevector);

  return;
}


/*----------------------------------------------------------------------*
 | output position, radius and sign                     rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntRK::OutputState
(
  bool& datawritten
)
{
  // call base class routine
  TimInt::OutputState(datawritten);

  // output particle sign vector if applicable
  if(sign_ != Teuchos::null)
  {
    output_->WriteVector("sign", sign_, IO::nodevector);
    output_->ClearMapCache();
  }

  return;
}


/*----------------------------------------------------------------------*
 | read and set restart state                           rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntRK::ReadRestartState()
{
  // call base class routine
  TimInt::ReadRestartState();

  // read particle sign vector if applicable
  if(sign_ != Teuchos::null)
  {
    IO::DiscretizationReader reader(discret_, step_);
    reader.ReadVector(sign_, "sign");
  }

  return;
}

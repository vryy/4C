/*----------------------------------------------------------------------*/
/*! \file
\brief Particle time integration with Runge-Kutta time integration
       scheme of 2nd/4th order (explicit),

\level 3
\maintainer Sebastian Fuchs
            fuchs@lnm.mw.tum.de
            http://www.lnm.mw.tum.de

*/

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_rk.H"
#include "particle_timint_strategy.H"
#include "particle_contact.H"
#include "particle_algorithm.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"

/*----------------------------------------------------------------------*
 | constructor                                          rasthofer 01/14 |
 *----------------------------------------------------------------------*/
PARTICLE::TimIntRK::TimIntRK(const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& particledynparams, const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<IO::DiscretizationWriter> output)
    : PARTICLE::TimIntExpl(ioparams, particledynparams, xparams, actdis, output),
      rk_scheme_(DRT::INPUT::IntegralValue<INPAR::PARTICLEOLD::DynamicType>(
          particledynparams, "DYNAMICTYP"))
{
  return;
}


/*----------------------------------------------------------------------*
 | time integration                                     rasthofer 01/14 |
 *----------------------------------------------------------------------*/
int PARTICLE::TimIntRK::IntegrateStep()
{
  // get time integration scheme
  if (rk_scheme_ == INPAR::PARTICLEOLD::dyna_rk2)
    // use Runge-Kutta scheme 2nd order
    Integrate_RK_Second();
  else if (rk_scheme_ == INPAR::PARTICLEOLD::dyna_rk4)
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
  const double dthalf = dt / 2.0;

  // perform time integration for pure particle simulations
  {
    // update particle connectivity
    particle_algorithm_->UpdateConnectivity();

    // initialize vectors for contact forces and moments
    const Teuchos::RCP<Epetra_Vector> f_contact = LINALG::CreateVector(*(discret_->DofRowMap()));
    const Teuchos::RCP<Epetra_Vector> m_contact = LINALG::CreateVector(*(discret_->DofRowMap()));

    // initialize variable for total internal energy, i.e., the sum of elastic and potential
    // energies
    intergy_ = 0.;

    // evaluate particle-particle and particle-wall collisions at time t_n
    if (collhandler_ != Teuchos::null)
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
    strategy_->PredictOrCorrectAngularVelocity(dthalf, *m_contact, orient);

    //
    // corrector step
    //

    // compute displacements at time t_{n+1}
    disn_->Update(dt, *veln_, 1.);

    // compute velocities at time t_{n+1}
    veln_->Update(1., *(*vel_)(0), dt, *accn_, 0.);

    // restore orientations at time t_n
    orient_->Update(1., orient, 0.);

    // compute orientations at time t_{n+1}
    strategy_->RotateOrientVector(dt);

    // restore angular velocities at time t_n
    angVeln_->Update(1., *(*angVel_)(0), 0.);

    // compute angular velocities at time t_{n+1}
    strategy_->PredictOrCorrectAngularVelocity(dt, *m_contact, orient);
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
 | output displacement                                       fang 10/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntRK::OutputDisplacement() const
{
  output_->WriteVector("displacement", disn_);

  return;
}

/*----------------------------------------------------------------------*/
/*!
\file particle_timint_rk.cpp
\brief Particle time integration with Runge-Kutta time integration
       scheme of 2nd/4th order (explicit),

\level 3
<pre>
\maintainer Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_rk.H"
#include "scatra_particle_coupling.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"


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
  Init();
  return;
}


/*----------------------------------------------------------------------*
 | initalization of time integration                    rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntRK::Init()
{
  // allocate vectors
  // displacements D_{n+1} at t_{n+1}
  disn_ = LINALG::CreateVector(*DofRowMapView(), true);
  // radii
  radius_  = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, NodeRowMapView(), true));
  // signs
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
  // -------------------------------------------------------------------
  //                predictor
  // -------------------------------------------------------------------

  // set velocity of particles: using vel
  Teuchos::RCP<Epetra_Vector> vel = Teuchos::rcp_dynamic_cast<PARTICLE::ScatraParticleCoupling>(particle_algorithm_)->GetVelocity(0.0);

  // compute position at time n+1/2
  const double dt = (*dt_)[0];
  const double dthalf = dt/2.0;
  // displacements at time_{n}
  Teuchos::RCP<Epetra_Vector> dis = Teuchos::rcp(new Epetra_Vector(*disn_));
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

  // transfer only available state vectors
  Teuchos::RCP<Epetra_Vector> old;

  if (disn_ != Teuchos::null)
  {
    old = disn_;
    disn_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
    LINALG::Export(*old, *disn_);
  }

  if (radius_ != Teuchos::null && (*radius_)(0) != Teuchos::null)
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

  return;
}


/*----------------------------------------------------------------------*
 | update step                                          rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntRK::UpdateStepState()
{
  // update of state vectors is not needed for RK schemes
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
  // Yes, we are going to write...
  datawritten = true;

  // mesh is written to disc
  output_->ParticleOutput(step_, (*time_)[0], true);
  output_->NewStep(step_, (*time_)[0]);
  output_->WriteVector("displacement", disn_);
  output_->WriteVector("radius", (*radius_)(0), IO::nodevector);
  output_->WriteVector("sign", sign_, IO::nodevector);

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


/*----------------------------------------------------------------------*
 | output position, radius and sign                     rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntRK::OutputState
(
  bool& datawritten
)
{
  // Yes, we are going to write...
  datawritten = true;

  // mesh is not written to disc, only maximum node id is important for output
  output_->ParticleOutput(step_, (*time_)[0], false);
  output_->NewStep(step_, (*time_)[0]);
  output_->WriteVector("displacement", disn_);

  output_->WriteVector("radius", (*radius_)(0), IO::nodevector);
  output_->WriteVector("sign", sign_, IO::nodevector);

  // maps are rebuild in every step so that reuse is not possible
  // keeps memory usage bounded
  output_->ClearMapCache();

  return;
}


/*----------------------------------------------------------------------*
 | add restart information to OutputState               rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntRK::AddRestartToOutputState()
{
  return;
}


/*----------------------------------------------------------------------*
 | read and set restart state                           rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntRK::ReadRestartState()
{
  IO::DiscretizationReader reader(discret_, step_);
  // maps need to be adapted to restarted discretization
  UpdateStatesAfterParticleTransfer();

  // now, state vectors an be read in
  reader.ReadVector(disn_, "displacement");
  reader.ReadVector((*radius_)(0), "radius");
  reader.ReadVector(sign_, "sign");

  return;
}

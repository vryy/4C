/*----------------------------------------------------------------------*/
/*!
\file particle_timint_rk.cpp
\brief Particle time integration with Runge-Kutta time integration
       scheme of 2nd/4th order (explicit),

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_rk.H"
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
  ) : PARTICLE::TimInt(ioparams,particledynparams,xparams,actdis,output),
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
  radius_  = Teuchos::rcp(new Epetra_Vector(*discret_->NodeRowMap(),true));
  sign_  =  Teuchos::rcp(new Epetra_Vector(*discret_->NodeRowMap(),true));
  
  return;
}


/*----------------------------------------------------------------------*
 | time integration                                     rasthofer 01/14 |
 *----------------------------------------------------------------------*/
int PARTICLE::TimIntRK::IntegrateStep()
{
  // get time integration scheme
  if (rk_scheme_ == INPAR::PARTICLE::dyna_rungekutta2)
    // use Runge-Kutta scheme 2nd order
    Integrate_RK_Second();
  else if (rk_scheme_ == INPAR::PARTICLE::dyna_rungekutta4)
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
  dserror("Runge-Kutta scheme of 2nd order not yet implemented!");

# if 0
  // TODO: get right velocity
  // -------------------------------------------------------------------
  //                predictor
  // -------------------------------------------------------------------
  // set velocity of particles: using veln
  SetVelocity();
  // get pointers to state vectors
  Teuchos::RCP<Epetra_Vector> velnp = particles_->ExtractVelnp();
  Teuchos::RCP<Epetra_Vector> disnp = particles_->ExtractDispn();

  // compute position at time n+1/2
  const double dthalf = dta_/2.0;
  disnp->Update(1.0, *particles_->ExtractDispn(), dthalf, *velnp, 0.0);

  // transfer particles to their correct bin, if required
  TransferParticles(false);
  // transfer state vectors to potentially new maps
  Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->UpdateStatesAfterParticleTransfer();

  // UpdateStatesAfterParticleTransfer() changes memory location, we thus have to adapt the pointers
  velnp = particles_->ExtractVelnp();
  disnp = particles_->ExtractDispnp();

  // -------------------------------------------------------------------
  //                corrector
  // -------------------------------------------------------------------
  // get velocity at intermediate position computed in predictor step: using veln to approximate velnp0.5
  SetVelocity();
  disnp->Update(1.0, *particles_->ExtractDispn(), dta_, *velnp, 0.0);

  // transfer particles to their correct bin, if required
  TransferParticles(false);
  // transfer state vectors to potentially new maps
  Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->UpdateStatesAfterParticleTransfer();

  // disnp = particles_->ExtractDispnp();
  //std::cout << "positions are " << std::setprecision(12) << *disnp << endl;
#endif

  return;
}


/*----------------------------------------------------------------------*
 | Runge-Kutta scheme 4th order                         rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntRK::Integrate_RK_Fourth()
{
  dserror("Runge-Kutta scheme of 4th order not yet implemented!");
#if 0 // ORIGINAL DATONG: VORSICHT!!!!!!!!!!!!!!
  // -------------------------------------------------------------------
  //                predictor
  // -------------------------------------------------------------------
  // set velocity of particles: using veln
  SetVelocity();
  // get pointers to state vectors
  Teuchos::RCP<Epetra_Vector> veln = particles_->ExtractVelnp();
  Teuchos::RCP<Epetra_Vector> disnp = particles_->ExtractDispnp();

  // compute position at time n+1/2
  const double dthalf = dta_/2.0;

  Teuchos::RCP<Epetra_Vector> disnp_halfpredictor = Teuchos::rcp(new Epetra_Vector(*(particledis_->DofRowMap()),true));
  Teuchos::RCP<Epetra_Vector> veln   = Teuchos::rcp(new Epetra_Vector(*(particledis_->DofRowMap()),true));

  Teuchos::RCP<Epetra_Vector> velnp_temp_one   = Teuchos::rcp(new Epetra_Vector(*(particledis_->DofRowMap()),true));
  Teuchos::RCP<Epetra_Vector> velnp_temp_two   = Teuchos::rcp(new Epetra_Vector(*(particledis_->DofRowMap()),true));
  Teuchos::RCP<Epetra_Vector> velnp_temp_three = Teuchos::rcp(new Epetra_Vector(*(particledis_->DofRowMap()),true));

  Teuchos::RCP<Epetra_Vector> disnp_temp_zero        = Teuchos::rcp(new Epetra_Vector(*(particledis_->DofRowMap()),true));
  Teuchos::RCP<Epetra_Vector> disnp_temp_part_one    = Teuchos::rcp(new Epetra_Vector(*(particledis_->DofRowMap()),true));
  Teuchos::RCP<Epetra_Vector> disnp_temp_part_two    = Teuchos::rcp(new Epetra_Vector(*(particledis_->DofRowMap()),true));
  Teuchos::RCP<Epetra_Vector> disnp_temp_one_and_two = Teuchos::rcp(new Epetra_Vector(*(particledis_->DofRowMap()),true));




  //1st predictor
  disnp_temp_zero->Update(1.0,*(particles_->Dispnp()),0.0);  //saves positions @ n
  velnp_temp_one->Update(1.0,*(particles_->ExtractVelnp()),0.0);  //save velocities @ n
  disnp->Update(dthalf, *velnp, 1.0); //positions @ n+1/2*
  TransferParticles(false);
  Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->UpdateStatesAfterParticleTransfer();
  // UpdateStatesAfterParticleTransfer() changes memory location, we thus have to adapt the pointers
  velnp = particles_->ExtractVelnp();
  disnp = particles_->ExtractDispnp();


  //2nd predictor
  SetVelocity(); //update velocities to state @ n+0.5*
  velnp_temp_two->Update(1.0,*(particles_->ExtractVelnp()),0.0); //save velocities @ n+0.5*
  disnp->Update(dthalf, *velnp, 1.0, *disnp_temp_zero, 0.0); //positions @ n+1/2 **
  TransferParticles(false);
  Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->UpdateStatesAfterParticleTransfer();
  // UpdateStatesAfterParticleTransfer() changes memory location, we thus have to adapt the pointers
  velnp = particles_->ExtractVelnp();
  disnp = particles_->ExtractDispnp();


  //3rd predictor
  SetVelocity(); //velocities @ n+0.5**
  velnp_temp_three->Update(1.0,*(particles_->ExtractVelnp()),0.0);  //save velocities @ n+0.5**
  disnp->Update(1.0, *disnp_temp_zero, dt, *velnp, 0.0); //positions @ n+1*
  TransferParticles(false);
  Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->UpdateStatesAfterParticleTransfer();
  // UpdateStatesAfterParticleTransfer() changes memory location, we thus have to adapt the pointers
  velnp = particles_->ExtractVelnp();
  disnp = particles_->ExtractDispnp();


  //4th step: correction
  SetVelocity();  //velocities at n+1*
  disnp_temp_part_one->Update(dt/6.0, *velnp_temp_one, dt/3.0, *velnp_temp_two, 0.0);
  disnp_temp_part_two->Update(dt/3.0, *velnp_temp_three, dt/6.0, *velnp, 0.0);
  disnp_temp_one_and_two->Update(1.0, *disnp_temp_part_one, 1.0, *disnp_temp_part_two, 0.0);
  disnp->Update(1.0, *disnp_temp_one_and_two, 1.0, *disnp_temp_zero, 0.0);
  TransferParticles(false);
  Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->UpdateStatesAfterParticleTransfer();
  // UpdateStatesAfterParticleTransfer() changes memory location, we thus have to adapt the pointers
  velnp = particles_->ExtractVelnp();
  disnp = particles_->ExtractDispnp();

  std::cout << "positions are " << std::setprecision(12) << *disnp << endl;

#endif
  return;
}


/*----------------------------------------------------------------------*
 | update state vectors according to new particle distribution          |
 |                                                      rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntRK::UpdateStatesAfterParticleTransfer()
{
  // first, transfer state vectors of basic class
  PARTICLE::TimInt::UpdateStatesAfterParticleTransfer();

  // then transfer additional state vectors for particle level-set method
  Teuchos::RCP<Epetra_Vector> old;
  if (sign_ != Teuchos::null)
  {
    old = sign_;
    radius_ = LINALG::CreateVector(*discret_->NodeRowMap(),true);
    LINALG::Export(*old, *sign_);
  }

  return;
}


/*----------------------------------------------------------------------*
 | update step                                          rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::TimIntRK::UpdateStepState()
{
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1}, D_{n-1} := D_{n}
  dis_->UpdateSteps(*disn_);
  // new velocities at t_{n+1} -> t_n
  //    V_{n} := V_{n+1}, V_{n-1} := V_{n}
  vel_->UpdateSteps(*veln_);

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
  output_->WriteVector("displacement", (*dis_)(0));
  output_->WriteVector("radius", radius_, output_->nodevector);
  output_->WriteVector("sign", sign_, output_->nodevector);

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
  output_->WriteVector("displacement", (*dis_)(0));

  output_->WriteVector("radius", radius_, output_->nodevector);
  output_->WriteVector("sign", sign_, output_->nodevector);

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
  dis_->UpdateSteps(*disn_);

  reader.ReadVector(radius_, "radius");
  reader.ReadVector(sign_, "sign");

  return;
}

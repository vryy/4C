/*----------------------------------------------------------------------*/
/*!
\file particle_timint_expl.cpp
\brief Particle time integration with explicit time integration

\level 2

<pre>
\maintainer Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/


/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_expl.H"
#include "particle_algorithm.H"
#include "particle_contact.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_timestepping/timintmstep.H"


/*----------------------------------------------------------------------*/
/* Constructor */
PARTICLE::TimIntExpl::TimIntExpl(
    const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& particledynparams,
    const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<IO::DiscretizationWriter> output
  ) : PARTICLE::TimInt
  (
    ioparams,
    particledynparams,
    xparams,
    actdis,
    output
  )
{
  return;
}


/*----------------------------------------------------------------------*/
/* mostly init of collision handling  */
void PARTICLE::TimIntExpl::Init()
{
  // call base class init
  TimInt::Init();

  // decide whether there is particle contact
  const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
  INPAR::PARTICLE::ParticleInteractions contact_strategy = DRT::INPUT::IntegralValue<INPAR::PARTICLE::ParticleInteractions>(particleparams,"PARTICLE_INTERACTION");

  if(contact_strategy != INPAR::PARTICLE::None)
  {
    // allocate vectors
    inertia_  = LINALG::CreateVector(*discret_->NodeRowMap(), true);

    ang_vel_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
    ang_acc_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));

    ang_veln_ = LINALG::CreateVector(*DofRowMapView(),true);
    ang_accn_ = LINALG::CreateVector(*DofRowMapView(),true);

    if(writeorientation_)
    {
      // initialize orientation-vector for visualization
      orient_ = LINALG::CreateVector(*DofRowMapView(),true);
      InitializeOrientVector();
    }

    // initialize inertia
    for(int n=0; n<discret_->NumMyRowNodes(); ++n)
    {
      double r_p = (*radius_)[n];

      //inertia-vector: sphere: I = 2/5 * m * r^2
      (*inertia_)[n] = 0.4 * (*mass_)[n] * r_p * r_p;
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/* update step */
void PARTICLE::TimIntExpl::UpdateStepState()
{
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1}, D_{n-1} := D_{n}
  dis_->UpdateSteps(*disn_);
  // new velocities at t_{n+1} -> t_n
  //    V_{n} := V_{n+1}, V_{n-1} := V_{n}
  vel_->UpdateSteps(*veln_);
  // new accelerations at t_{n+1} -> t_n
  //    A_{n} := A_{n+1}, A_{n-1} := A_{n}
  acc_->UpdateSteps(*accn_);
  //    T_{n} := T_{n+1}, T_{n-1} := T_{n}
  if (particle_algorithm_->trg_Temperature()) temperature_->UpdateSteps(*temperaturen_);

  if(collhandler_ != Teuchos::null)
  {
    // new angular-velocities at t_{n+1} -> t_n
    //    ang_V_{n} := ang_V_{n+1}, ang_V_{n-1} := ang_V_{n}
    ang_vel_->UpdateSteps(*ang_veln_);
    // new angular-accelerations at t_{n+1} -> t_n
    //    ang_A_{n} := ang_A_{n+1}, ang_A_{n-1} := ang_A_{n}
    ang_acc_->UpdateSteps(*ang_accn_);
  }

  return;
}


/*----------------------------------------------------------------------*/
/* states are given to the collision handler */
void PARTICLE::TimIntExpl::SetStatesForCollision()
{
  // dof based vectors
  discret_->SetState("bubblepos", disn_);
  discret_->SetState("bubblevel", veln_);
  discret_->SetState("bubbleangvel", ang_veln_);

  collhandler_->SetState(radius_, mass_);

  return;
}


/*----------------------------------------------------------------------*/
/* initialization of vector for visualization of the particle orientation */
void PARTICLE::TimIntExpl::InitializeOrientVector()
{
  int numrownodes = discret_->NodeRowMap()->NumMyElements();
  for(int i=0; i<numrownodes; ++i)
  {
    (*orient_)[i*3] = 0.0;
    (*orient_)[i*3+1] = 0.0;
    (*orient_)[i*3+2] = 1.0;
  }

  return;
}


/*----------------------------------------------------------------------*/
/* update of vector for visualization of the particle orientation */
void PARTICLE::TimIntExpl::RotateOrientVector(double dt)
{
  int numrownodes = discret_->NodeRowMap()->NumMyElements();
  for(int i=0; i<numrownodes; ++i)
  {
    double ang_vel[3];
    double r[3];

    for(int dim=0; dim<3; ++dim)
    {
      ang_vel[dim] = (*ang_veln_)[i*3+dim];
      r[dim] = (*orient_)[i*3+dim];
    }

    // visualization only valid for 2D when rotating around y-axis

    // simplified/linearized orient vector - just for visualization
    // delta_r = \Delta t * (ang_vel x r)
    (*orient_)[i*3]   += dt * (ang_vel[1] * r[2] - ang_vel[2] * r[1]);
    (*orient_)[i*3+1] += dt * (ang_vel[2] * r[0] - ang_vel[0] * r[2]);
    (*orient_)[i*3+2] += dt * (ang_vel[0] * r[1] - ang_vel[1] * r[0]);
    //--------------------------------------------------------------

    //more exactly------------------------------------------------
//    double d[3];
//    double norm_ang_vel=0.0;
//    //norm ang_vel
//    for(int dim=0; dim<3; ++dim)
//    {
//      norm_ang_vel += ang_vel[dim]*ang_vel[dim];
//    }
//    norm_ang_vel = sqrt(norm_ang_vel);
//
//    if(norm_ang_vel > 1E-10)
//    {
//      double scalar = 0.0;
//      double invnorm_ang_vel = 1.0 / norm_ang_vel;
//      for(int dim=0; dim<3; ++dim)
//      {
//        d[dim] = invnorm_ang_vel * ang_vel[dim];
//        scalar += d[dim] * r[dim];
//      }
//
//      for(int dim=0; dim<3; ++dim)
//      {
//        (*orient_)[i*3+dim] += (1-cos(norm_ang_vel*dt)) * scalar * d[dim] - (1-cos(norm_ang_vel*dt)) * r[dim]
//              + sin(norm_ang_vel*dt) * ( d[(dim+1)%3] * r[(dim+2)%3] - d[(dim+2)%3] * r[(dim+1)%3] );
//      }
//    }
    //--------------------------------------------------------------

  }

  return;
}

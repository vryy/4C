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

  switch (particle_algorithm_->ParticleInteractionType())
  {
  case INPAR::PARTICLE::MeshFree :
  case INPAR::PARTICLE::Normal_DEM_thermo :
  {
    //    R_{n} := R_{n+1}, R_{n-1} := R_{n}
    radius_->UpdateSteps(*radiusn_);
    //    D_{n} := D_{n+1}, D_{n-1} := D_{n}
    density_->UpdateSteps(*densityn_);
    //    H_{n} := H_{n+1}, H_{n-1} := H_{n}
    specEnthalpy_->UpdateSteps(*specEnthalpyn_);
    break;
  }
  default : // do nothing
    break;
  }

  if(collhandler_ != Teuchos::null)
  {
    // new angular-velocities at t_{n+1} -> t_n
    //    ang_V_{n} := ang_V_{n+1}, ang_V_{n-1} := ang_V_{n}
    angVel_->UpdateSteps(*angVeln_);
    // new angular-accelerations at t_{n+1} -> t_n
    //    ang_A_{n} := ang_A_{n+1}, ang_A_{n-1} := ang_A_{n}
    angAcc_->UpdateSteps(*angAccn_);
  }

  return;
}


/*----------------------------------------------------------------------*/
/* states are given to the collision handler */
void PARTICLE::TimIntExpl::SetStatesForCollision()
{
  switch (particle_algorithm_->ParticleInteractionType())
  {
  case INPAR::PARTICLE::MeshFree :
  case INPAR::PARTICLE::Normal_DEM_thermo :
  {
    collhandler_->Init(disn_, veln_, angVeln_, radiusn_, mass_);
    break;
  }
  default :
  {
    collhandler_->Init(disn_, veln_, angVeln_, (*radius_)(0), mass_);
    break;
  }
  }
}

/*----------------------------------------------------------------------*/
/* update of vector for visualization of the particle orientation */
void PARTICLE::TimIntExpl::RotateOrientVector(double dt)
{
  int numrownodes = discret_->NodeRowMap()->NumMyElements();
  for(int i=0; i<numrownodes; ++i)
  {
    double angVel[3];
    double r[3];

    for(int dim=0; dim<3; ++dim)
    {
      angVel[dim] = (*angVeln_)[i*3+dim];
      r[dim] = (*orient_)[i*3+dim];
    }

    // visualization only valid for 2D when rotating around y-axis

    // simplified/linearized orient vector - just for visualization
    // delta_r = \Delta t * (ang_vel x r)
    (*orient_)[i*3]   += dt * (angVel[1] * r[2] - angVel[2] * r[1]);
    (*orient_)[i*3+1] += dt * (angVel[2] * r[0] - angVel[0] * r[2]);
    (*orient_)[i*3+2] += dt * (angVel[0] * r[1] - angVel[1] * r[0]);
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

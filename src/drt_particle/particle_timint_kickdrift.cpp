/*----------------------------------------------------------------------*/
/*!
\file particle_timint_kickdrift.cpp

\brief SPH 2nd order (explicit) time integration scheme denoted as kick-drift-kick scheme as presented in Adami 2013, Eqs (14)-(17)

\level 2

<pre>
\maintainer Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/


/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_kickdrift.H"
#include "particle_algorithm.H"
#include "particle_contact.H"
#include "particle_utils.H"
#include "particleMeshFree_interaction.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/extparticle_mat.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "particle_algorithm.H"
#include "particle_heatSource.H"

#include "../drt_io/io_control.H"

/*----------------------------------------------------------------------*/
/* Constructor */
PARTICLE::TimIntKickDrift::TimIntKickDrift(
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
  return;
}


/*----------------------------------------------------------------------*/
/* mostly init of collision handling  */
void PARTICLE::TimIntKickDrift::Init()
{

  if (particle_algorithm_->ExtParticleMat()->surfaceVoidTension_ != 0)
    dserror("No surface tension possible in kick-drift scheme!");

  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"DENSITY_SUMMATION")==false)
    dserror("No density integration possible in kick-drift scheme!");

  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"SOLVE_THERMAL_PROBLEM"))
    dserror("No thermal problem possible in kick-drift scheme!");

  const INPAR::PARTICLE::WallInteractionType wallInteractionType=DRT::INPUT::IntegralValue<INPAR::PARTICLE::WallInteractionType>(DRT::Problem::Instance()->ParticleParams(),"WALL_INTERACTION_TYPE");
  if(wallInteractionType==INPAR::PARTICLE::Mirror or wallInteractionType==INPAR::PARTICLE::Custom or wallInteractionType==INPAR::PARTICLE::InitParticle)
    dserror("the wall interaction types Mirror, Custom and InitParticle not possible in kick-drift scheme!");

  // call base class init
  PARTICLE::TimIntExpl::Init();

  // decide whether there is particle contact
  const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();

  if(particle_algorithm_->ParticleInteractionType()==INPAR::PARTICLE::MeshFree)
  {
    interHandler_ = Teuchos::rcp(new PARTICLE::ParticleMeshFreeInteractionHandler(discret_, particle_algorithm_, particleparams, initDensity_, restDensity_, refdensfac_));
  }
  else
    dserror("The kick-drift-kick time integrator can currently exclusively be used for SPH/Meshfree applications!");
}

/*----------------------------------------------------------------------*/
/* Integrate step */
int PARTICLE::TimIntKickDrift::IntegrateStep()
{
  //Kick-drift-kick integration scheme with modified particle convection \tilde{v} due to modified acceleration contribution \tilde{a},
  //see e.g. Adami2013, Eqs. (14)-(17)
  /*
  v_{n+1/2}=v_{n}+dt/2*a_{n-1/2}                            with              a_{n-1/2}=a(r_{n},v_{n-1/2})
  \tilde{v}_{n+1/2}=v_{n+1/2}+dt/2*\tilde{a}_{n-1/2}        with      \tilde{a}_{n-1/2}=\tilde{a}(r_{n})
  r_{n+1}=r_{n}+dt*\tilde{v}_{n+1/2}
  v_{n+1}=v_{n+1/2}+dt/2*a_{n+1/2}                          with              a_{n+1/2}=a(r_{n+1},v_{n+1/2})
  */

  //Attention: At the end of IntegrateStep(), positions disn_=r_{n+1} and velocities veln_=v_{n+1} at the end time t_{n+1} are stored,
  //while accelerations accn_=a_{n+1/2} are stored at the mid time t_{n+1/2}!!!

  const double dt = (*dt_)[0];   // \f$\Delta t_{n}\f$
  const double dthalf = dt/2.0;  // \f$\Delta t_{n+1/2}\f$

  //v_{n+1/2}=v_{n}+dt/2*a_{n-1/2}, with *(*acc_)(0)=a_{n-1/2}
  veln_->Update(dthalf, *(*acc_)(0), 1.0);
  //\tilde{v}_{n+1/2}=v_{n+1/2}+dt/2*\tilde{a}_{n-1/2}, with \tilde{a}_{n-1/2}=*(*accmod_)(0), \tilde{v}_{n+1/2}=velConv
  Teuchos::RCP<Epetra_Vector> velConv = Teuchos::rcp(new Epetra_Vector(*veln_));
  velConv->Update(1.0,*veln_,0.0);

  //Apply modified convection velocity if required
  if (DRT::Problem::Instance()->ParticleParams().get<double>("BACKGROUND_PRESSURE")>0.0)
    velConv->Update(dthalf,*(*accmod_)(0),1.0);

  //r_{n+1}=r_{n}+dt*\tilde{v}_{n+1/2}=r_{n}+dt*\tilde{v}_{n+1/2}, with \tilde{v}_{n+1/2}=velConv
  disn_->Update(dt, *velConv, 1.0);

  //Apply Dirichlet BCs
  //(Accelerations accn_ at Dirichlet DoFs are required for boundary particles due to pressure averaging as suggested in Adami2012)
  ApplyDirichletBC(timen_, disn_, Teuchos::null, accn_, false);
  ApplyDirichletBC(timen_-dthalf, Teuchos::null, veln_, Teuchos::null, false);

  //Transfer particles and heat sources into their correct bins
  particle_algorithm_->UpdateConnectivity();

  //Determine density as well as physical and modified accelerations a_{n+1/2}=a(r_{n+1},v_{n+1/2})
  //and \tilde{a}_{n+1/2}=\tilde{a}(r_{n+1},v_{n+1/2}) based on r_{n+1} and v_{n+1/2} as well as external forces at t_{n+1}
  DetermineMeshfreeDensAndAcc(accn_,accmodn_,velConv,timen_);

  //Update of end-velocities v_{n+1}=v_{n+1/2}+dt/2*a_{n+1/2}
  veln_->Update(dthalf, *accn_, 1.0);

  //Apply Dirichlet BCs at t_{n+1} for v_{n+1}=veln_
  ApplyDirichletBC(timen_, Teuchos::null, veln_, Teuchos::null, false);
  //Apply Dirichlet BCs at t_{n+1/2} for a_{n+1/2}=accn_
  ApplyDirichletBC(timen_-dthalf, Teuchos::null, Teuchos::null, accn_, false);

  return 0;
}

/*----------------------------------------------------------------------*/
/* Read and set restart state */
void PARTICLE::TimIntKickDrift::ReadRestartState()
{
  // call the base class
  TimInt::ReadRestartState();

  IO::DiscretizationReader reader(discret_, step_);

  // read densityDot
  reader.ReadVector(densityDotn_, "densityDot");
  densityDot_->UpdateSteps(*densityDotn_);

  bool solve_thermal_problem=DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"SOLVE_THERMAL_PROBLEM");
  PARTICLE::Utils::Density2Pressure(restDensity_,refdensfac_,densityn_,specEnthalpyn_,pressure_,particle_algorithm_->ExtParticleMat(),true,solve_thermal_problem);
}

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
  ) :
  PARTICLE::TimIntExpl
  (
    ioparams,
    particledynparams,
    xparams,
    actdis,
    output
  ),
  alltimemin_pvp_dist_(1.0e10),
  alltimemin_pvw_dist_(1.0e10),
  alltimemin_pressure_(1.0e10),
  alltimemax_pressure_(-1.0e10)
{
  return;
}


/*----------------------------------------------------------------------*/
/* mostly init of collision handling  */
void PARTICLE::TimIntKickDrift::Init()
{
  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"SOLVE_THERMAL_PROBLEM"))
    dserror("No thermal problem possible in kick-drift scheme!");

  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"NO_VELDIFF_TERM")==true and DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"TRANSPORT_VELOCITY")==false)
      dserror("The parameter NO_VELDIFF_TERM only makes sense in combination with TRANSPORT_VELOCITY!");

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
  //Apply modified convection velocity if required
  if (DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"TRANSPORT_VELOCITY")==true)
  {
    //\tilde{v}_{n+1/2}=v_{n+1/2}+dt/2*\tilde{a}_{n-1/2}, with \tilde{a}_{n-1/2}=*(*accmod_)(0), \tilde{v}_{n+1/2}=velmodn_
    velmodn_->Update(1.0, *veln_, 0.0);
    velmodn_->Update(dthalf,*(*accmod_)(0),1.0);

    //r_{n+1}=r_{n}+dt*\tilde{v}_{n+1/2}, with \tilde{v}_{n+1/2}=velmodn_
    disn_->Update(dt, *velmodn_, 1.0);
  }
  else
  {
    //r_{n+1}=r_{n}+dt*v_{n+1/2}, with v_{n+1/2}=veln_
    disn_->Update(dt, *veln_, 1.0);
  }

  //Apply Dirichlet BCs
  //(Accelerations accn_ at Dirichlet DoFs are required for boundary particles due to pressure averaging as suggested in Adami2012)
  ApplyDirichletBC(timen_, disn_, Teuchos::null, accn_, false);
  ApplyDirichletBC(timen_-dthalf, Teuchos::null, veln_, Teuchos::null, false);

  if (DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"TRANSPORT_VELOCITY")==true)
    ApplyDirichletBC(timen_-dthalf, Teuchos::null, velmodn_, Teuchos::null, false);

  //Transfer particles and heat sources into their correct bins
  particle_algorithm_->UpdateConnectivity();

  Teuchos::RCP<Epetra_Vector> acc_A = Teuchos::rcp(new Epetra_Vector(*veln_));
  acc_A->PutScalar(0.0);

  //Determine density as well as physical and modified accelerations a_{n+1/2}=a(r_{n+1},v_{n+1/2})
  //and \tilde{a}_{n+1/2}=\tilde{a}(r_{n+1},v_{n+1/2}) based on r_{n+1} and v_{n+1/2} as well as external forces at t_{n+1}
  DetermineMeshfreeDensAndAcc(accn_,accmodn_,velmodn_,acc_A,timen_,dt);

  //Update of end-velocities v_{n+1}=v_{n+1/2}+dt/2*a_{n+1/2}
  veln_->Update(dthalf, *accn_, 1.0);

  //Apply Dirichlet BCs at t_{n+1} for veln_ and accn_
  ApplyDirichletBC(timen_, Teuchos::null, veln_, accn_, false);

  //Determination (and output) of extreme values
  GetExtremeValues();

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

  // read radius
  reader.ReadVector(radiusn_, "radius");
  radius_->UpdateSteps(*radiusn_);
  // read density
  reader.ReadVector(densityn_, "density");
  density_->UpdateSteps(*densityn_);
  // read specEnthalpy
  reader.ReadVector(specEnthalpyn_, "specEnthalpy");
  specEnthalpy_->UpdateSteps(*specEnthalpyn_);
  // read specEnthalpyDot
  reader.ReadVector(specEnthalpyDotn_, "specEnthalpyDot");
  specEnthalpyDot_->UpdateSteps(*specEnthalpyDotn_);
}

/* Determination (and output) of extreme values */
void PARTICLE::TimIntKickDrift::GetExtremeValues()
{
  double min_pvp_dist=0.0;
  double min_pvw_dist=0.0;
  double min_pressure=0.0;
  double max_pressure=0.0;
  /// Get min / max pressure etc.
  interHandler_->GetExtremeValues(min_pvp_dist,min_pvw_dist,min_pressure,max_pressure);

  if(min_pvp_dist<alltimemin_pvp_dist_)
    alltimemin_pvp_dist_=min_pvp_dist;
  if(min_pvw_dist<alltimemin_pvw_dist_)
    alltimemin_pvw_dist_=min_pvw_dist;
  if(min_pressure<alltimemin_pressure_)
    alltimemin_pressure_=min_pressure;
  if(max_pressure>alltimemax_pressure_)
    alltimemax_pressure_=max_pressure;

  double allprocalltime_min_pvp_dist=0.0;
  double allprocalltime_min_pvw_dist=0.0;
  double allprocalltime_min_pressure=0.0;
  double allprocalltime_max_pressure=0.0;
  double maxradius=0.0;

  discret_->Comm().MinAll(&alltimemin_pvp_dist_,&allprocalltime_min_pvp_dist,1);
  discret_->Comm().MinAll(&alltimemin_pvw_dist_,&allprocalltime_min_pvw_dist,1);
  discret_->Comm().MinAll(&alltimemin_pressure_,&allprocalltime_min_pressure,1);
  discret_->Comm().MaxAll(&alltimemax_pressure_,&allprocalltime_max_pressure,1);
  radiusn_->NormInf(&maxradius);

  #ifdef PARTICLE_MINMAXOUTPUT_EVRY
  if((discret_->Comm().MyPID()==0) and (stepn_ % 100==0))
  {
    const double c = particle_algorithm_->ExtParticleMat()->SpeedOfSoundL();
    std::cout << "allprocalltime_min_pvp_dist: " << allprocalltime_min_pvp_dist << std::endl;
    std::cout << "allprocalltime_min_pvw_dist: " << allprocalltime_min_pvw_dist << std::endl;
    std::cout << "maximal particle radius: " << maxradius << std::endl;
    std::cout << "allprocalltime_min_pressure: " << allprocalltime_min_pressure << std::endl;
    std::cout << "allprocalltime_max_pressure: " << allprocalltime_max_pressure << std::endl;
    std::cout << "reference pressure: " <<  c * c * restDensity_ << std::endl;
    std::cout << "background pressure: " << DRT::Problem::Instance()->ParticleParams().get<double>("BACKGROUND_PRESSURE") << std::endl;
  }
  #endif
}

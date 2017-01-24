/*----------------------------------------------------------------------*/
/*!
\file particle_timint_genalpha.cpp
\brief Implicit particle time integration scheme (backward Euler?)

\level 2

<pre>
\maintainer Alessandro Cattabiani
</pre>
*/


/*----------------------------------------------------------------------*/
/* headers */
#include "particle_timint_genalpha.H"

#include "particle_algorithm.H"
#include "particleMeshFree_interaction.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/extparticle_mat.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "particle_algorithm.H"
#include "particle_heatSource.H"
/*----------------------------------------------------------------------*/
/* Constructor */
PARTICLE::TimIntGenAlpha::TimIntGenAlpha(
    const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& particledynparams,
    const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<IO::DiscretizationWriter> output
  ) : PARTICLE::TimIntImpl(ioparams, particledynparams, xparams, actdis, output),
  midavg_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::MidAverageEnum>(particledynparams.sublist("GENALPHA"),"GENAVG")),
  beta_(particledynparams.sublist("GENALPHA").get<double>("BETA")),
  gamma_(particledynparams.sublist("GENALPHA").get<double>("GAMMA")),
  alphaf_(particledynparams.sublist("GENALPHA").get<double>("ALPHA_F")),
  alpham_(particledynparams.sublist("GENALPHA").get<double>("ALPHA_M")),
  rho_inf_(particledynparams.sublist("GENALPHA").get<double>("RHO_INF"))
{

  return;
}


/*----------------------------------------------------------------------*/
/* mostly init of collision handling  */
void PARTICLE::TimIntGenAlpha::Init()
{
  // check
  if (particle_algorithm_->ParticleInteractionType() != INPAR::PARTICLE::MeshFree)
  {
    dserror("this interaction model is not combined with the hybrid divergence free scheme");
  }

  // call base class init
  PARTICLE::TimIntImpl::Init();

  // initialize the additional state vectors
  mGradW_ = LINALG::CreateVector(*DofRowMapView(), true);
  mHessW_ = Teuchos::rcp(new Epetra_MultiVector(*DofRowMapView(), 3, true));
  accm_ = LINALG::CreateVector(*DofRowMapView(), true);
  velm_ = LINALG::CreateVector(*DofRowMapView(), true);
  dism_ = LINALG::CreateVector(*DofRowMapView(), true);

  // initialize the interaction handler
  const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
  interHandler_ = Teuchos::rcp(new PARTICLE::ParticleMeshFreeInteractionHandler(discret_, particle_algorithm_, particleparams, restDensity_));

  // set up the genalpha parameters...
  CalcCoeff();
  // ... and check them
  VerifyCoeff();
  // check that there are no convection or diffusion forces. genalpha is still not able to handle them
  if (particle_algorithm_->ExtParticleMat()->bulkViscosity_ != 0 || particle_algorithm_->ExtParticleMat()->dynamicViscosity_ != 0)
  {
    dserror("Genalpha is still unable to address the diffusion and convection terms");
  }

}

/*----------------------------------------------------------------------*/
/* Integrate step */
int PARTICLE::TimIntGenAlpha::IntegrateStep()
{


  interHandler_->Clear(step_,1);
  dserror("IntegrateStep is still in the todo list");

  return 0;
}

/*----------------------------------------------------------------------*/
// overload to determine the initial accelerations
void PARTICLE::TimIntGenAlpha::DetermineMassDampConsistAccel()
{
  // create the wall discretization
  particle_algorithm_->SetUpWallDiscret();

  // set up connections and the local state vectors in the interaction handler
  interHandler_->Init(disn_, veln_, radiusn_, mass_, densityn_, specEnthalpyn_, pressure_, temperature_,  stepn_);

  // compute and replace the density
  interHandler_->MF_mW(densityn_);
  interHandler_->SetStateVector(densityn_, PARTICLE::Density);

  interHandler_->Inter_pvp_acc(accn_);
  interHandler_->Inter_pvw_acc(accn_);
  acc_->UpdateSteps(*accn_);

  return;
}

/*----------------------------------------------------------------------*/
/* evaluate mid-state vectors by averaging end-point vectors */
void PARTICLE::TimIntGenAlpha::EvaluateMidState()
{
  // mid-displacements D_{n+1-alpha_f} (dism)
  //    D_{n+1-alpha_f} := (1.-alphaf) * D_{n+1} + alpha_f * D_{n}
  dism_->Update(1.-alphaf_, *disn_, alphaf_, (*dis_)[0], 0.0);

  // mid-velocities V_{n+1-alpha_f} (velm)
  //    V_{n+1-alpha_f} := (1.-alphaf) * V_{n+1} + alpha_f * V_{n}
  velm_->Update(1.-alphaf_, *veln_, alphaf_, (*vel_)[0], 0.0);

  // mid-accelerations A_{n+1-alpha_m} (accm)
  //    A_{n+1-alpha_m} := (1.-alpha_m) * A_{n+1} + alpha_m * A_{n}
  accm_->Update(1.-alpham_, *accn_, alpham_, (*acc_)[0], 0.0);

  // jump
  return;
}

/*----------------------------------------------------------------------*/
void PARTICLE::TimIntGenAlpha::VerifyCoeff()
{
  // beta
  if ( (beta_ <= 0.0) or (beta_ > 0.5) )
    dserror("beta out of range (0.0,0.5]");
  else
    std::cout << "   beta = " << beta_ << std::endl;
  // gamma
  if ( (gamma_ <= 0.0) or (gamma_ > 1.0) )
    dserror("gamma out of range (0.0,1.0]");
  else
    std::cout << "   gamma = " << gamma_ << std::endl;
  // alpha_f
  if ( (alphaf_ < 0.0) or (alphaf_ >= 1.0) )
    dserror("alpha_f out of range [0.0,1.0)");
  else
    std::cout << "   alpha_f = " << alphaf_ << std::endl;
  // alpha_m
  if ( (alpham_ < -1.0) or (alpham_ >= 1.0) )
    dserror("alpha_m out of range [-1.0,1.0)");
  else
    std::cout << "   alpha_m = " << alpham_ << std::endl;

  // mid-averaging type
  // In principle, there exist two mid-averaging possibilities, TR-like and IMR-like,
  // where TR-like means trapezoidal rule and IMR-like means implicit mid-point rule.
  // We used to maintain implementations of both variants, but due to its significantly
  // higher complexity, the IMR-like version has been deleted (popp 02/2013). The nice
  // thing about TR-like mid-averaging is that all element (and thus also material) calls
  // are exclusively(!) carried out at the end-point t_{n+1} of each time interval, but
  // never explicitly at some generalized midpoint, such as t_{n+1-\alpha_f}. Thus, any
  // cumbersome extrapolation of history variables, etc. becomes obsolete.
  if (midavg_ != INPAR::PARTICLE::midavg_trlike)
    dserror("mid-averaging of internal forces only implemented TR-like");
  else
    std::cout << "   midavg = " << INPAR::PARTICLE::MidAverageString(midavg_)<<std::endl;

  // done
  return;
}

/*----------------------------------------------------------------------*/
void PARTICLE::TimIntGenAlpha::CalcCoeff()
{
  // rho_inf specified --> calculate optimal parameters
  if (rho_inf_!=-1.)
  {
    if ( (rho_inf_ < 0.0) or (rho_inf_ > 1.0) )
      dserror("rho_inf out of range [0.0,1.0]");
    if ( (beta_!=0.25) or (gamma_!=0.5) or (alpham_!=0.5) or (alphaf_!=0.5) )
      dserror("you may only specify RHO_INF or the other four parameters");
    alpham_ = (2.0*rho_inf_-1.0)/(rho_inf_+1.0);
    alphaf_ = rho_inf_/(rho_inf_+1.0);
    beta_   = 0.25*(1.0-alpham_+alphaf_)*(1.0-alpham_+alphaf_);
    gamma_  = 0.5-alpham_+alphaf_;
  }
}


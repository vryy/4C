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
  rho_inf_(particledynparams.sublist("GENALPHA").get<double>("RHO_INF")),
  tol_(particledynparams.sublist("GENALPHA").get<double>("TOL")),
  maxIt_(particledynparams.sublist("GENALPHA").get<int>("MAXIT")),
  mGradW_(Teuchos::null),
  mHessW_(Teuchos::null),
  dism_(Teuchos::null),
  velm_(Teuchos::null),
  accm_(Teuchos::null),
  resAcc_(Teuchos::null)
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
  SetupStateVectors();

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
  // we do not need these vectors
  densityDotn_ = Teuchos::null;
  densityDot_ = Teuchos::null;

}

/*----------------------------------------------------------------------*/
/* Integrate step */
int PARTICLE::TimIntGenAlpha::IntegrateStep()
{
  const double dt = (*dt_)[0];   // \f$\Delta t_{n}\f$

  // update the interaction handler
  interHandler_->Init(stepn_, disn_, veln_, radiusn_, mass_, specEnthalpyn_);

  // heat balance
  interHandler_->Inter_pvp_specEnthalpyDot(specEnthalpyDotn_);
  interHandler_->Inter_pvhs_specEnthalpyDot(specEnthalpyDotn_);
  specEnthalpyn_->Update(dt, *specEnthalpyDotn_, 1.0);

  // compute densities and pressures
  interHandler_->MF_mW(densityn_);
  interHandler_->SetStateVector(densityn_, PARTICLE::Density);
  UpdatePressure();
  interHandler_->SetStateVector(pressure_, PARTICLE::Pressure);

  //-- genAlpha can start --//

  // predict new state
  PredictNewState(true);
  // predict mid state
  PredictMidState();
  // compute the residual

  ResAcc();

  // set up the error
  double resAccNorm2;

  // newton - rhapson iteration
  int nri = 0;
  bool notConverged = true;
  for (; nri<maxIt_; ++ nri)
  {
    // Is convergence reached?
    resAcc_->Norm2(&resAccNorm2);
    std::cout << "Iteration: " << nri << "/" << maxIt_ << std::endl;
    std::cout << "Residual: " << resAccNorm2 << std::endl;
    std::cout << "Required residual: " << tol_ << std::endl;
    if (resAccNorm2 <= tol_)
    {
      std::cout << "--- Converged! ---\n";
      notConverged = false;
      break;
    }


    // update gradient and hessian (necessary for gradResAcc)
    interHandler_->MF_mGradW(mGradW_);
    interHandler_->SetStateVector(mGradW_, PARTICLE::mGradW);
    interHandler_->MF_mHessW(mHessW_);
    interHandler_->SetStateVector(mHessW_);

    // compute the gradient of the residual acceleration
    GradResAcc(dt);

    // compute deltaDis and update disn - Newton-Rhapson iteration
    CorrectDis();

    // update the interaction handler
    interHandler_->SetStateVector(disn_, PARTICLE::Dis);
    // update weights
    interHandler_->UpdateWeights(step_);

    interHandler_->MF_mW(densityn_);
    interHandler_->SetStateVector(densityn_, PARTICLE::Density);
    UpdatePressure();
    interHandler_->SetStateVector(pressure_, PARTICLE::Pressure);
    // predict new state
    PredictNewState();
    // predict mid state
    PredictMidState();
    // compute the residual
    ResAcc();
  }

  if (notConverged)
  {
    std::cout << "Warning! The required tolerance was not reached!\n";
  }


  // erase the handler, information are outdated
  interHandler_->Clear();
  //dserror("IntegrateStep is still in the todo list");

  return 0;
}

/*----------------------------------------------------------------------*/
// overload to determine the initial accelerations
void PARTICLE::TimIntGenAlpha::DetermineMassDampConsistAccel()
{
  // create the wall discretization
  particle_algorithm_->SetUpWallDiscret();

  // set up connections and the local state vectors in the interaction handler
  interHandler_->Init(stepn_, disn_, veln_, radiusn_, mass_, specEnthalpyn_);

  // compute and replace density and pressure
  interHandler_->MF_mW(densityn_);
  interHandler_->SetStateVector(densityn_, PARTICLE::Density);
  UpdatePressure();
  interHandler_->SetStateVector(pressure_, PARTICLE::Pressure);

  interHandler_->Inter_pvp_acc(accn_);
  interHandler_->Inter_pvw_acc(accn_);
  acc_->UpdateSteps(*accn_);

  interHandler_->Clear();
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate mid-state vectors by averaging end-point vectors */
void PARTICLE::TimIntGenAlpha::PredictMidState()
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

/*----------------------------------------------------------------------*/
/* Consistent predictor with constant displacements
 * and consistent velocities and displacements */
void PARTICLE::TimIntGenAlpha::PredictNewState(const bool disAreEqual)
{
  if (disAreEqual)
  {
    veln_->PutScalar(0.0);
    accn_->PutScalar(0.0);
  }
  else
  {
    veln_->Update(1.0, *disn_, -1.0, *(*dis_)(0), 0.0);
    accn_->Update(1.0, *disn_, -1.0, *(*dis_)(0), 0.0);
  }

  // consistent velocities following Newmark formulas
  veln_->Update((beta_-gamma_)/beta_, *(*vel_)(0),
                (2.*beta_-gamma_)*(*dt_)[0]/(2.*beta_), *(*acc_)(0),
                gamma_/(beta_*(*dt_)[0]));

  // consistent accelerations following Newmark formulas
  accn_->Update(-1./(beta_*(*dt_)[0]), *(*vel_)(0),
                (2.*beta_-1.)/(2.*beta_), *(*acc_)(0),
                1./(beta_*(*dt_)[0]*(*dt_)[0]));

  // watch out
  return;
}


/*----------------------------------------------------------------------*/
/* Read and set restart state - overload */
void PARTICLE::TimIntGenAlpha::ReadRestartState()
{
  // call the base function
  TimInt::ReadRestartState();

  // reset to zero the additional state vectors
  SetupStateVectors();

}

/*----------------------------------------------------------------------*/
/* Read and set restart state - overload */
void PARTICLE::TimIntGenAlpha::SetupStateVectors()
{
  mGradW_ = LINALG::CreateVector(*DofRowMapView(), true);
  mHessW_ = Teuchos::rcp(new Epetra_MultiVector(*DofRowMapView(), 3, true));
  accm_ = LINALG::CreateVector(*DofRowMapView(), true);
  velm_ = LINALG::CreateVector(*DofRowMapView(), true);
  dism_ = LINALG::CreateVector(*DofRowMapView(), true);
  resAcc_ = LINALG::CreateVector(*DofRowMapView(), true);
  gradResAcc_ = Teuchos::rcp(new Epetra_MultiVector(*DofRowMapView(), 3, true));
}


/*----------------------------------------------------------------------*/
/* State vectors are updated according to the new distribution of particles */
void PARTICLE::TimIntGenAlpha::UpdateStatesAfterParticleTransfer()
{
  std::cout << "puppa\n";

  // call base function
  TimInt::UpdateStatesAfterParticleTransfer();

  std::cout << "puppa\n";

  UpdateStateVectorMap(mGradW_);

  std::cout << "puppa\n";
  UpdateStateVectorMap(mHessW_);
  std::cout << "puppa\n";
  UpdateStateVectorMap(dism_);
  std::cout << "puppa\n";
  UpdateStateVectorMap(velm_);
  std::cout << "puppa\n";
  UpdateStateVectorMap(accm_);
  std::cout << "puppa\n";
  UpdateStateVectorMap(resAcc_);
  std::cout << "puppa\n";
  UpdateStateVectorMap(gradResAcc_);
  std::cout << "puppa\n";
}


/*----------------------------------------------------------------------*/
/* Compute the residual (acceleration) */
void PARTICLE::TimIntGenAlpha::ResAcc()
{
  // erase the vector
  resAcc_->PutScalar(0.0);
  // build -resAcc
  GravityAcc(resAcc_);
  interHandler_->Inter_pvp_acc(resAcc_);
  interHandler_->Inter_pvw_acc(resAcc_);
  resAcc_->Update(-1.0, *accm_, 1.0);
  // reverse it
  resAcc_->Scale(-1.0);
}


/*----------------------------------------------------------------------*/
/* Consistent predictor with constant displacements
 * and consistent velocities and displacements */
void PARTICLE::TimIntGenAlpha::GradResAcc(const double dt)
{
  // erase the vector
  gradResAcc_->PutScalar(0.0);
  // build
  interHandler_->Inter_pvp_gradAccP(gradResAcc_, restDensity_);
  interHandler_->Inter_pvw_gradAccP(gradResAcc_, restDensity_);

  // add the acceleration part (in a collocation method is quite straight forward)
  const double accCoeff = (1 - alpham_) / (beta_ * dt * dt);

  for (int ii = 0; ii < discret_->DofRowMap()->NumMyElements(); ++ii)
  {
    ((*gradResAcc_)[ii%3])[ii] += accCoeff;
  }

}

/*----------------------------------------------------------------------*/
/* Consistent predictor with constant displacements
 * and consistent velocities and displacements */
void PARTICLE::TimIntGenAlpha::CorrectDis()
{
  for (int lidRowNode = 0; lidRowNode < discret_->NodeRowMap()->NumMyElements(); ++lidRowNode)
  {
    // indexes
    DRT::Node *particle = discret_->lRowNode(lidRowNode);
    std::vector<int> lm;
    lm.reserve(3);
    discret_->Dof(particle, lm);


    LINALG::Matrix<3,3> gradResAcc_i;
    PARTICLE::Utils::ExtractMyValues(*gradResAcc_, gradResAcc_i, lm);
    LINALG::Matrix<3,1> resAcc_i;
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*resAcc_, resAcc_i, lm);
    gradResAcc_i.Invert();

    LINALG::Matrix<3,1> deltaDis;
    deltaDis.MultiplyNN(gradResAcc_i, resAcc_i);
    deltaDis.Scale(-1.0);

    LINALG::Assemble(*disn_, deltaDis, lm, myrank_);
  }
}

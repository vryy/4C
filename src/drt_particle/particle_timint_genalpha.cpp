/*----------------------------------------------------------------------*/
/*!
\file particle_timint_genalpha.cpp

\brief Implicit particle time integration scheme

\level 2

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/
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
#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------*/
/* Constructor */
PARTICLE::TimIntGenAlpha::TimIntGenAlpha(
    const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& particledynparams,
    const Teuchos::ParameterList& xparams,
    Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<LINALG::Solver>& solver,
    Teuchos::RCP<IO::DiscretizationWriter> output
  ) : PARTICLE::TimIntImpl(ioparams, particledynparams, xparams, actdis, output),

  solver_(solver),
  midavg_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::MidAverageEnum>(particledynparams.sublist("GENALPHA"),"GENAVG")),\
  gradResAccApproxType_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::GAapproxType>(particledynparams.sublist("GENALPHA"),"GRADRES_APPROX")),
  beta_(particledynparams.sublist("GENALPHA").get<double>("BETA")),
  gamma_(particledynparams.sublist("GENALPHA").get<double>("GAMMA")),
  alphaf_(particledynparams.sublist("GENALPHA").get<double>("ALPHA_F")),
  alpham_(particledynparams.sublist("GENALPHA").get<double>("ALPHA_M")),
  rho_inf_(particledynparams.sublist("GENALPHA").get<double>("RHO_INF")),
  tol_(particledynparams.sublist("GENALPHA").get<double>("TOL")),
  maxIt_(particledynparams.sublist("GENALPHA").get<int>("MAXIT")),
  mGradW_(Teuchos::null),
  //mHessW_(Teuchos::null),
  dism_(Teuchos::null),
  velm_(Teuchos::null),
  accm_(Teuchos::null),
  resAcc_(Teuchos::null),
  gradResAcc_(Teuchos::null)
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
  interHandler_ = Teuchos::rcp(new PARTICLE::ParticleMeshFreeInteractionHandler(discret_, particle_algorithm_, particleparams, initDensity_, restDensity_, refdensfac_));

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

  //Perform some compatibility checks:
  const INPAR::PARTICLE::WallInteractionType wallInteractionType=DRT::INPUT::IntegralValue<INPAR::PARTICLE::WallInteractionType>(DRT::Problem::Instance()->ParticleParams(),"WALL_INTERACTION_TYPE");
  if(wallInteractionType==INPAR::PARTICLE::BoundarParticle_NoSlip or wallInteractionType==INPAR::PARTICLE::BoundarParticle_FreeSlip)
    dserror("Boundary particles are not possible in particle_timint_genalpha so far!");

  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"SOLVE_THERMAL_PROBLEM")==false)
    dserror("The avoidance of the thermal problem has not been tested in particle_timint_genalpha so far!");

  if(DRT::Problem::Instance()->ParticleParams().get<double>("VISCOUS_DAMPING")>0)
    dserror("The application of a viscous damping force has not been tested in particle_timint_genalpha so far!");


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

  // compute and set displacement-related quantities
  ComputeAndSetDisRelatedStateVectors();

  //-- genAlpha can start --//

  // predict new state
  PredictNewState(disn_, true);
  // predict mid state
  PredictMidState(disn_);

  // compute the residual
  ResAcc();

  // set up the error
  double resAccAvgNorm2;

  // newton - rhapson iteration
  bool notConverged = true;
  int nri = 0;
  for (; nri<maxIt_; ++ nri)
  {
    // Is convergence reached?
    resAcc_->Norm2(&resAccAvgNorm2);

    resAccAvgNorm2 /= discret_->NodeRowMap()->NumGlobalElements();
    if (myrank_ == 0)
    {
      std::cout << "Iteration: " << nri << "/" << maxIt_ << " and residual: " << resAccAvgNorm2 << "\n";
    }

    if (resAccAvgNorm2 <= tol_)
    {
      if (myrank_ == 0)
      {
        std::cout << "--- Converged! ---\n";
      }
      notConverged = false;
      break;
    }

    // compute the gradient of the residual acceleration
    GradResAcc();

    // predict new state
    PredictNewState(disn_);
    // predict mid state
    PredictMidState(disn_);


    // compute deltaDis and update disn - Newton-Rhapson iteration. Check dis every nn iterations
    CorrectDis((10 * nri) % maxIt_ == 0);

    // update the interaction handler
    interHandler_->SetStateVector(disn_, PARTICLE::Dis);
    // update weights
    interHandler_->UpdateWeights(step_);

    // update the displacement-related quantities
    ComputeAndSetDisRelatedStateVectors();

    // predict new state
    PredictNewState(disn_);

    // predict mid state
    PredictMidState(disn_);

    // compute the residual
    ResAcc();
  }

  if (notConverged)
  {
    if (myrank_ == 0)
    {
      double maxValue, minValue;
      resAcc_->MaxValue(&maxValue);
      resAcc_->MinValue(&minValue);
      std::cout << "\nWarning! The required tolerance was not reached! Final residual: " << resAccAvgNorm2 << "\n";
      std::cout << "Final residual extremes: " << maxValue << " " << minValue << "\n";
    }
    std::cin.get();
  }

  // erase the handler, information are outdated
  interHandler_->Clear();

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

  interHandler_->Inter_pvp_acc_var1(accn_);
  interHandler_->Inter_pvw_acc(accn_);

  acc_->UpdateSteps(*accn_);

  interHandler_->Clear();
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate mid-state vectors by averaging end-point vectors */
void PARTICLE::TimIntGenAlpha::PredictMidState(const Teuchos::RCP<Epetra_Vector> disn)
{
  // mid-displacements D_{n+1-alpha_f} (dism)
  //    D_{n+1-alpha_f} := (1.-alphaf) * D_{n+1} + alpha_f * D_{n}
  dism_->Update(1.-alphaf_, *disn, alphaf_, (*dis_)[0], 0.0);

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
void PARTICLE::TimIntGenAlpha::PredictNewState(const Teuchos::RCP<Epetra_Vector> disn, const bool disAreEqual)
{
  if (disAreEqual)
  {
    veln_->PutScalar(0.0);
    accn_->PutScalar(0.0);
  }
  else
  {
    veln_->Update(1.0, *disn, -1.0, *(*dis_)(0), 0.0);
    accn_->Update(1.0, *disn, -1.0, *(*dis_)(0), 0.0);
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

  IO::DiscretizationReader reader(discret_, step_);
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

  // set up the pressure
  PARTICLE::Utils::Density2Pressure(restDensity_,refdensfac_,densityn_,specEnthalpyn_,pressure_,particle_algorithm_->ExtParticleMat(),true);

  // reset to zero the additional state vectors
  SetupStateVectors();

}

/*----------------------------------------------------------------------*/
/* Read and set restart state - overload */
void PARTICLE::TimIntGenAlpha::SetupStateVectors()
{
  mGradW_ = LINALG::CreateVector(*DofRowMapView(), true);
  //mHessW_ = Teuchos::rcp(new Epetra_MultiVector(*DofRowMapView(), 3, true));
  accm_ = LINALG::CreateVector(*DofRowMapView(), true);
  velm_ = LINALG::CreateVector(*DofRowMapView(), true);
  dism_ = LINALG::CreateVector(*DofRowMapView(), true);
  resAcc_ = LINALG::CreateVector(*DofRowMapView(), true);
}


/*----------------------------------------------------------------------*/
/* State vectors are updated according to the new distribution of particles */
void PARTICLE::TimIntGenAlpha::UpdateStatesAfterParticleTransfer()
{
  // call base function
  TimInt::UpdateStatesAfterParticleTransfer();

  UpdateStateVectorMap(mGradW_);
  //UpdateStateVectorMap(mHessW_);
  UpdateStateVectorMap(dism_);
  UpdateStateVectorMap(velm_);
  UpdateStateVectorMap(accm_);
  UpdateStateVectorMap(resAcc_);
}


/*----------------------------------------------------------------------*/
/* Compute the residual (acceleration) */
void PARTICLE::TimIntGenAlpha::ResAcc()
{
  // erase the vector
  resAcc_->PutScalar(0.0);

  // build
  // F_ext - gravity
  //TODO: Time ramp for gravity forces not considered here so far!
  GravityAcc(resAcc_, -1.0, -1.0);
    // F_int - P
  interHandler_->Inter_pvp_acc_var1(resAcc_, Teuchos::null, -1.0);
  interHandler_->Inter_pvw_acc(resAcc_, -1.0);
    // Acc - Am
  resAcc_->Update(1.0, *accm_, 1.0);
}


/*----------------------------------------------------------------------*/
/* Consistent predictor with constant displacements
 * and consistent velocities and displacements */
void PARTICLE::TimIntGenAlpha::GradResAcc()
{
  const double dt = (*dt_)[0];
  // erase the vector
  gradResAcc_ = Teuchos::rcp(new LINALG::SparseMatrix(*(discret_->DofRowMap()), 0, false, false));

  // build
    // F_int - P
  const double accPmulti = 1 - alphaf_;

  if (gradResAccApproxType_ == INPAR::PARTICLE::gaapprox_full)
  {
    interHandler_->Inter_pvp_gradAccP(gradResAcc_, restDensity_, -accPmulti);
    interHandler_->Inter_pvw_gradAccP(gradResAcc_, restDensity_, -accPmulti);
  }
  else
  {
    interHandler_->Inter_pvp_gradAccPapproxOnlyHess(gradResAcc_, restDensity_, -accPmulti);
    interHandler_->Inter_pvw_gradAccPapproxOnlyHess(gradResAcc_, restDensity_, -accPmulti);
  }
    // Acc -Am
  double accCoeff = (1 - alpham_) / (beta_ * dt * dt);
  for (int lidDofRow = 0; lidDofRow < discret_->DofRowMap()->NumMyElements(); ++lidDofRow)
  {
    int gidDof = discret_->DofRowMap()->GID(lidDofRow);
    gradResAcc_->Assemble(accCoeff, gidDof, gidDof);
  }

  // the matrix has been filled
  gradResAcc_->Complete();

}

/*----------------------------------------------------------------------*/
/* Consistent predictor with constant displacements
 * and consistent velocities and displacements */
void PARTICLE::TimIntGenAlpha::CorrectDis(bool checkDis)
{

  Teuchos::RCP<Epetra_Vector> deltaDisNorms2 = Teuchos::null;

  if (checkDis)
  {
    deltaDisNorms2 = LINALG::CreateVector(*NodeRowMapView(), true);
  }

  gradResAcc_->Scale(-1.0);

  Teuchos::RCP<Epetra_Vector> deltaDis = LINALG::CreateVector(*DofRowMapView(), true);
  solver_->Solve(gradResAcc_->EpetraOperator(), deltaDis, resAcc_, true, true);

  disn_->Update(1.0, *deltaDis, 1.0);

  if (checkDis)
  {


    double maxDeltaDisNorms2 = 0;
    double minRadius = 0;
    deltaDisNorms2->MaxValue(&maxDeltaDisNorms2);
    radiusn_->MinValue(&minRadius);
    if (2 * maxDeltaDisNorms2 > minRadius)
    {
      dserror("Particles go too fast! They can bypass each other. Reduce the time step");
    }
  }
}

/*----------------------------------------------------------------------*/
/* Compute and set dis-related state vectors */
void PARTICLE::TimIntGenAlpha::ComputeAndSetDisRelatedStateVectors()
{
  interHandler_->MF_mW(densityn_);
  interHandler_->SetStateVector(densityn_, PARTICLE::Density);
  UpdatePressure();
  interHandler_->SetStateVector(pressure_, PARTICLE::Pressure);
  if (gradResAccApproxType_ == INPAR::PARTICLE::gaapprox_full)
  {
    interHandler_->MF_mGradW(mGradW_);
    interHandler_->SetStateVector(mGradW_, PARTICLE::mGradW);
  }
  // these lines can be useful in the future
  //interHandler_->MF_mHessW(mHessW_);
  //interHandler_->SetStateVector(mHessW_);
}

/*----------------------------------------------------------------------*/
// Compute the gradGenAcc with the finite differences
Teuchos::RCP<LINALG::SparseMatrix> PARTICLE::TimIntGenAlpha::FDGradResAcc()
{
  const double dt = (*dt_)[0];

  // compute the residual difference
  Teuchos::RCP<Epetra_Vector> resAcc0 = Teuchos::rcp(new Epetra_Vector(*resAcc_));
  Teuchos::RCP<LINALG::SparseMatrix> gradResAcc = Teuchos::rcp(new LINALG::SparseMatrix(*(discret_->DofRowMap()), 0, false, false));

  for (int jj = 0; jj<discret_->DofRowMap()->NumMyElements(); ++jj)
  {
    Teuchos::RCP<Epetra_Vector> newDis = Teuchos::rcp(new Epetra_Vector(*disn_));
    const double incr = dt * dt * std::abs((*resAcc0)[jj]) + dt * dt;

    (*newDis)[jj] += incr;

    interHandler_->SetStateVector(newDis, PARTICLE::Dis);
    interHandler_->UpdateWeights(step_);

    interHandler_->MF_mW(densityn_);
    interHandler_->SetStateVector(densityn_, PARTICLE::Density);
    UpdatePressure();
    interHandler_->SetStateVector(pressure_, PARTICLE::Pressure);
    interHandler_->MF_mGradW(mGradW_);
    interHandler_->SetStateVector(mGradW_, PARTICLE::mGradW);
    PredictNewState(newDis);
    PredictMidState(newDis);
    ResAcc();

    Teuchos::RCP<Epetra_Vector> deltaResAcc = Teuchos::rcp(new Epetra_Vector(*resAcc0));

    deltaResAcc->Update(1.0, *resAcc_, -1.0);
    deltaResAcc->Scale(1/incr);

    int gidj = discret_->DofRowMap()->GID(jj);
    for (int ii = 0; ii<discret_->DofRowMap()->NumMyElements(); ++ii)
    {
      int gidi = discret_->DofRowMap()->GID(ii);
      gradResAcc->Assemble((*deltaResAcc)[ii], gidi, gidj);
    }
  }

  gradResAcc->Complete();

  resAcc_ = resAcc0;

  return gradResAcc;
}

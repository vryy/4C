/*----------------------------------------------------------------------*/
/*!
\file strutimint_genalpha.cpp
\brief Structural time integration with generalised-alpha

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "strutimint_genalpha.H"

/*----------------------------------------------------------------------*/
/* convert input string to enum for mid-average type */
enum StruTimIntGenAlpha::MidAverageEnum StruTimIntGenAlpha::MapMidAvgStringToEnum
(
  const std::string name
)
{
  if (name == "Vague")
  {
    return midavg_vague;
  }
  else if (name == "ImrLike")
  {    
    return midavg_imrlike;
  }
  else if (name == "TrLike")
  {
    return midavg_trlike;
  }
  else
  {
    return midavg_vague;
  }
}

/*----------------------------------------------------------------------*/
/* constructor */
StruTimIntGenAlpha::StruTimIntGenAlpha
(
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& sdynparams,
  const Teuchos::ParameterList& xparams,
  const Teuchos::ParameterList& genalphaparams,
  DRT::Discretization& actis,
  LINALG::Solver& solver,
  IO::DiscretizationWriter& output
)
: StruTimIntImpl
  (
    ioparams,
    sdynparams,
    xparams,
    actis,
    solver,
    output
  ),
  midavg_(MapMidAvgStringToEnum(genalphaparams.get<string>("GENAVG"))),
  iterupditer_(false),
  beta_(genalphaparams.get<double>("BETA")),
  gamma_(genalphaparams.get<double>("GAMMA")),
  alphaf_(genalphaparams.get<double>("ALPHA_F")),
  alpham_(genalphaparams.get<double>("ALPHA_M")),
  dism_(Teuchos::null),
  velm_(Teuchos::null),
  accm_(Teuchos::null),
  fint_(Teuchos::null),
  fintm_(Teuchos::null),
  fintn_(Teuchos::null),
  fext_(Teuchos::null),
  fextm_(Teuchos::null),
  fextn_(Teuchos::null),
  finertm_(Teuchos::null),
  fviscm_(Teuchos::null),
  frobin_(Teuchos::null)
{
  // info to user : OST --- your oriental scheme
  if (myrank_ == 0)
  {
    std::cout << "with generalised-alpha ("
              << "beta=" << beta_
              << ", gamma=" << gamma_
              << ", alpha_f=" << alphaf_
              << ", alpha_m=" << alpham_
              << ")" << std::endl << std::endl;
  }

  // create state vectors
  
  // mid-displacements
  dism_ = LINALG::CreateVector(*dofrowmap_, true);
  // mid-velocities
  velm_ = LINALG::CreateVector(*dofrowmap_, true);
  // mid-accelerations
  accm_ = LINALG::CreateVector(*dofrowmap_, true);

  // create force vectors

  // internal forces
  if (midavg_ == midavg_trlike)
  {
    // internal force vector F_{int;n} at last time
    fint_ = LINALG::CreateVector(*dofrowmap_, true);
    // internal force vector F_{int;n+1} at new time
    fintn_ = LINALG::CreateVector(*dofrowmap_, true);
    // set initial internal force vector
    ApplyForceStiffInternal(time_, dis_, zeros_, fint_, Teuchos::null);
  } 
  else if (midavg_ == midavg_imrlike)
  {
    // internal force vector F_{int;m} at mid-time
    fintm_ = LINALG::CreateVector(*dofrowmap_, true);
  }

  // external force vector F_ext at last times
  fext_ = LINALG::CreateVector(*dofrowmap_, true);
  // external mid-force vector F_{ext;n+1-alpha_f}
  fextm_ = LINALG::CreateVector(*dofrowmap_, true);
  // external force vector F_{n+1} at new time
  fextn_ = LINALG::CreateVector(*dofrowmap_, true);
  // set initial external force vector
  ApplyForceExternal(time_, dis_, fext_);

  // inertial mid-point force vector F_inert
  finertm_ = LINALG::CreateVector(*dofrowmap_, true);
  // viscous mid-point force vector F_visc
  fviscm_ = LINALG::CreateVector(*dofrowmap_, true);

  // external pseudo force due to RobinBC
  frobin_ = LINALG::CreateVector(*dofrowmap_, true);

  // have a nice day
  return;
}

/*----------------------------------------------------------------------*/
/* Consistent predictor with constant displacements
 * and consistent velocities and displacements */
void StruTimIntGenAlpha::PredictConstDisConsistVelAcc()
{
  // constant predictor : displacement in domain
  disn_->Update(1.0, *dis_, 0.0);
  /*
  // apply Dirichlet condition onto this
  // This is not absolutely necessary, as the DBCs are set
  // by StruTimIntImpl::Predict. 
  // However, here they are
  // applied to achieve a more better guess for the remaining
  // DOFs.
  ApplyDirichletBC(timen_, disn_, Teuchos::null, Teuchos::null);
  */

  // consistent velocities
  veln_->Update(1.0, *disn_, -1.0, *dis_, 0.0);
  veln_->Update((beta_-gamma_)/beta_, *vel_,
                (2.*beta_-gamma_)*dt_/(2.*beta_), *acc_,
                gamma_/(beta_*dt_));

  // consistent accelerations
  accn_->Update(1.0, *disn_, -1.0, *dis_, 0.0);
  accn_->Update(-1./(beta_*dt_), *vel_,
                (2.*beta_-1.)/(2.*beta_), *acc_,
                1./(beta_*dt_*dt_));
  
  // watch out
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate residual force and its stiffness, ie derivative
 * with respect to end-point displacements \f$D_{n+1}\f$ */
void StruTimIntGenAlpha::EvaluateForceStiffResidual()
{
  // build by last converged state and predicted target state
  // the predicted mid-state
  EvaluateMidState();

  // build new external forces
  fextn_->PutScalar(0.0);
  ApplyForceExternal(timen_, dis_, fextn_);

  // external mid-forces F_{ext;n+1-alpha_f} (fextm)
  //    F_{ext;n+1-alpha_f} := (1.-alphaf) * F_{ext;n+1}
  //                         + alpha_f * F_{ext;n}
  fextm_->Update(1.-alphaf_, *fextn_, alphaf_, *fext_, 0.0);

  // initialise internal forces
  if (midavg_ == midavg_trlike)
  {
    fintn_->PutScalar(0.0);
  } 
  else if (midavg_ == midavg_imrlike)
  {
    fintm_->PutScalar(0.0);
  }

  // initialise stiffness matrix to zero
  stiff_->Zero();

  // ordinary internal force and stiffness
  if (midavg_ == midavg_trlike)
  {
    ApplyForceStiffInternal(timen_, disn_, disi_,  fintn_, stiff_);
  } 
  else if (midavg_ == midavg_imrlike)
  {
    disi_->Scale(1.-alphaf_);
    ApplyForceStiffInternal(timen_, dism_, disi_,  fintm_, stiff_);
  }

  // surface stress force
  if (midavg_ == midavg_trlike)
  {
    ApplyForceStiffSurfstress(disn_, fintn_, stiff_);
  } 
  else if (midavg_ == midavg_imrlike)
  {
    ApplyForceStiffSurfstress(dism_, fintm_, stiff_);
  }

  // potential forces
  if (midavg_ == midavg_trlike)
  {
    ApplyForceStiffPotential(disn_, fintn_, stiff_);
  } 
  else if (midavg_ == midavg_imrlike)
  {
    ApplyForceStiffPotential(dism_, fintm_, stiff_);
  }

  // inertial forces #finertm_
  mass_->Multiply(false, *accm_, *finertm_);

  // viscous forces due Rayleigh damping
  if (damping_)
  {
    damp_->Multiply(false, *velm_, *fviscm_);
  }

  // build negative residual
  //    Res = -( M . A_{n+1-alpha_m}
  //             + C . V_{n+1-alpha_f}
  //             + F_{int;m}
  //             - F_{ext;n+1-alpha_f} )
  fres_->Update(1.0, *fextm_, 0.0);
  if (midavg_ == midavg_trlike)
  {
    fres_->Update(-(1.-alphaf_), *fintn_, -alphaf_, *fint_, 1.0);
  }
  else if (midavg_ == midavg_imrlike)
  {
    fres_->Update(-1.0, *fintm_, 1.0);
  }
  if (damping_)
  {
    fres_->Update(-1.0, *fviscm_, 1.0);
  }
  fres_->Update(-1.0, *finertm_, 1.0);
  
  // build tangent matrix : effective dynamic stiffness matrix
  //    K_{Teffdyn} = (1 - alpha_m)/(beta*dt^2) M
  //                + (1 - alpha_f)*y/(beta*dt) C     
  //                + (1 - alpha_f) K_{T}
  stiff_->Add(*mass_, false, (1.-alpham_)/(beta_*dt_*dt_), 1.-alphaf_);
  if (damping_)
  {
    stiff_->Add(*damp_, false, (1.-alphaf_)*gamma_/(beta_*dt_), 1.0);
  }
  stiff_->Complete();  // close stiffness matrix

  // hallelujah
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate mid-state vectors by averaging end-point vectors */
void StruTimIntGenAlpha::EvaluateMidState()
{
  // mid-displacements D_{n+1-alpha_f} (dism)
  //    D_{n+1-alpha_f} := (1.-alphaf) * D_{n+1} + alpha_f * D_{n}
  dism_->Update(1.-alphaf_, *disn_, alphaf_, *dis_, 0.0);
  
  // mid-velocities V_{n+1-alpha_f} (velm)
  //    V_{n+1-alpha_f} := (1.-alphaf) * V_{n+1} + alpha_f * V_{n}
  velm_->Update(1.-alphaf_, *veln_, alphaf_, *vel_, 0.0);
  
  // mid-accelerations A_{n+1-alpha_m} (accm)
  //    A_{n+1-alpha_m} := (1.-alpha_m) * A_{n+1} + alpha_m * A_{n}
  accm_->Update(1.-alpham_, *accn_, alpham_, *acc_, 0.0);

  // jump
  return;
}

/*----------------------------------------------------------------------*/
/* calculate characteristic/reference norms for displacements
 * originally by lw */
double StruTimIntGenAlpha::CalcRefNormDisplacement()
{
  // The reference norms are used to scale the calculated iterative
  // displacement norm and/or the residual force norm. For this
  // purpose we only need the right order of magnitude, so we don't
  // mind evaluating the corresponding norms at possibly different
  // points within the timestep (end point, generalized midpoint).

  double charnormdis = 0.0;
  dis_->Norm2(&charnormdis);

  // rise your hat
  return charnormdis;
}

/*----------------------------------------------------------------------*/
/* calculate characteristic/reference norms for forces
 * originally by lw */
double StruTimIntGenAlpha::CalcRefNormForce()
{
  // The reference norms are used to scale the calculated iterative
  // displacement norm and/or the residual force norm. For this
  // purpose we only need the right order of magnitude, so we don't
  // mind evaluating the corresponding norms at possibly different
  // points within the timestep (end point, generalized midpoint).

  // norm of the internal forces
  double fintnorm = 0.0;
  if (midavg_ == midavg_trlike)
  {
    fintn_->Norm2(&fintnorm);
  }
  else if (midavg_ == midavg_imrlike)
  {
    fintm_->Norm2(&fintnorm);
  }

  // norm of the external forces
  double fextnorm = 0.0;
  fextm_->Norm2(&fextnorm);

  // norm of the inertial forces
  double finertnorm = 0.0;
  finertm_->Norm2(&finertnorm);

  // norm of viscous forces
  double fviscnorm = 0.0;
  if (damping_)
  {
    fviscm_->Norm2(&fviscnorm);
  }

  // determine worst value ==> charactersitic norm
  return max(fviscnorm, max(finertnorm, max(fintnorm, fextnorm)));
}

/*----------------------------------------------------------------------*/
/* iteration update of state */
void StruTimIntGenAlpha::UpdateIter()
{
  // iterative update method
  if (iterupditer_)
  {
    // new end-point displacements
    // D_{n+1}^{<k+1>} := D_{n+1}^{<k>} + IncD_{n+1}^{<k>}
    disn_->Update(1.0, *disi_, 1.0);
    
    // new end-point velocities
    veln_->Update(gamma_/(beta_*dt_), *disi_, 1.0);
    
    // new end-point accelerations
    accn_->Update(1.0/(beta_*dt_*dt_), *disi_, 1.0);
  }
  // incremental update method;
  else
  {
    // new end-point displacements
    // D_{n+1}^{<k+1>} := D_{n+1}^{<k>} + IncD_{n+1}^{<k>}
    disn_->Update(1.0, *disi_, 1.0);

    // new end-point velocities
    veln_->Update(1.0, *disn_, -1.0, *dis_, 0.0);
    veln_->Update((beta_-gamma_)/beta_, *vel_,
                  (2.0*beta_-gamma_)*dt_/(2.0*beta_), *acc_,
                  gamma_/(beta_*dt_));
    
    // new end-point accelerations
    accn_->Update(1.0, *disn_, -1.0, *dis_, 0.0);
    accn_->Update(-1.0/(beta_*dt_), *vel_,
                  (2.0*beta_-1.0)/(2.0*beta_), *acc_,
                  1.0/(beta_*dt_*dt_));
  }

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/* update after time step */
void StruTimIntGenAlpha::UpdateStep()
{
  // update state
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1}
  dis_->Update(1.0, *disn_, 0.0);
  // new velocities at t_{n+1} -> t_n
  //    V_{n} := V_{n+1}
  vel_->Update(1.0, *veln_, 0.0);
  // new accelerations at t_{n+1} -> t_n
  //    A_{n} := A_{n+1}
  acc_->Update(1.0, *accn_, 0.0);

  // update new external force
  //    F_{ext;n} := F_{ext;n+1}
  fext_->Update(1.0, *fextn_, 0.0);

  // update new internal force
  //    F_{int;n} := F_{int;n+1}
  if (midavg_ == midavg_trlike)
  {
    fint_->Update(1.0, *fintn_, 0.0);
  }

  // update anything that needs to be updated at the element level
  {
    // create the parameters for the discretization
    ParameterList p;
    // other parameters that might be needed by the elements
    p.set("total time", timen_);
    p.set("delta time", dt_);
    p.set("alpha f", alphaf_);
    // action for elements
    if (midavg_ == midavg_trlike) 
    {
      p.set("action", "calc_struct_update_istep");    
    }
    else if (midavg_ == midavg_imrlike)
    {
      p.set("action", "calc_struct_update_genalpha_imrlike");
    }
    // go to elements
    discret_.Evaluate(p, null, null, null, null, null);
  }

  // update surface stress
  UpdateStepSurfstress();

  // update potential forces
  UpdateStepPotential();

  // look out
  return;
}


/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET

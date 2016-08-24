 /*----------------------------------------------------------------------*/
/*!
\file strtimint_genalpha.cpp
\brief Structural time integration with generalised-alpha

<pre>
\level 1

\maintainer Alexander Popp
popp@lnm.mw.tum.de
http://www.lnm.mw.tum.de
089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint_genalpha.H"
#include "stru_aux.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_locsys.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io_pstream.H"
#include "../drt_crack/crackUtils.H"
#include "../drt_plastic_ssn/plastic_ssn_manager.H"

/*----------------------------------------------------------------------*/
void STR::TimIntGenAlpha::CalcCoeff()
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
void STR::TimIntGenAlpha::VerifyCoeff()
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
  if (midavg_ != INPAR::STR::midavg_trlike)
    dserror("mid-averaging of internal forces only implemented TR-like");
  else
    std::cout << "   midavg = " << INPAR::STR::MidAverageString(midavg_)<<std::endl;

  // done
  return;
}

/*----------------------------------------------------------------------*/
/* constructor */
STR::TimIntGenAlpha::TimIntGenAlpha
(
  const Teuchos::ParameterList& timeparams,
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& sdynparams,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization> actdis,
  Teuchos::RCP<LINALG::Solver> solver,
  Teuchos::RCP<LINALG::Solver> contactsolver,
  Teuchos::RCP<IO::DiscretizationWriter> output
)
: TimIntImpl
  (
    timeparams,
    ioparams,
    sdynparams,
    xparams,
    actdis,
    solver,
    contactsolver,
    output
  ),
  midavg_(DRT::INPUT::IntegralValue<INPAR::STR::MidAverageEnum>(sdynparams.sublist("GENALPHA"),"GENAVG")),
  beta_(sdynparams.sublist("GENALPHA").get<double>("BETA")),
  gamma_(sdynparams.sublist("GENALPHA").get<double>("GAMMA")),
  alphaf_(sdynparams.sublist("GENALPHA").get<double>("ALPHA_F")),
  alpham_(sdynparams.sublist("GENALPHA").get<double>("ALPHA_M")),
  rho_inf_(sdynparams.sublist("GENALPHA").get<double>("RHO_INF")),
  dism_(Teuchos::null),
  velm_(Teuchos::null),
  accm_(Teuchos::null),
  fint_(Teuchos::null),
  fintm_(Teuchos::null),
  fintn_(Teuchos::null),
  fext_(Teuchos::null),
  fextm_(Teuchos::null),
  fextn_(Teuchos::null),
  finert_(Teuchos::null),
  finertm_(Teuchos::null),
  finertn_(Teuchos::null),
  fviscm_(Teuchos::null),
  fint_str_(Teuchos::null)
{
  // calculate time integration parameters
  CalcCoeff();

  // info to user about current time integration scheme and its parametrization
  if (myrank_ == 0)
  {
    IO::cout << "with generalised-alpha" << IO::endl;
    VerifyCoeff();

    std::cout << "   p_dis = " << MethodOrderOfAccuracyDis() << std::endl
              << "   p_vel = " << MethodOrderOfAccuracyVel() << std::endl
              << std::endl;
  }

  if (!HaveNonlinearMass())
  {
    // determine mass, damping and initial accelerations
    DetermineMassDampConsistAccel();
  }
  else
  {
    /* the case of nonlinear inertia terms works so far only for examples with
     * vanishing initial accelerations, i.e. the initial external
     * forces and initial velocities have to be chosen consistently!!!
     */
    (*acc_)(0)->PutScalar(0.0);
  }

  // create state vectors

  // mid-displacements
  dism_ = LINALG::CreateVector(*DofRowMapView(), true);
  // mid-velocities
  velm_ = LINALG::CreateVector(*DofRowMapView(), true);
  // mid-accelerations
  accm_ = LINALG::CreateVector(*DofRowMapView(), true);

  // create force vectors

  // internal force vector F_{int;n} at last time
  fint_ = LINALG::CreateVector(*DofRowMapView(), true);
  // internal mid-force vector F_{int;n+1-alpha_f}
  fintm_ = LINALG::CreateVector(*DofRowMapView(), true);
  // internal force vector F_{int;n+1} at new time
  fintn_ = LINALG::CreateVector(*DofRowMapView(), true);

  // external force vector F_ext at last times
  fext_ = LINALG::CreateVector(*DofRowMapView(), true);
  // external mid-force vector F_{ext;n+1-alpha_f}
  fextm_ = LINALG::CreateVector(*DofRowMapView(), true);
  // external force vector F_{n+1} at new time
  fextn_ = LINALG::CreateVector(*DofRowMapView(), true);
  // set initial external force vector
  ApplyForceExternal((*time_)[0], (*dis_)(0), disn_, (*vel_)(0), fext_);

  // inertial force vector F_{int;n} at last time
  finert_ = LINALG::CreateVector(*DofRowMapView(), true);
  // inertial mid-force vector F_{int;n+1-alpha_f}
  finertm_ = LINALG::CreateVector(*DofRowMapView(), true);
  // inertial force vector F_{int;n+1} at new time
  finertn_ = LINALG::CreateVector(*DofRowMapView(), true);

  // viscous mid-point force vector F_visc
  fviscm_ = LINALG::CreateVector(*DofRowMapView(), true);

  // structural rhs for newton line search
  if (fresn_str_!=Teuchos::null)
    fint_str_=LINALG::CreateVector(*DofRowMapView(), true);

  // create parameter list
  Teuchos::ParameterList params;

  // for line search
  if (fintn_str_!=Teuchos::null)
  {
    params.set("cond_rhs_norm",0.);
    params.set("MyPID",myrank_);
  }

  // add initial forces due to 0D cardiovascular coupling conditions - needed in case of initial ventricular pressure!
  Teuchos::ParameterList pwindk;
  pwindk.set("scale_timint", 1.0);
  pwindk.set("time_step_size", (*dt_)[0]);
  ApplyForceStiffCardiovascular0D((*time_)[0], (*dis_)(0), fint_, stiff_, pwindk);

  if (HaveNonlinearMass()==INPAR::STR::ml_none)
  {
    // set initial internal force vector
    ApplyForceStiffInternal((*time_)[0], (*dt_)[0], (*dis_)(0), zeros_, (*vel_)(0), fint_, stiff_,params);
  }
  else
  {
    double timeintfac_dis=beta_*(*dt_)[0]*(*dt_)[0];
    double timeintfac_vel=gamma_*(*dt_)[0];

    // Check, if initial residuum really vanishes for acc_ = 0
    ApplyForceStiffInternalAndInertial((*time_)[0], (*dt_)[0], timeintfac_dis, timeintfac_vel, (*dis_)(0), zeros_, (*vel_)(0), (*acc_)(0), fint_, finert_, stiff_, mass_,params,beta_,gamma_,alphaf_,alpham_);

    NonlinearMassSanityCheck(fext_, (*dis_)(0), (*vel_)(0), (*acc_)(0), &sdynparams);

    if(HaveNonlinearMass() == INPAR::STR::ml_rotations and !SolelyBeam3Elements(actdis))
    {
      dserror("Multiplicative Gen-Alpha time integration scheme only implemented for beam elements so far!");
    }
  }

  // init old time step value
  if (fintn_str_!=Teuchos::null)
    fint_str_->Update(1.,*fintn_str_,0.);

  // have a nice day
  return;
}

/*----------------------------------------------------------------------*/
/* Consistent predictor with constant displacements
 * and consistent velocities and displacements */
void STR::TimIntGenAlpha::PredictConstDisConsistVelAcc()
{
  // constant predictor : displacement in domain
  disn_->Update(1.0, *(*dis_)(0), 0.0);

  // consistent velocities following Newmark formulas
  veln_->Update(1.0, *disn_, -1.0, *(*dis_)(0), 0.0);
  veln_->Update((beta_-gamma_)/beta_, *(*vel_)(0),
                (2.*beta_-gamma_)*(*dt_)[0]/(2.*beta_), *(*acc_)(0),
                gamma_/(beta_*(*dt_)[0]));

  // consistent accelerations following Newmark formulas
  accn_->Update(1.0, *disn_, -1.0, *(*dis_)(0), 0.0);
  accn_->Update(-1./(beta_*(*dt_)[0]), *(*vel_)(0),
                (2.*beta_-1.)/(2.*beta_), *(*acc_)(0),
                1./(beta_*(*dt_)[0]*(*dt_)[0]));

  // watch out
  return;
}

/*----------------------------------------------------------------------*/
/* Consistent predictor with constant velocities,
 * extrapolated displacements and consistent accelerations */
void STR::TimIntGenAlpha::PredictConstVelConsistAcc()
{
  // extrapolated displacements based upon constant velocities
  // d_{n+1} = d_{n} + dt * v_{n}
  disn_->Update(1.0, (*dis_)[0], (*dt_)[0], (*vel_)[0], 0.0);

  // consistent velocities following Newmark formulas
  veln_->Update(1.0, *disn_, -1.0, *(*dis_)(0), 0.0);
  veln_->Update((beta_-gamma_)/beta_, *(*vel_)(0),
                (2.*beta_-gamma_)*(*dt_)[0]/(2.*beta_), *(*acc_)(0),
                gamma_/(beta_*(*dt_)[0]));

  // consistent accelerations following Newmark formulas
  accn_->Update(1.0, *disn_, -1.0, *(*dis_)(0), 0.0);
  accn_->Update(-1./(beta_*(*dt_)[0]), *(*vel_)(0),
                (2.*beta_-1.)/(2.*beta_), *(*acc_)(0),
                1./(beta_*(*dt_)[0]*(*dt_)[0]));

  // That's it!
  return;
}

/*----------------------------------------------------------------------*/
/* Consistent predictor with constant accelerations
 * and extrapolated velocities and displacements */
void STR::TimIntGenAlpha::PredictConstAcc()
{
  // extrapolated displacements based upon constant accelerations
  // d_{n+1} = d_{n} + dt * v_{n} + dt^2 / 2 * a_{n}
  disn_->Update(1.0, (*dis_)[0], (*dt_)[0], (*vel_)[0], 0.0);
  disn_->Update((*dt_)[0] * (*dt_)[0] /2., (*acc_)[0], 1.0);

  // extrapolated velocities (equal to consistent velocities)
  // v_{n+1} = v_{n} + dt * a_{n}
  veln_->Update(1.0, (*vel_)[0], (*dt_)[0], (*acc_)[0],  0.0);

  // constant accelerations (equal to consistent accelerations)
  accn_->Update(1.0, (*acc_)[0], 0.0);

  // That's it!
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate residual force and its stiffness, i.e. derivative
 * with respect to end-point displacements \f$D_{n+1}\f$ */
void STR::TimIntGenAlpha::EvaluateForceStiffResidual(Teuchos::ParameterList& params)
{
  // get info about prediction step from parameter list
  bool predict = false;
  if(params.isParameter("predict"))
    predict = params.get<bool>("predict");

  // initialise stiffness matrix to zero
  stiff_->Zero();

  // build predicted mid-state by last converged state and predicted target state
  EvaluateMidState();


  // ************************** (1) EXTERNAL FORCES ***************************

  // build new external forces
  fextn_->PutScalar(0.0);
  ApplyForceStiffExternal(timen_, (*dis_)(0), disn_, (*vel_)(0), fextn_, stiff_);

  // additional external forces are added (e.g. interface forces)
  fextn_->Update(1.0, *fifc_, 1.0);

  // external mid-forces F_{ext;n+1-alpha_f} ----> TR-like
  // F_{ext;n+1-alpha_f} := (1.-alphaf) * F_{ext;n+1} + alpha_f * F_{ext;n}
  fextm_->Update(1.-alphaf_, *fextn_, alphaf_, *fext_, 0.0);

  // set time integration info for plastic TSI problem
  if (HaveSemiSmoothPlasticity())
    if (plastman_->TSI())
    {
      plastman_->SetData().scale_timint_=beta_/gamma_;
      plastman_->SetData().dt_=(*dt_)[0];
    }

  // ************************** (2) INTERNAL FORCES ***************************

  fintn_->PutScalar(0.0);
  // build new internal forces and stiffness
  if (HaveNonlinearMass() == INPAR::STR::ml_none)
  {
    ApplyForceStiffInternal(timen_, (*dt_)[0], disn_, disi_, veln_, fintn_, stiff_,params);
  }
  else
  {

    if (pred_ != INPAR::STR::pred_constdis)
      dserror("Only the predictor PredictConstDisConsistVelAcc() allowed for dynamic beam3r simulations!!!");

    //If we have nonlinear inertia forces, the corresponding contributions are computed together with the internal forces
    finertn_->PutScalar(0.0);
    mass_->Zero();

    // In general the nonlinear inertia force can depend on displacements, velocities and accelerations,
    // i.e     finertn_=finertn_(disn_, veln_, accn_):
    //
    //    LIN finertn_ = [ d(finertn_)/d(disn_) + gamma_/(beta_*dt_)*d(finertn_)/d(veln_)
    //                 + 1/(beta_*dt_*dt_)*d(finertn_)/d(accn_) ]*disi_
    //
    //    LIN finertm_ = (1-alpha_m)/(beta_*dt_*dt_)[ (beta_*dt_*dt_)*d(finertn_)/d(disn_)
    //                 + (gamma_*dt_)*d(finertn_)/d(veln_) + d(finertn_)/d(accn_)]*disi_
    //
    // While the factor (1-alpha_m/(beta_*dt_*dt_) is applied later on in strtimint_genalpha.cpp the
    // factors timintfac_dis=(beta_*dt_*dt_) and timeintfac_vel=(gamma_*dt_) have directly to be applied
    // on element level before the three contributions of the linearization are summed up in mass_.

    double timintfac_dis=beta_*(*dt_)[0]*(*dt_)[0];
    double timintfac_vel=gamma_*(*dt_)[0];
    ApplyForceStiffInternalAndInertial(timen_, (*dt_)[0], timintfac_dis, timintfac_vel, disn_, disi_, veln_, accn_, fintn_, finertn_, stiff_, mass_,params,beta_,gamma_,alphaf_,alpham_);
  }

  // add forces and stiffness due to constraints
  // (for TR scale constraint matrix with the same value fintn_ is scaled with)
  Teuchos::ParameterList pcon;
  pcon.set("scaleConstrMat", (1.0-alphaf_));
  ApplyForceStiffConstraint(timen_, (*dis_)(0), disn_, fintn_, stiff_, pcon);

  // add forces and stiffness due to 0D cardiovascular coupling conditions
  Teuchos::ParameterList pwindk;
  pwindk.set("scale_timint", (1.-alphaf_));
  pwindk.set("time_step_size", (*dt_)[0]);
  ApplyForceStiffCardiovascular0D(timen_, disn_, fintn_, stiff_, pwindk);

  // add surface stress force
  ApplyForceStiffSurfstress(timen_, (*dt_)[0], disn_, fintn_, stiff_);

  // add potential forces
  ApplyForceStiffPotential(timen_, disn_, fintn_, stiff_);
  TestForceStiffPotential(timen_, disn_, step_);

  // add forces and stiffness due to spring dashpot condition
  Teuchos::ParameterList psprdash;
  psprdash.set("time_fac", gamma_/(beta_*(*dt_)[0]));
  psprdash.set("dt", (*dt_)[0]); // needed only for cursurfnormal option!!
  ApplyForceStiffSpringDashpot(stiff_,fintn_,disn_,veln_,predict,psprdash);

  // total internal mid-forces F_{int;n+1-alpha_f} ----> TR-like
  // F_{int;n+1-alpha_f} := (1.-alphaf) * F_{int;n+1} + alpha_f * F_{int;n}
  fintm_->Update(1.-alphaf_, *fintn_, alphaf_, *fint_, 0.0);

  // ************************** (3) INERTIA FORCES ***************************

  // build new inertia forces and stiffness
  if (HaveNonlinearMass() == INPAR::STR::ml_none)
  {
    // build new inertia forces and stiffness
    finertm_->PutScalar(0.0);
    // inertia forces #finertm_
    mass_->Multiply(false, *accm_, *finertm_);
  }
  else
  {
    // total inertia mid-forces F_{inert;n+1-alpha_m} ----> TR-like
    // F_{inert;n+1-alpha_m} := (1.-alpham) * F_{inert;n+1} + alpha_m * F_{inert;n}
    finertm_->Update(1.-alpham_, *finertn_, alpham_, *finert_, 0.0);
  }

  // ************************** (4) DAMPING FORCES ****************************

  // viscous forces due to Rayleigh damping
  if (damping_ == INPAR::STR::damp_rayleigh)
  {
    damp_->Multiply(false, *velm_, *fviscm_);
  }

  // ******************** Finally, put everything together ********************

  //build residual and tangent matrix for standard case
  if (HaveNonlinearMass()!=INPAR::STR::ml_rotations)
  {
    // build residual
    //    Res = M . A_{n+1-alpha_m}
    //        + C . V_{n+1-alpha_f}
    //        + F_{int;n+1-alpha_f}
    //        - F_{ext;n+1-alpha_f}
    fres_->Update(-1.0, *fextm_, 0.0);
    fres_->Update( 1.0, *fintm_, 1.0);
    fres_->Update( 1.0, *finertm_, 1.0);
    if (damping_ == INPAR::STR::damp_rayleigh)
    {
      fres_->Update(1.0, *fviscm_, 1.0);
    }

    // build tangent matrix : effective dynamic stiffness matrix
    //    K_{Teffdyn} = (1 - alpha_m)/(beta*dt^2) M
    //                + (1 - alpha_f)*y/(beta*dt) C
    //                + (1 - alpha_f) K_{T}

    stiff_->Add(*mass_, false, (1.-alpham_)/(beta_*(*dt_)[0]*(*dt_)[0]), 1.-alphaf_);
    if (damping_ == INPAR::STR::damp_rayleigh)
    {
      stiff_->Add(*damp_, false, (1.-alphaf_)*gamma_/(beta_*(*dt_)[0]), 1.0);
    }
  }
  //build residual vector and tangent matrix if a multiplicative Gen-Alpha scheme for rotations is applied
  else
  {
    BuildResStiffNLMassRot(fres_, fextn_, fintn_, finertn_, stiff_, mass_);
  }

  // apply forces and stiffness due to beam contact ----> TR-like
  // F_{c;n+1-alpha_f} := (1.-alphaf) * F_{c;n+1} + alpha_f * F_{c;n}
  ApplyForceStiffBeamContact(stiff_,fres_,disn_,predict);

  // apply forces and stiffness due to contact / meshtying ----> TR-like
  // F_{c;n+1-alpha_f} := (1.-alphaf) * F_{c;n+1} + alpha_f * F_{c;n}
  ApplyForceStiffContactMeshtying(stiff_,fres_,disn_,predict);

  // calculate RHS without local condensations (for NewtonLs)
  if (fresn_str_!=Teuchos::null)
  {
    // total internal mid-forces F_{int;n+1-alpha_f} ----> TR-like
    // F_{int;n+1-alpha_f} := (1.-alphaf) * F_{int;n+1} + alpha_f * F_{int;n}
    fresn_str_->Update(1.,*fintn_str_,0.);
    fresn_str_->Update(alphaf_, *fint_str_, 1.-alphaf_);
    fresn_str_->Update(-1.0, *fextm_, 1.0);
    fresn_str_->Update( 1.0, *finertm_, 1.0);
    if (damping_ == INPAR::STR::damp_rayleigh)
      fresn_str_->Update(1.0, *fviscm_, 1.0);
    LINALG::ApplyDirichlettoSystem(fresn_str_,zeros_,*(dbcmaps_->CondMap()));
  }

  // close stiffness matrix
  stiff_->Complete();

  // hallelujah
  return;
}

/*----------------------------------------------------------------------*/
/* Evaluate/define the residual force vector #fres_ for
 * relaxation solution with SolveRelaxationLinear */
void STR::TimIntGenAlpha::EvaluateForceStiffResidualRelax(Teuchos::ParameterList& params)
{
  // compute residual forces #fres_ and stiffness #stiff_
  EvaluateForceStiffResidual(params);

  // overwrite the residual forces #fres_ with interface load
  if (HaveNonlinearMass()!=INPAR::STR::ml_rotations)
  {
    //standard case
    fres_->Update(-1+alphaf_, *fifc_, 0.0);
  }
  else
  {
    //Remark: In the case of an multiplicative Gen-Alpha time integration scheme, all forces are evaluated at the end point n+1.
    fres_->Update(-1.0, *fifc_, 0.0);
  }

  // oh gosh
  return;
}

/*----------------------------------------------------------------------*/
/* Evaluate residual */
void STR::TimIntGenAlpha::EvaluateForceResidual()
{
  // build predicted mid-state by last converged state and predicted target state
  EvaluateMidState();

  // ************************** (1) EXTERNAL FORCES ***************************

  // build new external forces
  fextn_->PutScalar(0.0);
  ApplyForceExternal(timen_, (*dis_)(0), disn_, (*vel_)(0), fextn_);

  // additional external forces are added (e.g. interface forces)
  fextn_->Update(1.0, *fifc_, 1.0);

  // external mid-forces F_{ext;n+1-alpha_f} ----> TR-like
  // F_{ext;n+1-alpha_f} := (1.-alphaf) * F_{ext;n+1} + alpha_f * F_{ext;n}
  fextm_->Update(1.-alphaf_, *fextn_, alphaf_, *fext_, 0.0);

  // ************************** (2) INTERNAL FORCES ***************************

  fintn_->PutScalar(0.0);

  // build new internal forces and stiffness
  if (HaveNonlinearMass() == INPAR::STR::ml_none)
  {
    ApplyForceInternal(timen_, (*dt_)[0], disn_, disi_, veln_, fintn_);
  }
  else
  {
    dserror("Not implemented, yet.");
  }

  // total internal mid-forces F_{int;n+1-alpha_f} ----> TR-like
  // F_{int;n+1-alpha_f} := (1.-alphaf) * F_{int;n+1} + alpha_f * F_{int;n}
  fintm_->Update(1.-alphaf_, *fintn_, alphaf_, *fint_, 0.0);

  // ************************** (3) INERTIAL FORCES ***************************

  // build new inertia forces and stiffness
  if (HaveNonlinearMass() == INPAR::STR::ml_none)
  {
    // build new inertia forces and stiffness
    finertm_->PutScalar(0.0);
    // inertia forces #finertm_
    mass_->Multiply(false, *accm_, *finertm_);
  }
  else
  {
    dserror("Not implemented, yet.");
  }

  // ************************** (4) DAMPING FORCES ****************************

  // viscous forces due to Rayleigh damping
  if (damping_ == INPAR::STR::damp_rayleigh)
  {
    damp_->Multiply(false, *velm_, *fviscm_);
  }

  // ******************** Finally, put everything together ********************

  // build residual and tangent matrix for standard case
  if (HaveNonlinearMass() != INPAR::STR::ml_rotations)
  {
    // build residual
    //    Res = M . A_{n+1-alpha_m}
    //        + C . V_{n+1-alpha_f}
    //        + F_{int;n+1-alpha_f}
    //        - F_{ext;n+1-alpha_f}
    fres_->Update(-1.0, *fextm_, 0.0);
    fres_->Update( 1.0, *fintm_, 1.0);
    fres_->Update( 1.0, *finertm_, 1.0);
    if (damping_ == INPAR::STR::damp_rayleigh)
    {
      fres_->Update(1.0, *fviscm_, 1.0);
    }
  }
  else /* build residual vector and tangent matrix if a multiplicative Gen-Alpha
          scheme for rotations is applied */
  {
    dserror("Not implemented, yet.");
  }

  // calculate RHS without local condensations (for NewtonLs)
  if (fresn_str_ != Teuchos::null)
  {
    // total internal mid-forces F_{int;n+1-alpha_f} ----> TR-like
    // F_{int;n+1-alpha_f} := (1.-alphaf) * F_{int;n+1} + alpha_f * F_{int;n}
    fresn_str_->Update(1., *fintn_str_, 0.);
    fresn_str_->Update(alphaf_, *fint_str_, 1.-alphaf_);
    fresn_str_->Update(-1.0, *fextm_, 1.0);
    fresn_str_->Update( 1.0, *finertm_, 1.0);
    if (damping_ == INPAR::STR::damp_rayleigh)
      fresn_str_->Update(1.0, *fviscm_, 1.0);

    LINALG::ApplyDirichlettoSystem(fresn_str_, zeros_, *(dbcmaps_->CondMap()));
  }

  return;
}

/*----------------------------------------------------------------------*/
/* evaluate mid-state vectors by averaging end-point vectors */
void STR::TimIntGenAlpha::EvaluateMidState()
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
/* calculate characteristic/reference norms for forces
 * originally by lw */
double STR::TimIntGenAlpha::CalcRefNormForce()
{
  // The reference norms are used to scale the calculated iterative
  // displacement norm and/or the residual force norm. For this
  // purpose we only need the right order of magnitude, so we don't
  // mind evaluating the corresponding norms at possibly different
  // points within the timestep (end point, generalized midpoint).

  // norm of the internal forces
  double fintnorm = 0.0;
  fintnorm = STR::AUX::CalculateVectorNorm(iternorm_, fintm_);

  // norm of the external forces
  double fextnorm = 0.0;
  fextnorm = STR::AUX::CalculateVectorNorm(iternorm_, fextm_);

  // norm of the inertial forces
  double finertnorm = 0.0;
  finertnorm = STR::AUX::CalculateVectorNorm(iternorm_, finertm_);

  // norm of viscous forces
  double fviscnorm = 0.0;
  if (damping_ == INPAR::STR::damp_rayleigh)
  {
    fviscnorm = STR::AUX::CalculateVectorNorm(iternorm_, fviscm_);
  }

  // norm of reaction forces
  double freactnorm = 0.0;
  freactnorm = STR::AUX::CalculateVectorNorm(iternorm_, freact_);

  // determine worst value ==> charactersitic norm
  return std::max(fviscnorm, std::max(finertnorm, std::max(fintnorm, std::max(fextnorm, freactnorm))));
}

/*----------------------------------------------------------------------*/
/* incremental iteration update of state */
void STR::TimIntGenAlpha::UpdateIterIncrementally()
{
  // auxiliary global vectors
  Teuchos::RCP<Epetra_Vector> aux
    = LINALG::CreateVector(*DofRowMapView(), true);

  // further auxiliary variables
  const double dt = (*dt_)[0];  // step size \f$\Delta t_{n}\f$

  // new end-point displacements
  // D_{n+1}^{<k+1>} := D_{n+1}^{<k>} + IncD_{n+1}^{<k>}
  disn_->Update(1.0, *disi_, 1.0);

  // new end-point velocities
  aux->Update(1.0, *disn_, -1.0, (*dis_)[0], 0.0);
  aux->Update((beta_-gamma_)/beta_, (*vel_)[0],
              (2.0*beta_-gamma_)*dt/(2.0*beta_), (*acc_)[0],
              gamma_/(beta_*dt));
  // put only to free/non-DBC DOFs
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(aux), veln_);

  // new end-point accelerations
  aux->Update(1.0, *disn_, -1.0, (*dis_)[0], 0.0);
  aux->Update(-1.0/(beta_*dt), (*vel_)[0],
              (2.0*beta_-1.0)/(2.0*beta_), (*acc_)[0],
              1.0/(beta_*dt*dt));
  // put only to free/non-DBC DOFs
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(aux), accn_);

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/* iterative iteration update of state */
void STR::TimIntGenAlpha::UpdateIterIteratively()
{
  // new end-point displacements
  // D_{n+1}^{<k+1>} := D_{n+1}^{<k>} + IncD_{n+1}^{<k>}
  disn_->Update(1.0, *disi_, 1.0);

  // new end-point velocities
  veln_->Update(gamma_/(beta_*(*dt_)[0]), *disi_, 1.0);

  // new end-point accelerations
  accn_->Update(1.0/(beta_*(*dt_)[0]*(*dt_)[0]), *disi_, 1.0);

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/* update after time step */
void STR::TimIntGenAlpha::UpdateStepState()
{
  // velocity update for contact
  // (must be called BEFORE the following update steps)
  UpdateStepContactVUM();

  // update all old state at t_{n-1} etc
  // important for step size adaptivity
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1}, etc
  dis_->UpdateSteps(*disn_);

  // material displacements (struct ale)
  if( (dismatn_!=Teuchos::null))
    dismat_->UpdateSteps(*dismatn_);

  // new velocities at t_{n+1} -> t_n
  //    V_{n} := V_{n+1}, etc
  vel_->UpdateSteps(*veln_);
  // new accelerations at t_{n+1} -> t_n
  //    A_{n} := A_{n+1}, etc
  acc_->UpdateSteps(*accn_);

  // update new external force
  //    F_{ext;n} := F_{ext;n+1}
  fext_->Update(1.0, *fextn_, 0.0);

  // update new internal force
  //    F_{int;n} := F_{int;n+1}
  fint_->Update(1.0, *fintn_, 0.0);

  // update new inertial force
  //    F_{inert;n} := F_{inert;n+1}
  finert_->Update(1.0, *finertn_, 0.0);

  // update residual force vector for NewtonLS
  if (fresn_str_!=Teuchos::null)
    fint_str_->Update(1.,*fintn_str_,0.);

  // update surface stress
  UpdateStepSurfstress();

  // update constraints
  UpdateStepConstraint();

  // update Cardiovascular0D
  UpdateStepCardiovascular0D();

  // update constraints
  UpdateStepSpringDashpot();

  // update contact / meshtying
  UpdateStepContactMeshtying();

  // update beam contact
  UpdateStepBeamContact();

  // look out
  return;
}

/*----------------------------------------------------------------------*/
/* update after time step after output on element level*/
// update anything that needs to be updated at the element level
void STR::TimIntGenAlpha::UpdateStepElement()
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // other parameters that might be needed by the elements
  p.set("total time", timen_);
  p.set("delta time", (*dt_)[0]);
  // action for elements
  p.set("action", "calc_struct_update_istep");

  // go to elements
  discret_->ClearState();
  discret_->SetState("displacement",(*dis_)(0));

  // Set material displacement state for ale-wear formulation
  if( (dismat_!=Teuchos::null))
    discret_->SetState("material_displacement",(*dismat_)(0));

  if (!HaveNonlinearMass())
  {
    discret_->Evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  }
  else
  {
    /* In the NonlinearMass-case its possible to make an update of
     * displacements, velocities and accelerations at the end of time step
     * (currently only necessary for Kirchhoff beams). An corresponding update
     * rule has to be implemented in the element, otherwise displacements,
     * velocities and accelerations remain unchanged.
     */
    discret_->SetState("velocity",(*vel_)(0));
    discret_->SetState("acceleration",(*acc_)(0));

    Teuchos::RCP<Epetra_Vector> update_disp;
    update_disp = LINALG::CreateVector(*DofRowMapView(), true);

    Teuchos::RCP<Epetra_Vector> update_vel;
    update_vel = LINALG::CreateVector(*DofRowMapView(), true);

    Teuchos::RCP<Epetra_Vector> update_acc;
    update_acc = LINALG::CreateVector(*DofRowMapView(), true);


    discret_->Evaluate(p, Teuchos::null, Teuchos::null, update_disp, update_vel, update_acc);

    disn_->Update(1.0,*update_disp,1.0);
    (*dis_)(0)->Update(1.0,*update_disp,1.0);
    veln_->Update(1.0,*update_vel,1.0);
    (*vel_)(0)->Update(1.0,*update_vel,1.0);
    accn_->Update(1.0,*update_acc,1.0);
    (*acc_)(0)->Update(1.0,*update_acc,1.0);

  }

  discret_->ClearState();
}

/*----------------------------------------------------------------------*/
/* read and/or calculate forces for restart */
void STR::TimIntGenAlpha::ReadRestartForce()
{
  IO::DiscretizationReader reader(discret_, step_);
  reader.ReadVector(fext_, "fexternal");
  reader.ReadVector(fint_, "fint");
  reader.ReadVector(finert_, "finert");

  return;
}

/*----------------------------------------------------------------------*/
/* write internal and external forces for restart */
void STR::TimIntGenAlpha::WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output)
{
  output->WriteVector("fexternal", fext_);
  output->WriteVector("fint",fint_);
  output->WriteVector("finert",finert_);
  return;
}

/*-----------------------------------------------------------------------------*
 * Update all field vectors defined specific for this method,     sudhakar 12/13
 *  to take into account of the new nodes introduced by crack propagation
 *----------------------------------------------------------------------------*/
void STR::TimIntGenAlpha::updateMethodSpecificEpetraCrack( std::map<int,int>& oldnew )
{
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( discret_, dism_, oldnew );
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( discret_, velm_, oldnew );
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( discret_, accm_, oldnew );

  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( discret_, fint_, oldnew );
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( discret_, fintm_, oldnew );
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( discret_, fintn_, oldnew );

  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( discret_, fext_, oldnew );
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( discret_, fextm_, oldnew );
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( discret_, fextn_, oldnew );

  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( discret_, finert_, oldnew );
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( discret_, finertm_, oldnew );
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( discret_, finertn_, oldnew );
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( discret_, fviscm_, oldnew );

  // set initial external force vector
  ApplyForceExternal((*time_)[0], (*dis_)(0), disn_, (*vel_)(0), fext_);

  if (!HaveNonlinearMass())
  {
    // determine mass, damping and initial accelerations
    DetermineMassDampConsistAccel();
  }
  else
  {
    /* the case of nonlinear inertia terms works so far only for examples with
     * vanishing initial accelerations, i.e. the initial external forces and
     * initial velocities have to be chosen consistently!!!
     */
    (*acc_)(0)->PutScalar(0.0);
  }

  ApplyForceStiffExternal((*time_)[0], (*dis_)(0), disn_, (*vel_)(0), fext_, stiff_);

    // create parameter list
    Teuchos::ParameterList params;

  if (!HaveNonlinearMass())
  {
    // set initial internal force vector
    ApplyForceStiffInternal((*time_)[0], (*dt_)[0], (*dis_)(0), zeros_, (*vel_)(0), fint_, stiff_, params);
  }
  else
  {
    double timeintfac_dis=beta_*(*dt_)[0]*(*dt_)[0];
    double timeintfac_vel=gamma_*(*dt_)[0];

    // Check, if initial residuum really vanishes for acc_ = 0
    ApplyForceStiffInternalAndInertial((*time_)[0], (*dt_)[0], timeintfac_dis, timeintfac_vel, (*dis_)(0), zeros_, (*vel_)(0), (*acc_)(0), fint_, finert_, stiff_, mass_,params,beta_,gamma_,alphaf_,alpham_);

    NonlinearMassSanityCheck(fext_, (*dis_)(0), (*vel_)(0), (*acc_)(0));
  }
}

/*-----------------------------------------------------------------------------*
 * Build total residual vector and effective tangential stiffness    meier 05/14
 * matrix in case of nonlinear, rotational inertia effects
 *----------------------------------------------------------------------------*/
void STR::TimIntGenAlpha::BuildResStiffNLMassRot(
    Teuchos::RCP<Epetra_Vector> fres_,
    Teuchos::RCP<Epetra_Vector> fextn_,
    Teuchos::RCP<Epetra_Vector> fintn_,
    Teuchos::RCP<Epetra_Vector> finertn_,
    Teuchos::RCP<LINALG::SparseOperator> stiff_,
    Teuchos::RCP<LINALG::SparseOperator> mass_
    )
{
  /* build residual
   *    Res = F_{inert;n+1}
   *        + F_{int;n+1}
   *        - F_{ext;n+1}
   * Remark: In the case of an multiplicative Gen-Alpha time integration scheme,
   * all forces are evaluated at the end point n+1.
   */
  fres_->Update(-1.0, *fextn_, 0.0);
  fres_->Update( 1.0, *fintn_, 1.0);
  fres_->Update( 1.0, *finertn_, 1.0);

  /* build tangent matrix : effective dynamic stiffness matrix
   *    K_{Teffdyn} = M
   *                + K_{T}
   * Remark: So far, all time integration pre-factors (only necessary for the
   * mass matrix since internal forces are evaluated at n+1) are already
   * considered at element level (see, e.g., beam3r_evaluate.cpp). Therefore,
   * we don't have to apply them here.
   */
  stiff_->Add(*mass_, false, 1.0 , 1.0);
}

/*-----------------------------------------------------------------------------*
 * Check, if there are solely beam elements in the whole             meier 05/14
 * discretization
 *----------------------------------------------------------------------------*/
bool STR::TimIntGenAlpha::SolelyBeam3Elements(Teuchos::RCP<DRT::Discretization> actdis)
{
  bool solelybeameles=true;

  for (int i=0;i<actdis->NumMyRowElements();i++)
  {
    DRT::Element* element = actdis->lColElement(i);
    DRT::Node* node = (element->Nodes())[0];
    int numdof = actdis->NumDof(node);

    //So far we simply check, if we have at least 6 DoFs per node, which is only true for beam elements
    if (numdof < 6)
      solelybeameles=false;
  }

  return solelybeameles;
}

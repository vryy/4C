/*----------------------------------------------------------------------*/
/*!
\file strtimint_gemm.cpp
\brief Structural time integration with generalised energy-momentum method

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
#include "strtimint_gemm.H"
#include "stru_aux.H"
#include "../drt_lib/drt_locsys.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*/
/* constructor */
STR::TimIntGEMM::TimIntGEMM
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
  beta_(sdynparams.sublist("GEMM").get<double>("BETA")),
  gamma_(sdynparams.sublist("GEMM").get<double>("GAMMA")),
  alphaf_(sdynparams.sublist("GEMM").get<double>("ALPHA_F")),
  alpham_(sdynparams.sublist("GEMM").get<double>("ALPHA_M")),
  xi_(sdynparams.sublist("GEMM").get<double>("XI")),
  dism_(Teuchos::null),
  velm_(Teuchos::null),
  accm_(Teuchos::null),
  fintm_(Teuchos::null),
  fext_(Teuchos::null),
  fextm_(Teuchos::null),
  fextn_(Teuchos::null),
  finertm_(Teuchos::null),
  fviscm_(Teuchos::null)
{
  // info to user about current time integration scheme and its parametrization
  if (myrank_ == 0)
  {
    std::cout << "with generalised energy-momentum method" << std::endl
              << "   alpha_f = " << alphaf_ << std::endl
              << "   alpha_m = " << alpham_ << std::endl
              << "   xi = " << xi_ << std::endl
              << "   p_dis = " << MethodOrderOfAccuracyDis() << std::endl
              << "   p_vel = " << MethodOrderOfAccuracyVel() << std::endl
              << std::endl;
  }

  // determine mass, damping and initial accelerations
  DetermineMassDampConsistAccel();

  // create state vectors

  // mid-displacements
  dism_ = LINALG::CreateVector(*DofRowMapView(), true);
  // mid-velocities
  velm_ = LINALG::CreateVector(*DofRowMapView(), true);
  // mid-accelerations
  accm_ = LINALG::CreateVector(*DofRowMapView(), true);

  // create force vectors

  // internal force vector F_{int;m} at mid-time
  fintm_ = LINALG::CreateVector(*DofRowMapView(), true);

  // external force vector F_ext at last times
  fext_ = LINALG::CreateVector(*DofRowMapView(), true);
  // external mid-force vector F_{ext;n+1-alpha_f}
  fextm_ = LINALG::CreateVector(*DofRowMapView(), true);
  // external force vector F_{n+1} at new time
  fextn_ = LINALG::CreateVector(*DofRowMapView(), true);
  // set initial external force vector
  ApplyForceExternal((*time_)[0], (*dis_)(0), disn_, (*vel_)(0), fext_);

  // inertia mid-point force vector F_inert
  finertm_ = LINALG::CreateVector(*DofRowMapView(), true);
  // viscous mid-point force vector F_visc
  fviscm_ = LINALG::CreateVector(*DofRowMapView(), true);

  // GEMM time integrator cannot handle nonlinear inertia forces
  if (HaveNonlinearMass())
    dserror("Gemm time integrator cannot handle nonlinear inertia forces "
        "(flag: MASSLIN)");

  // have a nice day
  return;
}

/*----------------------------------------------------------------------*/
/* Consistent predictor with constant displacements
 * and consistent velocities and displacements */
void STR::TimIntGEMM::PredictConstDisConsistVelAcc()
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
void STR::TimIntGEMM::PredictConstVelConsistAcc()
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
void STR::TimIntGEMM::PredictConstAcc()
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
/* evaluate residual force and its stiffness, ie derivative
 * with respect to end-point displacements \f$D_{n+1}\f$ */
void STR::TimIntGEMM::EvaluateForceStiffResidual(Teuchos::ParameterList& params)
{
  // get info about prediction step from parameter list
  bool predict = false;
  if(params.isParameter("predict"))
    predict = params.get<bool>("predict");

  // build by last converged state and predicted target state
  // the predicted mid-state
  EvaluateMidState();

  /* add forces and stiffness due to Windkessel bcs
   * necessarily has to be done BEFORE fextm_ is built, since the Windkessel
   * manager calls an EvaluateNeumann function and thus the correct application
   * and linearization of the follower load is needed !!! (mhv 11/2013)
   */
  Teuchos::ParameterList pwindk;
  pwindk.set("scale_timint", (1.-alphaf_));
  pwindk.set("time_step_size", (*dt_)[0]);
  ApplyForceStiffWindkessel(timen_, (*dis_)(0), disn_, pwindk);

  // initialise stiffness matrix to zero
  stiff_->Zero();

  // ************************** (1) EXTERNAL FORCES ***************************

  // build new external forces
  fextn_->PutScalar(0.0);
  ApplyForceStiffExternal(timen_, (*dis_)(0), disn_, (*vel_)(0), fextn_, stiff_);

  // additional external forces are added (e.g. interface forces)
  fextn_->Update(1.0, *fifc_, 1.0);

  // external mid-forces F_{ext;n+1-alpha_f} (fextm)
  //    F_{ext;n+1-alpha_f} := (1.-alphaf) * F_{ext;n+1}
  //                         + alpha_f * F_{ext;n}
  fextm_->Update(1.-alphaf_, *fextn_, alphaf_, *fext_, 0.0);

  // ************************** (2) INTERNAL FORCES ***************************

  // initialise internal forces
  fintm_->PutScalar(0.0);

  // ordinary internal force and stiffness
  disi_->Scale(1.-alphaf_);  // CHECK THIS
  ApplyForceStiffInternalMid(timen_, (*dt_)[0], (*dis_)(0), disn_, disi_, veln_,
                             fintm_, stiff_);

  // apply forces and stiffness due to constraints
  Teuchos::ParameterList pcon; //apply empty parameterlist, no scaling necessary
  ApplyForceStiffConstraint(timen_, (*dis_)(0), disn_, fintm_, stiff_, pcon);

  // potential forces
  ApplyForceStiffPotential(timen_, dism_, fintm_, stiff_);

  // add forces and stiffness due to spring dashpot condition
  Teuchos::ParameterList psprdash;
  psprdash.set("scale_gamma", gamma_);
  psprdash.set("scale_beta", beta_);
  psprdash.set("time_step_size", (*dt_)[0]);
  ApplyForceStiffSpringDashpot(stiff_,fintm_,disn_,veln_,predict,psprdash);

  // ************************** (3) INERTIAL FORCES ***************************

  // inertial forces #finertm_
  mass_->Multiply(false, *accm_, *finertm_);

  // ************************** (4) DAMPING FORCES ****************************

  // viscous forces due Rayleigh damping
  if (damping_ == INPAR::STR::damp_rayleigh)
  {
    damp_->Multiply(false, *velm_, *fviscm_);
  }

  // ******************** Finally, put everything together ********************

  // build residual
  //    Res = M . A_{n+1-alpha_m}
  //        + C . V_{n+1-alpha_f}
  //        + F_{int;m}
  //        - F_{ext;n+1-alpha_f}
  fres_->Update(-1.0, *fextm_, 0.0);
  fres_->Update(1.0, *fintm_, 1.0);
  if (damping_ == INPAR::STR::damp_rayleigh)
  {
    fres_->Update(1.0, *fviscm_, 1.0);
  }
  fres_->Update(1.0, *finertm_, 1.0);

  // build tangent matrix : effective dynamic stiffness matrix
  //    K_{Teffdyn} = (1 - alpha_m)/(beta*dt^2) M
  //                + (1 - alpha_f)*y/(beta*dt) C
  //                + K_{T;m}
  stiff_->Add(*mass_, false, (1.-alpham_)/(beta_*(*dt_)[0]*(*dt_)[0]), 1.0);
  if (damping_ == INPAR::STR::damp_rayleigh)
  {
    stiff_->Add(*damp_, false, (1.-alphaf_)*gamma_/(beta_*(*dt_)[0]), 1.0);
  }

  // apply forces and stiffness due to contact / meshtying
  // Note that we ALWAYS use a TR-like approach to compute the interface
  // forces. This means we never explicitly compute fc at the generalized
  // mid-point n+1-alphaf, but use a linear combination of the old end-
  // point n and the new end-point n+1 instead:
  // F_{c;n+1-alpha_f} := (1-alphaf) * F_{c;n+1} +  alpha_f * F_{c;n}
  ApplyForceStiffContactMeshtying(stiff_,fres_,disn_,predict);

  // close stiffness matrix
  stiff_->Complete();

  // hallelujah
  return;
}

/*----------------------------------------------------------------------*/
/* Evaluate/define the residual force vector #fres_ for
 * relaxation solution with SolveRelaxationLinear */
void STR::TimIntGEMM::EvaluateForceStiffResidualRelax(Teuchos::ParameterList& params)
{
  // compute residual forces #fres_ and stiffness #stiff_
  EvaluateForceStiffResidual(params);

  // overwrite the residual forces #fres_ with interface load
  fres_->Update(-(1.0-alphaf_), *fifc_, 0.0);
}

/*----------------------------------------------------------------------*/
/* Evaluate residual */
void STR::TimIntGEMM::EvaluateForceResidual()
{
  // build predicted mid-state by last converged state and predicted target state
  EvaluateMidState();

  // ************************** (1) EXTERNAL FORCES ***************************

  // build new external forces
  fextn_->PutScalar(0.0);
  ApplyForceExternal(timen_, (*dis_)(0), disn_, (*vel_)(0), fextn_);

  // additional external forces are added (e.g. interface forces)
  fextn_->Update(1.0, *fifc_, 1.0);

  // external mid-forces F_{ext;n+1-alpha_f} (fextm)
  //    F_{ext;n+1-alpha_f} := (1.-alphaf) * F_{ext;n+1}
  //                         + alpha_f * F_{ext;n}
  fextm_->Update(1.-alphaf_, *fextn_, alphaf_, *fext_, 0.0);

  // ************************** (2) INTERNAL FORCES ***************************

  // initialise internal forces
  fintm_->PutScalar(0.0);

  // ordinary internal force and stiffness
  disi_->Scale(1.-alphaf_);  // CHECK THIS
  ApplyForceInternalMid(timen_, (*dt_)[0], (*dis_)(0), disn_, disi_, veln_,
      fintm_);

  // ************************** (3) INERTIAL FORCES ***************************

  // inertial forces #finertm_
  mass_->Multiply(false, *accm_, *finertm_);

  // ************************** (4) DAMPING FORCES ****************************

  // viscous forces due Rayleigh damping
  if (damping_ == INPAR::STR::damp_rayleigh)
  {
    damp_->Multiply(false, *velm_, *fviscm_);
  }

  // ******************** Finally, put everything together ********************

  // build residual
  //    Res = M . A_{n+1-alpha_m}
  //        + C . V_{n+1-alpha_f}
  //        + F_{int;m}
  //        - F_{ext;n+1-alpha_f}
  fres_->Update(-1.0, *fextm_, 0.0);
  fres_->Update(1.0, *fintm_, 1.0);
  if (damping_ == INPAR::STR::damp_rayleigh)
  {
    fres_->Update(1.0, *fviscm_, 1.0);
  }
  fres_->Update(1.0, *finertm_, 1.0);

  return;
}

/*----------------------------------------------------------------------*/
/* evaluate mid-state vectors by averaging end-point vectors */
void STR::TimIntGEMM::EvaluateMidState()
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
double STR::TimIntGEMM::CalcRefNormForce()
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
void STR::TimIntGEMM::UpdateIterIncrementally()
{
  // auxiliary global vectors
  Teuchos::RCP<Epetra_Vector> aux
      = LINALG::CreateVector(*DofRowMapView(), false);

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
  // put new velocities only on non-DBC/free DOFs
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(aux), veln_);


  // new end-point accelerations
  aux->Update(1.0, *disn_, -1.0, (*dis_)[0], 0.0);
  aux->Update(-1.0/(beta_*dt), (*vel_)[0],
              (2.0*beta_-1.0)/(2.0*beta_), (*acc_)[0],
              1.0/(beta_*dt*dt));
  // put new accelerations only on free DOFs
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(aux), accn_);

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/* iterative iteration update of state */
void STR::TimIntGEMM::UpdateIterIteratively()
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
void STR::TimIntGEMM::UpdateStepState()
{
  // velocity update for contact
  // (must be called BEFORE the following update steps)
  UpdateStepContactVUM();

  // update all old state at t_{n-1} etc
  // important for step size adaptivity
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1}, etc
  dis_->UpdateSteps(*disn_);
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
  // nothing to be done

  // update surface stress
  UpdateStepSurfstress();

  // update constraints
  UpdateStepConstraint();

  // update Windkessel
  UpdateStepWindkessel();

  // update constraints
  UpdateStepSpringDashpot();

  // update contact  /meshtying
  UpdateStepContactMeshtying();

  // look out
  return;
}

/*----------------------------------------------------------------------*/
/* update after time step after output on element level*/
// update anything that needs to be updated at the element level
void STR::TimIntGEMM::UpdateStepElement()
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // other parameters that might be needed by the elements
  p.set("total time", timen_);
  //p.set("delta time", (*dt_)[0]);
  // action for elements
  //p.set("alpha f", alphaf_);
  p.set("action", "calc_struct_update_istep");
  // go to elements
  discret_->Evaluate(p, Teuchos::null, Teuchos::null,
                     Teuchos::null, Teuchos::null, Teuchos::null);
}

/*----------------------------------------------------------------------*/
/* evaluate ordinary internal force, its stiffness at mid-state */
void STR::TimIntGEMM::ApplyForceStiffInternalMid
(
  const double time,
  const double dt,
  const Teuchos::RCP<Epetra_Vector> dis,  // displacement state at t_n
  const Teuchos::RCP<Epetra_Vector> disn,  // displacement state at t_{n+1}
  const Teuchos::RCP<Epetra_Vector> disi,  // residual displacements
  const Teuchos::RCP<Epetra_Vector> vel,  // velocity state
  Teuchos::RCP<Epetra_Vector> fint,  // internal force
  Teuchos::RCP<LINALG::SparseOperator> stiff  // stiffness matrix
)
{
  // *********** time measurement ***********
  double dtcpu = timer_->WallTime();
  // *********** time measurement ***********

  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // action for elements
  const std::string action = "calc_struct_nlnstiff_gemm";
  p.set("action", action);
  // other parameters that might be needed by the elements
  p.set("total time", time);
  p.set("delta time", dt);
  p.set("alpha f", alphaf_);
  p.set("xi", xi_);
  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("old displacement", dis);
  discret_->SetState("displacement", disn);
  discret_->SetState("residual displacement", disi);
  if (damping_ == INPAR::STR::damp_material) discret_->SetState("velocity", vel);
  //fintn_->PutScalar(0.0);  // initialise internal force vector
  discret_->Evaluate(p, stiff, Teuchos::null,
                     fint, Teuchos::null, Teuchos::null);
  discret_->ClearState();

  // *********** time measurement ***********
  dtele_ = timer_->WallTime() - dtcpu;
  // *********** time measurement ***********

  return;
}

/*----------------------------------------------------------------------*/
/* evaluate ordinary internal force at mid-state */
void STR::TimIntGEMM::ApplyForceInternalMid
(
  const double time,
  const double dt,
  const Teuchos::RCP<Epetra_Vector> dis,
  const Teuchos::RCP<Epetra_Vector> disn,
  const Teuchos::RCP<Epetra_Vector> disi,
  const Teuchos::RCP<Epetra_Vector> vel,
  Teuchos::RCP<Epetra_Vector> fint
)
{
  // *********** time measurement ***********
  double dtcpu = timer_->WallTime();
  // *********** time measurement ***********

  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // action for elements
  const std::string action = "calc_struct_nlnstiff_gemm";
  p.set("action", action);
  // other parameters that might be needed by the elements
  p.set("total time", time);
  p.set("delta time", dt);
  p.set("alpha f", alphaf_);
  p.set("xi", xi_);
  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("old displacement", dis);
  discret_->SetState("displacement", disn);
  discret_->SetState("residual displacement", disi);
  if (damping_ == INPAR::STR::damp_material) discret_->SetState("velocity", vel);
  //fintn_->PutScalar(0.0);  // initialise internal force vector
  discret_->Evaluate(p, Teuchos::null, Teuchos::null,
                     fint, Teuchos::null, Teuchos::null);
  discret_->ClearState();

  // *********** time measurement ***********
  dtele_ = timer_->WallTime() - dtcpu;
  // *********** time measurement ***********

  return;
}

/*----------------------------------------------------------------------*/
/* read restart forces */
void STR::TimIntGEMM::ReadRestartForce()
{
  IO::DiscretizationReader reader(discret_, step_);
  reader.ReadVector(fext_, "fexternal");
  return;
}

/*----------------------------------------------------------------------*/
/* write external forces for restart */
void STR::TimIntGEMM::WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output)
{
  output->WriteVector("fexternal", fext_);
  return;
}

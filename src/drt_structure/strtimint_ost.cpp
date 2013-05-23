/*----------------------------------------------------------------------*/
/*!
\file strtimint_ost.cpp
\brief Structural time integration with one-step-theta

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint_ost.H"
#include "stru_aux.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_pstream.H"

/*----------------------------------------------------------------------*/
void STR::TimIntOneStepTheta::VerifyCoeff()
{
  // beta
  if ( (theta_ <= 0.0) or (theta_ > 1.0) )
    dserror("theta out of range (0.0,1.0]");

  // done
  return;
}

/*======================================================================*/
/* constructor */
STR::TimIntOneStepTheta::TimIntOneStepTheta
(
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
    ioparams,
    sdynparams,
    xparams,
    actdis,
    solver,
    contactsolver,
    output
  ),
  theta_(sdynparams.sublist("ONESTEPTHETA").get<double>("THETA")),
  dist_(Teuchos::null),
  velt_(Teuchos::null),
  acct_(Teuchos::null),
  fint_(Teuchos::null),
  fintn_(Teuchos::null),
  fext_(Teuchos::null),
  fextn_(Teuchos::null),
  finertt_(Teuchos::null),
  fvisct_(Teuchos::null)
{
  // info to user
  if (myrank_ == 0)
  {
    VerifyCoeff();
    IO::cout << "with one-step-theta" << IO::endl
             << "   theta = " << theta_ << IO::endl
             << IO::endl;
  }

  // determine mass, damping and initial accelerations
  DetermineMassDampConsistAccel();

  // create state vectors

  // mid-displacements
  dist_ = LINALG::CreateVector(*dofrowmap_, true);
  // mid-velocities
  velt_ = LINALG::CreateVector(*dofrowmap_, true);
  // mid-accelerations
  acct_ = LINALG::CreateVector(*dofrowmap_, true);

  // create force vectors

  // internal force vector F_{int;n} at last time
  fint_ = LINALG::CreateVector(*dofrowmap_, true);
  // internal force vector F_{int;n+1} at new time
  fintn_ = LINALG::CreateVector(*dofrowmap_, true);
  // set initial internal force vector
  ApplyForceStiffInternal((*time_)[0], (*dt_)[0], (*dis_)(0), zeros_, (*vel_)(0),
                          fint_, stiff_);

  // external force vector F_ext at last times
  fext_ = LINALG::CreateVector(*dofrowmap_, true);
  // external force vector F_{n+1} at new time
  fextn_ = LINALG::CreateVector(*dofrowmap_, true);
  // set initial external force vector
  ApplyForceExternal((*time_)[0], (*dis_)(0), disn_, (*vel_)(0), fext_, stiff_);

  // inertial mid-point force vector F_inert
  finertt_ = LINALG::CreateVector(*dofrowmap_, true);
  // viscous mid-point force vector F_visc
  fvisct_ = LINALG::CreateVector(*dofrowmap_, true);

  // have a nice day
  return;
}

/*----------------------------------------------------------------------*/
/* Consistent predictor with constant displacements
 * and consistent velocities and displacements */
void STR::TimIntOneStepTheta::PredictConstDisConsistVelAcc()
{
  // time step size
  const double dt = (*dt_)[0];

  // constant predictor : displacement in domain
  disn_->Update(1.0, *(*dis_)(0), 0.0);

  // new end-point velocities
  veln_->Update(1.0/(theta_*dt), *disn_,
                -1.0/(theta_*dt), *(*dis_)(0),
                0.0);
  veln_->Update(-(1.0-theta_)/theta_, *(*vel_)(0),
                1.0);

  // new end-point accelerations
  accn_->Update(1.0/(theta_*theta_*dt*dt), *disn_,
                -1.0/(theta_*theta_*dt*dt), *(*dis_)(0),
                0.0);
  accn_->Update(-1.0/(theta_*theta_*dt), *(*vel_)(0),
                -(1.0-theta_)/theta_, *(*acc_)(0),
                1.0);

  // watch out
  return;
}

/*----------------------------------------------------------------------*/
/* Consistent predictor with constant velocities,
 * extrapolated displacements and consistent accelerations */
void STR::TimIntOneStepTheta::PredictConstVelConsistAcc()
{
  // time step size
  const double dt = (*dt_)[0];

  // extrapolated displacements based upon constant velocities
  // d_{n+1} = d_{n} + dt * v_{n}
  disn_->Update(1.0, (*dis_)[0], dt, (*vel_)[0], 0.0);

  // new end-point velocities
  veln_->Update(1.0/(theta_*dt), *disn_,
                -1.0/(theta_*dt), *(*dis_)(0),
                0.0);
  veln_->Update(-(1.0-theta_)/theta_, *(*vel_)(0),
                1.0);

  // new end-point accelerations
  accn_->Update(1.0/(theta_*theta_*dt*dt), *disn_,
                -1.0/(theta_*theta_*dt*dt), *(*dis_)(0),
                0.0);
  accn_->Update(-1.0/(theta_*theta_*dt), *(*vel_)(0),
                -(1.0-theta_)/theta_, *(*acc_)(0),
                1.0);
  // That's it!
  return;
}

/*----------------------------------------------------------------------*/
/* Consistent predictor with constant accelerations
 * and extrapolated velocities and displacements */
void STR::TimIntOneStepTheta::PredictConstAcc()
{
  // time step size
  const double dt = (*dt_)[0];

  // extrapolated displacements based upon constant accelerations
  // d_{n+1} = d_{n} + dt * v_{n} + dt^2 / 2 * a_{n}
  disn_->Update(1.0, (*dis_)[0], dt, (*vel_)[0], 0.0);
  disn_->Update(dt * dt / 2., (*acc_)[0], 1.0);

  // new end-point velocities
  veln_->Update(1.0/(theta_*dt), *disn_,
                -1.0/(theta_*dt), *(*dis_)(0),
                0.0);
  veln_->Update(-(1.0-theta_)/theta_, *(*vel_)(0),
                1.0);

  // new end-point accelerations
  accn_->Update(1.0/(theta_*theta_*dt*dt), *disn_,
                -1.0/(theta_*theta_*dt*dt), *(*dis_)(0),
                0.0);
  accn_->Update(-1.0/(theta_*theta_*dt), *(*vel_)(0),
                -(1.0-theta_)/theta_, *(*acc_)(0),
                1.0);

  // That's it!
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate residual force and its stiffness, ie derivative
 * with respect to end-point displacements \f$D_{n+1}\f$ */
void STR::TimIntOneStepTheta::EvaluateForceStiffResidual(bool predict)
{
  // initialise stiffness matrix to zero
  stiff_->Zero();
  if(damping_ == INPAR::STR::damp_material)
    damp_->Zero();

  // theta-interpolate state vectors
  EvaluateMidState();

  // build new external forces
  fextn_->PutScalar(0.0);
  ApplyForceExternal(timen_, (*dis_)(0), disn_, (*vel_)(0), fextn_, stiff_);

  // additional external forces are added (e.g. interface forces)
  fextn_->Update(1.0, *fifc_, 1.0);

  // initialise internal forces
  fintn_->PutScalar(0.0);

  // ordinary internal force and stiffness
  ApplyForceStiffInternal(timen_, (*dt_)[0], disn_, disi_, veln_, fintn_, stiff_,damp_);

  // apply forces and stiffness due to constraints
  ParameterList pcon;
  //constraint matrix has to be scaled with the same value fintn_ is scaled with
  pcon.set("scaleConstrMat",theta_);
  ApplyForceStiffConstraint(timen_, (*dis_)(0), disn_, fintn_, stiff_, pcon);

  // potential forces
  ApplyForceStiffPotential(timen_, disn_, fintn_, stiff_);

  // inertial forces #finertt_
  mass_->Multiply(false, *acct_, *finertt_);

  // viscous forces due Rayleigh damping
  if (damping_ == INPAR::STR::damp_rayleigh)
  {
    damp_->Multiply(false, *velt_, *fvisct_);
  }

  // build residual  Res = M . A_{n+theta}
  //                     + C . V_{n+theta}
  //                     + F_{int;n+theta}
  //                     - F_{ext;n+theta}
  fres_->Update(-theta_, *fextn_, -(1.0-theta_), *fext_, 0.0);
  fres_->Update(theta_, *fintn_, (1.0-theta_), *fint_, 1.0);
  if (damping_ == INPAR::STR::damp_rayleigh)
  {
    fres_->Update(1.0, *fvisct_, 1.0);
  }
  fres_->Update(1.0, *finertt_, 1.0);

  //cout << STR::AUX::CalculateVectorNorm(vectornorm_l2, fextn_) << endl;

  // build tangent matrix : effective dynamic stiffness matrix
  //    K_{Teffdyn} = 1/(theta*dt^2) M
  //                + 1/dt C
  //                + theta K_{T}
  stiff_->Add(*mass_, false, 1.0/(theta_*(*dt_)[0]*(*dt_)[0]), theta_);
  if (damping_ != INPAR::STR::damp_none)
  {
    if(damping_ == INPAR::STR::damp_material)
      damp_->Complete();
    stiff_->Add(*damp_, false, 1.0/(*dt_)[0], 1.0);
  }

  // apply forces and stiffness due to contact / meshtying
  // Note that we ALWAYS use a TR-like approach to compute the interface
  // forces. This means we never explicitly compute fc at the generalized
  // mid-point n+theta, but use a linear combination of the old end-
  // point n and the new end-point n+1 instead:
  // F_{c;n+theta} := theta * F_{c;n+1} +  (1-theta) * F_{c;n}
  ApplyForceStiffContactMeshtying(stiff_,fres_,disn_,predict);

  // apply forces and stiffness due to beam contact
  ApplyForceStiffBeamContact(stiff_,fres_,disn_,predict);

  // close stiffness matrix
  stiff_->Complete();

  // hallelujah
  return;
}

/*----------------------------------------------------------------------*/
/* Evaluate/define the residual force vector #fres_ for
 * relaxation solution with SolveRelaxationLinear */
void STR::TimIntOneStepTheta::EvaluateForceStiffResidualRelax()
{
  // compute residual forces #fres_ and stiffness #stiff_
  EvaluateForceStiffResidual();

  // overwrite the residual forces #fres_ with interface load
  fres_->Update(-theta_, *fifc_, 0.0);
}

/*----------------------------------------------------------------------*/
/* evaluate theta-state vectors by averaging end-point vectors */
void STR::TimIntOneStepTheta::EvaluateMidState()
{
  // mid-displacements D_{n+1-alpha_f} (dism)
  //    D_{n+theta} := theta * D_{n+1} + (1-theta) * D_{n}
  dist_->Update(theta_, *disn_, 1.0-theta_, *(*dis_)(0), 0.0);

  // mid-velocities V_{n+1-alpha_f} (velm)
  //    V_{n+theta} := theta * V_{n+1} + (1-theta) * V_{n}
  velt_->Update(theta_, *veln_, 1.0-theta_, *(*vel_)(0), 0.0);

  // mid-accelerations A_{n+1-alpha_m} (accm)
  //    A_{n+theta} := theta * A_{n+1} + (1-theta) * A_{n}
  acct_->Update(theta_, *accn_, 1.0-theta_, *(*acc_)(0), 0.0);

  // jump
  return;
}

/*----------------------------------------------------------------------*/
/* calculate characteristic/reference norms for displacements
 * originally by lw */
double STR::TimIntOneStepTheta::CalcRefNormDisplacement()
{
  // The reference norms are used to scale the calculated iterative
  // displacement norm and/or the residual force norm. For this
  // purpose we only need the right order of magnitude, so we don't
  // mind evaluating the corresponding norms at possibly different
  // points within the timestep (end point, generalized midpoint).

  double charnormdis = 0.0;
  if (pressure_ != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> disp = pressure_->ExtractOtherVector((*dis_)(0));
    charnormdis = STR::AUX::CalculateVectorNorm(iternorm_, disp);
  }
  else
    charnormdis = STR::AUX::CalculateVectorNorm(iternorm_, (*dis_)(0));

  // rise your hat
  return charnormdis;
}

/*----------------------------------------------------------------------*/
/* calculate characteristic/reference norms for forces
 * originally by lw */
double STR::TimIntOneStepTheta::CalcRefNormForce()
{
  // The reference norms are used to scale the calculated iterative
  // displacement norm and/or the residual force norm. For this
  // purpose we only need the right order of magnitude, so we don't
  // mind evaluating the corresponding norms at possibly different
  // points within the timestep (end point, generalized midpoint).

  // norm of the internal forces
  double fintnorm = 0.0;
  fintnorm = STR::AUX::CalculateVectorNorm(iternorm_, fintn_);

  // norm of the external forces
  double fextnorm = 0.0;
  fextnorm = STR::AUX::CalculateVectorNorm(iternorm_, fextn_);

  // norm of the inertial forces
  double finertnorm = 0.0;
  finertnorm = STR::AUX::CalculateVectorNorm(iternorm_, finertt_);

  // norm of viscous forces
  double fviscnorm = 0.0;
  if (damping_ == INPAR::STR::damp_rayleigh)
  {
    fviscnorm = STR::AUX::CalculateVectorNorm(iternorm_, fvisct_);
  }

  // norm of reaction forces
  double freactnorm = 0.0;
  freactnorm = STR::AUX::CalculateVectorNorm(iternorm_, freact_);

  // return char norm
  return max(fviscnorm, max(finertnorm, max(fintnorm, max(fextnorm, freactnorm))));
}

/*----------------------------------------------------------------------*/
/* incremental iteration update of state */
void STR::TimIntOneStepTheta::UpdateIterIncrementally()
{
  // Auxiliary vector holding new velocities and accelerations
  // by extrapolation/scheme on __all__ DOFs. This includes
  // the Dirichlet DOFs as well. Thus we need to protect those
  // DOFs of overwriting; they already hold the
  // correctly 'predicted', final values.
  Teuchos::RCP<Epetra_Vector> aux
      = LINALG::CreateVector(*dofrowmap_, false);

  // new end-point displacements
  // D_{n+1}^{<k+1>} := D_{n+1}^{<k>} + IncD_{n+1}^{<k>}
  disn_->Update(1.0, *disi_, 1.0);

  // new end-point velocities
  aux->Update(1.0/(theta_*(*dt_)[0]), *disn_,
               -1.0/(theta_*(*dt_)[0]), *(*dis_)(0),
               0.0);
  aux->Update(-(1.0-theta_)/theta_, *(*vel_)(0), 1.0);
  // put only to free/non-DBC DOFs
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(aux), veln_);

  // new end-point accelerations
  aux->Update(1.0/(theta_*theta_*(*dt_)[0]*(*dt_)[0]), *disn_,
              -1.0/(theta_*theta_*(*dt_)[0]*(*dt_)[0]), *(*dis_)(0),
              0.0);
  aux->Update(-1.0/(theta_*theta_*(*dt_)[0]), *(*vel_)(0),
              -(1.0-theta_)/theta_, *(*acc_)(0),
              1.0);
  // put only to free/non-DBC DOFs
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(aux), accn_);

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/* iterative iteration update of state */
void STR::TimIntOneStepTheta::UpdateIterIteratively()
{
  // new end-point displacements
  // D_{n+1}^{<k+1>} := D_{n+1}^{<k>} + IncD_{n+1}^{<k>}
  disn_->Update(1.0, *disi_, 1.0);

  // new end-point velocities
  veln_->Update(1.0/(theta_*(*dt_)[0]), *disi_, 1.0);

  // new end-point accelerations
  accn_->Update(1.0/((*dt_)[0]*(*dt_)[0]*theta_*theta_), *disi_, 1.0);

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/* update after time step */
void STR::TimIntOneStepTheta::UpdateStepState()
{
  // update state
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1}
  dis_->UpdateSteps(*disn_);
  // new velocities at t_{n+1} -> t_n
  //    V_{n} := V_{n+1}
  vel_->UpdateSteps(*veln_);
  // new accelerations at t_{n+1} -> t_n
  //    A_{n} := A_{n+1}
  acc_->UpdateSteps(*accn_);

  // update new external force
  //    F_{ext;n} := F_{ext;n+1}
  fext_->Update(1.0, *fextn_, 0.0);

  // update new internal force
  //    F_{int;n} := F_{int;n+1}
  fint_->Update(1.0, *fintn_, 0.0);

  // update surface stress
  UpdateStepSurfstress();

  // update constraints
  UpdateStepConstraint();

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
void STR::TimIntOneStepTheta::UpdateStepElement()
{
  // create the parameters for the discretization
  ParameterList p;
  // other parameters that might be needed by the elements
  p.set("total time", timen_);
  p.set("delta time", (*dt_)[0]);
  //p.set("alpha f", theta_);
  // action for elements
  p.set("action", "calc_struct_update_istep");
  // go to elements
  discret_->ClearState();
  discret_->SetState("displacement",(*dis_)(0));
  discret_->Evaluate(p, Teuchos::null, Teuchos::null,
                     Teuchos::null, Teuchos::null, Teuchos::null);
  discret_->ClearState();
}

/*----------------------------------------------------------------------*/
/* read restart forces */
void STR::TimIntOneStepTheta::ReadRestartForce()
{
  IO::DiscretizationReader reader(discret_, step_);
  reader.ReadVector(fext_, "fexternal");
  reader.ReadVector(fint_, "fint");

  return;
}

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* write internal and external forces for restart */
void STR::TimIntOneStepTheta::WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output)
{
  output->WriteVector("fexternal",fext_);
  output->WriteVector("fint",fint_);
  return;
}

/*----------------------------------------------------------------------*/
/*!
\file thrtimint_statics.cpp
\brief Statics analysis

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
#include "thrtimint_statics.H"

/*======================================================================*/
/* constructor */
THR::TimIntStatics::TimIntStatics
(
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& tdynparams,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization> actdis,
  Teuchos::RCP<LINALG::Solver> solver,
  Teuchos::RCP<IO::DiscretizationWriter> output
)
: TimIntImpl
  (
    ioparams,
    tdynparams,
    xparams,
    actdis,
    solver,
    output
  ),
  fint_(Teuchos::null),
  fintn_(Teuchos::null),
  fext_(Teuchos::null),
  fextn_(Teuchos::null)
{
  // info to user
  if (myrank_ == 0)
  {
    std::cout << "with statics" << std::endl;
  }

  // create force vectors

  // internal force vector F_{int;n} at last time
  fint_ = LINALG::CreateVector(*dofrowmap_, true);
  // internal force vector F_{int;n+1} at new time
  fintn_ = LINALG::CreateVector(*dofrowmap_, true);
  // set initial internal force vector
  ApplyForceTangInternal((*time_)[0], (*dt_)[0], (*temp_)(0), zeros_,
                          fint_, tang_);

  // external force vector F_ext at last times
  fext_ = LINALG::CreateVector(*dofrowmap_, true);
  // external force vector F_{n+1} at new time
  fextn_ = LINALG::CreateVector(*dofrowmap_, true);
  // set initial external force vector
  ApplyForceExternal((*time_)[0], (*temp_)(0), fext_);

  // have a nice day
  return;
}

/*----------------------------------------------------------------------*/
/* Consistent predictor with constant temperatures
 * and consistent temperature rates and temperatures */
void THR::TimIntStatics::PredictConstTempConsistRate()
{
  // constant predictor : temperature in domain
  tempn_->Update(1.0, *(*temp_)(0), 0.0);

  // new end-point temperature rates, these stay zero in static calculation
  raten_->PutScalar(0.0);

  // watch out
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate residual force and its tangent, ie derivative
 * with respect to end-point temperatures \f$T_{n+1}\f$ */
void THR::TimIntStatics::EvaluateRhsTangResidual()
{
  // build new external forces
  fextn_->PutScalar(0.0);
  ApplyForceExternal(timen_, (*temp_)(0), fextn_);

  // initialise internal forces
  fintn_->PutScalar(0.0);

  // initialise tangent matrix to zero
  tang_->Zero();

  // ordinary internal force and tangent
  ApplyForceTangInternal(timen_, (*dt_)[0], tempn_, tempi_, fintn_, tang_);

  // potential forces
  //ApplyForceTangPotential(tempn_, fintn_, tang_);

  // build residual  Res = F_{int;n+1}
  //                     - F_{ext;n+1}
  fres_->Update(-1.0, *fextn_, 0.0);
  fres_->Update(1.0, *fintn_, 1.0);
//  cout << *fintn_ << endl;

  // build tangent matrix : effective dynamic tangent matrix
  //    K_{Teffdyn} = K_{T}
  // i.e. do nothing here
  tang_->Complete();  // close tangent matrix

  // hallelujah
  return;
}

/*----------------------------------------------------------------------*/
/* calculate characteristic/reference norms for temperatures
 * originally by lw */
double THR::TimIntStatics::CalcRefNormTemperature()
{
  // The reference norms are used to scale the calculated iterative
  // temperature norm and/or the residual force norm. For this
  // purpose we only need the right order of magnitude, so we don't
  // mind evaluating the corresponding norms at possibly different
  // points within the timestep (end point, generalized midpoint).

  double charnormtemp = 0.0;
    charnormtemp = THR::AUX::CalculateVectorNorm(iternorm_, (*temp_)(0));

  // rise your hat
  return charnormtemp;
}

/*----------------------------------------------------------------------*/
/* calculate characteristic/reference norms for forces
 * originally by lw */
double THR::TimIntStatics::CalcRefNormForce()
{
  // The reference norms are used to scale the calculated iterative
  // temperature norm and/or the residual force norm. For this
  // purpose we only need the right order of magnitude, so we don't
  // mind evaluating the corresponding norms at possibly different
  // points within the timestep (end point, generalized midpoint).

  // norm of the internal forces
  double fintnorm = 0.0;
  fintnorm = THR::AUX::CalculateVectorNorm(iternorm_, fintn_);

  // norm of the external forces
  double fextnorm = 0.0;
  fextnorm = THR::AUX::CalculateVectorNorm(iternorm_, fextn_);

  // norm of reaction forces
  double freactnorm = 0.0;
  freactnorm = THR::AUX::CalculateVectorNorm(iternorm_, freact_);

  // return char norm
  return max(fintnorm, max(fextnorm, freactnorm));
}

/*----------------------------------------------------------------------*/
/* incremental iteration update of state */
void THR::TimIntStatics::UpdateIterIncrementally()
{
  // new end-point temperatures
  // T_{n+1}^{<k+1>} := T_{n+1}^{<k>} + IncT_{n+1}^{<k>}
  tempn_->Update(1.0, *tempi_, 1.0);

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/* iterative iteration update of state */
void THR::TimIntStatics::UpdateIterIteratively()
{
  // new end-point temperatures
  // T_{n+1}^{<k+1>} := T_{n+1}^{<k>} + IncT_{n+1}^{<k>}
  tempn_->Update(1.0, *tempi_, 1.0);

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/* update after time step */
void THR::TimIntStatics::UpdateStepState()
{
  // update state
  // new temperatures at t_{n+1} -> t_n
  //    T_{n} := T_{n+1}
  temp_->UpdateSteps(*tempn_);
  // new temperature rates at t_{n+1} -> t_n
  //    T'_{n} := T'_{n+1}
  rate_->UpdateSteps(*raten_);  // this simply copies zero vectors

  // update new external force
  //    F_{ext;n} := F_{ext;n+1}
  fext_->Update(1.0, *fextn_, 0.0);

  // update new internal force
  //    F_{int;n} := F_{int;n+1}
  fint_->Update(1.0, *fintn_, 0.0);

  // update anything that needs to be updated at the element level
  {
    // create the parameters for the discretization
    ParameterList p;
    // other parameters that might be needed by the elements
    p.set("total time", timen_);
    p.set("delta time", (*dt_)[0]);
    //   //p.set("alpha f", theta_);
    // action for elements
    p.set("action", "thermo_update_istep");
    // go to elements
    discret_->Evaluate(p, null, null, null, null, null);
  }

  // look out
  return;
}

/*----------------------------------------------------------------------*/
/* read restart forces */
void THR::TimIntStatics::ReadRestartForce()
{
  IO::DiscretizationReader reader(discret_, step_);
  // set 'initial' external force
  reader.ReadVector(fext_, "fexternal");
  // set 'initial' internal force vector
  // Set dt to 0, since we do not propagate in time.
  ApplyForceInternal((*time_)[0], 0.0, (*temp_)(0), zeros_, fint_);

  return;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET

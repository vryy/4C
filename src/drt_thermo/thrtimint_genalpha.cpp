 /*----------------------------------------------------------------------*/
/*!
\file thrtimint_genalpha.cpp
\brief Thermal time integration with generalized-alpha

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 | headers                                                   dano 10/09 |
 *----------------------------------------------------------------------*/
#include "thrtimint_genalpha.H"

/*----------------------------------------------------------------------*
 |                                                           dano 10/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntGenAlpha::VerifyCoeff()
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
  if ( (alpham_ < 0.0) or (alpham_ >= 1.0) )
    dserror("alpha_m out of range [0.0,1.0)");
  else
    std::cout << "   alpha_m = " << alpham_ << std::endl;

  // mid-averaging type
  // In principle, there exits two common possibilities, namely TR-like and IMR-like,
  // where TR-like means trapezoidal rule and IMR-like means implicit mid-point rule.
  // We use to maintain implementations of both variants, but due to its significantly
  // higher complexity, the IMR-like version has been deleted (popp 02/2013). The nice
  // thing about TR-like mid-averaging is that all element (and thus also material) calls
  // are exclusively(!) carried out at the end-point t_{n+1} of each time interval, but
  // never explicitly at some generalized midpoint, such as t_{n+1-\alpha_f}.
  if (midavg_ != INPAR::THR::midavg_trlike)
    dserror("mid-averaging of internal forces only implemented TR-like");
  else
    std::cout << "   midavg = " << INPAR::THR::MidAverageString(midavg_)<<std::endl;

  // done
  return;
}

/*----------------------------------------------------------------------*
 | constructor                                               dano 10/09 |
 *----------------------------------------------------------------------*/
THR::TimIntGenAlpha::TimIntGenAlpha
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

  midavg_(DRT::INPUT::IntegralValue<INPAR::THR::MidAverageEnum>(tdynparams.sublist("GENALPHA"),"GENAVG")),
  /* iterupditer_(false), */
  beta_(tdynparams.sublist("GENALPHA").get<double>("BETA")),
  gamma_(tdynparams.sublist("GENALPHA").get<double>("GAMMA")),
  alphaf_(tdynparams.sublist("GENALPHA").get<double>("ALPHA_F")),
  alpham_(tdynparams.sublist("GENALPHA").get<double>("ALPHA_M")),
  tempm_(Teuchos::null),
  ratem_(Teuchos::null),
  fint_(Teuchos::null),
  fintm_(Teuchos::null),
  fintn_(Teuchos::null),
  fext_(Teuchos::null),
  fextm_(Teuchos::null),
  fextn_(Teuchos::null)
{
  // so far framework available but can not be used
  dserror("Not yet ready to use! In the meanwhile use STATICS or OneStepTheta!");

  // info to user
  if (myrank_ == 0)
  {
    std::cout << "with generalized-alpha" << std::endl;
    VerifyCoeff();
//    std::cout << "   p_temp = " << MethodOrderOfAccuracyTemp() << std::endl
//              << "   p_rate = " << MethodOrderOfAccuracyRate() << std::endl

//    std::cout << "   p_temp = " << MethodOrderOfAccuracy() << std::endl

//              << std::endl;
  }

  // determine capacity and initial temperature rates
  DetermineCapaConsistTempRate();

  // create state vectors

  // mid-temperatures
  tempm_ = LINALG::CreateVector(*dofrowmap_, true);
  // mid-temperature rates
  ratem_ = LINALG::CreateVector(*dofrowmap_, true);

  // create force vectors

  // internal force vector F_{int;n} at last time
  fint_ = LINALG::CreateVector(*dofrowmap_, true);
  // internal mid-force vector F_{int;n+1-alpha_f}
  fintm_ = LINALG::CreateVector(*dofrowmap_, true);
  // internal force vector F_{int;n+1} at new time
  fintn_ = LINALG::CreateVector(*dofrowmap_, true);
  // set initial internal force vector
  ApplyForceTangInternal((*time_)[0], (*dt_)[0], (*temp_)(0), zeros_,fint_, tang_);

  // external force vector F_ext at last times
  fext_ = LINALG::CreateVector(*dofrowmap_, true);
  // external mid-force vector F_{ext;n+1-alpha_f}
  fextm_ = LINALG::CreateVector(*dofrowmap_, true);
  // external force vector F_{n+1} at new time
  fextn_ = LINALG::CreateVector(*dofrowmap_, true);
  // set initial external force vector
  ApplyForceExternal((*time_)[0], (*temp_)(0), fext_);

  // have a nice day
  return;
}

/*----------------------------------------------------------------------*
 |  Consistent predictor with constant temperatures          dano 10/09 |
 |  and consistent temperature rates and temperatures                   |
 *----------------------------------------------------------------------*/
void THR::TimIntGenAlpha::PredictConstTempConsistRate()
{
  // constant predictor : temperature in domain
  tempn_->Update(1.0, *(*temp_)(0), 0.0);

  // consistent temperature rates
  //! Ueberpruefen, da acc_ noch drin!!!!!!!!!!!!!!!!! 01.10.2009

  raten_->Update(1.0, *tempn_, -1.0, *(*temp_)(0), 0.0);
//  raten_->Update((beta_-gamma_)/beta_, *(*rate_)(0),
//                (2.*beta_-gamma_)*(*dt_)[0]/(2.*beta_), *(*acc_)(0),
//                gamma_/(beta_*(*dt_)[0]));
  raten_->Update((beta_-gamma_)/beta_, *(*rate_)(0),
                 gamma_/(beta_*(*dt_)[0]));

  // watch out
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate residual force and its tangent, ie derivative   dano 10/09 |
 |  with respect to end-point temperatures \f$T_{n+1}\f$                |
 *----------------------------------------------------------------------*/
void THR::TimIntGenAlpha::EvaluateRhsTangResidual()
{
  // build by last converged state and predicted target state
  // the predicted mid-state
  EvaluateMidState();

  // build new external forces
  fextn_->PutScalar(0.0);
//  ApplyForceExternal(timen_, (*temp_)(0), (*rate_)(0), fextn_);
  ApplyForceExternal(timen_, (*temp_)(0), fextn_);

  // external mid-forces F_{ext;n+1-alpha_f} (fextm)
  //    F_{ext;n+1-alpha_f} := (1.-alphaf) * F_{ext;n+1}
  //                         + alpha_f * F_{ext;n}
  fextm_->Update(1.-alphaf_, *fextn_, alphaf_, *fext_, 0.0);

  // initialise internal forces
  fintn_->PutScalar(0.0);

  // initialise tangent matrix to zero
  tang_->Zero();

  // ordinary internal force and tangent

////    ApplyForceTangInternal(timen_, (*dt_)[0], tempn_, tempi_,  raten_,
////                            fintn_, tang_);

  // 02.10.09 nach dem Vorbild von thrtimint_ost.cpp
    ApplyForceTangInternal(timen_, (*dt_)[0], tempn_, tempi_, fintn_, tang_);

  // build residual
  //    Res = C . R_{n+1-alpha_f}
  //        + F_{int;m}
  //        - F_{ext;n+1-alpha_f}
  fres_->Update(-1.0, *fextm_, 0.0);
  fres_->Update((1.-alphaf_), *fintn_, alphaf_, *fint_, 1.0);

  // build tangent matrix : effective dynamic tangent matrix
  //    K_{Teffdyn} = (1 - alpha_m)/(beta*dt^2) M
  //                + (1 - alpha_f)*y/(beta*dt) C
  //                + (1 - alpha_f) K_{T}

  // vgl. thrtimint_ost.cpp: Linie auskommentiert!!! 02.10.09
//  tang_->Add(*capa_, false, (1.-alpham_)/(beta_*(*dt_)[0]*(*dt_)[0]), 1.-alphaf_);

  tang_->Complete();  // close tangent matrix

  // hallelujah
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate mid-state vectors by averaging end-point        dano 10/09 |
 |  vectors                                                             |
 *----------------------------------------------------------------------*/
void THR::TimIntGenAlpha::EvaluateMidState()
{
  // mid-temperatures T_{n+1-alpha_f} (tempm)
  //    T_{n+1-alpha_f} := (1.-alphaf) * T_{n+1} + alpha_f * T_{n}
  tempm_->Update(1.-alphaf_, *tempn_, alphaf_, (*temp_)[0], 0.0);

  // mid-temperature rates R_{n+1-alpha_f} (ratem)
  //    R_{n+1-alpha_f} := (1.-alphaf) * R_{n+1} + alpha_f * R_{n}
  ratem_->Update(1.-alphaf_, *raten_, alphaf_, (*rate_)[0], 0.0);

  // jump
  return;
}

/*----------------------------------------------------------------------*
 |  calculate characteristic/reference norms for             dano 10/09 |
 |  temperatures originally by lw                                       |
 *----------------------------------------------------------------------*/
double THR::TimIntGenAlpha::CalcRefNormTemperature()
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

/*----------------------------------------------------------------------*
 |  calculate characteristic/reference norms for forces      dano 10/09 |
 |  originally by lw                                                    |
 *----------------------------------------------------------------------*/
double THR::TimIntGenAlpha::CalcRefNormForce()
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
  fextnorm = THR::AUX::CalculateVectorNorm(iternorm_, fextm_);

  // norm of reaction forces
  double freactnorm = 0.0;
  freactnorm = THR::AUX::CalculateVectorNorm(iternorm_, freact_);

  // determine worst value ==> charactersitic norm
  return std::max(fintnorm, std::max(fextnorm, freactnorm));
}

/*----------------------------------------------------------------------*
 |  update after time step                                   dano 10/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntGenAlpha::UpdateIterIncrementally()
{
  // auxiliar global vectors
  Teuchos::RCP<Epetra_Vector> aux
    = LINALG::CreateVector(*dofrowmap_, true);
  // further auxiliar variables
  const double dt = (*dt_)[0];  // step size \f$\Delta t_{n}\f$

  // new end-point temperatures
  // T_{n+1}^{<k+1>} := T_{n+1}^{<k>} + IncT_{n+1}^{<k>}
  tempn_->Update(1.0, *tempi_, 1.0);

  // new end-point temperature rates
  aux->Update(1.0, *tempn_, -1.0, (*temp_)[0], 0.0);
//  aux->Update((beta_-gamma_)/beta_, (*rate_)[0],
//              (2.0*beta_-gamma_)*dt/(2.0*beta_), (*acc_)[0],
//              gamma_/(beta_*dt));
  aux->Update((beta_-gamma_)/beta_, (*rate_)[0],
              gamma_/(beta_*dt));

  // put only to free/non-DBC DOFs
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(aux), raten_);

  // bye
  return;
}

/*----------------------------------------------------------------------*
 |  iterative iteration update of state                      dano 10/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntGenAlpha::UpdateIterIteratively()
{
  // new end-point temperatures
  // T_{n+1}^{<k+1>} := T_{n+1}^{<k>} + IncT_{n+1}^{<k>}
  tempn_->Update(1.0, *tempi_, 1.0);

  // new end-point temperature rates
  raten_->Update(gamma_/(beta_*(*dt_)[0]), *tempi_, 1.0);

  // bye
  return;
}

/*----------------------------------------------------------------------*
 |  update after time step                                   dano 10/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntGenAlpha::UpdateStepState()
{
  // update all old state at t_{n-1} etc
  // important for step size adaptivity
  // new temperatures at t_{n+1} -> t_n
  //    T_{n} := T_{n+1}, etc
  temp_->UpdateSteps(*tempn_);
  // new temperature rates at t_{n+1} -> t_n
  //    R_{n} := R_{n+1}, etc
  rate_->UpdateSteps(*raten_);

  // update new external force
  //    F_{ext;n} := F_{ext;n+1}
  fext_->Update(1.0, *fextn_, 0.0);

  // update new internal force
  //    F_{int;n} := F_{int;n+1}
  fint_->Update(1.0, *fintn_, 0.0);

  // update anything that needs to be updated at the element level
  {
    // create the parameters for the discretization
    Teuchos::ParameterList p;
    // other parameters that might be needed by the elements
    p.set("total time", timen_);
    p.set("delta time", (*dt_)[0]);

//  // action for elements
//  p.set("action", "calc_therm_update_istep");

    // go to elements
    discret_->Evaluate(p, Teuchos::null, Teuchos::null,
                       Teuchos::null, Teuchos::null, Teuchos::null);
  }

  // look out
  return;
}

/*----------------------------------------------------------------------*
 |  read restart forces                                      dano 10/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntGenAlpha::ReadRestartForce()
{
  IO::DiscretizationReader reader(discret_, step_);
  // external force
  reader.ReadVector(fext_, "fexternal");
  // determine internal force
  fint_->PutScalar(0.0);
  // Set dt to 0, since we do not propagate in time.
  // No time integration on material level
  ApplyForceInternal((*time_)[0], 0.0, (*temp_)(0), zeros_, fint_);

  // bye
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the internal force and the tangent              dano 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntGenAlpha::ApplyForceTangInternal(
  const double time,  //!< evaluation time
  const double dt,  //!< step size
  const Teuchos::RCP<Epetra_Vector> temp,  //!< temperature state
  const Teuchos::RCP<Epetra_Vector> tempi,  //!< residual temperatures
  Teuchos::RCP<Epetra_Vector> fint,  //!< internal force
  Teuchos::RCP<LINALG::SparseMatrix> tang  //!< tangent matrix
)
{
  //! create the parameters for the discretization
  Teuchos::ParameterList p;
  //! set parameters
  // p.set();
  //! call the base function
  TimInt::ApplyForceTangInternal(p,time,dt,temp,tempi,fint,tang);
  //! finish
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the internal force                              dano 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntGenAlpha::ApplyForceInternal(
  const double time,  //!< evaluation time
  const double dt,  //!< step size
  const Teuchos::RCP<Epetra_Vector> temp,  //!< temperature state
  const Teuchos::RCP<Epetra_Vector> tempi,  //!< incremental temperatures
  Teuchos::RCP<Epetra_Vector> fint  //!< internal force
)
{
  //! create the parameters for the discretization
  Teuchos::ParameterList p;
  //! set parameters
//  p.set();
  //! call the base function
  TimInt::ApplyForceInternal(p,time,dt,temp,tempi,fint);
  //! finish
  return;
}


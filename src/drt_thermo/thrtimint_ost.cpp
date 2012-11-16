/*----------------------------------------------------------------------*/
/*!
\file thrtimint_ost.cpp
\brief Thermal time integration with one-step-theta

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237

</pre>
*/


/*----------------------------------------------------------------------*
 |  headers                                                  dano 08/09 |
 *----------------------------------------------------------------------*/
#include "thrtimint_ost.H"

/*----------------------------------------------------------------------*
 |  constructor                                              dano 08/09 |
 *----------------------------------------------------------------------*/
THR::TimIntOneStepTheta::TimIntOneStepTheta(
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& tdynparams,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization> actdis,
  Teuchos::RCP<LINALG::Solver> solver,
  Teuchos::RCP<IO::DiscretizationWriter> output
  )
: TimIntImpl(
    ioparams,
    tdynparams,
    xparams,
    actdis,
    solver,
    output
    ),
  theta_(tdynparams.sublist("ONESTEPTHETA").get<double>("THETA")),
  tempt_(Teuchos::null),
  ratet_(Teuchos::null),
  fint_(Teuchos::null),
  fintn_(Teuchos::null),
  fcap_(Teuchos::null),
  fcapn_(Teuchos::null),
  fext_(Teuchos::null),
  fextn_(Teuchos::null)
{
  // info to user
  if (myrank_ == 0)
  {
    std::cout << "with one-step-theta" << std::endl
              << "   theta = " << theta_ << std::endl
              << std::endl;
  }

  // determine capacity
  DetermineCapaConsistTempRate();

  // create state vectors
  // mid-temperatures
  tempt_ = LINALG::CreateVector(*dofrowmap_, true);
  // mid-temperature rates
  ratet_ = LINALG::CreateVector(*dofrowmap_, true);

  // create force vectors
  // internal force vector F_{int;n} at last time
  fint_ = LINALG::CreateVector(*dofrowmap_, true);
  // internal force vector F_{int;n+1} at new time
  fintn_ = LINALG::CreateVector(*dofrowmap_, true);
  // stored force vector F_{transient;n} at last time
  fcap_ = LINALG::CreateVector(*dofrowmap_, true);
  // stored force vector F_{transient;n+1} at new time
  fcapn_ = LINALG::CreateVector(*dofrowmap_, true);
  // set initial internal force vector
  ApplyForceTangInternal((*time_)[0], (*dt_)[0], (*temp_)(0),zeros_,fcap_,
                         fint_, tang_);
  // external force vector F_ext at last times
  fext_ = LINALG::CreateVector(*dofrowmap_, true);
  // external force vector F_{n+1} at new time
  fextn_ = LINALG::CreateVector(*dofrowmap_, true);

  // set initial external force vector
  ApplyForceExternal((*time_)[0], (*temp_)(0), fext_);

  // set initial external force vector of convective heat transfer boundary
  // conditions
  ApplyForceExternalConv((*time_)[0], (*temp_)(0), (*temp_)(0), fext_, tang_);

  // have a nice day
  return;
}

/*----------------------------------------------------------------------*
 |  Consistent predictor with constant temperatures          dano 08/09 |
 |  and consistent temperature rates and temperatures                   |
 *----------------------------------------------------------------------*/
void THR::TimIntOneStepTheta::PredictConstTempConsistRate()
{
  // time step size
  const double dt = (*dt_)[0];

  // constant predictor : temperature in domain
  tempn_->Update(1.0, *(*temp_)(0), 0.0);

  // new end-point temperature rates
  raten_->Update(1.0/(theta_*dt), *tempn_,
                -1.0/(theta_*dt), *(*temp_)(0),
                0.0);
  raten_->Update(-(1.0-theta_)/theta_, *(*rate_)(0),
                1.0);

  // watch out
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate residual force and its tangent, ie derivative   dano 08/09 |
 |  with respect to end-point temperatures \f$T_{n+1}\f$                |
 *----------------------------------------------------------------------*/
void THR::TimIntOneStepTheta::EvaluateRhsTangResidual()
{
  // theta-interpolate state vectors
  EvaluateMidState();

  // build new external forces
  fextn_->PutScalar(0.0);

  // initialise tangent matrix to zero
  tang_->Zero();

  // set initial external force vector of convective heat transfer boundary
  // conditions

  // if the boundary condition shall be dependent on the current temperature
  // solution T_n+1 --> linearisation must be uncommented
  // --> use tempn_

  // if the old temperature T_n  is sufficient --> no linearisation needed!
  // --> use (*temp_)(0)
  ApplyForceExternalConv(timen_, (*temp_)(0), tempn_, fextn_, tang_);

  ApplyForceExternal(timen_, (*temp_)(0), fextn_);

  // interface forces to external forces
  fextn_->Update(1.0, *fifc_, 1.0);

  // initialise internal forces
  fintn_->PutScalar(0.0);
  fcapn_->PutScalar(0.0);

  // ordinary internal force and tangent
  ApplyForceTangInternal(timen_, (*dt_)[0], tempn_, tempi_, fcapn_, fintn_,
                         tang_);

  // build residual  Res = C . R_{n+theta}
  //                     + F_{int;n+theta}
  //                     - F_{ext;n+theta}
  // with R_{n+theta}     = C . ( T_{n+1} - T_n ) / dt
  //      F_{int;n+theta} = theta * F_{int;n+1} + (1 - theta) * F_{int;n}
  //      F_{ext;n+theta} = - theta * F_{ext;n+1} - (1 - theta) * F_{ext;n}
  fres_->Update(1.0, *fcapn_, -1.0, *fcap_, 0.0);
  // here the time derivative is introduced needed for fcap depending on T'!
  fres_->Scale(1.0/(*dt_)[0]);
  fres_->Update(theta_, *fintn_, (1.0-theta_), *fint_, 1.0);
  // here is the negative sign for the external forces (heatfluxes)
  fres_->Update(-theta_, *fextn_, -(1.0-theta_), *fext_, 1.0);

  // apply modifications due to thermal contact
  if (thermcontman_!= Teuchos::null)
  {
    fres_->Scale(-1);
    thermcontman_->ApplyThermoContact(tang_,fres_,tempn_,(*dt_)[0]);
    fres_->Scale(-1);
  }

  tang_->Complete();  // close tangent matrix

  // hallelujah
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate theta-state vectors by averaging                dano 08/09 |
 |  end-point vector                                                    |
 *----------------------------------------------------------------------*/
void THR::TimIntOneStepTheta::EvaluateMidState()
{
  // mid-temperatures T_{n+1-alpha_f} (tempm)
  //    T_{n+theta} := theta * T_{n+1} + (1-theta) * T_{n}
  tempt_->Update(theta_, *tempn_, 1.0-theta_, *(*temp_)(0), 0.0);

  // mid-temperature rates R_{n+1-alpha_f} (ratem)
  //    R_{n+theta} := theta * R_{n+1} + (1-theta) * R_{n}
  ratet_->Update(theta_, *raten_, 1.0-theta_, *(*rate_)(0), 0.0);

  // jump
  return;
}

/*----------------------------------------------------------------------*
 |  calculate characteristic/reference norms for             dano 08/09 |
 |  temperatures originally by lw                                       |
 *----------------------------------------------------------------------*/
double THR::TimIntOneStepTheta::CalcRefNormTemperature()
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
 |  calculate characteristic/reference norms for forces      dano 08/09 |
 |  originally by lw                                                    |
 *----------------------------------------------------------------------*/
double THR::TimIntOneStepTheta::CalcRefNormForce()
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

  // norm of stored forces
  double fcapnorm = 0.0;
  fcapnorm = THR::AUX::CalculateVectorNorm(iternorm_, fcap_);

  // return char norm
  return max(fcapnorm, max(fintnorm, max(fextnorm, freactnorm)));
}

/*----------------------------------------------------------------------*
 |  incremental iteration update of state                    dano 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntOneStepTheta::UpdateIterIncrementally()
{
  // Auxiliar vector holding new temperature rates
  // by extrapolation/scheme on __all__ DOFs. This includes
  // the Dirichlet DOFs as well. Thus we need to protect those
  // DOFs of overwriting; they already hold the
  // correctly 'predicted', final values.
  Teuchos::RCP<Epetra_Vector> aux
      = LINALG::CreateVector(*dofrowmap_, false);

  // new end-point temperatures
  // T_{n+1}^{<k+1>} := T_{n+1}^{<k>} + IncT_{n+1}^{<k>}
  tempn_->Update(1.0, *tempi_, 1.0);

  // new end-point temperature rates
  // aux = 1/(theta . Dt) T_n+1^k+1 - 1/(theta . Dt) T_n+1^k - (1-theta)/theta T'_n+1^k
  aux->Update(1.0/(theta_*(*dt_)[0]), *tempn_,
               -1.0/(theta_*(*dt_)[0]), *(*temp_)(0),
               0.0);
  aux->Update(-(1.0-theta_)/theta_, *(*rate_)(0), 1.0);
  // put only to free/non-DBC DOFs
  dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(aux), raten_);

  // bye
  return;
}

/*----------------------------------------------------------------------*
 |  iterative iteration update of state                      dano 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntOneStepTheta::UpdateIterIteratively()
{
  // new end-point temperatures
  // T_{n+1}^{<k+1>} := T_{n+1}^{<k>} + IncT_{n+1}^{<k>}
  tempn_->Update(1.0, *tempi_, 1.0);

  // new end-point temperature rates
  raten_->Update(1.0/(theta_*(*dt_)[0]), *tempi_, 1.0);

  // bye
  return;
}

/*----------------------------------------------------------------------*
 |  update after time step                                   dano 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntOneStepTheta::UpdateStepState()
{
  // update state
  // new temperatures at t_{n+1} -> t_n
  //    T_{n} := T_{n+1}
  temp_->UpdateSteps(*tempn_);
  // new temperature rates at t_{n+1} -> t_n
  //    R_{n} := R_{n+1}
  rate_->UpdateSteps(*raten_);

  // update new external force
  //    F_{ext;n} := F_{ext;n+1}
  fext_->Update(1.0, *fextn_, 0.0);

  // update new internal force
  //    F_{int;n} := F_{int;n+1}
  fint_->Update(1.0, *fintn_, 0.0);

  // update new stored force
  //    F_{transient;n} := F_{transient;n+1}
  fcap_->Update(1.0, *fcapn_, 0.0);

  // update anything that needs to be updated at the element level
  {
    // create the parameters for the discretization
    Teuchos::ParameterList p;
    // other parameters that might be needed by the elements
    p.set("total time", timen_);
    p.set("delta time", (*dt_)[0]);
    //p.set("alpha f", theta_);
    // action for elements
    p.set("action", "calc_thermo_update_istep");
    // go to elements
    discret_->Evaluate(p, Teuchos::null, Teuchos::null,
                       Teuchos::null, Teuchos::null, Teuchos::null);
  }

  // look out
  return;
}

/*----------------------------------------------------------------------*
 |  read restart forces                                      dano 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntOneStepTheta::ReadRestartForce()
{
  IO::DiscretizationReader reader(discret_, step_);
  // set 'initial' external force
  reader.ReadVector(fext_, "fexternal");
  tang_->Zero();
  fint_->PutScalar(0.0);
  fcap_->PutScalar(0.0);
  // set 'initial' internal force vector and capacity
  // Set dt to 0, since we do not propagate in time.
  ApplyForceTangInternal(timen_, 0.0, tempn_, zeros_, fcap_, fint_,
                         tang_);
  tang_->Reset();
  
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the internal force and the tangent              dano 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntOneStepTheta::ApplyForceTangInternal(
  const double time,  //!< evaluation time
  const double dt,  //!< step size
  const Teuchos::RCP<Epetra_Vector> temp,  //!< temperature state
  const Teuchos::RCP<Epetra_Vector> tempi,  //!< residual temperatures
  Teuchos::RCP<Epetra_Vector> fcap,  //!< stored force (like inertial force in STR)
  Teuchos::RCP<Epetra_Vector> fint,  //!< internal force
  Teuchos::RCP<LINALG::SparseMatrix> tang  //!< tangent matrix
  )
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // set parameters
  p.set<double>("theta", theta_);
  // call the base function
  TimInt::ApplyForceTangInternal(p,time,dt,temp,tempi,fcap,fint,tang);
  // finish
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the internal force                              dano 08/09 |
 *----------------------------------------------------------------------*/
void THR::TimIntOneStepTheta::ApplyForceInternal(
  const double time,  //!< evaluation time
  const double dt,  //!< step size
  const Teuchos::RCP<Epetra_Vector> temp,  //!< temperature state
  const Teuchos::RCP<Epetra_Vector> tempi,  //!< incremental temperatures
  Teuchos::RCP<Epetra_Vector> fint  //!< internal force
  )
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // set parameters
  p.set("theta", theta_);
  // call the base function
  TimInt::ApplyForceInternal(p,time,dt,temp,tempi,fint);
  // finish
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the convective boundary condition               dano 12/10 |
 *----------------------------------------------------------------------*/
void THR::TimIntOneStepTheta::ApplyForceExternalConv(
  const double time,  //!< evaluation time
  const Teuchos::RCP<Epetra_Vector> tempn,  //!< old temperature state T_n
  const Teuchos::RCP<Epetra_Vector> temp,  //!< temperature state T_n+1
  Teuchos::RCP<Epetra_Vector> fext,  //!< external force
  Teuchos::RCP<LINALG::SparseMatrix> tang  //!< tangent matrix
  )
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // set parameters
  p.set<double>("theta", theta_);
  // call the base function
  TimInt::ApplyForceExternalConv(p,time,tempn,temp,fext,tang);
  // finish
  return;
}

/*----------------------------------------------------------------------*/

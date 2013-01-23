/*----------------------------------------------------------------------*/
/*!
\file strtimint_statics.cpp
\brief Statics analysis

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint_statics.H"
#include "stru_aux.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../linalg/linalg_utils.H"

/*======================================================================*/
/* constructor */
STR::TimIntStatics::TimIntStatics
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
  fint_(Teuchos::null),
  fintn_(Teuchos::null),
  fext_(Teuchos::null),
  fextn_(Teuchos::null)
{

  INPAR::STR::PreStress pstype = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdynparams,"PRESTRESS");
  INPAR::STR::DynamicType dyntype = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams,"DYNAMICTYP");

  if ( pstype != INPAR::STR::prestress_none && dyntype != INPAR::STR::dyna_statics ) 
  {
    dserror("Paranoia Error: PRESTRESS is only allowed in combinations with DYNAMICTYPE Statics!!");
  }

  // info to user
  if (myrank_ == 0)
  {
    // check if we are in prestressin mode
    if (pstype == INPAR::STR::prestress_mulf) IO::cout << "with static MULF prestress" << IO::endl;
    else if (pstype == INPAR::STR::prestress_id) IO::cout << "with static INVERSE DESIGN prestress" << IO::endl;
    else IO::cout << "with statics" << IO::endl;
  }
  

  // check if predictor is admissible for statics
  if ( (pred_ == INPAR::STR::pred_constvel) or (pred_ == INPAR::STR::pred_constacc))
    dserror("Predictor does not make sense for statics, only predictor with constant displacements.");

  // create force vectors

  // internal force vector F_{int;n} at last time
  fint_ = LINALG::CreateVector(*dofrowmap_, true);
  // internal force vector F_{int;n+1} at new time
  fintn_ = LINALG::CreateVector(*dofrowmap_, true);

  // external force vector F_ext at last times
  fext_ = LINALG::CreateVector(*dofrowmap_, true);
  // external force vector F_{n+1} at new time
  fextn_ = LINALG::CreateVector(*dofrowmap_, true);

  // initial resize of multi-step quantities
  if (xparams.isParameter("MSTEPEVRY") )
  {
    if (xparams.get<int>("MSTEPEVRY"))
      ResizeMStep();
  }

  // have a nice day
  return;
}

/*----------------------------------------------------------------------*/
/* Resizing of multi-step quantities */
void STR::TimIntStatics::ResizeMStep()
{
  // resize time and stepsize acc. to mstepevery_
  time_->Resize(-(stepmax_-1), 0, (*time_)[0]);
  dt_->Resize(-(stepmax_-1), 0, (*dt_)[0]);

  // resize state vectors according to mstepevery_
  dis_->Resize(-(stepmax_-1), 0, dofrowmap_, true);
  vel_->Resize(-(stepmax_-1), 0, dofrowmap_, true);
  acc_->Resize(-(stepmax_-1), 0, dofrowmap_, true);

}

/*----------------------------------------------------------------------*/
/* Consistent predictor with constant displacements
 * and consistent velocities and displacements */
void STR::TimIntStatics::PredictConstDisConsistVelAcc()
{
  // constant predictor : displacement in domain
  disn_->Update(1.0, *(*dis_)(0), 0.0);

  // new end-point velocities, these stay zero in static calculation
  veln_->PutScalar(0.0);

  // new end-point accelerations, these stay zero in static calculation
  accn_->PutScalar(0.0);

  // watch out
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate residual force and its stiffness, ie derivative
 * with respect to end-point displacements \f$D_{n+1}\f$ */
void STR::TimIntStatics::EvaluateForceStiffResidual(bool predict)
{
  // initialize stiffness matrix to zero
  stiff_->Zero();

  // build new external forces
  fextn_->PutScalar(0.0);
  ApplyForceExternal(timen_, (*dis_)(0), disn_, (*vel_)(0), fextn_, stiff_);

  // additional external forces are added (e.g. interface forces)
  fextn_->Update(1.0, *fifc_, 1.0);

  // initialize internal forces
  fintn_->PutScalar(0.0);

  // ordinary internal force and stiffness
  ApplyForceStiffInternal(timen_, (*dt_)[0], disn_, disi_, veln_, fintn_, stiff_);

  // apply forces and stiffness due to constraints
  Teuchos::ParameterList pcon; //apply empty parameterlist, no scaling necessary
  ApplyForceStiffConstraint(timen_, (*dis_)(0), disn_, fintn_, stiff_, pcon);

  // potential forces
  ApplyForceStiffPotential(timen_, disn_, fintn_, stiff_);
  TestForceStiffPotential(timen_, disn_, step_);

  // apply forces and stiffness due to embedding tissue condition
  ApplyForceStiffEmbedTissue(stiff_,fintn_,disn_,predict);

  // build residual  Res = F_{int;n+1}
  //                     - F_{ext;n+1}
  fres_->Update(-1.0, *fextn_, 0.0);
  fres_->Update(1.0, *fintn_, 1.0);
//  cout << *fintn_ << endl;

  // build tangent matrix : effective dynamic stiffness matrix
  //    K_{Teffdyn} = K_{T}
  // i.e. do nothing here

  // apply forces and stiffness due to contact / meshtying
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
void STR::TimIntStatics::EvaluateForceStiffResidualRelax()
{
  // compute residual forces #fres_ and stiffness #stiff_
  EvaluateForceStiffResidual();

  // overwrite the residual forces #fres_ with interface load
  fres_->Update(-1.0, *fifc_, 0.0);
}

/*----------------------------------------------------------------------*/
/* calculate characteristic/reference norms for displacements
 * originally by lw */
double STR::TimIntStatics::CalcRefNormDisplacement()
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
double STR::TimIntStatics::CalcRefNormForce()
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

  // norm of reaction forces
  double freactnorm = 0.0;
  freactnorm = STR::AUX::CalculateVectorNorm(iternorm_, freact_);

  // return char norm
  return max(fintnorm, max(fextnorm, freactnorm));
}

/*----------------------------------------------------------------------*/
/* incremental iteration update of state */
void STR::TimIntStatics::UpdateIterIncrementally()
{
  // new end-point displacements
  // D_{n+1}^{<k+1>} := D_{n+1}^{<k>} + IncD_{n+1}^{<k>}
  disn_->Update(1.0, *disi_, 1.0);

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/* iterative iteration update of state */
void STR::TimIntStatics::UpdateIterIteratively()
{
  // new end-point displacements
  // D_{n+1}^{<k+1>} := D_{n+1}^{<k>} + IncD_{n+1}^{<k>}
  disn_->Update(1.0, *disi_, 1.0);

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/* update after time step */
void STR::TimIntStatics::UpdateStepState()
{
  // update state
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1}
  dis_->UpdateSteps(*disn_);
  // new velocities at t_{n+1} -> t_n
  //    V_{n} := V_{n+1}
  vel_->UpdateSteps(*veln_);  // this simply copies zero vectors
  // new accelerations at t_{n+1} -> t_n
  //    A_{n} := A_{n+1}
  acc_->UpdateSteps(*accn_);  // this simply copies zero vectors

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
void STR::TimIntStatics::UpdateStepElement()
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
  discret_->Evaluate(p, null, null, null, null,Teuchos::null);
  discret_->ClearState();
}

/*----------------------------------------------------------------------*/
/* read restart forces */
void STR::TimIntStatics::ReadRestartForce()
{
  IO::DiscretizationReader reader(discret_, step_);
  // set 'initial' external force
  reader.ReadVector(fext_, "fexternal");
  // set 'initial' internal force vector
  // Set dt to 0, since we do not propagate in time.
  ApplyForceInternal((*time_)[0], 0.0, (*dis_)(0), zeros_, (*vel_)(0), fint_);

  ParameterList pcon; //no scaling necessary
  ApplyForceStiffConstraint((*time_)[0], (*dis_)(0), (*dis_)(0), fint_, stiff_, pcon);
  return;
}

/*----------------------------------------------------------------------*/


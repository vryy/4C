/*----------------------------------------------------------------------*/
/*!
\file strtimint_statics.cpp
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
#include "strtimint_statics.H"

/*======================================================================*/
/* constructor */
STR::TimIntStatics::TimIntStatics
(
  const Teuchos::ParameterList& ioparams,
  const Teuchos::ParameterList& sdynparams,
  const Teuchos::ParameterList& xparams,
  Teuchos::RCP<DRT::Discretization> actdis,
  Teuchos::RCP<LINALG::Solver> solver,
  Teuchos::RCP<IO::DiscretizationWriter> output
)
: TimIntImpl
  (
    ioparams,
    sdynparams,
    xparams,
    actdis,
    solver,
    output
  ),
  fint_(Teuchos::null),
  fintn_(Teuchos::null),
  fext_(Teuchos::null),
  fextn_(Teuchos::null),
  frobin_(Teuchos::null)
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
  ApplyForceStiffInternal((*time_)[0], (*dt_)[0], (*dis_)(0), zeros_, (*vel_)(0), 
                          fint_, stiff_);

  // external force vector F_ext at last times
  fext_ = LINALG::CreateVector(*dofrowmap_, true);
  // external force vector F_{n+1} at new time
  fextn_ = LINALG::CreateVector(*dofrowmap_, true);
  // set initial external force vector
  ApplyForceExternal((*time_)[0], (*dis_)(0), (*vel_)(0), fext_);

  // external pseudo force due to RobinBC
  frobin_ = LINALG::CreateVector(*dofrowmap_, true);

  // have a nice day
  return;
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
void STR::TimIntStatics::EvaluateForceStiffResidual()
{
  // build new external forces
  fextn_->PutScalar(0.0);
  ApplyForceExternal(timen_, disn_, veln_, fextn_);

  // interface forces to external forces
  if (fsisurface_)
  {
    fextn_->Update(1.0, *fifc_, 1.0);  
  }

  // initialise internal forces
  fintn_->PutScalar(0.0);

  // initialise stiffness matrix to zero
  stiff_->Zero();

  // ordinary internal force and stiffness
  ApplyForceStiffInternal(timen_, (*dt_)[0], disn_, disi_, veln_, fintn_, stiff_);

  // apply forces and stiffness due to constraints
  Teuchos::ParameterList pcon; //apply empty parameterlist, no scaling necessary
  ApplyForceStiffConstraint(timen_, (*dis_)(0), disn_, fint_, stiff_, pcon);

  // surface stress force
  ApplyForceStiffSurfstress(disn_, fintn_, stiff_);
  
  // potential forces
  ApplyForceStiffPotential(disn_, fintn_, stiff_);

  // build residual  Res = F_{int;n+1}
  //                     - F_{ext;n+1}
  fres_->Update(-1.0, *fextn_, 0.0);
  fres_->Update(1.0, *fintn_, 1.0);

  // build tangent matrix : effective dynamic stiffness matrix
  //    K_{Teffdyn} = K_{T}
  // i.e. do nothing here
  stiff_->Complete();  // close stiffness matrix

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
  (*dis_)(0)->Norm2(&charnormdis);

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
  fintn_->Norm2(&fintnorm);

  // norm of the external forces
  double fextnorm = 0.0;
  fextn_->Norm2(&fextnorm);

  // return char norm
  return max(fintnorm, fextnorm);
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
void STR::TimIntStatics::UpdateStep()
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

  // update anything that needs to be updated at the element level
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
    discret_->Evaluate(p, null, null, null, null, null);
  }

  // update surface stress
  UpdateStepSurfstress();

  // update potential forces
  UpdateStepPotential();

  // look out
  return;
}

/*----------------------------------------------------------------------*/
/* read restart forces */
void STR::TimIntStatics::ReadRestartForce()
{
  IO::DiscretizationReader reader(discret_, step_);
  // set 'initial' external force
  reader.ReadVector(fext_, "fexternal");
  // set 'initial' internal force vector
  ApplyForceStiffInternal((*time_)[0], (*dt_)[0], (*dis_)(0), zeros_, (*vel_)(0), 
                          fint_, stiff_);
  return;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET

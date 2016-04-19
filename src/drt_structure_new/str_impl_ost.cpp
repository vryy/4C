/*-----------------------------------------------------------*/
/*!
\file str_impl_ost.cpp

\maintainer Philipp Farah

\date Dec 14, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_impl_ost.H"

#include "str_dbc.H"
#include "str_model_evaluator.H"
#include "str_timint_base.H"
#include "str_timint_basedatasdyn.H"
#include "str_utils.H"

#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_sparsematrix.H"

#include <Epetra_Vector.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::IMPLICIT::OneStepTheta::OneStepTheta()
:    theta_(-1.0),
     fvisconp_ptr_(Teuchos::null),
     fviscon_ptr_(Teuchos::null),
     const_vel_acc_update_ptr_(Teuchos::null),
     finertian_ptr_(Teuchos::null),
     finertianp_ptr_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::Setup()
{
  CheckInit();
  // Call the Setup() of the abstract base class first.
  Generic::Setup();

  //cast to one-step-theta data
  const STR::TIMINT::OneStepThetaDataSDyn& onesteptheta_sdyn =
      dynamic_cast<const STR::TIMINT::OneStepThetaDataSDyn&>(TimInt().GetDataSDyn());

  // ---------------------------------------------------------------------------
  // setup time integration parameters
  // ---------------------------------------------------------------------------
  // get a copy of the input parameters
  theta_ = onesteptheta_sdyn.GetTheta();

  // sanity checks and some screen output
  if (GlobalState().GetMyRank()==0)
  {
    if ( (theta_ <= 0.0) or (theta_ > 1.0) )
      dserror("theta out of range (0.0,1.0]");
    else
      std::cout << "   theta = " << theta_ << std::endl;
  }

  // ---------------------------------------------------------------------------
  // setup mid-point vectors
  // ---------------------------------------------------------------------------
  const_vel_acc_update_ptr_ =
      Teuchos::rcp(new Epetra_MultiVector(*GlobalState().DofRowMapView(),2,true));

  // ---------------------------------------------------------------------------
  // setup pointers to the force vectors of the global state data container
  // ---------------------------------------------------------------------------
  finertian_ptr_ = GlobalState().GetMutableFinertialN();
  finertianp_ptr_ = GlobalState().GetMutableFinertialNp();

  fviscon_ptr_ = GlobalState().GetMutableFviscoN();
  fvisconp_ptr_ = GlobalState().GetMutableFviscoNp();

  // ---------------------------------------------------------------------------
  // initialize vectors and matrices
  // ---------------------------------------------------------------------------
  // set the constant parameters for the element evaluation
  if (TimInt().GetDataSDyn().GetMassLinType()==INPAR::STR::ml_rotations)
  {
    dserror("INPAR::STR::ml_rotations is currently unsupported!");
  }

  issetup_ = true;

  PostSetup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::PostSetup()
{
  CheckInitSetup();
  EquilibriateInitialState();
  // Has to be called after the equilibriate state routine!
  UpdateConstantStateContributions();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::SetState(const Epetra_Vector& x)
{
  CheckInitSetup();

  if (IsPredictorState())
    return;

  const double& dt = (*GlobalState().GetDeltaTime())[0];
  // ---------------------------------------------------------------------------
  // new end-point displacements
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> disnp_ptr = GlobalState().ExportDisplEntries(x);
  GlobalState().GetMutableDisNp()->Scale(1.0,*disnp_ptr);

  // ---------------------------------------------------------------------------
  // new end-point velocities
  // ---------------------------------------------------------------------------
  GlobalState().GetMutableVelNp()->Update(1.0,*(*const_vel_acc_update_ptr_)(0),
      1.0/(theta_*dt),*disnp_ptr,0.0);

  // ---------------------------------------------------------------------------
  // new end-point accelerations
  // ---------------------------------------------------------------------------
  GlobalState().GetMutableAccNp()->Update(1.0,*(*const_vel_acc_update_ptr_)(1),
      1.0/(theta_*theta_*dt*dt),*disnp_ptr,0.0);

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::UpdateConstantStateContributions()
{
  const double& dt = (*GlobalState().GetDeltaTime())[0];
  // copy the last converged state into the epetra_multivector
  // ---------------------------------------------------------------------------
  // create a copy of the last converged displacement solution and apply the
  // current DBC values --> D_{n+1}^{0}
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> disnp0_ptr =
      Teuchos::rcp(new Epetra_Vector(*GlobalState().GetDisN()));
  // ---------------------------------------------------------------------------
  // (1) constant velocity update contribution
  // ---------------------------------------------------------------------------
  (*const_vel_acc_update_ptr_)(0)->Scale(1.0/(theta_*dt),
      *GlobalState().GetVelN());
  // ---------------------------------------------------------------------------
  // (2) constant acceleration update contribution
  // ---------------------------------------------------------------------------
  (*const_vel_acc_update_ptr_)(1)->Scale(-(1.0-theta_)/theta_,
      *GlobalState().GetAccN());
  (*const_vel_acc_update_ptr_)(1)->Update(-1.0/(theta_*theta_*dt),
      *GlobalState().GetVelN(),1.0);
  // ---------------------------------------------------------------------------
  // set the current DBC values at the DBC-DoFs of the constant update vectors
  // ---------------------------------------------------------------------------
  Dbc().ApplyDirichletBC(GlobalState().GetTimeNp(),
      disnp0_ptr,                                          // displ.
      Teuchos::rcp((*const_vel_acc_update_ptr_)(0),false),    // veloc.
      Teuchos::rcp((*const_vel_acc_update_ptr_)(1),false),    // accel.
      false);
  /* Finally add the D_{n+1}^{0} vector to the constant velocity and
   * acceleration update vectors: */
  // ---------------------------------------------------------------------------
  // (1) velocity: V_{n+1}^{0} - 1.0/(theta * dt) * D_{n+1}^{0}
  // ---------------------------------------------------------------------------
  (*const_vel_acc_update_ptr_)(0)->Update(-1.0/(theta_*dt),*disnp0_ptr,1.0);
  // ---------------------------------------------------------------------------
  // (2) acceleration: A_{n+1}^{0} - 1.0/(theta^{2} * dt^{2}) * D_{n+1}^{0}
  // ---------------------------------------------------------------------------
  (*const_vel_acc_update_ptr_)(1)->Update(-1.0/(theta_*theta_*dt*dt),*disnp0_ptr,1.0);

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::OneStepTheta::ApplyForce(const Epetra_Vector& x,
        Epetra_Vector& f)
{
  CheckInitSetup();

  // ---------------------------------------------------------------------------
  // evaluate the different model types (static case) at t_{n+1}^{i}
  // ---------------------------------------------------------------------------
  // set the time step dependent parameters for the element evaluation
  ResetEvalParams();
  bool ok = ModelEval().ApplyForce(x,f);
  if (not ok) return ok;

  // ---------------------------------------------------------------------------
  // evaluate the mid state at t_{n+theta}^{i}
  // ---------------------------------------------------------------------------
  EvaluateMidStateForce(f);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::OneStepTheta::ApplyStiff(
    const Epetra_Vector& x,
    LINALG::SparseOperator& jac)
{
  CheckInitSetup();

  // ---------------------------------------------------------------------------
  // evaluate the different model types (static case) at t_{n+1}^{i}
  // ---------------------------------------------------------------------------
  // set the time step dependent parameters for the element evaluation
  ResetEvalParams();
  bool ok = ModelEval().ApplyStiff(x,jac);
  if (not ok) return ok;

  // ---------------------------------------------------------------------------
  // evaluate the mid state at t_{n+theta}^{i}
  // ---------------------------------------------------------------------------
  EvaluateMidStateJacobian(jac);

  jac.Complete();

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::OneStepTheta::ApplyForceStiff(
    const Epetra_Vector& x,
    Epetra_Vector& f,
    LINALG::SparseOperator& jac)
{
  CheckInitSetup();
  // ---------------------------------------------------------------------------
  // evaluate the different model types (static case) at t_{n+1}^{i}
  // ---------------------------------------------------------------------------
  // set the time step dependent parameters for the element evaluation
  ResetEvalParams();
  bool ok = ModelEval().ApplyForceStiff(x,f,jac);
  if (not ok) return ok;

  // ---------------------------------------------------------------------------
  // evaluate the mid state at t_{n+1-theta}^{i}
  // ---------------------------------------------------------------------------
  EvaluateMidStateForce(f);
  EvaluateMidStateJacobian(jac);

  jac.Complete();

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::EvaluateMidStateForce(Epetra_Vector& f) const
{
  // structural dofs of the right-hand-side vector at t_{n} (read-only)
  Teuchos::RCP<const Epetra_Vector> fstructn_ptr =
      GlobalState().GetFstructureN();

  // structural dofs of the right-hand-side vector at t_{n}
  STR::AssembleVector(theta_,f,1.0-theta_,*fstructn_ptr);

  // viscous damping forces at t_{n}
  STR::AssembleVector(1.0,f,1.0-theta_,*fviscon_ptr_);
  // viscous damping forces at t_{n+1}
  STR::AssembleVector(1.0,f,theta_,*fvisconp_ptr_);

  // inertial forces at t_{n}
  STR::AssembleVector(1.0,f,(1.0-theta_),*finertian_ptr_);
  // inertial forces at t_{n+1}
  STR::AssembleVector(1.0,f,(theta_),*finertianp_ptr_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::EvaluateMidStateJacobian(
    LINALG::SparseOperator& jac) const
{
  Teuchos::RCP<LINALG::SparseMatrix> stiff_ptr =
      GlobalState().ExtractDisplBlock(jac);
  const double& dt = (*GlobalState().GetDeltaTime())[0];
  // add inertial contributions and scale the structural stiffness block
  stiff_ptr->Add(*GlobalState().GetMassMatrix(),false,
      1.0/(theta_*dt*dt), theta_);
  // add Rayleigh damping contributions
  if (TimInt().GetDataSDyn().GetDampingType()==INPAR::STR::damp_rayleigh)
    stiff_ptr->Add(*GlobalState().GetDampMatrix(),false,1.0/dt,1.0);

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::IMPLICIT::OneStepTheta::CalcRefNormForce(
    const enum NOX::Abstract::Vector::NormType& type)
{
  dserror("Not yet implemented! (see the Statics integration for an example)");
  return -1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::IMPLICIT::OneStepTheta::GetIntParam() const
{
  dserror("Set the time integration parameter as return value!");
  return -1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::UpdateStepState()
{
  CheckInitSetup();
  // reset structural right hand side.
  //  It will be filled within the model evaluators
  GlobalState().GetMutableFstructureN()->PutScalar(0.0);
  // ---------------------------------------------------------------------------
  // dynamic effects
  // ---------------------------------------------------------------------------
  // new at t_{n+1} -> t_n
  //    finertial_{n} := finertial_{n+1}
  finertian_ptr_->Scale(1.0,*GlobalState().GetFinertialNp());
  // new at t_{n+1} -> t_n
  //    fviscous_{n} := fviscous_{n+1}
  fviscon_ptr_->Scale(1.0,*fvisconp_ptr_);

  // ---------------------------------------------------------------------------
  // update model specific variables
  // ---------------------------------------------------------------------------
  ModelEval().UpdateStepState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::UpdateStepElement()
{
  CheckInitSetup();
  ModelEval().UpdateStepElement();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::PostUpdate()
{
  UpdateConstantStateContributions();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::PredictConstDisConsistVelAcc(
    Epetra_Vector& disnp,
    Epetra_Vector& velnp,
    Epetra_Vector& accnp) const
{
  CheckInitSetup();
  Teuchos::RCP<const Epetra_Vector> disn = GlobalState().GetDisN();
  Teuchos::RCP<const Epetra_Vector> veln = GlobalState().GetVelN();
  Teuchos::RCP<const Epetra_Vector> accn = GlobalState().GetAccN();
  const double& dt = (*GlobalState().GetDeltaTime())[0];

  // constant predictor: displacement in domain
  disnp.Update(1.0,*disn,0.0);

  // consistent velocities
  velnp.Update(1.0,disnp,-1.0,*disn,0.0);
  velnp.Update(-(1.0-theta_)/theta_, *veln,
               1.0/(theta_*dt));

  // consistent accelerations
  accnp.Update(1.0,disnp,-1.0,*disn,0.0);
  accnp.Update(-1.0/(theta_*theta_*dt),*veln,
               -(1.0-theta_)/theta_,*accn,
                1.0/(theta_*theta_*dt*dt));


  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::OneStepTheta::PredictConstVelConsistAcc(
    Epetra_Vector& disnp,
    Epetra_Vector& velnp,
    Epetra_Vector& accnp) const
{
  CheckInitSetup();

  Teuchos::RCP<const Epetra_Vector> disn = GlobalState().GetDisN();
  Teuchos::RCP<const Epetra_Vector> veln = GlobalState().GetVelN();
  Teuchos::RCP<const Epetra_Vector> accn = GlobalState().GetAccN();
  const double& dt = (*GlobalState().GetDeltaTime())[0];

  // extrapolated displacements based upon constant velocities
  // d_{n+1} = d_{n} + dt * v_{n}
  disnp.Update(1.0,*disn,dt,*veln,0.0);

  // consistent velocities
  velnp.Update(1.0,disnp,-1.0,*disn,0.0);
  velnp.Update(-(1.0-theta_)/theta_,*veln,
                1.0/(theta_*dt));

  // consistent accelerations
  accnp.Update(1.0,disnp,-1.0,*disn,0.0);
  accnp.Update(-1.0/(theta_*theta_*dt),*veln,
               -(1.0-theta_)/theta_,*accn,
                1.0/(theta_*theta_*dt*dt));

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::OneStepTheta::PredictConstAcc(
    Epetra_Vector& disnp,
    Epetra_Vector& velnp,
    Epetra_Vector& accnp) const
{
  CheckInitSetup();

  Teuchos::RCP<const Epetra_Vector> disn = GlobalState().GetDisN();
  Teuchos::RCP<const Epetra_Vector> veln = GlobalState().GetVelN();
  Teuchos::RCP<const Epetra_Vector> accn = GlobalState().GetAccN();
  const double& dt = (*GlobalState().GetDeltaTime())[0];

  // extrapolated displacements based upon constant accelerations
  // d_{n+1} = d_{n} + dt * v_{n} + dt^2 / 2 * a_{n}
  disnp.Update(1.0,*disn,dt,*veln, 0.0);
  disnp.Update(0.5*dt*dt,*accn,1.0);

  // consistent velocities
  velnp.Update(1.0,disnp,-1.0,*disn,0.0);
  velnp.Update(-(1.0-theta_)/theta_,*veln,
                1.0/(theta_*dt));

  // consistent accelerations
  accnp.Update(1.0,disnp,-1.0,*disn,0.0);
  accnp.Update(-1.0/(theta_*theta_*dt),*veln,
               -(1.0-theta_)/theta_,*accn,
                1.0/(theta_*theta_*dt*dt));

  return true;
}

/*-----------------------------------------------------------*/
/*!
\file str_impl_ost.cpp

\brief One step theta time integrator.

\maintainer Philipp Farah

\date Dec 14, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_impl_ost.H"

#include "str_dbc.H"
#include "str_model_evaluator.H"
#include "str_model_evaluator_data.H"
#include "str_timint_base.H"
#include "str_timint_basedatasdyn.H"
#include "str_utils.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_io/io.H"
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

  issetup_ = true;

  PostSetup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::PostSetup()
{
  CheckInitSetup();

  if (SDyn().GetMassLinType() != INPAR::STR::ml_rotations and
      !SDyn().NeglectInertia() )
  {
    /* we can use this method for all elements with additive DoFs,
     * but it won't work like this for non-additive rotation vector DoFs */
    EquilibrateInitialState();
  }
  else
  {
    /* If we are restarting the simulation, we get the acceleration state from the
     * restart file. So we are already done at this point. */
    if (TimInt().IsRestarting())
      return;

    // so far, we are restricted to vanishing initial accelerations
    Teuchos::RCP<Epetra_Vector> accnp_ptr = GlobalState().GetMutableAccNp();
    accnp_ptr->PutScalar(0.0);

    // sanity check whether assumption is fulfilled
    /* ToDo tolerance value is experience and based on following consideration:
     * epsilon = O(1e-15) scaled with EA = O(1e8) yields residual contributions in
     * initial, stress free state of order 1e-8 */
    if(not CurrentStateIsEquilibrium(1.0e-6))
      dserror("OneStepTheta only supports initially vanishing acceleration state "
              "in case of ml_rotation = true,\ni.e. an initial state where the system "
              "is equilibrated");

    // call update routines to copy states from t_{n+1} to t_{n}
    // note that the time step is not incremented
    PreUpdate();
    UpdateStepState();
    UpdateStepElement();
    PostUpdate();
  }

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
  Teuchos::RCP<Epetra_Vector> disnp_ptr = GlobalState().ExtractDisplEntries(x);
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
  (*const_vel_acc_update_ptr_)(0)->Scale(-(1.0-theta_)/theta_,
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
  bool ok = ModelEval().ApplyForce(x,f,theta_);
  if (not ok) return ok;

  // ---------------------------------------------------------------------------
  // evaluate the mid state at t_{n+theta}^{i}
  // ---------------------------------------------------------------------------
  AddViscoMassContributions(f);

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
  bool ok = ModelEval().ApplyStiff(x,jac,theta_);
  if (not ok) return ok;

  // ---------------------------------------------------------------------------
  // evaluate the mid state at t_{n+theta}^{i}
  // ---------------------------------------------------------------------------
  AddViscoMassContributions(jac);

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
  bool ok = ModelEval().ApplyForceStiff(x,f,jac,theta_);
  if (not ok) return ok;

  // ---------------------------------------------------------------------------
  // add the visco and mass contributions
  // ---------------------------------------------------------------------------
  AddViscoMassContributions(f);
  AddViscoMassContributions(jac);

  jac.Complete();

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::AddViscoMassContributions(Epetra_Vector& f) const
{
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
void STR::IMPLICIT::OneStepTheta::AddViscoMassContributions(
    LINALG::SparseOperator& jac) const
{
  Teuchos::RCP<LINALG::SparseMatrix> stiff_ptr =
      GlobalState().ExtractDisplBlock(jac);
  const double& dt = (*GlobalState().GetDeltaTime())[0];
  // add inertial contributions and scale the structural stiffness block
  stiff_ptr->Add(*GlobalState().GetMassMatrix(),false,
      1.0/(theta_*dt*dt), 1.0);
  // add Rayleigh damping contributions
  if (TimInt().GetDataSDyn().GetDampingType()==INPAR::STR::damp_rayleigh)
    stiff_ptr->Add(*GlobalState().GetDampMatrix(),false,1.0/dt,1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::WriteRestart(
    IO::DiscretizationWriter& iowriter,
    const bool& forced_writerestart) const
{
  CheckInitSetup();
  // write dynamic forces
  iowriter.WriteVector("finert",finertian_ptr_);
  iowriter.WriteVector("fvisco",fviscon_ptr_);

  ModelEval().WriteRestart(iowriter,forced_writerestart);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::ReadRestart(IO::DiscretizationReader& ioreader)
{

  CheckInitSetup();
  ioreader.ReadVector(finertian_ptr_,"finert");
  ioreader.ReadVector(fviscon_ptr_,"fvisco");

  ModelEval().ReadRestart(ioreader);
  UpdateConstantStateContributions();
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
  return (1.0-theta_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::UpdateStepState()
{
  CheckInitSetup();
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
  ModelEval().UpdateStepState(1-theta_);
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
  disnp.Scale(1.0,*disn);

  // consistent velocities
  /* Since disnp and disn are equal we can skip the current
   * update part and have to consider only the old state at t_{n}.
   *           disnp-disn = 0.0                                 */
  velnp.Scale(-(1.0-theta_)/theta_, *veln);

  // consistent accelerations
  /* Since disnp and disn are equal we can skip the current
   * update part and have to consider only the old state at t_{n}.
   *           disnp-disn = 0.0                                 */
  accnp.Update(-1.0/(theta_*theta_*dt),*veln,
      -(1.0-theta_)/theta_,*accn,
      0.0);


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

  /* extrapolated displacements based upon constant velocities
   * d_{n+1} = d_{n} + dt * v_{n} */
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

  /* extrapolated displacements based upon constant accelerations
   * d_{n+1} = d_{n} + dt * v_{n} + dt^2 / 2 * a_{n} */
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

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::ResetEvalParams()
{
  // call base class
  STR::IMPLICIT::Generic::ResetEvalParams();

  // set the time step dependent parameters for the element evaluation
  const double& dt = (*GlobalState().GetDeltaTime())[0];
  double timeintfac_dis = theta_*theta_*dt*dt;
  double timeintfac_vel = theta_*dt;

  EvalData().SetTimIntFactorDisp(timeintfac_dis);
  EvalData().SetTimIntFactorVel(timeintfac_vel);
}

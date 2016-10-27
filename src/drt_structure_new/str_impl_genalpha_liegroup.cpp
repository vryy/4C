/*-----------------------------------------------------------*/
/*!
\file str_impl_genalpha_liegroup.cpp

\brief Modified Generalized Alpha time integrator for Lie groups
      (can handle non-additive rotational pseudo-vector DoFs)

\maintainer Maximilian Grill

\date Jun 29, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_impl_genalpha_liegroup.H"
#include "str_dbc.H"
#include "str_utils.H"
#include "str_model_evaluator.H"
#include "str_model_evaluator_data.H"
#include "str_timint_base.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_sparsematrix.H"

#include <Epetra_Vector.h>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::IMPLICIT::GenAlphaLieGroup::GenAlphaLieGroup()
    : accn_mod_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlphaLieGroup::Setup()
{
  CheckInit();

  // ---------------------------------------------------------------------------
  // setup additional state vectors of modified acceleration
  // ---------------------------------------------------------------------------
  accn_mod_ = Teuchos::rcp(new Epetra_Vector(*GlobalState().DofRowMapView(),true));

  // Call the Setup() of the parent GenAlpha class
  GenAlpha::Setup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlphaLieGroup::PostSetup()
{
  CheckInitSetup();

  if (SDyn().GetMassLinType() != INPAR::STR::ml_rotations and
      !SDyn().NeglectInertia() )
  {
    /* we can use this method for all elements with additive DoFs,
     * but it won't work like this for non-additive rotation vector DoFs */
    EquilibriateInitialState();
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
      dserror("Lie group GenAlpha only supports initially vanishing acceleration state,"
          " i.e. an initial state where the system is equilibrated");

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
void STR::IMPLICIT::GenAlphaLieGroup::SetState(const Epetra_Vector& x)
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

  /* ToDo in case we want to handle rotation vector DoFs correctly on time
   *      integrator level, the update procedure needs to be adapted here;
   *      use GlobalState().ExportAdditiveEntries() and ExportRotvecEntries() for
   *      this */

  // ---------------------------------------------------------------------------
  // new end-point velocities
  // ---------------------------------------------------------------------------
  GlobalState().GetMutableVelNp()->Update(1.0,*(*const_vel_acc_update_ptr_)(0),
      gamma_/(beta_*dt),*disnp_ptr,0.0);

  // ---------------------------------------------------------------------------
  // new end-point accelerations
  // ---------------------------------------------------------------------------
  GlobalState().GetMutableAccNp()->Update(1.0,*(*const_vel_acc_update_ptr_)(1),
      (1.0-alpham_)/(beta_*dt*dt*(1.0-alphaf_)),*disnp_ptr,0.0);

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlphaLieGroup::WriteRestart(
    IO::DiscretizationWriter& iowriter,
    const bool& forced_writerestart) const
{
  CheckInitSetup();
  // write modified acceleration vector
  iowriter.WriteVector("accn_mod",accn_mod_);

  GenAlpha::WriteRestart(iowriter,forced_writerestart);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlphaLieGroup::ReadRestart(IO::DiscretizationReader& ioreader)
{
  CheckInitSetup();
  ioreader.ReadVector(accn_mod_,"accn_mod");

  GenAlpha::ReadRestart(ioreader);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlphaLieGroup::UpdateStepState()
{
  CheckInitSetup();

  // new at t_{n+1} -> t_n
  //    acc_mod_{n} := -alpha_m/(1-alpha_m) * acc_mod_{n}
  accn_mod_->Scale(-alpham_/(1.0-alpham_));
  accn_mod_->Update(alphaf_/(1.0-alpham_),*GlobalState().GetAccN(),1.0);
  accn_mod_->Update((1.0-alphaf_)/(1.0-alpham_),*GlobalState().GetAccNp(),1.0);

  // ---------------------------------------------------------------------------
  // update model specific variables
  // ---------------------------------------------------------------------------
  ModelEval().UpdateStepState(0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlphaLieGroup::UpdateConstantStateContributions()
{
  const double& dt = (*GlobalState().GetDeltaTime())[0];

  // copy the last converged state into the epetra_multivector
  // ---------------------------------------------------------------------------
  // create a copy of the last converged displacement solution and apply the
  // current DBC values --> D_{n+1}^{0}
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> disnp0_ptr =
      Teuchos::rcp(new Epetra_Vector(*GlobalState().GetDisN()));

  /* ToDo in case we want to handle rotation vector DoFs correctly on time
   *      integrator level, the update procedure needs to be adapted here;
   *      use GlobalState().ExportAdditiveEntries() and ExportRotvecEntries() for
   *      this */

  // ---------------------------------------------------------------------------
  // (1) constant velocity update contribution
  // ---------------------------------------------------------------------------
  (*const_vel_acc_update_ptr_)(0)->Scale(1.0-gamma_/beta_,
      *GlobalState().GetVelN());
  (*const_vel_acc_update_ptr_)(0)->Update((1.0-gamma_/(2.0*beta_))*dt,
      *accn_mod_,1.0);

  // ---------------------------------------------------------------------------
  // (2) constant acceleration update contribution
  // ---------------------------------------------------------------------------
  (*const_vel_acc_update_ptr_)(1)->Scale(alphaf_/(alphaf_-1.0),
      *GlobalState().GetAccN());
  (*const_vel_acc_update_ptr_)(1)->Update(-(1.0-alpham_)/(beta_*dt*(1.0-alphaf_)),
      *GlobalState().GetVelN(),1.0);

  (*const_vel_acc_update_ptr_)(1)->Update(
      alpham_/(1.0-alphaf_) - (1.0-alpham_)*(0.5-beta_)/(beta_*(1.0-alphaf_)),
      *accn_mod_,1.0);

  // ---------------------------------------------------------------------------
  // set the current DBC values at the DBC-DoFs of the constant update vectors
  // ---------------------------------------------------------------------------
  Dbc().ApplyDirichletBC(GlobalState().GetTimeNp(),
      disnp0_ptr,                                             // displ.
      Teuchos::rcp((*const_vel_acc_update_ptr_)(0),false),    // veloc.
      Teuchos::rcp((*const_vel_acc_update_ptr_)(1),false),    // accel.
      false);

  /* Finally add the D_{n+1}^{0} vector to the constant velocity and
   * acceleration update vectors: */
  // ---------------------------------------------------------------------------
  // (1) velocity: V_{n+1}^{0} - gamma/(beta * dt) * D_{n+1}^{0}
  // ---------------------------------------------------------------------------
  (*const_vel_acc_update_ptr_)(0)->Update(-gamma_/(beta_*dt),*disnp0_ptr,1.0);
  // ---------------------------------------------------------------------------
  // (2) acceleration: A_{n+1}^{0} - 1.0/(beta * dt^{2}) * D_{n+1}^{0}
  // ---------------------------------------------------------------------------
  (*const_vel_acc_update_ptr_)(1)->Update(-(1.0-alpham_)/(beta_*dt*dt*(1.0-alphaf_)),
      *disnp0_ptr,1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::IMPLICIT::GenAlphaLieGroup::GetIntParam() const
{
  CheckInitSetup();
  return 0.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlphaLieGroup::AddViscoMassContributions(
    Epetra_Vector& f) const
{
  // viscous damping forces at t_{n+1}
  STR::AssembleVector(1.0,f,1.0,*fvisconp_ptr_);
  // inertia forces at t_{n+1}
  STR::AssembleVector(1.0,f,1.0,*finertianp_ptr_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlphaLieGroup::AddViscoMassContributions(
    LINALG::SparseOperator& jac) const
{
  Teuchos::RCP<LINALG::SparseMatrix> stiff_ptr =
      GlobalState().ExtractDisplBlock(jac);
  const double& dt = (*GlobalState().GetDeltaTime())[0];
  // add inertial contributions to structural stiffness block
  stiff_ptr->Add(*GlobalState().GetMassMatrix(),false,
      (1.0-alpham_)/(beta_*dt*dt*(1.0-alphaf_)),1.0);
  // add Rayleigh damping contributions
  if (TimInt().GetDataSDyn().GetDampingType()==INPAR::STR::damp_rayleigh)
    stiff_ptr->Add(*GlobalState().GetDampMatrix(),false,gamma_/(beta_*dt),1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlphaLieGroup::PredictConstDisConsistVelAcc(
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

  // consistent velocities following Newmark formulas
  velnp.Update(1.0,disnp,-1.0,*disn,0.0);
  velnp.Update((beta_-gamma_)/beta_, *veln,
      (2.0*beta_-gamma_)*dt/(2.0*beta_), *accn_mod_,
      gamma_/(beta_*dt));

  // consistent accelerations following Newmark formulas
  accnp.Update(1.0,disnp,-1.0,*disn,0.0);
  accnp.Update(-(1.0-alpham_)/(beta_*dt*(1-alphaf_)),*veln,
                -alphaf_/(1.0-alphaf_),*accn,
                (1.0-alpham_)/(beta_*dt*dt*(1.0-alphaf_)) );
  accnp.Update(alpham_/(1.0-alphaf_) - (1.0-alpham_)*(0.5-beta_)/(beta_*(1.0-alphaf_)),
                *accn_mod_,1.0);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::GenAlphaLieGroup::PredictConstVelConsistAcc(
    Epetra_Vector& disnp,
    Epetra_Vector& velnp,
    Epetra_Vector& accnp) const
{
  CheckInitSetup();

  dserror("Predictor ConstVelConsistAcc is not supported in Lie group GenAlpha so"
      " far! Use ConstDisConsistVelAcc!");

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::GenAlphaLieGroup::PredictConstAcc(
    Epetra_Vector& disnp,
    Epetra_Vector& velnp,
    Epetra_Vector& accnp) const
{
  CheckInitSetup();

  dserror("Predictor ConstAcc is not supported in Lie group GenAlpha so far! Use "
      "ConstDisConsistVelAcc!");

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlphaLieGroup::ResetEvalParams()
{
  // set the time step dependent parameters for the element evaluation
  GenAlpha::ResetEvalParams();

  /* in case we have non-additive rotation (pseudo-)vector DOFs, we need to pass
   * the GenAlpha parameters to the beam elements via beam parameter interface */
  if (TimInt().GetDataSDyn().GetMassLinType()==INPAR::STR::ml_rotations)
  {
    EvalData().GetBeamData().SetBeta(beta_);
    EvalData().GetBeamData().SetGamma(gamma_);
    EvalData().GetBeamData().SetAlphaf(alphaf_);
    EvalData().GetBeamData().SetAlpham(alpham_);
  }
}

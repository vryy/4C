/*-----------------------------------------------------------*/
/*! \file

\brief One step theta time integrator.


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_impl_ost.hpp"

#include "4C_io.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_structure_new_dbc.hpp"
#include "4C_structure_new_model_evaluator.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_timint_basedatasdyn.hpp"
#include "4C_structure_new_utils.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::IMPLICIT::OneStepTheta::OneStepTheta()
    : theta_(-1.0),
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
  check_init();
  // Call the Setup() of the abstract base class first.
  Generic::Setup();

  // ---------------------------------------------------------------------------
  // setup time integration parameters
  // ---------------------------------------------------------------------------
  // get a copy of the input parameters
  theta_ = get_theta();

  // sanity checks and some screen output
  if (global_state().GetMyRank() == 0)
  {
    if ((theta_ <= 0.0) or (theta_ > 1.0))
      FOUR_C_THROW("theta out of range (0.0,1.0]");
    else
      std::cout << "   theta = " << theta_ << std::endl;
  }

  // ---------------------------------------------------------------------------
  // setup mid-point vectors
  // ---------------------------------------------------------------------------
  const_vel_acc_update_ptr_ =
      Teuchos::rcp(new Epetra_MultiVector(*global_state().DofRowMapView(), 2, true));

  // ---------------------------------------------------------------------------
  // setup pointers to the force vectors of the global state data container
  // ---------------------------------------------------------------------------
  finertian_ptr_ = global_state().GetFinertialN();
  finertianp_ptr_ = global_state().GetFinertialNp();

  fviscon_ptr_ = global_state().GetFviscoN();
  fvisconp_ptr_ = global_state().GetFviscoNp();

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::post_setup()
{
  check_init_setup();

  if (s_dyn().GetMassLinType() != INPAR::STR::ml_rotations and !s_dyn().NeglectInertia())
  {
    /* we can use this method for all elements with additive DoFs,
     * but it won't work like this for non-additive rotation vector DoFs */
    equilibrate_initial_state();
  }
  else
  {
    /* If we are restarting the simulation, we get the acceleration state from the
     * restart file. So we are already done at this point. */
    if (tim_int().IsRestarting()) return;

    // so far, we are restricted to vanishing initial accelerations
    Teuchos::RCP<Epetra_Vector> accnp_ptr = global_state().GetAccNp();
    accnp_ptr->PutScalar(0.0);

    // sanity check whether assumption is fulfilled
    /* ToDo tolerance value is experience and based on following consideration:
     * epsilon = O(1e-15) scaled with EA = O(1e8) yields residual contributions in
     * initial, stress free state of order 1e-8 */
    if (not current_state_is_equilibrium(1.0e-6) and global_state().GetMyRank() == 0)
      std::cout << "\nSERIOUS WARNING: Initially non vanishing acceleration states "
                   "in case of ml_rotation = true,\ni.e. an initial state where the system "
                   "is not equilibrated, cannot yet be computed correctly.\nThis means your "
                   "results in the beginning are not physically correct\n"
                << std::endl;

    // call update routines to copy states from t_{n+1} to t_{n}
    // note that the time step is not incremented
    PreUpdate();
    UpdateStepState();
    UpdateStepElement();
    post_update();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::IMPLICIT::OneStepTheta::get_theta() const
{
  if (is_init() and is_setup()) return theta_;

  const STR::TIMINT::OneStepThetaDataSDyn& onesteptheta_sdyn =
      dynamic_cast<const STR::TIMINT::OneStepThetaDataSDyn&>(tim_int().GetDataSDyn());

  return onesteptheta_sdyn.get_theta();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::set_state(const Epetra_Vector& x)
{
  check_init_setup();

  if (IsPredictorState()) return;

  update_constant_state_contributions();

  const double& dt = (*global_state().GetDeltaTime())[0];
  // ---------------------------------------------------------------------------
  // new end-point displacements
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> disnp_ptr = global_state().ExtractDisplEntries(x);
  global_state().GetDisNp()->Scale(1.0, *disnp_ptr);

  // ---------------------------------------------------------------------------
  // new end-point velocities
  // ---------------------------------------------------------------------------
  global_state().GetVelNp()->Update(
      1.0, *(*const_vel_acc_update_ptr_)(0), 1.0 / (theta_ * dt), *disnp_ptr, 0.0);

  // ---------------------------------------------------------------------------
  // new end-point accelerations
  // ---------------------------------------------------------------------------
  global_state().GetAccNp()->Update(
      1.0, *(*const_vel_acc_update_ptr_)(1), 1.0 / (theta_ * theta_ * dt * dt), *disnp_ptr, 0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::update_constant_state_contributions()
{
  const double& dt = (*global_state().GetDeltaTime())[0];

  // ---------------------------------------------------------------------------
  // velocity
  // ---------------------------------------------------------------------------
  (*const_vel_acc_update_ptr_)(0)->Scale(-(1.0 - theta_) / theta_, *global_state().GetVelN());
  (*const_vel_acc_update_ptr_)(0)->Update(-1.0 / (theta_ * dt), *global_state().GetDisN(), 1.0);

  // ---------------------------------------------------------------------------
  // acceleration
  // ---------------------------------------------------------------------------
  (*const_vel_acc_update_ptr_)(1)->Scale(-(1.0 - theta_) / theta_, *global_state().GetAccN());
  (*const_vel_acc_update_ptr_)(1)->Update(
      -1.0 / (theta_ * theta_ * dt), *global_state().GetVelN(), 1.0);
  (*const_vel_acc_update_ptr_)(1)->Update(
      -1.0 / (theta_ * theta_ * dt * dt), *global_state().GetDisN(), 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::OneStepTheta::ApplyForce(const Epetra_Vector& x, Epetra_Vector& f)
{
  check_init_setup();

  // ---------------------------------------------------------------------------
  // evaluate the different model types (static case) at t_{n+1}^{i}
  // ---------------------------------------------------------------------------
  // set the time step dependent parameters for the element evaluation
  reset_eval_params();
  return ModelEval().ApplyForce(x, f, theta_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::OneStepTheta::ApplyStiff(
    const Epetra_Vector& x, CORE::LINALG::SparseOperator& jac)
{
  check_init_setup();

  // ---------------------------------------------------------------------------
  // evaluate the different model types (static case) at t_{n+1}^{i}
  // ---------------------------------------------------------------------------
  // set the time step dependent parameters for the element evaluation
  reset_eval_params();
  const bool ok = ModelEval().ApplyStiff(x, jac, theta_);

  if (not ok) return ok;

  jac.Complete();

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::OneStepTheta::ApplyForceStiff(
    const Epetra_Vector& x, Epetra_Vector& f, CORE::LINALG::SparseOperator& jac)
{
  check_init_setup();
  // ---------------------------------------------------------------------------
  // evaluate the different model types (static case) at t_{n+1}^{i}
  // ---------------------------------------------------------------------------
  // set the time step dependent parameters for the element evaluation
  reset_eval_params();
  const bool ok = ModelEval().ApplyForceStiff(x, f, jac, theta_);

  if (not ok) return ok;

  jac.Complete();

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::OneStepTheta::assemble_force(
    Epetra_Vector& f, const std::vector<INPAR::STR::ModelType>* without_these_models) const
{
  return ModelEval().assemble_force(theta_, f, without_these_models);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::add_visco_mass_contributions(Epetra_Vector& f) const
{
  // viscous damping forces at t_{n}
  CORE::LINALG::AssembleMyVector(1.0, f, 1.0 - theta_, *fviscon_ptr_);
  // viscous damping forces at t_{n+1}
  CORE::LINALG::AssembleMyVector(1.0, f, theta_, *fvisconp_ptr_);

  // inertial forces at t_{n}
  CORE::LINALG::AssembleMyVector(1.0, f, (1.0 - theta_), *finertian_ptr_);
  // inertial forces at t_{n+1}
  CORE::LINALG::AssembleMyVector(1.0, f, (theta_), *finertianp_ptr_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::add_visco_mass_contributions(
    CORE::LINALG::SparseOperator& jac) const
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> stiff_ptr = global_state().ExtractDisplBlock(jac);
  const double& dt = (*global_state().GetDeltaTime())[0];
  // add inertial contributions and scale the structural stiffness block
  stiff_ptr->Add(*global_state().GetMassMatrix(), false, 1.0 / (theta_ * dt * dt), 1.0);
  // add Rayleigh damping contributions
  if (tim_int().GetDataSDyn().GetDampingType() == INPAR::STR::damp_rayleigh)
    stiff_ptr->Add(*global_state().GetDampMatrix(), false, 1.0 / dt, 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::write_restart(
    IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  check_init_setup();
  // write dynamic forces
  iowriter.WriteVector("finert", finertian_ptr_);
  iowriter.WriteVector("fvisco", fviscon_ptr_);

  ModelEval().write_restart(iowriter, forced_writerestart);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::read_restart(IO::DiscretizationReader& ioreader)
{
  check_init_setup();
  ioreader.ReadVector(finertian_ptr_, "finert");
  ioreader.ReadVector(fviscon_ptr_, "fvisco");

  ModelEval().read_restart(ioreader);
  update_constant_state_contributions();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::IMPLICIT::OneStepTheta::CalcRefNormForce(
    const enum ::NOX::Abstract::Vector::NormType& type) const
{
  FOUR_C_THROW("Not yet implemented! (see the Statics integration for an example)");
  return -1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::IMPLICIT::OneStepTheta::GetIntParam() const { return (1.0 - get_theta()); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::UpdateStepState()
{
  check_init_setup();
  // ---------------------------------------------------------------------------
  // dynamic effects
  // ---------------------------------------------------------------------------
  // new at t_{n+1} -> t_n
  //    finertial_{n} := finertial_{n+1}
  finertian_ptr_->Scale(1.0, *global_state().GetFinertialNp());
  // new at t_{n+1} -> t_n
  //    fviscous_{n} := fviscous_{n+1}
  fviscon_ptr_->Scale(1.0, *fvisconp_ptr_);

  // ---------------------------------------------------------------------------
  // update model specific variables
  // ---------------------------------------------------------------------------
  ModelEval().UpdateStepState(1 - theta_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::UpdateStepElement()
{
  check_init_setup();
  ModelEval().UpdateStepElement();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::post_update() { update_constant_state_contributions(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::predict_const_dis_consist_vel_acc(
    Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const
{
  check_init_setup();
  Teuchos::RCP<const Epetra_Vector> disn = global_state().GetDisN();
  Teuchos::RCP<const Epetra_Vector> veln = global_state().GetVelN();
  Teuchos::RCP<const Epetra_Vector> accn = global_state().GetAccN();
  const double& dt = (*global_state().GetDeltaTime())[0];

  // constant predictor: displacement in domain
  disnp.Scale(1.0, *disn);

  // consistent velocities
  /* Since disnp and disn are equal we can skip the current
   * update part and have to consider only the old state at t_{n}.
   *           disnp-disn = 0.0                                 */
  velnp.Scale(-(1.0 - theta_) / theta_, *veln);

  // consistent accelerations
  /* Since disnp and disn are equal we can skip the current
   * update part and have to consider only the old state at t_{n}.
   *           disnp-disn = 0.0                                 */
  accnp.Update(-1.0 / (theta_ * theta_ * dt), *veln, -(1.0 - theta_) / theta_, *accn, 0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::OneStepTheta::predict_const_vel_consist_acc(
    Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const
{
  check_init_setup();

  Teuchos::RCP<const Epetra_Vector> disn = global_state().GetDisN();
  Teuchos::RCP<const Epetra_Vector> veln = global_state().GetVelN();
  Teuchos::RCP<const Epetra_Vector> accn = global_state().GetAccN();
  const double& dt = (*global_state().GetDeltaTime())[0];

  /* extrapolated displacements based upon constant velocities
   * d_{n+1} = d_{n} + dt * v_{n} */
  disnp.Update(1.0, *disn, dt, *veln, 0.0);

  // consistent velocities
  velnp.Update(1.0, disnp, -1.0, *disn, 0.0);
  velnp.Update(-(1.0 - theta_) / theta_, *veln, 1.0 / (theta_ * dt));

  // consistent accelerations
  accnp.Update(1.0, disnp, -1.0, *disn, 0.0);
  accnp.Update(-1.0 / (theta_ * theta_ * dt), *veln, -(1.0 - theta_) / theta_, *accn,
      1.0 / (theta_ * theta_ * dt * dt));

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::OneStepTheta::PredictConstAcc(
    Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const
{
  check_init_setup();

  Teuchos::RCP<const Epetra_Vector> disn = global_state().GetDisN();
  Teuchos::RCP<const Epetra_Vector> veln = global_state().GetVelN();
  Teuchos::RCP<const Epetra_Vector> accn = global_state().GetAccN();
  const double& dt = (*global_state().GetDeltaTime())[0];

  /* extrapolated displacements based upon constant accelerations
   * d_{n+1} = d_{n} + dt * v_{n} + dt^2 / 2 * a_{n} */
  disnp.Update(1.0, *disn, dt, *veln, 0.0);
  disnp.Update(0.5 * dt * dt, *accn, 1.0);

  // consistent velocities
  velnp.Update(1.0, disnp, -1.0, *disn, 0.0);
  velnp.Update(-(1.0 - theta_) / theta_, *veln, 1.0 / (theta_ * dt));

  // consistent accelerations
  accnp.Update(1.0, disnp, -1.0, *disn, 0.0);
  accnp.Update(-1.0 / (theta_ * theta_ * dt), *veln, -(1.0 - theta_) / theta_, *accn,
      1.0 / (theta_ * theta_ * dt * dt));

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::reset_eval_params()
{
  // call base class
  STR::IMPLICIT::Generic::reset_eval_params();

  // set the time step dependent parameters for the element evaluation
  const double& dt = (*global_state().GetDeltaTime())[0];
  double timeintfac_dis = theta_ * theta_ * dt * dt;
  double timeintfac_vel = theta_ * dt;

  EvalData().SetTimIntFactorDisp(timeintfac_dis);
  EvalData().SetTimIntFactorVel(timeintfac_vel);
}

FOUR_C_NAMESPACE_CLOSE

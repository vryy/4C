/*---------------------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all 0D cardiovascular model terms


\level 3

*/
/*---------------------------------------------------------------------*/

#include "4C_cardiovascular0d_structure_new_model_evaluator.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_structure_new_integrator.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_parameter_list.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Solid::MODELEVALUATOR::Cardiovascular0D::Cardiovascular0D()
    : disnp_ptr_(Teuchos::null),
      stiff_cardio_ptr_(Teuchos::null),
      fstructcardio_np_ptr_(Teuchos::null)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::MODELEVALUATOR::Cardiovascular0D::setup()
{
  check_init();

  Teuchos::RCP<Core::FE::Discretization> dis = discret_ptr();

  // setup the displacement pointer
  disnp_ptr_ = global_state().get_dis_np();

  // contributions of 0D model to structural rhs and stiffness
  fstructcardio_np_ptr_ = Teuchos::rcp(new Epetra_Vector(*global_state().dof_row_map_view()));
  stiff_cardio_ptr_ = Teuchos::rcp(
      new Core::LinAlg::SparseMatrix(*global_state().dof_row_map_view(), 81, true, true));

  Teuchos::ParameterList solvparams;
  Core::UTILS::AddEnumClassToParameterList<Core::LinearSolver::SolverType>(
      "SOLVER", Core::LinearSolver::SolverType::umfpack, solvparams);
  Teuchos::RCP<Core::LinAlg::Solver> dummysolver(new Core::LinAlg::Solver(
      solvparams, disnp_ptr_->Comm(), nullptr, Core::IO::Verbositylevel::standard));

  // ToDo: we do not want to hand in the structural dynamics parameter list
  // to the manager in the future! -> get rid of it as soon as old
  // time-integration dies ...
  // initialize 0D cardiovascular manager
  cardvasc0dman_ = Teuchos::rcp(new UTILS::Cardiovascular0DManager(dis, disnp_ptr_,
      Global::Problem::instance()->structural_dynamic_params(),
      Global::Problem::instance()->cardiovascular0_d_structural_params(), *dummysolver,
      Teuchos::null));

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::MODELEVALUATOR::Cardiovascular0D::reset(const Epetra_Vector& x)
{
  check_init_setup();

  // update the structural displacement vector
  disnp_ptr_ = global_state().get_dis_np();

  fstructcardio_np_ptr_->PutScalar(0.0);
  stiff_cardio_ptr_->zero();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Solid::MODELEVALUATOR::Cardiovascular0D::evaluate_force()
{
  check_init_setup();

  double time_np = global_state().get_time_np();
  Teuchos::ParameterList pcardvasc0d;
  pcardvasc0d.set("time_step_size", (*global_state().get_delta_time())[0]);

  // only forces are evaluated!
  cardvasc0dman_->evaluate_force_stiff(
      time_np, disnp_ptr_, fstructcardio_np_ptr_, Teuchos::null, pcardvasc0d);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Solid::MODELEVALUATOR::Cardiovascular0D::evaluate_stiff()
{
  check_init_setup();

  double time_np = global_state().get_time_np();
  Teuchos::ParameterList pcardvasc0d;
  pcardvasc0d.set("time_step_size", (*global_state().get_delta_time())[0]);

  // only stiffnesses are evaluated!
  cardvasc0dman_->evaluate_force_stiff(
      time_np, disnp_ptr_, Teuchos::null, stiff_cardio_ptr_, pcardvasc0d);

  if (not stiff_cardio_ptr_->filled()) stiff_cardio_ptr_->complete();

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Solid::MODELEVALUATOR::Cardiovascular0D::evaluate_force_stiff()
{
  check_init_setup();

  double time_np = global_state().get_time_np();
  Teuchos::ParameterList pcardvasc0d;
  pcardvasc0d.set("time_step_size", (*global_state().get_delta_time())[0]);

  cardvasc0dman_->evaluate_force_stiff(
      time_np, disnp_ptr_, fstructcardio_np_ptr_, stiff_cardio_ptr_, pcardvasc0d);

  if (not stiff_cardio_ptr_->filled()) stiff_cardio_ptr_->complete();

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Solid::MODELEVALUATOR::Cardiovascular0D::assemble_force(
    Epetra_Vector& f, const double& timefac_np) const
{
  Teuchos::RCP<const Epetra_Vector> block_vec_ptr = Teuchos::null;

  // assemble and scale with str time-integrator dependent value
  Core::LinAlg::AssembleMyVector(1.0, f, timefac_np, *fstructcardio_np_ptr_);

  // assemble 0D model rhs - already at the generalized mid-point t_{n+theta} !
  block_vec_ptr = cardvasc0dman_->get_cardiovascular0_drhs();

  if (block_vec_ptr.is_null())
    FOUR_C_THROW(
        "The 0D cardiovascular model vector is a nullptr pointer, although \n"
        "the structural part indicates, that 0D cardiovascular model contributions \n"
        "are present!");

  const int elements_f = f.Map().NumGlobalElements();
  const int max_gid = get_block_dof_row_map_ptr()->MaxAllGID();
  // only call when f is the full rhs of the coupled problem (not for structural
  // equilibriate initial state call)
  if (elements_f == max_gid + 1) Core::LinAlg::AssembleMyVector(1.0, f, 1.0, *block_vec_ptr);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Solid::MODELEVALUATOR::Cardiovascular0D::assemble_jacobian(
    Core::LinAlg::SparseOperator& jac, const double& timefac_np) const
{
  Teuchos::RCP<Core::LinAlg::SparseMatrix> block_ptr = Teuchos::null;

  // --- Kdd - block - scale with str time-integrator dependent value---
  Teuchos::RCP<Core::LinAlg::SparseMatrix> jac_dd_ptr = global_state().extract_displ_block(jac);
  jac_dd_ptr->add(*stiff_cardio_ptr_, false, timefac_np, 1.0);
  // no need to keep it
  stiff_cardio_ptr_->zero();

  // --- Kdz - block ---------------------------------------------------
  block_ptr = cardvasc0dman_->get_mat_dstruct_dcv0ddof();
  // scale with str time-integrator dependent value
  block_ptr->scale(timefac_np);
  global_state().assign_model_block(jac, *block_ptr, type(), MatBlockType::displ_lm);
  // reset the block pointer, just to be on the safe side
  block_ptr = Teuchos::null;

  // --- Kzd - block - already scaled correctly by 0D model !-----------
  block_ptr = cardvasc0dman_->get_mat_dcardvasc0d_dd()->transpose();
  global_state().assign_model_block(jac, *block_ptr, type(), MatBlockType::lm_displ);
  // reset the block pointer, just to be on the safe side
  block_ptr = Teuchos::null;

  // --- Kzz - block - already scaled with 0D theta by 0D model !-------
  block_ptr = cardvasc0dman_->get_cardiovascular0_d_stiffness();
  global_state().assign_model_block(jac, *block_ptr, type(), MatBlockType::lm_lm);
  // reset the block pointer, just to be on the safe side
  block_ptr = Teuchos::null;

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::MODELEVALUATOR::Cardiovascular0D::write_restart(
    Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  iowriter.write_vector("cv0d_df_np", cardvasc0dman_->get0_d_df_np());
  iowriter.write_vector("cv0d_f_np", cardvasc0dman_->get0_d_f_np());

  iowriter.write_vector("cv0d_dof_np", cardvasc0dman_->get0_d_dof_np());
  iowriter.write_vector("vol_np", cardvasc0dman_->get0_d_vol_np());

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::MODELEVALUATOR::Cardiovascular0D::read_restart(Core::IO::DiscretizationReader& ioreader)
{
  double time_n = global_state().get_time_n();
  cardvasc0dman_->read_restart(ioreader, time_n);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::MODELEVALUATOR::Cardiovascular0D::run_post_compute_x(
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  check_init_setup();

  Teuchos::RCP<Epetra_Vector> cv0d_incr =
      global_state().extract_model_entries(Inpar::Solid::model_cardiovascular0d, dir);

  cardvasc0dman_->update_cv0_d_dof(cv0d_incr);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::MODELEVALUATOR::Cardiovascular0D::update_step_state(const double& timefac_n)
{
  // only update 0D model at the end of the time step!
  if (eval_data().get_total_time() == global_state().get_time_np())
    cardvasc0dman_->update_time_step();

  // only print state variables after a finished time step, not when we're
  // in the equilibriate initial state routine
  if (eval_data().get_total_time() == global_state().get_time_np())
    cardvasc0dman_->print_pres_flux(false);

  // add the 0D cardiovascular force contributions to the old structural
  // residual state vector
  if (not fstructcardio_np_ptr_.is_null())
  {
    Teuchos::RCP<Epetra_Vector>& fstructold_ptr = global_state().get_fstructure_old();
    fstructold_ptr->Update(timefac_n, *fstructcardio_np_ptr_, 1.0);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::MODELEVALUATOR::Cardiovascular0D::update_step_element()
{
  // nothing to do
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::MODELEVALUATOR::Cardiovascular0D::determine_stress_strain()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::MODELEVALUATOR::Cardiovascular0D::determine_energy()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::MODELEVALUATOR::Cardiovascular0D::determine_optional_quantity()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::MODELEVALUATOR::Cardiovascular0D::output_step_state(
    Core::IO::DiscretizationWriter& iowriter) const
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::MODELEVALUATOR::Cardiovascular0D::reset_step_state()
{
  check_init_setup();

  FOUR_C_THROW("Not yet implemented");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> Solid::MODELEVALUATOR::Cardiovascular0D::get_block_dof_row_map_ptr()
    const
{
  check_init_setup();

  return cardvasc0dman_->get_cardiovascular0_d_map();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
Solid::MODELEVALUATOR::Cardiovascular0D::get_current_solution_ptr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
Solid::MODELEVALUATOR::Cardiovascular0D::get_last_time_step_solution_ptr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Solid::MODELEVALUATOR::Cardiovascular0D::post_output()
{
  check_init_setup();
  // empty

  return;
}  // post_output()

FOUR_C_NAMESPACE_CLOSE

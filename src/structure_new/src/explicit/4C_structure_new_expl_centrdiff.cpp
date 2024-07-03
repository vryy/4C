/*-----------------------------------------------------------*/
/*! \file

\brief Central differences time integration for solid dynamics


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_expl_centrdiff.hpp"

#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_structure_new_model_evaluator.hpp"
#include "4C_structure_new_model_evaluator_generic.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Solid::EXPLICIT::CentrDiff::CentrDiff()
    : fvisconp_ptr_(Teuchos::null),
      fviscon_ptr_(Teuchos::null),
      finertianp_ptr_(Teuchos::null),
      finertian_ptr_(Teuchos::null)
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::EXPLICIT::CentrDiff::setup()
{
  check_init();

  // Call the setup() of the abstract base class first.
  Generic::setup();

  // ---------------------------------------------------------------------------
  // setup pointers to the force vectors of the global state data container
  // ---------------------------------------------------------------------------
  finertian_ptr_ = global_state().get_finertial_n();
  finertianp_ptr_ = global_state().get_finertial_np();

  fviscon_ptr_ = global_state().get_fvisco_n();
  fvisconp_ptr_ = global_state().get_fvisco_np();

  // -------------------------------------------------------------------
  // set initial displacement
  // -------------------------------------------------------------------
  set_initial_displacement(
      tim_int().get_data_sdyn().get_initial_disp(), tim_int().get_data_sdyn().StartFuncNo());

  // Has to be set before the post_setup() routine is called!
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::EXPLICIT::CentrDiff::post_setup()
{
  check_init_setup();
  equilibrate_initial_state();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::EXPLICIT::CentrDiff::set_state(const Epetra_Vector& x)
{
  check_init_setup();

  const double dt = (*global_state().get_delta_time())[0];
  const double dthalf = dt / 2.0;

  model_eval().reset_step_state();

  // ---------------------------------------------------------------------------
  // new end-point acceleration
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> accnp_ptr = global_state().extract_displ_entries(x);
  global_state().get_acc_np()->Scale(1.0, *accnp_ptr);

  // ---------------------------------------------------------------------------
  // new half-point velocities
  // ---------------------------------------------------------------------------
  global_state().get_vel_np()->Update(1.0, *global_state().get_vel_n(), 0.0);
  global_state().get_vel_np()->Update(dthalf, *global_state().get_acc_n(), 1.0);

  // ---------------------------------------------------------------------------
  // new end-point displacements
  // ---------------------------------------------------------------------------
  global_state().get_dis_np()->Update(1.0, *global_state().get_dis_n(), 0.0);
  global_state().get_dis_np()->Update(dt, *global_state().get_vel_np(), 1.0);

  // ---------------------------------------------------------------------------
  // update the elemental state
  // ---------------------------------------------------------------------------
  model_eval().update_residual();
  model_eval().run_recover();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::EXPLICIT::CentrDiff::add_visco_mass_contributions(Epetra_Vector& f) const
{
  // viscous damping forces at t_{n+1}
  Core::LinAlg::AssembleMyVector(1.0, f, 1.0, *fvisconp_ptr_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::EXPLICIT::CentrDiff::add_visco_mass_contributions(
    Core::LinAlg::SparseOperator& jac) const
{
  Teuchos::RCP<Core::LinAlg::SparseMatrix> stiff_ptr = global_state().extract_displ_block(jac);
  // set mass matrix
  stiff_ptr->Add(*global_state().get_mass_matrix(), false, 1.0, 0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::EXPLICIT::CentrDiff::write_restart(
    Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  check_init_setup();
  // write dynamic forces
  iowriter.write_vector("finert", finertian_ptr_);
  iowriter.write_vector("fvisco", fviscon_ptr_);

  model_eval().write_restart(iowriter, forced_writerestart);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::EXPLICIT::CentrDiff::read_restart(Core::IO::DiscretizationReader& ioreader)
{
  check_init_setup();
  ioreader.read_vector(finertian_ptr_, "finert");
  ioreader.read_vector(fviscon_ptr_, "fvisco");

  model_eval().read_restart(ioreader);
  update_constant_state_contributions();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::EXPLICIT::CentrDiff::update_step_state()
{
  check_init_setup();

  const double dt = (*global_state().get_delta_time())[0];
  const double dthalf = dt / 2.0;

  // ---------------------------------------------------------------------------
  // dynamic effects
  // ---------------------------------------------------------------------------
  // new at t_{n+1} -> t_n
  //    finertial_{n} := finertial_{n+1}
  finertian_ptr_->Scale(1.0, *finertianp_ptr_);
  // new at t_{n+1} -> t_n
  //    fviscous_{n} := fviscous_{n+1}
  fviscon_ptr_->Scale(1.0, *fvisconp_ptr_);

  // recompute the velocity to account for new acceleration
  global_state().get_vel_np()->Update(1.0, *global_state().get_vel_n(), 0.0);
  global_state().get_vel_np()->Update(dthalf, *global_state().get_acc_n(), 1.0);
  global_state().get_vel_np()->Update(dthalf, *global_state().get_acc_np(), 1.0);

  // ---------------------------------------------------------------------------
  // update model specific variables
  // ---------------------------------------------------------------------------
  model_eval().update_step_state(0.0);
}

FOUR_C_NAMESPACE_CLOSE

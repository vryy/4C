/*-----------------------------------------------------------*/
/*! \file

\brief Adams-Bashforth-2 time integration for solid dynamics


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_expl_ab2.hpp"

#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_structure_new_model_evaluator.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::EXPLICIT::AdamsBashforth2::AdamsBashforth2()
    : fvisconp_ptr_(Teuchos::null),
      fviscon_ptr_(Teuchos::null),
      finertianp_ptr_(Teuchos::null),
      finertian_ptr_(Teuchos::null)
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::AdamsBashforth2::Setup()
{
  check_init();

  // Call the Setup() of the abstract base class first.
  Generic::Setup();

  // ---------------------------------------------------------------------------
  // setup pointers to the force vectors of the global state data container
  // ---------------------------------------------------------------------------
  finertian_ptr_ = global_state().GetFinertialN();
  finertianp_ptr_ = global_state().GetFinertialNp();

  fviscon_ptr_ = global_state().GetFviscoN();
  fvisconp_ptr_ = global_state().GetFviscoNp();

  // ---------------------------------------------------------------------------
  // resizing of multi-step quantities
  // ---------------------------------------------------------------------------
  global_state().GetMultiTime()->Resize(-1, 0, true);
  global_state().GetDeltaTime()->Resize(-1, 0, true);
  global_state().GetMultiDis()->Resize(-1, 0, global_state().DofRowMapView(), true);
  global_state().GetMultiVel()->Resize(-1, 0, global_state().DofRowMapView(), true);
  global_state().GetMultiAcc()->Resize(-1, 0, global_state().DofRowMapView(), true);

  // here we initialized the dt of previous steps in the database, since a resize is performed
  const double dt = (*global_state().GetDeltaTime())[0];
  global_state().GetDeltaTime()->UpdateSteps(dt);

  // -------------------------------------------------------------------
  // set initial displacement
  // -------------------------------------------------------------------
  set_initial_displacement(
      tim_int().GetDataSDyn().GetInitialDisp(), tim_int().GetDataSDyn().StartFuncNo());

  // Has to be set before the post_setup() routine is called!
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::AdamsBashforth2::post_setup()
{
  check_init_setup();
  equilibrate_initial_state();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::AdamsBashforth2::set_state(const Epetra_Vector& x)
{
  check_init_setup();

  const double dt = (*global_state().GetDeltaTime())[0];
  const double dto = (*global_state().GetDeltaTime())[-1];
  const double dta = (2.0 * dt * dto + dt * dt) / (2.0 * dto);
  const double dtb = -(dt * dt) / (2.0 * dto);

  // ---------------------------------------------------------------------------
  // new end-point acceleration
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> accnp_ptr = global_state().ExtractDisplEntries(x);
  global_state().GetAccNp()->Scale(1.0, *accnp_ptr);

  // ---------------------------------------------------------------------------
  // new end-point velocities
  // ---------------------------------------------------------------------------
  global_state().GetVelNp()->Update(1.0, (*(global_state().GetMultiVel()))[0], 0.0);
  global_state().GetVelNp()->Update(
      dta, (*(global_state().GetMultiAcc()))[0], dtb, (*(global_state().GetMultiAcc()))[-1], 1.0);

  // ---------------------------------------------------------------------------
  // new end-point displacements
  // ---------------------------------------------------------------------------
  global_state().GetDisNp()->Update(1.0, (*(global_state().GetMultiDis()))[0], 0.0);
  global_state().GetDisNp()->Update(
      dta, (*(global_state().GetMultiVel()))[0], dtb, (*(global_state().GetMultiVel()))[-1], 1.0);

  // ---------------------------------------------------------------------------
  // update the elemental state
  // ---------------------------------------------------------------------------
  ModelEval().update_residual();
  ModelEval().RunRecover();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::AdamsBashforth2::add_visco_mass_contributions(Epetra_Vector& f) const
{
  // viscous damping forces at t_{n+1}
  CORE::LINALG::AssembleMyVector(1.0, f, 1.0, *fvisconp_ptr_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::AdamsBashforth2::add_visco_mass_contributions(
    CORE::LINALG::SparseOperator& jac) const
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> stiff_ptr = global_state().ExtractDisplBlock(jac);
  // set mass matrix
  stiff_ptr->Add(*global_state().GetMassMatrix(), false, 1.0, 0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::AdamsBashforth2::write_restart(
    CORE::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  check_init_setup();
  // write dynamic forces
  iowriter.WriteVector("finert", finertian_ptr_);
  iowriter.WriteVector("fvisco", fviscon_ptr_);

  ModelEval().write_restart(iowriter, forced_writerestart);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::AdamsBashforth2::read_restart(CORE::IO::DiscretizationReader& ioreader)
{
  check_init_setup();
  ioreader.ReadVector(finertian_ptr_, "finert");
  ioreader.ReadVector(fviscon_ptr_, "fvisco");

  ModelEval().read_restart(ioreader);
  update_constant_state_contributions();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::AdamsBashforth2::UpdateStepState()
{
  check_init_setup();

  // ---------------------------------------------------------------------------
  // dynamic effects
  // ---------------------------------------------------------------------------
  // new at t_{n+1} -> t_n
  //    finertial_{n} := finertial_{n+1}
  finertian_ptr_->Scale(1.0, *finertianp_ptr_);
  // new at t_{n+1} -> t_n
  //    fviscous_{n} := fviscous_{n+1}
  fviscon_ptr_->Scale(1.0, *fvisconp_ptr_);

  // ---------------------------------------------------------------------------
  // update model specific variables
  // ---------------------------------------------------------------------------
  ModelEval().UpdateStepState(0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::EXPLICIT::AdamsBashforth2::method_lin_err_coeff_dis() const
{
  const double dt = (*global_state().GetDeltaTime())[0];
  const double dto = (*global_state().GetDeltaTime())[-1];
  return (2. * dt + 3. * dto) / (12. * dt);
}

FOUR_C_NAMESPACE_CLOSE

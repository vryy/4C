/*-----------------------------------------------------------*/
/*! \file

\brief Forward Euler time integration for solid dynamics


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_expl_forwardeuler.hpp"

#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_structure_new_model_evaluator.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::EXPLICIT::ForwardEuler::ForwardEuler()
    : modexpleuler_(true),
      fvisconp_ptr_(Teuchos::null),
      fviscon_ptr_(Teuchos::null),
      finertianp_ptr_(Teuchos::null),
      finertian_ptr_(Teuchos::null)
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::ForwardEuler::Setup()
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

  // -------------------------------------------------------------------
  // set initial displacement
  // -------------------------------------------------------------------
  set_initial_displacement(
      tim_int().GetDataSDyn().GetInitialDisp(), tim_int().GetDataSDyn().StartFuncNo());

  // mode of Forward Euler interpolation
  modexpleuler_ = dynamic_cast<const STR::TIMINT::ExplEulerDataSDyn&>(tim_int().GetDataSDyn())
                      .get_modified_forward_euler();

  // Has to be set before the post_setup() routine is called!
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::ForwardEuler::post_setup()
{
  check_init_setup();
  equilibrate_initial_state();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::ForwardEuler::set_state(const Epetra_Vector& x)
{
  check_init_setup();

  const double dt = (*global_state().GetDeltaTime())[0];

  // ---------------------------------------------------------------------------
  // new end-point acceleration
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> accnp_ptr = global_state().ExtractDisplEntries(x);
  global_state().GetAccNp()->Scale(1.0, *accnp_ptr);

  // ---------------------------------------------------------------------------
  // new end-point velocities
  // ---------------------------------------------------------------------------
  global_state().GetVelNp()->Update(1.0, *global_state().GetVelN(), 0.0);
  global_state().GetVelNp()->Update(dt, *global_state().GetAccN(), 1.0);

  // ---------------------------------------------------------------------------
  // new end-point displacements
  // ---------------------------------------------------------------------------
  global_state().GetDisNp()->Update(1.0, *global_state().GetDisN(), 0.0);
  if (modexpleuler_ == true)
    global_state().GetDisNp()->Update(dt, *global_state().GetVelNp(), 1.0);
  else
    global_state().GetDisNp()->Update(dt, *global_state().GetVelN(), 1.0);

  // ---------------------------------------------------------------------------
  // update the elemental state
  // ---------------------------------------------------------------------------
  ModelEval().update_residual();
  ModelEval().RunRecover();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::ForwardEuler::add_visco_mass_contributions(Epetra_Vector& f) const
{
  // viscous damping forces at t_{n+1}
  CORE::LINALG::AssembleMyVector(1.0, f, 1.0, *fvisconp_ptr_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::ForwardEuler::add_visco_mass_contributions(
    CORE::LINALG::SparseOperator& jac) const
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> stiff_ptr = global_state().ExtractDisplBlock(jac);
  // set mass matrix
  stiff_ptr->Add(*global_state().GetMassMatrix(), false, 1.0, 0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::ForwardEuler::write_restart(
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
void STR::EXPLICIT::ForwardEuler::read_restart(IO::DiscretizationReader& ioreader)
{
  check_init_setup();
  ioreader.ReadVector(finertian_ptr_, "finert");
  ioreader.ReadVector(fviscon_ptr_, "fvisco");

  ModelEval().read_restart(ioreader);
  update_constant_state_contributions();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::ForwardEuler::UpdateStepState()
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

FOUR_C_NAMESPACE_CLOSE

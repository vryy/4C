/*-----------------------------------------------------------*/
/*! \file

\brief Adams-Bashforth-2 time integration for solid dynamics


\level 3

*/
/*-----------------------------------------------------------*/

#include "baci_structure_new_expl_ab2.H"

#include "baci_io.H"
#include "baci_lib_globalproblem.H"
#include "baci_linalg_utils_sparse_algebra_assemble.H"
#include "baci_structure_new_model_evaluator.H"
#include "baci_structure_new_timint_base.H"
#include "baci_structure_new_timint_basedataglobalstate.H"

BACI_NAMESPACE_OPEN

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
  CheckInit();

  // Call the Setup() of the abstract base class first.
  Generic::Setup();

  // ---------------------------------------------------------------------------
  // setup pointers to the force vectors of the global state data container
  // ---------------------------------------------------------------------------
  finertian_ptr_ = GlobalState().GetFinertialN();
  finertianp_ptr_ = GlobalState().GetFinertialNp();

  fviscon_ptr_ = GlobalState().GetFviscoN();
  fvisconp_ptr_ = GlobalState().GetFviscoNp();

  // ---------------------------------------------------------------------------
  // resizing of multi-step quantities
  // ---------------------------------------------------------------------------
  GlobalState().GetMultiTime()->Resize(-1, 0, true);
  GlobalState().GetDeltaTime()->Resize(-1, 0, true);
  GlobalState().GetMultiDis()->Resize(-1, 0, GlobalState().DofRowMapView(), true);
  GlobalState().GetMultiVel()->Resize(-1, 0, GlobalState().DofRowMapView(), true);
  GlobalState().GetMultiAcc()->Resize(-1, 0, GlobalState().DofRowMapView(), true);

  // here we initialized the dt of previous steps in the database, since a resize is performed
  const double dt = (*GlobalState().GetDeltaTime())[0];
  GlobalState().GetDeltaTime()->UpdateSteps(dt);

  // -------------------------------------------------------------------
  // set initial displacement
  // -------------------------------------------------------------------
  SetInitialDisplacement(
      TimInt().GetDataSDyn().GetInitialDisp(), TimInt().GetDataSDyn().StartFuncNo());

  // Has to be set before the PostSetup() routine is called!
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::AdamsBashforth2::PostSetup()
{
  CheckInitSetup();
  EquilibrateInitialState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::AdamsBashforth2::SetState(const Epetra_Vector& x)
{
  CheckInitSetup();

  const double dt = (*GlobalState().GetDeltaTime())[0];
  const double dto = (*GlobalState().GetDeltaTime())[-1];
  const double dta = (2.0 * dt * dto + dt * dt) / (2.0 * dto);
  const double dtb = -(dt * dt) / (2.0 * dto);

  // ---------------------------------------------------------------------------
  // new end-point acceleration
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> accnp_ptr = GlobalState().ExtractDisplEntries(x);
  GlobalState().GetAccNp()->Scale(1.0, *accnp_ptr);

  // ---------------------------------------------------------------------------
  // new end-point velocities
  // ---------------------------------------------------------------------------
  GlobalState().GetVelNp()->Update(1.0, (*(GlobalState().GetMultiVel()))[0], 0.0);
  GlobalState().GetVelNp()->Update(
      dta, (*(GlobalState().GetMultiAcc()))[0], dtb, (*(GlobalState().GetMultiAcc()))[-1], 1.0);

  // ---------------------------------------------------------------------------
  // new end-point displacements
  // ---------------------------------------------------------------------------
  GlobalState().GetDisNp()->Update(1.0, (*(GlobalState().GetMultiDis()))[0], 0.0);
  GlobalState().GetDisNp()->Update(
      dta, (*(GlobalState().GetMultiVel()))[0], dtb, (*(GlobalState().GetMultiVel()))[-1], 1.0);

  // ---------------------------------------------------------------------------
  // update the elemental state
  // ---------------------------------------------------------------------------
  ModelEval().UpdateResidual();
  ModelEval().RunRecover();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::AdamsBashforth2::AddViscoMassContributions(Epetra_Vector& f) const
{
  // viscous damping forces at t_{n+1}
  CORE::LINALG::AssembleMyVector(1.0, f, 1.0, *fvisconp_ptr_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::AdamsBashforth2::AddViscoMassContributions(
    CORE::LINALG::SparseOperator& jac) const
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> stiff_ptr = GlobalState().ExtractDisplBlock(jac);
  // set mass matrix
  stiff_ptr->Add(*GlobalState().GetMassMatrix(), false, 1.0, 0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::AdamsBashforth2::WriteRestart(
    IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  CheckInitSetup();
  // write dynamic forces
  iowriter.WriteVector("finert", finertian_ptr_);
  iowriter.WriteVector("fvisco", fviscon_ptr_);

  ModelEval().WriteRestart(iowriter, forced_writerestart);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::AdamsBashforth2::ReadRestart(IO::DiscretizationReader& ioreader)
{
  CheckInitSetup();
  ioreader.ReadVector(finertian_ptr_, "finert");
  ioreader.ReadVector(fviscon_ptr_, "fvisco");

  ModelEval().ReadRestart(ioreader);
  UpdateConstantStateContributions();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::AdamsBashforth2::UpdateStepState()
{
  CheckInitSetup();

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
double STR::EXPLICIT::AdamsBashforth2::MethodLinErrCoeffDis() const
{
  const double dt = (*GlobalState().GetDeltaTime())[0];
  const double dto = (*GlobalState().GetDeltaTime())[-1];
  return (2. * dt + 3. * dto) / (12. * dt);
}

BACI_NAMESPACE_CLOSE

/*-----------------------------------------------------------*/
/*! \file

\brief High order Adams-Bashforth time integration for solid dynamics


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_expl_abx.hpp"

#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_structure_new_model_evaluator.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <int TOrder>
STR::EXPLICIT::AdamsBashforthX<TOrder>::AdamsBashforthX()
    : fvisconp_ptr_(Teuchos::null),
      fviscon_ptr_(Teuchos::null),
      finertianp_ptr_(Teuchos::null),
      finertian_ptr_(Teuchos::null),
      compute_phase_(0)
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <int TOrder>
void STR::EXPLICIT::AdamsBashforthX<TOrder>::Setup()
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
  constexpr int nhist = TOrder - 1;
  GlobalState().GetMultiTime()->Resize(-nhist, 0, true);
  GlobalState().GetDeltaTime()->Resize(-nhist, 0, true);
  GlobalState().GetMultiDis()->Resize(-nhist, 0, GlobalState().DofRowMapView(), true);
  GlobalState().GetMultiVel()->Resize(-nhist, 0, GlobalState().DofRowMapView(), true);
  GlobalState().GetMultiAcc()->Resize(-nhist, 0, GlobalState().DofRowMapView(), true);

  // here we initialized the dt of previous steps in the database, since a resize is performed
  const double dt = (*GlobalState().GetDeltaTime())[0];
  for (int i = 0; i < nhist; ++i) GlobalState().GetDeltaTime()->UpdateSteps(dt);

  // -------------------------------------------------------------------
  // set initial displacement
  // -------------------------------------------------------------------
  set_initial_displacement(
      TimInt().GetDataSDyn().GetInitialDisp(), TimInt().GetDataSDyn().StartFuncNo());

  // Has to be set before the PostSetup() routine is called!
  issetup_ = true;

  // Set the compute phase to compute the initial value
  compute_phase_ = 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <int TOrder>
void STR::EXPLICIT::AdamsBashforthX<TOrder>::PostSetup()
{
  CheckInitSetup();
  equilibrate_initial_state();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <int TOrder>
void STR::EXPLICIT::AdamsBashforthX<TOrder>::SetState(const Epetra_Vector& x)
{
  CheckInitSetup();

  // ---------------------------------------------------------------------------
  // new end-point acceleration
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> accnp_ptr = GlobalState().ExtractDisplEntries(x);
  GlobalState().GetAccNp()->Scale(1.0, *accnp_ptr);
  if (compute_phase_ < TOrder)
  {
    const double dt = (*GlobalState().GetDeltaTime())[0];

    // ---------------------------------------------------------------------------
    // new end-point velocities
    // ---------------------------------------------------------------------------
    GlobalState().GetVelNp()->Update(1.0, *GlobalState().GetVelN(), 0.0);
    GlobalState().GetVelNp()->Update(dt, *GlobalState().GetAccN(), 1.0);

    // ---------------------------------------------------------------------------
    // new end-point displacements
    // ---------------------------------------------------------------------------
    GlobalState().GetDisNp()->Update(1.0, *GlobalState().GetDisN(), 0.0);
    GlobalState().GetDisNp()->Update(dt, *GlobalState().GetVelNp(), 1.0);
  }
  else
  {
    constexpr int nhist = TOrder - 1;

    const double dt = (*GlobalState().GetDeltaTime())[0];

    // At present, a variable step size for high order Adams-Bashforth is not supported due to a
    // good reference is not yet been found. The time coefficient shall be adapted for a variable
    // step size approach.
    double test = 0.0, dti = dt;
    for (int i = 0; i < nhist; ++i)
    {
      const double dti1 = (*GlobalState().GetDeltaTime())[-i - 1];
      test += std::abs(dti - dti1);
      dti = dti1;
    }

    if (test > 1.0e-13)
      FOUR_C_THROW("High Order AdamsBashforth does not currently support the variable step size.");

    // ---------------------------------------------------------------------------
    // new end-point velocities
    // ---------------------------------------------------------------------------
    GlobalState().GetVelNp()->Update(1.0, (*(GlobalState().GetMultiVel()))[0], 0.0);
    for (int i = 0; i < TOrder; ++i)
    {
      double c = AdamsBashforthHelper<TOrder>::exc[i];
      GlobalState().GetVelNp()->Update(c * dt, (*(GlobalState().GetMultiAcc()))[-i], 1.0);
    }

    // ---------------------------------------------------------------------------
    // new end-point displacements
    // ---------------------------------------------------------------------------
    GlobalState().GetDisNp()->Update(1.0, (*(GlobalState().GetMultiDis()))[0], 0.0);
    for (int i = 0; i < TOrder; ++i)
    {
      double c = AdamsBashforthHelper<TOrder>::exc[i];
      GlobalState().GetDisNp()->Update(c * dt, (*(GlobalState().GetMultiVel()))[-i], 1.0);
    }
  }

  // ---------------------------------------------------------------------------
  // update the elemental state
  // ---------------------------------------------------------------------------
  ModelEval().UpdateResidual();
  ModelEval().RunRecover();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <int TOrder>
void STR::EXPLICIT::AdamsBashforthX<TOrder>::add_visco_mass_contributions(Epetra_Vector& f) const
{
  // viscous damping forces at t_{n+1}
  CORE::LINALG::AssembleMyVector(1.0, f, 1.0, *fvisconp_ptr_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <int TOrder>
void STR::EXPLICIT::AdamsBashforthX<TOrder>::add_visco_mass_contributions(
    CORE::LINALG::SparseOperator& jac) const
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> stiff_ptr = GlobalState().ExtractDisplBlock(jac);
  // set mass matrix
  stiff_ptr->Add(*GlobalState().GetMassMatrix(), false, 1.0, 0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <int TOrder>
void STR::EXPLICIT::AdamsBashforthX<TOrder>::WriteRestart(
    IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  CheckInitSetup();
  // write dynamic forces
  iowriter.WriteVector("finert", finertian_ptr_);
  iowriter.WriteVector("fvisco", fviscon_ptr_);

  // write compute phase
  iowriter.WriteInt("compute_phase", compute_phase_);

  // write velocities and accelerations
  if (compute_phase_ >= TOrder)
  {
    for (int i = 0; i < TOrder; ++i)
    {
      std::stringstream velname;
      velname << "histvel_" << i;
      Teuchos::RCP<const Epetra_Vector> vel_ptr_ =
          Teuchos::rcpFromRef<const Epetra_Vector>((*(GlobalState().GetMultiVel()))[-i]);
      iowriter.WriteVector(velname.str(), vel_ptr_);

      std::stringstream accname;
      accname << "histacc_" << i;
      Teuchos::RCP<const Epetra_Vector> acc_ptr_ =
          Teuchos::rcpFromRef<const Epetra_Vector>((*(GlobalState().GetMultiAcc()))[-i]);
      iowriter.WriteVector(accname.str(), acc_ptr_);
    }
  }

  ModelEval().WriteRestart(iowriter, forced_writerestart);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <int TOrder>
void STR::EXPLICIT::AdamsBashforthX<TOrder>::ReadRestart(IO::DiscretizationReader& ioreader)
{
  CheckInitSetup();
  // read dynamic forces
  ioreader.ReadVector(finertian_ptr_, "finert");
  ioreader.ReadVector(fviscon_ptr_, "fvisco");

  // read compute phase
  if (ioreader.HasInt("compute_phase"))
  {
    compute_phase_ = ioreader.ReadInt("compute_phase");
  }
  else
  {
    compute_phase_ = 0;
  }

  // read velocities and accelerations
  if (compute_phase_ >= TOrder)
  {
    for (int i = TOrder - 1; i >= 0; --i)
    {
      std::stringstream velname;
      velname << "histvel_" << i;
      Teuchos::RCP<Epetra_Vector> vel_ptr =
          Teuchos::rcp(new Epetra_Vector(*GlobalState().GetVelN()));
      ioreader.ReadVector(vel_ptr, velname.str());
      GlobalState().GetMultiVel()->UpdateSteps(*vel_ptr);

      std::stringstream accname;
      accname << "histacc_" << i;
      Teuchos::RCP<Epetra_Vector> acc_ptr =
          Teuchos::rcp(new Epetra_Vector(*GlobalState().GetAccN()));
      ioreader.ReadVector(acc_ptr, accname.str());
      GlobalState().GetMultiAcc()->UpdateSteps(*acc_ptr);
    }
  }

  ModelEval().ReadRestart(ioreader);
  update_constant_state_contributions();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <int TOrder>
void STR::EXPLICIT::AdamsBashforthX<TOrder>::UpdateStepState()
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

  // update the compute phase step flag
  if (compute_phase_ < TOrder) ++compute_phase_;
}

/*----------------------------------------------------------------------------*
 | Template instantiation for supported order                                 |
 *----------------------------------------------------------------------------*/
template class STR::EXPLICIT::AdamsBashforthX<2>;
// template class STR::EXPLICIT::AdamsBashforthX<3>;
template class STR::EXPLICIT::AdamsBashforthX<4>;

FOUR_C_NAMESPACE_CLOSE

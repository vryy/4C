/*-----------------------------------------------------------*/
/*! \file

\brief Central differences time integration for solid dynamics


\level 3

*/
/*-----------------------------------------------------------*/

#include "baci_structure_new_expl_centrdiff.H"

#include "baci_io.H"
#include "baci_lib_globalproblem.H"
#include "baci_linalg_utils_sparse_algebra_assemble.H"
#include "baci_structure_new_model_evaluator.H"
#include "baci_structure_new_model_evaluator_generic.H"
#include "baci_structure_new_timint_base.H"
#include "baci_structure_new_timint_basedataglobalstate.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::EXPLICIT::CentrDiff::CentrDiff()
    : fvisconp_ptr_(Teuchos::null),
      fviscon_ptr_(Teuchos::null),
      finertianp_ptr_(Teuchos::null),
      finertian_ptr_(Teuchos::null)
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::CentrDiff::Setup()
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

  // -------------------------------------------------------------------
  // set initial displacement
  // -------------------------------------------------------------------
  const Teuchos::ParameterList& params = DRT::Problem::Instance()->StructuralDynamicParams();
  SetInitialDisplacement(DRT::INPUT::IntegralValue<INPAR::STR::InitialDisp>(params, "INITIALDISP"),
      params.get<int>("STARTFUNCNO"));

  // Has to be set before the PostSetup() routine is called!
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::CentrDiff::PostSetup()
{
  CheckInitSetup();
  EquilibrateInitialState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::CentrDiff::SetState(const Epetra_Vector& x)
{
  CheckInitSetup();

  const double dt = (*GlobalState().GetDeltaTime())[0];
  const double dthalf = dt / 2.0;

  ModelEval().ResetStepState();

  // ---------------------------------------------------------------------------
  // new end-point acceleration
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> accnp_ptr = GlobalState().ExtractDisplEntries(x);
  GlobalState().GetAccNp()->Scale(1.0, *accnp_ptr);

  // ---------------------------------------------------------------------------
  // new half-point velocities
  // ---------------------------------------------------------------------------
  GlobalState().GetVelNp()->Update(1.0, *GlobalState().GetVelN(), 0.0);
  GlobalState().GetVelNp()->Update(dthalf, *GlobalState().GetAccN(), 1.0);

  // ---------------------------------------------------------------------------
  // new end-point displacements
  // ---------------------------------------------------------------------------
  GlobalState().GetDisNp()->Update(1.0, *GlobalState().GetDisN(), 0.0);
  GlobalState().GetDisNp()->Update(dt, *GlobalState().GetVelNp(), 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::CentrDiff::AddViscoMassContributions(Epetra_Vector& f) const
{
  // viscous damping forces at t_{n+1}
  CORE::LINALG::AssembleMyVector(1.0, f, 1.0, *fvisconp_ptr_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::CentrDiff::AddViscoMassContributions(CORE::LINALG::SparseOperator& jac) const
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> stiff_ptr = GlobalState().ExtractDisplBlock(jac);
  // set mass matrix
  stiff_ptr->Add(*GlobalState().GetMassMatrix(), false, 1.0, 0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::CentrDiff::WriteRestart(
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
void STR::EXPLICIT::CentrDiff::ReadRestart(IO::DiscretizationReader& ioreader)
{
  CheckInitSetup();
  ioreader.ReadVector(finertian_ptr_, "finert");
  ioreader.ReadVector(fviscon_ptr_, "fvisco");

  ModelEval().ReadRestart(ioreader);
  UpdateConstantStateContributions();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::EXPLICIT::CentrDiff::UpdateStepState()
{
  CheckInitSetup();

  const double dt = (*GlobalState().GetDeltaTime())[0];
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
  GlobalState().GetVelNp()->Update(1.0, *GlobalState().GetVelN(), 0.0);
  GlobalState().GetVelNp()->Update(dthalf, *GlobalState().GetAccN(), 1.0);
  GlobalState().GetVelNp()->Update(dthalf, *GlobalState().GetAccNp(), 1.0);

  // ---------------------------------------------------------------------------
  // update model specific variables
  // ---------------------------------------------------------------------------
  ModelEval().UpdateStepState(0.0);
}

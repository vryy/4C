/*---------------------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all constraint terms


\level 3

*/
/*---------------------------------------------------------------------*/

#include "baci_structure_new_model_evaluator_lagpenconstraint.H"

#include "baci_constraint_lagpenconstraint_noxinterface.H"
#include "baci_io.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_linalg_sparsematrix.H"
#include "baci_linalg_sparseoperator.H"
#include "baci_linalg_utils_sparse_algebra_assemble.H"
#include "baci_linalg_utils_sparse_algebra_manipulation.H"
#include "baci_structure_new_timint_base.H"
#include "baci_structure_new_utils.H"
#include "baci_utils_exceptions.H"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::LagPenConstraint::LagPenConstraint()
    : disnp_ptr_(Teuchos::null),
      stiff_constr_ptr_(Teuchos::null),
      fstrconstr_np_ptr_(Teuchos::null),
      noxinterface_ptr_(Teuchos::null),
      noxinterface_prec_ptr_(Teuchos::null)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::Setup()
{
  CheckInit();

  // build the NOX::NLN::CONSTRAINT::Interface::Required object
  noxinterface_ptr_ = Teuchos::rcp(new LAGPENCONSTRAINT::NoxInterface());
  noxinterface_ptr_->Init(GStatePtr());
  noxinterface_ptr_->Setup();

  // build the NOX::NLN::CONSTRAINT::Interface::Preconditioner object
  noxinterface_prec_ptr_ = Teuchos::rcp(new LAGPENCONSTRAINT::NoxInterfacePrec());
  noxinterface_prec_ptr_->Init(GStatePtr());
  noxinterface_prec_ptr_->Setup();

  Teuchos::RCP<DRT::Discretization> dis = DiscretPtr();

  // setup the displacement pointer
  disnp_ptr_ = GState().GetDisNp();

  // contributions of constraints to structural rhs and stiffness
  fstrconstr_np_ptr_ = Teuchos::rcp(new Epetra_Vector(*GState().DofRowMapView()));
  stiff_constr_ptr_ =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*GState().DofRowMapView(), 81, true, true));

  // ToDo: we do not want to hand in the structural dynamics parameter list
  // to the manager in the future! -> get rid of it as soon as old
  // time-integration dies ...
  // initialize constraint manager
  constrman_ = Teuchos::rcp(new UTILS::ConstrManager());
  constrman_->Init(dis, DRT::Problem::Instance()->StructuralDynamicParams());
  constrman_->Setup(disnp_ptr_, DRT::Problem::Instance()->StructuralDynamicParams());

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::Reset(const Epetra_Vector& x)
{
  CheckInitSetup();

  // update the structural displacement vector
  disnp_ptr_ = GState().GetDisNp();

  fstrconstr_np_ptr_->PutScalar(0.0);
  stiff_constr_ptr_->Zero();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::LagPenConstraint::EvaluateForce()
{
  CheckInitSetup();

  double time_np = GState().GetTimeNp();
  Teuchos::ParameterList pcon;  // empty parameter list
  Teuchos::RCP<const Epetra_Vector> disn = GState().GetDisN();

  // only forces are evaluated!
  constrman_->EvaluateForceStiff(
      time_np, disn, disnp_ptr_, fstrconstr_np_ptr_, Teuchos::null, pcon);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::LagPenConstraint::EvaluateStiff()
{
  CheckInitSetup();

  double time_np = GState().GetTimeNp();
  Teuchos::ParameterList pcon;  // empty parameter list
  Teuchos::RCP<const Epetra_Vector> disn = GState().GetDisN();

  // only stiffnesses are evaluated!
  constrman_->EvaluateForceStiff(time_np, disn, disnp_ptr_, Teuchos::null, stiff_constr_ptr_, pcon);

  if (not stiff_constr_ptr_->Filled()) stiff_constr_ptr_->Complete();

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::LagPenConstraint::EvaluateForceStiff()
{
  CheckInitSetup();

  double time_np = GState().GetTimeNp();
  Teuchos::ParameterList pcon;  // empty parameter list
  Teuchos::RCP<const Epetra_Vector> disn = GState().GetDisN();

  constrman_->EvaluateForceStiff(
      time_np, disn, disnp_ptr_, fstrconstr_np_ptr_, stiff_constr_ptr_, pcon);

  if (not stiff_constr_ptr_->Filled()) stiff_constr_ptr_->Complete();

  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::LagPenConstraint::AssembleForce(
    Epetra_Vector& f, const double& timefac_np) const
{
  Teuchos::RCP<const Epetra_Vector> block_vec_ptr = Teuchos::null;

  CORE::LINALG::AssembleMyVector(1.0, f, timefac_np, *fstrconstr_np_ptr_);

  if (noxinterface_prec_ptr_->IsSaddlePointSystem())
  {
    // assemble constraint rhs
    block_vec_ptr = constrman_->GetError();

    if (block_vec_ptr.is_null())
      dserror(
          "The constraint model vector is a nullptr pointer, although \n"
          "the structural part indicates, that constraint contributions \n"
          "are present!");

    const int elements_f = f.Map().NumGlobalElements();
    const int max_gid = GetBlockDofRowMapPtr()->MaxAllGID();
    // only call when f is the rhs of the full problem (not for structural
    // equilibriate initial state call)
    if (elements_f == max_gid + 1)
      CORE::LINALG::AssembleMyVector(1.0, f, timefac_np, *block_vec_ptr);
  }

  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::LagPenConstraint::AssembleJacobian(
    CORE::LINALG::SparseOperator& jac, const double& timefac_np) const
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> block_ptr = Teuchos::null;

  // --- Kdd - block ---------------------------------------------------
  Teuchos::RCP<CORE::LINALG::SparseMatrix> jac_dd_ptr = GState().ExtractDisplBlock(jac);
  jac_dd_ptr->Add(*stiff_constr_ptr_, false, timefac_np, 1.0);
  // no need to keep it
  stiff_constr_ptr_->Zero();

  if (noxinterface_prec_ptr_->IsSaddlePointSystem())
  {
    // --- Kdz - block - scale with time-integrator dependent value!-----
    block_ptr = (Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(
        constrman_->GetConstrMatrix(), true));
    block_ptr->Scale(timefac_np);
    GState().AssignModelBlock(jac, *block_ptr, Type(), DRT::UTILS::MatBlockType::displ_lm);
    // reset the block pointer, just to be on the safe side
    block_ptr = Teuchos::null;

    // --- Kzd - block - no scaling of this block (cf. diss Kloeppel p78)
    block_ptr =
        (Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(constrman_->GetConstrMatrix(), true))
            ->Transpose();
    GState().AssignModelBlock(jac, *block_ptr, Type(), DRT::UTILS::MatBlockType::lm_displ);
    // reset the block pointer, just to be on the safe side
    block_ptr = Teuchos::null;
  }

  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::WriteRestart(
    IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  iowriter.WriteVector("lagrmultiplier", constrman_->GetLagrMultVector());
  iowriter.WriteVector("refconval", constrman_->GetRefBaseValues());

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::ReadRestart(IO::DiscretizationReader& ioreader)
{
  double time_n = GState().GetTimeN();
  constrman_->ReadRestart(ioreader, time_n);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::RunPostComputeX(
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  CheckInitSetup();

  Teuchos::RCP<Epetra_Vector> lagmult_incr =
      Teuchos::rcp(new Epetra_Vector(*GetBlockDofRowMapPtr()));

  CORE::LINALG::Export(dir, *lagmult_incr);

  constrman_->UpdateLagrMult(lagmult_incr);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::UpdateStepState(const double& timefac_n)
{
  constrman_->Update();

  // add the constraint force contributions to the old structural
  // residual state vector
  if (not fstrconstr_np_ptr_.is_null())
  {
    Teuchos::RCP<Epetra_Vector>& fstructold_ptr = GState().GetFstructureOld();
    fstructold_ptr->Update(timefac_n, *fstrconstr_np_ptr_, 1.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::UpdateStepElement()
{
  // empty
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::DetermineStressStrain()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::DetermineEnergy()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::DetermineOptionalQuantity()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::OutputStepState(
    IO::DiscretizationWriter& iowriter) const
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::ResetStepState()
{
  CheckInitSetup();

  dserror("Not yet implemented");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Teuchos::RCP<LAGPENCONSTRAINT::NoxInterface>&
STR::MODELEVALUATOR::LagPenConstraint::NoxInterfacePtr()
{
  CheckInitSetup();

  return noxinterface_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Teuchos::RCP<LAGPENCONSTRAINT::NoxInterfacePrec>&
STR::MODELEVALUATOR::LagPenConstraint::NoxInterfacePrecPtr()
{
  CheckInitSetup();

  return noxinterface_prec_ptr_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::LagPenConstraint::GetBlockDofRowMapPtr() const
{
  CheckInitSetup();

  if (noxinterface_prec_ptr_->IsSaddlePointSystem())
  {
    return constrman_->GetConstraintMap();
  }
  else
  {
    return GState().DofRowMap();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::LagPenConstraint::GetCurrentSolutionPtr()
    const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
STR::MODELEVALUATOR::LagPenConstraint::GetLastTimeStepSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::PostOutput()
{
  CheckInitSetup();
  // empty

  return;
}  // PostOutput()

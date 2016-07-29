/*---------------------------------------------------------------------*/
/*!
\file str_model_evaluator_lagpenconstraint.cpp

\brief Evaluation and assembly of all constraint terms

\maintainer Marc Hirschvogel

\date Jun 29, 2016

\level 3

*/
/*---------------------------------------------------------------------*/

#include "str_model_evaluator_lagpenconstraint.H"

#include "str_timint_base.H"
#include "str_utils.H"

#include <Epetra_Vector.h>
#include <Epetra_Time.h>
#include <Teuchos_ParameterList.hpp>

#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_io/io.H"

#include "../drt_constraint/lagpenconstraint_noxinterface.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::LagPenConstraint::LagPenConstraint()
    : disnp_ptr_(Teuchos::null),
      stiff_constr_ptr_(Teuchos::null),
      fstrconstr_np_ptr_(Teuchos::null)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::Setup()
{

  CheckInit();

  Teuchos::RCP<DRT::Discretization> dis = GState().GetMutableDiscret();

  // setup the displacement pointer
  disnp_ptr_ = GState().GetMutableDisNp();
  fstrconstr_np_ptr_ =
      Teuchos::rcp(new Epetra_Vector(*GState().DofRowMapView()));
  stiff_constr_ptr_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *GState().DofRowMapView(), 81, true, true));

  // initialize constraint manager
  constrman_ = Teuchos::rcp(new UTILS::ConstrManager(dis,
      disnp_ptr_,
      DRT::Problem::Instance()->StructuralDynamicParams()));

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::Reset(const Epetra_Vector& x)
{
  CheckInitSetup();

  fstrconstr_np_ptr_->Scale(0.0);
  stiff_constr_ptr_->Zero();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::LagPenConstraint::EvaluateForce()
{
  CheckInitSetup();

  double time_np = GState().GetTimeNp();
  Teuchos::ParameterList pcon;
  pcon.set("total time", time_np);

  Teuchos::RCP<const Epetra_Vector> disn = GState().GetDisN();
//  Teuchos::RCP<const Epetra_Vector> disnp = GState().GetDisNp();

  constrman_->EvaluateForceStiff(time_np, disn, disnp_ptr_,
      fstrconstr_np_ptr_, Teuchos::null, pcon);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::LagPenConstraint::EvaluateStiff()
{

  CheckInitSetup();

  double time_np = GState().GetTimeNp();
  Teuchos::ParameterList pcon;
  pcon.set("total time", time_np);

  Teuchos::RCP<const Epetra_Vector> disn = GState().GetDisN();
//  Teuchos::RCP<const Epetra_Vector> disnp = GState().GetDisNp();

  constrman_->EvaluateForceStiff(time_np, disn, disnp_ptr_, Teuchos::null,
      stiff_constr_ptr_, pcon);

  if (not stiff_constr_ptr_->Filled())
    stiff_constr_ptr_->Complete();

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::LagPenConstraint::EvaluateForceStiff()
{

  CheckInitSetup();

  double time_np = GState().GetTimeNp();
  Teuchos::ParameterList pcon;
  pcon.set("total time", time_np);

  Teuchos::RCP<const Epetra_Vector> disn = GState().GetDisN();
//  Teuchos::RCP<const Epetra_Vector> disnp = GState().GetDisNp();

  constrman_->EvaluateForceStiff(time_np, disn,
      disnp_ptr_, fstrconstr_np_ptr_, stiff_constr_ptr_, pcon);

  if (not stiff_constr_ptr_->Filled())
    stiff_constr_ptr_->Complete();

  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::LagPenConstraint::AssembleForce(Epetra_Vector& f,
    const double & timefac_np) const
{
  if (timefac_np!=1.0)
    dserror("You have to consider the corresponding time integration"
      " factors.");

  Teuchos::RCP<const Epetra_Vector> block_vec_ptr = Teuchos::null;

  // assemble 0D model rhs
  block_vec_ptr = constrman_->GetError();

  if (block_vec_ptr.is_null())
    dserror("The constraint model vector is a NULL pointer, although \n"
        "the structural part indicates, that constraint contributions \n"
        "are present!");

  STR::AssembleVector(1.0,f,1.0,*fstrconstr_np_ptr_);
  STR::AssembleVector(1.0,f,1.0,*block_vec_ptr);

  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::LagPenConstraint::AssembleJacobian(
    LINALG::SparseOperator& jac,
    const double & timefac_np) const
{
  if (timefac_np!=1.0)
    dserror("You have to consider the corresponding time integration"
      " factors.");

  Teuchos::RCP<const LINALG::SparseMatrix> block_ptr = Teuchos::null;

  // --- Kdd - block ---------------------------------------------------
  Teuchos::RCP<LINALG::SparseMatrix> jac_dd_ptr =
      GState().ExtractDisplBlock(jac);
  jac_dd_ptr->Add(*stiff_constr_ptr_,false,timefac_np,1.0);
  // no need to keep it
  stiff_constr_ptr_->Zero();
  // --- Kdz - block ---------------------------------------------------
  block_ptr = (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(
      constrman_->GetConstrMatrix(),true));

  GState().AssignModelBlock(jac,*block_ptr,Type(),
      STR::block_displ_lm);
  // reset the block pointer, just to be on the safe side
  block_ptr = Teuchos::null;
  // --- Kzd - block ---------------------------------------------------
  block_ptr = (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(
      constrman_->GetConstrMatrix()))->Transpose();

  GState().AssignModelBlock(jac,*block_ptr,Type(),
      STR::block_lm_displ);
  // reset the block pointer, just to be on the safe side
  block_ptr = Teuchos::null;

  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::WriteRestart(
        IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const
{

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::ReadRestart(
    IO::DiscretizationReader& ioreader)
{

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::RecoverState(
    const Epetra_Vector& xold,
    const Epetra_Vector& dir,
    const Epetra_Vector& xnew)
{
  Teuchos::RCP<Epetra_Vector> lagmult_incr =
      Teuchos::rcp(new Epetra_Vector(*GetBlockDofRowMapPtr()));

  LINALG::Export(dir,*lagmult_incr);

  constrman_->UpdateLagrMult(lagmult_incr);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::UpdateStepState(
    const double& timefac_n)
{
  // empty
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

  // build the NOX::NLN::CONSTRAINT::Interface::Required object
  noxinterface_ptr_ = Teuchos::rcp(new LAGPENCONSTRAINT::NoxInterface());
  noxinterface_ptr_->Init(GStatePtr());
  noxinterface_ptr_->Setup();

  return noxinterface_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Teuchos::RCP<LAGPENCONSTRAINT::NoxInterfacePrec>&
    STR::MODELEVALUATOR::LagPenConstraint::NoxInterfacePtrPrec()
{
  CheckInitSetup();

  // build the NOX::NLN::CONSTRAINT::Interface::Preconditioner object
  noxinterface_ptr_prec_ = Teuchos::rcp(new LAGPENCONSTRAINT::NoxInterfacePrec());
  noxinterface_ptr_prec_->Init();
  noxinterface_ptr_prec_->Setup();

  return noxinterface_ptr_prec_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::LagPenConstraint::
    GetBlockDofRowMapPtr() const
{
  CheckInitSetup();
  return constrman_->GetConstraintMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::LagPenConstraint::
    GetCurrentSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::LagPenConstraint::
    GetLastTimeStepSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

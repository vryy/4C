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
      stiff_ptr_(Teuchos::null)
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

  // initialize constraint manager
  constrman_ = Teuchos::rcp(new UTILS::ConstrManager(dis,
                                                 disnp_ptr_,
                                                 DRT::Problem::Instance()->StructuralDynamicParams()));

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::LagPenConstraint::ApplyForce(
    const Epetra_Vector& x,
    Epetra_Vector& f)
{

  CheckInitSetup();
  Reset(x);

  // get structural displacement DOFs
  Teuchos::RCP<Epetra_Vector> block_vec_ptr =
        Teuchos::rcp(new Epetra_Vector(*GState().DofRowMap(),true));

  // dummy matrix
  Teuchos::RCP<LINALG::SparseMatrix> dummat = Teuchos::null;

  double time_np = GState().GetTimeNp();
  Teuchos::ParameterList pcon;
  pcon.set("total time", time_np);

  Teuchos::RCP<const Epetra_Vector> disn = GState().GetDisN();
//  Teuchos::RCP<const Epetra_Vector> disnp = GState().GetDisNp();

  constrman_->EvaluateForceStiff(time_np, disn, disnp_ptr_, block_vec_ptr, dummat, pcon);

  // assemble constraint contribution to structural rhs
  STR::AssembleVector(1.0,f,1.0,*block_vec_ptr);

  // reset the block pointer, just to be on the safe side
  block_vec_ptr = Teuchos::null;

  // --- assemble right-hand-side ---------------------------------------
  AssembleRhs(f);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::LagPenConstraint::ApplyStiff(
    const Epetra_Vector& x,
    LINALG::SparseOperator& jac)
{

  CheckInitSetup();
  Reset(x, jac);

  // dummy vector
  Teuchos::RCP<Epetra_Vector> dumvec = Teuchos::null;

  // get structural stiffness matrix
  Teuchos::RCP<LINALG::SparseMatrix> jac_dd =
      GState().ExtractModelBlock(jac,INPAR::STR::model_lag_pen_constraint,
          STR::block_displ_displ);

  double time_np = GState().GetTimeNp();
  Teuchos::ParameterList pcon;
  pcon.set("total time", time_np);

  Teuchos::RCP<const Epetra_Vector> disn = GState().GetDisN();
//  Teuchos::RCP<const Epetra_Vector> disnp = GState().GetDisNp();

  constrman_->EvaluateForceStiff(time_np, disn, disnp_ptr_, dumvec, jac_dd, pcon);

  // --- assemble jacobian matrix ---------------------------------------
  AssembleJacobian(jac);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::LagPenConstraint::ApplyForceStiff(
    const Epetra_Vector& x,
    Epetra_Vector& f,
    LINALG::SparseOperator& jac)
{

  CheckInitSetup();
  Reset(x, jac);

  // get structural displacement DOFs
  Teuchos::RCP<Epetra_Vector> block_vec_ptr =
        Teuchos::rcp(new Epetra_Vector(*GState().DofRowMap(),true));

  // get structural stiffness matrix
  Teuchos::RCP<LINALG::SparseMatrix> jac_dd =
      GState().ExtractModelBlock(jac,INPAR::STR::model_lag_pen_constraint,
          STR::block_displ_displ);

  double time_np = GState().GetTimeNp();
  Teuchos::ParameterList pcon;
  pcon.set("total time", time_np);

  Teuchos::RCP<const Epetra_Vector> disn = GState().GetDisN();
//  Teuchos::RCP<const Epetra_Vector> disnp = GState().GetDisNp();

  constrman_->EvaluateForceStiff(time_np, disn, disnp_ptr_, block_vec_ptr, jac_dd, pcon);

  // assemble constraint contribution to structural rhs
  STR::AssembleVector(1.0,f,1.0,*block_vec_ptr);

  // reset the block pointer, just to be on the safe side
  block_vec_ptr = Teuchos::null;

  // --- assemble right-hand-side ---------------------------------------
  AssembleRhs(f);
  // --- assemble jacobian matrix ---------------------------------------
  AssembleJacobian(jac);

  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::AssembleRhs(Epetra_Vector& f) const
{

  Teuchos::RCP<const Epetra_Vector> block_vec_ptr = Teuchos::null;

  // assemble 0D model rhs
  block_vec_ptr = constrman_->GetError();

  if (block_vec_ptr.is_null())
    dserror("The constraint model vector is a NULL pointer, although \n"
        "the structural part indicates, that constraint contributions \n"
        "are present!");

  STR::AssembleVector(1.0,f,1.0,*block_vec_ptr);

}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::AssembleJacobian(
    LINALG::SparseOperator& jac) const
{

  Teuchos::RCP<const LINALG::SparseMatrix> block_ptr = Teuchos::null;

  // --- Kdz - block ---------------------------------------------------
  block_ptr = (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(constrman_->GetConstrMatrix()));

  GState().AssignModelBlock(jac,*block_ptr,INPAR::STR::model_lag_pen_constraint,
      STR::block_displ_lm);
  // reset the block pointer, just to be on the safe side
  block_ptr = Teuchos::null;
  // --- Kzd - block ---------------------------------------------------
  block_ptr = (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(constrman_->GetConstrMatrix()))->Transpose();

  GState().AssignModelBlock(jac,*block_ptr,INPAR::STR::model_lag_pen_constraint,
      STR::block_lm_displ);
  // reset the block pointer, just to be on the safe side
  block_ptr = Teuchos::null;

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
void STR::MODELEVALUATOR::LagPenConstraint::UpdateStepState()
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
void STR::MODELEVALUATOR::LagPenConstraint::Reset(
    const Epetra_Vector& x,
    LINALG::SparseOperator& jac)
{
  CheckInitSetup();
  Reset(x);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::LagPenConstraint::Reset(const Epetra_Vector& x)
{
  CheckInitSetup();

  // update the structural displacement vector
  //disnp_ptr_ = GState().GetDisNp();

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

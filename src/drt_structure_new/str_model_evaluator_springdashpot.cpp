/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_springdashpot.cpp

\maintainer Martin Pfaller

\date Feb 29, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_model_evaluator_springdashpot.H"
#include "str_timint_base.H"
#include "str_utils.H"

#include <Epetra_Vector.h>
#include <Epetra_Time.h>
#include <Teuchos_ParameterList.hpp>

#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_io/io.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::SpringDashpot::SpringDashpot()
    : n_conds_(0),
      disnp_ptr_(Teuchos::null),
      stiff_ptr_(Teuchos::null)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::Setup()
{
  if (not IsInit())
    dserror("Init() has not been called, yet!");

  // get all spring dashpot conditions
  std::vector<Teuchos::RCP<DRT::Condition> > springdashpots;
  discret_ptr_->GetCondition("SpringDashpot",springdashpots);

  // number of spring dashpot conditions
  n_conds_ = (int)springdashpots.size();;

  // new instance of spring dashpot BC for each condition
  for (int i=0; i<n_conds_; ++i)
    springs_.push_back(Teuchos::rcp(new UTILS::SpringDashpotNew(discret_ptr_, springdashpots[i])));

  // setup the displacement pointer
  disnp_ptr_ = gstate_ptr_->GetMutableDisNp();

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::SpringDashpot::ApplyForce(
    const Epetra_Vector& x,
    Epetra_Vector& f)
{
  CheckInitSetup();
  Reset(x);

  // loop over all spring dashpot conditions and evaluate them
  for (int i=0; i<n_conds_; ++i)
    springs_[i]->EvaluateForce(f, disnp_ptr_);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::SpringDashpot::ApplyStiff(
    const Epetra_Vector& x,
    LINALG::SparseOperator& jac)
{
  CheckInitSetup();
  Reset(x, jac);

  // dummy vector
  Teuchos::RCP<Epetra_Vector> f_disp =
        Teuchos::rcp(new Epetra_Vector(*gstate_ptr_->DofRowMap(),true));

  // loop over all spring dashpot conditions and evaluate them
  for (int i=0; i<n_conds_; ++i)
    springs_[i]->EvaluateForceStiff(*stiff_ptr_, *f_disp, disnp_ptr_);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::SpringDashpot::ApplyForceStiff(
    const Epetra_Vector& x,
    Epetra_Vector& f,
    LINALG::SparseOperator& jac)
{
  CheckInitSetup();
  Reset(x, jac);

  // get displacement DOFs
  Teuchos::RCP<Epetra_Vector> f_disp =
      Teuchos::rcp(new Epetra_Vector(*gstate_ptr_->DofRowMap(),true));

  // loop over all spring dashpot conditions and evaluate them
  for (int i=0; i<n_conds_; ++i)
    springs_[i]->EvaluateForceStiff(*stiff_ptr_, *f_disp, disnp_ptr_);

  STR::AssembleForce(f,*f_disp,1.0);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::UpdateStepState()
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::UpdateStepElement()
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::OutputStepState()
{
  // row maps for export
  Teuchos::RCP<Epetra_Vector> gap = Teuchos::rcp(new Epetra_Vector(*(discret_ptr_->NodeRowMap()),true));
  Teuchos::RCP<Epetra_MultiVector> normals = Teuchos::rcp(new Epetra_MultiVector(*(discret_ptr_->NodeRowMap()),3,true));
  Teuchos::RCP<Epetra_MultiVector> springstress = Teuchos::rcp(new Epetra_MultiVector(*(discret_ptr_->NodeRowMap()),3,true));

  // collect outputs from all spring dashpot conditions
  bool found_cursurfnormal = false;
  for (int i=0; i<n_conds_; ++i)
  {
    // get spring type from current condition
    const UTILS::SpringDashpotNew::SpringType stype = springs_[i]->GetSpringType();
    if(stype == UTILS::SpringDashpotNew::cursurfnormal)
    {
      springs_[i]->OutputGapNormal(gap, normals, springstress);
      found_cursurfnormal = true;
    }
  }

  // write vectors to output
  if (found_cursurfnormal)
  {
    gio_ptr_->GetMutableOutputPtr()->WriteVector("gap", gap);
    gio_ptr_->GetMutableOutputPtr()->WriteVector("curnormals", normals);
    gio_ptr_->GetMutableOutputPtr()->WriteVector("springstress", springstress);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::Reset(
    const Epetra_Vector& x,
    LINALG::SparseOperator& jac)
{
  stiff_ptr_ = gstate_ptr_->ExtractDisplBlock(jac);

  CheckInitSetup();

  Reset(x);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::Reset(const Epetra_Vector& x)
{
  CheckInitSetup();

  // loop over all spring dashpot conditions and reset them
  for (int i=0; i<n_conds_; ++i)
    springs_[i]->ResetNewton();

  // update the structural displacement vector
  disnp_ptr_ = gstate_ptr_->GetDisNp();
}

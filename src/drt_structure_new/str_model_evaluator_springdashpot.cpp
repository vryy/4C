/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_springdashpot.cpp

\brief Evaluation and assembly of all spring dashpot terms

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
  Discret().GetCondition("SpringDashpot",springdashpots);

  // number of spring dashpot conditions
  n_conds_ = (int)springdashpots.size();;

  // new instance of spring dashpot BC for each condition
  for (int i=0; i<n_conds_; ++i)
    springs_.push_back(Teuchos::rcp(new UTILS::SpringDashpotNew(DiscretPtr(), springdashpots[i])));

  // setup the displacement pointer
  disnp_ptr_ = GState().GetMutableDisNp();

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
        Teuchos::rcp(new Epetra_Vector(*GState().DofRowMap(),true));

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
      Teuchos::rcp(new Epetra_Vector(*GState().DofRowMap(),true));

  // loop over all spring dashpot conditions and evaluate them
  for (int i=0; i<n_conds_; ++i)
    springs_[i]->EvaluateForceStiff(*stiff_ptr_, *f_disp, disnp_ptr_);

  STR::AssembleVector(1.0,f,1.0,*f_disp);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::WriteRestart(
        IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const
{
  // row maps for export
  Teuchos::RCP<Epetra_MultiVector> springoffsetprestr =
      Teuchos::rcp(new Epetra_MultiVector(*(Discret().NodeRowMap()),3,true));

  // collect outputs from all spring dashpot conditions
  for (int i=0; i<n_conds_; ++i)
    springs_[i]->OutputPrestrOffset(springoffsetprestr);

  // write vector to output for restart
  iowriter.WriteVector("springoffsetprestr", springoffsetprestr);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::ReadRestart(
    IO::DiscretizationReader& ioreader)
{
  Teuchos::RCP<Epetra_MultiVector> tempvec =
      Teuchos::rcp(new Epetra_MultiVector(*(Discret().NodeRowMap()),3,true));

  ioreader.ReadMultiVector(tempvec, "springoffsetprestr");
  // loop over all spring dashpot conditions and reset them
  for (int i=0; i<n_conds_; ++i)
    springs_[i]->SetRestart(tempvec);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::RecoverState(
    const Epetra_Vector& xold,
    const Epetra_Vector& dir,
    const Epetra_Vector& xnew)
{
  // empty
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
void STR::MODELEVALUATOR::SpringDashpot::DetermineStressStrain()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::DetermineEnergy()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::OutputStepState(
    IO::DiscretizationWriter& iowriter) const
{
  // row maps for export
  Teuchos::RCP<Epetra_Vector> gap =
      Teuchos::rcp(new Epetra_Vector(*(Discret().NodeRowMap()),true));
  Teuchos::RCP<Epetra_MultiVector> normals =
      Teuchos::rcp(new Epetra_MultiVector(*(Discret().NodeRowMap()),3,true));
  Teuchos::RCP<Epetra_MultiVector> springstress =
      Teuchos::rcp(new Epetra_MultiVector(*(Discret().NodeRowMap()),3,true));

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
    iowriter.WriteVector("gap", gap);
    iowriter.WriteVector("curnormals", normals);
    iowriter.WriteVector("springstress", springstress);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::Reset(
    const Epetra_Vector& x,
    LINALG::SparseOperator& jac)
{
  stiff_ptr_ = GState().ExtractDisplBlock(jac);

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
  disnp_ptr_ = GState().GetDisNp();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::SpringDashpot::
    GetBlockDofRowMapPtr() const
{
  CheckInitSetup();
  return GState().DofRowMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::SpringDashpot::
    GetCurrentSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::SpringDashpot::
    GetLastTimeStepSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

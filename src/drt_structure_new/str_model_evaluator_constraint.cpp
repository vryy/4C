/*---------------------------------------------------------------------*/
/*!
\file str_model_evaluator_constraint.cpp

\brief Evaluation and assembly of all constraint terms

\maintainer Marc Hirschvogel

\date Jun 29, 2016

\level 3

*/
/*---------------------------------------------------------------------*/

#include "str_model_evaluator_constraint.H"

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
STR::MODELEVALUATOR::Constraint::Constraint()
    : n_conds_(0),
      disnp_ptr_(Teuchos::null),
      stiff_ptr_(Teuchos::null)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraint::Setup()
{
  if (not IsInit())
    dserror("Init() has not been called, yet!");

  // get all constraint conditions
  std::vector<Teuchos::RCP<DRT::Condition> > constraints;
  Discret().GetCondition("Constraint",constraints);

  // number of spring dashpot conditions
  n_conds_ = (int)constraints.size();;



  // setup the displacement pointer
  disnp_ptr_ = GState().GetMutableDisNp();

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Constraint::ApplyForce(
    const Epetra_Vector& x,
    Epetra_Vector& f)
{
  CheckInitSetup();
  Reset(x);

//  // loop over all spring dashpot conditions and evaluate them
//  for (int i=0; i<n_conds_; ++i)
//    constraints_[i]->EvaluateForce(f, disnp_ptr_);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Constraint::ApplyStiff(
    const Epetra_Vector& x,
    LINALG::SparseOperator& jac)
{
  CheckInitSetup();
  Reset(x, jac);

  // dummy vector
  Teuchos::RCP<Epetra_Vector> f_disp =
        Teuchos::rcp(new Epetra_Vector(*GState().DofRowMap(),true));

  //StiffnessAndInternalForces(time, dis, disn, pwindk);



  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Constraint::ApplyForceStiff(
    const Epetra_Vector& x,
    Epetra_Vector& f,
    LINALG::SparseOperator& jac)
{
  CheckInitSetup();
  Reset(x, jac);

  // get displacement DOFs
  Teuchos::RCP<Epetra_Vector> f_disp =
      Teuchos::rcp(new Epetra_Vector(*GState().DofRowMap(),true));

  //StiffnessAndInternalForces(time, dis, disn, pwindk);

  STR::AssembleVector(1.0,f,1.0,*f_disp);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraint::WriteRestart(
        IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const
{
  // row maps for export
  Teuchos::RCP<Epetra_MultiVector> springoffsetprestr =
      Teuchos::rcp(new Epetra_MultiVector(*(Discret().NodeRowMap()),3,true));

//  // collect outputs from all spring dashpot conditions
//  for (int i=0; i<n_conds_; ++i)
//    constraints_[i]->OutputPrestrOffset(springoffsetprestr);

  // write vector to output for restart
  iowriter.WriteVector("springoffsetprestr", springoffsetprestr);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraint::ReadRestart(
    IO::DiscretizationReader& ioreader)
{
  Teuchos::RCP<Epetra_MultiVector> tempvec =
      Teuchos::rcp(new Epetra_MultiVector(*(Discret().NodeRowMap()),3,true));

  ioreader.ReadMultiVector(tempvec, "springoffsetprestr");
  // loop over all spring dashpot conditions and reset them
//  for (int i=0; i<n_conds_; ++i)
//    constraints_[i]->SetRestart(tempvec);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraint::RecoverState(
    const Epetra_Vector& xold,
    const Epetra_Vector& dir,
    const Epetra_Vector& xnew)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraint::UpdateStepState()
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraint::UpdateStepElement()
{
  // empty
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraint::DetermineStressStrain()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraint::DetermineEnergy()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraint::OutputStepState(
    IO::DiscretizationWriter& iowriter) const
{
  // row maps for export
  Teuchos::RCP<Epetra_Vector> gap =
      Teuchos::rcp(new Epetra_Vector(*(Discret().NodeRowMap()),true));
  Teuchos::RCP<Epetra_MultiVector> normals =
      Teuchos::rcp(new Epetra_MultiVector(*(Discret().NodeRowMap()),3,true));
  Teuchos::RCP<Epetra_MultiVector> constraintstress =
      Teuchos::rcp(new Epetra_MultiVector(*(Discret().NodeRowMap()),3,true));

//  // collect outputs from all spring dashpot conditions
//  bool found_cursurfnormal = false;
//  for (int i=0; i<n_conds_; ++i)
//  {
//    // get spring type from current condition
//    const UTILS::Constraint::SpringType stype = constraints_[i]->GetSpringType();
//    if(stype == UTILS::Constraint::cursurfnormal)
//    {
//      constraints_[i]->OutputGapNormal(gap, normals, constraintstress);
//      found_cursurfnormal = true;
//    }
//  }

//  // write vectors to output
//  if (found_cursurfnormal)
//  {
//    iowriter.WriteVector("gap", gap);
//    iowriter.WriteVector("curnormals", normals);
//    iowriter.WriteVector("constraintstress", constraintstress);
//  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraint::ResetStepState()
{
  CheckInitSetup();

  dserror("Not yet implemented");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraint::Reset(
    const Epetra_Vector& x,
    LINALG::SparseOperator& jac)
{
  stiff_ptr_ = GState().ExtractDisplBlock(jac);

  CheckInitSetup();

  Reset(x);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Constraint::Reset(const Epetra_Vector& x)
{
  CheckInitSetup();

//  // loop over all spring dashpot conditions and reset them
//  for (int i=0; i<n_conds_; ++i)
//    constraints_[i]->ResetNewton();

  // update the structural displacement vector
  disnp_ptr_ = GState().GetDisNp();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::Constraint::
    GetBlockDofRowMapPtr() const
{
  CheckInitSetup();
  return GState().DofRowMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Constraint::
    GetCurrentSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Constraint::
    GetLastTimeStepSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

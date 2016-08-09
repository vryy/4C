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

#include "../drt_structure_new/str_model_evaluator_data.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::SpringDashpot::SpringDashpot()
    : n_conds_(0),
      disnp_ptr_(Teuchos::null),
      velnp_ptr_(Teuchos::null),
      stiff_spring_ptr_(Teuchos::null),
      fspring_np_ptr_(Teuchos::null)
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
  n_conds_ = (int)springdashpots.size();

  // new instance of spring dashpot BC for each condition
  for (int i=0; i<n_conds_; ++i)
    springs_.push_back(Teuchos::rcp(new UTILS::SpringDashpotNew(DiscretPtr(), springdashpots[i])));

  // setup the displacement pointer
  disnp_ptr_ = GState().GetMutableDisNp();
  velnp_ptr_ = GState().GetMutableVelNp();
  stiff_spring_ptr_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *GState().DofRowMapView(), 81, true, true));

  // set flag
  issetup_ = true;
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

  // update the structural displacement vector
  velnp_ptr_ = GState().GetVelNp();

  // Zero out the stiffness contributions matrix
  stiff_spring_ptr_->Zero();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::SpringDashpot::EvaluateForce()
{
  CheckInitSetup();

  // loop over all spring dashpot conditions and evaluate them
  fspring_np_ptr_ = Teuchos::rcp(new Epetra_Vector(*GState().DofRowMapView()));
  for (int i=0; i<n_conds_; ++i)
    springs_[i]->EvaluateForce(*fspring_np_ptr_, disnp_ptr_, velnp_ptr_);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::SpringDashpot::EvaluateStiff()
{
  CheckInitSetup();

  fspring_np_ptr_ =
        Teuchos::rcp(new Epetra_Vector(*GState().DofRowMapView(),true));

  // factors from time-integrator for derivative of d(v_{n+1}) / d(d_{n+1})
  // needed for stiffness contribution from dashpot
  const double fac_vel = EvalData().GetTimIntFactorVel();
  const double fac_disp = EvalData().GetTimIntFactorDisp();
  const double time_fac = fac_vel / fac_disp;
  Teuchos::ParameterList springdashpotparams;
  springdashpotparams.set("time_fac", time_fac);

  // loop over all spring dashpot conditions and evaluate them
  for (int i=0; i<n_conds_; ++i)
    springs_[i]->EvaluateForceStiff(*stiff_spring_ptr_, *fspring_np_ptr_,
        disnp_ptr_, velnp_ptr_, springdashpotparams);

  if (not stiff_spring_ptr_->Filled())
    stiff_spring_ptr_->Complete();

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::SpringDashpot::EvaluateForceStiff()
{
  CheckInitSetup();

  // get displacement DOFs
  fspring_np_ptr_ =
      Teuchos::rcp(new Epetra_Vector(*GState().DofRowMap(),true));

  // factors from time-integrator for derivative of d(v_{n+1}) / d(d_{n+1})
  // needed for stiffness contribution from dashpot
  const double fac_vel = EvalData().GetTimIntFactorVel();
  const double fac_disp = EvalData().GetTimIntFactorDisp();
  const double time_fac = fac_vel / fac_disp;
  Teuchos::ParameterList springdashpotparams;
  springdashpotparams.set("time_fac", time_fac);

  // loop over all spring dashpot conditions and evaluate them
  for (int i=0; i<n_conds_; ++i)
    springs_[i]->EvaluateForceStiff(*stiff_spring_ptr_, *fspring_np_ptr_,
        disnp_ptr_, velnp_ptr_, springdashpotparams);

  if (not stiff_spring_ptr_->Filled())
    stiff_spring_ptr_->Complete();

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::SpringDashpot::AssembleForce(Epetra_Vector& f,
    const double & timefac_np) const
{
  STR::AssembleVector(1.0,f,timefac_np,*fspring_np_ptr_);
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::SpringDashpot::AssembleJacobian(
    LINALG::SparseOperator& jac,
    const double & timefac_np) const
{
  Teuchos::RCP<LINALG::SparseMatrix> jac_dd_ptr =
      GState().ExtractDisplBlock(jac);
  jac_dd_ptr->Add(*stiff_spring_ptr_,false,timefac_np,1.0);
  // no need to keep it
  stiff_spring_ptr_->Zero();
  // nothing to do
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
void STR::MODELEVALUATOR::SpringDashpot::UpdateStepState(
    const double& timefac_n)
{
  // add the old time factor scaled contributions to the residual
  Teuchos::RCP<Epetra_Vector>& fstructold_ptr =
      GState().GetMutableFstructureOld();
  fstructold_ptr->Update(timefac_n,*fspring_np_ptr_,1.0);
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
void STR::MODELEVALUATOR::SpringDashpot::ResetStepState()
{
  CheckInitSetup();

  dserror("Not yet implemented");

  return;
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

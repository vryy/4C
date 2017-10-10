/*---------------------------------------------------------------------*/
/*!
\file str_model_evaluator_cardiovascular0d.cpp

\brief Evaluation and assembly of all 0D cardiovascular model terms

\maintainer Marc Hirschvogel

\date Jun 29, 2016

\level 3

*/
/*---------------------------------------------------------------------*/

#include "str_model_evaluator_cardiovascular0d.H"
#include "str_model_evaluator_data.H"

#include "str_timint_base.H"
#include "str_utils.H"
#include "str_integrator.H"

#include <Epetra_Vector.h>
#include <Epetra_Time.h>
#include <Teuchos_ParameterList.hpp>

#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_io/io.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::Cardiovascular0D::Cardiovascular0D()
    : disnp_ptr_(Teuchos::null),
      stiff_cardio_ptr_(Teuchos::null),
      fstructcardio_np_ptr_(Teuchos::null)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::Setup()
{
  CheckInit();

  Teuchos::RCP<DRT::Discretization> dis =
      Teuchos::rcp_dynamic_cast<DRT::Discretization>(DiscretPtr(),true);

  // setup the displacement pointer
  disnp_ptr_ = GState().GetDisNp();

  // contributions of 0D model to structural rhs and stiffness
  fstructcardio_np_ptr_ =
      Teuchos::rcp(new Epetra_Vector(*GState().DofRowMapView()));
  stiff_cardio_ptr_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *GState().DofRowMapView(), 81, true, true));

  Teuchos::RCP<LINALG::Solver> dummysolver(new LINALG::Solver(disnp_ptr_->Comm()));

  // ToDo: we do not want to hand in the structural dynamics parameter list
  // to the manager in the future! -> get rid of it as soon as old
  // time-integration dies ...
  // initialize 0D cardiovascular manager
  cardvasc0dman_ = Teuchos::rcp(new UTILS::Cardiovascular0DManager(dis,
      disnp_ptr_,
      DRT::Problem::Instance()->StructuralDynamicParams(),
      DRT::Problem::Instance()->Cardiovascular0DStructuralParams(),
      *dummysolver,
      Teuchos::null));

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::Reset(const Epetra_Vector& x)
{
  CheckInitSetup();

  // update the structural displacement vector
  disnp_ptr_ = GState().GetDisNp();

  fstructcardio_np_ptr_->Scale(0.0);
  stiff_cardio_ptr_->Zero();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Cardiovascular0D::EvaluateForce()
{
  CheckInitSetup();

  double time_np = GState().GetTimeNp();
  Teuchos::ParameterList pcardvasc0d;
  pcardvasc0d.set("time_step_size", (*GState().GetDeltaTime())[0]);

  // only forces are evaluated!
  cardvasc0dman_->EvaluateForceStiff(time_np, disnp_ptr_,
      fstructcardio_np_ptr_, Teuchos::null, pcardvasc0d);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Cardiovascular0D::EvaluateStiff()
{
  CheckInitSetup();

  double time_np = GState().GetTimeNp();
  Teuchos::ParameterList pcardvasc0d;
  pcardvasc0d.set("time_step_size", (*GState().GetDeltaTime())[0]);

  // only stiffnesses are evaluated!
  cardvasc0dman_->EvaluateForceStiff(time_np, disnp_ptr_, Teuchos::null,
      stiff_cardio_ptr_, pcardvasc0d);

  if (not stiff_cardio_ptr_->Filled())
    stiff_cardio_ptr_->Complete();

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Cardiovascular0D::EvaluateForceStiff()
{
  CheckInitSetup();

  double time_np = GState().GetTimeNp();
  Teuchos::ParameterList pcardvasc0d;
  pcardvasc0d.set("time_step_size", (*GState().GetDeltaTime())[0]);

  cardvasc0dman_->EvaluateForceStiff(time_np, disnp_ptr_,
      fstructcardio_np_ptr_, stiff_cardio_ptr_, pcardvasc0d);

  if (not stiff_cardio_ptr_->Filled())
    stiff_cardio_ptr_->Complete();

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Cardiovascular0D::AssembleForce(
    Epetra_Vector& f, const double& timefac_np)
    const
{
  Teuchos::RCP<const Epetra_Vector> block_vec_ptr = Teuchos::null;

  // assemble and scale with str time-integrator dependent value
  LINALG::AssembleMyVector(1.0,f,timefac_np,*fstructcardio_np_ptr_);

  // assemble 0D model rhs - already at the generalized mid-point t_{n+theta} !
  block_vec_ptr = cardvasc0dman_->GetCardiovascular0DRHS();

  if (block_vec_ptr.is_null())
    dserror("The 0D cardiovascular model vector is a NULL pointer, although \n"
        "the structural part indicates, that 0D cardiovascular model contributions \n"
        "are present!");

  const int elements_f = f.Map().NumGlobalElements();
  const int max_gid = GetBlockDofRowMapPtr()->MaxAllGID();
  // only call when f is the full rhs of the coupled problem (not for structural
  // equilibriate initial state call)
  if (elements_f==max_gid+1)
    LINALG::AssembleMyVector(1.0,f,1.0,*block_vec_ptr);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Cardiovascular0D::AssembleJacobian(
    LINALG::SparseOperator& jac, const double& timefac_np) const
{
  Teuchos::RCP<LINALG::SparseMatrix> block_ptr = Teuchos::null;

  // --- Kdd - block - scale with str time-integrator dependent value---
  Teuchos::RCP<LINALG::SparseMatrix> jac_dd_ptr =
      GState().ExtractDisplBlock(jac);
  jac_dd_ptr->Add(*stiff_cardio_ptr_,false,timefac_np,1.0);
  // no need to keep it
  stiff_cardio_ptr_->Zero();

  // --- Kdz - block ---------------------------------------------------
  block_ptr = cardvasc0dman_->GetMatDstructDcv0ddof();
  // scale with str time-integrator dependent value
  block_ptr->Scale(timefac_np);
  GState().AssignModelBlock(jac,*block_ptr,Type(),
      DRT::UTILS::block_displ_lm);
  // reset the block pointer, just to be on the safe side
  block_ptr = Teuchos::null;

  // --- Kzd - block - already scaled correctly by 0D model !-----------
  block_ptr = cardvasc0dman_->GetMatDcardvasc0dDd()->Transpose();
  GState().AssignModelBlock(jac,*block_ptr,Type(),
      DRT::UTILS::block_lm_displ);
  // reset the block pointer, just to be on the safe side
  block_ptr = Teuchos::null;

  // --- Kzz - block - already scaled with 0D theta by 0D model !-------
  block_ptr = cardvasc0dman_->GetCardiovascular0DStiffness();
  GState().AssignModelBlock(jac,*block_ptr,Type(),
      DRT::UTILS::block_lm_lm);
  // reset the block pointer, just to be on the safe side
  block_ptr = Teuchos::null;

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::WriteRestart(
        IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const
{

  iowriter.WriteVector("cv0d_df_np",
                        cardvasc0dman_->Get0D_df_np());
  iowriter.WriteVector("cv0d_f_np",
                        cardvasc0dman_->Get0D_f_np());

  iowriter.WriteVector("cv0d_dof_np",
                        cardvasc0dman_->Get0D_dof_np());
  iowriter.WriteVector("vol_np",
                        cardvasc0dman_->Get0D_vol_np());

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::ReadRestart(
    IO::DiscretizationReader& ioreader)
{

  double time_n = GState().GetTimeN();
  cardvasc0dman_->ReadRestart(ioreader,time_n);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::RecoverState(
    const Epetra_Vector& xold,
    const Epetra_Vector& dir,
    const Epetra_Vector& xnew)
{
  CheckInitSetup();

  Teuchos::RCP<Epetra_Vector> cv0d_incr = GState().ExtractModelEntries(
      INPAR::STR::model_cardiovascular0d,dir);

  cardvasc0dman_->UpdateCv0DDof(cv0d_incr);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::UpdateStepState(
    const double& timefac_n)
{
  // only update 0D model at the end of the time step!
  if (EvalData().GetTotalTime()==GState().GetTimeNp())
    cardvasc0dman_->UpdateTimeStep();

  // only print state variables after a finished time step, not when we're
  // in the equilibriate initial state routine
  if (EvalData().GetTotalTime()==GState().GetTimeNp())
    cardvasc0dman_->PrintPresFlux(false);

  // add the 0D cardiovascular force contributions to the old structural
  // residual state vector
  if (not fstructcardio_np_ptr_.is_null())
  {
    Teuchos::RCP<Epetra_Vector>& fstructold_ptr =
        GState().GetMutableFstructureOld();
    fstructold_ptr->Update(timefac_n,*fstructcardio_np_ptr_,1.0);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::UpdateStepElement()
{
  // nothing to do
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::DetermineStressStrain()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::DetermineEnergy()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::OutputStepState(
    IO::DiscretizationWriter& iowriter) const
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::ResetStepState()
{
  CheckInitSetup();

  dserror("Not yet implemented");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::Cardiovascular0D::
    GetBlockDofRowMapPtr() const
{
  CheckInitSetup();

  return cardvasc0dman_->GetCardiovascular0DMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Cardiovascular0D::
    GetCurrentSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Cardiovascular0D::
    GetLastTimeStepSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::PostOutput()
{
  CheckInitSetup();
 // empty

  return;
} // PostOutput()

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

  Teuchos::RCP<DRT::Discretization> dis = GState().GetMutableDiscret();

  // setup the displacement pointer
  disnp_ptr_ = GState().GetMutableDisNp();

  // contributions of 0D model to structural rhs and stiffness
  fstructcardio_np_ptr_ =
      Teuchos::rcp(new Epetra_Vector(*GState().DofRowMapView()));
  stiff_cardio_ptr_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *GState().DofRowMapView(), 81, true, true));

  Teuchos::RCP<LINALG::Solver> dummysolver;

  // initialize 0D cardiovascular manager
  cardvasc0dman_ = Teuchos::rcp(new UTILS::Cardiovascular0DManager(dis,
      disnp_ptr_,
      DRT::Problem::Instance()->StructuralDynamicParams(),
      DRT::Problem::Instance()->Cardiovascular0DStructuralParams(),
      *dummysolver));

//  cardvasc0dman_->PrintNewtonHeader();

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::Reset(const Epetra_Vector& x)
{
  CheckInitSetup();

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

  cardvasc0dman_->EvaluateForceStiff(time_np, disnp_ptr_, fstructcardio_np_ptr_,
      stiff_cardio_ptr_, pcardvasc0d);

  if (not stiff_cardio_ptr_->Filled())
    stiff_cardio_ptr_->Complete();

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Cardiovascular0D::AssembleForce(
    Epetra_Vector& f,const double & timefac_np)
    const
{

  STR::AssembleVector(1.0,f,1.0,*fstructcardio_np_ptr_);

  Teuchos::RCP<const Epetra_Vector> block_vec_ptr = Teuchos::null;
  if (timefac_np!=1.0)
    dserror("You have to consider the corresponding time integration"
        " factors.");

  // assemble 0D model rhs
  block_vec_ptr = cardvasc0dman_->GetCardiovascular0DRHS();

  if (block_vec_ptr.is_null())
    dserror("The 0D cardiovascular model vector is a NULL pointer, although \n"
        "the structural part indicates, that 0D cardiovascular model contributions \n"
        "are present!");

  STR::AssembleVector(1.0,f,1.0,*block_vec_ptr);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Cardiovascular0D::AssembleJacobian(
    LINALG::SparseOperator& jac, const double & timefac_np) const
{
  Teuchos::RCP<const LINALG::SparseMatrix> block_ptr = Teuchos::null;
  if (timefac_np!=1.0)
    dserror("You have to consider the corresponding time integration"
        " factors.");

  // --- Kdd - block ---------------------------------------------------
  Teuchos::RCP<LINALG::SparseMatrix> jac_dd_ptr =
      GState().ExtractDisplBlock(jac);
  jac_dd_ptr->Add(*stiff_cardio_ptr_,false,1.0,1.0);
  // no need to keep it
  stiff_cardio_ptr_->Zero();
  // --- Kdz - block ---------------------------------------------------
  block_ptr = cardvasc0dman_->GetMatDstructDcv0ddof();

  GState().AssignModelBlock(jac,*block_ptr,Type(),
      STR::block_displ_lm);
  // reset the block pointer, just to be on the safe side
  block_ptr = Teuchos::null;
  // --- Kzd - block ---------------------------------------------------
  block_ptr = cardvasc0dman_->GetMatDcardvasc0dDd()->Transpose();

  GState().AssignModelBlock(jac,*block_ptr,Type(),
      STR::block_lm_displ);
  // reset the block pointer, just to be on the safe side
  block_ptr = Teuchos::null;
  // --- Kzz - block ---------------------------------------------------
  block_ptr = cardvasc0dman_->GetCardiovascular0DStiffness();

  GState().AssignModelBlock(jac,*block_ptr,Type(),
      STR::block_lm_lm);
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

  iowriter.WriteVector("cv0ddof",
                        cardvasc0dman_->Get0DDofVector());
  iowriter.WriteVector("refvolval",
                        cardvasc0dman_->GetRefVolValue());
  iowriter.WriteVector("reffluxval",
                        cardvasc0dman_->GetRefFluxValue());
  iowriter.WriteVector("refdfluxval",
                        cardvasc0dman_->GetRefDFluxValue());
  iowriter.WriteVector("refddfluxval",
                        cardvasc0dman_->GetRefDDFluxValue());

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

  Teuchos::RCP<Epetra_Vector> cv0d_incr = GState().ExtractModelEntries(
      INPAR::STR::model_cardiovascular0d,dir);
  Teuchos::RCP<Epetra_Vector> dis_incr = GState().ExtractModelEntries(
      INPAR::STR::model_structure,dir);

  cardvasc0dman_->UpdateCv0DDof(cv0d_incr);

  // store for manager-internal monitoring...
  cardvasc0dman_->StoreCv0dDofIncrement(cv0d_incr);
  cardvasc0dman_->StoreStructuralDisplIncrement(dis_incr);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Cardiovascular0D::UpdateStepState(
    const double& timefac_n)
{
  cardvasc0dman_->UpdateTimeStep();
  cardvasc0dman_->PrintPresFlux(false);

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

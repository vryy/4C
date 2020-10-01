/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all spring dashpot terms


\level 3

*/
/*-----------------------------------------------------------*/

#include "str_model_evaluator_springdashpot.H"
#include "str_timint_base.H"
#include "str_utils.H"

#include <Epetra_Vector.h>
#include <Epetra_Time.h>
#include <Teuchos_ParameterList.hpp>

#include "../drt_inpar/inpar_structure.H"

#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/prestress_service.H"

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
  if (not IsInit()) dserror("Init() has not been called, yet!");

  // get all spring dashpot conditions
  std::vector<Teuchos::RCP<DRT::Condition>> springdashpots;
  Discret().GetCondition("RobinSpringDashpot", springdashpots);

  // number of spring dashpot conditions
  n_conds_ = (int)springdashpots.size();

  // new instance of spring dashpot BC for each condition
  Teuchos::RCP<DRT::Discretization> discret_ptr =
      Teuchos::rcp_dynamic_cast<DRT::Discretization>(DiscretPtr(), true);
  for (int i = 0; i < n_conds_; ++i)
    springs_.push_back(Teuchos::rcp(new UTILS::SpringDashpot(discret_ptr, springdashpots[i])));

  // setup the displacement pointer
  disnp_ptr_ = GState().GetMutableDisNp();
  velnp_ptr_ = GState().GetMutableVelNp();

  fspring_np_ptr_ = Teuchos::rcp(new Epetra_Vector(*GState().DofRowMapView()));
  stiff_spring_ptr_ =
      Teuchos::rcp(new LINALG::SparseMatrix(*GState().DofRowMapView(), 81, true, true));

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::Reset(const Epetra_Vector& x)
{
  CheckInitSetup();

  // loop over all spring dashpot conditions and reset them
  for (int i = 0; i < n_conds_; ++i) springs_[i]->ResetNewton();

  // update the structural displacement vector
  disnp_ptr_ = GState().GetDisNp();

  // update the structural displacement vector
  velnp_ptr_ = GState().GetVelNp();

  fspring_np_ptr_->Scale(0.0);
  stiff_spring_ptr_->Zero();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::SpringDashpot::EvaluateForce()
{
  CheckInitSetup();

  Teuchos::ParameterList springdashpotparams;
  // loop over all spring dashpot conditions and evaluate them
  fspring_np_ptr_ = Teuchos::rcp(new Epetra_Vector(*GState().DofRowMapView()));
  for (int i = 0; i < n_conds_; ++i)
  {
    const UTILS::SpringDashpot::SpringType stype = springs_[i]->GetSpringType();

    if (stype == UTILS::SpringDashpot::xyz or stype == UTILS::SpringDashpot::refsurfnormal)
    {
      springdashpotparams.set("total time", GState().GetTimeNp());
      springs_[i]->EvaluateRobin(
          Teuchos::null, fspring_np_ptr_, disnp_ptr_, velnp_ptr_, springdashpotparams);
    }
    if (stype == UTILS::SpringDashpot::cursurfnormal)
    {
      springdashpotparams.set("dt", (*GState().GetDeltaTime())[0]);
      springs_[i]->EvaluateForce(*fspring_np_ptr_, disnp_ptr_, velnp_ptr_, springdashpotparams);
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::SpringDashpot::EvaluateStiff()
{
  CheckInitSetup();

  fspring_np_ptr_ = Teuchos::rcp(new Epetra_Vector(*GState().DofRowMapView(), true));

  // factors from time-integrator for derivative of d(v_{n+1}) / d(d_{n+1})
  // needed for stiffness contribution from dashpot
  const double fac_vel = EvalData().GetTimIntFactorVel();
  const double fac_disp = EvalData().GetTimIntFactorDisp();
  const double time_fac = fac_vel / fac_disp;
  Teuchos::ParameterList springdashpotparams;
  if (fac_vel > 0.0) springdashpotparams.set("time_fac", time_fac);

  // loop over all spring dashpot conditions and evaluate them
  for (int i = 0; i < n_conds_; ++i)
  {
    const UTILS::SpringDashpot::SpringType stype = springs_[i]->GetSpringType();

    if (stype == UTILS::SpringDashpot::xyz or stype == UTILS::SpringDashpot::refsurfnormal)
    {
      springdashpotparams.set("total time", GState().GetTimeNp());
      springs_[i]->EvaluateRobin(
          stiff_spring_ptr_, Teuchos::null, disnp_ptr_, velnp_ptr_, springdashpotparams);
    }
    if (stype == UTILS::SpringDashpot::cursurfnormal)
    {
      springdashpotparams.set("dt", (*GState().GetDeltaTime())[0]);
      springs_[i]->EvaluateForceStiff(
          *stiff_spring_ptr_, *fspring_np_ptr_, disnp_ptr_, velnp_ptr_, springdashpotparams);
    }
  }

  if (not stiff_spring_ptr_->Filled()) stiff_spring_ptr_->Complete();

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::SpringDashpot::EvaluateForceStiff()
{
  CheckInitSetup();

  // get displacement DOFs
  fspring_np_ptr_ = Teuchos::rcp(new Epetra_Vector(*GState().DofRowMap(), true));

  // factors from time-integrator for derivative of d(v_{n+1}) / d(d_{n+1})
  // needed for stiffness contribution from dashpot
  const double fac_vel = EvalData().GetTimIntFactorVel();
  const double fac_disp = EvalData().GetTimIntFactorDisp();
  const double time_fac = fac_vel / fac_disp;
  Teuchos::ParameterList springdashpotparams;

  if (fac_vel > 0.0) springdashpotparams.set("time_fac", time_fac);

  // loop over all spring dashpot conditions and evaluate them
  for (int i = 0; i < n_conds_; ++i)
  {
    const UTILS::SpringDashpot::SpringType stype = springs_[i]->GetSpringType();

    if (stype == UTILS::SpringDashpot::xyz or stype == UTILS::SpringDashpot::refsurfnormal)
    {
      springdashpotparams.set("total time", GState().GetTimeNp());
      springs_[i]->EvaluateRobin(
          stiff_spring_ptr_, fspring_np_ptr_, disnp_ptr_, velnp_ptr_, springdashpotparams);
    }
    if (stype == UTILS::SpringDashpot::cursurfnormal)
    {
      springdashpotparams.set("dt", (*GState().GetDeltaTime())[0]);
      springs_[i]->EvaluateForceStiff(
          *stiff_spring_ptr_, *fspring_np_ptr_, disnp_ptr_, velnp_ptr_, springdashpotparams);
    }
  }

  if (not stiff_spring_ptr_->Filled()) stiff_spring_ptr_->Complete();

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::SpringDashpot::AssembleForce(
    Epetra_Vector& f, const double& timefac_np) const
{
  LINALG::AssembleMyVector(1.0, f, timefac_np, *fspring_np_ptr_);
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::SpringDashpot::AssembleJacobian(
    LINALG::SparseOperator& jac, const double& timefac_np) const
{
  Teuchos::RCP<LINALG::SparseMatrix> jac_dd_ptr = GState().ExtractDisplBlock(jac);
  jac_dd_ptr->Add(*stiff_spring_ptr_, false, timefac_np, 1.0);
  // no need to keep it
  stiff_spring_ptr_->Zero();
  // nothing to do
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::WriteRestart(
    IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  // row maps for export
  Teuchos::RCP<Epetra_Vector> springoffsetprestr =
      Teuchos::rcp(new Epetra_Vector(*Discret().DofRowMap()));
  Teuchos::RCP<Epetra_MultiVector> springoffsetprestr_old =
      Teuchos::rcp(new Epetra_MultiVector(*(Discret().NodeRowMap()), 3, true));

  // collect outputs from all spring dashpot conditions
  for (int i = 0; i < n_conds_; ++i)
  {
    // get spring type from current condition
    const UTILS::SpringDashpot::SpringType stype = springs_[i]->GetSpringType();

    if (stype == UTILS::SpringDashpot::xyz or stype == UTILS::SpringDashpot::refsurfnormal)
      springs_[i]->OutputPrestrOffset(springoffsetprestr);
    if (stype == UTILS::SpringDashpot::cursurfnormal)
      springs_[i]->OutputPrestrOffsetOld(springoffsetprestr_old);
  }

  // write vector to output for restart
  iowriter.WriteVector("springoffsetprestr", springoffsetprestr);
  // write vector to output for restart
  iowriter.WriteVector("springoffsetprestr_old", springoffsetprestr_old);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::ReadRestart(IO::DiscretizationReader& ioreader)
{
  Teuchos::RCP<Epetra_Vector> tempvec = Teuchos::rcp(new Epetra_Vector(*Discret().DofRowMap()));
  Teuchos::RCP<Epetra_MultiVector> tempvecold =
      Teuchos::rcp(new Epetra_MultiVector(*(Discret().NodeRowMap()), 3, true));

  ioreader.ReadVector(tempvec, "springoffsetprestr");
  ioreader.ReadMultiVector(tempvecold, "springoffsetprestr_old");

  // loop over all spring dashpot conditions and set restart
  for (int i = 0; i < n_conds_; ++i)
  {
    // get spring type from current condition
    const UTILS::SpringDashpot::SpringType stype = springs_[i]->GetSpringType();

    if (stype == UTILS::SpringDashpot::xyz or stype == UTILS::SpringDashpot::refsurfnormal)
      springs_[i]->SetRestart(tempvec);
    if (stype == UTILS::SpringDashpot::cursurfnormal) springs_[i]->SetRestartOld(tempvecold);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::RunPostComputeX(
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::UpdateStepState(const double& timefac_n)
{
  // add the old time factor scaled contributions to the residual
  Teuchos::RCP<Epetra_Vector>& fstructold_ptr = GState().GetMutableFstructureOld();
  fstructold_ptr->Update(timefac_n, *fspring_np_ptr_, 1.0);

  // check for prestressing and reset if necessary
  const INPAR::STR::PreStress prestress_type = TimInt().GetDataSDyn().GetPreStressType();
  const double prestress_time = TimInt().GetDataSDyn().GetPreStressTime();

  if (::UTILS::PRESTRESS::IsActive(GState().GetTimeNp(), prestress_type, prestress_time))
  {
    switch (prestress_type)
    {
      case INPAR::STR::PreStress::mulf:
      case INPAR::STR::PreStress::material_iterative:
        for (auto spring : springs_) spring->ResetPrestress(GState().GetDisNp());
      default:
        break;
    }
  }
  for (int i = 0; i < n_conds_; ++i) springs_[i]->Update();
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
void STR::MODELEVALUATOR::SpringDashpot::DetermineOptionalQuantity()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::OutputStepState(IO::DiscretizationWriter& iowriter) const
{
  // row maps for export
  Teuchos::RCP<Epetra_Vector> gap =
      Teuchos::rcp(new Epetra_Vector(*(Discret().NodeRowMap()), true));
  Teuchos::RCP<Epetra_MultiVector> normals =
      Teuchos::rcp(new Epetra_MultiVector(*(Discret().NodeRowMap()), 3, true));
  Teuchos::RCP<Epetra_MultiVector> springstress =
      Teuchos::rcp(new Epetra_MultiVector(*(Discret().NodeRowMap()), 3, true));

  // collect outputs from all spring dashpot conditions
  bool found_cursurfnormal = false;
  for (int i = 0; i < n_conds_; ++i)
  {
    springs_[i]->OutputGapNormal(gap, normals, springstress);

    // get spring type from current condition
    const UTILS::SpringDashpot::SpringType stype = springs_[i]->GetSpringType();
    if (stype == UTILS::SpringDashpot::cursurfnormal) found_cursurfnormal = true;
  }

  // write vectors to output
  if (found_cursurfnormal)
  {
    iowriter.WriteVector("gap", gap);
    iowriter.WriteVector("curnormals", normals);
  }

  // write spring stress if defined in io-flag
  if (DRT::INPUT::IntegralValue<bool>(DRT::Problem::Instance()->IOParams(), "OUTPUT_SPRING") ==
      true)
    iowriter.WriteVector("springstress", springstress);

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
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::SpringDashpot::GetBlockDofRowMapPtr() const
{
  CheckInitSetup();
  return GState().DofRowMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::SpringDashpot::GetCurrentSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::SpringDashpot::GetLastTimeStepSolutionPtr()
    const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::PostOutput()
{
  CheckInitSetup();
  // empty

  return;
}  // PostOutput()

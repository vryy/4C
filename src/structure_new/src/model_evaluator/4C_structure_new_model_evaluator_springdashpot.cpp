/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all spring dashpot terms


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_model_evaluator_springdashpot.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_io.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_utils.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::SpringDashpot::SpringDashpot()
    : disnp_ptr_(Teuchos::null),
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
  FOUR_C_ASSERT(is_init(), "Init() has not been called, yet!");

  // get all spring dashpot conditions
  std::vector<Teuchos::RCP<Core::Conditions::Condition>> springdashpots;
  discret().GetCondition("RobinSpringDashpot", springdashpots);

  // new instance of spring dashpot BC for each condition
  for (auto& springdashpot : springdashpots)
    springs_.emplace_back(
        Teuchos::rcp(new CONSTRAINTS::SpringDashpot(discret_ptr(), springdashpot)));

  // setup the displacement pointer
  disnp_ptr_ = g_state().GetDisNp();
  velnp_ptr_ = g_state().GetVelNp();

  fspring_np_ptr_ = Teuchos::rcp(new Epetra_Vector(*g_state().DofRowMapView()));
  stiff_spring_ptr_ =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*g_state().DofRowMapView(), 81, true, true));

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::Reset(const Epetra_Vector& x)
{
  check_init_setup();

  // loop over all spring dashpot conditions and reset them
  for (const auto& spring : springs_) spring->ResetNewton();

  // update the structural displacement vector
  disnp_ptr_ = g_state().GetDisNp();

  // update the structural displacement vector
  velnp_ptr_ = g_state().GetVelNp();

  fspring_np_ptr_->PutScalar(0.0);
  stiff_spring_ptr_->Zero();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::SpringDashpot::evaluate_force()
{
  check_init_setup();

  Teuchos::ParameterList springdashpotparams;
  // loop over all spring dashpot conditions and evaluate them
  fspring_np_ptr_ = Teuchos::rcp(new Epetra_Vector(*g_state().DofRowMapView()));
  for (const auto& spring : springs_)
  {
    const CONSTRAINTS::SpringDashpot::SpringType stype = spring->GetSpringType();

    if (stype == CONSTRAINTS::SpringDashpot::xyz or
        stype == CONSTRAINTS::SpringDashpot::refsurfnormal)
    {
      springdashpotparams.set("total time", g_state().GetTimeNp());
      spring->EvaluateRobin(
          Teuchos::null, fspring_np_ptr_, disnp_ptr_, velnp_ptr_, springdashpotparams);
    }
    if (stype == CONSTRAINTS::SpringDashpot::cursurfnormal)
    {
      springdashpotparams.set("dt", (*g_state().GetDeltaTime())[0]);
      spring->evaluate_force(*fspring_np_ptr_, disnp_ptr_, velnp_ptr_, springdashpotparams);
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::SpringDashpot::evaluate_stiff()
{
  check_init_setup();

  fspring_np_ptr_ = Teuchos::rcp(new Epetra_Vector(*g_state().DofRowMapView(), true));

  // factors from time-integrator for derivative of d(v_{n+1}) / d(d_{n+1})
  // needed for stiffness contribution from dashpot
  const double fac_vel = eval_data().GetTimIntFactorVel();
  const double fac_disp = eval_data().GetTimIntFactorDisp();
  const double time_fac = fac_vel / fac_disp;
  Teuchos::ParameterList springdashpotparams;
  if (fac_vel > 0.0) springdashpotparams.set("time_fac", time_fac);

  // loop over all spring dashpot conditions and evaluate them
  for (const auto& spring : springs_)
  {
    const CONSTRAINTS::SpringDashpot::SpringType stype = spring->GetSpringType();

    if (stype == CONSTRAINTS::SpringDashpot::xyz or
        stype == CONSTRAINTS::SpringDashpot::refsurfnormal)
    {
      springdashpotparams.set("total time", g_state().GetTimeNp());
      spring->EvaluateRobin(
          stiff_spring_ptr_, Teuchos::null, disnp_ptr_, velnp_ptr_, springdashpotparams);
    }
    if (stype == CONSTRAINTS::SpringDashpot::cursurfnormal)
    {
      springdashpotparams.set("dt", (*g_state().GetDeltaTime())[0]);
      spring->evaluate_force_stiff(
          *stiff_spring_ptr_, *fspring_np_ptr_, disnp_ptr_, velnp_ptr_, springdashpotparams);
    }
  }

  if (not stiff_spring_ptr_->Filled()) stiff_spring_ptr_->Complete();

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::SpringDashpot::evaluate_force_stiff()
{
  check_init_setup();

  // get displacement DOFs
  fspring_np_ptr_ = Teuchos::rcp(new Epetra_Vector(*g_state().dof_row_map(), true));

  // factors from time-integrator for derivative of d(v_{n+1}) / d(d_{n+1})
  // needed for stiffness contribution from dashpot
  const double fac_vel = eval_data().GetTimIntFactorVel();
  const double fac_disp = eval_data().GetTimIntFactorDisp();
  const double time_fac = fac_vel / fac_disp;
  Teuchos::ParameterList springdashpotparams;

  if (fac_vel > 0.0) springdashpotparams.set("time_fac", time_fac);

  // loop over all spring dashpot conditions and evaluate them
  for (const auto& spring : springs_)
  {
    const CONSTRAINTS::SpringDashpot::SpringType stype = spring->GetSpringType();

    if (stype == CONSTRAINTS::SpringDashpot::xyz or
        stype == CONSTRAINTS::SpringDashpot::refsurfnormal)
    {
      springdashpotparams.set("total time", g_state().GetTimeNp());
      spring->EvaluateRobin(
          stiff_spring_ptr_, fspring_np_ptr_, disnp_ptr_, velnp_ptr_, springdashpotparams);
    }
    if (stype == CONSTRAINTS::SpringDashpot::cursurfnormal)
    {
      springdashpotparams.set("dt", (*g_state().GetDeltaTime())[0]);
      spring->evaluate_force_stiff(
          *stiff_spring_ptr_, *fspring_np_ptr_, disnp_ptr_, velnp_ptr_, springdashpotparams);
    }
  }

  if (not stiff_spring_ptr_->Filled()) stiff_spring_ptr_->Complete();

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::SpringDashpot::assemble_force(
    Epetra_Vector& f, const double& timefac_np) const
{
  Core::LinAlg::AssembleMyVector(1.0, f, timefac_np, *fspring_np_ptr_);
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::SpringDashpot::assemble_jacobian(
    Core::LinAlg::SparseOperator& jac, const double& timefac_np) const
{
  Teuchos::RCP<Core::LinAlg::SparseMatrix> jac_dd_ptr = GState().ExtractDisplBlock(jac);
  jac_dd_ptr->Add(*stiff_spring_ptr_, false, timefac_np, 1.0);
  // no need to keep it
  stiff_spring_ptr_->Zero();
  // nothing to do
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::write_restart(
    Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  // row maps for export
  Teuchos::RCP<Epetra_Vector> springoffsetprestr =
      Teuchos::rcp(new Epetra_Vector(*Discret().dof_row_map()));
  Teuchos::RCP<Epetra_MultiVector> springoffsetprestr_old =
      Teuchos::rcp(new Epetra_MultiVector(*(Discret().NodeRowMap()), 3, true));

  // collect outputs from all spring dashpot conditions
  for (const auto& spring : springs_)
  {
    // get spring type from current condition
    const CONSTRAINTS::SpringDashpot::SpringType stype = spring->GetSpringType();

    if (stype == CONSTRAINTS::SpringDashpot::xyz or
        stype == CONSTRAINTS::SpringDashpot::refsurfnormal)
      spring->OutputPrestrOffset(springoffsetprestr);
    if (stype == CONSTRAINTS::SpringDashpot::cursurfnormal)
      spring->output_prestr_offset_old(springoffsetprestr_old);
  }

  // write vector to output for restart
  iowriter.WriteVector("springoffsetprestr", springoffsetprestr);
  // write vector to output for restart
  iowriter.WriteVector("springoffsetprestr_old", springoffsetprestr_old);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::read_restart(Core::IO::DiscretizationReader& ioreader)
{
  Teuchos::RCP<Epetra_Vector> tempvec = Teuchos::rcp(new Epetra_Vector(*discret().dof_row_map()));
  Teuchos::RCP<Epetra_MultiVector> tempvecold =
      Teuchos::rcp(new Epetra_MultiVector(*(discret().NodeRowMap()), 3, true));

  ioreader.ReadVector(tempvec, "springoffsetprestr");
  ioreader.ReadMultiVector(tempvecold, "springoffsetprestr_old");

  // loop over all spring dashpot conditions and set restart
  for (const auto& spring : springs_)
  {
    // get spring type from current condition
    const CONSTRAINTS::SpringDashpot::SpringType stype = spring->GetSpringType();

    if (stype == CONSTRAINTS::SpringDashpot::xyz or
        stype == CONSTRAINTS::SpringDashpot::refsurfnormal)
      spring->SetRestart(tempvec);
    if (stype == CONSTRAINTS::SpringDashpot::cursurfnormal) spring->SetRestartOld(tempvecold);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::UpdateStepState(const double& timefac_n)
{
  // add the old time factor scaled contributions to the residual
  Teuchos::RCP<Epetra_Vector>& fstructold_ptr = g_state().GetFstructureOld();
  fstructold_ptr->Update(timefac_n, *fspring_np_ptr_, 1.0);

  // check for prestressing and reset if necessary
  const Inpar::STR::PreStress prestress_type = TimInt().GetDataSDyn().GetPreStressType();
  const double prestress_time = TimInt().GetDataSDyn().GetPreStressTime();

  if (prestress_type != Inpar::STR::PreStress::none &&
      g_state().GetTimeNp() <= prestress_time + 1.0e-15)
  {
    switch (prestress_type)
    {
      case Inpar::STR::PreStress::mulf:
      case Inpar::STR::PreStress::material_iterative:
        for (const auto& spring : springs_) spring->ResetPrestress(g_state().GetDisNp());
      default:
        break;
    }
  }
  for (const auto& spring : springs_) spring->Update();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::OutputStepState(
    Core::IO::DiscretizationWriter& iowriter) const
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
  for (const auto& spring : springs_)
  {
    spring->OutputGapNormal(gap, normals, springstress);

    // get spring type from current condition
    const CONSTRAINTS::SpringDashpot::SpringType stype = spring->GetSpringType();
    if (stype == CONSTRAINTS::SpringDashpot::cursurfnormal) found_cursurfnormal = true;
  }

  // write vectors to output
  if (found_cursurfnormal)
  {
    iowriter.WriteVector("gap", gap);
    iowriter.WriteVector("curnormals", normals);
  }

  // write spring stress if defined in io-flag
  if (Core::UTILS::IntegralValue<bool>(Global::Problem::Instance()->IOParams(), "OUTPUT_SPRING"))
    iowriter.WriteVector("springstress", springstress);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::ResetStepState()
{
  check_init_setup();

  for (auto& spring : springs_)
  {
    spring->ResetStepState();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::SpringDashpot::get_block_dof_row_map_ptr() const
{
  check_init_setup();
  return GState().dof_row_map();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::SpringDashpot::get_current_solution_ptr()
    const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
STR::MODELEVALUATOR::SpringDashpot::get_last_time_step_solution_ptr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::SpringDashpot::PostOutput() { check_init_setup(); }

FOUR_C_NAMESPACE_CLOSE

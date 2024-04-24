/*----------------------------------------------------------------------------*/
/** \file
  \brief Lagrange multiplier function: solve a least squares problem to compute
  the Lagrange multiplier value dependent on the current displacement state

  \level 3
*/
/*----------------------------------------------------------------------------*/

#include "baci_contact_aug_lagrange_multiplier_function.hpp"

#include "baci_contact_aug_interface.hpp"
#include "baci_contact_aug_strategy.hpp"
#include "baci_contact_paramsinterface.hpp"
#include "baci_global_data.hpp"
#include "baci_io_control.hpp"
#include "baci_linalg_multiply.hpp"
#include "baci_linalg_sparsematrix.hpp"
#include "baci_linalg_utils_sparse_algebra_create.hpp"
#include "baci_linalg_utils_sparse_algebra_manipulation.hpp"
#include "baci_linear_solver_method_linalg.hpp"
#include "baci_mortar_matrix_transform.hpp"
#include "baci_structure_new_model_evaluator_contact.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::LagrangeMultiplierFunction::LagrangeMultiplierFunction()
    : isinit_(false),
      issetup_(false),
      strategy_(nullptr),
      interfaces_(0),
      data_(Teuchos::null),
      lin_solver_type_(CORE::LINEAR_SOLVER::SolverType::undefined),
      lin_solver_(Teuchos::null),
      bmat_(Teuchos::null)
{
  /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::LagrangeMultiplierFunction::Init(
    const Strategy* const strategy, CONTACT::AUG::DataContainer& data)
{
  issetup_ = false;

  strategy_ = strategy;
  interfaces_.reserve(strategy->ContactInterfaces().size());
  std::copy(strategy->ContactInterfaces().begin(), strategy->ContactInterfaces().end(),
      std::back_inserter(interfaces_));

  data_ = Teuchos::rcpFromRef(data);

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::LagrangeMultiplierFunction::Setup()
{
  CheckInit();

  const Teuchos::ParameterList& p_lm_func =
      strategy_->Params().sublist("AUGMENTED").sublist("LAGRANGE_MULTIPLIER_FUNCTION");

  lin_solver_ =
      CreateLinearSolver(p_lm_func.get<int>("LINEAR_SOLVER"), strategy_->Comm(), lin_solver_type_);

  Redistribute();

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::LagrangeMultiplierFunction::Redistribute()
{
  const Epetra_Map& slMaDofRowMap = *data_->GSlMaDofRowMapPtr();
  bmat_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(slMaDofRowMap, 100, false, false));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::Solver> CONTACT::AUG::LagrangeMultiplierFunction::CreateLinearSolver(
    const int lin_sol_id, const Epetra_Comm& comm,
    enum CORE::LINEAR_SOLVER::SolverType& solver_type) const
{
  if (lin_sol_id == -1) FOUR_C_THROW("You must specify a meaningful LINEAR_SOLVER!");

  // get solver parameter list of linear solver
  const Teuchos::ParameterList& solverparams =
      GLOBAL::Problem::Instance()->SolverParams(lin_sol_id);
  solver_type = Teuchos::getIntegralValue<CORE::LINEAR_SOLVER::SolverType>(solverparams, "SOLVER");

  Teuchos::RCP<CORE::LINALG::Solver> solver =
      Teuchos::rcp(new CORE::LINALG::Solver(solverparams, comm));

  if (solver_type != CORE::LINEAR_SOLVER::SolverType::umfpack and
      solver_type != CORE::LINEAR_SOLVER::SolverType::superlu)
    FOUR_C_THROW("Currently only direct linear solvers are supported!");

  return solver;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::LagrangeMultiplierFunction::LinSolve(
    CORE::LINALG::SparseOperator& mat, Epetra_MultiVector& rhs, Epetra_MultiVector& sol)
{
  if (rhs.NumVectors() > 1 or sol.NumVectors() > 1)
    FOUR_C_THROW("MultiVector support is not yet implemented!");

  CORE::LINALG::SolverParams solver_params;
  solver_params.refactor = true;
  solver_params.reset = true;
  int err = lin_solver_->Solve(
      mat.EpetraOperator(), Teuchos::rcpFromRef(sol), Teuchos::rcpFromRef(rhs), solver_params);

  if (err) FOUR_C_THROW("LinSolve failed with err = %d", err);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> CONTACT::AUG::LagrangeMultiplierFunction::Compute(
    const CONTACT::ParamsInterface& cparams)
{
  TEUCHOS_FUNC_TIME_MONITOR(CONTACT_FUNC_NAME);

  CheckInitSetup();

  Teuchos::RCP<Epetra_Vector> lmn_vec = Teuchos::rcp(new Epetra_Vector(data_->GActiveNDofRowMap()));

  Teuchos::RCP<Epetra_Vector> str_gradient = GetStructureGradient(cparams);

  CreateBMatrix();

  Epetra_Vector str_gradient_exp(*data_->GSlMaDofRowMapPtr(), true);
  CORE::LINALG::Export(*str_gradient, str_gradient_exp);

  Epetra_Vector rhs(data_->GActiveNDofRowMap(), true);
  bmat_->Multiply(true, str_gradient_exp, rhs);

  Teuchos::RCP<CORE::LINALG::SparseMatrix> bbmat =
      CORE::LINALG::MLMultiply(*bmat_, true, *bmat_, false, false, false, true);

  LinSolve(*bbmat, rhs, *lmn_vec);

  return lmn_vec;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> CONTACT::AUG::LagrangeMultiplierFunction::GetStructureGradient(
    const CONTACT::ParamsInterface& cparams) const
{
  const STR::MODELEVALUATOR::Generic& model = cparams.GetModelEvaluator();
  const STR::MODELEVALUATOR::Contact& cmodel =
      dynamic_cast<const STR::MODELEVALUATOR::Contact&>(model);

  const std::vector<INPAR::STR::ModelType> without_contact_model(1, model.Type());

  Teuchos::RCP<Epetra_Vector> str_gradient =
      cmodel.AssembleForceOfModels(&without_contact_model, true);

  return str_gradient;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::LagrangeMultiplierFunction::CreateBMatrix()
{
  bmat_->Reset();

  bmat_->Add(data_->DMatrix(), true, 1.0, 0.0);
  bmat_->Add(data_->MMatrix(), true, 1.0, 1.0);
  //  bmat_->Add( data_->DLmNWGapLinMatrix(), true, 1.0, 0.0 );

  bmat_->Complete(data_->GActiveNDofRowMap(), *data_->GSlMaDofRowMapPtr());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> CONTACT::AUG::LagrangeMultiplierFunction::FirstOrderDirection(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& dincr)
{
  TEUCHOS_FUNC_TIME_MONITOR(CONTACT_FUNC_NAME);

  Epetra_Vector rhs(data_->GActiveNDofRowMap(), true);

  const STR::MODELEVALUATOR::Generic& model = cparams.GetModelEvaluator();
  const STR::MODELEVALUATOR::Contact& cmodel =
      dynamic_cast<const STR::MODELEVALUATOR::Contact&>(model);

  // access the full stiffness matrix
  CORE::LINALG::SparseMatrix full_stiff(
      *cmodel.GetJacobianBlock(STR::MatBlockType::displ_displ), CORE::LINALG::Copy);

  Teuchos::RCP<CORE::LINALG::SparseMatrix> kdd_ptr =
      strategy_->GetMatrixBlockPtr(CONTACT::MatBlockType::displ_displ);

  // undo matrix contributions
  full_stiff.Add(*kdd_ptr, false, -1.0, 1.0);

  // --- first summand
  Teuchos::RCP<Epetra_Vector> tmp_vec = CORE::LINALG::CreateVector(full_stiff.RangeMap(), true);

  int err = full_stiff.Multiply(false, dincr, *tmp_vec);
  if (err) FOUR_C_THROW("Multiply failed with err = %d", err);

  Teuchos::RCP<Epetra_Vector> tmp_vec_exp =
      CORE::LINALG::CreateVector(*data_->GSlMaDofRowMapPtr(), true);

  // build necessary exporter
  Epetra_Export exporter(tmp_vec_exp->Map(), tmp_vec->Map());
  err = tmp_vec_exp->Import(*tmp_vec, exporter, Insert);
  if (err) FOUR_C_THROW("Import failed with err = %d", err);

  CreateBMatrix();

  bmat_->Multiply(true, *tmp_vec_exp, rhs);

  // --- second summand
  tmp_vec = GetStructureGradient(cparams);
  tmp_vec_exp->PutScalar(0.0);
  tmp_vec_exp->Import(*tmp_vec, exporter, Insert);

  Teuchos::RCP<Epetra_Vector> dincr_exp =
      CORE::LINALG::CreateVector(*data_->GSlMaDofRowMapPtr(), true);
  err = dincr_exp->Import(dincr, exporter, Insert);
  if (err) FOUR_C_THROW("Import failed with err = %d", err);

  AssembleGradientBMatrixContribution(*dincr_exp, *tmp_vec_exp, rhs);

  // --- 3rd summand
  AssembleGradientBBMatrixContribution(*dincr_exp, data_->LmN(), rhs);

  Teuchos::RCP<Epetra_Vector> lmincr = CORE::LINALG::CreateVector(data_->GActiveNDofRowMap(), true);

  Teuchos::RCP<CORE::LINALG::SparseMatrix> bbmat =
      CORE::LINALG::MLMultiply(*bmat_, true, *bmat_, false, false, false, true);

  LinSolve(*bbmat, rhs, *lmincr);

  return lmincr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::LagrangeMultiplierFunction::AssembleGradientBBMatrixContribution(
    const Epetra_Vector& dincr, const Epetra_Vector& lm, Epetra_Vector& lmincr) const
{
  if (not lmincr.Map().SameAs(lm.Map())) FOUR_C_THROW("The maps must be identical!");

  for (plain_interface_set::const_iterator cit = interfaces_.begin(); cit != interfaces_.end();
       ++cit)
  {
    const CONTACT::AUG::Interface& interface = dynamic_cast<const CONTACT::AUG::Interface&>(**cit);

    interface.AssembleGradientBBMatrixContribution(dincr, lm, lmincr);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::LagrangeMultiplierFunction::AssembleGradientBMatrixContribution(
    const Epetra_Vector& dincr, const Epetra_Vector& str_grad, Epetra_Vector& lmincr) const
{
  if (not dincr.Map().SameAs(str_grad.Map())) FOUR_C_THROW("The maps must be identical!");

  for (plain_interface_set::const_iterator cit = interfaces_.begin(); cit != interfaces_.end();
       ++cit)
  {
    const CONTACT::AUG::Interface& interface = dynamic_cast<const CONTACT::AUG::Interface&>(**cit);

    interface.AssembleGradientBMatrixContribution(dincr, str_grad, lmincr);
  }
}

FOUR_C_NAMESPACE_CLOSE

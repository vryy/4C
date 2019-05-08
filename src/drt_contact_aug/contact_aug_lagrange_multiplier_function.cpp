/*----------------------------------------------------------------------------*/
/**
  \file contact_aug_lagrange_multiplier_function.cpp

  \brief Lagrange multiplier function: solve a least squares problem to compute
  the Lagrange multiplier value dependent on the current displacement state

  \level 3
  \maintainer Matthias Mayr
  \date Apr 3, 2017
*/
/*----------------------------------------------------------------------------*/

#include "contact_aug_lagrange_multiplier_function.H"
#include "contact_augmented_strategy.H"
#include "contact_augmented_interface.H"

#include "../drt_mortar/mortar_matrix_transform.H"

#include "../drt_contact/contact_paramsinterface.H"

#include "../drt_structure_new/str_model_evaluator_contact.H"

#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_multiply.H"

#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::LagrangeMultiplierFunction::LagrangeMultiplierFunction()
    : isinit_(false),
      issetup_(false),
      strategy_(NULL),
      interfaces_(0),
      data_(Teuchos::null),
      lin_solver_type_(INPAR::SOLVER::undefined),
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
  bmat_ = Teuchos::rcp(new LINALG::SparseMatrix(slMaDofRowMap, 100, false, false));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::Solver> CONTACT::AUG::LagrangeMultiplierFunction::CreateLinearSolver(
    const int lin_sol_id, const Epetra_Comm& comm,
    enum INPAR::SOLVER::SolverType& solver_type) const
{
  if (lin_sol_id == -1) dserror("You must specify a meaningful LINEAR_SOLVER!");

  // get solver parameter list of linear solver
  const Teuchos::ParameterList& solverparams = DRT::Problem::Instance()->SolverParams(lin_sol_id);
  solver_type = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(solverparams, "SOLVER");

  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(
      new LINALG::Solver(solverparams, comm, DRT::Problem::Instance()->ErrorFile()->Handle()));

  if (solver_type != INPAR::SOLVER::umfpack and solver_type != INPAR::SOLVER::superlu)
    dserror("Currently only direct linear solvers are supported!");

  return solver;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::LagrangeMultiplierFunction::LinSolve(
    LINALG::SparseOperator& mat, Epetra_MultiVector& rhs, Epetra_MultiVector& sol)
{
  if (rhs.NumVectors() > 1 or sol.NumVectors() > 1)
    dserror("MultiVector support is not yet implemented!");

  int err = lin_solver_->Solve(
      mat.EpetraOperator(), Teuchos::rcpFromRef(sol), Teuchos::rcpFromRef(rhs), true, true);

  if (err) dserror("LinSolve failed with err = %d", err);
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
  LINALG::Export(*str_gradient, str_gradient_exp);

  Epetra_Vector rhs(data_->GActiveNDofRowMap(), true);
  bmat_->Multiply(true, str_gradient_exp, rhs);

  Teuchos::RCP<LINALG::SparseMatrix> bbmat =
      LINALG::MLMultiply(*bmat_, true, *bmat_, false, false, false, true);

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
  LINALG::SparseMatrix full_stiff(
      *cmodel.GetJacobianBlock(DRT::UTILS::block_displ_displ), LINALG::Copy);

  Teuchos::RCP<LINALG::SparseMatrix> kdd_ptr =
      strategy_->GetMatrixBlockPtr(DRT::UTILS::block_displ_displ);

  // undo matrix contributions
  full_stiff.Add(*kdd_ptr, false, -1.0, 1.0);

  // --- first summand
  Teuchos::RCP<Epetra_Vector> tmp_vec = LINALG::CreateVector(full_stiff.RangeMap(), true);

  int err = full_stiff.Multiply(false, dincr, *tmp_vec);
  if (err) dserror("Multiply failed with err = %d", err);

  Teuchos::RCP<Epetra_Vector> tmp_vec_exp = LINALG::CreateVector(*data_->GSlMaDofRowMapPtr(), true);

  // build necessary exporter
  Epetra_Export exporter(tmp_vec_exp->Map(), tmp_vec->Map());
  err = tmp_vec_exp->Import(*tmp_vec, exporter, Insert);
  if (err) dserror("Import failed with err = %d", err);

  CreateBMatrix();

  bmat_->Multiply(true, *tmp_vec_exp, rhs);

  // --- second summand
  tmp_vec = GetStructureGradient(cparams);
  tmp_vec_exp->Scale(0.0);
  tmp_vec_exp->Import(*tmp_vec, exporter, Insert);

  Teuchos::RCP<Epetra_Vector> dincr_exp = LINALG::CreateVector(*data_->GSlMaDofRowMapPtr(), true);
  err = dincr_exp->Import(dincr, exporter, Insert);
  if (err) dserror("Import failed with err = %d", err);

  AssembleGradientBMatrixContribution(*dincr_exp, *tmp_vec_exp, rhs);

  // --- 3rd summand
  AssembleGradientBBMatrixContribution(*dincr_exp, data_->LmN(), rhs);

  Teuchos::RCP<Epetra_Vector> lmincr = LINALG::CreateVector(data_->GActiveNDofRowMap(), true);

  Teuchos::RCP<LINALG::SparseMatrix> bbmat =
      LINALG::MLMultiply(*bmat_, true, *bmat_, false, false, false, true);

  LinSolve(*bbmat, rhs, *lmincr);

  return lmincr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::LagrangeMultiplierFunction::AssembleGradientBBMatrixContribution(
    const Epetra_Vector& dincr, const Epetra_Vector& lm, Epetra_Vector& lmincr) const
{
  if (not lmincr.Map().SameAs(lm.Map())) dserror("The maps must be identical!");

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
  if (not dincr.Map().SameAs(str_grad.Map())) dserror("The maps must be identical!");

  for (plain_interface_set::const_iterator cit = interfaces_.begin(); cit != interfaces_.end();
       ++cit)
  {
    const CONTACT::AUG::Interface& interface = dynamic_cast<const CONTACT::AUG::Interface&>(**cit);

    interface.AssembleGradientBMatrixContribution(dincr, str_grad, lmincr);
  }
}

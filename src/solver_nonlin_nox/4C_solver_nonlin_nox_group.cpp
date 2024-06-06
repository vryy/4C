/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN implementation of a %::NOX::Epetra::Group
       to handle unconstrained problems.



\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_solver_nonlin_nox_group.hpp"

#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_solver_nonlin_nox_group_prepostoperator.hpp"
#include "4C_solver_nonlin_nox_interface_jacobian.hpp"
#include "4C_solver_nonlin_nox_interface_required.hpp"
#include "4C_solver_nonlin_nox_linearsystem.hpp"
#include "4C_solver_nonlin_nox_solver_ptc.hpp"

#include <NOX_StatusTest_NormF.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Group::Group(Teuchos::ParameterList& printParams, Teuchos::ParameterList& grpOptionParams,
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& i, const ::NOX::Epetra::Vector& x,
    const Teuchos::RCP<::NOX::Epetra::LinearSystem>& linSys)
    : ::NOX::Epetra::Group(printParams, i, x, linSys),
      skipUpdateX_(false),
      corr_type_(NOX::Nln::CorrectionType::vague),
      prePostOperatorPtr_(Teuchos::rcp(new NOX::Nln::GROUP::PrePostOperator(grpOptionParams)))
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Group::Group(const NOX::Nln::Group& source, ::NOX::CopyType type)
    : ::NOX::Epetra::Group(source, type),
      skipUpdateX_(false),
      corr_type_(NOX::Nln::CorrectionType::vague),
      prePostOperatorPtr_(source.prePostOperatorPtr_)
{
  switch (type)
  {
    case ::NOX::DeepCopy:
    {
      skipUpdateX_ = source.skipUpdateX_;
      corr_type_ = source.corr_type_;
      ev_ = source.ev_;
      break;
    }
    default:
      break;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Abstract::Group> NOX::Nln::Group::clone(::NOX::CopyType type) const
{
  Teuchos::RCP<::NOX::Abstract::Group> newgrp = Teuchos::rcp(new NOX::Nln::Group(*this, type));
  return newgrp;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group& NOX::Nln::Group::operator=(const ::NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const ::NOX::Epetra::Group&>(source));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group& NOX::Nln::Group::operator=(const ::NOX::Epetra::Group& source)
{
  ::NOX::Epetra::Group::operator=(source);
  const NOX::Nln::Group& nln_src = dynamic_cast<const NOX::Nln::Group&>(source);

  this->skipUpdateX_ = nln_src.skipUpdateX_;
  this->corr_type_ = nln_src.corr_type_;
  this->ev_ = nln_src.ev_;

  return *this;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Group::resetIsValid()
{
  ::NOX::Epetra::Group::resetIsValid();
  ev_.isvalid_ = false;
  corr_type_ = NOX::Nln::CorrectionType::vague;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Group::computeX(
    const ::NOX::Abstract::Group& grp, const ::NOX::Abstract::Vector& d, double step)
{
  // Cast to appropriate type, then call the "native" computeX
  const NOX::Nln::Group* nlngrp = dynamic_cast<const NOX::Nln::Group*>(&grp);
  if (nlngrp == nullptr) throw_error("computeX", "dyn_cast to nox_nln_group failed!");
  const ::NOX::Epetra::Vector& epetrad = dynamic_cast<const ::NOX::Epetra::Vector&>(d);

  computeX(*nlngrp, epetrad, step);
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Group::computeX(
    const NOX::Nln::Group& grp, const ::NOX::Epetra::Vector& d, double step)
{
  skipUpdateX_ = false;
  prePostOperatorPtr_->runPreComputeX(grp, d.getEpetraVector(), step, *this);

  if (isPreconditioner()) sharedLinearSystem.getObject(this)->destroyPreconditioner();

  resetIsValid();

  if (not skipUpdateX_) xVector.update(1.0, grp.xVector, step, d);

  prePostOperatorPtr_->runPostComputeX(grp, d.getEpetraVector(), step, *this);
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group::ReturnType NOX::Nln::Group::setF(Teuchos::RCP<::NOX::Epetra::Vector> Fptr)
{
  if (Fptr == Teuchos::null or Fptr->getEpetraVector().Map().NumGlobalElements() == 0)
    return ::NOX::Abstract::Group::BadDependency;

  RHSVector = *Fptr;
  isValidRHS = true;

  return ::NOX::Abstract::Group::Ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group::ReturnType NOX::Nln::Group::setJacobianOperator(
    const Teuchos::RCP<const Epetra_Operator> jacOperator)
{
  if (jacOperator == Teuchos::null or jacOperator->OperatorRangeMap().NumGlobalElements() == 0)
    return ::NOX::Abstract::Group::BadDependency;

  sharedLinearSystem.getObject(this)->setJacobianOperatorForSolve(jacOperator);
  isValidJacobian = true;

  return ::NOX::Abstract::Group::Ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Group::resetX() { xVector.init(0.0); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Group::setSkipUpdateX(bool skipUpdateX) { skipUpdateX_ = skipUpdateX; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group::ReturnType NOX::Nln::Group::computeF()
{
  prePostOperatorPtr_->runPreComputeF(RHSVector.getEpetraVector(), *this);

  if (isF()) return ::NOX::Abstract::Group::Ok;

  const bool success = userInterfacePtr->computeF(xVector.getEpetraVector(),
      RHSVector.getEpetraVector(), ::NOX::Epetra::Interface::Required::Residual);

  if (not success)
  {
    throw "NOX::Nln::Group::computeF() - fill failed";
  }

  isValidRHS = true;

  prePostOperatorPtr_->runPostComputeF(RHSVector.getEpetraVector(), *this);
  return ::NOX::Abstract::Group::Ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group::ReturnType NOX::Nln::Group::computeFandJacobian()
{
  // initialize the return type
  ::NOX::Abstract::Group::ReturnType ret = ::NOX::Abstract::Group::Failed;

  // update right hand side vector
  if (!isF() and isJacobian())
  {
    ret = computeF();
  }
  // update right hand side vector and jacobian
  else if (!isJacobian())
  {
    isValidRHS = false;
    prePostOperatorPtr_->runPreComputeF(RHSVector.getEpetraVector(), *this);
    bool status = false;
    Teuchos::RCP<NOX::Nln::LinearSystem> nlnSharedLinearSystem =
        Teuchos::rcp_dynamic_cast<NOX::Nln::LinearSystem>(sharedLinearSystem.getObject(this));

    if (nlnSharedLinearSystem.is_null())
      throw_error("computeFandJacobian", "Dynamic cast of the shared linear system failed!");

    status = nlnSharedLinearSystem->computeFandJacobian(xVector, RHSVector);
    if (!status) throw_error("computeFandJacobian", "evaluation failed!");

    isValidRHS = true;
    isValidJacobian = true;

    ret = ::NOX::Abstract::Group::Ok;
    prePostOperatorPtr_->runPostComputeF(RHSVector.getEpetraVector(), *this);
  }
  // nothing to do, because all quantities are up-to-date
  else
  {
    prePostOperatorPtr_->runPreComputeF(RHSVector.getEpetraVector(), *this);
    ret = ::NOX::Abstract::Group::Ok;
  }

  return ret;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group::ReturnType NOX::Nln::Group::compute_element_volumes(
    Teuchos::RCP<Epetra_Vector>& ele_vols) const
{
  const bool success =
      get_nln_req_interface_ptr()->compute_element_volumes(xVector.getEpetraVector(), ele_vols);

  return (success ? Ok : Failed);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group::ReturnType NOX::Nln::Group::compute_trial_element_volumes(
    Teuchos::RCP<Epetra_Vector>& ele_vols, const ::NOX::Abstract::Vector& dir, double step)
{
  if (tmpVectorPtr.is_null() or !tmpVectorPtr->Map().SameAs(xVector.getEpetraVector().Map()) or
      tmpVectorPtr.get() == &xVector.getEpetraVector())
    tmpVectorPtr = Teuchos::rcp(new Epetra_Vector(xVector.getEpetraVector()));
  else
    tmpVectorPtr->Scale(1.0, xVector.getEpetraVector());

  const ::NOX::Epetra::Vector& dir_epetra = dynamic_cast<const ::NOX::Epetra::Vector&>(dir);
  tmpVectorPtr->Update(step, dir_epetra.getEpetraVector(), 1.0);

  const bool success =
      get_nln_req_interface_ptr()->compute_element_volumes(*tmpVectorPtr, ele_vols);

  return (success ? Ok : Failed);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group::ReturnType NOX::Nln::Group::applyJacobianInverse(Teuchos::ParameterList& p,
    const ::NOX::Epetra::Vector& input, ::NOX::Epetra::Vector& result) const
{
  prePostOperatorPtr_->run_pre_apply_jacobian_inverse(input, result, xVector, *this);

  ::NOX::Abstract::Group::ReturnType status =
      ::NOX::Epetra::Group::applyJacobianInverse(p, input, result);

  prePostOperatorPtr_->run_post_apply_jacobian_inverse(input, result, xVector, *this);

  return status;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group::ReturnType NOX::Nln::Group::applyJacobianBlock(
    const ::NOX::Epetra::Vector& input, Teuchos::RCP<::NOX::Epetra::Vector>& result, unsigned rbid,
    unsigned cbid) const
{
  if (not isJacobian())
    FOUR_C_THROW(
        "It is not possible to access the Jacobian since it has not yet "
        "been evaluated.");

  Teuchos::RCP<NOX::Nln::LinearSystem> nlnSharedLinearSystem =
      Teuchos::rcp_dynamic_cast<NOX::Nln::LinearSystem>(sharedLinearSystem.getObject(this), true);

  const bool success = nlnSharedLinearSystem->applyJacobianBlock(input, result, rbid, cbid);

  return (success ? ::NOX::Abstract::Group::Ok : ::NOX::Abstract::Group::Failed);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group::ReturnType NOX::Nln::Group::compute_correction_system(
    const enum NOX::Nln::CorrectionType type)
{
  prePostOperatorPtr_->runPreComputeF(RHSVector.getEpetraVector(), *this);

  Teuchos::RCP<NOX::Nln::LinearSystem> nlnSharedLinearSystem =
      Teuchos::rcp_dynamic_cast<NOX::Nln::LinearSystem>(sharedLinearSystem.getObject(this));

  if (nlnSharedLinearSystem.is_null())
    FOUR_C_THROW("Dynamic cast of the shared linear system failed!");

  isValidRHS = false;
  isValidJacobian = false;

  const bool success =
      nlnSharedLinearSystem->compute_correction_system(type, *this, xVector, RHSVector);

  if (not success) FOUR_C_THROW("compute_correction_system failed!");

  isValidRHS = true;
  isValidJacobian = true;
  corr_type_ = type;

  prePostOperatorPtr_->runPostComputeF(RHSVector.getEpetraVector(), *this);

  return ::NOX::Abstract::Group::Ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const NOX::Nln::Interface::Required> NOX::Nln::Group::get_nln_req_interface_ptr() const
{
  Teuchos::RCP<NOX::Nln::Interface::Required> userInterfaceNlnPtr =
      Teuchos::rcp_dynamic_cast<NOX::Nln::Interface::Required>(userInterfacePtr);

  if (userInterfaceNlnPtr.is_null())
    throw_error("get_nln_req_interface_ptr",
        "Dynamic cast of the userInterfacePtr to NOX::Nln::Interface::Required failed!");

  return userInterfaceNlnPtr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const std::vector<double>> NOX::Nln::Group::GetRHSNorms(
    const std::vector<::NOX::Abstract::Vector::NormType>& type,
    const std::vector<NOX::Nln::StatusTest::QuantityType>& chQ,
    Teuchos::RCP<const std::vector<::NOX::StatusTest::NormF::ScaleType>> scale) const
{
  if (scale.is_null())
    scale = Teuchos::rcp(new std::vector<::NOX::StatusTest::NormF::ScaleType>(
        chQ.size(), ::NOX::StatusTest::NormF::Unscaled));

  Teuchos::RCP<std::vector<double>> norms = Teuchos::rcp(new std::vector<double>(0));

  double rval = -1.0;
  for (std::size_t i = 0; i < chQ.size(); ++i)
  {
    rval = get_nln_req_interface_ptr()->GetPrimaryRHSNorms(RHSVector.getEpetraVector(), chQ[i],
        type[i], (*scale)[i] == ::NOX::StatusTest::NormF::Scaled);
    if (rval >= 0.0)
    {
      norms->push_back(rval);
    }
    else
    {
      std::ostringstream msg;
      msg << "The desired quantity"
             " for the \"NormF\" Status Test could not be found! (enum="
          << chQ[i] << " | " << NOX::Nln::StatusTest::QuantityType2String(chQ[i]) << ")"
          << std::endl;
      throw_error("GetRHSNorms", msg.str());
    }
  }

  return norms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<double>> NOX::Nln::Group::get_solution_update_rms(
    const ::NOX::Abstract::Vector& xOld, const std::vector<double>& aTol,
    const std::vector<double>& rTol, const std::vector<NOX::Nln::StatusTest::QuantityType>& chQ,
    const std::vector<bool>& disable_implicit_weighting) const
{
  const ::NOX::Epetra::Vector& xOldEpetra = dynamic_cast<const ::NOX::Epetra::Vector&>(xOld);
  Teuchos::RCP<std::vector<double>> rms = Teuchos::rcp(new std::vector<double>(0));

  double rval = -1.0;
  for (std::size_t i = 0; i < chQ.size(); ++i)
  {
    rval = get_nln_req_interface_ptr()->get_primary_solution_update_rms(xVector.getEpetraVector(),
        xOldEpetra.getEpetraVector(), aTol[i], rTol[i], chQ[i], disable_implicit_weighting[i]);
    if (rval >= 0.0)
    {
      rms->push_back(rval);
    }
    else
    {
      std::ostringstream msg;
      msg << "The desired quantity"
             " for the \"NormWRMS\" Status Test could not be found! (enum="
          << chQ[i] << " | " << NOX::Nln::StatusTest::QuantityType2String(chQ[i]) << ")"
          << std::endl;
      throw_error("get_solution_update_rms", msg.str());
    }
  }

  return rms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::Nln::Group::GetTrialUpdateNorm(const ::NOX::Abstract::Vector& dir,
    const ::NOX::Abstract::Vector::NormType normtype, const StatusTest::QuantityType quantity,
    const StatusTest::NormUpdate::ScaleType scale) const
{
  const std::vector<::NOX::Abstract::Vector::NormType> normtypes(1, normtype);
  const std::vector<StatusTest::QuantityType> quantities(1, quantity);
  const std::vector<StatusTest::NormUpdate::ScaleType> scales(1, scale);

  if (tmpVectorPtr.is_null() or !tmpVectorPtr->Map().SameAs(xVector.getEpetraVector().Map()) or
      tmpVectorPtr.get() == &xVector.getEpetraVector())
    tmpVectorPtr = Teuchos::rcp(new Epetra_Vector(xVector.getEpetraVector()));
  else
    tmpVectorPtr->Scale(1.0, xVector.getEpetraVector());

  // change the internally stored x-vector for the norm evaluation
  ::NOX::Epetra::Vector& x_mutable = const_cast<::NOX::Epetra::Vector&>(xVector);
  x_mutable.update(1.0, dir, 1.0);

  ::NOX::Epetra::Vector xold(tmpVectorPtr, ::NOX::Epetra::Vector::CreateView);

  const double rval =
      get_solution_update_norms(xold, normtypes, quantities, Teuchos::rcpFromRef(scales))->at(0);

  // un-do the changes to the x-vector
  x_mutable.getEpetraVector().Scale(1.0, *tmpVectorPtr);

  return rval;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<double>> NOX::Nln::Group::get_solution_update_norms(
    const ::NOX::Abstract::Vector& xOld, const std::vector<::NOX::Abstract::Vector::NormType>& type,
    const std::vector<StatusTest::QuantityType>& chQ,
    Teuchos::RCP<const std::vector<StatusTest::NormUpdate::ScaleType>> scale) const
{
  const ::NOX::Epetra::Vector& xOldEpetra = dynamic_cast<const ::NOX::Epetra::Vector&>(xOld);
  if (scale.is_null())
    scale = Teuchos::rcp(new std::vector<StatusTest::NormUpdate::ScaleType>(
        chQ.size(), StatusTest::NormUpdate::Unscaled));

  Teuchos::RCP<std::vector<double>> norms = Teuchos::rcp(new std::vector<double>(0));

  double rval = -1.0;
  for (std::size_t i = 0; i < chQ.size(); ++i)
  {
    rval = get_nln_req_interface_ptr()->get_primary_solution_update_norms(xVector.getEpetraVector(),
        xOldEpetra.getEpetraVector(), chQ[i], type[i],
        (*scale)[i] == StatusTest::NormUpdate::Scaled);
    if (rval >= 0.0)
    {
      norms->push_back(rval);
    }
    else
    {
      std::ostringstream msg;
      msg << "The desired quantity"
             " for the \"NormIncr\" Status Test could not be found! (enum="
          << chQ[i] << " | " << NOX::Nln::StatusTest::QuantityType2String(chQ[i]) << ")"
          << std::endl;
      throw_error("get_solution_update_norms", msg.str());
    }
  }

  return norms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<double>> NOX::Nln::Group::get_previous_solution_norms(
    const ::NOX::Abstract::Vector& xOld, const std::vector<::NOX::Abstract::Vector::NormType>& type,
    const std::vector<StatusTest::QuantityType>& chQ,
    Teuchos::RCP<const std::vector<StatusTest::NormUpdate::ScaleType>> scale) const
{
  const ::NOX::Epetra::Vector& xOldEpetra = dynamic_cast<const ::NOX::Epetra::Vector&>(xOld);
  if (scale.is_null())
    scale = Teuchos::rcp(new std::vector<StatusTest::NormUpdate::ScaleType>(
        chQ.size(), StatusTest::NormUpdate::Unscaled));

  Teuchos::RCP<std::vector<double>> norms = Teuchos::rcp(new std::vector<double>(0));

  double rval = -1.0;
  for (std::size_t i = 0; i < chQ.size(); ++i)
  {
    rval = get_nln_req_interface_ptr()->get_previous_primary_solution_norms(
        xOldEpetra.getEpetraVector(), chQ[i], type[i],
        (*scale)[i] == StatusTest::NormUpdate::Scaled);
    if (rval >= 0.0)
    {
      norms->push_back(rval);
    }
    else
    {
      std::ostringstream msg;
      msg << "The desired quantity"
             " for the \"NormUpdate\" Status Test could not be found! (enum="
          << chQ[i] << " | " << NOX::Nln::StatusTest::QuantityType2String(chQ[i]) << ")"
          << std::endl;
      throw_error("get_previous_solution_norms", msg.str());
    }
  }

  return norms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::Nln::Group::GetModelValue(const enum MeritFunction::MeritFctName merit_func_type) const
{
  return get_nln_req_interface_ptr()->GetModelValue(
      xVector.getEpetraVector(), RHSVector.getEpetraVector(), merit_func_type);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::Nln::Group::get_linearized_model_terms(const ::NOX::Abstract::Vector& dir,
    const enum NOX::Nln::MeritFunction::MeritFctName mf_type,
    const enum NOX::Nln::MeritFunction::LinOrder linorder,
    const enum NOX::Nln::MeritFunction::LinType lintype) const
{
  const ::NOX::Epetra::Vector& dir_nox_epetra = dynamic_cast<const ::NOX::Epetra::Vector&>(dir);
  const Epetra_Vector& dir_epetra = dir_nox_epetra.getEpetraVector();

  return get_nln_req_interface_ptr()->get_linearized_model_terms(
      this, dir_epetra, mf_type, linorder, lintype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Group::reset_pre_post_operator(
    Teuchos::ParameterList& grpOptionParams, const bool& resetIsValidFlags)
{
  if (resetIsValidFlags) resetIsValid();

  prePostOperatorPtr_->reset(grpOptionParams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Group::reset_lin_sys_pre_post_operator(
    Teuchos::ParameterList& linearSolverParams, const bool& resetIsValidFlags)
{
  if (resetIsValidFlags) resetIsValid();

  Teuchos::RCP<NOX::Nln::LinearSystem> nlnLinsysPtr =
      Teuchos::rcp_dynamic_cast<NOX::Nln::LinearSystem>(getLinearSystem());

  if (nlnLinsysPtr.is_null())
    throw_error("reset_lin_sys_pre_post_operator", "The linear system cast failed!");

  nlnLinsysPtr->reset_pre_post_operator(linearSolverParams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Group::adjust_pseudo_time_step(double& delta, const double& stepSize,
    const ::NOX::Abstract::Vector& dir, const NOX::Nln::Solver::PseudoTransient& ptcsolver)
{
  const ::NOX::Epetra::Vector& dirEpetra = dynamic_cast<const ::NOX::Epetra::Vector&>(dir);
  adjust_pseudo_time_step(delta, stepSize, dirEpetra, ptcsolver);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Group::adjust_pseudo_time_step(double& delta, const double& stepSize,
    const ::NOX::Epetra::Vector& dir, const NOX::Nln::Solver::PseudoTransient& ptcsolver)
{
  if (!isF() or !isJacobian())
    throw_error("AdjustPseudoTimeStep", "F and/or the jacobian are not evaluated!");

  Teuchos::RCP<NOX::Nln::LinearSystem> nlnSharedLinearSystem =
      Teuchos::rcp_dynamic_cast<NOX::Nln::LinearSystem>(sharedLinearSystem.getObject(this));

  if (nlnSharedLinearSystem.is_null())
    throw_error("AdjustPseudoTimeStep()", "Dynamic cast of the shared linear system failed!");

  nlnSharedLinearSystem->adjust_pseudo_time_step(delta, stepSize, dir, RHSVector, ptcsolver);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> NOX::Nln::Group::get_lumped_mass_matrix_ptr() const
{
  return Teuchos::rcp_dynamic_cast<NOX::Nln::Interface::Required>(userInterfacePtr)
      ->get_lumped_mass_matrix_ptr();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseMatrix> NOX::Nln::Group::get_contributions_from_element_level()
{
  return Teuchos::rcp_dynamic_cast<NOX::Nln::Interface::Jacobian>(userInterfacePtr)
      ->calc_jacobian_contributions_from_element_level_for_ptc();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Group::DestroyJacobian()
{
  Teuchos::RCP<NOX::Nln::LinearSystem> nlnSharedLinearSystem =
      Teuchos::rcp_dynamic_cast<NOX::Nln::LinearSystem>(sharedLinearSystem.getObject(this));

  isValidJacobian = false;

  return nlnSharedLinearSystem->DestroyJacobian();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Group::isJacobian() const
{
  if (isValidJacobian and not sharedLinearSystem.isOwner(this)) sharedLinearSystem.getObject(this);

  return ::NOX::Epetra::Group::isJacobian();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Group::throw_error(
    const std::string& functionName, const std::string& errorMsg) const
{
  std::ostringstream msg;
  msg << "ERROR - NOX::Nln::Group::" << functionName << " - " << errorMsg << std::endl;

  FOUR_C_THROW(msg.str());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_BlockMap& NOX::Nln::Group::getDofMap() const
{
  return xVector.getEpetraVector().Map();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Group::CreateBackupState(const ::NOX::Abstract::Vector& dir) const
{
  Teuchos::RCP<NOX::Nln::Interface::Required> nln_required =
      Teuchos::rcp_dynamic_cast<NOX::Nln::Interface::Required>(userInterfacePtr, true);

  const ::NOX::Epetra::Vector* epetra_dir = dynamic_cast<const ::NOX::Epetra::Vector*>(&dir);
  if (not epetra_dir) FOUR_C_THROW("Dynamic cast failed.");

  nln_required->CreateBackupState(epetra_dir->getEpetraVector());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Group::recover_from_backup_state()
{
  Teuchos::RCP<NOX::Nln::Interface::Required> nln_required =
      Teuchos::rcp_dynamic_cast<NOX::Nln::Interface::Required>(userInterfacePtr);

  if (nln_required.is_null()) FOUR_C_THROW("Dynamic cast failed.");

  nln_required->recover_from_backup_state();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map& NOX::Nln::Group::getJacobianRangeMap(unsigned rbid, unsigned cbid) const
{
  if (not isJacobian())
    FOUR_C_THROW(
        "It is not possible to access the Jacobian since it has not yet "
        "been evaluated.");

  Teuchos::RCP<NOX::Nln::LinearSystem> nlnSharedLinearSystem =
      Teuchos::rcp_dynamic_cast<NOX::Nln::LinearSystem>(sharedLinearSystem.getObject(this), true);

  return nlnSharedLinearSystem->getJacobianRangeMap(rbid, cbid);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> NOX::Nln::Group::get_diagonal_of_jacobian(unsigned diag_bid) const
{
  if (not isJacobian())
    FOUR_C_THROW(
        "It is not possible to access the Jacobian since it has not yet "
        "been evaluated.");

  Teuchos::RCP<NOX::Nln::LinearSystem> nlnSharedLinearSystem =
      Teuchos::rcp_dynamic_cast<NOX::Nln::LinearSystem>(sharedLinearSystem.getObject(this), true);

  return nlnSharedLinearSystem->get_diagonal_of_jacobian(diag_bid);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Group::replace_diagonal_of_jacobian(
    const Epetra_Vector& new_diag, unsigned diag_bid) const
{
  if (not isJacobian())
    FOUR_C_THROW(
        "It is not possible to access the Jacobian since it has not yet "
        "been evaluated.");

  Teuchos::RCP<NOX::Nln::LinearSystem> nlnSharedLinearSystem =
      Teuchos::rcp_dynamic_cast<NOX::Nln::LinearSystem>(sharedLinearSystem.getObject(this), true);

  return nlnSharedLinearSystem->replace_diagonal_of_jacobian(new_diag, diag_bid);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group::ReturnType NOX::Nln::Group::compute_serial_jacobian_condition_number(
    const NOX::Nln::LinSystem::ConditionNumber condnum_type, bool printOutput)
{
  if (isConditionNumber()) return ::NOX::Abstract::Group::Ok;

  if (xVector.getEpetraVector().Comm().NumProc() != 1)
    FOUR_C_THROW("Only serial mode is supported!");

  if (!isJacobian()) FOUR_C_THROW("Jacobian is invalid wrt the solution.");

  switch (condnum_type)
  {
    case LinSystem::ConditionNumber::max_min_ev_ratio:
    {
      compute_serial_jacobian_eigenvalues(printOutput);

      if (ev_.real_min_ == 0)
      {
        utils.out(::NOX::Utils::Warning) << "Jacobian is singular!" << std::endl;
        conditionNumber = std::numeric_limits<double>::infinity();
      }
      else
        conditionNumber = ev_.real_max_ / ev_.real_min_;

      break;
    }
    case LinSystem::ConditionNumber::one_norm:
    case LinSystem::ConditionNumber::inf_norm:
    {
      Teuchos::RCP<NOX::Nln::LinearSystem> nlnSharedLinearSystem =
          Teuchos::rcp_dynamic_cast<NOX::Nln::LinearSystem>(
              sharedLinearSystem.getObject(this), true);

      conditionNumber =
          nlnSharedLinearSystem->compute_serial_condition_number_of_jacobian(condnum_type);

      break;
    }
    default:
      FOUR_C_THROW("Unknown LinSystem::ConditionNumber type!");
      exit(EXIT_FAILURE);
  }

  if (printOutput)
    utils.out() << "Condition number = " << ::NOX::Utils::sciformat(conditionNumber, 5)
                << std::endl;

  isValidConditionNumber = true;

  return ::NOX::Abstract::Group::Ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group::ReturnType NOX::Nln::Group::compute_serial_jacobian_eigenvalues(
    bool printOutput)
{
  if (xVector.getEpetraVector().Comm().NumProc() != 1) FOUR_C_THROW("Works only in serial mode!");

  if (ev_.isvalid_) return ::NOX::Abstract::Group::Ok;

  Teuchos::RCP<NOX::Nln::LinearSystem> nlnSharedLinearSystem =
      Teuchos::rcp_dynamic_cast<NOX::Nln::LinearSystem>(sharedLinearSystem.getObject(this), true);

  nlnSharedLinearSystem->compute_serial_eigenvalues_of_jacobian(ev_.realpart_, ev_.imaginarypart_);

  const int length = ev_.realpart_.length();
  ev_.real_max_ = *std::max_element(ev_.realpart_.values(), ev_.realpart_.values() + length);
  ev_.real_min_ = *std::min_element(ev_.realpart_.values(), ev_.realpart_.values() + length);

  if (printOutput)
  {
    utils.out() << "maximal eigenvalue = " << utils.sciformat(ev_.real_max_, 5) << "\n";
    utils.out() << "minimal eigenvalue = " << utils.sciformat(ev_.real_min_, 5) << "\n";
  }

  ev_.isvalid_ = true;

  return ::NOX::Abstract::Group::Ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Group::isEigenvalues() const { return ev_.isvalid_; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Core::LinAlg::SerialDenseVector& NOX::Nln::Group::get_jacobian_real_eigenvalues() const
{
  if (not isEigenvalues()) FOUR_C_THROW("The eigenvalues has not yet been computed!");
  return ev_.realpart_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Core::LinAlg::SerialDenseVector& NOX::Nln::Group::get_jacobian_imaginary_eigenvalues() const
{
  if (not isEigenvalues()) FOUR_C_THROW("The eigenvalues has not yet been computed!");
  return ev_.imaginarypart_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::Nln::Group::get_jacobian_max_real_eigenvalue() const
{
  if (not isEigenvalues()) FOUR_C_THROW("The eigenvalues has not yet been computed!");
  return ev_.real_max_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::Nln::Group::get_jacobian_min_real_eigenvalue() const
{
  if (not isEigenvalues()) FOUR_C_THROW("The eigenvalues has not yet been computed!");
  return ev_.real_min_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Group::Eigenvalues& NOX::Nln::Group::Eigenvalues::operator=(const Eigenvalues& src)
{
  if (not src.isvalid_)
  {
    this->isvalid_ = false;
    return *this;
  }

  this->realpart_.resize(src.realpart_.length());
  this->realpart_ = src.realpart_;

  this->imaginarypart_.resize(src.imaginarypart_.length());
  this->imaginarypart_ = src.imaginarypart_;

  this->real_max_ = src.real_max_;
  this->real_min_ = src.real_min_;

  this->isvalid_ = src.isvalid_;

  return *this;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Group::getDofsFromElements(
    const std::vector<int>& my_ele_gids, std::set<int>& my_ele_dofs) const
{
  get_nln_req_interface_ptr()->getDofsFromElements(my_ele_gids, my_ele_dofs);
}

FOUR_C_NAMESPACE_CLOSE

/*-----------------------------------------------------------*/
/*!
\file nox_nln_group.cpp

\brief %NOX::NLN implementation of a %NOX::Epetra::Group
       to handle unconstrained problems.

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_group.H"
#include "nox_nln_interface_required.H"
#include "nox_nln_interface_jacobian.H"
#include "nox_nln_linearsystem.H"
#include "nox_nln_group_prepostoperator.H"
#include "nox_nln_solver_ptc.H"

#include <NOX_StatusTest_NormF.H>

#include "../linalg/linalg_utils.H"
#include "../solver/solver_aztecoo_conditionnumber.H"
#include <az_aztec_defs.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Group::Group(Teuchos::ParameterList& printParams, Teuchos::ParameterList& grpOptionParams,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& i, const NOX::Epetra::Vector& x,
    const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys)
    : NOX::Epetra::Group(printParams, i, x, linSys),
      skipUpdateX_(false),
      corr_type_(NOX::NLN::CorrectionType::vague),
      prePostOperatorPtr_(Teuchos::rcp(new NOX::NLN::GROUP::PrePostOperator(grpOptionParams)))
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Group::Group(const NOX::NLN::Group& source, CopyType type)
    : NOX::Epetra::Group(source, type),
      skipUpdateX_(false),
      corr_type_(NOX::NLN::CorrectionType::vague),
      prePostOperatorPtr_(source.prePostOperatorPtr_)
{
  switch (type)
  {
    case DeepCopy:
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
Teuchos::RCP<NOX::Abstract::Group> NOX::NLN::Group::clone(CopyType type) const
{
  Teuchos::RCP<NOX::Abstract::Group> newgrp = Teuchos::rcp(new NOX::NLN::Group(*this, type));
  return newgrp;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group& NOX::NLN::Group::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const NOX::Epetra::Group&>(source));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group& NOX::NLN::Group::operator=(const NOX::Epetra::Group& source)
{
  NOX::Epetra::Group::operator=(source);
  const NOX::NLN::Group& nln_src = dynamic_cast<const NOX::NLN::Group&>(source);

  this->skipUpdateX_ = nln_src.skipUpdateX_;
  this->corr_type_ = nln_src.corr_type_;
  this->ev_ = nln_src.ev_;

  return *this;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Group::resetIsValid()
{
  NOX::Epetra::Group::resetIsValid();
  ev_.isvalid_ = false;
  corr_type_ = NOX::NLN::CorrectionType::vague;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Group::computeX(
    const NOX::Abstract::Group& grp, const NOX::Abstract::Vector& d, double step)
{
  // Cast to appropriate type, then call the "native" computeX
  const NOX::NLN::Group* nlngrp = dynamic_cast<const NOX::NLN::Group*>(&grp);
  if (nlngrp == NULL) throwError("computeX", "dyn_cast to nox_nln_group failed!");
  const NOX::Epetra::Vector& epetrad = dynamic_cast<const NOX::Epetra::Vector&>(d);

  computeX(*nlngrp, epetrad, step);
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Group::computeX(
    const NOX::NLN::Group& grp, const NOX::Epetra::Vector& d, double step)
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
NOX::Abstract::Group::ReturnType NOX::NLN::Group::setF(Teuchos::RCP<NOX::Epetra::Vector> Fptr)
{
  if (Fptr == Teuchos::null or Fptr->getEpetraVector().Map().NumGlobalElements() == 0)
    return NOX::Abstract::Group::BadDependency;

  RHSVector = *Fptr;
  isValidRHS = true;

  return NOX::Abstract::Group::Ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::NLN::Group::setJacobianOperator(
    const Teuchos::RCP<const Epetra_Operator> jacOperator)
{
  if (jacOperator == Teuchos::null or jacOperator->OperatorRangeMap().NumGlobalElements() == 0)
    return NOX::Abstract::Group::BadDependency;

  sharedLinearSystem.getObject(this)->setJacobianOperatorForSolve(jacOperator);
  isValidJacobian = true;

  return NOX::Abstract::Group::Ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Group::setSkipUpdateX(bool skipUpdateX) { skipUpdateX_ = skipUpdateX; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::NLN::Group::computeF()
{
  prePostOperatorPtr_->runPreComputeF(RHSVector.getEpetraVector(), *this);

  if (isF()) return NOX::Abstract::Group::Ok;

  const bool success = userInterfacePtr->computeF(xVector.getEpetraVector(),
      RHSVector.getEpetraVector(), NOX::Epetra::Interface::Required::Residual);

  if (not success)
  {
    throw "NOX::NLN::Group::computeF() - fill failed";
  }

  isValidRHS = true;

  prePostOperatorPtr_->runPostComputeF(RHSVector.getEpetraVector(), *this);
  return NOX::Abstract::Group::Ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::NLN::Group::computeFandJacobian()
{
  // initialize the return type
  NOX::Abstract::Group::ReturnType ret = NOX::Abstract::Group::Failed;

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
    Teuchos::RCP<NOX::NLN::LinearSystem> nlnSharedLinearSystem =
        Teuchos::rcp_dynamic_cast<NOX::NLN::LinearSystem>(sharedLinearSystem.getObject(this));

    if (nlnSharedLinearSystem.is_null())
      throwError("computeFandJacobian", "Dynamic cast of the shared linear system failed!");

    status = nlnSharedLinearSystem->computeFandJacobian(xVector, RHSVector);
    if (!status) throwError("computeFandJacobian", "evaluation failed!");

    isValidRHS = true;
    isValidJacobian = true;

    ret = NOX::Abstract::Group::Ok;
    prePostOperatorPtr_->runPostComputeF(RHSVector.getEpetraVector(), *this);
  }
  // nothing to do, because all quantities are up-to-date
  else
  {
    prePostOperatorPtr_->runPreComputeF(RHSVector.getEpetraVector(), *this);
    ret = NOX::Abstract::Group::Ok;
  }

  return ret;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::NLN::Group::computeElementVolumes(
    Teuchos::RCP<Epetra_Vector>& ele_vols) const
{
  const bool success =
      GetNlnReqInterfacePtr()->computeElementVolumes(xVector.getEpetraVector(), ele_vols);

  return (success ? Ok : Failed);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::NLN::Group::computeTrialElementVolumes(
    Teuchos::RCP<Epetra_Vector>& ele_vols, const NOX::Abstract::Vector& dir, double step)
{
  if (tmpVectorPtr.is_null() or !tmpVectorPtr->Map().SameAs(xVector.getEpetraVector().Map()) or
      tmpVectorPtr.get() == &xVector.getEpetraVector())
    tmpVectorPtr = Teuchos::rcp(new Epetra_Vector(xVector.getEpetraVector()));
  else
    tmpVectorPtr->Scale(1.0, xVector.getEpetraVector());

  const NOX::Epetra::Vector& dir_epetra = dynamic_cast<const NOX::Epetra::Vector&>(dir);
  tmpVectorPtr->Update(step, dir_epetra.getEpetraVector(), 1.0);

  const bool success = GetNlnReqInterfacePtr()->computeElementVolumes(*tmpVectorPtr, ele_vols);

  return (success ? Ok : Failed);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::NLN::Group::applyJacobianInverse(
    Teuchos::ParameterList& p, const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const
{
  prePostOperatorPtr_->runPreApplyJacobianInverse(input, result, xVector, *this);

  NOX::Abstract::Group::ReturnType status =
      NOX::Epetra::Group::applyJacobianInverse(p, input, result);

  prePostOperatorPtr_->runPostApplyJacobianInverse(input, result, xVector, *this);

  return status;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::NLN::Group::applyJacobianBlock(
    const NOX::Epetra::Vector& input, Teuchos::RCP<NOX::Epetra::Vector>& result, unsigned rbid,
    unsigned cbid) const
{
  if (not isJacobian())
    dserror(
        "It is not possible to access the Jacobian since it has not yet "
        "been evaluated.");

  Teuchos::RCP<NOX::NLN::LinearSystem> nlnSharedLinearSystem =
      Teuchos::rcp_dynamic_cast<NOX::NLN::LinearSystem>(sharedLinearSystem.getObject(this), true);

  const bool success = nlnSharedLinearSystem->applyJacobianBlock(input, result, rbid, cbid);

  return (success ? NOX::Abstract::Group::Ok : NOX::Abstract::Group::Failed);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::NLN::Group::computeCorrectionSystem(
    const enum NOX::NLN::CorrectionType type)
{
  prePostOperatorPtr_->runPreComputeF(RHSVector.getEpetraVector(), *this);

  Teuchos::RCP<NOX::NLN::LinearSystem> nlnSharedLinearSystem =
      Teuchos::rcp_dynamic_cast<NOX::NLN::LinearSystem>(sharedLinearSystem.getObject(this));

  if (nlnSharedLinearSystem.is_null()) dserror("Dynamic cast of the shared linear system failed!");

  isValidRHS = false;
  isValidJacobian = false;

  const bool success =
      nlnSharedLinearSystem->computeCorrectionSystem(type, *this, xVector, RHSVector);

  if (not success) run_time_error("computeCorrectionSystem failed!");

  isValidRHS = true;
  isValidJacobian = true;
  corr_type_ = type;

  prePostOperatorPtr_->runPostComputeF(RHSVector.getEpetraVector(), *this);

  return NOX::Abstract::Group::Ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const NOX::NLN::Interface::Required> NOX::NLN::Group::GetNlnReqInterfacePtr() const
{
  Teuchos::RCP<NOX::NLN::Interface::Required> userInterfaceNlnPtr =
      Teuchos::rcp_dynamic_cast<NOX::NLN::Interface::Required>(userInterfacePtr);

  if (userInterfaceNlnPtr.is_null())
    throwError("GetNlnReqInterfacePtr",
        "Dynamic cast of the userInterfacePtr to NOX::NLN::Interface::Required failed!");

  return userInterfaceNlnPtr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const std::vector<double>> NOX::NLN::Group::GetRHSNorms(
    const std::vector<NOX::Abstract::Vector::NormType>& type,
    const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
    Teuchos::RCP<const std::vector<NOX::StatusTest::NormF::ScaleType>> scale) const
{
  if (scale.is_null())
    scale = Teuchos::rcp(new std::vector<NOX::StatusTest::NormF::ScaleType>(
        chQ.size(), NOX::StatusTest::NormF::Unscaled));

  Teuchos::RCP<std::vector<double>> norms = Teuchos::rcp(new std::vector<double>(0));

  double rval = -1.0;
  for (std::size_t i = 0; i < chQ.size(); ++i)
  {
    rval = GetNlnReqInterfacePtr()->GetPrimaryRHSNorms(RHSVector.getEpetraVector(), chQ[i], type[i],
        (*scale)[i] == NOX::StatusTest::NormF::Scaled);
    if (rval >= 0.0)
    {
      norms->push_back(rval);
    }
    else
    {
      std::ostringstream msg;
      msg << "The desired quantity"
             " for the \"NormF\" Status Test could not be found! (enum="
          << chQ[i] << " | " << NOX::NLN::StatusTest::QuantityType2String(chQ[i]) << ")"
          << std::endl;
      throwError("GetRHSNorms", msg.str());
    }
  }

  return norms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<double>> NOX::NLN::Group::GetSolutionUpdateRMS(
    const NOX::Abstract::Vector& xOld, const std::vector<double>& aTol,
    const std::vector<double>& rTol, const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
    const std::vector<bool>& disable_implicit_weighting) const
{
  const NOX::Epetra::Vector& xOldEpetra = dynamic_cast<const NOX::Epetra::Vector&>(xOld);
  Teuchos::RCP<std::vector<double>> rms = Teuchos::rcp(new std::vector<double>(0));

  double rval = -1.0;
  for (std::size_t i = 0; i < chQ.size(); ++i)
  {
    rval = GetNlnReqInterfacePtr()->GetPrimarySolutionUpdateRMS(xVector.getEpetraVector(),
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
          << chQ[i] << " | " << NOX::NLN::StatusTest::QuantityType2String(chQ[i]) << ")"
          << std::endl;
      throwError("GetSolutionUpdateRMS", msg.str());
    }
  }

  return rms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::Group::GetTrialUpdateNorm(const NOX::Abstract::Vector& dir,
    const NOX::Abstract::Vector::NormType normtype, const StatusTest::QuantityType quantity,
    const StatusTest::NormUpdate::ScaleType scale) const
{
  const std::vector<NOX::Abstract::Vector::NormType> normtypes(1, normtype);
  const std::vector<StatusTest::QuantityType> quantities(1, quantity);
  const std::vector<StatusTest::NormUpdate::ScaleType> scales(1, scale);

  if (tmpVectorPtr.is_null() or !tmpVectorPtr->Map().SameAs(xVector.getEpetraVector().Map()) or
      tmpVectorPtr.get() == &xVector.getEpetraVector())
    tmpVectorPtr = Teuchos::rcp(new Epetra_Vector(xVector.getEpetraVector()));
  else
    tmpVectorPtr->Scale(1.0, xVector.getEpetraVector());

  // change the internally stored x-vector for the norm evaluation
  NOX::Epetra::Vector& x_mutable = const_cast<NOX::Epetra::Vector&>(xVector);
  x_mutable.update(1.0, dir, 1.0);

  NOX::Epetra::Vector xold(tmpVectorPtr, NOX::Epetra::Vector::CreateView);

  const double rval =
      GetSolutionUpdateNorms(xold, normtypes, quantities, Teuchos::rcpFromRef(scales))->at(0);

  // un-do the changes to the x-vector
  x_mutable.getEpetraVector().Scale(1.0, *tmpVectorPtr);

  return rval;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<double>> NOX::NLN::Group::GetSolutionUpdateNorms(
    const NOX::Abstract::Vector& xOld, const std::vector<NOX::Abstract::Vector::NormType>& type,
    const std::vector<StatusTest::QuantityType>& chQ,
    Teuchos::RCP<const std::vector<StatusTest::NormUpdate::ScaleType>> scale) const
{
  const NOX::Epetra::Vector& xOldEpetra = dynamic_cast<const NOX::Epetra::Vector&>(xOld);
  if (scale.is_null())
    scale = Teuchos::rcp(new std::vector<StatusTest::NormUpdate::ScaleType>(
        chQ.size(), StatusTest::NormUpdate::Unscaled));

  Teuchos::RCP<std::vector<double>> norms = Teuchos::rcp(new std::vector<double>(0));

  double rval = -1.0;
  for (std::size_t i = 0; i < chQ.size(); ++i)
  {
    rval = GetNlnReqInterfacePtr()->GetPrimarySolutionUpdateNorms(xVector.getEpetraVector(),
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
          << chQ[i] << " | " << NOX::NLN::StatusTest::QuantityType2String(chQ[i]) << ")"
          << std::endl;
      throwError("GetSolutionUpdateNorms", msg.str());
    }
  }

  return norms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<double>> NOX::NLN::Group::GetPreviousSolutionNorms(
    const NOX::Abstract::Vector& xOld, const std::vector<NOX::Abstract::Vector::NormType>& type,
    const std::vector<StatusTest::QuantityType>& chQ,
    Teuchos::RCP<const std::vector<StatusTest::NormUpdate::ScaleType>> scale) const
{
  const NOX::Epetra::Vector& xOldEpetra = dynamic_cast<const NOX::Epetra::Vector&>(xOld);
  if (scale.is_null())
    scale = Teuchos::rcp(new std::vector<StatusTest::NormUpdate::ScaleType>(
        chQ.size(), StatusTest::NormUpdate::Unscaled));

  Teuchos::RCP<std::vector<double>> norms = Teuchos::rcp(new std::vector<double>(0));

  double rval = -1.0;
  for (std::size_t i = 0; i < chQ.size(); ++i)
  {
    rval = GetNlnReqInterfacePtr()->GetPreviousPrimarySolutionNorms(xOldEpetra.getEpetraVector(),
        chQ[i], type[i], (*scale)[i] == StatusTest::NormUpdate::Scaled);
    if (rval >= 0.0)
    {
      norms->push_back(rval);
    }
    else
    {
      std::ostringstream msg;
      msg << "The desired quantity"
             " for the \"NormUpdate\" Status Test could not be found! (enum="
          << chQ[i] << " | " << NOX::NLN::StatusTest::QuantityType2String(chQ[i]) << ")"
          << std::endl;
      throwError("GetPreviousSolutionNorms", msg.str());
    }
  }

  return norms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::Group::GetModelValue(const enum MeritFunction::MeritFctName merit_func_type) const
{
  return GetNlnReqInterfacePtr()->GetModelValue(
      xVector.getEpetraVector(), RHSVector.getEpetraVector(), merit_func_type);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::Group::GetLinearizedModelTerms(const NOX::Abstract::Vector& dir,
    const enum NOX::NLN::MeritFunction::MeritFctName mf_type,
    const enum NOX::NLN::MeritFunction::LinOrder linorder,
    const enum NOX::NLN::MeritFunction::LinType lintype) const
{
  const NOX::Epetra::Vector& dir_nox_epetra = dynamic_cast<const NOX::Epetra::Vector&>(dir);
  const Epetra_Vector& dir_epetra = dir_nox_epetra.getEpetraVector();

  return GetNlnReqInterfacePtr()->GetLinearizedModelTerms(
      this, dir_epetra, mf_type, linorder, lintype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Group::ResetPrePostOperator(
    Teuchos::ParameterList& grpOptionParams, const bool& resetIsValidFlags)
{
  if (resetIsValidFlags) resetIsValid();

  prePostOperatorPtr_->reset(grpOptionParams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Group::ResetLinSysPrePostOperator(
    Teuchos::ParameterList& linearSolverParams, const bool& resetIsValidFlags)
{
  if (resetIsValidFlags) resetIsValid();

  Teuchos::RCP<NOX::NLN::LinearSystem> nlnLinsysPtr =
      Teuchos::rcp_dynamic_cast<NOX::NLN::LinearSystem>(getLinearSystem());

  if (nlnLinsysPtr.is_null())
    throwError("ResetLinSysPrePostOperator", "The linear system cast failed!");

  nlnLinsysPtr->resetPrePostOperator(linearSolverParams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Group::adjustPseudoTimeStep(double& delta, const double& stepSize,
    const NOX::Abstract::Vector& dir, const NOX::NLN::Solver::PseudoTransient& ptcsolver)
{
  const NOX::Epetra::Vector& dirEpetra = dynamic_cast<const NOX::Epetra::Vector&>(dir);
  adjustPseudoTimeStep(delta, stepSize, dirEpetra, ptcsolver);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Group::adjustPseudoTimeStep(double& delta, const double& stepSize,
    const NOX::Epetra::Vector& dir, const NOX::NLN::Solver::PseudoTransient& ptcsolver)
{
  if (!isF() or !isJacobian())
    throwError("AdjustPseudoTimeStep", "F and/or the jacobian are not evaluated!");

  Teuchos::RCP<NOX::NLN::LinearSystem> nlnSharedLinearSystem =
      Teuchos::rcp_dynamic_cast<NOX::NLN::LinearSystem>(sharedLinearSystem.getObject(this));

  if (nlnSharedLinearSystem.is_null())
    throwError("AdjustPseudoTimeStep()", "Dynamic cast of the shared linear system failed!");

  nlnSharedLinearSystem->adjustPseudoTimeStep(delta, stepSize, dir, RHSVector, ptcsolver);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> NOX::NLN::Group::GetLumpedMassMatrixPtr() const
{
  return Teuchos::rcp_dynamic_cast<NOX::NLN::Interface::Required>(userInterfacePtr)
      ->GetLumpedMassMatrixPtr();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> NOX::NLN::Group::GetContributionsFromElementLevel()
{
  return Teuchos::rcp_dynamic_cast<NOX::NLN::Interface::Jacobian>(userInterfacePtr)
      ->CalcJacobianContributionsFromElementLevelForPTC();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Group::DestroyJacobian()
{
  Teuchos::RCP<NOX::NLN::LinearSystem> nlnSharedLinearSystem =
      Teuchos::rcp_dynamic_cast<NOX::NLN::LinearSystem>(sharedLinearSystem.getObject(this));

  isValidJacobian = false;

  return nlnSharedLinearSystem->DestroyJacobian();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Group::isJacobian() const
{
  if (isValidJacobian and not sharedLinearSystem.isOwner(this)) sharedLinearSystem.getObject(this);

  return NOX::Epetra::Group::isJacobian();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Group::throwError(const std::string& functionName, const std::string& errorMsg) const
{
  std::ostringstream msg;
  msg << "ERROR - NOX::NLN::Group::" << functionName << " - " << errorMsg << std::endl;

  dserror(msg.str());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_BlockMap& NOX::NLN::Group::getDofMap() const
{
  return xVector.getEpetraVector().Map();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Group::CreateBackupState(const NOX::Abstract::Vector& dir) const
{
  Teuchos::RCP<NOX::NLN::Interface::Required> nln_required =
      Teuchos::rcp_dynamic_cast<NOX::NLN::Interface::Required>(userInterfacePtr, true);

  const NOX::Epetra::Vector* epetra_dir = dynamic_cast<const NOX::Epetra::Vector*>(&dir);
  if (not epetra_dir) dserror("Dynamic cast failed.");

  nln_required->CreateBackupState(epetra_dir->getEpetraVector());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Group::RecoverFromBackupState()
{
  Teuchos::RCP<NOX::NLN::Interface::Required> nln_required =
      Teuchos::rcp_dynamic_cast<NOX::NLN::Interface::Required>(userInterfacePtr);

  if (nln_required.is_null()) dserror("Dynamic cast failed.");

  nln_required->RecoverFromBackupState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map& NOX::NLN::Group::getJacobianRangeMap(unsigned rbid, unsigned cbid) const
{
  if (not isJacobian())
    dserror(
        "It is not possible to access the Jacobian since it has not yet "
        "been evaluated.");

  Teuchos::RCP<NOX::NLN::LinearSystem> nlnSharedLinearSystem =
      Teuchos::rcp_dynamic_cast<NOX::NLN::LinearSystem>(sharedLinearSystem.getObject(this), true);

  return nlnSharedLinearSystem->getJacobianRangeMap(rbid, cbid);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> NOX::NLN::Group::getDiagonalOfJacobian(unsigned diag_bid) const
{
  if (not isJacobian())
    dserror(
        "It is not possible to access the Jacobian since it has not yet "
        "been evaluated.");

  Teuchos::RCP<NOX::NLN::LinearSystem> nlnSharedLinearSystem =
      Teuchos::rcp_dynamic_cast<NOX::NLN::LinearSystem>(sharedLinearSystem.getObject(this), true);

  return nlnSharedLinearSystem->getDiagonalOfJacobian(diag_bid);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Group::replaceDiagonalOfJacobian(
    const Epetra_Vector& new_diag, unsigned diag_bid) const
{
  if (not isJacobian())
    dserror(
        "It is not possible to access the Jacobian since it has not yet "
        "been evaluated.");

  Teuchos::RCP<NOX::NLN::LinearSystem> nlnSharedLinearSystem =
      Teuchos::rcp_dynamic_cast<NOX::NLN::LinearSystem>(sharedLinearSystem.getObject(this), true);

  return nlnSharedLinearSystem->replaceDiagonalOfJacobian(new_diag, diag_bid);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::NLN::Group::computeJacobianConditionNumber(
    int maxIters, double tolerance, int krylovSubspaceSize, bool printOutput)
{
  if (isConditionNumber()) return NOX::Abstract::Group::Ok;

  if (maxIters <= 0)
    dserror(
        "The direct computation of the condition number via LAPACK is "
        "in parallel not possible. Please try the GMRES based variant by "
        "providing a meaningful set of input parameters.");

  if (azconditionnumberptr_.is_null())
    azconditionnumberptr_ = Teuchos::rcp(new LINALG::AztecOOConditionNumber);
  azConditionNumberPtr = Teuchos::rcpFromRef(*azconditionnumberptr_);

  NOX::Abstract::Group::ReturnType rtype = NOX::Epetra::Group::computeJacobianConditionNumber(
      maxIters, tolerance, krylovSubspaceSize, printOutput);

  ev_.setAztecEstimates(*azconditionnumberptr_);

  return rtype;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::NLN::Group::computeSerialJacobianConditionNumber(
    const NOX::NLN::LinSystem::ConditionNumber condnum_type, bool printOutput)
{
  if (isConditionNumber()) return NOX::Abstract::Group::Ok;

  if (xVector.getEpetraVector().Comm().NumProc() != 1) dserror("Only serial mode is supported!");

  if (!isJacobian()) dserror("Jacobian is invalid wrt the solution.");

  switch (condnum_type)
  {
    case LinSystem::ConditionNumber::max_min_ev_ratio:
    {
      computeSerialJacobianEigenvalues(printOutput);

      if (ev_.real_min_ == 0)
      {
        utils.out(NOX::Utils::Warning) << "Jacobian is singular!" << std::endl;
        conditionNumber = std::numeric_limits<double>::infinity();
      }
      else
        conditionNumber = ev_.real_max_ / ev_.real_min_;

      break;
    }
    case LinSystem::ConditionNumber::one_norm:
    case LinSystem::ConditionNumber::inf_norm:
    {
      Teuchos::RCP<NOX::NLN::LinearSystem> nlnSharedLinearSystem =
          Teuchos::rcp_dynamic_cast<NOX::NLN::LinearSystem>(
              sharedLinearSystem.getObject(this), true);

      conditionNumber = nlnSharedLinearSystem->computeSerialConditionNumberOfJacobian(condnum_type);

      break;
    }
    default:
      dserror("Unknown LinSystem::ConditionNumber type!");
      exit(EXIT_FAILURE);
  }

  if (printOutput)
    utils.out() << "Condition number = " << NOX::Utils::sciformat(conditionNumber, 5) << std::endl;

  isValidConditionNumber = true;

  return NOX::Abstract::Group::Ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::NLN::Group::computeSerialJacobianEigenvalues(bool printOutput)
{
  if (xVector.getEpetraVector().Comm().NumProc() != 1) dserror("Works only in serial mode!");

  if (ev_.isvalid_) return NOX::Abstract::Group::Ok;

  Teuchos::RCP<NOX::NLN::LinearSystem> nlnSharedLinearSystem =
      Teuchos::rcp_dynamic_cast<NOX::NLN::LinearSystem>(sharedLinearSystem.getObject(this), true);

  nlnSharedLinearSystem->computeSerialEigenvaluesOfJacobian(ev_.realpart_, ev_.imaginarypart_);

  const int length = ev_.realpart_.Length();
  ev_.real_max_ = *std::max_element(ev_.realpart_.A(), ev_.realpart_.A() + length);
  ev_.real_min_ = *std::min_element(ev_.realpart_.A(), ev_.realpart_.A() + length);

  if (printOutput)
  {
    utils.out() << "maximal eigenvalue = " << utils.sciformat(ev_.real_max_, 5) << "\n";
    utils.out() << "minimal eigenvalue = " << utils.sciformat(ev_.real_min_, 5) << "\n";
  }

  ev_.isvalid_ = true;

  return NOX::Abstract::Group::Ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Group::isEigenvalues() const { return ev_.isvalid_; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const LINALG::SerialDenseVector& NOX::NLN::Group::getJacoianRealEigenvalues() const
{
  if (not isEigenvalues()) dserror("The eigenvalues has not yet been computed!");
  return ev_.realpart_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const LINALG::SerialDenseVector& NOX::NLN::Group::getJacoianImaginaryEigenvalues() const
{
  if (not isEigenvalues()) dserror("The eigenvalues has not yet been computed!");
  return ev_.imaginarypart_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::Group::getJacobianMaxRealEigenvalue() const
{
  if (not isEigenvalues()) dserror("The eigenvalues has not yet been computed!");
  return ev_.real_max_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::Group::getJacobianMinRealEigenvalue() const
{
  if (not isEigenvalues()) dserror("The eigenvalues has not yet been computed!");
  return ev_.real_min_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Group::Eigenvalues& NOX::NLN::Group::Eigenvalues::operator=(const Eigenvalues& src)
{
  if (not src.isvalid_)
  {
    this->isvalid_ = false;
    return *this;
  }

  this->realpart_.Resize(src.realpart_.Length());
  this->realpart_ = src.realpart_;

  this->imaginarypart_.Resize(src.imaginarypart_.Length());
  this->imaginarypart_ = src.imaginarypart_;

  this->real_max_ = src.real_max_;
  this->real_min_ = src.real_min_;

  this->isvalid_ = src.isvalid_;

  return *this;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Group::Eigenvalues::setAztecEstimates(
    LINALG::AztecOOConditionNumber& azconditionnumber)
{
  if (isvalid_) return;

  real_max_ = azconditionnumber.getStatus(AZ_lambda_real_max);
  real_min_ = azconditionnumber.getStatus(AZ_lambda_real_min);

  realpart_.Size(2);
  realpart_(0) = azconditionnumber.getStatus(AZ_lambda_real_max);
  realpart_(1) = azconditionnumber.getStatus(AZ_lambda_real_min);

  imaginarypart_.Size(2);
  imaginarypart_(0) = azconditionnumber.getStatus(AZ_lambda_imag_max);
  imaginarypart_(1) = azconditionnumber.getStatus(AZ_lambda_imag_min);

  isvalid_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Group::getDofsFromElements(
    const std::vector<int>& my_ele_gids, std::set<int>& my_ele_dofs) const
{
  GetNlnReqInterfacePtr()->getDofsFromElements(my_ele_gids, my_ele_dofs);
}

/*-----------------------------------------------------------*/
/*! \file

\brief Implementation of LOCA::MultiContinuation::ConstraintInterfaceMVDX for arc-length
continuation.


\level 3

*/
/*-----------------------------------------------------------*/

#include "loca_str_constraint_arclength.H"
#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LOCA::STR::MultiContinuation::ArcLengthConstraint::ArcLengthConstraint()
    : isinit_(false),
      issetup_(false),
      isvalid_constraints_(false),
      isvalid_dx_(false),
      numconstraints_(0),
      constraints_(0, 0),
      loca_param_vec_(),
      gstate_ptr_(Teuchos::null),
      xvector_ptr_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LOCA::STR::MultiContinuation::ArcLengthConstraint::ArcLengthConstraint(
    const ArcLengthConstraint& source, NOX::CopyType type)
    : isinit_(source.isinit_),
      issetup_(source.issetup_),
      numconstraints_(source.numconstraints_),
      constraints_(source.constraints_),
      loca_param_vec_(source.loca_param_vec_),
      gstate_ptr_(source.gstate_ptr_),
      xvector_ptr_(source.xvector_ptr_->clone(type))
{
  if (type == NOX::DeepCopy)
  {
    isvalid_constraints_ = source.isvalid_constraints_;
    isvalid_dx_ = source.isvalid_dx_;
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LOCA::STR::MultiContinuation::ArcLengthConstraint::Init(
    const Teuchos::RCP<const ::STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
    const Teuchos::RCP<LOCA::ParameterVector>& loca_param_vec_ptr,
    const Teuchos::RCP<NOX::Abstract::Vector>& xvector_ptr)
{
  issetup_ = false;

  loca_param_vec_ = *loca_param_vec_ptr;
  gstate_ptr_ = gstate_ptr;
  xvector_ptr_ = xvector_ptr->createMultiVector(1);

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LOCA::STR::MultiContinuation::ArcLengthConstraint::Setup()
{
  CheckInit();

  numconstraints_ = 1;
  constraints_ = NOX::Abstract::MultiVector::DenseMatrix(1, 1);

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LOCA::STR::MultiContinuation::ArcLengthConstraint::CheckInitSetup() const
{
  if (!IsInit() or !IsSetup()) dserror("Call first Init() and Setup()!");

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LOCA::STR::MultiContinuation::ArcLengthConstraint::CheckInit() const
{
  if (!IsInit()) dserror("Call first Init()!");

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LOCA::STR::MultiContinuation::ArcLengthConstraint::copy(const ConstraintInterface& src)
{
  const ArcLengthConstraint& source = dynamic_cast<const ArcLengthConstraint&>(src);

  if (this == &source) return;

  isinit_ = source.isinit_;
  issetup_ = source.issetup_;
  isvalid_constraints_ = source.isvalid_constraints_;
  isvalid_dx_ = source.isvalid_dx_;
  numconstraints_ = source.numconstraints_;
  constraints_ = source.constraints_;
  loca_param_vec_ = source.loca_param_vec_;
  gstate_ptr_ = source.gstate_ptr_;
  xvector_ptr_ = source.xvector_ptr_;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
LOCA::STR::MultiContinuation::ArcLengthConstraint::clone(NOX::CopyType type) const
{
  return ArcLengthConstraint(*this, type);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LOCA::STR::MultiContinuation::ArcLengthConstraint::resetIsValid()
{
  isvalid_constraints_ = false;
  isvalid_dx_ = false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int LOCA::STR::MultiContinuation::ArcLengthConstraint::numConstraints() const
{
  return numconstraints_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LOCA::STR::MultiContinuation::ArcLengthConstraint::setX(const NOX::Abstract::Vector& y)
{
  (*xvector_ptr_)[0] = y;
  resetIsValid();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LOCA::STR::MultiContinuation::ArcLengthConstraint::setParam(int paramID, double val)
{
  loca_param_vec_[paramID] = val;
  resetIsValid();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LOCA::STR::MultiContinuation::ArcLengthConstraint::setParams(
    const std::vector<int>& paramIDs, const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  for (std::size_t i = 0; i < paramIDs.size(); ++i) loca_param_vec_[paramIDs[i]] = vals(i, 0);
  resetIsValid();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType
LOCA::STR::MultiContinuation::ArcLengthConstraint::computeConstraints()
{
  return NOX::Abstract::Group::Ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool LOCA::STR::MultiContinuation::ArcLengthConstraint::isConstraints() const
{
  return isvalid_constraints_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool LOCA::STR::MultiContinuation::ArcLengthConstraint::isDX() const { return isvalid_dx_; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LOCA::STR::MultiContinuation::ArcLengthConstraint::preProcessContinuationStep(
    LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  /* currently empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LOCA::STR::MultiContinuation::ArcLengthConstraint::postProcessContinuationStep(
    LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  /* currently empty */
}

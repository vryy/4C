/*----------------------------------------------------------------------*/
/*! \file

\brief A family of abstract nice algebra operations (ANA)

\level 1

*----------------------------------------------------------------------*/
#include "4C_linear_solver_preconditioner_linalg_ana.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
  implicit operator product Operator*Operator
 *----------------------------------------------------------------------*/
int Core::LinAlg::Ana::OperatorProduct::Apply(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (!usetransposed_)  // A*B*x
  {
    Core::LinAlg::Ana::Vector tmp(right_->OperatorRangeMap(), false);
    int err1 = right_->Apply(X, tmp);
    if (err1) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err1);
    int err2 = left_->Apply(tmp, Y);
    if (err2) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err2);
    return err1 - err2;
  }
  else  // (A*B)^T * x = B^T * A^T * x
  {
    Core::LinAlg::Ana::Vector tmp(left_->OperatorDomainMap(), false);
    const_cast<LightWeightOperatorBase&>(*left_).SetUseTranspose(!left_->UseTranspose());
    int err1 = left_->Apply(X, tmp);
    if (err1) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err1);
    const_cast<LightWeightOperatorBase&>(*left_).SetUseTranspose(!left_->UseTranspose());
    const_cast<LightWeightOperatorBase&>(*right_).SetUseTranspose(!right_->UseTranspose());
    int err2 = right_->Apply(tmp, Y);
    if (err2) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err2);
    const_cast<LightWeightOperatorBase&>(*right_).SetUseTranspose(!right_->UseTranspose());
    return err1 - err2;
  }
  return 0;
}
int Core::LinAlg::Ana::OperatorProduct::ApplyInverse(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (!usetransposed_)  // (A*B)^{-1} * x = B^{-1} * A^{-1} * x
  {
    Core::LinAlg::Ana::Vector tmp(left_->OperatorDomainMap(), false);
    int err1 = left_->ApplyInverse(X, tmp);
    if (err1) FOUR_C_THROW("LightWeightOperatorBase::ApplyInverse returned err=%d", err1);
    int err2 = right_->ApplyInverse(tmp, Y);
    if (err2) FOUR_C_THROW("LightWeightOperatorBase::ApplyInverse returned err=%d", err2);
    return err1 - err2;
  }
  else  // [ (A*B)^T ]^{-1} * x = A^{T,-1} * B^{T,-1} * x
  {
    Core::LinAlg::Ana::Vector tmp(right_->OperatorRangeMap(), false);
    const_cast<LightWeightOperatorBase&>(*right_).SetUseTranspose(!right_->UseTranspose());
    int err1 = right_->ApplyInverse(X, tmp);
    if (err1) FOUR_C_THROW("LightWeightOperatorBase::ApplyInverse returned err=%d", err1);
    const_cast<LightWeightOperatorBase&>(*right_).SetUseTranspose(!right_->UseTranspose());
    const_cast<LightWeightOperatorBase&>(*left_).SetUseTranspose(!left_->UseTranspose());
    int err2 = left_->ApplyInverse(tmp, Y);
    if (err2) FOUR_C_THROW("LightWeightOperatorBase::ApplyInverse returned err=%d", err2);
    const_cast<LightWeightOperatorBase&>(*left_).SetUseTranspose(!left_->UseTranspose());
    return err1 - err2;
  }
  return -1;
}



/*----------------------------------------------------------------------*
  implicit operator sum (Operator+Operator)*x
 *----------------------------------------------------------------------*/
int Core::LinAlg::Ana::OperatorSum::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (!usetransposed_)  // (A+-B)*x = Ax +- Bx
  {
    Core::LinAlg::Ana::Vector tmp(right_->OperatorRangeMap(), false);
    int err1 = right_->Apply(X, tmp);
    if (err1) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err1);
    int err2 = left_->Apply(X, Y);
    if (err2) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err2);
    Y.Update((double)sign_, tmp, 1.0);
    return err1 - err2;
  }
  else  // (A+-B)^T * x = A^T x +- B^T * x
  {
    Core::LinAlg::Ana::Vector tmp(right_->OperatorDomainMap(), false);
    const_cast<LightWeightOperatorBase&>(*right_).SetUseTranspose(!right_->UseTranspose());
    int err1 = right_->Apply(X, tmp);
    if (err1) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err1);
    const_cast<LightWeightOperatorBase&>(*right_).SetUseTranspose(!right_->UseTranspose());
    const_cast<LightWeightOperatorBase&>(*left_).SetUseTranspose(!left_->UseTranspose());
    int err2 = left_->Apply(X, Y);
    if (err2) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err2);
    const_cast<LightWeightOperatorBase&>(*left_).SetUseTranspose(!left_->UseTranspose());
    Y.Update((double)sign_, tmp, 1.0);
    return err1 - err2;
  }
  return 0;
}


/*----------------------------------------------------------------------*
  apply inverse of an operator
 *----------------------------------------------------------------------*/
int Core::LinAlg::Ana::OperatorInverse::Apply(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // wrap column 0 of in and output vectors to Epetra_Vector
  const Epetra_Vector* invec = X(0);
  Teuchos::RCP<Epetra_Vector> in = Teuchos::rcp(const_cast<Epetra_Vector*>(invec), false);
  Epetra_Vector* outvec = Y(0);
  Teuchos::RCP<Epetra_Vector> out = Teuchos::rcp(outvec, false);

  // wrap underlying operator
  Teuchos::RCP<Epetra_Operator> rcpop = Teuchos::rcp(const_cast<Epetra_Operator*>(&op_), false);

  out->PutScalar(0.0);
  Core::LinAlg::SolverParams solver_params;
  solver_params.refactor = true;
  solver_params.reset = reset_;
  solver_.Solve(rcpop, out, in, solver_params);

  return 0;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::Ana::LcVecPointwiseLc::Update(
    Core::LinAlg::Ana::Vector& v, const double& scale) const
{
#if DEBUGGING_ANA
  cout << "Core::LinAlg::Ana::LC_vec_pointwise_lc::Update(Core::LinAlg::Ana::Vector& v, const "
          "double& scale)"
       << endl;
  fflush(stdout);
#endif
  // we can not avoid extra memory in this case
  Core::LinAlg::Ana::Vector tmp(right_.RangeMap(), false);
  right_.set(tmp, 1.0);
  v.Multiply(scale, vec_, tmp, 1.0);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::Ana::LcLcPointwiseLc::set(
    Core::LinAlg::Ana::Vector& v, const double& scale) const
{
#if DEBUGGING_ANA
  cout << "Core::LinAlg::Ana::LC_lc_pointwise_lc::set(Core::LinAlg::Ana::Vector& v, const double& "
          "scale)"
       << endl;
  fflush(stdout);
#endif
  left_.set(v, 1.0);
  Core::LinAlg::Ana::Vector tmp(right_.RangeMap(), false);
  right_.set(tmp, 1.0);
  v.Multiply(scale, v, tmp, 0.0);
}
void Core::LinAlg::Ana::LcLcPointwiseLc::Update(
    Core::LinAlg::Ana::Vector& v, const double& scale) const
{
#if DEBUGGING_ANA
  cout << "Core::LinAlg::Ana::LC_lc_pointwise_lc::Update(Core::LinAlg::Ana::Vector& v, const "
          "double& scale)"
       << endl;
  fflush(stdout);
#endif
  Core::LinAlg::Ana::Vector tmp1(left_.RangeMap(), false);
  left_.set(tmp1, 1.0);
  Core::LinAlg::Ana::Vector tmp2(right_.RangeMap(), false);
  right_.set(tmp2, 1.0);
  v.Multiply(scale, tmp1, tmp2, 1.0);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::Ana::LcOperatorTimesLcsv::set(
    Core::LinAlg::Ana::Vector& v, const double& scale) const
{
#if DEBUGGING_ANA
  cout << "Core::LinAlg::Ana::LC_Operator_times_lcsv::set(Core::LinAlg::Ana::Vector& v, const "
          "double& scale)"
       << endl;
  fflush(stdout);
#endif
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (!v.Map().SameAs(op_->OperatorRangeMap())) FOUR_C_THROW("Range maps don't match - fatal");
  if (!right_.Vector().Map().SameAs(op_->OperatorDomainMap()))
    FOUR_C_THROW("Domain maps don't match - fatal");
#endif
  int err = op_->Apply(right_.Vector(), v);
  if (err) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err);
  if (scale * right_.Scalar() != 1.0) v.Scale(scale * right_.Scalar());
}
void Core::LinAlg::Ana::LcOperatorTimesLcsv::Update(
    Core::LinAlg::Ana::Vector& v, const double& scale) const
{
#if DEBUGGING_ANA
  cout << "Core::LinAlg::Ana::LC_Operator_times_lcsv::Update(Core::LinAlg::Ana::Vector& v, const "
          "double& scale)"
       << endl;
  fflush(stdout);
#endif
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (!v.Map().SameAs(op_->OperatorRangeMap())) FOUR_C_THROW("Range maps don't match - fatal");
  if (!right_.Vector().Map().SameAs(op_->OperatorDomainMap()))
    FOUR_C_THROW("Domain maps don't match - fatal");
#endif
  Core::LinAlg::Ana::Vector tmp(op_->OperatorRangeMap(), false);
  int err = op_->Apply(right_.Vector(), tmp);
  if (err) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err);
  v.Update(scale * right_.Scalar(), tmp, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::Ana::LcOperatorTimesLc::set(
    Core::LinAlg::Ana::Vector& v, const double& scale) const
{
#if DEBUGGING_ANA
  cout << "Core::LinAlg::Ana::LC_Operator_times_lc::set(Core::LinAlg::Ana::Vector& v, const "
          "double& scale)"
       << endl;
  fflush(stdout);
#endif
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (!v.Map().SameAs(op_->OperatorRangeMap())) FOUR_C_THROW("Range maps don't match - fatal");
  if (!right_.RangeMap().SameAs(op_->OperatorDomainMap()))
    FOUR_C_THROW("Domain maps don't match - fatal");
#endif
  Core::LinAlg::Ana::Vector tmp(right_.RangeMap(), false);
  right_.set(tmp, scale);
  int err = op_->Apply(tmp, v);
  if (err) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err);
}
void Core::LinAlg::Ana::LcOperatorTimesLc::Update(
    Core::LinAlg::Ana::Vector& v, const double& scale) const
{
#if DEBUGGING_ANA
  cout << "Core::LinAlg::Ana::LC_Operator_times_lc::Update(Core::LinAlg::Ana::Vector& v, const "
          "double& scale)"
       << endl;
  fflush(stdout);
#endif
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (!v.Map().SameAs(op_->OperatorRangeMap())) FOUR_C_THROW("Range maps don't match - fatal");
  if (!right_.RangeMap().SameAs(op_->OperatorDomainMap()))
    FOUR_C_THROW("Domain maps don't match - fatal");
#endif
  Core::LinAlg::Ana::Vector tmp1(right_.RangeMap(), false);
  right_.set(tmp1, 1.0);
  Core::LinAlg::Ana::Vector tmp2(op_->OperatorRangeMap(), false);
  int err = op_->Apply(tmp1, tmp2);
  if (err) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err);
  v.Update(scale, tmp2, 1.0);
}



/*----------------------------------------------------------------------*
   vec dot vec (result is scalar)
 *----------------------------------------------------------------------*/
double Core::LinAlg::Ana::operator*(
    const Core::LinAlg::Ana::Vector& vec1, const Core::LinAlg::Ana::Vector& vec2)
{
#if DEBUGGING_ANA
  cout << "double operator* (const Core::LinAlg::Ana::Vector& vec1, const "
          "Core::LinAlg::Ana::Vector& vec2)"
       << endl;
  fflush(stdout);
#endif
  double result;
  vec1.Dot(vec2, &result);
  return result;
}

/*----------------------------------------------------------------------*
   lc dot lc (result is scalar)
 *----------------------------------------------------------------------*/
double Core::LinAlg::Ana::operator*(
    const Core::LinAlg::Ana::LCBase& left, const Core::LinAlg::Ana::LCBase& right)
{
#if DEBUGGING_ANA
  cout << "double operator* (const Core::LinAlg::Ana::LCBase& left, const "
          "Core::LinAlg::Ana::LCBase& right)"
       << endl;
  fflush(stdout);
#endif
  Core::LinAlg::Ana::Vector tmpleft(left.RangeMap(), false);
  left.set(tmpleft, 1.0);
  Core::LinAlg::Ana::Vector tmpright(right.RangeMap(), false);
  right.set(tmpright, 1.0);
  double result;
  tmpleft.Dot(tmpright, &result);
  return result;
}

/*----------------------------------------------------------------------*
   vec dot lc (result is scalar)
 *----------------------------------------------------------------------*/
double Core::LinAlg::Ana::operator*(
    const Core::LinAlg::Ana::Vector& vec1, const Core::LinAlg::Ana::LCBase& right)
{
#if DEBUGGING_ANA
  cout << "double operator* (const Core::LinAlg::Ana::Vector& vec1, const "
          "Core::LinAlg::Ana::LCBase& right)"
       << endl;
  fflush(stdout);
#endif
  Core::LinAlg::Ana::Vector tmp(right.RangeMap(), false);
  right.set(tmp, 1.0);
  double result;
  vec1.Dot(tmp, &result);
  return result;
}

/*----------------------------------------------------------------------*
   vec dot lcsv (result is scalar) (specialization)
 *----------------------------------------------------------------------*/
double Core::LinAlg::Ana::operator*(
    const Core::LinAlg::Ana::Vector& vec1, const Core::LinAlg::Ana::LCSTimesVec& right)
{
#if DEBUGGING_ANA
  cout << "double operator* (const Core::LinAlg::Ana::Vector& vec1, const "
          "Core::LinAlg::Ana::LCSTimesVec& "
          "right)"
       << endl;
  fflush(stdout);
#endif
  double result;
  vec1.Dot(right.Vector(), &result);
  result *= right.Scalar();
  return result;
}

/*----------------------------------------------------------------------*
   lcsv dot lcsv (result is scalar) (specialization)
 *----------------------------------------------------------------------*/
double Core::LinAlg::Ana::operator*(
    const Core::LinAlg::Ana::LCSTimesVec& left, const Core::LinAlg::Ana::LCSTimesVec& right)
{
#if DEBUGGING_ANA
  cout << "double operator* (const Core::LinAlg::Ana::LCSTimesVec& left, const "
          "Core::LinAlg::Ana::LCSTimesVec& right)"
       << endl;
  fflush(stdout);
#endif
  double result;
  left.Vector().Dot(right.Vector(), &result);
  return (result * left.Scalar() * right.Scalar());
}

/*----------------------------------------------------------------------*
   norms
 *----------------------------------------------------------------------*/
double Core::LinAlg::Ana::norm2(const Core::LinAlg::Ana::LCBase& lc)
{
  Core::LinAlg::Ana::Vector tmp(lc.RangeMap(), false);
  lc.set(tmp, 1.0);
  return norm2(tmp);
}
double Core::LinAlg::Ana::norm1(const Core::LinAlg::Ana::LCBase& lc)
{
  Core::LinAlg::Ana::Vector tmp(lc.RangeMap(), false);
  lc.set(tmp, 1.0);
  return norm1(tmp);
}
double Core::LinAlg::Ana::norminf(const Core::LinAlg::Ana::LCBase& lc)
{
  Core::LinAlg::Ana::Vector tmp(lc.RangeMap(), false);
  lc.set(tmp, 1.0);
  return norminf(tmp);
}

FOUR_C_NAMESPACE_CLOSE

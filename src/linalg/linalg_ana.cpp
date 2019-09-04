/*----------------------------------------------------------------------*/
/*! \file

\brief A family of abstract nice algebra operations (ANA)

\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*----------------------------------------------------------------------*/
#include "linalg_ana.H"

/*----------------------------------------------------------------------*
  implicit operator product Operator*Operator
 *----------------------------------------------------------------------*/
int LINALG::ANA::OperatorProduct::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (!usetransposed_)  // A*B*x
  {
    LINALG::ANA::Vector tmp(right_->OperatorRangeMap(), false);
    int err1 = right_->Apply(X, tmp);
    if (err1) dserror("LightWeightOperatorBase::Apply returned err=%d", err1);
    int err2 = left_->Apply(tmp, Y);
    if (err2) dserror("LightWeightOperatorBase::Apply returned err=%d", err2);
    return err1 - err2;
  }
  else  // (A*B)^T * x = B^T * A^T * x
  {
    LINALG::ANA::Vector tmp(left_->OperatorDomainMap(), false);
    const_cast<LightWeightOperatorBase&>(*left_).SetUseTranspose(!left_->UseTranspose());
    int err1 = left_->Apply(X, tmp);
    if (err1) dserror("LightWeightOperatorBase::Apply returned err=%d", err1);
    const_cast<LightWeightOperatorBase&>(*left_).SetUseTranspose(!left_->UseTranspose());
    const_cast<LightWeightOperatorBase&>(*right_).SetUseTranspose(!right_->UseTranspose());
    int err2 = right_->Apply(tmp, Y);
    if (err2) dserror("LightWeightOperatorBase::Apply returned err=%d", err2);
    const_cast<LightWeightOperatorBase&>(*right_).SetUseTranspose(!right_->UseTranspose());
    return err1 - err2;
  }
  return 0;
}
int LINALG::ANA::OperatorProduct::ApplyInverse(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (!usetransposed_)  // (A*B)^{-1} * x = B^{-1} * A^{-1} * x
  {
    LINALG::ANA::Vector tmp(left_->OperatorDomainMap(), false);
    int err1 = left_->ApplyInverse(X, tmp);
    if (err1) dserror("LightWeightOperatorBase::ApplyInverse returned err=%d", err1);
    int err2 = right_->ApplyInverse(tmp, Y);
    if (err2) dserror("LightWeightOperatorBase::ApplyInverse returned err=%d", err2);
    return err1 - err2;
  }
  else  // [ (A*B)^T ]^{-1} * x = A^{T,-1} * B^{T,-1} * x
  {
    LINALG::ANA::Vector tmp(right_->OperatorRangeMap(), false);
    const_cast<LightWeightOperatorBase&>(*right_).SetUseTranspose(!right_->UseTranspose());
    int err1 = right_->ApplyInverse(X, tmp);
    if (err1) dserror("LightWeightOperatorBase::ApplyInverse returned err=%d", err1);
    const_cast<LightWeightOperatorBase&>(*right_).SetUseTranspose(!right_->UseTranspose());
    const_cast<LightWeightOperatorBase&>(*left_).SetUseTranspose(!left_->UseTranspose());
    int err2 = left_->ApplyInverse(tmp, Y);
    if (err2) dserror("LightWeightOperatorBase::ApplyInverse returned err=%d", err2);
    const_cast<LightWeightOperatorBase&>(*left_).SetUseTranspose(!left_->UseTranspose());
    return err1 - err2;
  }
  return -1;
}



/*----------------------------------------------------------------------*
  implicit operator sum (Operator+Operator)*x
 *----------------------------------------------------------------------*/
int LINALG::ANA::OperatorSum::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (!usetransposed_)  // (A+-B)*x = Ax +- Bx
  {
    LINALG::ANA::Vector tmp(right_->OperatorRangeMap(), false);
    int err1 = right_->Apply(X, tmp);
    if (err1) dserror("LightWeightOperatorBase::Apply returned err=%d", err1);
    int err2 = left_->Apply(X, Y);
    if (err2) dserror("LightWeightOperatorBase::Apply returned err=%d", err2);
    Y.Update((double)sign_, tmp, 1.0);
    return err1 - err2;
  }
  else  // (A+-B)^T * x = A^T x +- B^T * x
  {
    LINALG::ANA::Vector tmp(right_->OperatorDomainMap(), false);
    const_cast<LightWeightOperatorBase&>(*right_).SetUseTranspose(!right_->UseTranspose());
    int err1 = right_->Apply(X, tmp);
    if (err1) dserror("LightWeightOperatorBase::Apply returned err=%d", err1);
    const_cast<LightWeightOperatorBase&>(*right_).SetUseTranspose(!right_->UseTranspose());
    const_cast<LightWeightOperatorBase&>(*left_).SetUseTranspose(!left_->UseTranspose());
    int err2 = left_->Apply(X, Y);
    if (err2) dserror("LightWeightOperatorBase::Apply returned err=%d", err2);
    const_cast<LightWeightOperatorBase&>(*left_).SetUseTranspose(!left_->UseTranspose());
    Y.Update((double)sign_, tmp, 1.0);
    return err1 - err2;
  }
  return 0;
}


/*----------------------------------------------------------------------*
  apply inverse of an operator
 *----------------------------------------------------------------------*/
int LINALG::ANA::OperatorInverse::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // wrap column 0 of in and output vectors to Epetra_Vector
  const Epetra_Vector* invec = X(0);
  Teuchos::RCP<Epetra_Vector> in = Teuchos::rcp(const_cast<Epetra_Vector*>(invec), false);
  Epetra_Vector* outvec = Y(0);
  Teuchos::RCP<Epetra_Vector> out = Teuchos::rcp(outvec, false);

  // wrap underlying operator
  Teuchos::RCP<Epetra_Operator> rcpop = Teuchos::rcp(const_cast<Epetra_Operator*>(&op_), false);

  out->PutScalar(0.0);
  solver_.Solve(rcpop, out, in, true, reset_);

  return 0;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::ANA::LC_vec_pointwise_lc::Update(LINALG::ANA::Vector& v, const double& scale) const
{
#if DEBUGGING_ANA
  cout << "LINALG::ANA::LC_vec_pointwise_lc::Update(LINALG::ANA::Vector& v, const double& scale)"
       << endl;
  fflush(stdout);
#endif
  // we can not avoid extra memory in this case
  LINALG::ANA::Vector tmp(right_.RangeMap(), false);
  right_.Set(tmp, 1.0);
  v.Multiply(scale, vec_, tmp, 1.0);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::ANA::LC_lc_pointwise_lc::Set(LINALG::ANA::Vector& v, const double& scale) const
{
#if DEBUGGING_ANA
  cout << "LINALG::ANA::LC_lc_pointwise_lc::Set(LINALG::ANA::Vector& v, const double& scale)"
       << endl;
  fflush(stdout);
#endif
  left_.Set(v, 1.0);
  LINALG::ANA::Vector tmp(right_.RangeMap(), false);
  right_.Set(tmp, 1.0);
  v.Multiply(scale, v, tmp, 0.0);
}
void LINALG::ANA::LC_lc_pointwise_lc::Update(LINALG::ANA::Vector& v, const double& scale) const
{
#if DEBUGGING_ANA
  cout << "LINALG::ANA::LC_lc_pointwise_lc::Update(LINALG::ANA::Vector& v, const double& scale)"
       << endl;
  fflush(stdout);
#endif
  LINALG::ANA::Vector tmp1(left_.RangeMap(), false);
  left_.Set(tmp1, 1.0);
  LINALG::ANA::Vector tmp2(right_.RangeMap(), false);
  right_.Set(tmp2, 1.0);
  v.Multiply(scale, tmp1, tmp2, 1.0);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::ANA::LC_Operator_times_lcsv::Set(LINALG::ANA::Vector& v, const double& scale) const
{
#if DEBUGGING_ANA
  cout << "LINALG::ANA::LC_Operator_times_lcsv::Set(LINALG::ANA::Vector& v, const double& scale)"
       << endl;
  fflush(stdout);
#endif
#ifdef DEBUG
  if (!v.Map().SameAs(op_->OperatorRangeMap())) dserror("Range maps don't match - fatal");
  if (!right_.Vector().Map().SameAs(op_->OperatorDomainMap()))
    dserror("Domain maps don't match - fatal");
#endif
  int err = op_->Apply(right_.Vector(), v);
  if (err) dserror("LightWeightOperatorBase::Apply returned err=%d", err);
  if (scale * right_.Scalar() != 1.0) v.Scale(scale * right_.Scalar());
}
void LINALG::ANA::LC_Operator_times_lcsv::Update(LINALG::ANA::Vector& v, const double& scale) const
{
#if DEBUGGING_ANA
  cout << "LINALG::ANA::LC_Operator_times_lcsv::Update(LINALG::ANA::Vector& v, const double& scale)"
       << endl;
  fflush(stdout);
#endif
#ifdef DEBUG
  if (!v.Map().SameAs(op_->OperatorRangeMap())) dserror("Range maps don't match - fatal");
  if (!right_.Vector().Map().SameAs(op_->OperatorDomainMap()))
    dserror("Domain maps don't match - fatal");
#endif
  LINALG::ANA::Vector tmp(op_->OperatorRangeMap(), false);
  int err = op_->Apply(right_.Vector(), tmp);
  if (err) dserror("LightWeightOperatorBase::Apply returned err=%d", err);
  v.Update(scale * right_.Scalar(), tmp, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::ANA::LC_Operator_times_lc::Set(LINALG::ANA::Vector& v, const double& scale) const
{
#if DEBUGGING_ANA
  cout << "LINALG::ANA::LC_Operator_times_lc::Set(LINALG::ANA::Vector& v, const double& scale)"
       << endl;
  fflush(stdout);
#endif
#ifdef DEBUG
  if (!v.Map().SameAs(op_->OperatorRangeMap())) dserror("Range maps don't match - fatal");
  if (!right_.RangeMap().SameAs(op_->OperatorDomainMap()))
    dserror("Domain maps don't match - fatal");
#endif
  LINALG::ANA::Vector tmp(right_.RangeMap(), false);
  right_.Set(tmp, scale);
  int err = op_->Apply(tmp, v);
  if (err) dserror("LightWeightOperatorBase::Apply returned err=%d", err);
}
void LINALG::ANA::LC_Operator_times_lc::Update(LINALG::ANA::Vector& v, const double& scale) const
{
#if DEBUGGING_ANA
  cout << "LINALG::ANA::LC_Operator_times_lc::Update(LINALG::ANA::Vector& v, const double& scale)"
       << endl;
  fflush(stdout);
#endif
#ifdef DEBUG
  if (!v.Map().SameAs(op_->OperatorRangeMap())) dserror("Range maps don't match - fatal");
  if (!right_.RangeMap().SameAs(op_->OperatorDomainMap()))
    dserror("Domain maps don't match - fatal");
#endif
  LINALG::ANA::Vector tmp1(right_.RangeMap(), false);
  right_.Set(tmp1, 1.0);
  LINALG::ANA::Vector tmp2(op_->OperatorRangeMap(), false);
  int err = op_->Apply(tmp1, tmp2);
  if (err) dserror("LightWeightOperatorBase::Apply returned err=%d", err);
  v.Update(scale, tmp2, 1.0);
}



/*----------------------------------------------------------------------*
   vec dot vec (result is scalar)
 *----------------------------------------------------------------------*/
double LINALG::ANA::operator*(const LINALG::ANA::Vector& vec1, const LINALG::ANA::Vector& vec2)
{
#if DEBUGGING_ANA
  cout << "double operator* (const LINALG::ANA::Vector& vec1, const LINALG::ANA::Vector& vec2)"
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
double LINALG::ANA::operator*(const LINALG::ANA::LCBase& left, const LINALG::ANA::LCBase& right)
{
#if DEBUGGING_ANA
  cout << "double operator* (const LINALG::ANA::LCBase& left, const LINALG::ANA::LCBase& right)"
       << endl;
  fflush(stdout);
#endif
  LINALG::ANA::Vector tmpleft(left.RangeMap(), false);
  left.Set(tmpleft, 1.0);
  LINALG::ANA::Vector tmpright(right.RangeMap(), false);
  right.Set(tmpright, 1.0);
  double result;
  tmpleft.Dot(tmpright, &result);
  return result;
}

/*----------------------------------------------------------------------*
   vec dot lc (result is scalar)
 *----------------------------------------------------------------------*/
double LINALG::ANA::operator*(const LINALG::ANA::Vector& vec1, const LINALG::ANA::LCBase& right)
{
#if DEBUGGING_ANA
  cout << "double operator* (const LINALG::ANA::Vector& vec1, const LINALG::ANA::LCBase& right)"
       << endl;
  fflush(stdout);
#endif
  LINALG::ANA::Vector tmp(right.RangeMap(), false);
  right.Set(tmp, 1.0);
  double result;
  vec1.Dot(tmp, &result);
  return result;
}

/*----------------------------------------------------------------------*
   vec dot lcsv (result is scalar) (specialization)
 *----------------------------------------------------------------------*/
double LINALG::ANA::operator*(
    const LINALG::ANA::Vector& vec1, const LINALG::ANA::LC_s_times_vec& right)
{
#if DEBUGGING_ANA
  cout << "double operator* (const LINALG::ANA::Vector& vec1, const LINALG::ANA::LC_s_times_vec& "
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
double LINALG::ANA::operator*(
    const LINALG::ANA::LC_s_times_vec& left, const LINALG::ANA::LC_s_times_vec& right)
{
#if DEBUGGING_ANA
  cout << "double operator* (const LINALG::ANA::LC_s_times_vec& left, const "
          "LINALG::ANA::LC_s_times_vec& right)"
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
double LINALG::ANA::norm2(const LINALG::ANA::LCBase& lc)
{
  LINALG::ANA::Vector tmp(lc.RangeMap(), false);
  lc.Set(tmp, 1.0);
  return norm2(tmp);
}
double LINALG::ANA::norm1(const LINALG::ANA::LCBase& lc)
{
  LINALG::ANA::Vector tmp(lc.RangeMap(), false);
  lc.Set(tmp, 1.0);
  return norm1(tmp);
}
double LINALG::ANA::norminf(const LINALG::ANA::LCBase& lc)
{
  LINALG::ANA::Vector tmp(lc.RangeMap(), false);
  lc.Set(tmp, 1.0);
  return norminf(tmp);
}

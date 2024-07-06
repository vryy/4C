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
int Core::LinAlg::Ana::OperatorProduct::apply(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (!usetransposed_)  // A*B*x
  {
    Core::LinAlg::Ana::Vector tmp(right_->operator_range_map(), false);
    int err1 = right_->apply(X, tmp);
    if (err1) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err1);
    int err2 = left_->apply(tmp, Y);
    if (err2) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err2);
    return err1 - err2;
  }
  else  // (A*B)^T * x = B^T * A^T * x
  {
    Core::LinAlg::Ana::Vector tmp(left_->operator_domain_map(), false);
    const_cast<LightWeightOperatorBase&>(*left_).set_use_transpose(!left_->use_transpose());
    int err1 = left_->apply(X, tmp);
    if (err1) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err1);
    const_cast<LightWeightOperatorBase&>(*left_).set_use_transpose(!left_->use_transpose());
    const_cast<LightWeightOperatorBase&>(*right_).set_use_transpose(!right_->use_transpose());
    int err2 = right_->apply(tmp, Y);
    if (err2) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err2);
    const_cast<LightWeightOperatorBase&>(*right_).set_use_transpose(!right_->use_transpose());
    return err1 - err2;
  }
  return 0;
}
int Core::LinAlg::Ana::OperatorProduct::apply_inverse(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (!usetransposed_)  // (A*B)^{-1} * x = B^{-1} * A^{-1} * x
  {
    Core::LinAlg::Ana::Vector tmp(left_->operator_domain_map(), false);
    int err1 = left_->apply_inverse(X, tmp);
    if (err1) FOUR_C_THROW("LightWeightOperatorBase::ApplyInverse returned err=%d", err1);
    int err2 = right_->apply_inverse(tmp, Y);
    if (err2) FOUR_C_THROW("LightWeightOperatorBase::ApplyInverse returned err=%d", err2);
    return err1 - err2;
  }
  else  // [ (A*B)^T ]^{-1} * x = A^{T,-1} * B^{T,-1} * x
  {
    Core::LinAlg::Ana::Vector tmp(right_->operator_range_map(), false);
    const_cast<LightWeightOperatorBase&>(*right_).set_use_transpose(!right_->use_transpose());
    int err1 = right_->apply_inverse(X, tmp);
    if (err1) FOUR_C_THROW("LightWeightOperatorBase::ApplyInverse returned err=%d", err1);
    const_cast<LightWeightOperatorBase&>(*right_).set_use_transpose(!right_->use_transpose());
    const_cast<LightWeightOperatorBase&>(*left_).set_use_transpose(!left_->use_transpose());
    int err2 = left_->apply_inverse(tmp, Y);
    if (err2) FOUR_C_THROW("LightWeightOperatorBase::ApplyInverse returned err=%d", err2);
    const_cast<LightWeightOperatorBase&>(*left_).set_use_transpose(!left_->use_transpose());
    return err1 - err2;
  }
  return -1;
}



/*----------------------------------------------------------------------*
  implicit operator sum (Operator+Operator)*x
 *----------------------------------------------------------------------*/
int Core::LinAlg::Ana::OperatorSum::apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (!usetransposed_)  // (A+-B)*x = Ax +- Bx
  {
    Core::LinAlg::Ana::Vector tmp(right_->operator_range_map(), false);
    int err1 = right_->apply(X, tmp);
    if (err1) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err1);
    int err2 = left_->apply(X, Y);
    if (err2) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err2);
    Y.Update((double)sign_, tmp, 1.0);
    return err1 - err2;
  }
  else  // (A+-B)^T * x = A^T x +- B^T * x
  {
    Core::LinAlg::Ana::Vector tmp(right_->operator_domain_map(), false);
    const_cast<LightWeightOperatorBase&>(*right_).set_use_transpose(!right_->use_transpose());
    int err1 = right_->apply(X, tmp);
    if (err1) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err1);
    const_cast<LightWeightOperatorBase&>(*right_).set_use_transpose(!right_->use_transpose());
    const_cast<LightWeightOperatorBase&>(*left_).set_use_transpose(!left_->use_transpose());
    int err2 = left_->apply(X, Y);
    if (err2) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err2);
    const_cast<LightWeightOperatorBase&>(*left_).set_use_transpose(!left_->use_transpose());
    Y.Update((double)sign_, tmp, 1.0);
    return err1 - err2;
  }
  return 0;
}


/*----------------------------------------------------------------------*
  apply inverse of an operator
 *----------------------------------------------------------------------*/
int Core::LinAlg::Ana::OperatorInverse::apply(
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
  solver_.solve(rcpop, out, in, solver_params);

  return 0;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::Ana::LcVecPointwiseLc::update(
    Core::LinAlg::Ana::Vector& v, const double& scale) const
{
#if DEBUGGING_ANA
  cout << "Core::LinAlg::Ana::LC_vec_pointwise_lc::Update(Core::LinAlg::Ana::Vector& v, const "
          "double& scale)"
       << endl;
  fflush(stdout);
#endif
  // we can not avoid extra memory in this case
  Core::LinAlg::Ana::Vector tmp(right_.range_map(), false);
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
  Core::LinAlg::Ana::Vector tmp(right_.range_map(), false);
  right_.set(tmp, 1.0);
  v.Multiply(scale, v, tmp, 0.0);
}
void Core::LinAlg::Ana::LcLcPointwiseLc::update(
    Core::LinAlg::Ana::Vector& v, const double& scale) const
{
#if DEBUGGING_ANA
  cout << "Core::LinAlg::Ana::LC_lc_pointwise_lc::Update(Core::LinAlg::Ana::Vector& v, const "
          "double& scale)"
       << endl;
  fflush(stdout);
#endif
  Core::LinAlg::Ana::Vector tmp1(left_.range_map(), false);
  left_.set(tmp1, 1.0);
  Core::LinAlg::Ana::Vector tmp2(right_.range_map(), false);
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
  if (!v.Map().SameAs(op_->operator_range_map())) FOUR_C_THROW("Range maps don't match - fatal");
  if (!right_.vector().Map().SameAs(op_->operator_domain_map()))
    FOUR_C_THROW("Domain maps don't match - fatal");
#endif
  int err = op_->apply(right_.vector(), v);
  if (err) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err);
  if (scale * right_.scalar() != 1.0) v.Scale(scale * right_.scalar());
}
void Core::LinAlg::Ana::LcOperatorTimesLcsv::update(
    Core::LinAlg::Ana::Vector& v, const double& scale) const
{
#if DEBUGGING_ANA
  cout << "Core::LinAlg::Ana::LC_Operator_times_lcsv::Update(Core::LinAlg::Ana::Vector& v, const "
          "double& scale)"
       << endl;
  fflush(stdout);
#endif
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (!v.Map().SameAs(op_->operator_range_map())) FOUR_C_THROW("Range maps don't match - fatal");
  if (!right_.vector().Map().SameAs(op_->operator_domain_map()))
    FOUR_C_THROW("Domain maps don't match - fatal");
#endif
  Core::LinAlg::Ana::Vector tmp(op_->operator_range_map(), false);
  int err = op_->apply(right_.vector(), tmp);
  if (err) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err);
  v.Update(scale * right_.scalar(), tmp, 1.0);
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
  if (!v.Map().SameAs(op_->operator_range_map())) FOUR_C_THROW("Range maps don't match - fatal");
  if (!right_.range_map().SameAs(op_->operator_domain_map()))
    FOUR_C_THROW("Domain maps don't match - fatal");
#endif
  Core::LinAlg::Ana::Vector tmp(right_.range_map(), false);
  right_.set(tmp, scale);
  int err = op_->apply(tmp, v);
  if (err) FOUR_C_THROW("LightWeightOperatorBase::Apply returned err=%d", err);
}
void Core::LinAlg::Ana::LcOperatorTimesLc::update(
    Core::LinAlg::Ana::Vector& v, const double& scale) const
{
#if DEBUGGING_ANA
  cout << "Core::LinAlg::Ana::LC_Operator_times_lc::Update(Core::LinAlg::Ana::Vector& v, const "
          "double& scale)"
       << endl;
  fflush(stdout);
#endif
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (!v.Map().SameAs(op_->operator_range_map())) FOUR_C_THROW("Range maps don't match - fatal");
  if (!right_.range_map().SameAs(op_->operator_domain_map()))
    FOUR_C_THROW("Domain maps don't match - fatal");
#endif
  Core::LinAlg::Ana::Vector tmp1(right_.range_map(), false);
  right_.set(tmp1, 1.0);
  Core::LinAlg::Ana::Vector tmp2(op_->operator_range_map(), false);
  int err = op_->apply(tmp1, tmp2);
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
  Core::LinAlg::Ana::Vector tmpleft(left.range_map(), false);
  left.set(tmpleft, 1.0);
  Core::LinAlg::Ana::Vector tmpright(right.range_map(), false);
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
  Core::LinAlg::Ana::Vector tmp(right.range_map(), false);
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
  vec1.Dot(right.vector(), &result);
  result *= right.scalar();
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
  left.vector().Dot(right.vector(), &result);
  return (result * left.scalar() * right.scalar());
}

/*----------------------------------------------------------------------*
   norms
 *----------------------------------------------------------------------*/
double Core::LinAlg::Ana::norm2(const Core::LinAlg::Ana::LCBase& lc)
{
  Core::LinAlg::Ana::Vector tmp(lc.range_map(), false);
  lc.set(tmp, 1.0);
  return norm2(tmp);
}
double Core::LinAlg::Ana::norm1(const Core::LinAlg::Ana::LCBase& lc)
{
  Core::LinAlg::Ana::Vector tmp(lc.range_map(), false);
  lc.set(tmp, 1.0);
  return norm1(tmp);
}
double Core::LinAlg::Ana::norminf(const Core::LinAlg::Ana::LCBase& lc)
{
  Core::LinAlg::Ana::Vector tmp(lc.range_map(), false);
  lc.set(tmp, 1.0);
  return norminf(tmp);
}

FOUR_C_NAMESPACE_CLOSE

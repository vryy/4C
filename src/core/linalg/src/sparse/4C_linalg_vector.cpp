/*----------------------------------------------------------------------*/
/*! \file

\brief A wrapper for the Epetra_Vector

\level 0
*----------------------------------------------------------------------*/
#include "4C_linalg_vector.hpp"

#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN
Core::LinAlg::Vector::Vector(const Epetra_BlockMap& Map, bool zeroOut)
    : vector_(Epetra_Vector(Map, zeroOut))
{
}

Core::LinAlg::Vector::Vector(const Epetra_Vector& Source) : vector_(Source) {}

int Core::LinAlg::Vector::Norm1(double* Result) const { return vector_.Norm1(Result); }

int Core::LinAlg::Vector::Norm2(double* Result) const { return vector_.Norm2(Result); }

int Core::LinAlg::Vector::NormInf(double* Result) const { return vector_.NormInf(Result); }

int Core::LinAlg::Vector::MinValue(double* Result) const { return vector_.MinValue(Result); }

int Core::LinAlg::Vector::MaxValue(double* Result) const { return vector_.MaxValue(Result); }

int Core::LinAlg::Vector::MeanValue(double* Result) const { return vector_.MeanValue(Result); }

int Core::LinAlg::Vector::Dot(const Epetra_MultiVector& A, double* Result) const
{
  return vector_.Dot(A, Result);
}

int Core::LinAlg::Vector::Abs(const Epetra_MultiVector& A) { return vector_.Abs(A); }

int Core::LinAlg::Vector::Scale(double ScalarA, const Epetra_MultiVector& A)
{
  return vector_.Scale(ScalarA, A);
}

int Core::LinAlg::Vector::Update(double ScalarA, const Epetra_MultiVector& A, double ScalarThis)
{
  return vector_.Update(ScalarA, A, ScalarThis);
}

int Core::LinAlg::Vector::Update(double ScalarA, const Epetra_MultiVector& A, double ScalarB,
    const Epetra_MultiVector& B, double ScalarThis)
{
  return vector_.Update(ScalarA, A, ScalarB, B, ScalarThis);
}


int Core::LinAlg::Vector::Dot(const Vector& A, double* Result) const
{
  return vector_.Dot(A.get_Epetra_Vector(), Result);
}

int Core::LinAlg::Vector::Abs(const Vector& A) { return vector_.Abs(A.get_Epetra_Vector()); }

int Core::LinAlg::Vector::Scale(double ScalarA, const Vector& A)
{
  return vector_.Scale(ScalarA, A.get_Epetra_Vector());
}

int Core::LinAlg::Vector::Update(double ScalarA, const Vector& A, double ScalarThis)
{
  return vector_.Update(ScalarA, A.get_Epetra_Vector(), ScalarThis);
}

int Core::LinAlg::Vector::Update(
    double ScalarA, const Vector& A, double ScalarB, const Vector& B, double ScalarThis)
{
  return vector_.Update(ScalarA, A.get_Epetra_Vector(), ScalarB, B.get_Epetra_Vector(), ScalarThis);
}

int Core::LinAlg::Vector::PutScalar(double ScalarConstant)
{
  return vector_.PutScalar(ScalarConstant);
}

FOUR_C_NAMESPACE_CLOSE
// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_multi_vector.hpp"

#include "4C_linalg_vector.hpp"
#include "4C_utils_exceptions.hpp"


// Do not lint the file for identifier names, since the naming of the Wrapper functions follow the
// naming of the Core::LinAlg::MultiVector<double>

// NOLINTBEGIN(readability-identifier-naming)

FOUR_C_NAMESPACE_OPEN
template <typename T>
Core::LinAlg::MultiVector<T>::MultiVector(const Epetra_BlockMap& Map, int num_columns, bool zeroOut)
    : vector_(std::make_shared<Epetra_MultiVector>(Map, num_columns, zeroOut))
{
}


template <typename T>
Core::LinAlg::MultiVector<T>::MultiVector(const Epetra_MultiVector& source)
    : vector_(std::make_shared<Epetra_MultiVector>(source))
{
}

template <typename T>
Core::LinAlg::MultiVector<T>::MultiVector(const MultiVector& other)
    : vector_(std::make_shared<Epetra_MultiVector>(*other.vector_))
{
}


template <typename T>
Core::LinAlg::MultiVector<T>& Core::LinAlg::MultiVector<T>::operator=(const MultiVector& other)
{
  *vector_ = *other.vector_;
  return *this;
}

template <typename T>
int Core::LinAlg::MultiVector<T>::Norm1(double* Result) const
{
  return vector_->Norm1(Result);
}

template <typename T>
int Core::LinAlg::MultiVector<T>::Norm2(double* Result) const
{
  return vector_->Norm2(Result);
}

template <typename T>
int Core::LinAlg::MultiVector<T>::NormInf(double* Result) const
{
  return vector_->NormInf(Result);
}

template <typename T>
int Core::LinAlg::MultiVector<T>::MinValue(double* Result) const
{
  return vector_->MinValue(Result);
}

template <typename T>
int Core::LinAlg::MultiVector<T>::MaxValue(double* Result) const
{
  return vector_->MaxValue(Result);
}

template <typename T>
int Core::LinAlg::MultiVector<T>::MeanValue(double* Result) const
{
  return vector_->MeanValue(Result);
}

template <typename T>
int Core::LinAlg::MultiVector<T>::Dot(const MultiVector& A, double* Result) const
{
  return vector_->Dot(A, Result);
}

template <typename T>
int Core::LinAlg::MultiVector<T>::Abs(const MultiVector& A)
{
  return vector_->Abs(A);
}

template <typename T>
int Core::LinAlg::MultiVector<T>::Scale(double ScalarA, const MultiVector& A)
{
  return vector_->Scale(ScalarA, A);
}

template <typename T>
int Core::LinAlg::MultiVector<T>::Update(double ScalarA, const MultiVector& A, double ScalarThis)
{
  return vector_->Update(ScalarA, A, ScalarThis);
}

template <typename T>
int Core::LinAlg::MultiVector<T>::Update(
    double ScalarA, const MultiVector& A, double ScalarB, const MultiVector& B, double ScalarThis)
{
  return vector_->Update(ScalarA, A, ScalarB, *B.vector_, ScalarThis);
}

template <typename T>
int Core::LinAlg::MultiVector<T>::PutScalar(double ScalarConstant)
{
  return vector_->PutScalar(ScalarConstant);
}

template <typename T>
int Core::LinAlg::MultiVector<T>::ReplaceMap(const Epetra_BlockMap& map)
{
  for (auto& view : column_vector_view_) view->ReplaceMap(map);
  return vector_->ReplaceMap(map);
}

template <typename T>
Core::LinAlg::Vector<double>& Core::LinAlg::MultiVector<T>::operator()(int i)
{
  sync_view();
  FOUR_C_ASSERT(
      column_vector_view_.size() == static_cast<std::size_t>(NumVectors()), "Internal error.");
  return *column_vector_view_[i];
}

template <typename T>
const Core::LinAlg::Vector<double>& Core::LinAlg::MultiVector<T>::operator()(int i) const
{
  sync_view();
  FOUR_C_ASSERT(
      column_vector_view_.size() == static_cast<std::size_t>(NumVectors()), "Internal error.");
  return *column_vector_view_[i];
}

template <typename T>
void Core::LinAlg::MultiVector<T>::sync_view() const
{
  // Ensure that this only runs once.
  if (column_vector_view_.empty())
  {
    column_vector_view_.reserve(NumVectors());
    for (int i = 0; i < NumVectors(); ++i)
      column_vector_view_.emplace_back(Vector<T>::create_view(*(*vector_)(i)));
  }
}


// explicit instantiation
template class Core::LinAlg::MultiVector<double>;

FOUR_C_NAMESPACE_CLOSE

// NOLINTEND(readability-identifier-naming)
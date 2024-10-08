/*----------------------------------------------------------------------*/
/*! \file

\brief A wrapper for the Epetra_Vector

\level 0
*----------------------------------------------------------------------*/
#include "4C_linalg_vector.hpp"

#include "4C_linalg_multi_vector.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Vector.h>

// Do not lint the file for identifier names, since the naming of the Wrapper functions follow the
// naming of the Epetra_Vector

// NOLINTBEGIN(readability-identifier-naming)

FOUR_C_NAMESPACE_OPEN

template <typename T>
Core::LinAlg::Vector<T>::Vector(const Epetra_BlockMap& Map, bool zeroOut)
    : vector_(Teuchos::make_rcp<Epetra_Vector>(Map, zeroOut))
{
}

template <typename T>
Core::LinAlg::Vector<T>::Vector(const Epetra_Vector& Source)
    : vector_(Teuchos::make_rcp<Epetra_Vector>(Source))
{
}

template <typename T>
Core::LinAlg::Vector<T>::Vector(const Epetra_FEVector& Source)
    : vector_(Teuchos::make_rcp<Epetra_Vector>(Epetra_DataAccess::Copy, Source, 0))
{
}

template <typename T>
Core::LinAlg::Vector<T>::Vector(const Vector& other)
    : vector_(Teuchos::make_rcp<Epetra_Vector>(other.get_ref_of_Epetra_Vector()))
{
}


template <typename T>
Core::LinAlg::Vector<T>& Core::LinAlg::Vector<T>::operator=(const Vector& other)
{
  *vector_ = other.get_ref_of_Epetra_Vector();
  return *this;
}


template <typename T>
Core::LinAlg::Vector<T>::operator const Core::LinAlg::MultiVector<T>&() const
{
  sync_view();
  FOUR_C_ASSERT(multi_vector_view_ != Teuchos::null, "Internal error.");
  return *multi_vector_view_;
}


template <typename T>
Core::LinAlg::Vector<T>::operator Core::LinAlg::MultiVector<T>&()
{
  sync_view();
  FOUR_C_ASSERT(multi_vector_view_ != Teuchos::null, "Internal error.");
  return *multi_vector_view_;
}


template <typename T>
int Core::LinAlg::Vector<T>::Norm1(double* Result) const
{
  return vector_->Norm1(Result);
}

template <typename T>
int Core::LinAlg::Vector<T>::Norm2(double* Result) const
{
  return vector_->Norm2(Result);
}

template <typename T>
int Core::LinAlg::Vector<T>::NormInf(double* Result) const
{
  return vector_->NormInf(Result);
}

template <typename T>
int Core::LinAlg::Vector<T>::MinValue(double* Result) const
{
  return vector_->MinValue(Result);
}

template <typename T>
int Core::LinAlg::Vector<T>::MaxValue(double* Result) const
{
  return vector_->MaxValue(Result);
}

template <typename T>
int Core::LinAlg::Vector<T>::MeanValue(double* Result) const
{
  return vector_->MeanValue(Result);
}

template <typename T>
int Core::LinAlg::Vector<T>::Dot(const Epetra_MultiVector& A, double* Result) const
{
  return vector_->Dot(A, Result);
}

template <typename T>
int Core::LinAlg::Vector<T>::Abs(const Epetra_MultiVector& A)
{
  return vector_->Abs(A);
}

template <typename T>
int Core::LinAlg::Vector<T>::Scale(double ScalarA, const Epetra_MultiVector& A)
{
  return vector_->Scale(ScalarA, A);
}

template <typename T>
int Core::LinAlg::Vector<T>::Update(double ScalarA, const Epetra_MultiVector& A, double ScalarThis)
{
  return vector_->Update(ScalarA, A, ScalarThis);
}

template <typename T>
int Core::LinAlg::Vector<T>::Update(double ScalarA, const Epetra_MultiVector& A, double ScalarB,
    const Epetra_MultiVector& B, double ScalarThis)
{
  return vector_->Update(ScalarA, A, ScalarB, B, ScalarThis);
}


template <typename T>
int Core::LinAlg::Vector<T>::Dot(const Vector& A, double* Result) const
{
  return vector_->Dot(A, Result);
}

template <typename T>
int Core::LinAlg::Vector<T>::Abs(const Vector& A)
{
  return vector_->Abs(A);
}

template <typename T>
int Core::LinAlg::Vector<T>::Scale(double ScalarA, const Vector& A)
{
  return vector_->Scale(ScalarA, A);
}

template <typename T>
int Core::LinAlg::Vector<T>::Update(double ScalarA, const Vector& A, double ScalarThis)
{
  return vector_->Update(ScalarA, A, ScalarThis);
}

template <typename T>
int Core::LinAlg::Vector<T>::Update(
    double ScalarA, const Vector& A, double ScalarB, const Vector& B, double ScalarThis)
{
  return vector_->Update(ScalarA, A, ScalarB, B.get_ref_of_Epetra_Vector(), ScalarThis);
}

template <typename T>
int Core::LinAlg::Vector<T>::PutScalar(double ScalarConstant)
{
  return vector_->PutScalar(ScalarConstant);
}

template <typename T>
void Core::LinAlg::Vector<T>::sync_view() const
{
  // Only do this once.
  if (!multi_vector_view_) multi_vector_view_ = MultiVector<T>::create_view(*vector_);
}

// explicit instantiation
template class Core::LinAlg::Vector<double>;



Core::LinAlg::Vector<int>::Vector(const Epetra_BlockMap& map, bool zeroOut)
    : vector_(Teuchos::make_rcp<Epetra_IntVector>(map, zeroOut))
{
}

Core::LinAlg::Vector<int>::Vector(const Epetra_BlockMap& map, int* values)
    : vector_(Teuchos::make_rcp<Epetra_IntVector>(Epetra_DataAccess::Copy, map, values))
{
}

Core::LinAlg::Vector<int>::Vector(const Vector& other)
    : vector_(Teuchos::make_rcp<Epetra_IntVector>(*other.vector_))
{
}

Core::LinAlg::Vector<int>::Vector(Vector&& other) noexcept : vector_(std::move(other.vector_)) {}

Core::LinAlg::Vector<int>& Core::LinAlg::Vector<int>::operator=(const Vector& other)
{
  vector_ = Teuchos::make_rcp<Epetra_IntVector>(*other.vector_);
  return *this;
}

Core::LinAlg::Vector<int>& Core::LinAlg::Vector<int>::operator=(Vector&& other) noexcept
{
  vector_ = std::move(other.vector_);
  return *this;
}


int Core::LinAlg::Vector<int>::PutValue(int Value) { return vector_->PutValue(Value); }

int Core::LinAlg::Vector<int>::MaxValue() { return vector_->MaxValue(); }

int Core::LinAlg::Vector<int>::MinValue() { return vector_->MinValue(); }

void Core::LinAlg::Vector<int>::Print(std::ostream& os) const { vector_->Print(os); }


FOUR_C_NAMESPACE_CLOSE

// NOLINTEND(readability-identifier-naming)
/*----------------------------------------------------------------------*/
/*! \file

\brief A wrapper for the Epetra_Vector

\level 0
*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINALG_VECTOR_HPP
#define FOUR_C_LINALG_VECTOR_HPP


#include "4C_config.hpp"

#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPDecl.hpp>

// Do not lint the file for identifier names, since the naming of the Wrapper functions follow the
// naming of the Epetra_Vector

// NOLINTBEGIN(readability-identifier-naming)

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  // Sparse Vector which will replace the Epetra_Vector
  class Vector
  {
   public:
    /// Epetra related methods

    /// Basic vector constructor to create vector based on a map and initialize memory with zeros
    Vector(const Epetra_BlockMap& Map, bool zeroOut = true);

    /// Copy constructor from epetra to vector
    Vector(const Epetra_Vector& Source);

    //! get epetra vector
    Epetra_Vector get_Epetra_Vector() const { return vector_; }

    //! get teuchos pointer of epetra vector
    Teuchos::RCP<Epetra_Vector> get_ptr_of_Epetra_Vector() { return Teuchos::rcpFromRef(vector_); }

    //! get epetra multi vector from vector
    Epetra_MultiVector get_Epetra_MultiVector() const
    {
      return static_cast<Epetra_MultiVector>(vector_);
    }

    //! get pointer of epetra multi vector
    Teuchos::RCP<Epetra_MultiVector> get_ptr_of_Epetra_MultiVector()
    {
      return Teuchos::rcpFromRef(vector_);
    }

    //! Computes dot product of each corresponding pair of vectors.
    int Dot(const Epetra_MultiVector& A, double* Result) const;

    //! Puts element-wise absolute values of input Multi-vector in target.
    int Abs(const Epetra_MultiVector& A);

    //! Replace multi-vector values with scaled values of A, \e this = ScalarA*A.
    int Scale(double ScalarA, const Epetra_MultiVector& A);

    //! Update multi-vector values with scaled values of A, \e this = ScalarThis*\e this +
    //! ScalarA*A.
    int Update(double ScalarA, const Epetra_MultiVector& A, double ScalarThis);

    //! Update multi-vector with scaled values of A and B, \e this = ScalarThis*\e this + ScalarA*A
    //! + ScalarB*B.
    int Update(double ScalarA, const Epetra_MultiVector& A, double ScalarB,
        const Epetra_MultiVector& B, double ScalarThis);


    ///

    //! Compute 1-norm of each vector
    int Norm1(double* Result) const;

    //! Compute 2-norm of each vector
    int Norm2(double* Result) const;

    //! Compute Inf-norm of each vector
    int NormInf(double* Result) const;

    //! Compute minimum value of each vector
    int MinValue(double* Result) const;

    //! Compute maximum value of each vector
    int MaxValue(double* Result) const;

    //! Compute mean (average) value of each vector
    int MeanValue(double* Result) const;

    //! Scale the current values of a multi-vector, \e this = ScalarValue*\e this.
    int Scale(double ScalarValue);

    //! Computes dot product of each corresponding pair of vectors.
    int Dot(const Vector& A, double* Result) const;

    //! Puts element-wise absolute values of input Multi-vector in target.
    int Abs(const Vector& A);

    //! Replace multi-vector values with scaled values of A, \e this = ScalarA*A.
    int Scale(double ScalarA, const Vector& A);

    //! Update multi-vector values with scaled values of A, \e this = ScalarThis*\e this +
    //! ScalarA*A.
    int Update(double ScalarA, const Vector& A, double ScalarThis);

    //! Update multi-vector with scaled values of A and B, \e this = ScalarThis*\e this + ScalarA*A
    //! + ScalarB*B.
    int Update(double ScalarA, const Vector& A, double ScalarB, const Vector& B, double ScalarThis);

    //! Initialize all values in a multi-vector with const value.
    int PutScalar(double ScalarConstant);

    //! Element access function
    double& operator[](int index) { return vector_[index]; }

    //! Returns the address of the Epetra_BlockMap for this multi-vector.
    const Epetra_BlockMap& Map() const { return (vector_.Map()); };

    //! Returns the address of the Epetra_Comm for this multi-vector.
    const Epetra_Comm& Comm() const { return (vector_.Comm()); };

    //! Returns true if this multi-vector is distributed global, i.e., not local replicated.
    bool DistributedGlobal() const { return (vector_.Map().DistributedGlobal()); };

    //! Print method
    void Print(std::ostream& os) const { vector_.Print(os); }

    //! Returns the number of vectors in the multi-vector.
    int NumVectors() const { return vector_.NumVectors(); }

    //! Returns the local vector length on the calling processor of vectors in the multi-vector.
    int MyLength() const { return vector_.MyLength(); }

    //! Returns the global vector length of vectors in the multi-vector.
    int GlobalLength() const { return vector_.GlobalLength(); }


   private:
    Epetra_Vector vector_;
  };



}  // namespace Core::LinAlg



FOUR_C_NAMESPACE_CLOSE

// NOLINTEND(readability-identifier-naming)

#endif
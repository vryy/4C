/*----------------------------------------------------------------------*/
/*! \file
* \file linalg_utils_sparse_algebra_math.H

\brief A collection of algebraic mathematical methods for namespace LINALG

\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*-----------------------------------------------------------------------*/
#ifndef SRC_LINALG_LINALG_UTILS_SPARSE_ALGEBRA_MATH_H_
#define SRC_LINALG_LINALG_UTILS_SPARSE_ALGEBRA_MATH_H_

#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_exporter.H"
#include "linalg_blocksparsematrix.H"

namespace LINALG
{
  /*!
    \brief Make Epetra_SerialDenseMatrix symmetric by averaging upper and lower traingular indices

    \param A (in/out): Matrix to be symmetrised
    */
  void SymmetriseMatrix(Epetra_SerialDenseMatrix& A);

  /*!
     \brief Add a (transposed) Epetra_CrsMatrix to another: B = B*scalarB + A(^T)*scalarA

     Add one matrix to another.

     The matrix B may or may not be completed. If B is completed, no new elements can be
     inserted and the addition only succeeds in case the sparsity pattern of B is a superset of
     the sparsity pattern of A (otherwise: dserror).

     Performance characterization: If B is filled (completed), this function is pretty fast,
     typically on the order of two to four matrix-vector products with B. The case where B is
     un-filled runs much slower (on the order of up to 100 matrix-vector products).

     Sparsity patterns of A and B need not match and A and B can be
     nonsymmetric in value and pattern.

     Row map of A has to be a processor-local subset of the row map of B.

     Note that this is a true parallel add, even in the transposed case!

     \param A          (in)     : Matrix to add to B (must have Filled()==true)
     \param transposeA (in)     : flag indicating whether transposed of A should be used
     \param scalarA    (in)     : scaling factor for A
     \param B          (in/out) : Matrix to be added to (must have Filled()==false)
     \param scalarB    (in)     : scaling factor for B
     */
  void Add(const Epetra_CrsMatrix& A, const bool transposeA, const double scalarA,
      Epetra_CrsMatrix& B, const double scalarB);

  /*!
   \brief Add a (transposed) Epetra_CrsMatrix to another: B = B*scalarB + A(^T)*scalarA

   Add one matrix to another.

   The matrix B may or may not be completed. If B is completed, no new elements can be
   inserted and the addition only succeeds in case the sparsity pattern of B is a superset of
   the sparsity pattern of A (otherwise: dserror).

   Performance characterization: If B is filled (completed), this function is pretty fast,
   typically on the order of two to four matrix-vector products with B. The case where B is
   un-filled runs much slower (on the order of up to 100 matrix-vector products).

   Sparsity patterns of A and B need not match and A and B can be
   nonsymmetric in value and pattern.

   Row map of A has to be a processor-local subset of the row map of B.


   Note that this is a true parallel add, even in the transposed case!
   This is the Teuchos::RCP wrapper of the above method.

   \param A          (in)     : Matrix to add to B (must have Filled()==true)
   \param transposeA (in)     : flag indicating whether transposed of A should be used
   \param scalarA    (in)     : scaling factor for A
   \param B          (in/out) : Matrix to be added to (must have Filled()==false)
   \param scalarB    (in)     : scaling factor for B
   */
  inline void Add(const Teuchos::RCP<Epetra_CrsMatrix> A, const bool transposeA,
      const double scalarA, Teuchos::RCP<Epetra_CrsMatrix> B, const double scalarB)
  {
    LINALG::Add(*A, transposeA, scalarA, *B, scalarB);
    return;
  }

  /*!
   \brief Add a (transposed) Epetra_CrsMatrix to a LINALG::SparseMatrix: B = B*scalarB +
   A(^T)*scalarA

   Add one matrix to another.

   As opposed to the other Add() functions, this method can handle both the case where
   matrix B is fill-completed (for performance reasons) but does not have to.
   If B is completed and new matrix elements are detected, the matrix is un-completed and
   rebuild internally (expensive).

   Sparsity patterns of A and B need not match and A and B can be
   nonsymmetric in value and pattern.

   Row map of A has to be a processor-local subset of the row map of B.

   Note that this is a true parallel add, even in the transposed case!

   \param A          (in)     : Matrix to add to B (must have Filled()==true)
   \param transposeA (in)     : flag indicating whether transposed of A should be used
   \param scalarA    (in)     : scaling factor for A
   \param B          (in/out) : Matrix to be added to (must have Filled()==false)
   \param scalarB    (in)     : scaling factor for B
   */
  void Add(const Epetra_CrsMatrix& A, const bool transposeA, const double scalarA,
      LINALG::SparseMatrixBase& B, const double scalarB);

  /*!
   \brief Compute transposed matrix of an Epetra_CrsMatrix explicitly

   Returns Teuchos::RCP to the transposed matrix of the input matrix A.

   \param A          (in)     : Matrix to transpose (must have Filled()==true)

   */
  Teuchos::RCP<Epetra_CrsMatrix> Transpose(const Epetra_CrsMatrix& A);

  /*!
   \brief Compute transposed matrix of an Epetra_CrsMatrix explicitly

   Returns Teuchos::RCP to the transposed matrix of the input matrix A.
   This is the Teuchos::RCP wrapper of the above method.

   \param A          (in)     : Matrix to transpose (must have Filled()==true)

   */
  inline Teuchos::RCP<Epetra_CrsMatrix> Transpose(const Teuchos::RCP<Epetra_CrsMatrix> A)
  {
    return LINALG::Transpose(*A);
  }

  /*!
   \brief Multiply a (transposed) Epetra_CrsMatrix with another (transposed): C = A(^T)*B(^T)

   Multiply one matrix with another. Both matrices must be completed. Sparsity
   Respective Range, Row and Domain maps of A(^T) and B(^T) have to match.

   Note that this is a true parallel multiplication, even in the transposed case!

   \param A          (in)     : Matrix to multiply with B (must have Filled()==true)
   \param transA     (in)     : flag indicating whether transposed of A should be used
   \param B          (in)     : Matrix to multiply with A (must have Filled()==true)
   \param transB     (in)     : flag indicating whether transposed of B should be used
   \param complete   (in)     : flag indicating whether FillComplete should be called on C upon
   exit, (defaults to true) \return Matrix product A(^T)*B(^T)
   */
  Teuchos::RCP<Epetra_CrsMatrix> Multiply(const Epetra_CrsMatrix& A, bool transA,
      const Epetra_CrsMatrix& B, bool transB, bool complete = true);

  /*!
   \brief Multiply a (transposed) Epetra_CrsMatrix with another (transposed): C = A(^T)*B(^T)

   Multiply one matrix with another. Both matrices must be completed. Sparsity
   Respective Range, Row and Domain maps of A(^T) and B(^T) have to match.

   Note that this is a true parallel multiplication, even in the transposed case!
   This is the Teuchos::RCP wrapper of the above method.

   \param A          (in)     : Matrix to multiply with B (must have Filled()==true)
   \param transA     (in)     : flag indicating whether transposed of A should be used
   \param B          (in)     : Matrix to multiply with A (must have Filled()==true)
   \param transB     (in)     : flag indicating whether transposed of B should be used
   \param complete   (in)     : flag indicating whether FillComplete should be called on C upon
   exit, (defaults to true) \return Matrix product A(^T)*B(^T)
   */
  inline Teuchos::RCP<Epetra_CrsMatrix> Multiply(const Teuchos::RCP<Epetra_CrsMatrix>& A,
      bool transA, const Teuchos::RCP<Epetra_CrsMatrix>& B, bool transB, bool complete = true)
  {
    return Multiply(*A, transA, *B, transB, complete);
  }

  /*!
   \brief Triple matrix product: D = A(^T)*B(^T)*C(^T)

   Multiply one matrix with another. All input matrices must be completed. Sparsity
   Respective Range, Row and Domain maps of A(^T) and B(^T) C(^T) have to match.

   Note that this is a true parallel multiplication, even in the transposed case!

   \param A          (in)     : Matrix to multiply with B (must have Filled()==true)
   \param transA     (in)     : flag indicating whether transposed of A should be used
   \param B          (in)     : Matrix to multiply with C (must have Filled()==true)
   \param transB     (in)     : flag indicating whether transposed of B should be used
   \param C          (in)     : Matrix C (must have Filled()==true)
   \param transC     (in)     : flag indicating whether transposed of C should be used
   \param complete   (in)     : flag indicating whether FillComplete should be called on C upon
   exit, (defaults to true) \return Matrix product A(^T)*B(^T)*C(^T)
   */
  Teuchos::RCP<Epetra_CrsMatrix> Multiply(const Epetra_CrsMatrix& A, bool transA,
      const Epetra_CrsMatrix& B, bool transB, const Epetra_CrsMatrix& C, bool transC,
      bool complete = true);

  /*!
   \brief Triple matrix product: D = A(^T)*B(^T)*C(^T)

   Multiply one matrix with another. All input matrices must be completed. Sparsity
   Respective Range, Row and Domain maps of A(^T) and B(^T) C(^T) have to match.

   Note that this is a true parallel multiplication, even in the transposed case!
   This is the Teuchos::RCP wrapper of the above method.

   \param A          (in)     : Matrix to multiply with B (must have Filled()==true)
   \param transA     (in)     : flag indicating whether transposed of A should be used
   \param B          (in)     : Matrix to multiply with C (must have Filled()==true)
   \param transB     (in)     : flag indicating whether transposed of B should be used
   \param C          (in)     : Matrix C (must have Filled()==true)
   \param transC     (in)     : flag indicating whether transposed of C should be used
   \param complete   (in)     : flag indicating whether FillComplete should be called on C upon
   exit, (defaults to true) \return Matrix product A(^T)*B(^T)*C(^T)
   */
  inline Teuchos::RCP<Epetra_CrsMatrix> Multiply(const Teuchos::RCP<Epetra_CrsMatrix>& A,
      bool transA, const Teuchos::RCP<Epetra_CrsMatrix>& B, bool transB,
      const Teuchos::RCP<Epetra_CrsMatrix>& C, bool transC, bool complete = true)
  {
    return Multiply(*A, transA, *B, transB, *C, transC, complete);
  }

  /// Returns a dim x dim identity matrix
  template <unsigned int dim>
  void IdentityMatrix(LINALG::Matrix<dim, dim>& identity_matrix)
  {
    for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
      {
        if (i == j)
          identity_matrix(i, j) = 1.0;
        else
          identity_matrix(i, j) = 0.0;
      }
    }

    return;
  }

}  // namespace LINALG

#endif /* SRC_LINALG_LINALG_UTILS_SPARSE_ALGEBRA_MATH_H_ */

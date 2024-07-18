/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of algebraic mathematical methods for namespace Core::LinAlg

\level 0
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINALG_UTILS_SPARSE_ALGEBRA_MATH_HPP
#define FOUR_C_LINALG_UTILS_SPARSE_ALGEBRA_MATH_HPP

#include "4C_config.hpp"

#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"

#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  /*!
   \brief Add a (transposed) Epetra_CrsMatrix to another: B = B*scalarB + A(^T)*scalarA

   Add one matrix to another.

   The matrix B may or may not be completed. If B is completed, no new elements can be
   inserted and the addition only succeeds in case the sparsity pattern of B is a superset of
   the sparsity pattern of A (otherwise: FOUR_C_THROW).

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
   the sparsity pattern of A (otherwise: FOUR_C_THROW).

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
    Core::LinAlg::Add(*A, transposeA, scalarA, *B, scalarB);
  }

  /*!
   \brief Add a (transposed) Epetra_CrsMatrix to a Core::LinAlg::SparseMatrix: B = B*scalarB +
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
      Core::LinAlg::SparseMatrixBase& B, const double scalarB);

  /*!
   \brief Multiply a (transposed) sparse matrix with another (transposed): C = A(^T)*B(^T)

   Multiply one matrix with another. Both matrices must be completed. Sparsity
   Respective Range, Row and Domain maps of A(^T) and B(^T) have to match.

   \note that this is a true parallel multiplication, even in the transposed case!

   \note Does call complete on C upon exit by default.

   \param A          (in) : Matrix to multiply with B (must have Filled()==true)
   \param transA     (in) : flag indicating whether transposed of A should be used
   \param B          (in) : Matrix to multiply with A (must have Filled()==true)
   \param transB     (in) : flag indicating whether transposed of B should be used
   \param complete   (in) : flag indicating whether fill_complete should be called on C upon
                            exit, (defaults to true)
   \return Matrix product A(^T)*B(^T)
   */
  Teuchos::RCP<SparseMatrix> MatrixMultiply(
      const SparseMatrix& A, bool transA, const SparseMatrix& B, bool transB, bool complete = true);

  /*!
   \brief Multiply a (transposed) sparse matrix with another (transposed): C = A(^T)*B(^T)

   Multiply one matrix with another. Both matrices must be completed. Sparsity
   Respective Range, Row and Domain maps of A(^T) and B(^T) have to match.

   \note that this is a true parallel multiplication, even in the transposed case!

   \note Does call complete on C upon exit by default.

   \note In this version the flags explicitdirichlet and savegraph must be handed in.
   Thus, they can be defined explicitly, while in the standard version of MatrixMultiply()
   above, result matrix C automatically inherits these flags from input matrix A.

   \param A                 (in) : Matrix to multiply with B (must have Filled()==true)
   \param transA            (in) : flag indicating whether transposed of A should be used
   \param B                 (in) : Matrix to multiply with A (must have Filled()==true)
   \param transB            (in) : flag indicating whether transposed of B should be used
   \param explicitdirichlet (in) : flag deciding on explicitdirichlet flag of C
   \param savegraph         (in) : flag deciding on savegraph flag of C
   \param complete          (in) : flag indicating whether fill_complete should be called on C upon
                                   exit, (defaults to true)
   \return Matrix product A(^T)*B(^T)
   */
  Teuchos::RCP<SparseMatrix> MatrixMultiply(const SparseMatrix& A, bool transA,
      const SparseMatrix& B, bool transB, bool explicitdirichlet, bool savegraph,
      bool complete = true);

}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif

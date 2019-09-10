/*----------------------------------------------------------------------*/
/*! \file

\brief matrix multiplication header

\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
*/
/*----------------------------------------------------------------------*/

#include <Epetra_CrsMatrix.h>

namespace LINALG
{
  // forward declarations

  class SparseMatrix;
  class SparseOperator;


  /*!
  \brief Matrix-Matrix Multiply C = A*B using ML

  Depending on the structure of your matrix, this method is potentially
  (significantly) faster than the LINALG::Multiply implementations.
  A factor of 5 - 10 in speed has been observed when comparing.
  It uses ML instead of EpetraExt for the multiplication kernel.

                           ***Warning***
           This method is not well tested (yet) and there is some risk of
           failure AND wrong results without notice. This depends
           on structural properties of the input matrices A and B.
           First:
           A,B should NOT contain columns in the
           column and/or domain map on a processor,
           where there is no nonzero entry in that column on a specific processor.
           This will definitely lead to failure but is a very rare case.
           Second:
           There should not be a row in row/rangemap with no
           nonzero entry in that row on any processor.
           Third:
           Very weirdo column orderings might lead to trouble, but
           scattered row/column ordering as usually the case in PDE
           problems should work fine.
           Fourth:
           If you have a choice how to formulate your problem:

  \note The Mat-Mat-Mult kernel used herein
        was optimized for algebraic multigrid. This means, the method
        will be significantly faster when B is sparser and/or smaller than A.

  \note There is no natural transpose multiply in ML. So if you need a
        transpose-Multiply you have to transpose outside (or write a wrapper
        for this method)

  \note ML wipes all exact zero entries in the product. So its perfectly
        ok if your product has less nonzeros than with the other
        LINALG::Multiply methods.

  \param Aorig (in)       : Matrix A for C=A*B
  \param Borig (in)       : Matrix B for C=A*B
  \param complete (in): flag indicating whether C shall be called FillComplete.
                        There intentionally is NO default value for this here.
  */
  Teuchos::RCP<SparseMatrix> MLMultiply(const Epetra_CrsMatrix& Aorig,
      const Epetra_CrsMatrix& Borig, bool explicitdirichlet, bool savegraph, bool complete);

  /*!
  \brief Matrix-Matrix Multiply C = A*B using ML

  \sa MLMultiply(const Epetra_CrsMatrix& A,const Epetra_CrsMatrix& B,bool complete);

  */
  Teuchos::RCP<SparseMatrix> MLMultiply(const SparseMatrix& A, const SparseMatrix& B,
      bool explicitdirichlet, bool savegraph, bool complete);

  /*!
  \brief Matrix-Matrix Multiply C = A*B using ML

  \sa MLMultiply(const Epetra_CrsMatrix& A,const Epetra_CrsMatrix& B,bool complete);

  */
  Teuchos::RCP<SparseMatrix> MLMultiply(
      const SparseMatrix& A, const SparseMatrix& B, bool complete);


  /// Multiply a (transposed) matrix with another (transposed): C = A(^T)*B(^T)
  /*!
    Multiply one matrix with another. Both matrices must be completed.
    Respective Range, Row and Domain maps of A(^T) and B(^T) have to match.

    \note This is a true parallel multiplication, even in the transposed case.

    \note Does call complete on C upon exit by default.

    \note Uses ML as multiplication kernel, not EpetraExt.

    \note In this version the flags explicitdirichlet and savegraph must be handed in.
          Thus, they can be defined explicitly, while in the standard version of Multipliy()
          above, result matrix C automatically inherits these flags from input matrix A

    \param A              (in)     : Matrix to multiply with B (must have Filled()==true)
    \param transA         (in)     : flag indicating whether transposed of A should be used
    \param B              (in)     : Matrix to multiply with A (must have Filled()==true)
    \param transB         (in)     : flag indicating whether transposed of B should be used
    \param explicitdirichlet (in)  : flag deciding on explicitdirichlet flag of C
    \param savegraph      (in)     : flag deciding on savegraph flag of C
    \param completeoutput (in)     : flag indicating whether Complete(...) shall be called on C upon
    output \return Matrix product A(^T)*B(^T)
  */
  Teuchos::RCP<SparseMatrix> MLMultiply(const SparseMatrix& A, bool transA, const SparseMatrix& B,
      bool transB, bool explicitdirichlet, bool savegraph, bool completeoutput);

}  // namespace LINALG

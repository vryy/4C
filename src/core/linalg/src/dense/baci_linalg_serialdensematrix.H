/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of serial dense matrix wrapper class

\level 0
*/
/*----------------------------------------------------------------------*/
#ifndef BACI_LINALG_SERIALDENSEMATRIX_H
#define BACI_LINALG_SERIALDENSEMATRIX_H


#include <Teuchos_SerialDenseMatrix.hpp>


namespace CORE::LINALG
{
  /*!
  \brief A class that wraps Teuchos::SerialDenseMatrix

       This is done in favor of typedef to allow forward declaration
  */
  class SerialDenseMatrix : public Teuchos::SerialDenseMatrix<int, double>
  {
   public:
    /// Base type definition
    using Base = Teuchos::SerialDenseMatrix<int, double>;

    /// Using the base class constructor
    using Base::SerialDenseMatrix;

    /*!
       \brief Standard Copy Constructor wraps
        Teuchos::SerialDenseMatrix(const SerialDenseMatrix& Source);
       This allows to use the CORE::LINALG::SerialDenseVector in place of
       CORE::LINALG::SerialDenseMatrix at the cost of a copy operation
    */
    SerialDenseMatrix(const Base& Source, Teuchos::ETransp trans = Teuchos::NO_TRANS)
        : Base(Source, trans)
    {
    }
  };

  /*!
    \brief Update matrix components with scaled values of A,
           B = alpha*A + beta*B
    */
  void Update(double alpha, const SerialDenseMatrix& A, double beta, SerialDenseMatrix& B);

  /*!
   * \brief Zero out first n elements of the matrix
   */
  void Zero(SerialDenseMatrix& mat, int Length);

  /*!
    \brief Determinant Computation using the Sarrus rule.
    Internal computation is based on the long double data type (80/128 bit depending on platform)
    this allows for higher precision when used on a bigger matrix. Long double is also used for
    output, therefore it has to be expicitly casted to a double if used in a "double only" context.
    */
  long double Det_long(const SerialDenseMatrix& mat);

  /*!
   * \brief Utility function to copy the data to SerialDenseMatrix
   */
  void copy(const double* vec, SerialDenseMatrix::Base& mat);

  // output stream operator
  inline std::ostream& operator<<(std::ostream& out, const SerialDenseMatrix& mat)
  {
    mat.print(out);
    return out;
  }

  /*!
   \brief General matrix-matrix multiplication C = A*B.
   It assumes that the dimension is compatible

   \param A (in):        Matrix A
   \param B (in):        Matrix B
   \param C (out):       Matrix C
   */
  template <typename Tmat1, typename Tmat2, typename Tmat3>
  int AxB(const Tmat1& A, const Tmat2& B, Tmat3& C)
  {
    for (int i = 0; i < C.numRows(); ++i)
    {
      for (int j = 0; j < C.numCols(); ++j)
      {
        double cij = 0.0;
        for (int k = 0; k < A.numCols(); ++k) cij += A(i, k) * B(k, j);
        C(i, j) = cij;
      }
    }
    return 0;
  }
}  // namespace CORE::LINALG

#endif  // BACI_LINALG_SERIALDENSEMATRIX_H

/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of serial dense matrix wrapper class

\level 0
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINALG_SERIALDENSEMATRIX_HPP
#define FOUR_C_LINALG_SERIALDENSEMATRIX_HPP


#include "4C_config.hpp"

#include <Teuchos_SerialDenseMatrix.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
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
       This allows to use the Core::LinAlg::SerialDenseVector in place of
       Core::LinAlg::SerialDenseMatrix at the cost of a copy operation
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
  void Zero(SerialDenseMatrix& mat, int length);

  /*!
    \brief Determinant Computation using the Sarrus rule.
    Internal computation is based on the long double data type (80/128 bit depending on platform)
    this allows for higher precision when used on a bigger matrix. Long double is also used for
    output, therefore it has to be expicitly casted to a double if used in a "double only" context.
    */
  long double Det_long(const SerialDenseMatrix& matrix);

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
}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif

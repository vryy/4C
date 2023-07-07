/*---------------------------------------------------------------------*/
/*! \file

\brief Implementation of least squares by Sudhakar for Moment-fitting

\level 3


*----------------------------------------------------------------------*/

#ifndef CUT_LEAST_SQUARES_H
#define CUT_LEAST_SQUARES_H

#include <vector>
#include <Epetra_SerialDenseSolver.h>
#include "linalg_utils_sparse_algebra_math.H"
#include "linalg_serialdensematrix.H"
#include "linalg_serialdensevector.H"
namespace CORE::GEO
{
  namespace CUT
  {
    /*!
    \brief This class solves the system of equations using linear least squares method
    */
    class LeastSquares
    {
     public:
      LeastSquares(std::vector<std::vector<double>> matri, Epetra_SerialDenseVector sourc)
          : matri_(matri), sourc_(sourc)
      {
      }
      ~LeastSquares(){};

      /*!
      \brief Performs the linear least squares. The resulting equations are solved using
      non-iterative method because the least squares produces full matrix
      */
      Epetra_SerialDenseVector linear_least_square();

     private:
      /*!
      \brief Construct the square matrix by multiplying the transpose of the original system matrix
      */
      Epetra_SerialDenseMatrix get_square_matrix(Epetra_SerialDenseVector &rhs);

      std::vector<std::vector<double>> matri_;
      Epetra_SerialDenseVector sourc_;
      Epetra_SerialDenseVector unknown_;
    };
  }  // namespace CUT
}  // namespace CORE::GEO

#endif

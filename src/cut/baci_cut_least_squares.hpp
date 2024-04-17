/*---------------------------------------------------------------------*/
/*! \file

\brief Implementation of least squares by Sudhakar for Moment-fitting

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_LEAST_SQUARES_HPP
#define FOUR_C_CUT_LEAST_SQUARES_HPP

#include "baci_config.hpp"

#include "baci_linalg_serialdensematrix.hpp"
#include "baci_linalg_serialdensevector.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"

#include <Teuchos_SerialDenseSolver.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN
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
      LeastSquares(std::vector<std::vector<double>> matri, CORE::LINALG::SerialDenseVector sourc)
          : matri_(matri), sourc_(sourc)
      {
      }

      /*!
      \brief Performs the linear least squares. The resulting equations are solved using
      non-iterative method because the least squares produces full matrix
      */
      CORE::LINALG::SerialDenseVector linear_least_square();

     private:
      /*!
      \brief Construct the square matrix by multiplying the transpose of the original system matrix
      */
      CORE::LINALG::SerialDenseMatrix get_square_matrix(CORE::LINALG::SerialDenseVector &rhs);

      std::vector<std::vector<double>> matri_;
      CORE::LINALG::SerialDenseVector sourc_;
      CORE::LINALG::SerialDenseVector unknown_;
    };
  }  // namespace CUT
}  // namespace CORE::GEO

FOUR_C_NAMESPACE_CLOSE

#endif

/*---------------------------------------------------------------------*/
/*! \file

\brief Implementation of least squares by Sudhakar for Moment-fitting

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_LEAST_SQUARES_HPP
#define FOUR_C_CUT_LEAST_SQUARES_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include <Teuchos_SerialDenseSolver.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN
namespace Core::Geo
{
  namespace Cut
  {
    /*!
    \brief This class solves the system of equations using linear least squares method
    */
    class LeastSquares
    {
     public:
      LeastSquares(std::vector<std::vector<double>> matri, Core::LinAlg::SerialDenseVector sourc)
          : matri_(matri), sourc_(sourc)
      {
      }

      /*!
      \brief Performs the linear least squares. The resulting equations are solved using
      non-iterative method because the least squares produces full matrix
      */
      Core::LinAlg::SerialDenseVector linear_least_square();

     private:
      /*!
      \brief Construct the square matrix by multiplying the transpose of the original system matrix
      */
      Core::LinAlg::SerialDenseMatrix get_square_matrix(Core::LinAlg::SerialDenseVector &rhs);

      std::vector<std::vector<double>> matri_;
      Core::LinAlg::SerialDenseVector sourc_;
      Core::LinAlg::SerialDenseVector unknown_;
    };
  }  // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif

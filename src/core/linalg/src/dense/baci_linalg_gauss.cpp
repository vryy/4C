/*----------------------------------------------------------------------*/
/*! \file

\brief Gauss elimination for small nxn systems

\level 1

*/
/*----------------------------------------------------------------------*/

#include "baci_linalg_gauss_templates.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CORE::LINALG
{
  template double gaussElimination<true, 1, double>(
      CORE::LINALG::Matrix<1, 1, double>& A,  ///< (in)    : system matrix
      CORE::LINALG::Matrix<1, 1, double>& b,  ///< (in)    : right-hand-side
      CORE::LINALG::Matrix<1, 1, double>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<false, 1, double>(
      CORE::LINALG::Matrix<1, 1>& A,  ///< (in)    : system matrix
      CORE::LINALG::Matrix<1, 1>& b,  ///< (in)    : right-hand-side
      CORE::LINALG::Matrix<1, 1>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<true, 2, double>(
      CORE::LINALG::Matrix<2, 2>& A,  ///< (in)    : system matrix
      CORE::LINALG::Matrix<2, 1>& b,  ///< (in)    : right-hand-side
      CORE::LINALG::Matrix<2, 1>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<false, 2, double>(
      CORE::LINALG::Matrix<2, 2>& A,  ///< (in)    : system matrix
      CORE::LINALG::Matrix<2, 1>& b,  ///< (in)    : right-hand-side
      CORE::LINALG::Matrix<2, 1>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<true, 3, double>(
      CORE::LINALG::Matrix<3, 3>& A,  ///< (in)    : system matrix
      CORE::LINALG::Matrix<3, 1>& b,  ///< (in)    : right-hand-side
      CORE::LINALG::Matrix<3, 1>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<false, 3, double>(
      CORE::LINALG::Matrix<3, 3>& A,  ///< (in)    : system matrix
      CORE::LINALG::Matrix<3, 1>& b,  ///< (in)    : right-hand-side
      CORE::LINALG::Matrix<3, 1>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<true, 4, double>(
      CORE::LINALG::Matrix<4, 4>& A,  ///< (in)    : system matrix
      CORE::LINALG::Matrix<4, 1>& b,  ///< (in)    : right-hand-side
      CORE::LINALG::Matrix<4, 1>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<false, 4, double>(
      CORE::LINALG::Matrix<4, 4>& A,  ///< (in)    : system matrix
      CORE::LINALG::Matrix<4, 1>& b,  ///< (in)    : right-hand-side
      CORE::LINALG::Matrix<4, 1>& x   ///< (out)   : solution vector
  );

  template double scaledGaussElimination<2>(
      CORE::LINALG::Matrix<2, 2>& A,  ///< (in)    : system matrix
      CORE::LINALG::Matrix<2, 1>& b,  ///< (in)    : right-hand-side
      CORE::LINALG::Matrix<2, 1>& x   ///< (out)   : solution vector
  );
  template double scaledGaussElimination<3>(
      CORE::LINALG::Matrix<3, 3>& A,  ///< (in)    : system matrix
      CORE::LINALG::Matrix<3, 1>& b,  ///< (in)    : right-hand-side
      CORE::LINALG::Matrix<3, 1>& x   ///< (out)   : solution vector
  );
  template double scaledGaussElimination<4>(
      CORE::LINALG::Matrix<4, 4>& A,  ///< (in)    : system matrix
      CORE::LINALG::Matrix<4, 1>& b,  ///< (in)    : right-hand-side
      CORE::LINALG::Matrix<4, 1>& x   ///< (out)   : solution vector
  );

}  // namespace CORE::LINALG

FOUR_C_NAMESPACE_CLOSE

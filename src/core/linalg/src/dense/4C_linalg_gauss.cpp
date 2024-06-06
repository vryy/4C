/*----------------------------------------------------------------------*/
/*! \file

\brief Gauss elimination for small nxn systems

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_linalg_gauss_templates.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  template double gaussElimination<true, 1, double>(
      Core::LinAlg::Matrix<1, 1, double>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<1, 1, double>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<1, 1, double>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<false, 1, double>(
      Core::LinAlg::Matrix<1, 1>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<1, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<1, 1>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<true, 2, double>(
      Core::LinAlg::Matrix<2, 2>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<2, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<2, 1>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<false, 2, double>(
      Core::LinAlg::Matrix<2, 2>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<2, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<2, 1>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<true, 3, double>(
      Core::LinAlg::Matrix<3, 3>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<3, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<3, 1>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<false, 3, double>(
      Core::LinAlg::Matrix<3, 3>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<3, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<3, 1>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<true, 4, double>(
      Core::LinAlg::Matrix<4, 4>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<4, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<4, 1>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<false, 4, double>(
      Core::LinAlg::Matrix<4, 4>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<4, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<4, 1>& x   ///< (out)   : solution vector
  );

  template double scaledGaussElimination<2>(
      Core::LinAlg::Matrix<2, 2>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<2, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<2, 1>& x   ///< (out)   : solution vector
  );
  template double scaledGaussElimination<3>(
      Core::LinAlg::Matrix<3, 3>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<3, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<3, 1>& x   ///< (out)   : solution vector
  );
  template double scaledGaussElimination<4>(
      Core::LinAlg::Matrix<4, 4>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<4, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<4, 1>& x   ///< (out)   : solution vector
  );

}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

/*----------------------------------------------------------------------*/
/*! \file

\brief Gauss elimination for small nxn systems

\level 1

*/
/*----------------------------------------------------------------------*/

#include "linalg_gauss_templates.H"

namespace LINALG
{
  template double gaussElimination<true, 1, double>(
      LINALG::Matrix<1, 1, double>& A,  ///< (in)    : system matrix
      LINALG::Matrix<1, 1, double>& b,  ///< (in)    : right-hand-side
      LINALG::Matrix<1, 1, double>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<false, 1, double>(
      LINALG::Matrix<1, 1>& A,  ///< (in)    : system matrix
      LINALG::Matrix<1, 1>& b,  ///< (in)    : right-hand-side
      LINALG::Matrix<1, 1>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<true, 2, double>(
      LINALG::Matrix<2, 2>& A,  ///< (in)    : system matrix
      LINALG::Matrix<2, 1>& b,  ///< (in)    : right-hand-side
      LINALG::Matrix<2, 1>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<false, 2, double>(
      LINALG::Matrix<2, 2>& A,  ///< (in)    : system matrix
      LINALG::Matrix<2, 1>& b,  ///< (in)    : right-hand-side
      LINALG::Matrix<2, 1>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<true, 3, double>(
      LINALG::Matrix<3, 3>& A,  ///< (in)    : system matrix
      LINALG::Matrix<3, 1>& b,  ///< (in)    : right-hand-side
      LINALG::Matrix<3, 1>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<false, 3, double>(
      LINALG::Matrix<3, 3>& A,  ///< (in)    : system matrix
      LINALG::Matrix<3, 1>& b,  ///< (in)    : right-hand-side
      LINALG::Matrix<3, 1>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<true, 4, double>(
      LINALG::Matrix<4, 4>& A,  ///< (in)    : system matrix
      LINALG::Matrix<4, 1>& b,  ///< (in)    : right-hand-side
      LINALG::Matrix<4, 1>& x   ///< (out)   : solution vector
  );
  template double gaussElimination<false, 4, double>(
      LINALG::Matrix<4, 4>& A,  ///< (in)    : system matrix
      LINALG::Matrix<4, 1>& b,  ///< (in)    : right-hand-side
      LINALG::Matrix<4, 1>& x   ///< (out)   : solution vector
  );

  template double scaledGaussElimination<2>(LINALG::Matrix<2, 2>& A,  ///< (in)    : system matrix
      LINALG::Matrix<2, 1>& b,                                        ///< (in)    : right-hand-side
      LINALG::Matrix<2, 1>& x                                         ///< (out)   : solution vector
  );
  template double scaledGaussElimination<3>(LINALG::Matrix<3, 3>& A,  ///< (in)    : system matrix
      LINALG::Matrix<3, 1>& b,                                        ///< (in)    : right-hand-side
      LINALG::Matrix<3, 1>& x                                         ///< (out)   : solution vector
  );
  template double scaledGaussElimination<4>(LINALG::Matrix<4, 4>& A,  ///< (in)    : system matrix
      LINALG::Matrix<4, 1>& b,                                        ///< (in)    : right-hand-side
      LINALG::Matrix<4, 1>& x                                         ///< (out)   : solution vector
  );

}  // namespace LINALG

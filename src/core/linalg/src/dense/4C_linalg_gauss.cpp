// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_gauss.templates.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  template double gauss_elimination<true, 1, double>(
      Core::LinAlg::Matrix<1, 1, double>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<1, 1, double>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<1, 1, double>& x   ///< (out)   : solution vector
  );
  template double gauss_elimination<false, 1, double>(
      Core::LinAlg::Matrix<1, 1>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<1, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<1, 1>& x   ///< (out)   : solution vector
  );
  template double gauss_elimination<true, 2, double>(
      Core::LinAlg::Matrix<2, 2>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<2, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<2, 1>& x   ///< (out)   : solution vector
  );
  template double gauss_elimination<false, 2, double>(
      Core::LinAlg::Matrix<2, 2>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<2, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<2, 1>& x   ///< (out)   : solution vector
  );
  template double gauss_elimination<true, 3, double>(
      Core::LinAlg::Matrix<3, 3>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<3, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<3, 1>& x   ///< (out)   : solution vector
  );
  template double gauss_elimination<false, 3, double>(
      Core::LinAlg::Matrix<3, 3>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<3, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<3, 1>& x   ///< (out)   : solution vector
  );
  template double gauss_elimination<true, 4, double>(
      Core::LinAlg::Matrix<4, 4>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<4, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<4, 1>& x   ///< (out)   : solution vector
  );
  template double gauss_elimination<false, 4, double>(
      Core::LinAlg::Matrix<4, 4>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<4, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<4, 1>& x   ///< (out)   : solution vector
  );

  template double scaled_gauss_elimination<2>(
      Core::LinAlg::Matrix<2, 2>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<2, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<2, 1>& x   ///< (out)   : solution vector
  );
  template double scaled_gauss_elimination<3>(
      Core::LinAlg::Matrix<3, 3>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<3, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<3, 1>& x   ///< (out)   : solution vector
  );
  template double scaled_gauss_elimination<4>(
      Core::LinAlg::Matrix<4, 4>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<4, 1>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<4, 1>& x   ///< (out)   : solution vector
  );

}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

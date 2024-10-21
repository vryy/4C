// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_serialdensevector.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Compute vector 2-norm                                               |
 *----------------------------------------------------------------------*/
double Core::LinAlg::norm2(const Core::LinAlg::SerialDenseVector& v) { return v.normFrobenius(); }

/*----------------------------------------------------------------------*
 |  b = alpha*a + beta*b                                                |
 *----------------------------------------------------------------------*/
void Core::LinAlg::update(double alpha, const Core::LinAlg::SerialDenseVector& a, double beta,
    Core::LinAlg::SerialDenseVector& b)
{
  b.scale(beta);
  Core::LinAlg::SerialDenseVector acopy(a);
  acopy.scale(alpha);
  b += acopy;
}

FOUR_C_NAMESPACE_CLOSE

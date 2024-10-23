// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_fixedsizematrix_tensor_derivatives.hpp"



FOUR_C_NAMESPACE_OPEN

void Core::LinAlg::Tensor::add_derivative_of_squared_tensor(Core::LinAlg::Matrix<6, 6>& C,
    double scalar_squared_dx, Core::LinAlg::Matrix<3, 3> X, double scalar_this)
{
  C(0, 0) = scalar_this * C(0, 0) + scalar_squared_dx * 2. * X(0, 0);  // C1111
  C(0, 1) = scalar_this * C(0, 1);                                     // C1122
  C(0, 2) = scalar_this * C(0, 2);                                     // C1133
  C(0, 3) = scalar_this * C(0, 3) + scalar_squared_dx * X(0, 1);       // C1112
  C(0, 4) = scalar_this * C(0, 4);                                     // C1123
  C(0, 5) = scalar_this * C(0, 5) + scalar_squared_dx * X(0, 2);       // C1113

  C(1, 0) = scalar_this * C(1, 0);                                     // C2211
  C(1, 1) = scalar_this * C(1, 1) + scalar_squared_dx * 2. * X(1, 1);  // C2222
  C(1, 2) = scalar_this * C(1, 2);                                     // C2233
  C(1, 3) = scalar_this * C(1, 3) + scalar_squared_dx * X(0, 1);       // C2212
  C(1, 4) = scalar_this * C(1, 4) + scalar_squared_dx * X(1, 2);       // C2223
  C(1, 5) = scalar_this * C(1, 5);                                     // C2213

  C(2, 0) = scalar_this * C(2, 0);                                     // C3311
  C(2, 1) = scalar_this * C(2, 1);                                     // C3322
  C(2, 2) = scalar_this * C(2, 2) + scalar_squared_dx * 2. * X(2, 2);  // C3333
  C(2, 3) = scalar_this * C(2, 3);                                     // C3312
  C(2, 4) = scalar_this * C(2, 4) + scalar_squared_dx * X(1, 2);       // C3323
  C(2, 5) = scalar_this * C(2, 5) + scalar_squared_dx * X(0, 2);       // C3313

  C(3, 0) = scalar_this * C(3, 0) + scalar_squared_dx * X(0, 1);                    // C1211
  C(3, 1) = scalar_this * C(3, 1) + scalar_squared_dx * X(0, 1);                    // C1222
  C(3, 2) = scalar_this * C(3, 2);                                                  // C1233
  C(3, 3) = scalar_this * C(3, 3) + scalar_squared_dx * 0.5 * (X(0, 0) + X(1, 1));  // C1212
  C(3, 4) = scalar_this * C(3, 4) + scalar_squared_dx * 0.5 * X(0, 2);              // C1223
  C(3, 5) = scalar_this * C(3, 5) + scalar_squared_dx * 0.5 * X(1, 2);              // C1213

  C(4, 0) = scalar_this * C(4, 0);                                                  // C2311
  C(4, 1) = scalar_this * C(4, 1) + scalar_squared_dx * X(1, 2);                    // C2322
  C(4, 2) = scalar_this * C(4, 2) + scalar_squared_dx * X(1, 2);                    // C2333
  C(4, 3) = scalar_this * C(4, 3) + scalar_squared_dx * 0.5 * X(0, 2);              // C2312
  C(4, 4) = scalar_this * C(4, 4) + scalar_squared_dx * 0.5 * (X(1, 1) + X(2, 2));  // C2323
  C(4, 5) = scalar_this * C(4, 5) + scalar_squared_dx * 0.5 * X(0, 1);              // C2313

  C(5, 0) = scalar_this * C(5, 0) + scalar_squared_dx * X(0, 2);                    // C1311
  C(5, 1) = scalar_this * C(5, 1);                                                  // C1322
  C(5, 2) = scalar_this * C(5, 2) + scalar_squared_dx * X(0, 2);                    // C1333
  C(5, 3) = scalar_this * C(5, 3) + scalar_squared_dx * 0.5 * X(1, 2);              // C1312
  C(5, 4) = scalar_this * C(5, 4) + scalar_squared_dx * 0.5 * X(0, 1);              // C1323
  C(5, 5) = scalar_this * C(5, 5) + scalar_squared_dx * 0.5 * (X(2, 2) + X(0, 0));  // C1313
}

void Core::LinAlg::Tensor::add_derivative_of_inva_b_inva_product(double const& fac,
    const Core::LinAlg::Matrix<6, 1>& invA, const Core::LinAlg::Matrix<6, 1>& invABinvA,
    Core::LinAlg::Matrix<6, 6>& out)
{
  out(0, 0) -= 2.0 * fac * (invA(0) * invABinvA(0));
  out(0, 1) -= 2.0 * fac * (invA(3) * invABinvA(3));
  out(0, 2) -= 2.0 * fac * (invA(5) * invABinvA(5));
  out(0, 3) -= fac * (invA(0) * invABinvA(3) + invA(3) * invABinvA(0));
  out(0, 4) -= fac * (invA(3) * invABinvA(5) + invA(5) * invABinvA(3));
  out(0, 5) -= fac * (invA(0) * invABinvA(5) + invA(5) * invABinvA(0));

  out(1, 0) -= 2.0 * fac * (invA(3) * invABinvA(3));
  out(1, 1) -= 2.0 * fac * (invA(1) * invABinvA(1));
  out(1, 2) -= 2.0 * fac * (invA(4) * invABinvA(4));
  out(1, 3) -= fac * (invA(3) * invABinvA(1) + invA(1) * invABinvA(3));
  out(1, 4) -= fac * (invA(1) * invABinvA(4) + invA(4) * invABinvA(1));
  out(1, 5) -= fac * (invA(3) * invABinvA(4) + invA(4) * invABinvA(3));

  out(2, 0) -= 2.0 * fac * (invA(5) * invABinvA(5));
  out(2, 1) -= 2.0 * fac * (invA(4) * invABinvA(4));
  out(2, 2) -= 2.0 * fac * (invA(2) * invABinvA(2));
  out(2, 3) -= fac * (invA(5) * invABinvA(4) + invA(4) * invABinvA(5));
  out(2, 4) -= fac * (invA(4) * invABinvA(2) + invA(2) * invABinvA(4));
  out(2, 5) -= fac * (invA(5) * invABinvA(2) + invA(2) * invABinvA(5));

  out(3, 0) -= fac * (invA(0) * invABinvA(3) + invA(3) * invABinvA(0));
  out(3, 1) -= fac * (invA(3) * invABinvA(1) + invA(1) * invABinvA(3));
  out(3, 2) -= fac * (invA(5) * invABinvA(4) + invA(4) * invABinvA(5));
  out(3, 3) -=
      fac * (0.5 * (invA(0) * invABinvA(1) + invA(1) * invABinvA(0)) + invA(3) * invABinvA(3));
  out(3, 4) -= 0.5 * fac *
               (invA(3) * invABinvA(4) + invA(5) * invABinvA(1) + invA(1) * invABinvA(5) +
                   invA(4) * invABinvA(3));
  out(3, 5) -= 0.5 * fac *
               (invA(0) * invABinvA(4) + invA(5) * invABinvA(3) + invA(3) * invABinvA(5) +
                   invA(4) * invABinvA(0));

  out(4, 0) -= fac * (invA(3) * invABinvA(5) + invA(5) * invABinvA(3));
  out(4, 1) -= fac * (invA(1) * invABinvA(4) + invA(4) * invABinvA(1));
  out(4, 2) -= fac * (invA(4) * invABinvA(2) + invA(2) * invABinvA(4));
  out(4, 3) -= 0.5 * fac *
               (invA(3) * invABinvA(4) + invA(1) * invABinvA(5) + invA(5) * invABinvA(1) +
                   invA(4) * invABinvA(3));
  out(4, 4) -=
      fac * (0.5 * (invA(1) * invABinvA(2) + invA(2) * invABinvA(1)) + invA(4) * invABinvA(4));
  out(4, 5) -= 0.5 * fac *
               (invA(3) * invABinvA(2) + invA(4) * invABinvA(5) + invA(5) * invABinvA(4) +
                   invA(2) * invABinvA(3));

  out(5, 0) -= fac * (invA(0) * invABinvA(5) + invA(5) * invABinvA(0));
  out(5, 1) -= fac * (invA(3) * invABinvA(4) + invA(4) * invABinvA(3));
  out(5, 2) -= fac * (invA(5) * invABinvA(2) + invA(2) * invABinvA(5));
  out(5, 3) -= 0.5 * fac *
               (invA(0) * invABinvA(4) + invA(3) * invABinvA(5) + invA(5) * invABinvA(3) +
                   invA(4) * invABinvA(0));
  out(5, 4) -= 0.5 * fac *
               (invA(3) * invABinvA(2) + invA(5) * invABinvA(4) + invA(4) * invABinvA(5) +
                   invA(2) * invABinvA(3));
  out(5, 5) -=
      fac * (0.5 * (invA(0) * invABinvA(2) + invA(2) * invABinvA(0)) + invA(5) * invABinvA(5));
}


FOUR_C_NAMESPACE_CLOSE

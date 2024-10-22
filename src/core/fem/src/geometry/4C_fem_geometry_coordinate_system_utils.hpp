// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_GEOMETRY_COORDINATE_SYSTEM_UTILS_HPP
#define FOUR_C_FEM_GEOMETRY_COORDINATE_SYSTEM_UTILS_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  /*!
   * create a right-handed orthornormal basis from a unit vector based on
   * https://doi.org/10.1080/10867651.1999.10487513
   * @param unitvec[i]    : unit vector based on which basis will be created (third basis vector)
   * @param basisvec_1[o] : first basis vector
   * @param basisvec_2[o] : second basis vector
   */
  void build_orthonormal_basis_from_unit_vector(const Core::LinAlg::Matrix<3, 1>& unitvec,
      Core::LinAlg::Matrix<3, 1>& basisvec_1, Core::LinAlg::Matrix<3, 1>& basisvec_2)
  {
    if (std::abs(unitvec.norm2() - 1.0) > 1.0e-14) FOUR_C_THROW("given vector not normalized!");

    if ((std::abs(unitvec(0)) <= std::abs(unitvec(1))) and
        (std::abs(unitvec(0)) <= std::abs(unitvec(2))))
    {
      basisvec_1(0) = 0.0;
      basisvec_1(1) = -unitvec(2);
      basisvec_1(2) = unitvec(1);
    }
    else if (std::abs(unitvec(1)) <= std::abs(unitvec(2)))
    {
      basisvec_1(0) = -unitvec(2);
      basisvec_1(1) = 0.0;
      basisvec_1(2) = unitvec(0);
    }
    else
    {
      basisvec_1(0) = -unitvec(1);
      basisvec_1(1) = unitvec(0);
      basisvec_1(2) = 0.0;
    }

    basisvec_1.scale(1.0 / basisvec_1.norm2());

    basisvec_2.cross_product(unitvec, basisvec_1);
  }
}  // namespace Core::Geo


FOUR_C_NAMESPACE_CLOSE

#endif

// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_MONOLITHIC_SOLID_SCALAR_MATERIAL_HPP
#define FOUR_C_MAT_MONOLITHIC_SOLID_SCALAR_MATERIAL_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  class MonolithicSolidScalarMaterial
  {
   public:
    virtual ~MonolithicSolidScalarMaterial() = default;

    /*!
     * @brief Evaluates the added derivatives of the stress w.r.t. all scalars
     *
     * @param defgrad (in) : Deformation gradient
     * @param glstrain (in) : Green-Lagrange strain
     * @param params (in) : ParameterList for additional parameters
     * @param gp (in) : Gauss points
     * @param eleGID (in) : global element id
     * @return std::vector<std::optional<Core::LinAlg::Matrix<6, 1>>>
     */
    virtual Core::LinAlg::SymmetricTensor<double, 3, 3> evaluate_d_stress_d_scalar(
        const Core::LinAlg::Tensor<double, 3, 3>& defgrad,
        const Core::LinAlg::SymmetricTensor<double, 3, 3>& glstrain,
        const Teuchos::ParameterList& params, const EvaluationContext<3>& context, int gp,
        int eleGID) = 0;
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
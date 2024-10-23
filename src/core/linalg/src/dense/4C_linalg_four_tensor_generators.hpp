// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_FOUR_TENSOR_GENERATORS_HPP
#define FOUR_C_LINALG_FOUR_TENSOR_GENERATORS_HPP


#include "4C_config.hpp"

#include "4C_linalg_four_tensor.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  /*!
   * @brief Compute the fourth order deviatoric tensor.
   *
   *  @note This tensor projects the stress tensor to its deviatoric variant via contraction
   * operator, i.e., \f$\mathbb{I}_{d} : \boldsymbol{sigma} = \mathbf{s} (\boldsymbol{\sigma})\f$.
   *
   * This tensor can be computed as
   *  \f$\mathbb{I}_{d} = \dfrac{1}{2} \delta_{ik} \delta_{jl} + \dfrac{1}{2} \delta_{il}
   * \delta_{jk} - \dfrac{1}{3} \delta_{ij} \delta_{kl}\f$
   *
   * Reference:
   *    Souze de Neto, Computational Methods for Plasticity, 2008, page 10.
   *
   * @param[out]  four_tensor   the resulting deviatoric tensor
   * @param[in]   scale         the scaling factor
   */
  template <unsigned int dim>
  Core::LinAlg::FourTensor<dim> setup_deviatoric_projection_tensor(const double scale = 1.0)
  {
    // declare output variable
    Core::LinAlg::FourTensor<dim> deviatoric_projection_tensor;

    const auto eye = [](int i, int j) { return i == j ? 1.0 : 0.0; };

    for (unsigned int i = 0; i < dim; ++i)
    {
      for (unsigned int j = 0; j < dim; ++j)
      {
        for (unsigned int k = 0; k < dim; ++k)
        {
          for (unsigned int l = 0; l < dim; ++l)
            deviatoric_projection_tensor(i, j, k, l) =
                scale * (0.5 * eye(i, k) * eye(j, l) + 0.5 * eye(i, l) * eye(j, k) -
                            1.0 / 3 * eye(i, j) * eye(k, l));
        }
      }
    }

    return deviatoric_projection_tensor;
  }
}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif

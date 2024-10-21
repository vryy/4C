// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MATELAST_ACTIVESUMMAND_HPP
#define FOUR_C_MATELAST_ACTIVESUMMAND_HPP
#include "4C_config.hpp"

#include "4C_matelast_summand.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  namespace Elastic
  {
    /*!
     * \brief This is a pure abstract extension of the Summand class to be used for active
     * materials.
     */
    class ActiveSummand : public Summand
    {
     public:
      /*!
       * \brief adds the active part of the summand to the stress and stiffness matrix
       *
       * \note This is ONLY the ACTIVE PART of the stress response
       *
       * \param CM Right Cauchy-Green deformation tensor in Tensor notation
       * \param cmat Material stiffness matrix
       * \param stress 2nd PK-stress
       * \param eleGID global element id
       */
      virtual void add_active_stress_cmat_aniso(
          const Core::LinAlg::Matrix<3, 3>& CM,  ///< right Cauchy Green tensor
          Core::LinAlg::Matrix<6, 6>& cmat,      ///< material stiffness matrix
          Core::LinAlg::Matrix<6, 1>& stress,    ///< 2nd PK-stress
          int gp,                                ///< Gauss point
          int eleGID) const = 0;                 ///< element GID

      /*!
       * \brief Evaluates the first derivative of active fiber potential w.r.t. the active fiber
       * stretch
       *
       * \return double First derivative of active fiber potential w.r.t. the active fiber stretch
       */
      virtual double get_derivative_aniso_active() const = 0;
    };
  }  // namespace Elastic

}  // namespace Mat


FOUR_C_NAMESPACE_CLOSE

#endif
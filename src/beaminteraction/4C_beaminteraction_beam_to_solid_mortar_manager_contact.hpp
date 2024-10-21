// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_MORTAR_MANAGER_CONTACT_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_MORTAR_MANAGER_CONTACT_HPP

#include "4C_config.hpp"

#include "4C_beaminteraction_beam_to_solid_mortar_manager.hpp"

FOUR_C_NAMESPACE_OPEN

// Forward declaration
namespace BEAMINTERACTION
{
  class BeamToSolidParamsBase;
}  // namespace BEAMINTERACTION


namespace BEAMINTERACTION
{
  /**
   * \brief This is a specialization of the mesh tying mortar manager for contact
   */
  class BeamToSolidMortarManagerContact : public BeamToSolidMortarManager
  {
   public:
    /**
     * \brief Standard Constructor
     *
     * @params discret (in) Pointer to the discretization.
     * @params params (in) Beam-to-solid parameters.
     * @params start_value_lambda_gid (in) Start value for the Lagrange multiplier global IDs.
     */
    BeamToSolidMortarManagerContact(const Teuchos::RCP<const Core::FE::Discretization>& discret,
        const Teuchos::RCP<const BEAMINTERACTION::BeamToSolidParamsBase>& params,
        int start_value_lambda_gid);

   protected:
    /**
     * \brief Get the penalty regularization of the constraint vector (derived)
     */
    [[nodiscard]] std::tuple<Teuchos::RCP<Core::LinAlg::Vector<double>>,
        Teuchos::RCP<Core::LinAlg::Vector<double>>, Teuchos::RCP<Core::LinAlg::Vector<double>>>
    get_penalty_regularization(const bool compute_linearization = false) const override;
  };
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif

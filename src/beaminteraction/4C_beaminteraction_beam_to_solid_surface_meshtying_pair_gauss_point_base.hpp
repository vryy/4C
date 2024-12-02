// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_GAUSS_POINT_BASE_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_GAUSS_POINT_BASE_HPP


#include "4C_config.hpp"

#include "4C_beaminteraction_beam_to_solid_surface_meshtying_pair_base.hpp"

FOUR_C_NAMESPACE_OPEN


namespace BeamInteraction
{
  /**
   * \brief Class for Gauss-point-to-segment beam to surface surface mesh tying.
   * @tparam scalar_type Scalar type for variables.
   * @tparam beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @tparam surface Type from GEOMETRYPAIR::ElementDiscretization... representing the surface.
   */
  template <typename ScalarType, typename Beam, typename Surface>
  class BeamToSolidSurfaceMeshtyingPairGaussPointBase
      : public BeamToSolidSurfaceMeshtyingPairBase<ScalarType, Beam, Surface>
  {
   private:
    //! Shortcut to the base class.
    using base_class = BeamToSolidSurfaceMeshtyingPairBase<ScalarType, Beam, Surface>;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidSurfaceMeshtyingPairGaussPointBase();


    /**
     * \brief Get energy for Gauss point to segment coupling.
     */
    double get_energy() const override;

   protected:
    /**
     * \brief Get the penalty potential for Gauss point to segment coupling.
     * @return Penalty potential.
     */
    ScalarType get_penalty_potential() const;
  };
}  // namespace BeamInteraction

FOUR_C_NAMESPACE_CLOSE

#endif

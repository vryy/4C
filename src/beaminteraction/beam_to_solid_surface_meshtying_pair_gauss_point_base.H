/*----------------------------------------------------------------------*/
/*! \file

\brief Base Gauss point to segment mesh tying element for between a 3D beam and a surface element.

\level 3
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_GAUSS_POINT_BASE_H_
#define BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_GAUSS_POINT_BASE_H_


#include "beam_to_solid_surface_meshtying_pair_base.H"


namespace BEAMINTERACTION
{
  /**
   * \brief Class for Gauss-point-to-segment beam to surface surface mesh tying.
   * @tparam scalar_type Scalar type for variables.
   * @tparam beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @tparam surface Type from GEOMETRYPAIR::ElementDiscretization... representing the surface.
   */
  template <typename scalar_type, typename beam, typename surface>
  class BeamToSolidSurfaceMeshtyingPairGaussPointBase
      : public BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam, surface>
  {
   private:
    //! Shortcut to the base class.
    using base_class = BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam, surface>;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidSurfaceMeshtyingPairGaussPointBase();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidSurfaceMeshtyingPairGaussPointBase(){};

    /**
     * \brief Get energy for Gauss point to segment coupling.
     */
    double GetEnergy() const override;

   protected:
    /**
     * \brief Get the penalty potential for Gauss point to segment coupling.
     * @return Penalty potential.
     */
    scalar_type GetPenaltyPotential() const;
  };
}  // namespace BEAMINTERACTION

#endif

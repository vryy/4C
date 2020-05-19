/*----------------------------------------------------------------------*/
/*! \file

\brief Mortar mesh tying element for between a 3D beam and a surface element.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_MORTAR_H_
#define BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_MORTAR_H_


#include "beam_to_solid_surface_meshtying_pair_base.H"

#include "../drt_geometry_pair/geometry_pair_scalar_types.H"


namespace BEAMINTERACTION
{
  /**
   * \brief Class for Mortar beam to surface surface mesh tying.
   * @tparam beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @tparam surface Type from GEOMETRYPAIR::ElementDiscretization... representing the surface.
   * @tparam mortar Type from BEAMINTERACTION::ElementDiscretization... representing the mortar
   * shape functions.
   */
  template <typename beam, typename surface, typename mortar>
  class BeamToSolidSurfaceMeshtyingPairMortar
      : public BeamToSolidSurfaceMeshtyingPairBase<
            GEOMETRYPAIR::line_to_surface_scalar_type<beam, surface>, beam, surface>
  {
   private:
    //! Type to be used for scalar AD variables.
    using scalar_type = GEOMETRYPAIR::line_to_surface_scalar_type<beam, surface>;

    //! Shortcut to the base class.
    using base_class = BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam, surface>;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidSurfaceMeshtyingPairMortar();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidSurfaceMeshtyingPairMortar(){};

    /**
     * \brief This pair has no terms that need to be directly assembled into the global matrices.
     */
    inline bool IsAssemblyDirect() const override { return false; };

    /**
     * \brief This pair has terms that need to be indirectly (mortar like) assembled into the global
     * matrices.
     */
    inline bool IsAssemblyIndirect() const override { return true; };
  };
}  // namespace BEAMINTERACTION

#endif

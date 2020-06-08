/*----------------------------------------------------------------------*/
/*! \file

\brief Base mortar mesh tying element for between a 3D beam and a surface element.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_MORTAR_BASE_H_
#define BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_MORTAR_BASE_H_


#include "beam_to_solid_surface_meshtying_pair_base.H"

#include "../drt_geometry_pair/geometry_pair_scalar_types.H"


namespace BEAMINTERACTION
{
  /**
   * \brief Base class for Mortar beam to surface surface mesh tying.
   * @tparam scalar_type Type for scalar variables.
   * @tparam beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @tparam surface Type from GEOMETRYPAIR::ElementDiscretization... representing the surface.
   * @tparam mortar Type from BEAMINTERACTION::ElementDiscretization... representing the mortar
   * shape functions.
   */
  template <typename scalar_type, typename beam, typename surface, typename mortar>
  class BeamToSolidSurfaceMeshtyingPairMortarBase
      : public BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam, surface>
  {
   private:
    //! Shortcut to the base class.
    using base_class = BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam, surface>;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidSurfaceMeshtyingPairMortarBase();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidSurfaceMeshtyingPairMortarBase(){};

    /**
     * \brief This pair has no terms that need to be directly assembled into the global matrices.
     */
    inline bool IsAssemblyDirect() const override { return false; };

    /**
     * \brief This pair has terms that need to be indirectly (mortar like) assembled into the global
     * matrices.
     */
    inline bool IsAssemblyIndirect() const override { return true; };

    /**
     * \brief Add the visualization of this pair to the beam to solid vtu output writer. This will
     * add mortar specific data to the output.
     * @param visualization_writer (out) Object that manages all visualization related data for beam
     * to solid pairs.
     * @param visualization_params (in) Parameter list (not used in this class).
     */
    void GetPairVisualization(Teuchos::RCP<BeamToSolidVtuOutputWriterBase> visualization_writer,
        Teuchos::ParameterList& visualization_params) const override;

    /**
     * \brief The mortar energy contribution will be calculated globally in the mortar manager.
     */
    double GetEnergy() const override { return 0.0; }
  };
}  // namespace BEAMINTERACTION

#endif

/*----------------------------------------------------------------------*/
/*! \file

\brief Base mortar mesh tying element for between a 3D beam and a surface element.

\level 3
*/
// End doxygen header.


#ifndef BACI_BEAMINTERACTION_BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_MORTAR_BASE_H
#define BACI_BEAMINTERACTION_BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_MORTAR_BASE_H


#include "baci_beaminteraction_beam_to_solid_surface_meshtying_pair_base.H"
#include "baci_geometry_pair_scalar_types.H"


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
     * \brief This pair enforces constraints via a mortar-type method, which requires an own
     * assembly method (provided by the mortar manager).
     */
    inline bool IsAssemblyDirect() const override { return false; };

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

   protected:
    //! Number of rotational Lagrange multipliers.
    unsigned int n_mortar_rot_;
  };
}  // namespace BEAMINTERACTION

#endif

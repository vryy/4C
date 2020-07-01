/*----------------------------------------------------------------------*/
/*! \file

\brief Mortar mesh tying element for between a 3D beam and a surface element.

\level 3
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_MORTAR_H_
#define BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_MORTAR_H_


#include "beam_to_solid_surface_meshtying_pair_mortar_base.H"

#include "../drt_geometry_pair/geometry_pair_scalar_types.H"


// Forward declaration.
namespace INPAR
{
  namespace BEAMTOSOLID
  {
    enum class BeamToSolidMortarShapefunctions;
  }
}  // namespace INPAR


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
      : public BeamToSolidSurfaceMeshtyingPairMortarBase<
            GEOMETRYPAIR::line_to_surface_scalar_type<beam, surface>, beam, surface, mortar>
  {
   private:
    //! Type to be used for scalar AD variables.
    using scalar_type = GEOMETRYPAIR::line_to_surface_scalar_type<beam, surface>;

    //! Shortcut to the base class.
    using base_class =
        BeamToSolidSurfaceMeshtyingPairMortarBase<scalar_type, beam, surface, mortar>;

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
     * \brief Evaluate the mortar matrices $D$ and $M$ for this contact element pair. (derived)
     */
    bool EvaluateDM(LINALG::SerialDenseMatrix& local_D, LINALG::SerialDenseMatrix& local_M,
        LINALG::SerialDenseVector& local_kappa,
        LINALG::SerialDenseVector& local_constraint) override;
  };

  /**
   * \brief Factory function for beam-to-solid mortar pairs.
   * @param surface_shape (in) Type of surface element.
   * @param mortar_shapefunction (in) Type of mortar shape function.
   * @return Pointer to the created pair.
   */
  Teuchos::RCP<BEAMINTERACTION::BeamContactPair> BeamToSolidSurfaceMeshtyingPairMortarFactory(
      const DRT::Element::DiscretizationType surface_shape,
      const INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions mortar_shapefunction);
}  // namespace BEAMINTERACTION

#endif

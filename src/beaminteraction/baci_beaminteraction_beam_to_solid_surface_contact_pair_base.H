/*----------------------------------------------------------------------*/
/*! \file

\brief Base contact element for contact between a 3D beam and a surface element.

\level 3
*/
// End doxygen header.


#ifndef BACI_BEAMINTERACTION_BEAM_TO_SOLID_SURFACE_CONTACT_PAIR_BASE_H
#define BACI_BEAMINTERACTION_BEAM_TO_SOLID_SURFACE_CONTACT_PAIR_BASE_H


#include "baci_beaminteraction_beam3contact_defines.H"
#include "baci_beaminteraction_beam_to_solid_pair_base.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_linalg_sparsematrix.H"


// Forward declarations.
namespace CORE::LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace CORE::LINALG
namespace GEOMETRYPAIR
{
  template <typename scalar_type, typename line, typename surface>
  class GeometryPairLineToSurface;

  class FaceElement;

  template <typename surface, typename scalar_type>
  class FaceElementTemplate;
}  // namespace GEOMETRYPAIR
namespace BEAMINTERACTION
{
  class BeamToSolidVtuOutputWriterVisualization;
}  // namespace BEAMINTERACTION


namespace BEAMINTERACTION
{
  /**
   * \brief Base class for beam to surface surface contact.
   * @tparam scalar_type Type for scalar DOF values.
   * @tparam beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @tparam surface Type from GEOMETRYPAIR::ElementDiscretization... representing the surface.
   */
  template <typename scalar_type, typename beam, typename surface>
  class BeamToSolidSurfaceContactPairBase
      : public BeamToSolidPairBase<scalar_type, scalar_type, beam, surface>
  {
   protected:
    //! Shortcut to the base class.
    using base_class = BeamToSolidPairBase<scalar_type, scalar_type, beam, surface>;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidSurfaceContactPairBase();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidSurfaceContactPairBase(){};

    /**
     * \brief Update state of translational nodal DoFs (absolute positions and tangents) of the beam
     * element. (derived)
     *
     * This function has to be overwritten here, since the size of FAD variables for surface
     * elements is not known at compile time and has to be set depending on the surface patch that
     * the surface element is part of.
     *
     * @param beam_centerline_dofvec
     * @param solid_nodal_dofvec
     */
    void ResetState(const std::vector<double>& beam_centerline_dofvec,
        const std::vector<double>& solid_nodal_dofvec) override;

    /**
     * \brief Things that need to be done in a separate loop before the actual evaluation loop over
     * the contact pairs.
     */
    void PreEvaluate() override;

    /**
     * \brief Create the geometry pair for this contact pair.
     * @param geometry_evaluation_data_ptr Evaluation data that will be linked to the pair.
     */
    void CreateGeometryPair(
        const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>& geometry_evaluation_data_ptr)
        override;

    /**
     * \brief Link the contact pair with the face element storing information on the averaged nodal
     * normals (derived).
     */
    void SetFaceElement(Teuchos::RCP<GEOMETRYPAIR::FaceElement>& face_element) override;

   protected:
    /**
     * \brief Return a cast of the geometry pair to the type for this contact pair.
     * @return RPC with the type of geometry pair for this beam contact pair.
     */
    Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, beam, surface>>
    CastGeometryPair() const;

    /**
     * \brief Get the GIDs of the pair, i.e. first the beam GIDs and then the pair GIDs.
     * @param discret (in) Discretization.
     * @return Vector with the GIDs of this pair.
     */
    std::vector<int> GetPairGID(const ::DRT::Discretization& discret) const;

   protected:
    //! Pointer to the face element object which manages the positions on the surface, including the
    //! averaged nodal normals.
    Teuchos::RCP<GEOMETRYPAIR::FaceElementTemplate<surface, scalar_type>> face_element_;
  };
}  // namespace BEAMINTERACTION

#endif

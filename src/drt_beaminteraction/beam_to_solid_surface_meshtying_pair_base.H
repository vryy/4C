/*----------------------------------------------------------------------*/
/*! \file

\brief Base meshtying element for meshtying between a 3D beam and a surface element.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_BASE_H_
#define BEAM_TO_SOLID_SURFACE_MESHTYING_PAIR_BASE_H_


#include "beam_to_solid_pair_base.H"

#include "beam3contact_defines.H"

#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_sparsematrix.H"


// Forward declarations.
namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace LINALG
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
   * \brief Class for beam to surface surface mesh tying.
   * @tparam scalar_type Type for scalar DOF values.
   * @tparam beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @tparam surface Type from GEOMETRYPAIR::ElementDiscretization... representing the surface.
   */
  template <typename scalar_type, typename beam, typename surface>
  class BeamToSolidSurfaceMeshtyingPairBase : public BeamToSolidPairBase<scalar_type, beam, surface>
  {
   protected:
    //! Shortcut to the base class.
    using base_class = BeamToSolidPairBase<scalar_type, beam, surface>;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidSurfaceMeshtyingPairBase();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidSurfaceMeshtyingPairBase(){};

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
     * \brief Add the visualization of this pair to the beam to solid vtu output writer.
     *
     * Create segmentation and integration points output.
     *
     * @param visualization_writer (out) Object that manages all visualization related data for beam
     * to solid pairs.
     * @param visualization_params (in) Parameter list (not used in this class).
     */
    void GetPairVisualization(Teuchos::RCP<BeamToSolidVtuOutputWriterBase> visualization_writer,
        Teuchos::ParameterList& visualization_params) const override;

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
    Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToSurface<double, beam, surface>> CastGeometryPair()
        const;

    /**
     * \brief Evaluate the coupling vector at an evaluation point.
     * @param evaluation_point (in) Projection point from 1D to 3D where the coupling should be
     * evaluated.
     * @return 3D vector with the difference in the positions / displacements at the projection
     * point.
     */
    LINALG::Matrix<3, 1, scalar_type> EvaluateCoupling(
        const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& evaluation_point) const;

    /**
     * \brief Get the GIDs of the pair, i.e. first the beam GIDs and then the pair GIDs.
     * @param discret (in) Discretization.
     * @return Vector with the GIDs of this pair.
     */
    std::vector<int> GetPairGID(const DRT::Discretization& discret) const;

   private:
    /**
     * \brief Add points on the beam element to an output writer.
     * @param visualization_writer (in/out) Output writer the points are appended to.
     * @param points (in) Vector with the projection points.
     */
    void AddVisualizationIntegrationPoints(
        const Teuchos::RCP<BeamToSolidVtuOutputWriterVisualization>& visualization_writer,
        const std::vector<GEOMETRYPAIR::ProjectionPoint1DTo3D<double>>& points) const;

   protected:
    //! Flag if the meshtying has been evaluated already.
    bool meshtying_is_evaluated_;

    //! Pointer to the face element object which manages the positions on the surface, including the
    //! averaged nodal normals.
    Teuchos::RCP<GEOMETRYPAIR::FaceElementTemplate<surface, scalar_type>> face_element_;
  };
}  // namespace BEAMINTERACTION

#endif

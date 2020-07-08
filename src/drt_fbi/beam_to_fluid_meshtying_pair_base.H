/*----------------------------------------------------------------------*/
/*! \file

\brief Base meshtying element for meshtying between a 1D beam and a 3D fluid element.

\level 2
*/

#ifndef BEAM_TO_FLUID_VOLUME_MESHTYING_PAIR_BASE_H_
#define BEAM_TO_FLUID_VOLUME_MESHTYING_PAIR_BASE_H_

#include "../drt_beaminteraction/beam_to_solid_volume_meshtying_pair_base.H"

// Forward declarations.
namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace LINALG

namespace GEOMETRYPAIR
{
  template <typename scalar_type, typename line, typename volume>
  class GeometryPairLineToVolume;
  class LineTo3DEvaluationData;
  class GeometryEvaluationDataBase;
}  // namespace GEOMETRYPAIR


namespace BEAMINTERACTION
{
  /**
   * \brief Class representing a pair of elements for beam to fluid meshtying
   *
   * \param[in] beam Type from GEOMETRYPAIR::ElementDiscretization representing the beam.
   * \param[in] fluid Type from GEOMETRYPAIR::ElementDiscretization representing the fluid.
   */
  template <typename beam, typename fluid>
  class BeamToFluidMeshtyingPairBase
      : public BeamToSolidVolumeMeshtyingPairBase<beam,
            fluid>  // todo can i just derice from beamtosolidmeshtyingpairbase?
  {
   protected:
    //! Shortcut to base class.
    using base_class = BeamToSolidVolumeMeshtyingPairBase<beam, fluid>;

    //! Scalar type for FAD variables.
    using scalar_type = typename base_class::scalar_type;

   public:
    /**
     * \brief empty Destructor.
     */
    virtual ~BeamToFluidMeshtyingPairBase(){};

    /**
     * \brief Setup the contact pair and set information on the current position of the elements in
     * the pair
     */
    void Setup() override;

    /**
     * \brief Things that need to be done in a separate loop before the actual evaluation loop over
     * all contact pairs.
     */
    void PreEvaluate() override;

    /**
     * \brief Update state of translational nodal DoFs (absolute positions, tangents and velocities)
     * of both elements.
     * @param beam_centerline_dofvec current nodal beam positions extracted from the element and
     * nodal velocities computed by the time integrator
     * @param fluid_nodal_dofvec current nodal fluid positions (only for ALE different from the
     * reference nodal values) and nodal velocities
     */
    void ResetState(const std::vector<double>& beam_centerline_dofvec,
        const std::vector<double>& fluid_nodal_dofvec) override;

    /**
     * \brief Print information about this beam contact element pair to screen.
     */
    void Print(std::ostream& out) const override;

    /**
     * \brief Print this beam contact element pair to screen.
     */
    void PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const override;

    /**
     * \brief Add the visualization of this pair to the beam to solid vtu output writer.
     *
     * This base class creates output of (if selected in the input file) the segmentation, the
     * integration points - and if implemented in the derived classes - the forces at the
     * integration points.
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

   protected:
    /** \brief You will have to use the FBI::PairFactory
     *
     */
    BeamToFluidMeshtyingPairBase();

    void EvaluateBeamPosition(const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& integration_point,
        LINALG::Matrix<3, 1, scalar_type>& r_beam, bool reference) const;

    //! Current nodal velocities of the two elements.
    LINALG::Matrix<beam::n_dof_, 1, scalar_type> ele1vel_;
    LINALG::Matrix<fluid::n_dof_, 1, scalar_type> ele2vel_;

    //! Current nodal positions (and tangents) of the two elements.
    LINALG::Matrix<beam::n_dof_, 1, double> ele1poscur_;
    LINALG::Matrix<fluid::n_dof_, 1, double> ele2poscur_;
  };
}  // namespace BEAMINTERACTION

#endif

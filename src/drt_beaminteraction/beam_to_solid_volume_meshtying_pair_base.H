/*----------------------------------------------------------------------*/
/*! \file

\brief Base meshtying element for meshtying between a 3D beam and a 3D solid element.

\level 3
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_BASE_H_
#define BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_BASE_H_


#include "beam_to_solid_pair_base.H"

#include "../drt_geometry_pair/geometry_pair_scalar_types.H"
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
  template <typename scalar_type, typename line, typename volume>
  class GeometryPairLineToVolume;
}  // namespace GEOMETRYPAIR


namespace BEAMINTERACTION
{
  /**
   * \brief Class for beam to solid meshtying.
   * @param beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @param solid Type from GEOMETRYPAIR::ElementDiscretization... representing the solid.
   */
  template <typename beam, typename solid>
  class BeamToSolidVolumeMeshtyingPairBase
      : public BeamToSolidPairBase<GEOMETRYPAIR::line_to_volume_scalar_type<beam, solid>, beam,
            solid>
  {
   protected:
    //! Type to be used for scalar AD variables.
    using scalar_type = GEOMETRYPAIR::line_to_volume_scalar_type<beam, solid>;

    //! Shortcut to the base class.
    using base_class = BeamToSolidPairBase<scalar_type, beam, solid>;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidVolumeMeshtyingPairBase();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidVolumeMeshtyingPairBase(){};

    /**
     * \brief Setup the contact pair (derived).
     *
     * This method sets the solid reference positions for this pair. This can not be done in the
     * base class, since the beam-to-surface (which derive from the base class) need a different
     * handling of the solid DOF.
     */
    void Setup() override;

    /**
     * \brief Things that need to be done in a separate loop before the actual evaluation loop over
     * all contact pairs.
     */
    void PreEvaluate() override;

    /**
     * \brief Update state of translational nodal DoFs (absolute positions and tangents) of both
     * elements.
     *
     * Update of the solid positions is performed in this class method, the beam positions are set
     * in the parent class method, which is called here.
     *
     * @param beam_centerline_dofvec
     * @param solid_nodal_dofvec
     */
    void ResetState(const std::vector<double>& beam_centerline_dofvec,
        const std::vector<double>& solid_nodal_dofvec) override;

    /**
     * \brief Set the restart displacement in this pair.
     *
     * If coupling interactions should be evaluated w.r.t the restart state, this method will set
     * them in the pair accordingly.
     *
     * @param centerline_restart_vec_ (in) Vector with the centerline displacements at the restart
     * step, for all contained elements (Vector of vector).
     */
    virtual void SetRestartDisplacement(
        const std::vector<std::vector<double>>& centerline_restart_vec_) override;

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

   protected:
    /**
     * \brief Return a cast of the geometry pair to the type for this contact pair.
     * @return RPC with the type of geometry pair for this beam contact pair.
     */
    inline Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToVolume<double, beam, solid>>
    CastGeometryPair() const
    {
      return Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::GeometryPairLineToVolume<double, beam, solid>>(
          this->geometry_pair_, true);
    };

    /**
     * \brief This function evaluates the penalty force from a given beam position and a given solid
     * position.
     *
     * This method is mainly used for visualization.
     *
     * @param r_beam (in) Position on the beam.
     * @param r_solid (in) Position on the solid.
     * @param force (out) Force acting on the beam (the negative force acts on the solid).
     */
    virtual void EvaluatePenaltyForce(const LINALG::Matrix<3, 1, scalar_type>& r_beam,
        const LINALG::Matrix<3, 1, scalar_type>& r_solid,
        LINALG::Matrix<3, 1, scalar_type>& force) const;

    /**
     * \brief Get the reference position to be used for the evaluation of the coupling terms.
     *
     * @param beam_coupling_ref (out) shifted reference position of the beam.
     * @param solid_coupling_ref (out) shifted reference position of the solid.
     */
    void GetCouplingReferencePosition(LINALG::Matrix<beam::n_dof_, 1, double>& beam_coupling_ref,
        LINALG::Matrix<solid::n_dof_, 1, double>& solid_coupling_ref) const;

   protected:
    //! Flag if the meshtying has been evaluated already.
    bool meshtying_is_evaluated_;

    //! Current nodal positions (and tangents) of the solid.
    LINALG::Matrix<solid::n_dof_, 1, scalar_type> ele2pos_;

    //! Reference nodal positions (and tangents) of the solid.
    LINALG::Matrix<solid::n_dof_, 1, double> ele2posref_;

    //! Offset of solid DOFs for coupling. This will be used when the state that should be coupled
    //! is not the undeformed reference position, i.e. in restart simulations where the restart
    //! state is coupled. This only makes sense for volume mesh tying, which is why we also define
    //! the beam restart DOFs here.
    LINALG::Matrix<beam::n_dof_, 1, double> ele1posref_offset_;
    LINALG::Matrix<solid::n_dof_, 1, double> ele2posref_offset_;
  };
}  // namespace BEAMINTERACTION

#endif

/*----------------------------------------------------------------------*/
/*! \file

\brief Mesh tying for 2D rotational coupling examples with integration over the beam surface.

\level 3
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_GAUSS_POINT_CROSS_SECTION_ROTAION_H_
#define BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_GAUSS_POINT_CROSS_SECTION_ROTAION_H_


#include "beam_to_solid_volume_meshtying_pair_base.H"


// Forward declaration.
namespace GEOMETRYPAIR
{
  template <typename scalar_type, typename line, typename volume>
  class GeometryPairLineToVolumeGaussPointProjectionCrossSection;
}  // namespace GEOMETRYPAIR


namespace BEAMINTERACTION
{
  /**
   * \brief Class for rotational beam to solid mesh tying using Gauss point projection on the
   * (circular) beam cross section.
   * @param beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @param solid Type from GEOMETRYPAIR::ElementDiscretization... representing the solid.
   */
  template <typename beam, typename solid>
  class BeamToSolidVolumeMeshtyingPairGaussPointCrossSectionRotation
      : public BeamToSolidVolumeMeshtyingPairBase<beam, solid>
  {
   private:
    //! Shortcut to the base class.
    using base_class = BeamToSolidVolumeMeshtyingPairBase<beam, solid>;

    //! Type to be used for scalar AD variables. This can not be inherited from the base class.
    using scalar_type = typename base_class::scalar_type;

    //! Number of rotational DOF for the SR beams;
    static const unsigned int n_dof_rot_ = 3;
    static const unsigned int n_dof_pair_ = n_dof_rot_ + beam::n_dof_ + solid::n_dof_;

    //! FAD type to evaluate the rotational coupling terms. The first 3 entries are the values of
    //! psi_beam, the following entries are the discrete solid DOFs.
    using scalar_type_rot_1st = typename Sacado::Fad::SLFad<double, n_dof_pair_>;
    using scalar_type_rot_2nd = typename FADUTILS::HigherOrderFadType<2, scalar_type_rot_1st>::type;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidVolumeMeshtyingPairGaussPointCrossSectionRotation();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidVolumeMeshtyingPairGaussPointCrossSectionRotation(){};

    /**
     * \brief Evaluate the pair and directly assemble it into the global force vector and stiffness
     * matrix (derived).
     *
     * Rotational coupling contributions will be added in this method.
     */
    void EvaluateAndAssemble(const Teuchos::RCP<const DRT::Discretization>& discret,
        const Teuchos::RCP<Epetra_FEVector>& force_vector,
        const Teuchos::RCP<LINALG::SparseMatrix>& stiffness_matrix,
        const Teuchos::RCP<const Epetra_Vector>& displacement_vector) override;

    /**
     * \brief Create the geometry pair for this contact pair. We overload that function because this
     * pair requires explicitly that a cross section projection pair is created.
     * @param geometry_evaluation_data_ptr Evaluation data that will be linked to the pair.
     */
    void CreateGeometryPair(
        const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>& geometry_evaluation_data_ptr)
        override;

   protected:
    /**
     * \brief Calculate the position on the beam, also taking into account parameter coordinates on
     * the cross section.
     * @param integration_point (in) Integration where the position should be evaluated.
     * @param r_beam (out) Position on the beam.
     * @param reference (in) True -> the reference position is calculated, False -> the current
     * position is calculated.
     */
    void EvaluateBeamPosition(const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& integration_point,
        LINALG::Matrix<3, 1, scalar_type>& r_beam, bool reference) const override;
  };
}  // namespace BEAMINTERACTION

#endif

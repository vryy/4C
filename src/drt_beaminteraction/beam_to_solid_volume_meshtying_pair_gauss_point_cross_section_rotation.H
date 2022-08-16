/*----------------------------------------------------------------------*/
/*! \file

\brief Mesh tying for 2D rotational coupling examples with integration over the beam surface.

\level 3
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_GAUSS_POINT_CROSS_SECTION_ROTAION_H_
#define BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_GAUSS_POINT_CROSS_SECTION_ROTAION_H_


#include "beam_to_solid_volume_meshtying_pair_gauss_point_cross_section_base.H"


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
      : public BeamToSolidVolumeMeshtyingPairGaussPointCrossSectionBase<beam, solid>
  {
   private:
    //! Shortcut to the base class.
    using base_class = BeamToSolidVolumeMeshtyingPairGaussPointCrossSectionBase<beam, solid>;

    //! Type to be used for scalar AD variables. This can not be inherited from the base class.
    using scalar_type = typename base_class::scalar_type;

    //! Rotation representation.
    static const unsigned int n_nodes_rot_ = 3;
    static const unsigned int n_dof_rot_ = n_nodes_rot_ * 3;

    //! Number of dependent DOFs for the pair.
    static const unsigned int n_dof_pair_ = beam::n_dof_ + solid::n_dof_ + n_dof_rot_;
    static const unsigned int n_dof_fad_ = 3 + beam::n_dof_ + solid::n_dof_ + n_dof_rot_;

    //! FAD type to evaluate the rotational coupling terms. The first ordering is as follows: beam
    //! DOFs, solid DOFs, 3 components of the rotation vector psi.
    using scalar_type_pair = typename Sacado::Fad::SLFad<double, n_dof_fad_>;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidVolumeMeshtyingPairGaussPointCrossSectionRotation();

    /*!
     *\brief things that need to be done in a separate loop before the actual evaluation loop
     *      over all contact pairs
     */
    virtual void PreEvaluate() override;

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
     * \brief Update state of rotational DoFs of both elements
     */
    void ResetRotationState(const DRT::Discretization& discret,
        const Teuchos::RCP<const Epetra_Vector>& ia_discolnp) override;

   protected:
    /**
     * \brief Get the triad of the beam at the parameter coordinate xi (derived)
     */
    void GetTriadAtXiDouble(
        const double xi, LINALG::Matrix<3, 3, double>& triad, const bool reference) const override;

   private:
    //! Reference triad interpolation in the beam element
    LARGEROTATIONS::TriadInterpolationLocalRotationVectors<3, double>
        triad_interpolation_scheme_ref_;

    //! Current triad interpolation in the beam element
    LARGEROTATIONS::TriadInterpolationLocalRotationVectors<3, double> triad_interpolation_scheme_;
  };
}  // namespace BEAMINTERACTION

#endif

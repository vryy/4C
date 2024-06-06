/*----------------------------------------------------------------------*/
/*! \file

\brief Class for full 2D-3D beam-to-solid volume mesh tying based on a Simo-Reissner beam element.

\level 3
*/
// End doxygen header.


#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_2D_3D_FULL_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_2D_3D_FULL_HPP


#include "4C_config.hpp"

#include "4C_beaminteraction_beam_to_solid_volume_meshtying_pair_2d-3d_base.hpp"

FOUR_C_NAMESPACE_OPEN


namespace BEAMINTERACTION
{
  /**
   * \brief Class for full 2D-3D beam-to-solid volume mesh tying based on a Simo-Reissner beam
   * element.
   * @param beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @param solid Type from GEOMETRYPAIR::ElementDiscretization... representing the solid.
   */
  template <typename beam, typename solid>
  class BeamToSolidVolumeMeshtyingPair2D3DFull
      : public BeamToSolidVolumeMeshtyingPair2D3DBase<beam, solid>
  {
   private:
    //! Shortcut to the base class.
    using base_class = BeamToSolidVolumeMeshtyingPair2D3DBase<beam, solid>;

    //! Type to be used for scalar AD variables. This can not be inherited from the base class.
    using scalar_type = typename base_class::scalar_type;

    //! Rotation representation. Each node has a rotation vector with 3 components.
    static constexpr unsigned int n_nodes_rot_ = 3;
    static constexpr unsigned int n_dof_rot_ = n_nodes_rot_ * 3;

    //! Number of DOFs for the pair.
    static constexpr unsigned int n_dof_pair_ = beam::n_dof_ + solid::n_dof_ + n_dof_rot_;

    //! Number of dependent variables for the pair. The ordering is as follows: first the beam DOFs,
    //! then the solid DOFs, then the components of the cross section rotation vector.
    static constexpr unsigned int n_dof_fad_ = beam::n_dof_ + solid::n_dof_ + n_dof_rot_ + 3;

    //! FAD type to evaluate the rotational coupling terms. The first ordering is as follows: beam
    //! DOFs, solid DOFs, 3 components of the rotation vector psi.
    using scalar_type_pair = typename Sacado::Fad::SLFad<double, n_dof_fad_>;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidVolumeMeshtyingPair2D3DFull() = default;

    /*!
     *\brief things that need to be done in a separate loop before the actual evaluation loop
     *      over all contact pairs
     */
    void pre_evaluate() override;

    /**
     * \brief Evaluate the pair and directly assemble it into the global force vector and stiffness
     * matrix (derived).
     *
     * Rotational coupling contributions will be added in this method.
     */
    void EvaluateAndAssemble(const Teuchos::RCP<const Discret::Discretization>& discret,
        const Teuchos::RCP<Epetra_FEVector>& force_vector,
        const Teuchos::RCP<Core::LinAlg::SparseMatrix>& stiffness_matrix,
        const Teuchos::RCP<const Epetra_Vector>& displacement_vector) override;

    /**
     * \brief Update state of rotational DoFs of both elements
     */
    void ResetRotationState(const Discret::Discretization& discret,
        const Teuchos::RCP<const Epetra_Vector>& ia_discolnp) override;

   protected:
    /**
     * \brief Get the triad of the beam at the parameter coordinate xi (derived)
     */
    void get_triad_at_xi_double(const double xi, Core::LinAlg::Matrix<3, 3, double>& triad,
        const bool reference) const override;

   private:
    //! Reference triad interpolation in the beam element
    LargeRotations::TriadInterpolationLocalRotationVectors<3, double>
        triad_interpolation_scheme_ref_;

    //! Current triad interpolation in the beam element
    LargeRotations::TriadInterpolationLocalRotationVectors<3, double> triad_interpolation_scheme_;
  };
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif

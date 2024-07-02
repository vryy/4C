/*----------------------------------------------------------------------*/
/*! \file

\brief Contact element for mortar contact between a 3D beam and a surface element.

\level 3
*/
// End doxygen header.


#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_SURFACE_CONTACT_PAIR_MORTAR_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_SURFACE_CONTACT_PAIR_MORTAR_HPP

#include "4C_config.hpp"

#include "4C_beaminteraction_beam_to_solid_surface_contact_pair_base.hpp"
#include "4C_inpar_beam_to_solid.hpp"
#include "4C_linalg_sparsematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace BEAMINTERACTION
{
  //! Forward declarations
  class BeamToSolidSurfaceContactParams;

  /**
   * \brief Class for mortar beam to surface surface contact based on a scalar Lagrange multiplier
   * interpolation
   * @tparam scalar_type Type for scalar DOF values.
   * @tparam beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @tparam surface Type from GEOMETRYPAIR::ElementDiscretization... representing the surface.
   * @tparam mortar Type from GEOMETRYPAIR::ElementDiscretization... representing the interpolation
   * of the Lagrange multiplier.
   */
  template <typename ScalarType, typename Beam, typename Surface, typename Mortar>
  class BeamToSolidSurfaceContactPairMortar
      : public BeamToSolidSurfaceContactPairBase<ScalarType, Beam, Surface>
  {
   protected:
    //! Shortcut to the base class.
    using base_class = BeamToSolidSurfaceContactPairBase<ScalarType, Beam, Surface>;

    // Type from GEOMETRYPAIR::ElementDiscretization... representing the interpolation of the
    // Lagrange multiplier variations. For now this is always equal to the primal interpolation.
    using mortar_trial = Mortar;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidSurfaceContactPairMortar();

    /**
     * \brief Destructor.
     */
    ~BeamToSolidSurfaceContactPairMortar() override{};

    /**
     * \brief Evaluate the global matrices and vectors resulting from mortar coupling. (derived)
     */
    void evaluate_and_assemble_mortar_contributions(const Core::FE::Discretization& discret,
        const BeamToSolidMortarManager* mortar_manager,
        Core::LinAlg::SparseMatrix& global_constraint_lin_beam,
        Core::LinAlg::SparseMatrix& global_constraint_lin_solid,
        Core::LinAlg::SparseMatrix& global_force_beam_lin_lambda,
        Core::LinAlg::SparseMatrix& global_force_solid_lin_lambda,
        Epetra_FEVector& global_constraint, Epetra_FEVector& global_kappa,
        Core::LinAlg::SparseMatrix& global_kappa_lin_beam,
        Core::LinAlg::SparseMatrix& global_kappa_lin_solid, Epetra_FEVector& global_lambda_active,
        const Teuchos::RCP<const Epetra_Vector>& displacement_vector) override;

    /**
     * \brief Evaluate the pair and directly assemble it into the global force vector and stiffness
     * matrix (derived).
     */
    void EvaluateAndAssemble(const Core::FE::Discretization& discret,
        const BeamToSolidMortarManager* mortar_manager,
        const Teuchos::RCP<Epetra_FEVector>& force_vector,
        const Teuchos::RCP<Core::LinAlg::SparseMatrix>& stiffness_matrix,
        const Epetra_Vector& global_lambda, const Epetra_Vector& displacement_vector) override;

   private:
    /**
     * @brief Get the Jacobian for the configuration the Lagrange multipliers are defined in
     */
    ScalarType get_jacobian_for_configuration(const ScalarType& eta,
        const Inpar::BeamToSolid::BeamToSolidSurfaceContactMortarDefinedIn mortar_configuration)
        const;

   private:
    //! Integral of the shape function matrix of the beam (transposed) multiplied with the normal
    //! vector multiplied with the shape function matrix of the Lagrange multipliers
    Core::LinAlg::Matrix<Beam::n_dof_, Mortar::n_dof_, ScalarType>
        beam_shape_times_normal_times_lambda_shape_;
    //! Integral of the shape function matrix of the surface (transposed) multiplied with the normal
    //! vector multiplied with the shape function matrix of the Lagrange multipliers
    Core::LinAlg::Matrix<Surface::n_dof_, Mortar::n_dof_, ScalarType>
        surface_shape_times_normal_times_lambda_shape_;
  };

  /**
   * \brief Factory function for beam-to-solid contact mortar pairs.
   */
  Teuchos::RCP<BEAMINTERACTION::BeamContactPair> BeamToSolidSurfaceContactPairMortarFactory(
      const Teuchos::RCP<const BeamToSolidSurfaceContactParams> beam_to_surface_contact_params,
      const Core::FE::CellType& surface_type, const bool beam_is_hermite);

}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif

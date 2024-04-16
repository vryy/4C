/*----------------------------------------------------------------------*/
/*! \file

\brief Meshtying element for meshtying between a 3D beam and a 3D solid element using mortar shape
functions for the traction.

\level 3
*/
// End doxygen header.


#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_MORTAR_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_VOLUME_MESHTYING_PAIR_MORTAR_HPP


#include "baci_config.hpp"

#include "baci_beaminteraction_beam_to_solid_volume_meshtying_pair_base.hpp"

// Forward declarations.
namespace
{
  class BeamToSolidVolumeMeshtyingPairMortarTest;
}

BACI_NAMESPACE_OPEN

namespace BEAMINTERACTION
{
  /**
   * \brief Class for beam to solid meshtying using mortar shape functions for the contact
   * tractions.
   * @param beam Type from GEOMETRYPAIR::ElementDiscretization... representing the beam.
   * @param solid Type from GEOMETRYPAIR::ElementDiscretization... representing the solid.
   * @param mortar Type from BEAMINTERACTION::ElementDiscretization... representing the mortar shape
   * functions.
   */
  template <typename beam, typename solid, typename mortar>
  class BeamToSolidVolumeMeshtyingPairMortar
      : public BeamToSolidVolumeMeshtyingPairBase<beam, solid>
  {
    //! Define the unit test class as friend so it can set up a valid pair state for the test cases.
    friend BeamToSolidVolumeMeshtyingPairMortarTest;

   protected:
    //! Shortcut to the base class.
    using base_class = BeamToSolidVolumeMeshtyingPairBase<beam, solid>;

    //! Type to be used for scalar AD variables.
    using scalar_type = typename base_class::scalar_type;

   public:
    /**
     * \brief Standard Constructor
     */
    BeamToSolidVolumeMeshtyingPairMortar();


    /**
     * \brief Evaluate the global matrices and vectors resulting from mortar coupling. (derived)
     */
    void EvaluateAndAssembleMortarContributions(const DRT::Discretization& discret,
        const BeamToSolidMortarManager* mortar_manager, CORE::LINALG::SparseMatrix& global_G_B,
        CORE::LINALG::SparseMatrix& global_G_S, CORE::LINALG::SparseMatrix& global_FB_L,
        CORE::LINALG::SparseMatrix& global_FS_L, Epetra_FEVector& global_constraint,
        Epetra_FEVector& global_kappa, Epetra_FEVector& global_lambda_active,
        const Teuchos::RCP<const Epetra_Vector>& displacement_vector) override;

    /**
     * \brief This pair enforces constraints via a mortar-type method, which requires an own
     * assembly method (provided by the mortar manager).
     */
    inline bool IsAssemblyDirect() const override { return false; };

    /**
     * \brief Add the visualization of this pair to the beam to solid visualization output writer.
     * This will add mortar specific data to the output.
     * @param visualization_writer (out) Object that manages all visualization related data for beam
     * to solid pairs.
     * @param visualization_params (in) Parameter list (not used in this class).
     */
    void GetPairVisualization(
        Teuchos::RCP<BeamToSolidVisualizationOutputWriterBase> visualization_writer,
        Teuchos::ParameterList& visualization_params) const override;

   protected:
    /**
     * \brief Evaluate the local mortar matrices for this contact element pair.
     */
    void EvaluateDM(CORE::LINALG::Matrix<mortar::n_dof_, beam::n_dof_, double>& local_D,
        CORE::LINALG::Matrix<mortar::n_dof_, solid::n_dof_, double>& local_M,
        CORE::LINALG::Matrix<mortar::n_dof_, 1, double>& local_kappa,
        CORE::LINALG::Matrix<mortar::n_dof_, 1, double>& local_constraint) const;

    /**
     * \brief For the mortar pairs it does not make sense to calculate forces at the integration
     * points.
     * @param r_beam (in) Position on the beam.
     * @param r_solid (in) Position on the solid.
     * @param force (out) Return 0 by default.
     */
    void EvaluatePenaltyForceDouble(const CORE::LINALG::Matrix<3, 1, double>& r_beam,
        const CORE::LINALG::Matrix<3, 1, double>& r_solid,
        CORE::LINALG::Matrix<3, 1, double>& force) const override;

   protected:
    //! Number of rotational Lagrange multiplier DOFS.
    unsigned int n_mortar_rot_;
  };
}  // namespace BEAMINTERACTION

BACI_NAMESPACE_CLOSE

#endif

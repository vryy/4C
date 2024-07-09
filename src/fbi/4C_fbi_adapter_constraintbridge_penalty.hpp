/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation connecting the penalty constraint enforcement technique with a discretization
approach for Fluid-beam interaction.

\level 2

*----------------------------------------------------------------------*/

#ifndef FOUR_C_FBI_ADAPTER_CONSTRAINTBRIDGE_PENALTY_HPP
#define FOUR_C_FBI_ADAPTER_CONSTRAINTBRIDGE_PENALTY_HPP

#include "4C_config.hpp"

#include "4C_fbi_adapter_constraintbridge.hpp"

#include <Epetra_FEVector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace BEAMINTERACTION
{
  class BeamToFluidMeshtyingVtkOutputWriter;
}

namespace Core::LinAlg
{
  class SparseMatrix;
  class SparseOperator;
}  // namespace Core::LinAlg
namespace Adapter
{
  /**
   *   \brief Implementation connecting the penalty constraint enforcement technique with a
   *   discretization approach for Fluid-beam interaction
   *
   */
  class FBIConstraintBridgePenalty : public FBIConstraintBridge
  {
    friend class ConstraintEnforcerFactory;
    friend class FBIPenaltyConstraintenforcer;
    friend class BEAMINTERACTION::BeamToFluidMeshtyingVtkOutputWriter;
    friend class FBIConstraintenforcer;

   public:
    /**
     * \brief Initializes all members of the class     *
     */
    void setup(const Epetra_Map* beam_map, const Epetra_Map* fluid_map,
        Teuchos::RCP<Core::LinAlg::SparseOperator> fluidmatrix, bool fluidmeshtying) override;

    /**
     * \brief Computes the coupling matrices
     *
     * This is where the magic happens. The global meshtying contributions are integrated using
     * information of the beam elements, the fluid elements and their position relative to each
     * other.
     *
     */

    void evaluate(Teuchos::RCP<const Core::FE::Discretization> discretization1,
        Teuchos::RCP<const Core::FE::Discretization> discretization2,
        Teuchos::RCP<const Epetra_Vector> fluid_vel,
        Teuchos::RCP<const Epetra_Vector> beam_vel) override;

    /// resets the matrices and vectors to zero
    void reset_bridge() override;

    void prepare_fluid_solve() override { set_weak_dirichlet_flag(); };

    /// Matrix containing only structure side contributions \f$C_{ss}\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> get_css() const override { return css_; };

    /// Matrix containing only fluid side contributions \f$C_{ff}\f$
    Teuchos::RCP<const Core::LinAlg::SparseOperator> get_cff() const override { return cff_; };

    /// Matrix containing mixed fluid side contributions \f$C_{fs}\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> get_cfs() const override { return cfs_; };

    /// Matrix containing mixed structure side contributions \f$C_{sf}\f$
    Teuchos::RCP<const Core::LinAlg::SparseMatrix> get_csf() const override { return csf_; };

    /// Negative RHS coupling contribution for the fluid partition \f$f_f\f$
    Teuchos::RCP<const Epetra_FEVector> get_fluid_coupling_residual() const override
    {
      return ff_;
    };

    /// Force vector acting on the structure side \f$f_s\f$
    Teuchos::RCP<const Epetra_FEVector> get_structure_coupling_residual() const override
    {
      return fs_;
    };

   protected:
    /** \brief You will have to use the Adapter::ConstraintEnforcerFactory
     *
     */
    FBIConstraintBridgePenalty()
        : css_(Teuchos::null),
          cff_(Teuchos::null),
          cfs_(Teuchos::null),
          csf_(Teuchos::null),
          ff_(Teuchos::null),
          fs_(Teuchos::null),
          fluid_scaled_(false),
          structure_scaled_(false){};

    /**
     * \brief Sets the flag to compute only force contributions from the beam
     *
     * This allows for a more efficient implementation for the assembly of weak dirichlet
     * contributions to the fluid field, since it avoids a global multiplication of the stiffness
     * $\f C_sf \f$ matrix with the structure velocity
     */
    void set_weak_dirichlet_flag();

    /// Sets the flag to compute force contributions from beam and fluid
    void unset_weak_dirichlet_flag();

    /// Scales all structure vectors and martrices with penalty
    void scale_penalty_structure_contributions();

    /// Scales all fluid vectors and martrices with penalty
    void scale_penalty_fluid_contributions();

   private:
    /// Coupling matrix containing only structure side contributions \f$C_ss\f$
    Teuchos::RCP<Core::LinAlg::SparseMatrix> css_;

    /// Coupling matrix containing only fluid side contributions \f$C_ff\f$
    Teuchos::RCP<Core::LinAlg::SparseOperator> cff_;

    /// Coupling matrix containing mixed fluid side contributions \f$C_fs\f$
    Teuchos::RCP<Core::LinAlg::SparseMatrix> cfs_;

    /// Coupling matrix containing mixed structure side contributions \f$C_sf\f$
    Teuchos::RCP<Core::LinAlg::SparseMatrix> csf_;

    /// Force vector acting on the fluid side \f$f_f\f$
    Teuchos::RCP<Epetra_FEVector> ff_;

    /// Force vector acting on the structure side \f$f_s\f$
    Teuchos::RCP<Epetra_FEVector> fs_;

    /// Bool to keep track if the fluid coupling contributions were already scaled with the penalty
    /// parameter
    bool fluid_scaled_;

    /// Bool to keep track if the structure coupling contributions were already scaled with the
    /// penalty parameter
    bool structure_scaled_;
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif

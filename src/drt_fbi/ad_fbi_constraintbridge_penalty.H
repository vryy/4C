/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation connecting the penalty constraint enforcement technique with a discretization
approach for Fluid-beam interaction.

\level 2

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/

#ifndef AD_FBI_CONSTRAINTBRIDGE_PENALTY_H_
#define AD_FBI_CONSTRAINTBRIDGE_PENALTY_H_

#include <Teuchos_RCP.hpp>
#include <Epetra_FEVector.h>
#include "ad_fbi_constraintbridge.H"

namespace BEAMINTERACTION
{
  class BeamToFluidMeshtyingVtkOutputWriter;
}

namespace LINALG
{
  class SparseMatrix;
}
namespace ADAPTER
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
    /// empty destructor
    virtual ~FBIConstraintBridgePenalty(){};

    /**
     * \brief Initializes all members of the class     *
     */
    virtual void Setup(const Epetra_Map* beam_map, const Epetra_Map* fluid_map);

    /**
     * \brief Computes the coupling matrices
     *
     * This is where the magic happens. The global meshtying contributions are integrated using
     * information of the beam elements, the fluid elements and their position relative to each
     * other.
     *
     */

    void Evaluate(Teuchos::RCP<const DRT::Discretization> discretization1,
        Teuchos::RCP<const DRT::Discretization> discretization2,
        Teuchos::RCP<const Epetra_Vector> fluid_vel,
        Teuchos::RCP<const Epetra_Vector> beam_vel) override;

    /// resets the matrices and vectors to zero
    void ResetBridge() override;

    void PrepareFluidSolve() override { SetWeakDirichletFlag(); };

    /// Matrix containing only structure side contributions \f$C_ss\f$
    Teuchos::RCP<const LINALG::SparseMatrix> GetCss() const override { return Css_; };

    /// Matrix containing only fluid side contributions \f$C_ff\f$
    Teuchos::RCP<const LINALG::SparseMatrix> GetCff() const override { return Cff_; };

    /// Matrix containing mixed fluid side contributions \f$C_fs\f$
    Teuchos::RCP<const LINALG::SparseMatrix> GetCfs() const override { return Cfs_; };

    /// Matrix containing mixed structure side contributions \f$C_sf\f$
    Teuchos::RCP<const LINALG::SparseMatrix> GetCsf() const override { return Csf_; };

    /// Force vector acting on the fluid side \f$f_f\f$
    Teuchos::RCP<const Epetra_FEVector> GetFluidCouplingResidual() const override { return ff_; };

    /// Force vector acting on the structure side \f$f_s\f$
    Teuchos::RCP<const Epetra_FEVector> GetStructureCouplingResidual() const override
    {
      return fs_;
    };

   protected:
    /** \brief You will have to use the ADAPTER::ConstraintEnforcerFactory
     *
     */
    FBIConstraintBridgePenalty()
        : Css_(Teuchos::null),
          Cff_(Teuchos::null),
          Cfs_(Teuchos::null),
          Csf_(Teuchos::null),
          ff_(Teuchos::null),
          fs_(Teuchos::null){};

    /**
     * \brief Sets the flag to compute only force contributions from the beam
     *
     * This allows for a more efficient implementation for the assembly of weak dirichlet
     * contributions to the fluid field, since it avoids a global multiplication of the stiffness
     * $\f C_sf \f$ matrix with the structure velocity
     */
    void SetWeakDirichletFlag();

    /// Sets the flag to compute force contributions from beam and fluid
    void UnsetWeakDirichletFlag();

    /// Scales all vectors and martrices with penalty
    void ScalePenaltyContributions();

   private:
    /// Coupling matrix containing only structure side contributions \f$C_ss\f$
    Teuchos::RCP<LINALG::SparseMatrix> Css_;

    /// Coupling matrix containing only fluid side contributions \f$C_ff\f$
    Teuchos::RCP<LINALG::SparseMatrix> Cff_;

    /// Coupling matrix containing mixed fluid side contributions \f$C_fs\f$
    Teuchos::RCP<LINALG::SparseMatrix> Cfs_;

    /// Coupling matrix containing mixed structure side contributions \f$C_sf\f$
    Teuchos::RCP<LINALG::SparseMatrix> Csf_;

    /// Force vector acting on the fluid side \f$f_f\f$
    Teuchos::RCP<Epetra_FEVector> ff_;

    /// Force vector acting on the structure side \f$f_s\f$
    Teuchos::RCP<Epetra_FEVector> fs_;
  };
}  // namespace ADAPTER

#endif /* AD_FBI_CONSTRAINTBRIDGE_PENALTY_H_ */

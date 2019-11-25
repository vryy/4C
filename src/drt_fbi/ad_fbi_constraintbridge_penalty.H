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
    friend class FBIPenaltyConstraintenforcer;                          // todo why?
    friend class BEAMINTERACTION::BeamToFluidMeshtyingVtkOutputWriter;  // todo why?

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
     * other. Depending on
     *
     */

    virtual void Evaluate(
        Teuchos::RCP<const std::vector<Teuchos::RCP<DRT::Discretization>>> discretizations)
        override;

    /// resets the matrices and vectors to zero
    virtual void ResetBridge() override;

    virtual void PrepareFluidSolve() override { SetWeakDirichletFlag(); };

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

    /// Matrix containing only structure side contributions \f$C_ss\f$
    Teuchos::RCP<LINALG::SparseMatrix> const GetCss() { return Css_; };

    /// Matrix containing only fluid side contributions \f$C_ff\f$
    Teuchos::RCP<LINALG::SparseMatrix> const GetCff() { return Cff_; };

    /// Matrix containing mixed fluid side contributions \f$C_fs\f$
    Teuchos::RCP<LINALG::SparseMatrix> const GetCfs() { return Cfs_; };

    /// Matrix containing mixed structure side contributions \f$C_sf\f$
    Teuchos::RCP<LINALG::SparseMatrix> const GetCsf() { return Csf_; };

    /// Force vector acting on the fluid side \f$f_f\f$
    Teuchos::RCP<Epetra_FEVector> const GetFf() { return ff_; };

    /// Force vector acting on the structure side \f$f_s\f$
    Teuchos::RCP<Epetra_FEVector> const GetFs() { return fs_; };

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

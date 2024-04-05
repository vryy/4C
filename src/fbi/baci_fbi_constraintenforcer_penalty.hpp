/*----------------------------------------------------------------------*/
/*! \file

\brief Implements the constraintenforcement technique of a penalty approach (Mortar and GPTS) (for
fluid-beam interaction).

\level 3

*----------------------------------------------------------------------*/
#ifndef FOUR_C_FBI_CONSTRAINTENFORCER_PENALTY_HPP
#define FOUR_C_FBI_CONSTRAINTENFORCER_PENALTY_HPP

#include "baci_config.hpp"

#include "baci_fbi_constraintenforcer.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

namespace ADAPTER
{
  class FSIStructureWrapper;
  class FluidMovingBoundary;

}  // namespace ADAPTER
namespace CORE::GEO
{
  class SearchTree;
}

namespace FBI
{
  class FBIGeometryCoupler;
}

namespace CORE::LINALG
{
  class SparseOperator;
  class SparseMatrix;
}  // namespace CORE::LINALG
namespace ADAPTER
{
  class ConstraintEnforcerFactory;

  /**
   * \brief Implements the constraint enforcement technique of a penalty approach (Mortar and GPTS)
   * (for fluid-beam interaction)
   *
   * This class decides which information to pass to the actual (partitioned) algorithm.
   * It assembles the stiffness and force contributions of a penalty constraintenforcement approach
   * coupling the velocities of to non-matching meshes
   *
   */
  class FBIPenaltyConstraintenforcer : public FBIConstraintenforcer
  {
    friend ConstraintEnforcerFactory;

   public:
    /**
     * \brief Sets up the constraint enforcer
     *
     *\param[in] structure wrapper for the structure solver
     *\param[in] fluid moving boundary wrapper for the fluid solver
     */
    void Setup(Teuchos::RCP<ADAPTER::FSIStructureWrapper> structure,
        Teuchos::RCP<ADAPTER::FluidMovingBoundary> fluid) override;

    /// Handle fbi specific output
    void Output(double time, int step) override;

   protected:
    /** \brief You will have to use the ADAPTER::ConstraintEnforcerFactory
     *
     * \param[in] bridge an object managing the pair contributins
     * \param[in] geometrycoupler an object managing the search, parallel communication, ect.
     */

    FBIPenaltyConstraintenforcer(Teuchos::RCP<ADAPTER::FBIConstraintBridge> bridge,
        Teuchos::RCP<FBI::FBIGeometryCoupler> geometrycoupler)
        : FBIConstraintenforcer(bridge, geometrycoupler){};

    /**
     * \brief Computes the contributions to the stiffness matrix of the structure field.
     *
     *
     * \returns coupling contributions to the structure system matrix
     */
    Teuchos::RCP<const CORE::LINALG::SparseMatrix> AssembleStructureCouplingMatrix() const override;

    /**
     * \brief Computes the contributions to the stiffness matrix of the fluid field.
     *
     * \returns coupling contributions to the fluid system matrix
     */
    Teuchos::RCP<const CORE::LINALG::SparseOperator> AssembleFluidCouplingMatrix() const override;

    /**
     * \brief Computes the contributions to the rhs of the structure field.
     *
     * \returns coupling contributions to the structure residual
     */
    Teuchos::RCP<Epetra_Vector> AssembleStructureCouplingResidual() const override;

    /**
     * \brief Computes the contributions to the rhs of the slave field.
     *
     * \returns coupling contributions to the fluid residual
     */
    Teuchos::RCP<Epetra_Vector> AssembleFluidCouplingResidual() const override;

    /// Interface to do preparations to solve the fluid
    void PrepareFluidSolve() override;

    /// Output the constraint violation
    virtual void PrintViolation(double time, int step);

   private:
    FBIPenaltyConstraintenforcer() = delete;
  };
}  // namespace ADAPTER
BACI_NAMESPACE_CLOSE

#endif

/*----------------------------------------------------------------------*/
/*! \file

\brief Implements the constraintenforcement technique of a penalty approach (Mortar and GPTS) (for
fluid-beam interaction).

\level 3

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/
#ifndef CONSTRAINTENFORCER_FBI_PENALTY_H
#define CONSTRAINTENFORCER_FBI_PENALTY_H

#include "Teuchos_RCP.hpp"
#include "constraintenforcer_fbi.H"
#include <Epetra_Vector.h>

namespace ADAPTER
{
  class FSIStructureWrapper;
  class FluidMovingBoundary;

}  // namespace ADAPTER
namespace GEO
{
  class SearchTree;
}

namespace FBI
{
  class FBIGeometryCoupler;
}

namespace LINALG
{
  class SparseMatrix;
}
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
    /// empty destructor
    virtual ~FBIPenaltyConstraintenforcer(){};

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
    Teuchos::RCP<const LINALG::SparseMatrix> AssembleStructureCouplingMatrix() const override;

    /**
     * \brief Computes the contributions to the stiffness matrix of the fluid field.
     *
     * \returns coupling contributions to the fluid system matrix
     */
    Teuchos::RCP<const LINALG::SparseMatrix> AssembleFluidCouplingMatrix() const override;

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
#endif

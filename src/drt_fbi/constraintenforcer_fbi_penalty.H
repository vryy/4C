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

   protected:
    /** \brief You will have to use the ADAPTER::ConstraintEnforcerFactory
     *
     * \params[in] bridge an object managing the pair contributins
     * \params[in] geometrycoupler an object managing the search, parallel communication, ect.
     */

    FBIPenaltyConstraintenforcer(Teuchos::RCP<ADAPTER::FBIConstraintBridge> bridge,
        Teuchos::RCP<FBI::FBIGeometryCoupler> geometrycoupler)
        : FBIConstraintenforcer(bridge, geometrycoupler){};

    /**
     * \brief Computes the contributions to the stiffness matrix of the structure field.
     *
     *
     * \param[inout] k system/stiffness matrix
     */
    virtual Teuchos::RCP<LINALG::SparseMatrix> AssembleStructureStiffness() const override;

    /**
     * \brief Computes the contributions to the stiffness matrix of the fluid field.
     *
     * \returns k system/stiffness matrix
     */
    virtual Teuchos::RCP<LINALG::SparseMatrix> AssembleFluidStiffness() const override;

    /**
     * \brief Computes the contributions to the rhs of the structure field.
     *
     * \returns f force/rhs matrix
     */
    virtual Teuchos::RCP<Epetra_Vector> AssembleStructureForce() const override;

    /**
     * \brief Computes the contributions to the rhs of the slave field.
     *
     * \returns f force/rhs matrix
     */
    virtual Teuchos::RCP<Epetra_Vector> AssembleFluidForce() const override;

    /// Interface to do preparations to solve the fluid
    virtual void PrepareFluidSolve() override;

   private:
    FBIPenaltyConstraintenforcer() = delete;
  };
}  // namespace ADAPTER
#endif

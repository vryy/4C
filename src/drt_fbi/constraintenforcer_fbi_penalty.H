/*----------------------------------------------------------------------*/
/*! \file
\file constraintenforcer_fbi_penalty.H

\brief Implements the constraintenforcement technique of a penalty approach (Mortar and GPTS) (for
fluid-beam interaction). todo only for FBI?

\level 3

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/
#ifndef CONSTRAINTENFORCER_FBI_PENALTY_H
#define CONSTRAINTENFORCER_FBI_PENALTY_H

#include "Teuchos_RCP.hpp"
#include "../linalg/linalg_utils.H"
#include "constraintenforcer_fbi.H"

namespace ADAPTER
{
  class FSIStructureWrapper;
  class FluidMovingBoundary;

}  // namespace ADAPTER
namespace GEO
{
  class SearchTree;
}

namespace ADAPTER
{
  class ConstraintEnforcerFactory;

  class FBIPenaltyConstraintenforcer : public FBIConstraintenforcer  // todo fix this
  {
    friend ConstraintEnforcerFactory;

   public:
    /// empty destructor
    virtual ~FBIPenaltyConstraintenforcer(){};

   protected:
    /// Standard constructor
    FBIPenaltyConstraintenforcer(Teuchos::RCP<ADAPTER::FBIConstraintBridge> bridge)
        : FBIConstraintenforcer(bridge){};

    // Do we need set-routines here?

    /**
     * \brief Computes the contributions to the stiffness matrix of the slave field.
     *
     * This has to be implemented differently depending on the concrete constraint enforcement
     * strategy.
     *
     * \param[inout] k system/stiffness matrix
     */
    virtual Teuchos::RCP<LINALG::SparseMatrix> AssembleMasterStiffness() const override;

    /**
     * \brief Computes the contributions to the stiffness matrix of the master field.
     *
     * This has to be implemented differently depending on the concrete constraint enforcement
     * strategy.
     *
     * \returns k system/stiffness matrix
     */
    virtual Teuchos::RCP<LINALG::SparseMatrix> AssembleSlaveStiffness() const override;

    /**
     * \brief Computes the contributions to the residuum of the master field.
     *
     * This has to be implemented differently depending on the concrete constraint enforcement
     * strategy.
     *
     * \returns f force/rhs matrix
     */
    virtual Teuchos::RCP<Epetra_Vector> AssembleMasterForce() const override;

    /**
     * \brief Computes the contributions to the residuum of the slave field.
     *
     * This has to be implemented differently depending on the concrete constraint enforcement
     * strategy.
     *
     * \returns f force/rhs matrix
     */
    virtual Teuchos::RCP<Epetra_Vector> AssembleSlaveForce() override;

    /// Interface to do preparations to solve the fluid
    virtual void PrepareFluidSolve() override;
  };
}  // namespace ADAPTER
#endif

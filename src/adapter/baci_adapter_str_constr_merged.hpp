/*----------------------------------------------------------------------*/
/*! \file

\brief Adapter Layer for Structures with Algebraic Constraints

\level 2

*/

/*----------------------------------------------------------------------*/
/* macros */


#ifndef FOUR_C_ADAPTER_STR_CONSTR_MERGED_HPP
#define FOUR_C_ADAPTER_STR_CONSTR_MERGED_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "baci_config.hpp"

#include "baci_adapter_str_fsiwrapper.hpp"

#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


// forward declarations
namespace STR
{
  namespace AUX
  {
    class MapExtractor;
  }
}  // namespace STR

namespace CORE::LINALG
{
  class MapExtractor;
}

/*----------------------------------------------------------------------*/
namespace ADAPTER
{
  /*====================================================================*/
  /*!
   * \brief Adapter to constrained structural time integration.
   * This class wraps one of the standard adapters for structural time
   * integration. The results are modified and/or merged to account for the
   * additional degrees of freedom of the lagrange multipliers.
   *
   * \date 11/08
   */
  class StructureConstrMerged : public FSIStructureWrapper
  {
   public:
    /// Constructor
    StructureConstrMerged(Teuchos::RCP<Structure> stru);

    /// setup this object
    void Setup() override;

    /// initial guess of Newton's method
    Teuchos::RCP<const Epetra_Vector> InitialGuess() override;

    /// right-hand-side of Newton's method
    Teuchos::RCP<const Epetra_Vector> RHS() override;

    /// unknown displacements at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> Dispnp() const override;

    /// known displacements at \f$t_{n}\f$
    Teuchos::RCP<const Epetra_Vector> Dispn() const override;

    /*! \brief known velocity at \f$t_{n}\f$
     *
     *  Lagrange multiplier does not have a time derivative. Though we need a map
     *  including the Lagrange multiplier, thus, we include it and set it to zero.
     */
    Teuchos::RCP<const Epetra_Vector> Veln() const override;

    /*! known acceleration at \f$t_{n}\f$
     *
     *  Lagrange multiplier does not have a time derivative. Though we need a map
     *  including the Lagrange multiplier, thus, we include it and set it to zero.
     */
    Teuchos::RCP<const Epetra_Vector> Accn() const override;

    /// dof map of vector of unknowns
    Teuchos::RCP<const Epetra_Map> DofRowMap() override;

    /// apply interface forces to structural solver
    ///
    /// This prepares a new solve of the structural field within one time
    /// step. The middle values are newly created.
    ///
    /// \note This is not yet the most efficient implementation.
    void ApplyInterfaceForcesTemporaryDeprecated(Teuchos::RCP<Epetra_Vector> iforce) override;

    /// direct access to system matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> SystemMatrix() override;

    /// direct access to system matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> BlockSystemMatrix() override;

    /// update displacement and evaluate elements
    void Evaluate(Teuchos::RCP<const Epetra_Vector>
            dispstepinc  ///< solution increment between time step n and n+1
        ) override;

    //! Return MapExtractor for Dirichlet boundary conditions
    Teuchos::RCP<const CORE::LINALG::MapExtractor> GetDBCMapExtractor() override
    {
      return structure_->GetDBCMapExtractor();
    };

    /// domain map of system matrix
    const Epetra_Map& DomainMap() const override;

    /// are there any algebraic constraints?
    bool HaveConstraint() override { return structure_->HaveConstraint(); };

    /// Return bool indicating if constraints are defined
    Teuchos::RCP<CONSTRAINTS::ConstrManager> GetConstraintManager() override
    {
      return structure_->GetConstraintManager();
    };

    INPAR::STR::StcScale GetSTCAlgo() override { return structure_->GetSTCAlgo(); };

    Teuchos::RCP<CORE::LINALG::SparseMatrix> GetSTCMat() override
    {
      dserror("FSI with merged structural constraints does not work in combination with STC!");
      return structure_->GetSTCMat();
    }

    //! Update iteration
    //! Add residual increment to Lagrange multipliers stored in Constraint manager
    void UpdateIterIncrConstr(
        Teuchos::RCP<Epetra_Vector> lagrincr  ///< Lagrange multiplier increment
        ) override
    {
      structure_->UpdateIterIncrConstr(lagrincr);
    }

    /// @name Apply interface forces

    //@}

    /// Integrate from t1 to t2
    int Integrate() override { return structure_->Integrate(); }

   private:
    /// the constraint map setup for full <-> stuct+constr transition
    Teuchos::RCP<CORE::LINALG::MapExtractor> conmerger_;

    /// the complete non-overlapping degree of freedom row map for structure and lagrange
    /// multipliers
    Teuchos::RCP<Epetra_Map> dofrowmap_;

    /// @name local copies of input parameters
    //{@
    Teuchos::RCP<DRT::Discretization> discret_;  ///< the discretization
    Teuchos::RCP<Teuchos::ParameterList>
        sdynparams_;  ///< dynamic control flags ... used, but could/should be circumvented
    Teuchos::RCP<Teuchos::ParameterList> xparams_;   ///< eXtra input parameters
    Teuchos::RCP<CORE::LINALG::Solver> solver_;      ///< the linear solver
    Teuchos::RCP<IO::DiscretizationWriter> output_;  ///< the output writer

    //@}

    /// flag indicating if Setup() was called
    bool issetup_;

  };  // class StructureConstrained

}  // namespace ADAPTER

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif

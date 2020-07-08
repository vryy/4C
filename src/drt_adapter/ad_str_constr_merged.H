/*----------------------------------------------------------------------*/
/*! \file

\brief Adapter Layer for Structures with Algebraic Constraints

\level 2

*/

/*----------------------------------------------------------------------*/
/* macros */


#ifndef ADAPTER_STRUCTURE_CONSTR_MERGED_H
#define ADAPTER_STRUCTURE_CONSTR_MERGED_H

/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "ad_str_fsiwrapper.H"


// forward declarations
namespace STR
{
  namespace AUX
  {
    class MapExtractor;
  }
}  // namespace STR

namespace LINALG
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
    virtual void Setup();

    /// initial guess of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess();

    /// right-hand-side of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> RHS();

    /// unknown displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp() const;

    /// known displacements at \f$t_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispn() const;

    /*! \brief known velocity at \f$t_{n}\f$
     *
     *  Lagrange multiplier does not have a time derivative. Though we need a map
     *  including the Lagrange multiplier, thus, we include it and set it to zero.
     */
    virtual Teuchos::RCP<const Epetra_Vector> Veln() const;

    /*! known acceleration at \f$t_{n}\f$
     *
     *  Lagrange multiplier does not have a time derivative. Though we need a map
     *  including the Lagrange multiplier, thus, we include it and set it to zero.
     */
    virtual Teuchos::RCP<const Epetra_Vector> Accn() const;

    /// dof map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap();

    /// apply interface forces to structural solver
    ///
    /// This prepares a new solve of the structural field within one time
    /// step. The middle values are newly created.
    ///
    /// \note This is not yet the most efficient implementation.
    virtual void ApplyInterfaceForcesTemporaryDeprecated(Teuchos::RCP<Epetra_Vector> iforce);

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

    /// update displacement and evaluate elements
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector>
            dispstepinc  ///< solution increment between time step n and n+1
    );

    //! Return MapExtractor for Dirichlet boundary conditions
    virtual const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
    {
      return structure_->GetDBCMapExtractor();
    };

    /// domain map of system matrix
    virtual const Epetra_Map& DomainMap() const;

    /// are there any algebraic constraints?
    virtual bool HaveConstraint() { return structure_->HaveConstraint(); };

    /// Return bool indicating if constraints are defined
    virtual const Teuchos::RCP<UTILS::ConstrManager> GetConstraintManager()
    {
      return structure_->GetConstraintManager();
    };

    virtual INPAR::STR::STC_Scale GetSTCAlgo() { return structure_->GetSTCAlgo(); };

    virtual Teuchos::RCP<LINALG::SparseMatrix> GetSTCMat()
    {
      dserror("FSI with merged structural constraints does not work in combination with STC!");
      return structure_->GetSTCMat();
    }

    //! Update iteration
    //! Add residual increment to Lagrange multipliers stored in Constraint manager
    virtual void UpdateIterIncrConstr(
        Teuchos::RCP<Epetra_Vector> lagrincr  ///< Lagrange multiplier increment
    )
    {
      structure_->UpdateIterIncrConstr(lagrincr);
    }

    /// @name Apply interface forces

    //@}

    /// Integrate from t1 to t2
    virtual int Integrate() { return structure_->Integrate(); }

   private:
    /// the constraint map setup for full <-> stuct+constr transition
    Teuchos::RCP<LINALG::MapExtractor> conmerger_;

    /// the complete non-overlapping degree of freedom row map for structure and lagrange
    /// multipliers
    Teuchos::RCP<Epetra_Map> dofrowmap_;

    /// @name local copies of input parameters
    //{@
    Teuchos::RCP<DRT::Discretization> discret_;  ///< the discretization
    Teuchos::RCP<Teuchos::ParameterList>
        sdynparams_;  ///< dynamic control flags ... used, but could/should be circumvented
    Teuchos::RCP<Teuchos::ParameterList> xparams_;   ///< eXtra input parameters
    Teuchos::RCP<LINALG::Solver> solver_;            ///< the linear solver
    Teuchos::RCP<IO::DiscretizationWriter> output_;  ///< the output writer

    //@}

    /// flag indicating if Setup() was called
    bool issetup_;

  };  // class StructureConstrained

}  // namespace ADAPTER

/*----------------------------------------------------------------------*/
#endif

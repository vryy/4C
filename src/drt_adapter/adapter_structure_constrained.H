/*----------------------------------------------------------------------*/
/*!
\file adapter_structure_constrained.H

\brief Adapter Layer for Structures with Algebraic Constraints

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

#ifndef ADAPTER_STRUCTURE_CONSTRAINED_H
#define ADAPTER_STRUCTURE_CONSTRAINED_H

/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_lib/drt_utils.H"
#include "../drt_structure/stru_resulttest.H"

#include "../drt_structure/strtimint_impl.H"
#include "adapter_structure.H"

#include "adapter_structure_wrapper.H"

using namespace Teuchos;
/*----------------------------------------------------------------------*/
/* adapting adapter */
namespace ADAPTER {

  /*====================================================================*/
  /*!
   * \brief Adapter to constrained structural time integration.
   * This class wraps one of the standard adapters for structural time
   * integration. The results are modified and/or merged to account for the
   * additional degrees of freedom of the lagrange multipliers.
   *
   * \date 11/08
   */
  class StructureConstrained : public StructureWrapper
  {
  public:

    /// Constructor
    StructureConstrained
    (
      RCP<Structure> stru
    );

    /// initial guess of Newton's method
    RCP<const Epetra_Vector> InitialGuess();

    /// right-hand-side of Newton's method
    RCP<const Epetra_Vector> RHS();

    /// unknown displacements at \f$t_{n+1}\f$
    RCP<const Epetra_Vector> Dispnp();

    /// unknown displacements at \f$t_{n}\f$
    RCP<const Epetra_Vector> Dispn();

    /// dof map of vector of unknowns
    Teuchos::RCP<const Epetra_Map> DofRowMap();

    /// direct access to system matrix
    RCP<LINALG::SparseMatrix> SystemMatrix();

    /// External force at t_{n+1}
    RCP<const Epetra_Vector> FExtn();

    /// IT'S ACTUALLY ALWAYS ZERO AND RATHER POINTLESS, OR?
    RCP<const Epetra_Vector> FRobin();

    /// evaluate elements with given displacement
    void Evaluate(
      RCP<const Epetra_Vector> disp  ///< current state
    );

    /// communication object at the interface
    virtual const LINALG::MapExtractor& Interface() const
    {
      return interface_;
    }

    //! Return MapExtractor for Dirichlet boundary conditions
    virtual const RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
    {
      return structure_->GetDBCMapExtractor();
    };

    /// domain map of system matrix (do we really need this?)
    virtual const Epetra_Map& DomainMap();

    /// are there any algebraic constraints?
     virtual bool HaveConstraint()
     {
       return structure_->HaveConstraint();
     };

    /// Return bool indicating if constraints are defined
    virtual const RCP<const UTILS::ConstrManager> GetConstraintManager()
    {
      return structure_->GetConstraintManager();
    };

    //! Update iteration
    //! Add residual increment to Lagrange multipliers stored in Constraint manager
    virtual void UpdateIterIncrConstr
    (
        Teuchos::RCP<Epetra_Vector> lagrincr ///< Lagrange multiplier increment
    )
    {
      structure_->UpdateIterIncrConstr(lagrincr);
    }

    /// @name Apply interface forces

    virtual void ApplyInterfaceRobinValue(
      RCP<Epetra_Vector> iforce,
      RCP<Epetra_Vector> ifluidvel
    );

    //@}

  private:

    /// the interface map setup for interface <-> full translation
    LINALG::MapExtractor interface_;

    /// the constraint map setup for full <-> stuct+constr transition
    LINALG::MapExtractor conmerger_;

    /// the complete non-overlapping degree of freedom row map for structure and lagrange multipliers
    RCP<Epetra_Map> dofrowmap_;

    /// @name local copies of input parameters
    //{@
    RCP<DRT::Discretization> discret_;  ///< the discretisation
    RCP<ParameterList> ioparams_;  ///< I/O flags ... not sure of really needed
    RCP<ParameterList> sdynparams_;  ///< dynamic control flags ... used,
                                                       ///< but could/should be circumvented
    RCP<ParameterList> xparams_;  ///< eXtra input parameters
    RCP<LINALG::Solver> solver_;  ///< the solver
    RCP<IO::DiscretizationWriter> output_;  ///< the output writer

    //@}

  };  // class StructureConstrained

}  // namespace ADAPTER

/*----------------------------------------------------------------------*/
#endif
#endif

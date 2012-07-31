/*----------------------------------------------------------------------*/
/*!
\file ad_str_constr_merged.H

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


#ifndef ADAPTER_STRUCTURE_CONSTR_MERGED_H
#define ADAPTER_STRUCTURE_CONSTR_MERGED_H

/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_structure/stru_resulttest.H"


#include "../drt_structure/strtimint_impl.H"
#include "ad_str_structure.H"

#include "ad_str_fsiwrapper.H"


// forward declarations
namespace STR
{
  namespace AUX
  {
    class MapExtractor;
  }
}

namespace LINALG
{
  class MapExtractor;
}

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
  class StructureConstrMerged : public FSIStructureWrapper
  {
  public:

    /// Constructor
    StructureConstrMerged
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

    /// direct access to system matrix
    RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

    /// External force at t_{n+1}
//    RCP<const Epetra_Vector> FExtn();

    /// update displacement and evaluate elements
    void Evaluate(
      RCP<const Epetra_Vector> dispstepinc  ///< solution increment between time step n and n+1
    );

    /// communication object at the interface
    virtual Teuchos::RCP<const STR::AUX::MapExtractor> Interface() const
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
    virtual const RCP< UTILS::ConstrManager> GetConstraintManager()
    {
      return structure_->GetConstraintManager();
    };

    virtual INPAR::STR::STC_Scale GetSTCAlgo()
    { return structure_->GetSTCAlgo(); };

    virtual Teuchos::RCP<LINALG::SparseMatrix> GetSTCMat()
    {
      dserror("FSI with merged structural constraints does not work in combination with STC!");
      return structure_->GetSTCMat();
    }

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

    //@}

    /// Integrate from t1 to t2
    virtual void Integrate() { structure_->Integrate(); }

  private:

    /// the interface map setup for interface <-> full translation
    Teuchos::RCP<STR::AUX::MapExtractor> interface_;

    /// the constraint map setup for full <-> stuct+constr transition
    Teuchos::RCP<LINALG::MapExtractor> conmerger_;

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

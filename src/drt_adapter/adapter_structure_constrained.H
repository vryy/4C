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
  class StructureConstrained : public Structure
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

    /// direct access to discretization
    RCP<DRT::Discretization> Discretization();

    /// External force at t_{n+1}
    RCP<const Epetra_Vector> FExtn();

    /// IT'S ACTUALLY ALWAYS ZERO AND RATHER POINTLESS, OR?
    RCP<const Epetra_Vector> FRobin();

    /// Get the current time
    virtual double GetTime() const { return structure_->GetTime(); }
    
    /// Get upper limit of time range of interest
    virtual double GetTimeEnd() const { return structure_->GetTimeEnd(); }

    /// Get time step size \f$\Delta t_n\f$
    virtual double GetTimeStepSize() const { return structure_->GetTimeStepSize(); }

    /// Get current time step \f$n\f$ (sometimes called \f$k\f$ [but not in Baci])
    virtual int GetTimeStep() const { return structure_->GetTimeStep(); }

    /// Get number of time steps
    virtual int GetTimeNumStep() const { return structure_->GetTimeNumStep(); }

    /// start new time step
    void PrepareTimeStep();

    /// evaluate elements with given displacement
    void Evaluate(
      RCP<const Epetra_Vector> disp  ///< current state
    );

    /// update at time step end
    void Update();

    /// output results
    void Output();

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

    /// read restart information for given time step
    virtual void ReadRestart
    (
      int step  ///< given time step
    );

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
     
    /// @name Solver calls
    //@{

    /// non-linear solve
    ///
    /// Do the nonlinear solve for the time step. All boundary conditions have
    /// been set.
    virtual void Solve();

    /// get the linear solver object used for this field
    virtual RCP<LINALG::Solver> LinearSolver() { return structure_->LinearSolver(); }

    /// linear structure solve with just a interface load
    ///
    ///  The very special solve done in steepest descent relaxation
    ///  calculation (and matrix free Newton Krylov).
    ///
    ///  \note Can only be called after a valid structural solve.
    virtual RCP<Epetra_Vector> RelaxationSolve(
      RCP<Epetra_Vector> iforce
    );

    //@}

    /// @name Extract interface values

    /// extract displacements
    virtual RCP<Epetra_Vector> ExtractInterfaceDispn();

    /// extract displacements
    virtual RCP<Epetra_Vector> ExtractInterfaceDispnp();

    /// extract interface forces (necessary for robin-BC)
    virtual RCP<Epetra_Vector> ExtractInterfaceForces();

    /// Predictor for interface displacements
    virtual RCP<Epetra_Vector> PredictInterfaceDispnp();

    //@}

    /// @name Apply interface forces

    /// apply interface forces to structural solver
    ///
    ///  This prepares a new solve of the structural field within one time
    ///  step. The middle values are newly created.
    ///
    ///  \note This is not yet the most efficient implementation.
    virtual void ApplyInterfaceForces(
      RCP<Epetra_Vector> iforce
    );

    virtual void ApplyInterfaceRobinValue(
      RCP<Epetra_Vector> iforce,
      RCP<Epetra_Vector> ifluidvel
    );

    //@}

    /// create result test for encapsulated structure algorithm
    virtual RCP<DRT::ResultTest> CreateFieldTest();

  private:

    /// the underlying structural adapter containing the actual structural algorithm
    RCP<ADAPTER::Structure> structure_;

    /// the interface map setup for interface <-> full translation
    LINALG::MapExtractor interface_;
    
    /// the constraint map setup for full <-> stuct+constr transition 
    LINALG::MapExtractor conmerger_;
    
    /// the complete non-overlapping degree of freedom row map for structure and lagrange multipliers 
    RCP<Epetra_Map> dofrowmap_;

    /// sum of displacement vectors already applied,
    /// i.e. the incremental displacements
    RCP<Epetra_Vector> disinc_;

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

/*----------------------------------------------------------------------*/
/*!
\file adapter_structure_timint.H

\brief Structural time integration field adapter

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

#ifndef ADAPTER_STRUCTURE_TIMINT_H
#define ADAPTER_STRUCTURE_TIMINT_H

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


/*----------------------------------------------------------------------*/
/* adapting adapter */
namespace ADAPTER {

  /*====================================================================*/
  /*!
   * \brief Adapter to general implicit structural time integration
   *
   * \date 08/08
   */
  class StructureTimInt : public Structure
  {
  public:

    /// Constructor
    StructureTimInt(
      Teuchos::RCP<Teuchos::ParameterList> ioparams,  ///< I/O flags
      Teuchos::RCP<Teuchos::ParameterList> sdynparams,  ///< input parameters
      Teuchos::RCP<Teuchos::ParameterList> xparams,  ///< extra flags
      Teuchos::RCP<DRT::Discretization> discret,  ///< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,  ///< the solver
      Teuchos::RCP<IO::DiscretizationWriter> output  ///< the output
    );

    /// initial guess of Newton's method
    Teuchos::RCP<const Epetra_Vector> InitialGuess();

    /// right-hand-side of Newton's method
    Teuchos::RCP<const Epetra_Vector> RHS();

    /// unknown displacements at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> Dispnp();

    /// unknown displacements at \f$t_{n}\f$
    Teuchos::RCP<const Epetra_Vector> Dispn();

    /// dof map of vector of unknowns
    Teuchos::RCP<const Epetra_Map> DofRowMap();

    /// direct access to system matrix
    Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

    /// direct access to discretization
    Teuchos::RCP<DRT::Discretization> Discretization();

    /// External force at t_{n+1}
    Teuchos::RCP<const Epetra_Vector> FExtn();

    /// IT'S ACTUALLY ALWAYS ZERO AND RATHER POINTLESS, OR?
    Teuchos::RCP<const Epetra_Vector> FRobin();

    /// @name Time step helpers
    //@{

    /// Get the current time
    virtual double GetTime() const { return structure_->GetTime(); }
    
    /// Get upper limit of time range of interest
    virtual double GetTimeEnd() const { return structure_->GetTimeEnd(); }

    /// Get time step size \f$\Delta t_n\f$
    virtual double GetTimeStepSize() const { return structure_->GetTimeStepSize(); }

    /// Get current time step \f$n\f$ (sometimes called \f$k\f$ [but not in Baci])
    virtual int GetTimeStep() const { return structure_->GetStep(); }

    /// Get number of time steps
    virtual int GetTimeNumStep() const { return structure_->GetTimeNumStep(); }

    /// start new time step
    void PrepareTimeStep();

    /// evaluate elements with given displacement
    void Evaluate(
      Teuchos::RCP<const Epetra_Vector> disp  ///< current state
    );

    /// update at time step end
    void Update();
    
    /// Update iteration
    /// Add residual increment to Lagrange multipliers stored in Constraint manager
    virtual void UpdateIterIncrConstr
    (
      Teuchos::RCP<Epetra_Vector> lagrincr  ///< Lagrange multiplier increment
    )
    {
      structure_->UpdateIterIncrConstr(lagrincr);
    }
    
    /// Return MapExtractor for Dirichlet boundary conditions 
    virtual const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
    {
      return structure_->GetDBCMapExtractor();
    }

    /// output results
    void Output();

    //@}

    /// communication object at the interface
    virtual const LINALG::MapExtractor& Interface() const
    { 
      return interface_;
    }

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
    
    /// get constraint manager defined in the structure
    virtual const RCP<const UTILS::ConstrManager> GetConstraintManager()
    {
      return structure_->GetConstraintManager();
    };
    
    /// @name Solver calls
    //@{

    /// non-linear solve
    ///
    /// Do the nonlinear solve for the time step. All boundary conditions have
    /// been set.
    virtual void Solve();

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

    /// linear structure solve with just a interface load
    ///
    /// The very special solve done in steepest descent relaxation
    /// calculation (and matrix free Newton Krylov).
    ///
    /// \note Can only be called after a valid structural solve.
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(
      Teuchos::RCP<Epetra_Vector> iforce
    );

    //@}

    /// @name Extract interface values

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispn();

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispnp();

    /// extract interface forces (necessary for robin-BC)
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();

    /// Predictor for interface displacements
    virtual Teuchos::RCP<Epetra_Vector> PredictInterfaceDispnp();

    //@}

    /// @name Apply interface forces

    /// apply interface forces to structural solver
    ///
    /// This prepares a new solve of the structural field within one time
    /// step. The middle values are newly created.
    ///
    /// \note This is not yet the most efficient implementation.
    virtual void ApplyInterfaceForces(
      Teuchos::RCP<Epetra_Vector> iforce
    );

    virtual void ApplyInterfaceRobinValue(
      Teuchos::RCP<Epetra_Vector> iforce,
      Teuchos::RCP<Epetra_Vector> ifluidvel
    );

    //@}

    /// create result test for encapsulated structure algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  private:

    /// the actual structure algorithm
    Teuchos::RCP<STR::TimIntImpl> structure_;

    /// the interface map setup for interface <-> full translation
    LINALG::MapExtractor interface_;

    /// sum of displacement vectors already applied,
    /// i.e. the incremental displacements
    Teuchos::RCP<Epetra_Vector> disinc_;

    /// @name local copies of input parameters
    //{@
    Teuchos::RCP<DRT::Discretization> discret_;  ///< the discretisation
    Teuchos::RCP<Teuchos::ParameterList> ioparams_;  ///< I/O flags ... not sure of really needed
    Teuchos::RCP<Teuchos::ParameterList> sdynparams_;  ///< dynamic control flags ... used, 
                                                       ///< but could/should be circumvented
    Teuchos::RCP<Teuchos::ParameterList> xparams_;  ///< eXtra input parameters
    Teuchos::RCP<LINALG::Solver> solver_;  ///< the solver
    Teuchos::RCP<IO::DiscretizationWriter> output_;  ///< the output writer
    //@}

  };  // class StructureTimInt

}  // namespace ADAPTER

/*----------------------------------------------------------------------*/
#endif
#endif

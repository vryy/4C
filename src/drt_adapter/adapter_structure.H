/*----------------------------------------------------------------------*/
/*!
\file adapter_structure.H

\brief Structure field adapter

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ADAPTER_STRUCTURE_H
#define ADAPTER_STRUCTURE_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_lib/drt_utils.H"
#include "../drt_structure/strugenalpha.H"
#include "../drt_structure/stru_resulttest.H"



namespace ADAPTER {

  /// general structural field interface
  /*!

  The point is to keep FSI as far apart from our field solvers as
  possible. Each structure field solver we want to use should get its own
  subclass of this. The FSI algorithm should be able to extract all the
  information from the structure field it needs using this interface.

  All FSI algorithms use this adapter to communicate with the structural
  field. There are different ways to use this adapter.

  In all cases you need to tell the structural algorithm about your time
  step. Therefore PrepareTimeStep(), Update() and Output() must be called at
  the appropriate position in the FSI algorithm.

  <h3>Dirichlet-Neumann coupled FSI</h3>

  A good starting displacement can be guessed with PredictInterfaceDispnp().

  Dirichlet-Neumann coupled FSI will need to Solve() the nonlinear
  structural problem for each time step after the fluid forces have been
  applied (ApplyInterfaceForces()). Solve() will be called many times for each
  time step until the interface equilibrium is reached. The structural
  algorithm has to preserve its state until Update() is called.

  After each Solve() you get the interface forces by ExtractInterfaceDispnp().

  A Dirichlet-Neumann FSI with steepest descent relaxation or matrix free
  Newton Krylov will want to solve the structural problem linearly without
  history and prescribed interface forces: RelaxationSolve().

  <h3>Monolithic FSI</h3>

  Monolithic FSI is based on Evaluate() of elements. This results in a new
  RHS() and a new SysMat(). Together with the InitialGuess() these form the
  building blocks for a block based Newton's method.

  \warning Further cleanup is still needed.

  \sa Fluid, Ale
  \author u.kue
  \date 11/07
  */
  class Structure
  {
  public:

    /// virtual to get polymorph destruction
    virtual ~Structure();

    //! @name Vector access

    /// initial guess of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() = 0;

    /// rhs of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> RHS() = 0;

    /// unknown displacements at t(n+1)
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp() = 0;

    /// unknown displacements at t(n)
    virtual Teuchos::RCP<const Epetra_Vector> Dispn() = 0;

    //@}

    //! @name Misc

    /// dof map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() = 0;

    /// domain map of system matrix (do we really need this?)
    virtual const Epetra_Map& DomainMap() = 0;

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() = 0;

    /// direct access to discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization() = 0;

    /// are there any algebraic constraints?
    virtual bool HaveConstraint() = 0;
    
    /// get constraint manager defined in the structure
    virtual const RCP<const UTILS::ConstrManager> GetConstraintManager() = 0;
    
    //! Return MapExtractor for Dirichlet boundary conditions 
    virtual const RCP<const LINALG::MapExtractor> GetDBCMapExtractor() = 0;

    /// structure part of the rhs to enforce coupling condition at the interface
    /*!
      structure part of
      \f$ \mathbf{u}^{n+1}_{\Gamma} \Delta t = \mathbf{d}^{n+1}_{\Gamma} - \mathbf{d}^{n}_{\Gamma} \f$
    */
    //virtual Teuchos::RCP<Epetra_Vector> FluidCondRHS() const = 0;


    /// structure part of the rhs to enforce structure -- mesh coupling at the interface
    /*!
      \f$ \mathbf{d}^{G,n+1}_{\Gamma} = \mathbf{d}^{n+1}_{\Gamma} \f$
    */
    //virtual Teuchos::RCP<Epetra_Vector> MeshCondRHS() const = 0;

    /// communication object at the interface
    virtual const LINALG::MapExtractor& Interface() const = 0;

    //@}

    //! @name Time step helpers
    //@{

    /// Get the current time
    virtual const double GetTime() const = 0;
    
    /// Get upper limit of time range of interest
    virtual const double GetTimeEnd() const = 0;

    /// Get time step size \f$\Delta t_n\f$
    virtual const double GetTimeStepSize() const = 0;

    /// Get current time step \f$n\f$ (sometimes called \f$k\f$ [but not in Baci])
    virtual const int GetTimeStep() const = 0;

    /// Get number of time steps
    virtual const int GetTimeNumStep() const = 0;

    /// Take the time and integrate (time loop)
    /// \date 11/08
    void Integrate();

    /// start new time step
    virtual void PrepareTimeStep() = 0;

    /// evaluate elements with given displacement
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disp) = 0;

    /// update at time step end
    virtual void Update() = 0;
    
    //! Update iteration
    //! Add residual increment to Lagrange multipliers stored in Constraint manager
    virtual void UpdateIterIncrConstr(Teuchos::RCP<Epetra_Vector> lagrincr) = 0;
 

    /// output results
    virtual void Output() = 0;

    /// read restart information for given time step
    virtual void ReadRestart(int step) = 0;

    //@}

    //! @name Solver calls

    /// nonlinear solve
    ///
    ///  Do the nonlinear solve for the time step. All boundary conditions have
    ///  been set.
    virtual void Solve() = 0;

    /// linear structure solve with just a interface load
    ///
    ///  The very special solve done in steepest descent relaxation
    ///  calculation (and matrix free Newton Krylov).
    ///
    /// \note Can only be called after a valid structural solve.
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> iforce) = 0;

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() = 0;

    //@}

    //! @name Extract interface values

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispn() = 0;

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispnp() = 0;

    /// Predictor for interface displacements
    virtual Teuchos::RCP<Epetra_Vector> PredictInterfaceDispnp() = 0;

    /// extract interface forces (necessary for robin-BC)
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces() = 0;

    //@}

    //! @name Apply interface forces

    /// apply interface forces to structural solver
    ///
    /// This prepares a new solve of the structural field within one time
    /// step. The middle values are newly created.
    ///
    /// \note This is not yet the most efficient implementation.
    virtual void ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce) = 0;

    /// apply robin enhanced interface forces to structural solver
    virtual void ApplyInterfaceRobinValue(Teuchos::RCP<Epetra_Vector> iforce,
                                          Teuchos::RCP<Epetra_Vector> ifluidvel) = 0;

    //@}

    /// create result test for encapulated structure algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest() = 0;
  };


  /// structure field solver
  class StructureBaseAlgorithm
  {
  public:

    /// constructor
    explicit StructureBaseAlgorithm(const Teuchos::ParameterList& prbdyn);

    /// virtual destructor to support polymorph destruction
    virtual ~StructureBaseAlgorithm();

    /// structural field solver
    Structure& StructureField() { return *structure_; }
    const Structure& StructureField() const { return *structure_; }

  private:

    /// setup structure algorithm
    void SetupStructure(const Teuchos::ParameterList& prbdyn);

    /// setup structure algorithm of StruGenAlpha type
    void SetupStruGenAlpha(const Teuchos::ParameterList& prbdyn);

    /// setup structure algorithm of STR::TimIntImpl type
    void SetupTimIntImpl(const Teuchos::ParameterList& prbdyn);

    /// structural field solver
    Teuchos::RCP<Structure> structure_;
  };

}

#endif
#endif

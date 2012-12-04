/*----------------------------------------------------------------------*/
/*!
\file ad_str_structure.H

\brief Structure field adapter

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef ADAPTER_STRUCTURE_H
#define ADAPTER_STRUCTURE_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_inpar/inpar_structure.H"
#include "../drt_lib/drt_dserror.H"


// forward declarations:
namespace DRT
{
  class Discretization;
  class ResultTest;
  namespace UTILS
  {
    class LocsysManager;
  }
}

namespace IO
{
  class DiscretizationWriter;
}

namespace LINALG
{
  class Solver;
  class SparseMatrix;
  class BlockSparseMatrixBase;
  class MapExtractor;
  class MultiMapExtractor;
}

namespace MORTAR
{
  class ManagerBase;
}

namespace UTILS
{
  class ConstrManager;
}


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

    /// DOF map of vector of unknowns for multiple dofsets
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds) = 0;

    /// domain map of system matrix (do we really need this?)
    virtual const Epetra_Map& DomainMap() = 0;

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() = 0;

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() = 0;

    /// switch structure field to block matrix
    virtual void UseBlockMatrix(Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps,
        Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps) = 0;

    /// recalculate structural matrices for tsi simulations
    virtual void TSIMatrix() = 0;

    /// return contact/meshtying manager
    virtual Teuchos::RCP<MORTAR::ManagerBase> ContactManager() = 0;

    // access to locsys manager
    virtual Teuchos::RCP<DRT::UTILS::LocsysManager> LocsysManager() = 0;

    /// direct access to discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization() = 0;

    /// are there any algebraic constraints?
    virtual bool HaveConstraint() = 0;

    /// get constraint manager defined in the structure
    virtual const Teuchos::RCP<UTILS::ConstrManager> GetConstraintManager() = 0;

    /// Get type of thickness scaling for thin shell structures
    virtual INPAR::STR::STC_Scale GetSTCAlgo() = 0;

    /// Access to scaling matrix for STC
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetSTCMat() = 0;

    /// Return MapExtractor for Dirichlet boundary conditions
    virtual const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() = 0;

    /*
    \brief Shorten the Dirichlet DOF set

    The method shortens the DOF set (map) which contains the DOFs
    subjected to Dirichlet boundary conditions. For instance, the method is
    called by the monolithic FSI where fluid boundary conditions on the
    FSI interface have to be removed if it is a fluid split.

    */
    virtual void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove) = 0;

    //@}

    //! @name Time step helpers
    //@{

    /// return time integration factor
    virtual double TimIntParam() = 0;

    /// Get the current time
    virtual double GetTime() const = 0;

    /// Get upper limit of time range of interest
    virtual double GetTimeEnd() const = 0;

    /// Get time step size \f$\Delta t_n\f$
    virtual double GetTimeStepSize() const = 0;

    /// Get current time step \f$n\f$ (sometimes called \f$k\f$ [but not in Baci])
    virtual int GetTimeStep() const = 0;

    /// Get number of time steps
    virtual int GetTimeNumStep() const = 0;

    /// Take the time and integrate (time loop)
    /// \date 11/08
    virtual void Integrate() = 0;

    /// integrate a single time step
    virtual void IntegrateStep() = 0;

    /// start new time step
    virtual void PrepareTimeStep() = 0;

    /// set time step size
    virtual void SetTimeStepSize(double timestepsize) = 0;

    /// set time step size and time
    virtual void SetInitialTimeStepAndTime(double timestepsize) = 0;

    /// update displacement and evaluate elements
    ///
    /// there are two displacement increments possible
    ///
    /// x^n+1_i+1 = x^n+1_i + disiterinc  (sometimes referred to as residual increment), and
    ///
    /// x^n+1_i+1 = x^n     + disstepinc
    ///
    /// with n and i being time and Newton iteration step
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> disiterinc ///< displacement increment between Newton iteration i and i+1
        ) = 0;

    /// update at time step end
    virtual void Update() = 0;

    /// Update iteration
    /// Add residual increment to Lagrange multipliers stored in Constraint manager
    virtual void UpdateIterIncrConstr(Teuchos::RCP<Epetra_Vector> lagrincr) = 0;

    /// Access to output object
    virtual Teuchos::RCP<IO::DiscretizationWriter> DiscWriter() = 0;

    /// prepare output (i.e. calculate stresses, strains, energies)
    virtual void PrepareOutput() = 0;

    /// output results
    virtual void Output() = 0;

    /// read restart information for given time step
    virtual void ReadRestart(int step) = 0;

    /// evaluate reference state (for frictional contact)
    virtual void EvaluateReferenceState() = 0;

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
//    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> iforce) = 0;

    //! Linear structure solve with just an interface load
    virtual Teuchos::RCP<Epetra_Vector> SolveRelaxationLinear() = 0;

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() = 0;

    //@}

    //! @name Extract methods
    //@{

    //! Extract coupling values (displacements)

    /// extract displacements, stress and velocities for TSI
    virtual Teuchos::RCP<Epetra_Vector> ExtractDispn() = 0;

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> ExtractDispnp() = 0;

    /// extract velocities
    virtual Teuchos::RCP<Epetra_Vector> ExtractVeln() = 0;

    /// extract velocities
    virtual Teuchos::RCP<Epetra_Vector> ExtractVelnp() = 0;

    /// extract velocities
    virtual Teuchos::RCP<Epetra_Vector> ExtractVelaf() = 0;

    /// extract velocities
    virtual Teuchos::RCP<Epetra_Vector> ExtractAccn() = 0;

    //@}

    //! @name TSI specific methods
    //@{

    //! Apply temperatures  for TSI

    /// apply temperatures to structural solver
    ///
    /// This prepares a new solve of the structural field within one time
    /// step. The middle values are newly created.
    ///
    /// \note This is not yet the most efficient implementation.
    virtual void ApplyTemperatures(Teuchos::RCP<const Epetra_Vector> temp) = 0;

    virtual void ApplyVelAndPress(Teuchos::RCP<const Epetra_Vector> velnp) = 0;

    /// Set forces due to interface with fluid, the force is expected external-force-like
    virtual void SetForceInterface(Teuchos::RCP<Epetra_Vector> iforce) = 0;

    //! specific method for iterative staggered partitioned TSI

    /// Identify residual
    /// This method does not predict the target solution but
    /// evaluates the residual and the stiffness matrix.
    /// In partitioned solution schemes, it is better to keep the current
    /// solution instead of evaluating the initial guess (as the predictor)
    /// does.
    virtual void PreparePartitionStep() = 0;

    //@}

    //! @name Structure with ale specific methods
    //@{

    /// material displacements (structure with ale)
    virtual Teuchos::RCP<Epetra_Vector> DispMat() = 0;

    /// apply material displacements to structure field (structure with ale)
    virtual void ApplyDisMat(Teuchos::RCP<Epetra_Vector> dismat) = 0;

    //@}

    /// create result test for encapsulated structure algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest() = 0;

    /// reset time and state vectors (needed for biofilm simulations)
    virtual void Reset() = 0;

  };


  /// structure field solver
  class StructureBaseAlgorithm
  {
  public:

    /// constructor
    explicit StructureBaseAlgorithm(const Teuchos::ParameterList& prbdyn, Teuchos::RCP<DRT::Discretization> actdis);

    /// virtual destructor to support polymorph destruction
    virtual ~StructureBaseAlgorithm();

    /// structural field solver
    Structure& StructureField() { return *structure_; }

    /// const version of structural field solver
    const Structure& StructureField() const { return *structure_; }

    /// Teuchos::rcp version of structural field solver
    Teuchos::RCP<Structure> StructureFieldrcp() { return structure_; }

  private:

    /// setup structure algorithm
    void SetupStructure(const Teuchos::ParameterList& prbdyn, Teuchos::RCP<DRT::Discretization> actdis);

    /// setup structure algorithm of STR::TimIntImpl type
    void SetupTimInt(const Teuchos::ParameterList& prbdyn, Teuchos::RCP<DRT::Discretization> actdis);

    /// create linear solver for contact/meshtying problems
    /** \param Teuchos::RCP<DRT::Discretization> the discretization class with the corresponding structural elements
     *
     * returns a contact solver object.
     *
     * Per default the CONTACT SOLVER block from the dat file is used for generating the solver object.
     * The idea is, that this linear solver object is used whenever there is contact between (two) structures.
     * Otherwise the standard structural solver block is used (generated by <tt>CreateLinearSolver</tt>. So we can
     * use highly optimized solvers for symmetric pure structural problems, but choose a different solver for the
     * hard nonsymmetric contact case. We automatically switch from CONTACT SOLVER (in case of contact) to STRUCT SOLVER
     * (pure structural problem, no contact) and back again.
     *
     * \note For contact/meshtying problems in the saddlepoint formulation (not condensed) this routines requires a block preconditioner
     * (eg <tt>CheapSIMPLE</tt> as preconditioner for the CONTACT SOLVER. The STRUCT SOLVER block is used for the primary
     * (structural) variables and the CONTACT CONSTRAINT SOLVER block is used for the Lagrange multiplier block.
     *
     * \note ATTENTION: condensed meshtying problems use the standard STRUCT SOLVER block (generated by <tt>CreateLinearSolver</tt>).
     * We assume that in contrary to contact problems the domain configuration is not changing for meshtying over the time.
     */
    Teuchos::RCP<LINALG::Solver> CreateContactMeshtyingSolver(Teuchos::RCP<DRT::Discretization>& actdis);

    /// create linear solver for pure structure problems
    /** returns linear solver object for pure structural problems
     *
     *  uses the parameters from the STRUCT SOLVER block in the BACI dat file. This solver is used for pure structural problems,
     *  whenever there is no contact. For contact/meshtying problems in saddlepoint formulation the STRUCT SOLVER parameters are
     *  also used for the block preconditioner (structural block).
     */
    Teuchos::RCP<LINALG::Solver> CreateLinearSolver(Teuchos::RCP<DRT::Discretization>& actdis);

    /// structural field solver
    Teuchos::RCP<Structure> structure_;
  };

}

#endif

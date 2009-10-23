/*----------------------------------------------------------------------*/
/*!
\file adapter_thermo.H

\brief Thermo field adapter

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253

</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ADAPTER_THERMO_H
#define ADAPTER_THERMO_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_lib/drt_utils.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_thermo/thr_resulttest.H"


namespace ADAPTER {

  /// general thermal field interface
  /*!

  The point is to keep TSI as far apart from our field solvers as
  possible. Each thermal field solver we want to use should get its own
  subclass of this. The TSI algorithm should be able to extract all the
  information from the thermal field it needs using this interface.

  All TSI algorithms use this adapter to communicate with the thermal
  field. There are different ways to use this adapter.

  In all cases you need to tell the thermal algorithm about your time
  step. Therefore PrepareTimeStep(), Update() and Output() must be called at
  the appropriate position in the TSI algorithm.

  <h3>Dirichlet-Neumann coupled TSI</h3>

  A good starting temperature can be guessed with PredictInterfaceTempnp().

  Dirichlet-Neumann coupled TSI will need to Solve() the nonlinear
  thermal problem for each time step after the structure forces have been
  applied (ApplyInterfaceForces()). Solve() will be called many times for each
  time step until the interface equilibrium is reached. The thermal
  algorithm has to preserve its state until Update() is called.

  After each Solve() you get the interface forces by ExtractInterfaceTempnp().

  <h3>Monolithic TSI</h3>

  Monolithic TSI is based on Evaluate() of elements. This results in a new
  RHS() and a new SysMat(). Together with the InitialGuess() these form the
  building blocks for a block based Newton's method.

  \warning Further cleanup is still needed.

  \sa Fluid, Ale, Structure
  \author cd, bborn
  \date 07/09
  */
  class Thermo
  {
  public:

    /// virtual to get polymorph destruction
    virtual ~Thermo();

    /// @name Vector access
    //@{

    /// initial guess of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() = 0;

    /// RHS of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> RHS() = 0;

    /// unknown temperatures at t(n+1)
    virtual Teuchos::RCP<const Epetra_Vector> Tempnp() = 0;

    /// unknown temperatures at t(n)
    virtual Teuchos::RCP<const Epetra_Vector> Tempn() = 0;

    //@}

    /// @name Misc
    //@{

    /// DOF map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() = 0;

    /// domain map of system matrix (do we really need this?)
    virtual const Epetra_Map& DomainMap() = 0;

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() = 0;

    /// direct access to discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization() = 0;

    /// Return MapExtractor for Dirichlet boundary conditions
    virtual Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() = 0;

    //@}

    //! @name Time step helpers
    //@{

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
    void Integrate();

    /// start new time step
    virtual void PrepareTimeStep() = 0;

    /// evaluate residual at given temperature (?)
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> temp) = 0;

    /// update at time step end
    virtual void Update() = 0;

    /// output results
    virtual void Output() = 0;

    /// read restart information for given time step
    virtual void ReadRestart(const int step) = 0;

    //@}

    //! @name Solver calls
    //@{

    /// nonlinear solve
    ///
    ///  Do the nonlinear solve for the time step. All boundary conditions have
    ///  been set.
    virtual void Solve() = 0;

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() = 0;

    //@}

    /// create result test for encapulated thermo algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest() = 0;
  };


  /// thermo field solver
  class ThermoBaseAlgorithm
  {
  public:

    /// constructor
    explicit ThermoBaseAlgorithm(const Teuchos::ParameterList& prbdyn);

    /// virtual destructor to support polymorph destruction
    virtual ~ThermoBaseAlgorithm();

    /// thermal field solver
    Thermo& ThermoField() { return *thermo_; }
    /// const version of thermal field solver
    const Thermo& ThermoField() const { return *thermo_; }

  private:

    /// setup thermo algorithm
    void SetupThermo(const Teuchos::ParameterList& prbdyn);

    /// setup thermo algorithm of THR::TimIntImpl type
    void SetupTimIntImpl(const Teuchos::ParameterList& prbdyn);

    /// thermal field solver
    Teuchos::RCP<Thermo> thermo_;

  };  // class ThermoBaseAlgorithm

}  // namespace ADAPTER

#endif // ADAPTER_THERMO_H
#endif // CCADISCRET

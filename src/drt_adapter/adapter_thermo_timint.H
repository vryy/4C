/*----------------------------------------------------------------------*/
/*!
\file adapter_thermo_timint.H

\brief Thermal time integration field adapter

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237

            changed to thermo 07.09 bb, cd

</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

#ifndef ADAPTER_THERMO_TIMINT_H
#define ADAPTER_THERMO_TIMINT_H

/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_lib/drt_utils.H"
//#include "../drt_thermo/thr_resulttest.H"

#include "../drt_thermo/thrtimint_impl.H"
#include "adapter_thermo.H"


/*----------------------------------------------------------------------*/
/* adapting adapter */
namespace ADAPTER {

  /*====================================================================*/
  /*!
   * \brief Adapter to general implicit thermal time integration
   *
   * \date 08/08
   */
  class ThermoTimInt : public Thermo
  {
  public:

    /// Constructor
    ThermoTimInt(
      Teuchos::RCP<Teuchos::ParameterList> ioparams,  ///< I/O flags
      Teuchos::RCP<Teuchos::ParameterList> tdynparams,  ///< input parameters
      Teuchos::RCP<Teuchos::ParameterList> xparams,  ///< extra flags
      Teuchos::RCP<DRT::Discretization> discret,  ///< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,  ///< the solver
      Teuchos::RCP<IO::DiscretizationWriter> output  ///< the output
    );

    /// initial guess of Newton's method
    Teuchos::RCP<const Epetra_Vector> InitialGuess();

    /// right-hand-side of Newton's method
    Teuchos::RCP<const Epetra_Vector> RHS();

    /// unknown temperature at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> Tempnp();

    /// unknown temperature at \f$t_{n}\f$
    Teuchos::RCP<const Epetra_Vector> Tempn();

    /// dof map of vector of unknowns
    Teuchos::RCP<const Epetra_Map> DofRowMap();

    /// domain map of system matrix (do we really need this?)
    const Epetra_Map& DomainMap();

    /// direct access to system matrix
    Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

    /// direct access to discretization
    Teuchos::RCP<DRT::Discretization> Discretization();

    /// Return MapExtractor for Dirichlet boundary conditions
    virtual Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
    {
      return thermo_->GetDBCMapExtractor();
    }

    /// External force at t_{n+1} (corresponds in structure to FExtn())
    Teuchos::RCP<const Epetra_Vector> FExtnp();

    /// @name Time step helpers
    //@{

    /// Get the current time
    virtual double GetTime() const { return thermo_->GetTime(); }

    /// Get upper limit of time range of interest
    virtual double GetTimeEnd() const { return thermo_->GetTimeEnd(); }

    /// Get time step size \f$\Delta t_n\f$
    virtual double GetTimeStepSize() const { return thermo_->GetTimeStepSize(); }

    /// Get current time step \f$n\f$ (sometimes called \f$k\f$ [but not in Baci])
    virtual int GetTimeStep() const { return thermo_->GetStep(); }

    /// Get number of time steps
    virtual int GetTimeNumStep() const { return thermo_->GetTimeNumStep(); }

    /// start new time step
    void PrepareTimeStep();

    /// evaluate elements with given temperature
    void Evaluate(
      Teuchos::RCP<const Epetra_Vector> temp  ///< current state
    );

    /// update at time step end
    void Update();

    /// output results
    void Output();

    //@}

    /// read restart information for given time step
    virtual void ReadRestart
    (
      const int step  ///< given time step
    );

    /// @name Solver calls
    //@{

    /// non-linear solve
    ///
    /// Do the nonlinear solve for the time step. All boundary conditions have
    /// been set.
    virtual void Solve();

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

    //@}

    /// create result test for encapsulated thermo algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  private:

    /// create implicit marching time integrator
    Teuchos::RCP<THR::TimIntImpl> Create(
      const Teuchos::ParameterList& ioflags,  //!< input-output-flags
      const Teuchos::ParameterList& tdyn,  //!< thermal dynamic flags
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization>& actdis,  //!< discretisation
      Teuchos::RCP<LINALG::Solver>& solver,  //!< the solver
      Teuchos::RCP<IO::DiscretizationWriter>& output  //!< output writer
    );

    /// the actual thermo algorithm
    Teuchos::RCP<THR::TimIntImpl> thermo_;

    /// sum of temperature vectors already applied,
    /// i.e. the incremental temperature
    Teuchos::RCP<Epetra_Vector> tempinc_;

    /// @name local copies of input parameters
    //{@
    Teuchos::RCP<DRT::Discretization> discret_;  ///< the discretisation
    Teuchos::RCP<Teuchos::ParameterList> ioparams_;  ///< I/O flags ... not sure of really needed
    Teuchos::RCP<Teuchos::ParameterList> tdynparams_;  ///< dynamic control flags ... used,
                                                       ///< but could/should be circumvented
    Teuchos::RCP<Teuchos::ParameterList> xparams_;  ///< eXtra input parameters
    Teuchos::RCP<LINALG::Solver> solver_;  ///< the solver
    Teuchos::RCP<IO::DiscretizationWriter> output_;  ///< the output writer
    //@}

  };  // class ThermoTimInt

}  // namespace ADAPTER

/*----------------------------------------------------------------------*/
#endif
#endif

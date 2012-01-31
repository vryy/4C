/*----------------------------------------------------------------------*/
/*!
\file adapter_thermo_timint_expl.H

\brief Thermal time integration field adapter

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253

</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                              dano 01/12 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ADAPTER_THERMO_TIMINT_EXPL_H
#define ADAPTER_THERMO_TIMINT_EXPL_H

/*----------------------------------------------------------------------*
 |  headers                                                  dano 01/12 |
 *----------------------------------------------------------------------*/
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_lib/drt_utils.H"

#include "../drt_thermo/thrtimint_expl.H"
#include "adapter_thermo.H"

/*----------------------------------------------------------------------*
 |  adapting adapter                                         dano 01/12 |
 *----------------------------------------------------------------------*/
namespace ADAPTER
{

  /*====================================================================*/
  /*!
   * \brief Adapter to general explicit thermal time integration
   *
   * \date 08/08
   */
  class ThermoTimIntExpl : public Thermo
  {
  public:

    /// Constructor
    ThermoTimIntExpl(
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

    /// dof map of vector of unknowns for multiple dofsets
    Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds);

    /// domain map of system matrix (do we really need this?)
    const Epetra_Map& DomainMap();

    /// direct access to system matrix
    Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

    /// recalculate thermal matrices for tsi simulations
    void TSIMatrix();

    /// direct access to discretization
    Teuchos::RCP<DRT::Discretization> Discretization();

    /// Return MapExtractor for Dirichlet boundary conditions
    virtual Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
    {
      return thermo_->GetDBCMapExtractor();
    }

    /// External force at t_{n+1} (corresponds in structure to FExtn())
    Teuchos::RCP<const Epetra_Vector> FExtnp();

    // access to thermal contact manager
    Teuchos::RCP<THR::ThermoContactMan> ThermoContactManager();

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

    /// Set time step size for the current step
    virtual void SetTimeStepSize(double timestepsize){ thermo_->SetTimeStepSize(timestepsize); }

    /// Set time step size for the first step
    virtual void SetInitialTimeStepAndTime(double timestepsize){ thermo_->SetInitialTimeStepAndTime(timestepsize); }

    /// start new time step
    void PrepareTimeStep();

    /// set the interface map extractor in the thermo field
    void SetSurfaceTFSI(
      Teuchos::RCP<const LINALG::MapExtractor> tfsisurface  //!< the TFSI surface
      );

    /// apply external interface loads to the thermo field
    void ApplyInterfaceForces(
      Teuchos::RCP<Epetra_Vector> ithermoload   ///< load due to the tfsi interface
      );

    /// evaluate element with given temperature increment
    void Evaluate(
      Teuchos::RCP<const Epetra_Vector> temp  ///< current state
      );

    /// evaluate element with given temperature increment
    void Evaluate();

    /// update temperature increment after Newton step
    void UpdateNewton(
      Teuchos::RCP<const Epetra_Vector> temp  ///< current temperature increment
      );

    /// update at time step end
    void Update();

    /// print info about finished time step
    void PrintStep();

    //! Access to output object
    Teuchos::RCP<IO::DiscretizationWriter> DiscWriter()
    {
      return output_;
    }

    /// output results
    void Output();

    //@}

    /// read restart information for given time step
    virtual void ReadRestart(
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

    //! @name Extract current temperature values needed for TSI
    //@{
    /// extract temperatures for inserting in structure field
    virtual Teuchos::RCP<Epetra_Vector> ExtractTempn();

    /// extract current temperatures for inserting in structure field
    virtual Teuchos::RCP<Epetra_Vector> ExtractTempnp();

    /// prepare thermal contact
    void PrepareThermoContact(Teuchos::RCP<MORTAR::ManagerBase> cmtman,
                              Teuchos::RCP<DRT::Discretization> discretstruct);

    //@}

    //! @name Apply current displacements and velocities (for TSI)
    //@{

    /// apply current displacements and velocities thermo field
    void ApplyStructVariables(
      Teuchos::RCP<const Epetra_Vector> disp,
      Teuchos::RCP<const Epetra_Vector> vel
      );

    //@}

    /// Identify residual
    /// This method does not predict the target solution but
    /// evaluates the residual and the stiffness matrix.
    /// In partitioned solution schemes, it is better to keep the current
    /// solution instead of evaluating the initial guess (as the predictor)
    /// does.
    void PreparePartitionStep();

    /// create result test for encapsulated thermo algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  private:

    /// create implicit marching time integrator
    Teuchos::RCP<THR::TimIntExpl> Create(
      const Teuchos::ParameterList& ioflags,  //!< input-output-flags
      const Teuchos::ParameterList& tdyn,  //!< thermal dynamic flags
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization>& actdis,  //!< discretisation
      Teuchos::RCP<LINALG::Solver>& solver,  //!< the solver
      Teuchos::RCP<IO::DiscretizationWriter>& output  //!< output writer
      );

    /// the actual explicit thermo algorithm
    Teuchos::RCP<THR::TimIntExpl> thermo_;

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

  };  // class ThermoTimIntExpl

}  // namespace ADAPTER

/*----------------------------------------------------------------------*/
#endif // ADAPTER_THERMO_TIMINT_EXPL_H
#endif // CCADISCRET

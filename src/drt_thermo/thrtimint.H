/*----------------------------------------------------------------------*/
/*! \file
\brief base class for all time integrators in thermal field
\level 1
\maintainer Sebastian Proell
*/


/*----------------------------------------------------------------------*
 | definitions                                              bborn 06/08 |
 *----------------------------------------------------------------------*/
#ifndef THRTIMINT_H
#define THRTIMINT_H


/*----------------------------------------------------------------------*
 | headers                                                  bborn 06/08 |
 *----------------------------------------------------------------------*/

#include "../drt_adapter/adapter_thermo.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_timestepping/timintmstep.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../drt_inpar/drt_validparameters.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"

// forward declaration
namespace CONTACT
{
  class CoNitscheStrategyTsi;
}

/*----------------------------------------------------------------------*
 | belongs to thermal dynamics namespace                    bborn 08/09 |
 *----------------------------------------------------------------------*/
namespace THR
{
  /*====================================================================*/
  /*!
   * \brief Front-end for thermal dynamics by integrating in time.
   *
   * <h3> Intention </h3>
   * This front-end for thermal dynamics defines an interface to call
   * several derived time integrators. Thus it describes a plethora of pure
   * virtual methods which have to be implemented at the derived integrators.
   * However, it also offers a few non-empty methods and stores associated
   * data. The most important method of this base time integrator object
   * is #Integrate().
   *
   * #Integrate() performs a time integration (time loop) with constant
   * time steps and other parameters as set by the user.
   *
   * Although #Integrate is the main interface, this base time integrator
   * allows the public to access a few of its datum objects, for instance
   * the tangent system matrix #tang_ by #Tang(). This selective access
   * is needed in environments in which a independent time loop is provided.
   * This happens e.g. in fluid-structure-interaction.
   *
   * <h3> Responsibilties </h3>
   * Most importantly the base integrator manages the system state vectors and
   * matrices. It also deals with the output to files and offers method to
   * determine forces and tangents.
   *
   * \author bborn
   * \date 06/08
   */
  class TimInt : public ::ADAPTER::Thermo
  {
   public:
    //! @name Life
    //@{

    //! Print thermo time logo
    void Logo();

    //! Constructor
    TimInt(const Teuchos::ParameterList& ioparams,     //!< ioflags
        const Teuchos::ParameterList& tdynparams,      //!< input parameters
        const Teuchos::ParameterList& xparams,         //!< extra flags
        Teuchos::RCP<DRT::Discretization> actdis,      //!< current discretisation
        Teuchos::RCP<LINALG::Solver> solver,           //!< the solver
        Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Empty constructor
    TimInt() { ; }

    //! Copy constructor
    TimInt(const TimInt& old) { ; }

    //! Destructor
    virtual ~TimInt() { DetachEnergyFile(); }

    //! Resize #TimIntMStep<T> multi-step quantities
    virtual void ResizeMStep() = 0;

    //@}

    //! @name Actions
    //@{

    //! Equilibrate the initial state by identifying the consistent
    //! initial accelerations and (if applicable) internal variables
    //! Make capacity matrix
    void DetermineCapaConsistTempRate();

    //! Apply Dirichlet boundary conditions on provided state vectors
    void ApplyDirichletBC(const double time,  //!< at time
        Teuchos::RCP<Epetra_Vector> temp,     //!< temperatures
                                              //!< (may be Teuchos::null)
        Teuchos::RCP<Epetra_Vector> rate,     //!< temperature rate
                                              //!< (may be Teuchos::null)
        bool recreatemap                      //!< recreate mapextractor/toggle-vector
                                              //!< which stores the DOF IDs subjected
                                              //!< to Dirichlet BCs
                                              //!< This needs to be true if the bounded DOFs
                                              //!< have been changed.
    );

    //! prepare thermal contact
    void SetNitscheContactStrategy(Teuchos::RCP<CONTACT::CoNitscheStrategyTsi> strategy)
    {
      contact_strategy_nitsche_ = strategy;
    }

    //! prepare time step
    virtual void PrepareTimeStep() = 0;

    //! Do time integration of single step
    virtual void IntegrateStep() = 0;

    /// tests if there are more time steps to do
    virtual bool NotFinished() const
    {
      return timen_ <= timemax_ + 1.0e-8 * (*dt_)[0] and stepn_ <= stepmax_;
    }

    //! non-linear solve
    //!
    //! Do the nonlinear solve, i.e. (multiple) corrector,
    //! for the time step. All boundary conditions have
    //! been set.
    virtual INPAR::THR::ConvergenceStatus Solve() = 0;

    //! build linear system tangent matrix, rhs/force residual
    //! Monolithic TSI accesses the linearised thermo problem
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> tempi) = 0;

    //! build linear system tangent matrix, rhs/force residual
    //! Monolithic TSI accesses the linearised thermo problem
    virtual void Evaluate() = 0;

    //! Update configuration after time step
    //!
    //! Thus the 'last' converged is lost and a reset of the time step
    //! becomes impossible. We are ready and keen awating the next
    //! time step.
    virtual void UpdateStepState() = 0;

    //! Update everything on element level after time step and after output
    //!
    //! Thus the 'last' converged is lost and a reset of the time step
    //! becomes impossible. We are ready and keen awating the next time step.
    virtual void UpdateStepElement() = 0;

    //! Update time and step counter
    void UpdateStepTime();

    //! update at time step end
    virtual void Update() = 0;

    //! update Newton step
    virtual void UpdateNewton(Teuchos::RCP<const Epetra_Vector> tempi) = 0;

    //! Reset configuration after time step
    //!
    //! Thus the last converged state is copied back on the predictor
    //! for current time step. This applies only to elemet-wise
    //! quantities
    virtual void ResetStep();

    //! set the initial thermal field
    void SetInitialField(const INPAR::THR::InitialField,  //!< type of initial field
        const int startfuncno                             //!< number of spatial function
    );

    //@}

    //! @name Output
    //@{

    //! print summary after step
    virtual void PrintStep() = 0;

    //! Output to file
    //! This routine always prints the last converged state, i.e.
    //! \f$T_{n}, R_{n}\f$. So, #UpdateIncrement should be called
    //! upon object prior to writing stuff here.
    //! \author mwgee (originally) \date 03/07
    void OutputStep(bool forced_writerestart);

    //! output
    virtual void Output(bool forced_writerestart = false) { OutputStep(forced_writerestart); }

    //! Write restart
    //! \author mwgee (originally) \date 03/07
    void OutputRestart(bool& datawritten  //!< (in/out) read and append if
                                          //!< it was written at this time step
    );

    //! Output temperatures, temperature rates
    //! and more system vectors
    //! \author mwgee (originally) \date 03/07
    void OutputState(bool& datawritten  //!< (in/out) read and append if
                                        //!< it was written at this time step
    );

    //! Add restart information to OutputState
    void AddRestartToOutputState();

    //! Heatflux & temperature gradient output
    //! \author lw (originally)
    void OutputHeatfluxTempgrad(bool& datawritten  //!< (in/out) read and append if
                                                   //!< it was written at this time step
    );

    //! Energy output
    void OutputEnergy();

    //! phase information in SLM
    void OutputPhase(bool& datawritten);

    //! postproc output of material parameters
    void OutputConductivityAndCapacity(bool& datawritten);

    //! Write internal and external forces (if necessary for restart)
    virtual void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output) = 0;

    //! Check wether energy output file is attached
    bool AttachedEnergyFile()
    {
      if (energyfile_)
        return true;
      else
        return false;
    }

    //! Attach file handle for energy file #energyfile_
    void AttachEnergyFile()
    {
      if (not energyfile_)
      {
        std::string energyname =
            DRT::Problem::Instance()->OutputControlFile()->FileName() + ".thermo.energy";
        energyfile_ = new std::ofstream(energyname.c_str());
        *energyfile_ << "# timestep time internal_energy" << std::endl;
      }
    }

    //! Detach file handle for energy file #energyfile_
    void DetachEnergyFile()
    {
      if (energyfile_) delete energyfile_;
    }

    //! Identify residual
    //! This method does not predict the target solution but
    //! evaluates the residual and the stiffness matrix.
    //! In partitioned solution schemes, it is better to keep the current
    //! solution instead of evaluating the initial guess (as the predictor)
    //! does.
    virtual void PreparePartitionStep() = 0;

    //! thermal result test
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

    //@}

    //! @name Forces and Tangent
    //@{

    //! Apply external force
    void ApplyForceExternal(const double time,   //!< evaluation time
        const Teuchos::RCP<Epetra_Vector> temp,  //!< temperature state
        Teuchos::RCP<Epetra_Vector>& fext        //!< external force
    );

    //! Apply convective boundary conditions force
    void ApplyForceExternalConv(Teuchos::ParameterList& p,
        const double time,                         //!< evaluation time
        const Teuchos::RCP<Epetra_Vector> tempn,   //!< old temperature state T_n
        const Teuchos::RCP<Epetra_Vector> temp,    //!< temperature state T_n+1
        Teuchos::RCP<Epetra_Vector>& fext,         //!< external force
        Teuchos::RCP<LINALG::SparseOperator> tang  //!< tangent at time n+1
    );

    //! Evaluate ordinary internal force, its tangent at state
    void ApplyForceTangInternal(
        Teuchos::ParameterList& p,                //!< parameter list handed down to elements
        const double time,                        //!< evaluation time
        const double dt,                          //!< step size
        const Teuchos::RCP<Epetra_Vector> temp,   //!< temperature state
        const Teuchos::RCP<Epetra_Vector> tempi,  //!< residual temperatures
        Teuchos::RCP<Epetra_Vector> fint,         //!< internal force
        Teuchos::RCP<LINALG::SparseMatrix> tang   //!< tangent matrix
    );

    //! Evaluate ordinary internal force, its tangent and the stored force at state
    void ApplyForceTangInternal(
        Teuchos::ParameterList& p,                //!< parameter list handed down to elements
        const double time,                        //!< evaluation time
        const double dt,                          //!< step size
        const Teuchos::RCP<Epetra_Vector> temp,   //!< temperature state
        const Teuchos::RCP<Epetra_Vector> tempi,  //!< residual temperatures
        Teuchos::RCP<Epetra_Vector> fcap,         //!< capacity force
        Teuchos::RCP<Epetra_Vector> fint,         //!< internal force
        Teuchos::RCP<LINALG::SparseMatrix> tang   //!< tangent matrix
    );

    //! Evaluate ordinary internal force
    //!
    //! We need incremental temperatures, because the internal variables,
    //! chiefly EAS parameters with an algebraic constraint, are treated
    //! as well. They are not treated perfectly, ie they are not iteratively
    //! equilibriated according to their (non-linear) constraint and
    //! the pre-determined temperatures -- we talk explicit time integration
    //! here, but they are applied in linearised manner. The linearised
    //! manner means the static condensation is applied once with
    //! residual temperatures replaced by the full-step temperature
    //! increment \f$T_{n+1}-T_{n}\f$.
    void ApplyForceInternal(Teuchos::ParameterList& p,  //!< parameter list handed down to elements
        const double time,                              //!< evaluation time
        const double dt,                                //!< step size
        const Teuchos::RCP<Epetra_Vector> temp,         //!< temperature state
        const Teuchos::RCP<Epetra_Vector> tempi,        //!< incremental temperatures
        Teuchos::RCP<Epetra_Vector> fint                //!< internal force
    );

    //@}

    //! @name Thermo-structure-interaction specific methods
    //@{

    //! Set external loads (heat flux) due to tfsi interface
    virtual void SetForceInterface(
        Teuchos::RCP<Epetra_Vector> ithermoload  //!< thermal interface load
    );

    //@}

    //! @name Attributes
    //@{

    //! Provide Name
    virtual enum INPAR::THR::DynamicType MethodName() const = 0;

    //! Provide title
    std::string MethodTitle() const { return INPAR::THR::DynamicTypeString(MethodName()); }

    //! Return true, if time integrator is implicit
    virtual bool MethodImplicit() = 0;

    //! Return true, if time integrator is explicit
    bool MethodExplicit() { return (not MethodImplicit()); }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a \f$m\f$-multistep method returns \f$m\f$
    virtual int MethodSteps() = 0;

    //! Give order of accuracy
    virtual int MethodOrderOfAccuracy() = 0;

    //! Return linear error coefficient of temperatures
    virtual double MethodLinErrCoeff() = 0;

    //@}

    //! @name Access methods
    //@{

    //! Access dicretisation
    Teuchos::RCP<DRT::Discretization> Discretization() { return discret_; }

    //! non-overlapping DOF map for multiple dofsets
    Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds)
    {
      const Epetra_Map* dofrowmap = discret_->DofRowMap(nds);
      return Teuchos::rcp(new Epetra_Map(*dofrowmap));
    }

    //! non-overlapping DOF map
    Teuchos::RCP<const Epetra_Map> DofRowMap()
    {
      const Epetra_Map* dofrowmap = discret_->DofRowMap();
      return Teuchos::rcp(new Epetra_Map(*dofrowmap));
    }

    //! Access solver
    Teuchos::RCP<LINALG::Solver> Solver() { return solver_; }

    //! get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

    //! Access output object
    Teuchos::RCP<IO::DiscretizationWriter> DiscWriter() { return output_; }

    //! prepare output (do nothing)
    virtual void PrepareOutput() { ; }

    //! Read restart values
    void ReadRestart(const int step  //!< restart step
    );

    //! Read and set restart state
    void ReadRestartState();

    //! Read and set restart forces
    virtual void ReadRestartForce() = 0;

    //! Return temperatures \f$T_{n}\f$
    Teuchos::RCP<Epetra_Vector> Temp() { return (*temp_)(0); }

    //! Return temperatures \f$T_{n}\f$
    Teuchos::RCP<Epetra_Vector> WriteAccessTempn() { return (*temp_)(0); }

    //! Return temperatures \f$T_{n}\f$
    Teuchos::RCP<const Epetra_Vector> Tempn() { return (*temp_)(0); }

    //! initial guess of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() = 0;

    //! Return temperatures \f$T_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> WriteAccessTempnp() { return tempn_; }

    //! Return temperatures \f$T_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> Tempnp() { return tempn_; }

    //! Return temperature rates \f$R_{n}\f$
    Teuchos::RCP<Epetra_Vector> Rate() { return (*rate_)(0); }

    //! Return temperature rates \f$R_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> RateNew() { return raten_; }

    //! Return external force \f$F_{ext,n}\f$
    virtual Teuchos::RCP<Epetra_Vector> Fext() = 0;

    //! Return reaction forces
    virtual Teuchos::RCP<Epetra_Vector> Freact() = 0;

    //! right-hand side alias the dynamic thermal residual
    virtual Teuchos::RCP<const Epetra_Vector> RHS() = 0;

    //! Return tangent, i.e. thermal residual differentiated by temperatures
    //! (SystemMatrix()/stiff_ in STR)
    Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return tang_; }

    //! Return domain map
    const Epetra_Map& GetDomainMap() { return tang_->DomainMap(); }

    //! Return domain map
    const Epetra_Map& DomainMap() { return tang_->DomainMap(); }

    //! Return current time \f$t_{n}\f$
    virtual double TimeOld() const { return (*time_)[0]; }

    //! Return target time \f$t_{n+1}\f$
    virtual double Time() const { return timen_; }

    //! Get upper limit of time range of interest
    virtual double GetTimeEnd() const { return timemax_; }

    //! Get time step size \f$\Delta t_n\f$
    virtual double Dt() const { return (*dt_)[0]; }

    //! Set time step size \f$\Delta t_n\f$
    virtual void SetDt(double timestepsize) { (*dt_)[0] = timestepsize; }

    //! Sets the target time \f$t_{n+1}\f$ of this time step
    virtual void SetTimen(const double time) { timen_ = time; }

    //! Return current step number $n$
    virtual int StepOld() const { return step_; }

    //! Return current step number $n+1$
    virtual int Step() const { return stepn_; }

    //! Get number of time steps
    virtual int NumStep() const { return stepmax_; }

    //! Update number of time steps (in adaptivity)
    virtual void SetNumStep(const int newNumStep) { stepmax_ = newNumStep; }

    //! Get communicator
    virtual inline const Epetra_Comm& Comm() const { return discret_->Comm(); }

    //! Return MapExtractor for Dirichlet boundary conditions
    const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() const { return dbcmaps_; }

    //! Return MapExtractor for Dirichlet boundary conditions
    Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() { return dbcmaps_; }

    //! Return thermal contact manager
    Teuchos::RCP<CONTACT::CoNitscheStrategyTsi> NitscheContactStrategy()
    {
      return contact_strategy_nitsche_;
    }

    //@}

   protected:
    //! evaluate error compared to analytical solution
    Teuchos::RCP<std::vector<double>> EvaluateErrorComparedToAnalyticalSol();

    //! @name General purpose algorithm members
    //@{

    Teuchos::RCP<DRT::Discretization> discret_;        //!< attached discretisation
    Teuchos::RCP<DRT::Discretization> discretstruct_;  //!< structural discretisation

    int myrank_;                                  //!< ID of actual processor in parallel
    Teuchos::RCP<LINALG::Solver> solver_;         //!< linear algebraic solver
    bool solveradapttol_;                         //!< adapt solver tolerance
    double solveradaptolbetter_;                  //!< tolerance to which is adpated ????
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;  //!< map extractor object
                                                  //!< containing non-overlapping
                                                  //!< map of global DOFs on Dirichlet
                                                  //!< boundary conditions
    //@}

    //! @name Printing and output
    //@{

    Teuchos::RCP<IO::DiscretizationWriter> output_;  //!< binary output
    bool printlogo_;                                 //!< true: enjoy your cuppa
    int printscreen_;                                //!< print infos to standard out every n steps
    FILE* errfile_;                                  //!< error file handle
    bool printerrfile_;                              //!< print infos to error file #errfile_
    bool printiter_;         //!< print intermediate iterations during solution
    int writerestartevery_;  //!< write restart every given step;
                             //!< if 0, restart is not written
    bool writeglob_;         //!< write state on/off
    int writeglobevery_;     //!< write state every given step
    INPAR::THR::HeatFluxType writeheatflux_;
    INPAR::THR::TempGradType writetempgrad_;
    INPAR::THR::PhaseType writephase_;
    int writeenergyevery_;             //!< write system energy every given step
    std::ofstream* energyfile_;        //!< outputfile for energy
    INPAR::THR::CalcError calcerror_;  //!< evaluate error compared to analytical solution
    int errorfunctno_;  //!< function number of analytical solution for error evaluation
    //@}

    //! @name General control parameters
    //@{

    Teuchos::RCP<TIMINT::TimIntMStep<double>> time_;  //!< time \f$t_{n}\f$ of last converged step
    double timen_;                                    //!< target time \f$t_{n+1}\f$
    Teuchos::RCP<TIMINT::TimIntMStep<double>> dt_;    //!< time step size \f$\Delta t\f$
    double timemax_;                                  //!< final time \f$t_\text{fin}\f$
    int stepmax_;                                     //!< final step \f$N\f$
    int step_;                                        //!< time step index \f$n\f$
    int stepn_;                                       //!< time step index \f$n+1\f$
    bool firstoutputofrun_;  //!< flag whether this output step is the first one (restarted or not)
    bool lumpcapa_;          //!< flag for lumping the capacity matrix, default: false

    //@}

    //! @name Global vectors
    //@{

    Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length

    //@}


    //! @name Global state vectors
    //@{

    //! global temperatures \f${T}_{n}, T_{n-1}, ...\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector>> temp_;
    //! global temperature rates \f${R}_{n}, R_{n-1}, ...\f$
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector>> rate_;
    Teuchos::RCP<Epetra_Vector> tempn_;  //!< global temperatures
                                         //!< \f${T}_{n+1}\f$
                                         //!< at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> raten_;  //!< global temperature rates
                                         //!< \f${R}_{n+1}\f$
                                         //!< at \f$t_{n+1}\f$
    //@}

    //! @name Interface stuff
    //@{

    Teuchos::RCP<Epetra_Vector> fifc_;  //!< external interface loads

    //@}

    //! @name System matrices
    //@{

    //! holds eventually effective tangent (STR: stiff_)
    Teuchos::RCP<LINALG::SparseMatrix> tang_;
    //! capacity matrix (constant)
    // Teuchos::RCP<LINALG::SparseMatrix> capa_;

    //@}

    //! @name Nitsche contact stuff
    //@{

    // thermo contact manager
    Teuchos::RCP<CONTACT::CoNitscheStrategyTsi> contact_strategy_nitsche_;

    //@}

  };  // TimInt

}  // namespace THR

/*----------------------------------------------------------------------*/
#endif  // #ifndef THRTIMINT_H

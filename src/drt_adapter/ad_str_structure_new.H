/*-----------------------------------------------------------*/
/*!
\file ad_str_structure_new.H

\brief Adapter for the new structural time integration framework.

\maintainer Michael Hiermeier

\date Sep 2, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef AD_STR_STRUCTURE_NEW_H_
#define AD_STR_STRUCTURE_NEW_H_

#include "../drt_inpar/inpar_structure.H"
#include "ad_str_structure.H"


namespace Teuchos{
  class ParameterList;
} // namespace Teuchos
namespace STR{
namespace TIMINT{
  class Base;
} // namespace TIMINT
} // namespace STR
namespace LINALG{
  class Solver;
} // namespace LINALG
namespace DRT{
  class Discretization;
} // namespace DRT

namespace ADAPTER {

  class StructureNew : public Structure
  {
  public:
    StructureNew() {};

    /// virtual to get polymorph destruction
    virtual ~StructureNew() {};

    /// @name Vector access
    ///@{
    /// initial guess of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() = 0;

    /// rhs of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> RHS() = 0;

    /// unknown displacements at \f$t_{n+1}\f$
    /// ToDo Replace the deprecated version with the new version
    virtual Teuchos::RCP<const Epetra_Vector> DispNp() const = 0;
    Teuchos::RCP<const Epetra_Vector> Dispnp() const
    { return DispNp(); }
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispNp() = 0;
    Teuchos::RCP<Epetra_Vector> WriteAccessDispnp()
    { return WriteAccessDispNp(); }

    /// known displacements at \f$t_{n}\f$
    /// ToDo Replace the deprecated version with the new version
    virtual Teuchos::RCP<const Epetra_Vector> DispN() const = 0;
    Teuchos::RCP<const Epetra_Vector> Dispn() const
    { return DispN(); }
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispN() = 0;
    Teuchos::RCP<Epetra_Vector> WriteAccessDispn()
    { return WriteAccessDispN(); }

    /// unknown velocity at \f$t_{n+1}\f$
    /// ToDo Replace the deprecated version with the new version
    virtual Teuchos::RCP<const Epetra_Vector> VelNp() const = 0;
    Teuchos::RCP<const Epetra_Vector> Velnp() const
    { return VelNp(); }
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessVelNp() = 0;
    Teuchos::RCP<Epetra_Vector> WriteAccessVelnp()
    { return WriteAccessVelNp(); }

    /// known velocity at \f$t_{n}\f$
    /// ToDo Replace the deprecated version with the new version
    virtual Teuchos::RCP<const Epetra_Vector> VelN() const = 0;
    Teuchos::RCP<const Epetra_Vector> Veln() const
    { return VelN(); }
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessVelN() = 0;
    Teuchos::RCP<Epetra_Vector> WriteAccessVeln()
    { return WriteAccessVelN(); }

    /// known velocity at \f$t_{n-1}\f$
    /// ToDo Replace the deprecated version with the new version
    virtual Teuchos::RCP<const Epetra_Vector> VelNm() const = 0;
    Teuchos::RCP<const Epetra_Vector> Velnm() const
    { return VelNm(); }

    /// unknown acceleration at \f$t_{n+1}\f$
    /// ToDo Replace the deprecated version with the new version
    virtual Teuchos::RCP<const Epetra_Vector> AccNp() const = 0;
    Teuchos::RCP<const Epetra_Vector> Accnp() const
    { return AccNp(); }

    /// known acceleration at \f$t_{n}\f$
    /// ToDo Replace the deprecated version with the new version
    virtual Teuchos::RCP<const Epetra_Vector> AccN() const = 0;
    Teuchos::RCP<const Epetra_Vector> Accn() const
    { return AccN(); }

    /// resize the multi step class vector
    virtual void ResizeMStepTimAda() = 0;
    ///@}

    /// @name Time step helpers
    ///@{

    /// return time integration factor
    virtual double TimIntParam() const = 0;

    /// Return current time \f$t_{n}\f$
    /// ToDo Replace the deprecated version with the new version
    virtual double GetTimeN() const = 0;
    double TimeOld() const
    { return GetTimeN(); }

    /// Sets the current time \f$t_{n}\f$
    /// ToDo Replace the deprecated version with the new version
    virtual void SetTimeN(const double time_n) = 0;
    void SetTime(const double time_n)
    { SetTimeN(time_n); }

    /// Return target time \f$t_{n+1}\f$
    /// ToDo Replace the deprecated version with the new version
    virtual double GetTimeNp() const = 0;
    double Time() const
    { return GetTimeNp(); }

    /// Get upper limit of time range of interest
    virtual double GetTimeEnd() const = 0;

    //! Set upper limit of time range of interest
    virtual void SetTimeEnd(double timemax) = 0;

    /// Sets the target time \f$t_{n+1}\f$ of this time step
    /// ToDo Replace the deprecated version with the new version
    virtual void SetTimeNp(const double time_np) = 0;
    void SetTimen(const double time_np)
    { SetTimeNp(time_np); }

    /// Get time step size \f$\Delta t_n\f$
    /// ToDo Replace the deprecated version with the new version
    virtual double GetDeltaTime() const = 0;
    double Dt() const
    { return GetDeltaTime(); }

    /// set time step size
    /// ToDo Replace the deprecated version with the new version
    virtual void SetDeltaTime(const double dt) = 0;
    void SetDt(const double dt)
    { SetDeltaTime(dt); }

    /// Return current step number $n$
    /// ToDo Replace the deprecated version with the new version
    virtual int GetStepN() const = 0;
    int StepOld() const
    { return GetStepN(); }

    /// Sets the current step \f$n\f$
    /// ToDo Replace the deprecated version with the new version
    virtual void SetStepN(int step_n) = 0;
    void SetStep(int step_n)
    { SetStepN(step_n); }

    /// Return current step number $n+1$
    /// ToDo Replace the deprecated version with the new version
    virtual int GetStepNp() const = 0;
    int Step() const
    { return GetStepNp(); }

    /// Sets the current step \f$n+1\f$
    /// ToDo Replace the deprecated version with the new version
    virtual void SetStepNp(int step_np) = 0;
    void SetStepn(int step_np)
    { SetStepNp(step_np); }

    /// Get number of time steps
    /// ToDo Replace the deprecated version with the new version
    virtual int GetStepEnd() const = 0;
    int NumStep() const
    { return GetStepEnd(); }

    /// Take the time and integrate (time loop)
    /// \date 11/08
    virtual int Integrate() = 0;

    /// do something in case nonlinear solution does not converge for some reason
    /// fixme: this can go when the old structure time integration is gone and PerformErrorAction is only called in STR::TIMINT::Implicit::Solve()
    /// and not on the structure in the adapter time loop
    INPAR::STR::ConvergenceStatus PerformErrorAction(INPAR::STR::ConvergenceStatus nonlinsoldiv) {dserror("You should not be here"); return nonlinsoldiv;};

    /// tests if there are more time steps to do
    virtual bool NotFinished() const = 0;

    /// start new time step
    virtual void PrepareTimeStep() = 0;

    /*!
    \brief update displacement and evaluate elements

    There are two displacement increments possible

    \f$x^n+1_i+1 = x^n+1_i + disiterinc\f$  (sometimes referred to as residual increment), and

    \f$x^n+1_i+1 = x^n     + disstepinc\f$

    with \f$n\f$ and \f$i\f$ being time and Newton iteration step

    Note: The structure expects an iteration increment.
    In case the StructureNOXCorrectionWrapper is applied, the step increment is expected
    which is then transformed into an iteration increment
    */
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> disiterinc ///< displacement increment between Newton iteration i and i+1
        ) = 0;

    /// update at time step end
    virtual void Update() = 0;

    /// update at time step end in case of FSI time adaptivity
    virtual void Update(double endtime) = 0;

    /// Update iteration
    /// Add residual increment to Lagrange multipliers stored in Constraint manager
    virtual void UpdateIterIncrConstr(Teuchos::RCP<Epetra_Vector> lagrincr) = 0;

    /// Update iteration
    /// Add residual increment to pressures stored in Cardiovascular0D manager
    virtual void UpdateIterIncrCardiovascular0D(Teuchos::RCP<Epetra_Vector> presincr) = 0;

    /// Access to output object
    virtual Teuchos::RCP<IO::DiscretizationWriter> DiscWriter() = 0;

    /// prepare output (i.e. calculate stresses, strains, energies)
    virtual void PrepareOutput() = 0;

    // Get restart data
    virtual void GetRestartData
    (
      Teuchos::RCP<int> step,
      Teuchos::RCP<double> time,
      Teuchos::RCP<Epetra_Vector> disn,
      Teuchos::RCP<Epetra_Vector> veln,
      Teuchos::RCP<Epetra_Vector> accn,
      Teuchos::RCP<std::vector<char> > elementdata,
      Teuchos::RCP<std::vector<char> > nodedata
    ) = 0;

    /// output results
    virtual void Output(bool forced_writerestart = false) = 0;

    /// output results to screen
    virtual void PrintStep() = 0;

    /// read restart information for given time step
    virtual void ReadRestart(const int step) = 0;

    /*!
    \brief Reset time step

    In case of time step size adaptivity, time steps might have to be repeated.
    Therefore, we need to reset the solution back to the initial solution of the
    time step.

    \author mayr.mt
    \date 08/2013
    */
    virtual void ResetStep () = 0;

    /// set restart information for parameter continuation
    virtual void SetRestart
    (
      int step,
      double time,
      Teuchos::RCP<Epetra_Vector> disn,
      Teuchos::RCP<Epetra_Vector> veln,
      Teuchos::RCP<Epetra_Vector> accn,
      Teuchos::RCP<std::vector<char> > elementdata,
      Teuchos::RCP<std::vector<char> > nodedata
    ) = 0;

    /// wrapper for things that should be done before PrepareTimeStep is called
    virtual void PrePredict() = 0;

    /// wrapper for things that should be done after PrepareTimeStep is called
    virtual void PostPredict() = 0;

    /// wrapper for things that should be done before solving the nonlinear iterations
    virtual void PreSolve() = 0;

    /// wrapper for things that should be done after solving the nonlinear iterations
    virtual void PostSolve() = 0;

    /// wrapper for things that should be done before updating
    virtual void PreUpdate() = 0;

    /// wrapper for things that should be done after solving the update
    virtual void PostUpdate() = 0;

    /// wrapper for things that should be done before the output
    virtual void PreOutput() = 0;

    /// wrapper for things that should be done after the output
    virtual void PostOutput() = 0;

    /// wrapper for things that should be done after the actual time loop is finished
    virtual void PostTimeLoop() = 0;

    ///@}

    //! @name Solver calls

    /*!
    \brief nonlinear solve

    Do the nonlinear solve, i.e. (multiple) corrector,
    for the time step. All boundary conditions have
    been set.
    */
    virtual INPAR::STR::ConvergenceStatus Solve() = 0;

    /*!
    \brief linear structure solve with just a interface load

    The very special solve done in steepest descent relaxation
    calculation (and matrix free Newton Krylov).

    \note Can only be called after a valid structural solve.
    */
    virtual Teuchos::RCP<Epetra_Vector> SolveRelaxationLinear() = 0;

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() = 0;

    //@}

    /// extract rhs (used to calculate reaction force for post-processing)
    virtual Teuchos::RCP<Epetra_Vector> Freact() = 0;


    //! @name volume coupled specific methods
    //@{

    /// Set forces due to interface with fluid, the force is expected external-force-like
    virtual void SetForceInterface(Teuchos::RCP<Epetra_MultiVector> iforce) = 0;

    //! specific method for iterative staggered partitioned TSI

    /// Identify residual
    /// This method does not predict the target solution but
    /// evaluates the residual and the stiffness matrix.
    /// In partitioned solution schemes, it is better to keep the current
    /// solution instead of evaluating the initial guess (as the predictor)
    /// does.
    virtual void PreparePartitionStep() = 0;

    //@}

    /// @name Structure with ale specific methods
    ///@{
    /// unknown material displacements at \f$t_{n+1}\f$
    /// ToDo Replace the deprecated version with the new version
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispMatNp() = 0;
    Teuchos::RCP<Epetra_Vector> DispMat()
    { return WriteAccessDispMatNp(); }

    /// set/apply material displacements to structure field (structure with ale)
    virtual void SetDispMatNp(Teuchos::RCP<Epetra_Vector> dispmatnp) = 0;
    void ApplyDisMat(Teuchos::RCP<Epetra_Vector> dismat)
    { SetDispMatNp(dismat); };
    ///@}

    /// @name Misc
    ///@{
    /// dof map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() = 0;

    /// DOF map of vector of unknowns for multiple dofsets
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds) = 0;

    /// DOF map view of vector of unknowns
    virtual const Epetra_Map* DofRowMapView() = 0;

    /// domain map of system matrix (do we really need this?)
    /// ToDo Replace the deprecated version with the new version
    virtual const Epetra_Map& GetMassDomainMap() const = 0;
    virtual const Epetra_Map& DomainMap() const
    { return GetMassDomainMap(); }

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() = 0;

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() = 0;

    /// switch structure field to block matrix
    virtual void UseBlockMatrix(Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps,
        Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps) = 0;

    /// return contact/meshtying bridge
    virtual Teuchos::RCP<CONTACT::MeshtyingContactBridge> MeshtyingContactBridge() = 0;

    // access to locsys manager
    virtual Teuchos::RCP<DRT::UTILS::LocsysManager> LocsysManager() = 0;

    /// direct access to discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization() = 0;

    /// are there any algebraic constraints?
    virtual bool HaveConstraint() = 0;

    /// do we have semi-smooth Newton-type plasticity algorithm
    virtual bool HaveSemiSmoothPlasticity() = 0;

    /// get constraint manager defined in the structure
    virtual const Teuchos::RCP<UTILS::ConstrManager> GetConstraintManager() = 0;

    /// get plasticity manager defined in the structure
    virtual const Teuchos::RCP<DRT::UTILS::PlastSsnManager> GetPlasticityManager() = 0;

    /// Get type of thickness scaling for thin shell structures
    virtual INPAR::STR::STC_Scale GetSTCAlgo() = 0;

    /// Access to scaling matrix for STC
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetSTCMat() = 0;

    /// Return MapExtractor for Dirichlet boundary conditions
    virtual const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() = 0;

    /// create result test for encapsulated structure algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest() = 0;

    /// reset time and state vectors (needed for biofilm growth simulations)
    virtual void Reset() = 0;

    /// set structure displacement vector due to biofilm growth
    virtual void SetStrGrDisp(Teuchos::RCP<Epetra_Vector> struct_growth_disp) = 0;

    /// bool indicating if micro material is used
    virtual bool HaveMicroMat() = 0;

    /// return statmech manager
    virtual Teuchos::RCP<STATMECH::StatMechManager> StatMechManager() = 0;

    //! return statmech managerbilayer
    virtual Teuchos::RCP<STATMECH::StatMechManagerBilayer> StatMechManagerBilayer() = 0;

    ///@}

    /// @name Crack propagation specific methods
    /// @{
    /// Propagate crack within the structure and accordingly modify the discretization
    virtual bool UpdateCrackInformation(Teuchos::RCP<const Epetra_Vector> displace) = 0;

    /// Write Gmsh output for structural field
    virtual void writeGmshStrucOutputStep() = 0;

    /// Get map of newly introduced node ids with their corresponding old node ids
    virtual std::map<int,int> getOldNewCrackNodes() = 0;

    /// Get new crack tip nodes after propagation
    virtual std::vector<int> GetCrackTipNodes() = 0;
    /// @}



    /// @name Currently unused functions, which will be deleted in the near future,
    /// if they stay unnecessary.
    ///@{
    /// are there any spring dashpot bcs?
    virtual bool HaveSpringDashpot()
    { dserror("This function seems to be unused!"); return false; }

    /// get SpringDashpot manager defined in the structure
    virtual const Teuchos::RCP<UTILS::SpringDashpotManager> GetSpringDashpotManager()
    { dserror("This function seems to be unused!"); return Teuchos::null; }

    ///@}


  };  // class StructureNew

  /// structure field solver
  class StructureBaseAlgorithmNew
  {
  public:

    /// constructor
    StructureBaseAlgorithmNew();

    /// virtual destructor to support polymorph destruction
    virtual ~StructureBaseAlgorithmNew() {};

    /// initialize all class internal variables
    virtual void Init(const Teuchos::ParameterList& prbdyn, Teuchos::ParameterList& sdyn, Teuchos::RCP<DRT::Discretization> actdis);

    /// setup
    virtual void Setup();

    /// structural field solver
    Teuchos::RCP<Structure> StructureField() { return str_wrapper_; }

  protected:

    inline const bool& IsInit() const { return isinit_; };

    inline const bool& IsSetup() const { return issetup_; };

  protected:
    /// setup structure algorithm of STR::TimInt::Implicit or STR::TimInt::Explicit type
    void SetupTimInt();

    /// Set all found model types. This is e.g. necessary for the model evaluator.
    void SetModelTypes(
        std::set<enum INPAR::STR::ModelType>& modeltypes) const;

    /// Set all found model types.
    void DetectElementTechnologies(
        std::set<enum INPAR::STR::EleTech>& eletechs) const;

    /// Set different time integrator specific parameters in the different parameter lists
    virtual void SetParams(
        Teuchos::ParameterList& ioflags,
        Teuchos::ParameterList& xparams,
        Teuchos::ParameterList& timeadaptivity
        );

    /// set the final structure time integrator object
    virtual void SetStructureWrapper(
        const Teuchos::ParameterList& ioflags,
        const Teuchos::ParameterList& sdyn,
        const Teuchos::ParameterList& xparams,
        const Teuchos::ParameterList& timeadaptivity,
        Teuchos::RCP< ::STR::TIMINT::Base> ti_strategy);

    /// create the adaptive time integrator wrapper
    void CreateAdaptiveWrapper(
        const Teuchos::ParameterList& ioflags,
        const Teuchos::ParameterList& sdyn,
        const Teuchos::ParameterList& xparams,
        const Teuchos::ParameterList& taflags,
        Teuchos::RCP< ::STR::TIMINT::Base> ti_strategy);

    /// create the time integrator wrapper
    void CreateWrapper(Teuchos::RCP< ::STR::TIMINT::Base> ti_strategy);

  protected:
    /// structural field solver
    Teuchos::RCP<Structure> str_wrapper_;

    /// parameter list of the problem dynamics (read only)
    Teuchos::RCP<const Teuchos::ParameterList> prbdyn_;

    /// parameter list of the structural dynamics (mutable)
    Teuchos::RCP<Teuchos::ParameterList> sdyn_;

    /// current discretization
    Teuchos::RCP<DRT::Discretization> actdis_;

    /// init flag
    bool isinit_;

    /// setup flag
    bool issetup_;
  };  // class StructureBaseAlgorithmNew
} // namespace ADAPTER

#endif /* AD_STR_STRUCTURE_NEW_H_ */

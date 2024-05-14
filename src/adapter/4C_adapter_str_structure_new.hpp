/*-----------------------------------------------------------*/
/*! \file

\brief Adapter for the new structural time integration framework.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_STR_STRUCTURE_NEW_HPP
#define FOUR_C_ADAPTER_STR_STRUCTURE_NEW_HPP

#include "4C_config.hpp"

#include "4C_adapter_str_structure.hpp"
#include "4C_inpar_structure.hpp"

#include <set>

// forward declarations
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

FOUR_C_NAMESPACE_OPEN

namespace STR
{
  namespace TIMINT
  {
    class Base;
    class BaseDataGlobalState;
    class BaseDataSDyn;
    class BaseDataIO;
  }  // namespace TIMINT
  namespace MODELEVALUATOR
  {
    class Generic;
  }  // namespace MODELEVALUATOR
}  // namespace STR

namespace CORE::LINALG
{
  class Solver;
}  // namespace CORE::LINALG

namespace DRT
{
  class Discretization;
}  // namespace DRT

namespace ADAPTER
{
  class StructureNew : public Structure
  {
   public:
    /// @name General methods
    ///@{
    /// Setup the structure integrator
    void Setup() override = 0;
    ///@}

    /// @name Vector access
    ///@{
    /// initial guess of Newton's method
    Teuchos::RCP<const Epetra_Vector> InitialGuess() override = 0;

    /// rhs of Newton's method
    Teuchos::RCP<const Epetra_Vector> RHS() override = 0;

    /// unknown displacements at \f$t_{n+1}\f$
    /// ToDo Replace the deprecated version with the new version
    [[nodiscard]] virtual Teuchos::RCP<const Epetra_Vector> DispNp() const = 0;
    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> Dispnp() const override { return DispNp(); }
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispNp() = 0;
    Teuchos::RCP<Epetra_Vector> WriteAccessDispnp() override { return WriteAccessDispNp(); }

    /// known displacements at \f$t_{n}\f$
    /// ToDo Replace the deprecated version with the new version
    [[nodiscard]] virtual Teuchos::RCP<const Epetra_Vector> DispN() const = 0;
    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> Dispn() const override { return DispN(); }
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispN() = 0;
    Teuchos::RCP<Epetra_Vector> WriteAccessDispn() override { return WriteAccessDispN(); }

    /// unknown velocity at \f$t_{n+1}\f$
    /// ToDo Replace the deprecated version with the new version
    [[nodiscard]] virtual Teuchos::RCP<const Epetra_Vector> VelNp() const = 0;
    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> Velnp() const override { return VelNp(); }
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessVelNp() = 0;
    Teuchos::RCP<Epetra_Vector> WriteAccessVelnp() override { return WriteAccessVelNp(); }

    /// known velocity at \f$t_{n}\f$
    /// ToDo Replace the deprecated version with the new version
    [[nodiscard]] virtual Teuchos::RCP<const Epetra_Vector> VelN() const = 0;
    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> Veln() const override { return VelN(); }
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessVelN() = 0;
    Teuchos::RCP<Epetra_Vector> WriteAccessVeln() override { return WriteAccessVelN(); }

    /// known velocity at \f$t_{n-1}\f$
    /// ToDo Replace the deprecated version with the new version
    [[nodiscard]] virtual Teuchos::RCP<const Epetra_Vector> VelNm() const = 0;
    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> Velnm() const override { return VelNm(); }

    /// unknown acceleration at \f$t_{n+1}\f$
    /// ToDo Replace the deprecated version with the new version
    [[nodiscard]] virtual Teuchos::RCP<const Epetra_Vector> AccNp() const = 0;
    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> Accnp() const override { return AccNp(); }

    /// known acceleration at \f$t_{n}\f$
    /// ToDo Replace the deprecated version with the new version
    [[nodiscard]] virtual Teuchos::RCP<const Epetra_Vector> AccN() const = 0;
    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> Accn() const override { return AccN(); }

    /// resize the multi step class vector
    void ResizeMStepTimAda() override = 0;
    ///@}

    /// @name Time step helpers
    ///@{

    /// return time integration factor
    [[nodiscard]] double TimIntParam() const override = 0;

    /// Return current time \f$t_{n}\f$
    /// ToDo Replace the deprecated version with the new version
    [[nodiscard]] virtual double GetTimeN() const = 0;
    [[nodiscard]] double TimeOld() const override { return GetTimeN(); }

    /// Sets the current time \f$t_{n}\f$
    /// ToDo Replace the deprecated version with the new version
    virtual void SetTimeN(const double time_n) = 0;
    void SetTime(const double time_n) override { SetTimeN(time_n); }

    /// Return target time \f$t_{n+1}\f$
    /// ToDo Replace the deprecated version with the new version
    [[nodiscard]] virtual double GetTimeNp() const = 0;
    [[nodiscard]] double Time() const override { return GetTimeNp(); }

    /// Get upper limit of time range of interest
    [[nodiscard]] double GetTimeEnd() const override = 0;

    //! Set upper limit of time range of interest
    void SetTimeEnd(double timemax) override = 0;

    /// Sets the target time \f$t_{n+1}\f$ of this time step
    /// ToDo Replace the deprecated version with the new version
    virtual void SetTimeNp(const double time_np) = 0;
    void SetTimen(const double time_np) override { SetTimeNp(time_np); }

    /// Get time step size \f$\Delta t_n\f$
    /// ToDo Replace the deprecated version with the new version
    [[nodiscard]] virtual double GetDeltaTime() const = 0;
    [[nodiscard]] double Dt() const override { return GetDeltaTime(); }

    /// set time step size
    /// ToDo Replace the deprecated version with the new version
    virtual void SetDeltaTime(const double dt) = 0;
    void SetDt(const double dt) override { SetDeltaTime(dt); }

    /// Return current step number $n$
    /// ToDo Replace the deprecated version with the new version
    [[nodiscard]] virtual int GetStepN() const = 0;
    [[nodiscard]] int StepOld() const override { return GetStepN(); }

    /// Sets the current step \f$n\f$
    /// ToDo Replace the deprecated version with the new version
    virtual void SetStepN(int step_n) = 0;
    void SetStep(int step_n) override { SetStepN(step_n); }

    /// Return current step number $n+1$
    /// ToDo Replace the deprecated version with the new version
    [[nodiscard]] virtual int GetStepNp() const = 0;
    [[nodiscard]] int Step() const override { return GetStepNp(); }

    /// Sets the current step \f$n+1\f$
    /// ToDo Replace the deprecated version with the new version
    virtual void SetStepNp(int step_np) = 0;
    void SetStepn(int step_np) override { SetStepNp(step_np); }

    /// Get number of time steps
    /// ToDo Replace the deprecated version with the new version
    [[nodiscard]] virtual int GetStepEnd() const = 0;
    [[nodiscard]] int NumStep() const override { return GetStepEnd(); }

    /// Sets number of time steps (in case of time adaptivity)
    virtual void SetStepEnd(int step_end) = 0;

    /// Take the time and integrate (time loop)
    /// \date 11/08
    int Integrate() override = 0;

    /// fixme: this can go when the old structure time integration is gone and PerformErrorAction is
    /// only called in STR::TIMINT::Implicit::Solve() and not on the structure in the adapter time
    /// loop
    INPAR::STR::ConvergenceStatus PerformErrorAction(
        INPAR::STR::ConvergenceStatus nonlinsoldiv) override
    {
      FOUR_C_THROW("You should not be here");
      return nonlinsoldiv;
    };

    /// tests if there are more time steps to do
    [[nodiscard]] bool NotFinished() const override = 0;

    /// start new time step
    void PrepareTimeStep() override = 0;

    /*!
     \brief update displacement
     There are two displacement increments possible

     \f$x^n+1_i+1 = x^n+1_i + disiterinc\f$  (sometimes referred to as residual increment), and

     \f$x^n+1_i+1 = x^n     + disstepinc\f$

     with \f$n\f$ and \f$i\f$ being time and Newton iteration step

     Note: The structure expects an iteration increment.
     In case the StructureNOXCorrectionWrapper is applied, the step increment is expected
     which is then transformed into an iteration increment
     */
    void UpdateStateIncrementally(Teuchos::RCP<const Epetra_Vector>
            disiterinc  ///< displacement increment between Newton iteration i and i+1
        ) override = 0;

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
    void Evaluate(Teuchos::RCP<const Epetra_Vector>
            disiterinc  ///< displacement increment between Newton iteration i and i+1
        ) override = 0;

    /// don't update displacement but evaluate elements (implicit only)
    void Evaluate() override = 0;

    /// update at time step end
    void Update() override = 0;

    /// update at time step end in case of FSI time adaptivity
    void Update(double endtime) override = 0;

    /// Update iteration
    /// Add residual increment to Lagrange multipliers stored in Constraint manager
    void UpdateIterIncrConstr(Teuchos::RCP<Epetra_Vector> lagrincr) override = 0;

    /// Update iteration
    /// Add residual increment to pressures stored in Cardiovascular0D manager
    void UpdateIterIncrCardiovascular0D(Teuchos::RCP<Epetra_Vector> presincr) override = 0;

    /// Access to output object
    Teuchos::RCP<IO::DiscretizationWriter> DiscWriter() override = 0;

    /// prepare output (i.e. calculate stresses, strains, energies)
    void PrepareOutput(bool force_prepare_timestep) override = 0;

    // Get restart data
    void GetRestartData(Teuchos::RCP<int> step, Teuchos::RCP<double> time,
        Teuchos::RCP<Epetra_Vector> disn, Teuchos::RCP<Epetra_Vector> veln,
        Teuchos::RCP<Epetra_Vector> accn, Teuchos::RCP<std::vector<char>> elementdata,
        Teuchos::RCP<std::vector<char>> nodedata) override = 0;

    /// output results
    void Output(bool forced_writerestart = false) override = 0;

    /// output results to screen
    void PrintStep() override = 0;

    /// read restart information for given time step
    void ReadRestart(const int step) override = 0;

    /*!
    \brief Reset time step

    In case of time step size adaptivity, time steps might have to be repeated.
    Therefore, we need to reset the solution back to the initial solution of the
    time step.

    \author mayr.mt
    \date 08/2013
    */
    void ResetStep() override = 0;

    /// set restart information for parameter continuation
    void SetRestart(int step, double time, Teuchos::RCP<Epetra_Vector> disn,
        Teuchos::RCP<Epetra_Vector> veln, Teuchos::RCP<Epetra_Vector> accn,
        Teuchos::RCP<std::vector<char>> elementdata,
        Teuchos::RCP<std::vector<char>> nodedata) override = 0;

    /// wrapper for things that should be done before PrepareTimeStep is called
    void PrePredict() override = 0;

    /// wrapper for things that should be done before solving the nonlinear iterations
    void PreSolve() override = 0;

    /// wrapper for things that should be done before updating
    void PreUpdate() override = 0;

    /// wrapper for things that should be done after solving the update
    void PostUpdate() override = 0;

    /// wrapper for things that should be done after the output
    void PostOutput() override = 0;

    /// wrapper for things that should be done after the actual time loop is finished
    void PostTimeLoop() override = 0;

    ///@}

    //! @name Solver calls

    /*!
    \brief nonlinear solve

    Do the nonlinear solve, i.e. (multiple) corrector,
    for the time step. All boundary conditions have
    been set.
    */
    INPAR::STR::ConvergenceStatus Solve() override = 0;

    /*!
    \brief linear structure solve with just a interface load

    The very special solve done in steepest descent relaxation
    calculation (and matrix free Newton Krylov).

    \note Can only be called after a valid structural solve.
    */
    Teuchos::RCP<Epetra_Vector> SolveRelaxationLinear() override
    {
      FOUR_C_THROW(
          "In the new structural timeintegration this method is"
          "no longer needed inside the structure. Since this is"
          "FSI specific, the functionality is shifted to the"
          "STR::MODELEVALUATOR::PartitionedFSI.");
      return Teuchos::null;
    };

    /// get the linear solver object used for this field
    Teuchos::RCP<CORE::LINALG::Solver> LinearSolver() override = 0;

    //@}

    /// extract rhs (used to calculate reaction force for post-processing)
    Teuchos::RCP<Epetra_Vector> Freact() override = 0;


    //! @name volume coupled specific methods
    //@{

    /// Set forces due to interface with fluid, the force is expected external-force-like
    void SetForceInterface(Teuchos::RCP<Epetra_MultiVector> iforce) override
    {
      FOUR_C_THROW(
          "This method is deprecated. In the new structural time integration"
          "this functionality is taken over by the problem specific model "
          "evaluators. Remove this method as soon as possible.");
    };

    //! specific method for iterative staggered partitioned TSI

    /// Identify residual
    /// This method does not predict the target solution but
    /// evaluates the residual and the stiffness matrix.
    /// In partitioned solution schemes, it is better to keep the current
    /// solution instead of evaluating the initial guess (as the predictor)
    /// does.
    void PreparePartitionStep() override = 0;

    //@}

    /// @name Structure with ale specific methods
    ///@{
    /// unknown material displacements at \f$t_{n+1}\f$
    /// ToDo Replace the deprecated version with the new version
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispMatNp() = 0;
    Teuchos::RCP<Epetra_Vector> DispMat() override { return WriteAccessDispMatNp(); }

    /// set/apply material displacements to structure field (structure with ale)
    virtual void SetDispMatNp(Teuchos::RCP<Epetra_Vector> dispmatnp) = 0;
    void ApplyDisMat(Teuchos::RCP<Epetra_Vector> dismat) override { SetDispMatNp(dismat); };
    ///@}

    /// @name Misc
    ///@{
    /// dof map of vector of unknowns
    Teuchos::RCP<const Epetra_Map> DofRowMap() override = 0;

    /// DOF map of vector of unknowns for multiple dofsets
    Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds) override = 0;

    /// DOF map view of vector of unknowns
    const Epetra_Map* DofRowMapView() override = 0;

    /// domain map of system matrix (do we really need this?)
    /// ToDo Replace the deprecated version with the new version
    [[nodiscard]] virtual const Epetra_Map& GetMassDomainMap() const = 0;
    [[nodiscard]] const Epetra_Map& DomainMap() const override { return GetMassDomainMap(); }

    /// direct access to system matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> SystemMatrix() override = 0;

    /// direct access to system matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> BlockSystemMatrix() override = 0;

    /// switch structure field to block matrix
    void UseBlockMatrix(Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> domainmaps,
        Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> rangemaps) override = 0;

    /// return contact/meshtying bridge
    Teuchos::RCP<CONTACT::MeshtyingContactBridge> MeshtyingContactBridge() override = 0;

    /// access to locsys manager
    Teuchos::RCP<DRT::UTILS::LocsysManager> LocsysManager() override = 0;

    /// access the desired model evaluator (read-only)
    [[nodiscard]] virtual const STR::MODELEVALUATOR::Generic& ModelEvaluator(
        INPAR::STR::ModelType mtype) const = 0;

    /// access the desired model evaluator (read and write)
    STR::MODELEVALUATOR::Generic& ModelEvaluator(INPAR::STR::ModelType mtype) override = 0;

    /// direct access to discretization
    Teuchos::RCP<DRT::Discretization> Discretization() override = 0;

    /// are there any algebraic constraints?
    bool HaveConstraint() override = 0;

    /// get constraint manager defined in the structure
    Teuchos::RCP<CONSTRAINTS::ConstrManager> GetConstraintManager() override = 0;

    /// Get type of thickness scaling for thin shell structures
    INPAR::STR::StcScale GetSTCAlgo() override = 0;

    /// Access to scaling matrix for STC
    Teuchos::RCP<CORE::LINALG::SparseMatrix> GetSTCMat() override = 0;

    /// Return MapExtractor for Dirichlet boundary conditions
    Teuchos::RCP<const CORE::LINALG::MapExtractor> GetDBCMapExtractor() override = 0;

    /// create result test for encapsulated structure algorithm
    Teuchos::RCP<CORE::UTILS::ResultTest> CreateFieldTest() override = 0;

    /// reset time and state vectors (needed for biofilm growth simulations)
    void Reset() override = 0;

    /// set structure displacement vector due to biofilm growth
    void SetStrGrDisp(Teuchos::RCP<Epetra_Vector> struct_growth_disp) override = 0;

    /// bool indicating if micro material is used
    bool HaveMicroMat() override = 0;

    ///@}

    /// @name Currently unused functions, which will be deleted in the near future,
    /// if they stay unnecessary.
    ///@{
    /// are there any spring dashpot bcs?
    bool HaveSpringDashpot() override
    {
      FOUR_C_THROW("This function seems to be unused!");
      return false;
    }

    /// get SpringDashpot manager defined in the structure
    Teuchos::RCP<CONSTRAINTS::SpringDashpotManager> GetSpringDashpotManager() override
    {
      FOUR_C_THROW("This function seems to be unused!");
      return Teuchos::null;
    }

    ///@}

    /// @name Multiphysics related stuff
    ///@{

    /** \brief Set the state of the NOX group and the global state data container.
     *
     * This method is needed because there are two parallel ways to handle the
     * global state in the 'new' structural time integration.
     *
     * 1) The current state is held in the global state data container:
     *    \ref STR::TIMINT::BaseDataGlobalState
     *
     * 2) Also the NOX group (that means the nonlinear solver) has its
     *    own state vector (called 'X').
     *
     * This method sets the provided state consistently in both objects.
     *
     * This is useful for multiphysics in case a manipulated state needs
     * to be set from outside:
     *
     *   For example, see the class:
     *   \ref IMMERSED::ImmersedPartitionedAdhesionTraction
     *
     *   There, see the method
     *   \ref IMMERSED::ImmersedPartitionedAdhesionTraction::ImmersedOp
     *
     *   First, the displacement state with write access is request from
     *   the time integrator via the corresponding adapter.
     *
     *   Then, the displacement state is manipulated.
     *
     *   Finally, the displacement state is handed over to this method in
     *   order to apply the new state to both (i) the global state object;
     *   and (ii) the NOX group.
     *
     *   \note velocities and accelerations are recalculated inside by invoking
     *   SetState(x) on the concrete time integrator (e.g. OST, GenAlpha, etc.)
     *   It never makes any sense to call velocities or displacements as WriteAccess
     *   variant from outside, because these vectors should always be consistent with
     *   our primary variable (i.e. the displacements).
     *
     *  \author rauch
     *  \date 10/17
     */
    void SetState(const Teuchos::RCP<Epetra_Vector>& x) override = 0;

    ///@}

  };  // class StructureNew

  /// structure field solver
  class StructureBaseAlgorithmNew
  {
   public:
    /// constructor
    StructureBaseAlgorithmNew();

    /// virtual destructor to support polymorph destruction
    virtual ~StructureBaseAlgorithmNew() = default;

    /// initialize all class internal variables
    virtual void Init(const Teuchos::ParameterList& prbdyn, Teuchos::ParameterList& sdyn,
        Teuchos::RCP<DRT::Discretization> actdis);

    /// setup
    virtual void Setup();

    /** \brief Register an externally created model evaluator.
     *
     *  This can be used e.g. by coupled problems.
     *
     *  \date 11/16 */
    void RegisterModelEvaluator(
        const std::string name, Teuchos::RCP<STR::MODELEVALUATOR::Generic> me);

    /// structural field solver
    Teuchos::RCP<Structure> StructureField() { return str_wrapper_; }

   public:
    [[nodiscard]] inline const bool& IsInit() const { return isinit_; };

    [[nodiscard]] inline const bool& IsSetup() const { return issetup_; };

   protected:
    /// setup structure algorithm of STR::TimInt::Implicit or STR::TimInt::Explicit type
    void SetupTimInt();

    /** \brief Set all model types. This is necessary for the model evaluation.
     *
     *  The inherent structural models are identified by the corresponding conditions and/or
     *  other unique criteria. If your intention is to solve a partitioned coupled problem and
     *  you need to modify the structural right-hand-side in any way, then you have to implement
     *  your own concrete implementation of a STR::MODELEVALUATOR::Generic class and register it
     *  as a Teuchos::RCP<STR::MODELEVALUATOR::Generic> pointer in your problem dynamic parameter-
     *  list. For partitioned problems you have to use the parameter-name
     *  \"Partitioned Coupling Model\".
     *
     *  For example: To create and use a coupling model evaluator for the partitioned FSI you
     *  have to insert the object as follows:
     *
     *  <ol>
     *
     *  <li> Create a model evaluator that derives from
     *  STR::MODELEVALUATOR::Generic. For example, the model evaluator
     *  \c FSI_Partitioned might be defined as shown below.
     *
     *  \code
     *  class FSI_Partitioned : public STR::MODELEVALUATOR::Generic
     *  {
     *  // Insert class definition here
     *  }
     *  \endcode
     *
     *  <li> Create the appropriate entries in your problem dynamic parameter-list \c prbdyn
     *  and initialize member variables, if desired (optional):
     *
     *  \code
     *  Teuchos::RCP<FSI_Partitioned> fsi_model_ptr = Teuchos::rcp(new FSI_Partitioned());
     *  // optional: call of your own 2-nd Init() method
     *  fsi_model_ptr->Init(stuff_you_need_inside_the_model_evaluator);
     *  prbdyn.set<Teuchos::RCP<STR::MODELEVALUATOR::Generic> >("Partitioned Coupling Model",
     *      fsi_model_ptr);
     *  \endcode
     *
     *  </ol>
     *
     *  \remark Please keep in mind, that the prescribed Generic::Init() and Generic::Setup()
     *  methods will be called automatically in the STR::ModelEvaluator::Setup() routine. If
     *  you need a different Init() method, just define a second Init() function with different
     *  input variables in your concrete class implementation and call it somewhere in your code
     *  (see upper example code).
     *  The constructor is supposed to stay empty. If you need a safety check, you can overload
     *  the Generic::CheckInit() and Generic::CheckInitSetup() routines, instead.
     *
     *  \author hiermeier
     *  \date 09/16 */
    void SetModelTypes(std::set<enum INPAR::STR::ModelType>& modeltypes) const;

    /// Set all found model types.
    void DetectElementTechnologies(std::set<enum INPAR::STR::EleTech>& eletechs) const;

    /// Set different time integrator specific parameters in the different parameter lists
    virtual void SetParams(Teuchos::ParameterList& ioflags, Teuchos::ParameterList& xparams,
        Teuchos::ParameterList& time_adaptivity_params);

    /// Create, initialize and setup the global state data container
    virtual void SetGlobalState(Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& dataglobalstate,
        const Teuchos::RCP<const STR::TIMINT::BaseDataSDyn>& datasdyn_ptr);

    /// Create, initialize and setup the time integration strategy object
    virtual void SetTimeIntegrationStrategy(Teuchos::RCP<STR::TIMINT::Base>& ti_strategy,
        const Teuchos::RCP<STR::TIMINT::BaseDataIO>& dataio,
        const Teuchos::RCP<STR::TIMINT::BaseDataSDyn>& datasdyn,
        const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& dataglobalstate, const int& restart);

    /// set the final structure time integrator object
    virtual void SetStructureWrapper(const Teuchos::ParameterList& ioflags,
        const Teuchos::ParameterList& sdyn, const Teuchos::ParameterList& xparams,
        const Teuchos::ParameterList& time_adaptivity_params,
        Teuchos::RCP<STR::TIMINT::Base> ti_strategy);

    /// create the time integrator wrapper
    void CreateWrapper(Teuchos::RCP<STR::TIMINT::Base> ti_strategy);

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
}  // namespace ADAPTER

FOUR_C_NAMESPACE_CLOSE

#endif

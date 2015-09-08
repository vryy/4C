/*
 * str_timint_base.H
 *
 *  Created on: Aug 12, 2015
 *      Author: farah
 */

#ifndef STR_TI_TIMINT_H_
#define STR_TI_TIMINT_H_

#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_lib/drt_utils_timintmstep.H"

// forward declaration
class Epetra_Vector;
class Epetra_Time;

namespace Teuchos {
  class ParameterList;
}
namespace LINALG {
  class SparseOperator;
}
namespace IO {
  class DiscretizationWriter;
}
namespace DRT {
  class Discretization;
}

namespace STR {
namespace TIMINT {
  class BaseDataIO
  {
  public:
    BaseDataIO();

    virtual ~BaseDataIO() {};

  public:

    void Init(const Teuchos::ParameterList& IOParams,
        const Teuchos::ParameterList& sDynParams,
        const Teuchos::ParameterList& xParams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    void Setup();

  protected:
    const bool& IsInit();

    const bool& IsSetup();

  protected:
    //! @name variables for internal use only
    //@{

    bool isInit_;

    bool isSetup_;
    //@}

  private:
    //! @name Printing and output
    //@{

    //! binary output
    Teuchos::RCP<IO::DiscretizationWriter> output_;

    //! outputfile for energy
    Teuchos::RCP<std::ofstream> energyfile_;

    //! error file handle
    FILE* errfile_;

    //! Is GMSH output of displacements required?
    bool gmsh_out_;

    //! print the logo (or not)?
    bool printlogo_;

    //! print infos to error file #errfile_
    bool printerrfile_;

    //! print intermediate iterations during solution
    bool printiter_;

    //! switch
    bool outputeveryiter_;

    //! write surfactant output
    bool writesurfactant_;

    //! write state on/off
    bool writestate_;

    //! write velocity and acceleration on/off
    bool writevelacc_;

    //! print infos to standard out every n steps
    int printscreen_;

    //! filename counter
    int oei_filecounter_;

    //! output counter for OutputEveryIter
    int outputcounter_;

    //! write restart every given step. if 0, restart is not written
    int writerestartevery_;

    //! Write reduced Restart Data i.e no mesh for mlmc purposes
    int writereducedrestart_;

    //! write state/stress/strain every given step
    int writeresultsevery_;

    //! write system energy every given step
    int writeenergyevery_;

    //! kinetic energy of system
    double kinergy_;

    //! internal energy of system
    double intergy_;

    //! external energy of system
    double extergy_;

    //! stress output type
    INPAR::STR::StressType writestress_;

    //! output type of coupling stress
    INPAR::STR::StressType writecouplstress_;

    //! strain output type
    INPAR::STR::StrainType writestrain_;

    //!< plastic strain output type
    INPAR::STR::StrainType writeplstrain_;

    //@}
  };  // BaseDataIO


  class BaseDataSDyn
  {
  public:
    BaseDataSDyn();

    virtual ~BaseDataSDyn() {};

    virtual void Init(
        const Teuchos::RCP<DRT::Discretization> discret,
        const Teuchos::ParameterList& sDynParams,
        const std::vector<const enum INPAR::STR::ModelType>& modeltypes,
        std::map<const enum INPAR::STR::ModelType,Teuchos::RCP<LINALG::Solver> >& linsolvers);

    virtual void Setup();

  protected:
    virtual const bool& IsInit() const { return isInit_; };

    virtual const bool& IsSetup() const { return isSetup_; };

    virtual const void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup())
        dserror("Call Init() and Setup() first!");
    }

  public:

    //! @name Get general control parameters (read only access)
    //@{

    //! Return final time \f$t_\text{fin}\f$
    const double& GetTimeMax() const
    {
      CheckInitSetup();
      return timemax_;
    };

    //! Return final time step \f$N\f$
    const int& GetStepMax() const
    {
      CheckInitSetup();
      return stepmax_;
    };

    //! Return timer for solution technique
    Teuchos::RCP<const Epetra_Time>  GetTimer() const
    {
      CheckInitSetup();
      return timer_;
    };

    //! Return dynamic type
    const enum INPAR::STR::DynamicType& GetDynamicType() const
    {
      CheckInitSetup();
      return dynType_;
    };

    //! Return prestress type
    const enum INPAR::STR::PreStress& GetPreStressType() const
    {
      CheckInitSetup();
      return preStressType_;
    };

    //! Return nonlinear solver type
    const enum INPAR::STR::NonlinSolTech& GetNlnSolverType() const
    {
      CheckInitSetup();
      return nlnSolverType_;
    };

    //@}

    //! @name Get damping control parameters (read only access)
    //@{
    //! Return damping type
    const enum INPAR::STR::DampKind& GetDampingType() const
    {
      CheckInitSetup();
      return damptype_;
    };

    //! Return damping factor for stiffness \f$c_\text{K}\f$
    const double& GetDampingStiffnessFactor() const
    {
      CheckInitSetup();
      return dampk_;
    };

    //! Return damping factor for mass \f$c_\text{M}\f$
    const double& GetDampingMassFactor() const
    {
      CheckInitSetup();
      return dampm_;
    };
    //@}

    //! @name Get mass and inertia control parameters (read only access)
    //@{
    //! Return mass linearization type
    const double& GetMassLinType() const
    {
      CheckInitSetup();
      return masslintype_;
    };
    //@}

    //! @name Get model evaluator control parameters (read only access)
    //@{
    //! Return types of the current models
    const std::vector<const enum INPAR::STR::ModelType> GetModelType() const
    {
      CheckInitSetup();
      return *modeltypes_;
    };

    //! Return map of the current models
    const std::map<const enum INPAR::STR::ModelType, Teuchos::RCP<STR::ModelEvaluator::Generic> >& GetModelEvaluators() const
    {
      CheckInitSetup();
      return modelevaluators_;
    };
    //@}

    //! @name Get linear solver parameters (read only access)
    //@{
    const std::map<const enum INPAR::STR::ModelType,Teuchos::RCP<LINALG::Solver> > GetLinSolvers() const
    {
      CheckInitSetup();
      return *linsolvers_;
    };
    //@}

    //! @name Get implicit time integrator parameters (read only access)
    //@{
    //! Return specific implicit time integrator pointer
    //! NOTE: Only filled for the implicit case.
    const STR::IMPLICIT::Generic& GetImplicitIntegrator() const
    {
      CheckInitSetup();
      if (implint_.is_null())
        dserror("Implicit Integrator is not initialized!");
      return *implint_;
    };

    //! non-linear solver pointer
    //! NOTE: Only filled for the implicit case
    const STR::NLN::SOLVER::Generic& GetNlnSolvers() const
    {
      CheckInitSetup();
      if (nlnsolvers_.is_null())
        dserror("Non-linear solvers are not initialized!");
      return *nlnsolvers_;
    };

    //@}

    //! @name Get mutable general control parameters (read and write access)
    //@{

    //! Return final time \f$t_\text{fin}\f$
    double& GetMutableTimeMax()
    {
      CheckInitSetup();
      return timemax_;
    };

    //! Return final time step \f$N\f$
    int& GetMutableStepMax()
    {
      CheckInitSetup();
      return stepmax_; };

    //! Return timer for solution technique
    Teuchos::RCP<Epetra_Time>& GetMutableTimer()
    {
      CheckInitSetup();
      return timer_;
    };

    //! Return prestress type
    enum INPAR::STR::PreStress& GetMutablePreStressType()
    {
      CheckInitSetup();
      return preStressType_;
    };

    //! Return dynamic type
    enum INPAR::STR::DynamicType& GetMutableDynamicType()
    {
      CheckInitSetup();
      return dynType_;
    };

    //! Return nonlinear solver type
    enum INPAR::STR::NonlinSolTech& GetMutableNlnSolverType()
    {
      CheckInitSetup();
      return nlnSolverType_;
    };

    //@}

    //! @name Get mutable damping control parameters (read only access)
    //@{
    //! Return damping type
    enum INPAR::STR::DampKind& GetMutableDampingType()
    {
      CheckInitSetup();
      return damptype_;
    };

    //! Return damping factor for stiffness \f$c_\text{K}\f$
    double& GetMutableDampingStiffnessFactor()
    {
      CheckInitSetup();
      return dampk_;
    };

    //! Return damping factor for mass \f$c_\text{M}\f$
    double& GetMutableDampingMassFactor()
    {
      CheckInitSetup();
      return dampm_;
    };
    //@}

    //! @name Get mutable linear solver parameters (read and write access)
    //@{
    std::map<const enum INPAR::STR::ModelType,Teuchos::RCP<LINALG::Solver> > GetMutableLinSolvers()
    {
      CheckInitSetup();
      return *linsolvers_;
    };
    //@}

    //! @name Get mutable mass and inertia control parameters (read and write access)
    //@{
    //! Return mass linearization type
    double& GetMutableMassLinType()
    {
      CheckInitSetup();
      return masslintype_;
    };
    //@}

    //! @name Get model evaluator control parameters (read and write access)
    //@{
    //! Return types of the current models
    //! Note: It is possible to add or remove a model type, but you are not allowed to
    //!       change an already existing model type
    std::vector<const enum INPAR::STR::ModelType>& GetMutableModelTypes()
    {
      CheckInitSetup();
      return *modeltypes_;
    };

    //! Return map of the current models
    std::map<const enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> >& GetMutableModelEvaluators()
    {
      CheckInitSetup();
      return modelevaluators_;
    };
    //@}

    //! @name Get implicit time integrator parameters (read and write access)
    //@{
    //! Return specific implicit time integrator pointer
    //! NOTE: Only filled for the implicit case.
    STR::IMPLICIT::Generic& GetMutableImplicitIntegrator()
    {
      CheckInitSetup();
      if (implint_.is_null())
        dserror("Implicit Integrator is not initialized!");
      return *implint_;
    };

    //! non-linear solver pointer
    //! NOTE: Only filled for the implicit case
    STR::NLN::SOLVER::Generic& GetMutableNlnSolvers()
    {
      CheckInitSetup();
      if (nlnsolvers_.is_null())
        dserror("Non-linear solvers are not initialized!");
      return *nlnsolvers_;
    };

    //@}

  protected:
    //! @name variables for internal use only
    //@{

    bool isInit_;

    bool isSetup_;
    //@}

  private:

    //! @name General control parameters
    //@{

    //! final time \f$t_\text{fin}\f$
    double timemax_;

    //! final time step \f$N\f$
    int stepmax_;

    //@}

    //! @name Time measurement
    //@{
    //! timer for solution technique
    Teuchos::RCP<Epetra_Time> timer_;

    //@}

    //! @name Damping control parameters
    //! Rayleigh damping means \f${C} = c_\text{K} {K} + c_\text{M} {M}\f$
    //@{
    //! damping type
    enum INPAR::STR::DampKind damptype_;

    //! damping factor for stiffness \f$c_\text{K}\f$
    double dampk_;

    //! damping factor for mass \f$c_\text{M}\f$
    double dampm_;
    //@}

    //! @name Mass and inertia control parameters
    //@{
    //! have inertia forces to be linearized?
    enum INPAR::STR::MassLin masslintype_;
    //@}

    //! @name Model evaluator control parameters
    //@{

    //! current active model types for the model evaluator
    Teuchos::RCP<const std::vector<const enum INPAR::STR::ModelType> > modeltypes_;

    //! model evaluators corresponding to the model types (see variable modeltypes_)
    //! NOTE: This variable is initialized in the ex-/implicit time integrator function!
    std::map<const enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> > modelevaluators_;
    //@}

    //! @name Linear solver parameters
    //@{

    //! map of all necessary linear solvers (corresponding to the active model types)
    Teuchos::RCP<std::map<const enum INPAR::STR::ModelType,Teuchos::RCP<LINALG::Solver> > > linsolvers_;
    //@}

    //! @name implicit and explicit time integrator parameters
    //@{
    //! dynamic type
    enum INPAR::STR::DynamicType dynType_;

    //@}

    //! @name implicit time integrator parameters
    //@{
    //! type of pre-stressing
    enum INPAR::STR::PreStress preStressType_;

    //! type of nonlinear solver
    enum INPAR::STR::NonlinSolTech nlnSolverType_;

    //! specific implicit time integrator pointer
    //! NOTE: Only filled for the implicit case.
    Teuchos::RCP<STR::IMPLICIT::Generic> implint_;

    //! non-linear solver pointer
    //! NOTE: Only filled for the implicit case
    Teuchos::RCP<STR::NLN::SOLVER::Generic> nlnsolvers_;

    //@}
  };  // BaseDataSDyn


  class BaseDataGlobalState
    {
    public:
      BaseDataGlobalState();

      virtual ~BaseDataGlobalState() {};

      void Init(const Teuchos::RCP<DRT::Discretization> discret,
                const Teuchos::RCP<const BaseDataSDyn> dataSDyn);

      void Setup();

    public:
      //! @name Get general purpose algorithm members (read only access)
      //@{

      //! attached discretisation
      Teuchos::RCP<const DRT::Discretization> GetDiscret()
      { return discret_; };

      //! communicator
      Teuchos::RCP<const Epetra_Comm> GetComm()
      { return comm_; };

      //! ID of actual processor in parallel
      const int GetMyRank()
      { return myRank_; };

      //@}

      //! @name Get general control parameters (read only access)
      //@{

      //! Return target time \f$t_{n+1}\f$
      const double& GetTimeNp() const
      { return timenp_; };

      //! Return time for lin solver
      const double& GetLinearSolverTime() const
      { return dtsolve_; };

      //! Return element evaluation time
      const double& GetElementEvaluationTime() const
      { return dtele_; };

      //! Return time \f$t_{n}\f$ of last converged step
      Teuchos::RCP<const DRT::UTILS::TimIntMStep<double> > GetTimeN() const
      { return timen_; };

      //! Return time step size \f$\Delta t\f$
      Teuchos::RCP<const DRT::UTILS::TimIntMStep<double> > GetDeltaTime() const
      { return dt_; };

      //! Return timer for solution technique
      Teuchos::RCP<const Epetra_Time>  GetTimer() const
      { return timer_; };

      //@}

      //! @name Get state variables (read only access)
      //@{

      //! Return displacements \f$D_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> GetDisNp() const
      { return disnp_; }

      //! Return displacements \f$D_{n}\f$
      Teuchos::RCP<const Epetra_Vector> GetDisN() const
      { return (*dis_)(0); }

      //! Return velocities \f$V_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> GetVelNp() const
      { return velnp_; }

      //! Return velocities \f$V_{n}\f$
      Teuchos::RCP<const Epetra_Vector> GetVelN() const
      { return (*vel_)(0); }

      //! Return accelerations \f$A_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> GetAccNp() const
      { return accnp_; }

      //! Return accelerations \f$A_{n}\f$
      Teuchos::RCP<const Epetra_Vector> GetAccN() const
      { return (*acc_)(0); }

      //@}

      //! @name Get mutable general control parameters (read and write access)
      //@{

      //! Return target time \f$t_{n+1}\f$
      double& GetMutableTimeNp()
      { return timenp_; };

      //! Return time for lin solver
      double& GetMutableLinearSolverTime()
      { return dtsolve_; };

      //! Return element evaluation time
      double& GetMutableElementEvaluationTime()
      { return dtele_; };

      //! Return time \f$t_{n}\f$ of last converged step
      Teuchos::RCP<DRT::UTILS::TimIntMStep<double> >& GetMutableTimeN()
      { return timen_; };

      //! Return time step size \f$\Delta t\f$
      Teuchos::RCP<DRT::UTILS::TimIntMStep<double> >& GetMutableDeltaTime()
      { return dt_; };

      //! Return timer for solution technique
      Teuchos::RCP<Epetra_Time>& GetMutableTimer()
      { return timer_; };

      //@}

      //! @name Get mutable state variables (read and write access)
      //@{

      //! Return displacements \f$D_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& GetMutableDisNp()
      { return disnp_; }

      //! Return displacements \f$D_{n}\f$
      Teuchos::RCP<Epetra_Vector> GetMutableDisN()
      { return (*dis_)(0); }

      //! Return velocities \f$V_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& GetMutableVelNp()
      { return velnp_; }

      //! Return velocities \f$V_{n}\f$
      Teuchos::RCP<Epetra_Vector> GetMutableVelN()
      { return (*vel_)(0); }

      //! Return accelerations \f$A_{n+1}\f$
      Teuchos::RCP<Epetra_Vector>& GetMutableAccNp()
      { return accnp_; }

      //! Return accelerations \f$A_{n}\f$
      Teuchos::RCP<Epetra_Vector> GetMutableAccN()
      { return (*acc_)(0); }

      //@}

    protected:
      const bool& IsInit();

      const bool& IsSetup();

    protected:
      //! @name variables for internal use only
      //@{

      bool isInit_;

      bool isSetup_;

      // read only access
      Teuchos::RCP<const BaseDataSDyn> dataSDyn_;
      //@}

    private:
      //! @name General purpose algorithm members
      //@{

      //! attached discretisation
      Teuchos::RCP<DRT::Discretization> discret_;

      //! communicator
      Teuchos::RCP<const Epetra_Comm> comm_;

      //! ID of actual processor in parallel
      int myRank_;

      //@}

      //! @name General control parameters
      //@{
      //! target time \f$t_{n+1}\f$
      double timenp_;

      //! time \f$t_{n}\f$ of last converged step
      Teuchos::RCP<DRT::UTILS::TimIntMStep<double> > timen_;

      //! time step size \f$\Delta t\f$
      Teuchos::RCP<DRT::UTILS::TimIntMStep<double> > dt_;

      //! time step index \f$n\f$
      int stepn_;

      //! time step index \f$n+1\f$
      int stepnp_;
      //@}

      //! @name Global state vectors
      //@{

      //! global displacements \f${D}_{n}, D_{n-1}, ...\f$
      Teuchos::RCP<DRT::UTILS::TimIntMStep<Epetra_Vector> > dis_;

      //! global velocities \f${V}_{n}, V_{n-1}, ...\f$
      Teuchos::RCP<DRT::UTILS::TimIntMStep<Epetra_Vector> > vel_;

      //! global accelerations \f${A}_{n}, A_{n-1}, ...\f$
      Teuchos::RCP<DRT::UTILS::TimIntMStep<Epetra_Vector> > acc_;

      //! global displacements \f${D}_{n+1}\f$ at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> disnp_;

      //! global velocities \f${V}_{n+1}\f$ at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> velnp_;

      //! global accelerations \f${A}_{n+1}\f$ at \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> accnp_;

      //! @name System matrices
      //@{
      //! holds eventually effective stiffness TODO
      Teuchos::RCP<LINALG::SparseOperator> stiff_;

      //! mass matrix (constant)
      Teuchos::RCP<LINALG::SparseOperator> mass_;

      //! damping matrix
      Teuchos::RCP<LINALG::SparseOperator> damp_;
      //@}

      //! @name Time measurement
      //@{
      //! timer for solution technique
      Teuchos::RCP<Epetra_Time> timer_;

      //! linear solver time
      double dtsolve_;

      //! element evaluation time
      double dtele_;
      //@}
    };  // BaseDataGlobalState


  class Base : public ::ADAPTER::StructureNew
  {
  public:
    //! constructor
    Base();

    //! destructor
    virtual ~Base() {};

    virtual void IntegrateStep() = 0;

    /// wrapper for things that should be done before PrepareTimeStep is called
    virtual void PrePredict() {};

    /// wrapper for things that should be done after PrepareTimeStep is called
    virtual void PostPredict() {};

    /// wrapper for things that should be done before solving the nonlinear iterations
    virtual void PreSolve() {};

    /// wrapper for things that should be done after solving the nonlinear iterations
    virtual void PostSolve();

    /// wrapper for things that should be done before updating
    virtual void PreUpdate() {};

    /// wrapper for things that should be done after solving the update
    virtual void PostUpdate() {};

    /// wrapper for things that should be done before the output
    virtual void PreOutput() {};

    /// wrapper for things that should be done after convergence of Newton scheme (like crack propagation)
    virtual void PostOutput() {};

    // Get TimIntBase data for global state quantities (read access)
    Teuchos::RCP<const BaseDataGlobalState> DataGlobalState() const;

    // Get TimIntBase data for io quantities (read access)
    Teuchos::RCP<const BaseDataIO> DataIO() const;

    // Get TimIntBase data or struct dynamics quantitites (read access)
    Teuchos::RCP<const BaseDataSDyn> DataSDyn() const;

    virtual void Init(
        const Teuchos::ParameterList& ioparams,
        const Teuchos::ParameterList& sdynparams,
        const Teuchos::ParameterList& xparams,
        const Teuchos::RCP<DRT::Discretization> discret,
        const std::vector<const enum INPAR::STR::ModelType>& modeltypes,
        std::map<const enum INPAR::STR::ModelType, Teuchos::RCP<LINALG::Solver> >& linsolvers,
        Teuchos::RCP<IO::DiscretizationWriter> output
    );

    virtual void Setup() = 0;

  protected:

    const bool& IsInit() const { return isInit_; };

    const bool& IsSetup() const { return isSetup_; };

    // Get internal TimIntBase data for global state quantities (read and write access)
    BaseDataGlobalState& DataGlobalState()
    { return *dataGlobalState_; };

    // Get internal TimIntBase data for io quantities (read and write access)
    BaseDataIO& DataIO()
    { return *dataIO_; };

    // Get internal TimIntBase data for struct dynamics quantitites (read and write access)
    BaseDataSDyn& DataSDyn()
    { return *dataSDyn_; };

  protected:
    bool isInit_;

    bool isSetup_;

  private:

    Teuchos::RCP<BaseDataIO>          dataIO_;
    Teuchos::RCP<BaseDataSDyn>        dataSDyn_;
    Teuchos::RCP<BaseDataGlobalState> dataGlobalState_;
  };  // class Base
} // namespace TIMINT
} // namespace STR


#endif /* STR_TIMINT_BASE_H_ */

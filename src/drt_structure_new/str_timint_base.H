/*-----------------------------------------------------------*/
/*!

\brief Base class for all structural time integration strategies.

\maintainer Anh-Tu Vuong

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_TIMINT_BASE_H_
#define STR_TIMINT_BASE_H_

#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_io/every_iteration_writer.H"

// structural data containers
#include "str_timint_basedataio.H"
#include "str_timint_basedataglobalstate.H"
#include "str_timint_basedatasdyn.H"


// forward declaration
class Epetra_Vector;
class Epetra_Map;

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos
namespace LINALG
{
  class BlockSparseMatrixBase;
}  // namespace LINALG
namespace STR
{
  class ModelEvaluator;
  class Dbc;
  class Integrator;
  namespace MODELEVALUATOR
  {
    class Generic;
  }  // namespace MODELEVALUATOR
  namespace TIMINT
  {
    /** \brief Abstract class for all time integration strategies
     *
     *  \author Michael Hiermeier */
    class Base : public ::ADAPTER::StructureNew, IO::EveryIterationWriterInterface
    {
     public:
      /// constructor
      Base();

      /// destructor
      virtual ~Base(){};

      /// initialize (all already existing) class variables
      virtual void Init(const Teuchos::RCP<STR::TIMINT::BaseDataIO> dataio,
          const Teuchos::RCP<STR::TIMINT::BaseDataSDyn> datasdyn,
          const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> dataglobalstate);

      /// setup of the new class variables
      virtual void Setup();

      /// tests if there are more time steps to do
      virtual bool NotFinished() const;

      /// reset everything (needed for biofilm simulations)
      virtual void Reset();

      /** \brief reset step configuration after time step
       *
       *  This function is supposed to reset all variables which are directly related
       *  to the current new step n+1. To be more precise all variables ending with "Np"
       *  have to be reseted. */
      virtual void ResetStep();

      /// wrapper for things that should be done before PrepareTimeStep is called
      virtual void PrePredict(){};

      /// wrapper for things that should be done after PrepareTimeStep is called
      virtual void PostPredict(){};

      /// wrapper for things that should be done before solving the nonlinear iterations
      virtual void PreSolve(){};

      /// wrapper for things that should be done after solving the nonlinear iterations
      virtual void PostSolve(){};

      /// wrapper for things that should be done before the output
      virtual void PreOutput(){};

      /// wrapper for things that should be done after convergence of Newton scheme
      virtual void PostOutput(){};

      /// things that should be done after the actual time loop is finished
      virtual void PostTimeLoop(){/* FixMe HaveMicroMat() is doing something here */};

      /// @name General access methods
      ///@{
      /// Access discretization (structure only)
      virtual Teuchos::RCP<DRT::Discretization> Discretization();

      /// Access generic discretization pointer
      virtual Teuchos::RCP<DRT::DiscretizationInterface> DiscretizationInterface();

      /// Access to pointer to DoF row map of the discretization (structure only)
      virtual const Epetra_Map* DofRowMapView()
      {
        CheckInit();
        return dataglobalstate_->DofRowMapView();
      }

      /// DoF map of structural vector of unknowns
      virtual Teuchos::RCP<const Epetra_Map> DofRowMap()
      {
        CheckInit();
        return dataglobalstate_->DofRowMap();
      }

      //! DoF map of vector of unknowns
      // Alternative method capable of multiple DoF sets
      virtual Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds)
      {
        CheckInit();
        return dataglobalstate_->DofRowMap(nds);
      }

      /// Access linear structural solver
      virtual Teuchos::RCP<LINALG::Solver> LinearSolver()
      {
        CheckInit();
        return datasdyn_->GetMutableLinSolvers()[INPAR::STR::model_structure];
      }

      /// Return MapExtractor for Dirichlet boundary conditions
      virtual const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor();
      Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() const;

      //! Return locsys manager
      virtual Teuchos::RCP<DRT::UTILS::LocsysManager> LocsysManager();

      //! Return the desired model evaluator (read-only)
      virtual const STR::MODELEVALUATOR::Generic& ModelEvaluator(INPAR::STR::ModelType mtype) const;

      //! Return the desired model evaluator (read and write)
      virtual STR::MODELEVALUATOR::Generic& ModelEvaluator(INPAR::STR::ModelType mtype);

      ///@}

      /// Return domain map of the mass matrix (implicit and explicit)
      virtual const Epetra_Map& GetMassDomainMap() const;

      /// @name Coupled problem routines
      /// @{
      /// wrapper for things that should be done before updating
      virtual void PreUpdate(){};

      /// Update routine for coupled problems with monolithic approach
      virtual void Update();

      /// Update time and step counter
      virtual void UpdateStepTime();

      /// wrapper for things that should be done after solving the update
      virtual void PostUpdate(){};
      /// @}

      /// @name Access global state from outside via adapter (needed for coupled problems)
      ///@{
      /// unknown displacements at \f$t_{n+1}\f$
      virtual Teuchos::RCP<const Epetra_Vector> DispNp() const
      {
        CheckInit();
        return dataglobalstate_->GetDisNp();
      }

      /* \brief write access to displacements at \f$t^{n+1}\f$
       *
       * Calling this method makes only sense if state is supposed
       * to be manipulated. We must not forget to synchronize the
       * manipulated state with the NOX group.
       * Otherwise, the manipulations will be overwritten by NOX.
       * Therefore, we set the flag state_is_insync_with_noxgroup_
       * to false.
       * This will be checked:
       * See \ref CheckStateInSyncWithNOXGroup
       *
       * See also \ref ADAPTER::StructureNew::SetState
       */
      virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispNp()
      {
        CheckInit();
        SetStateInSyncWithNOXGroup(false);
        return dataglobalstate_->GetMutableDisNp();
      }

      /// known displacements at \f$t_{n}\f$
      virtual Teuchos::RCP<const Epetra_Vector> DispN() const
      {
        CheckInit();
        return dataglobalstate_->GetDisN();
      }

      /// write access to displacements at \f$t^{n}\f$
      virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispN()
      {
        CheckInit();
        return dataglobalstate_->GetMutableDisN();
      }

      /// unknown velocities at \f$t_{n+1}\f$
      virtual Teuchos::RCP<const Epetra_Vector> VelNp() const
      {
        CheckInit();
        return dataglobalstate_->GetVelNp();
      }

      /// write access to velocities at \f$t^{n+1}\f$
      virtual Teuchos::RCP<Epetra_Vector> WriteAccessVelNp()
      {
        CheckInit();
        return dataglobalstate_->GetMutableVelNp();
      }

      /// unknown velocities at \f$t_{n}\f$
      virtual Teuchos::RCP<const Epetra_Vector> VelN() const
      {
        CheckInit();
        return dataglobalstate_->GetVelN();
      }

      /// write access to velocities at \f$t^{n}\f$
      virtual Teuchos::RCP<Epetra_Vector> WriteAccessVelN()
      {
        CheckInit();
        return dataglobalstate_->GetMutableVelN();
      }

      /// known velocities at \f$t_{n-1}\f$
      virtual Teuchos::RCP<const Epetra_Vector> VelNm() const
      {
        CheckInit();
        return dataglobalstate_->GetVelNm();
      }

      /// unknown accelerations at \f$t_{n+1}\f$
      virtual Teuchos::RCP<const Epetra_Vector> AccNp() const
      {
        CheckInit();
        return dataglobalstate_->GetAccNp();
      }

      //! known accelerations at \f$t_{n}\f$
      virtual Teuchos::RCP<const Epetra_Vector> AccN() const
      {
        CheckInit();
        return dataglobalstate_->GetAccN();
      }
      ///@}

      /// @name access and modify model evaluator stuff via adapter
      /// @{
      /// are there any algebraic constraints?
      virtual bool HaveConstraint()
      {
        CheckInitSetup();
        return datasdyn_->HaveModelType(INPAR::STR::model_lag_pen_constraint);
      }

      /// do we need a semi-smooth Newton-type plasticity algorithm
      virtual bool HaveSemiSmoothPlasticity()
      {
        CheckInitSetup();
        return datasdyn_->HaveEleTech(INPAR::STR::eletech_plasticity);
      }

      /// FixMe get constraint manager defined in the structure
      virtual const Teuchos::RCP<UTILS::ConstrManager> GetConstraintManager()
      {
        dserror("Not yet implemented!");
        return Teuchos::null;
      }

      /// FixMe get contact/meshtying manager
      virtual Teuchos::RCP<CONTACT::MeshtyingContactBridge> MeshtyingContactBridge()
      {
        dserror("Not yet implemented!");
        return Teuchos::null;
      }

      /// do we have this model
      virtual bool HaveModel(INPAR::STR::ModelType model)
      {
        return datasdyn_->HaveModelType(model);
      }

      /// Add residual increment to Lagrange multipliers stored in Constraint manager (derived)
      /// FixMe Different behavior for the implicit and explicit case!!!
      virtual void UpdateIterIncrConstr(Teuchos::RCP<Epetra_Vector> lagrincr)
      {
        dserror("Not yet implemented!");
        return;
      }

      /// Add residual increment to pressures stored in Cardiovascular0D manager (derived)
      /// FixMe Different behavior for the implicit and explicit case!!!
      virtual void UpdateIterIncrCardiovascular0D(Teuchos::RCP<Epetra_Vector> presincr)
      {
        dserror("Not yet implemented!");
        return;
      }
      /// @}

      /// @name Time step helpers
      ///@{
      /// Return current time \f$t_{n}\f$ (derived)
      virtual double GetTimeN() const
      {
        CheckInit();
        return dataglobalstate_->GetTimeN();
      }

      /// Sets the current time \f$t_{n}\f$ (derived)
      virtual void SetTimeN(const double time_n)
      {
        CheckInit();
        dataglobalstate_->GetMutableTimeN() = time_n;
      }

      /// Return target time \f$t_{n+1}\f$ (derived)
      virtual double GetTimeNp() const
      {
        CheckInit();
        return dataglobalstate_->GetTimeNp();
      }

      /// Sets the target time \f$t_{n+1}\f$ of this time step (derived)
      virtual void SetTimeNp(const double time_np)
      {
        CheckInit();
        dataglobalstate_->GetMutableTimeNp() = time_np;
      }

      /// Get upper limit of time range of interest (derived)
      virtual double GetTimeEnd() const
      {
        CheckInit();
        return datasdyn_->GetTimeMax();
      }

      /// Get upper limit of time range of interest (derived)
      virtual void SetTimeEnd(double timemax)
      {
        CheckInit();
        datasdyn_->GetMutableTimeMax() = timemax;
      }

      /// Get time step size \f$\Delta t_n\f$
      virtual double GetDeltaTime() const
      {
        CheckInit();
        return (*dataglobalstate_->GetDeltaTime())[0];
      }

      /// Set time step size \f$\Delta t_n\f$
      virtual void SetDeltaTime(const double dt)
      {
        CheckInit();
        (*dataglobalstate_->GetMutableDeltaTime())[0] = dt;
      }

      /// Return time integration factor
      virtual double TimIntParam() const;

      /// Return current step number \f$n\f$
      virtual int GetStepN() const
      {
        CheckInit();
        return dataglobalstate_->GetStepN();
      }

      /// Sets the current step \f$n\f$
      virtual void SetStepN(int step_n)
      {
        CheckInit();
        dataglobalstate_->GetMutableStepN() = step_n;
      }

      /// Return current step number $n+1$
      virtual int GetStepNp() const
      {
        CheckInit();
        return dataglobalstate_->GetStepNp();
      }

      /// Sets the current step number \f$n+1\f$
      virtual void SetStepNp(int step_np)
      {
        CheckInitSetup();
        dataglobalstate_->GetMutableStepNp() = step_np;
      }

      //! Get number of time steps
      virtual int GetStepEnd() const
      {
        CheckInit();
        return datasdyn_->GetStepMax();
      }

      /// Sets number of time steps
      virtual void SetStepEnd(int step_end)
      {
        CheckInitSetup();
        datasdyn_->GetMutableStepMax() = step_end;
      }

      //! Get divcont type
      virtual enum INPAR::STR::DivContAct GetDivergenceAction() const
      {
        CheckInitSetup();
        return datasdyn_->GetDivergenceAction();
      }

      //! Get number of times you want to halve your timestep in case nonlinear solver diverges
      virtual int GetMaxDivConRefineLevel() const
      {
        CheckInitSetup();
        return datasdyn_->GetMaxDivConRefineLevel();
      }

      //! Get random factor for time step adaption
      virtual double GetRandomTimeStepFactor() const
      {
        CheckInitSetup();
        return datasdyn_->GetRandomTimeStepFactor();
      }

      //! Set random factor for time step adaption
      virtual double SetRandomTimeStepFactor(double rand_tsfac)
      {
        CheckInitSetup();
        return datasdyn_->GetMutableRandomTimeStepFactor() = rand_tsfac;
      }

      //! Get random factor for time step adaption
      virtual int GetDivConRefineLevel() const
      {
        CheckInitSetup();
        return datasdyn_->GetDivConRefineLevel();
      }

      //! Set random factor for time step adaption
      virtual int SetDivConRefineLevel(int divconrefinementlevel)
      {
        CheckInitSetup();
        return datasdyn_->GetMutableDivConRefineLevel() = divconrefinementlevel;
      }

      //! Get random factor for time step adaption
      virtual int GetDivConNumFineStep() const
      {
        CheckInitSetup();
        return datasdyn_->GetDivConNumFineStep();
      }

      //! Set random factor for time step adaption
      virtual int SetDivConNumFineStep(int divconnumfinestep)
      {
        CheckInitSetup();
        return datasdyn_->GetMutableDivConNumFineStep() = divconnumfinestep;
      }

      /// set evaluation action
      virtual void SetActionType(const DRT::ELEMENTS::ActionType& action);

      // group id in nested parallelity
      int GroupId() const;
      ///@}

      /// @name Structure with ale specific methods
      ///@{
      /// FixMe set/apply material displacements to structure field (structure with ale)
      virtual void SetDispMatNp(Teuchos::RCP<Epetra_Vector> dispmatnp)
      {
        dserror("Not supported at the moment!");
      }

      /// FixMe write access to material displacements (strutcure with ale) at \f$t^{n+1}\f$
      virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispMatNp()
      {
        CheckInitSetup();
        dserror("Not yet supported!");
        return Teuchos::null;
      }
      ///@}


      /// Time adaptivity (derived pure virtual functionality)
      /// @{
      /// Resize MStep Object due to time adaptivity in FSI (derived)
      virtual void ResizeMStepTimAda();

      /// @}

      /// Output writer related routines (file and screen output)
      /// @{
      /// Access output object
      virtual Teuchos::RCP<IO::DiscretizationWriter> DiscWriter()
      {
        return DataIO().GetMutableOutputPtr();
      }

      /// Calculate all output quantities depending on the constitutive model
      /// (and, hence, on a potential material history)
      virtual void PrepareOutput();

      /// output results (implicit and explicit)
      virtual void Output() { Output(false); };
      virtual void Output(bool forced_writerestart);

      /// Write Gmsh output for structural field
      virtual void writeGmshStrucOutputStep();

      /// FixMe Check if there are any elements with the micro material definition.
      /// Maybe the detection can be moved to the element loop in the ad_str_structure_new.cpp.
      /// There is already one.
      virtual bool HaveMicroMat()
      {
        dserror("Not yet considered!");
        return false;
      };

      /// create result test for encapsulated structure algorithm
      virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

      /** \brief Get data that is written during restart
       *
       *  This routine is only for simple structure problems!
       *  \date 06/13
       *  \author biehler */
      virtual void GetRestartData(Teuchos::RCP<int> step, Teuchos::RCP<double> time,
          Teuchos::RCP<Epetra_Vector> disnp, Teuchos::RCP<Epetra_Vector> velnp,
          Teuchos::RCP<Epetra_Vector> accnp, Teuchos::RCP<std::vector<char>> elementdata,
          Teuchos::RCP<std::vector<char>> nodedata);

      /** Read restart values
       *
       * \param stepn (in): restart step at \f${n}\f$
       */
      virtual void ReadRestart(const int stepn);

      /// Set restart values (deprecated)
      virtual void SetRestart(int stepn,                //!< restart step at \f${n}\f$
          double timen,                                 //!< restart time at \f$t_{n}\f$
          Teuchos::RCP<Epetra_Vector> disn,             //!< restart displacements at \f$t_{n}\f$
          Teuchos::RCP<Epetra_Vector> veln,             //!< restart velocities at \f$t_{n}\f$
          Teuchos::RCP<Epetra_Vector> accn,             //!< restart accelerations at \f$t_{n}\f$
          Teuchos::RCP<std::vector<char>> elementdata,  //!< restart element data
          Teuchos::RCP<std::vector<char>> nodedata      //!< restart element data
      );
      /// @}

      /// Biofilm related stuff
      /// @{
      /// FixMe set structure displacement vector due to biofilm growth
      virtual void SetStrGrDisp(Teuchos::RCP<Epetra_Vector> struct_growth_disp)
      {
        dserror("Currently unsupported!");
      }
      /// @}

      /// @name Pure virtual adapter functions (have to be implemented in the derived classes)
      /// @{
      /// Take the time and integrate (time loop) (implicit and explicit)
      virtual int Integrate() = 0;
      /// integrate the current step (implicit and explicit)
      virtual int IntegrateStep() = 0;
      /// nonlinear solve (implicit and explicit)
      virtual INPAR::STR::ConvergenceStatus Solve() = 0;
      //! print summary after step (implicit and explicit)
      virtual void PrintStep() = 0;
      /// start new time step (implicit and explicit)
      virtual void PrepareTimeStep() = 0;
      /// specific method for iterative staggered partitioned TSI (implicit and explicit)
      virtual void PreparePartitionStep() = 0;
      /// set the state of the nox group and the global state data container
      /// see class \ref ADAPTER::StructureNew for detailed documentation.
      /// \note implemented only for implicit
      virtual void SetState(const Teuchos::RCP<Epetra_Vector>& x) = 0;
      /// update displacement and evaluate elements (implicit only)
      virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc) = 0;
      /// don't update displacement but evaluate elements (implicit only)
      virtual void Evaluate() = 0;
      /// Update routine for coupled problems with monolithic approach with time adaptivity
      /// (implicit only)
      virtual void Update(double endtime) = 0;
      /// Get type of thickness scaling for thin shell structures (implicit only)
      virtual INPAR::STR::STC_Scale GetSTCAlgo() = 0;
      /// Access to scaling matrix for STC (implicit only)
      virtual Teuchos::RCP<LINALG::SparseMatrix> GetSTCMat() = 0;
      /// initial guess of Newton's method (implicit only)
      virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() = 0;
      /// right-hand-side of Newton's method (implicit only)
      virtual Teuchos::RCP<const Epetra_Vector> RHS() { return GetF(); };
      virtual Teuchos::RCP<const Epetra_Vector> GetF() const = 0;
      /// Return reaction forces at \f$t_{n+1}\f$ (read and write) (implicit only)
      virtual Teuchos::RCP<Epetra_Vector> Freact() = 0;
      //! Return stiffness (implicit only)
      virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() = 0;
      /// Return stiffness (implicit only)
      virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() = 0;
      ///! Switch structure field to block matrix in fsi simulations (implicit only)
      virtual void UseBlockMatrix(Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps,
          Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps) = 0;
      /// @}

      /// @name XFEM/XSTR related methods
      ///@{

      /// access the xfield state pointer
      virtual Teuchos::RCP<XFEM::XFieldState> XFieldState();
      virtual Teuchos::RCP<const XFEM::XFieldState> XFieldState() const;

      /** \brief create a new XFieldState object by transferring the information to the new dofs
       *
       *  \author hiermeier \date 02/17 */
      virtual void CreateNewXFieldState(const Teuchos::RCP<XFEM::XFieldState>& new_xstate) const;

      /** \brief reset the non-standard / enrichment dofs of the xfield state object
       *
       *  \author hiermeier \date 02/17 */
      virtual void ResetXFieldNonStandardDofs();

      /** \brief destroy the current state
       *
       *  \author hiermeier \date 02/17 */
      virtual bool DestroyState();

      /** \brief destroy the nonlinear solver state
       *
       *  \author hiermeier \date 02/17 */
      virtual bool DestroyNoxState() = 0;

      ///@}

     public:
      /// @name External accessors for the class variables
      ///@{
      /// Get the indicator if we are currently restarting the simulation
      inline const bool& IsRestarting() const { return isrestarting_; };

      /// Get TimIntBase data for global state quantities (read access)
      Teuchos::RCP<const BaseDataGlobalState> GetDataGlobalStatePtr() const
      {
        CheckInit();
        return dataglobalstate_;
      };

      /// Get TimIntBase data for global state quantities (read & write access)
      Teuchos::RCP<BaseDataGlobalState>& GetMutableDataGlobalStatePtr()
      {
        CheckInit();
        return dataglobalstate_;
      };

      const BaseDataGlobalState& GetDataGlobalState() const
      {
        CheckInit();
        return *dataglobalstate_;
      };

      /// Get TimIntBase data for io quantities (read access)
      Teuchos::RCP<const BaseDataIO> GetDataIOPtr() const
      {
        CheckInit();
        return dataio_;
      }

      const BaseDataIO& GetDataIO() const
      {
        CheckInit();
        return *dataio_;
      }

      /// Get TimIntBase data or struct dynamics quantitites (read access)
      Teuchos::RCP<const BaseDataSDyn> GetDataSDynPtr() const
      {
        CheckInit();
        return datasdyn_;
      }

      const BaseDataSDyn& GetDataSDyn() const
      {
        CheckInit();
        return *datasdyn_;
      }

      /// return a reference to the Dirichlet Boundary Condition handler (read access)
      const STR::Dbc& GetDBC() const
      {
        CheckInitSetup();
        return *dbc_ptr_;
      }

      /// return a reference to the Dirichlet Boundary Condition handler (write access)
      STR::Dbc& GetMutableDBC()
      {
        CheckInitSetup();
        return *dbc_ptr_;
      }

      /// return a pointer to the Dirichlet Boundary Condition handler (read access)
      Teuchos::RCP<const STR::Dbc> GetDBCPtr() const
      {
        CheckInitSetup();
        return dbc_ptr_;
      }

      /// return the integrator (read-only)
      const STR::Integrator& Integrator() const
      {
        CheckInitSetup();
        return *int_ptr_;
      }

      /// Get internal TimIntBase data for structural dynamics quantities (read and write access)
      BaseDataSDyn& DataSDyn()
      {
        CheckInit();
        return *datasdyn_;
      };

      /// return a pointer to the Dirichlet Boundary Condition handler (read and write access)
      const Teuchos::RCP<STR::Dbc>& DBCPtr()
      {
        CheckInitSetup();
        return dbc_ptr_;
      }

      ///@}
     protected:
      /// get the indicator state
      inline const bool& IsInit() const { return isinit_; };

      /// get the indicator state
      inline const bool& IsSetup() const { return issetup_; };

      /// Check if Init() and Setup() have been called, yet.
      inline void CheckInitSetup() const
      {
        if (!IsInit() or !IsSetup()) dserror("Call Init() and Setup() first!");
      }

      /// Check if Init() has been called
      inline void CheckInit() const
      {
        if (not IsInit()) dserror("Call Init() first!");
      }

      BaseDataGlobalState& DataGlobalState()
      {
        CheckInit();
        return *dataglobalstate_;
      };

      const Teuchos::RCP<BaseDataGlobalState>& DataGlobalStatePtr()
      {
        CheckInit();
        return dataglobalstate_;
      };

      /// Get internal TimIntBase data for io quantities (read and write access)
      BaseDataIO& DataIO()
      {
        CheckInit();
        return *dataio_;
      };

      /// return a pointer to the input/output data container (read and write access)
      const Teuchos::RCP<BaseDataIO>& DataIOPtr()
      {
        CheckInit();
        return dataio_;
      };

      /// return a pointer to the structural dynamic data container (read and write access)
      const Teuchos::RCP<BaseDataSDyn>& DataSDynPtr()
      {
        CheckInit();
        return datasdyn_;
      };

      /// return a reference to the Dirichlet Boundary Condition handler (read and write access)
      STR::Dbc& DBC()
      {
        CheckInitSetup();
        return *dbc_ptr_;
      }

      /// return a reference to the integrator (read and write access)
      STR::Integrator& Integrator()
      {
        CheckInitSetup();
        return *int_ptr_;
      }

      /// return a pointer to the integrator (read and write access)
      const Teuchos::RCP<STR::Integrator>& IntegratorPtr()
      {
        CheckInitSetup();
        return int_ptr_;
      }

      /** \brief Output to file
       *
       *  This routine always prints the last converged state, i.e.
       *  \f$D_{n}, V_{n}, A_{n}\f$.
       *
       *  \date 03/07
       *  \author mwgee (originally) */
      void OutputStep(bool forced_writerestart);

     private:
      /*! \brief Create a new input/output step in the output writer
       *
       * New step is created only once per time step. This is controlled by \c datawritten.
       * Do nothing if data has already been written in this time step.
       *
       * \param[in,out] Indicator whether data has already been written in this time step (true) or
       *                not (false)
       */
      void NewIOStep(bool& datawritten);

      /// output of the current state
      void OutputState();

      /** \brief output of the current state */
      void OutputState(IO::DiscretizationWriter& iowriter, bool write_owner) const;

      /** \brief output of the debug state */
      void OutputDebugState(IO::DiscretizationWriter& iowriter, bool write_owner) const override;

      /// output during runtime
      void RuntimeOutputState();

      /// output element volumes
      void OutputElementVolume(IO::DiscretizationWriter& iowriter) const;

      /// output stress and/or strain state
      void OutputStressStrain();

      /// output energy
      void OutputEnergy() const;

      /// output optional quantity
      void OutputOptionalQuantity();

      /// write restart information
      void OutputRestart(bool& datawritten);

      /// add restart information to output state
      void AddRestartToOutputState();

      /** \brief set the number of nonlinear iterations of the last time step
       *
       *  \pre UpdateStepTime() must be called beforehand, otherwise the wrong
       *  step-id is considered.
       *
       *  \author hiermeier \date 11/17 */
      void SetNumberOfNonlinearIterations();

      /** \brief decide which contributions to the total system energy shall be
       *         computed and written during simulation
       *
       *  \author grill */
      void SelectEnergyTypesToBeWritten();

      /** \brief initialize file stream for energy values and write all the
       *         column headers for the previously selected energy contributions
       *         to be written separately
       *
       *  \author grill */
      void InitializeEnergyFileStreamAndWriteHeaders();

     protected:
      /// flag indicating if Init() has been called
      bool isinit_;

      /// flag indicating if Setup() has been called
      bool issetup_;

      /// flag indicating that the simulation is currently restarting
      bool isrestarting_;

      /// flag indicating that displacement state was manipulated
      /// but NOX group has not been informed.
      bool state_is_insync_with_noxgroup_;

     protected:
      inline void SetStateInSyncWithNOXGroup(const bool insync)
      {
        state_is_insync_with_noxgroup_ = insync;
      }

      inline void ThrowIfStateNotInSyncWithNOXGroup() const
      {
        if (!state_is_insync_with_noxgroup_)
          dserror(
              "Mutable state has been requested but the manipulated state has\n"
              "not been communicated to NOX.\n"
              "Manipulations made in the state vector will have no effect.\n"
              "Call SetState(x) to synchronize the states stored in the global\n"
              "state object and in the NOX group!");
      }

     private:
      /// pointer to the different data containers
      Teuchos::RCP<BaseDataIO> dataio_;
      Teuchos::RCP<BaseDataSDyn> datasdyn_;
      Teuchos::RCP<BaseDataGlobalState> dataglobalstate_;

      /// pointer to the integrator (implicit or explicit)
      Teuchos::RCP<STR::Integrator> int_ptr_;

      /// pointer to the dirichlet boundary condition handler
      Teuchos::RCP<STR::Dbc> dbc_ptr_;
    };  // class Base
  }     // namespace TIMINT
}  // namespace STR


#endif /* STR_TIMINT_BASE_H_ */

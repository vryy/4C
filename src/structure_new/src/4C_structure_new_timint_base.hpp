/*-----------------------------------------------------------*/
/*! \file

\brief Base class for all structural time integration strategies.


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_TIMINT_BASE_HPP
#define FOUR_C_STRUCTURE_NEW_TIMINT_BASE_HPP

#include "4C_config.hpp"

#include "4C_adapter_str_structure_new.hpp"
#include "4C_io_every_iteration_writer.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"
#include "4C_structure_new_timint_basedataio.hpp"
#include "4C_structure_new_timint_basedatasdyn.hpp"

// forward declaration
class Epetra_Vector;
class Epetra_Map;

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class BlockSparseMatrixBase;
}  // namespace Core::LinAlg
namespace Adapter
{
  class StructureTimeAda;
}
namespace STR
{
  class ModelEvaluator;
  class Dbc;
  class Integrator;
  namespace MODELEVALUATOR
  {
    class Generic;
  }  // namespace MODELEVALUATOR
  namespace TimeInt
  {
    /** \brief Abstract class for all time integration strategies
     *
     *  \author Michael Hiermeier */
    class Base : public Adapter::StructureNew, Core::IO::EveryIterationWriterInterface
    {
      friend class Adapter::StructureTimeAda;

     public:
      /// constructor
      Base();


      /// initialize (all already existing) class variables
      virtual void Init(const Teuchos::RCP<STR::TimeInt::BaseDataIO> dataio,
          const Teuchos::RCP<STR::TimeInt::BaseDataSDyn> datasdyn,
          const Teuchos::RCP<STR::TimeInt::BaseDataGlobalState> dataglobalstate);

      /// setup of the new class variables
      void Setup() override;

      /// tests if there are more time steps to do
      [[nodiscard]] bool NotFinished() const override;

      /// reset everything (needed for biofilm simulations)
      void Reset() override;

      /** \brief reset step configuration after time step
       *
       *  This function is supposed to reset all variables which are directly related
       *  to the current new step n+1. To be more precise all variables ending with "Np"
       *  have to be reset. */
      void reset_step() override;

      /// wrapper for things that should be done before prepare_time_step is called
      void PrePredict() override {}

      /// wrapper for things that should be done before solving the nonlinear iterations
      void PreSolve() override {}

      /// wrapper for things that should be done after convergence of Newton scheme
      void PostOutput() override {}

      /// things that should be done after the actual time loop is finished
      void PostTimeLoop() override;

      /// @name General access methods
      ///@{
      /// Access discretization (structure only)
      Teuchos::RCP<Discret::Discretization> discretization() override;

      /// Access to pointer to DoF row map of the discretization (structure only)
      const Epetra_Map* DofRowMapView() override
      {
        check_init();
        return dataglobalstate_->DofRowMapView();
      }

      /// DoF map of structural vector of unknowns
      Teuchos::RCP<const Epetra_Map> dof_row_map() override
      {
        check_init();
        return dataglobalstate_->dof_row_map();
      }

      //! DoF map of vector of unknowns
      //! Alternative method capable of multiple DoF sets
      Teuchos::RCP<const Epetra_Map> dof_row_map(unsigned nds) override
      {
        check_init();
        return dataglobalstate_->dof_row_map(nds);
      }

      /// Access linear structural solver
      Teuchos::RCP<Core::LinAlg::Solver> LinearSolver() override
      {
        check_init();
        return datasdyn_->GetLinSolvers()[Inpar::STR::model_structure];
      }

      /// Return MapExtractor for Dirichlet boundary conditions
      Teuchos::RCP<const Core::LinAlg::MapExtractor> GetDBCMapExtractor() override;
      [[nodiscard]] Teuchos::RCP<const Core::LinAlg::MapExtractor> GetDBCMapExtractor() const;

      //! Return locsys manager
      Teuchos::RCP<Core::Conditions::LocsysManager> LocsysManager() override;

      //! Return the desired model evaluator (read-only)
      [[nodiscard]] const STR::MODELEVALUATOR::Generic& ModelEvaluator(
          Inpar::STR::ModelType mtype) const override;

      //! Return the desired model evaluator (read and write)
      STR::MODELEVALUATOR::Generic& ModelEvaluator(Inpar::STR::ModelType mtype) override;

      ///@}

      /// Return domain map of the mass matrix (implicit and explicit)
      [[nodiscard]] const Epetra_Map& GetMassDomainMap() const override;

      /// @name Coupled problem routines
      /// @{
      /// wrapper for things that should be done before updating
      void PreUpdate() override {}

      /// Update routine for coupled problems with monolithic approach
      void Update() override;

      /// Update routine for coupled problems with monolithic approach with time adaptivity
      void Update(double endtime) override = 0;

      /// Update time and step counter
      virtual void UpdateStepTime();

      /// wrapper for things that should be done after solving the update
      void post_update() override;
      /// @}

      /// @name Access global state from outside via adapter (needed for coupled problems)
      ///@{
      /// unknown displacements at \f$t_{n+1}\f$
      [[nodiscard]] Teuchos::RCP<const Epetra_Vector> DispNp() const override
      {
        check_init();
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
       * See also \ref Adapter::StructureNew::set_state
       */
      Teuchos::RCP<Epetra_Vector> WriteAccessDispNp() override
      {
        check_init();
        set_state_in_sync_with_nox_group(false);
        return dataglobalstate_->GetDisNp();
      }

      /// known displacements at \f$t_{n}\f$
      [[nodiscard]] Teuchos::RCP<const Epetra_Vector> DispN() const override
      {
        check_init();
        return dataglobalstate_->GetDisN();
      }

      /// write access to displacements at \f$t^{n}\f$
      Teuchos::RCP<Epetra_Vector> WriteAccessDispN() override
      {
        check_init();
        return dataglobalstate_->GetDisN();
      }

      /// unknown velocities at \f$t_{n+1}\f$
      [[nodiscard]] Teuchos::RCP<const Epetra_Vector> VelNp() const override
      {
        check_init();
        return dataglobalstate_->GetVelNp();
      }

      /// write access to velocities at \f$t^{n+1}\f$
      Teuchos::RCP<Epetra_Vector> WriteAccessVelNp() override
      {
        check_init();
        return dataglobalstate_->GetVelNp();
      }

      /// unknown velocities at \f$t_{n}\f$
      [[nodiscard]] Teuchos::RCP<const Epetra_Vector> VelN() const override
      {
        check_init();
        return dataglobalstate_->GetVelN();
      }

      /// write access to velocities at \f$t^{n}\f$
      Teuchos::RCP<Epetra_Vector> WriteAccessVelN() override
      {
        check_init();
        return dataglobalstate_->GetVelN();
      }

      /// known velocities at \f$t_{n-1}\f$
      [[nodiscard]] Teuchos::RCP<const Epetra_Vector> VelNm() const override
      {
        check_init();
        return dataglobalstate_->GetVelNm();
      }

      /// unknown accelerations at \f$t_{n+1}\f$
      [[nodiscard]] Teuchos::RCP<const Epetra_Vector> AccNp() const override
      {
        check_init();
        return dataglobalstate_->GetAccNp();
      }

      //! known accelerations at \f$t_{n}\f$
      [[nodiscard]] Teuchos::RCP<const Epetra_Vector> AccN() const override
      {
        check_init();
        return dataglobalstate_->GetAccN();
      }
      ///@}

      /// @name access and modify model evaluator stuff via adapter
      /// @{
      /// are there any algebraic constraints?
      bool HaveConstraint() override
      {
        check_init_setup();
        return datasdyn_->HaveModelType(Inpar::STR::model_lag_pen_constraint);
      }

      /// do we need a semi-smooth Newton-type plasticity algorithm
      virtual bool have_semi_smooth_plasticity()
      {
        check_init_setup();
        return datasdyn_->HaveEleTech(Inpar::STR::EleTech::plasticity);
      }

      /// FixMe get constraint manager defined in the structure
      Teuchos::RCP<CONSTRAINTS::ConstrManager> get_constraint_manager() override
      {
        FOUR_C_THROW("Not yet implemented!");
        return Teuchos::null;
      }

      /// FixMe get contact/meshtying manager
      Teuchos::RCP<CONTACT::MeshtyingContactBridge> meshtying_contact_bridge() override
      {
        FOUR_C_THROW("Not yet implemented!");
        return Teuchos::null;
      }

      /// do we have this model
      bool HaveModel(Inpar::STR::ModelType model) override
      {
        return datasdyn_->HaveModelType(model);
      }

      /// Add residual increment to Lagrange multipliers stored in Constraint manager (derived)
      /// FixMe Different behavior for the implicit and explicit case!!!
      void update_iter_incr_constr(Teuchos::RCP<Epetra_Vector> lagrincr) override
      {
        FOUR_C_THROW("Not yet implemented!");
      }

      /// Add residual increment to pressures stored in Cardiovascular0D manager (derived)
      /// FixMe Different behavior for the implicit and explicit case!!!
      void update_iter_incr_cardiovascular0_d(Teuchos::RCP<Epetra_Vector> presincr) override
      {
        FOUR_C_THROW("Not yet implemented!");
      }
      /// @}

      /// @name Time step helpers
      ///@{
      /// Return current time \f$t_{n}\f$ (derived)
      [[nodiscard]] double GetTimeN() const override
      {
        check_init();
        return dataglobalstate_->GetTimeN();
      }

      /// Sets the current time \f$t_{n}\f$ (derived)
      void SetTimeN(const double time_n) override
      {
        check_init();
        dataglobalstate_->GetTimeN() = time_n;
      }

      /// Return target time \f$t_{n+1}\f$ (derived)
      [[nodiscard]] double GetTimeNp() const override
      {
        check_init();
        return dataglobalstate_->GetTimeNp();
      }

      /// Sets the target time \f$t_{n+1}\f$ of this time step (derived)
      void SetTimeNp(const double time_np) override
      {
        check_init();
        dataglobalstate_->GetTimeNp() = time_np;
      }

      /// Get upper limit of time range of interest (derived)
      [[nodiscard]] double GetTimeEnd() const override
      {
        check_init();
        return datasdyn_->GetTimeMax();
      }

      /// Get upper limit of time range of interest (derived)
      void SetTimeEnd(double timemax) override
      {
        check_init();
        datasdyn_->GetTimeMax() = timemax;
      }

      /// Get time step size \f$\Delta t_n\f$
      [[nodiscard]] double GetDeltaTime() const override
      {
        check_init();
        return (*dataglobalstate_->GetDeltaTime())[0];
      }

      /// Set time step size \f$\Delta t_n\f$
      void SetDeltaTime(const double dt) override
      {
        check_init();
        (*dataglobalstate_->GetDeltaTime())[0] = dt;
      }

      /// Return time integration factor
      [[nodiscard]] double TimIntParam() const override;

      /// Return current step number \f$n\f$
      [[nodiscard]] int GetStepN() const override
      {
        check_init();
        return dataglobalstate_->GetStepN();
      }

      /// Sets the current step \f$n\f$
      void SetStepN(int step_n) override
      {
        check_init();
        dataglobalstate_->GetStepN() = step_n;
      }

      /// Return current step number $n+1$
      [[nodiscard]] int GetStepNp() const override
      {
        check_init();
        return dataglobalstate_->GetStepNp();
      }

      /// Sets the current step number \f$n+1\f$
      void SetStepNp(int step_np) override
      {
        check_init_setup();
        dataglobalstate_->GetStepNp() = step_np;
      }

      //! Get number of time steps
      [[nodiscard]] int GetStepEnd() const override
      {
        check_init();
        return datasdyn_->GetStepMax();
      }

      /// Sets number of time steps
      void SetStepEnd(int step_end) override
      {
        check_init_setup();
        datasdyn_->GetStepMax() = step_end;
      }

      //! Get divcont type
      [[nodiscard]] virtual enum Inpar::STR::DivContAct GetDivergenceAction() const
      {
        check_init_setup();
        return datasdyn_->GetDivergenceAction();
      }

      //! Get number of times you want to halve your timestep in case nonlinear solver diverges
      [[nodiscard]] virtual int get_max_div_con_refine_level() const
      {
        check_init_setup();
        return datasdyn_->get_max_div_con_refine_level();
      }

      //! Get random factor for time step adaption
      [[nodiscard]] virtual double get_random_time_step_factor() const
      {
        check_init_setup();
        return datasdyn_->get_random_time_step_factor();
      }

      //! Set random factor for time step adaption
      virtual double set_random_time_step_factor(double rand_tsfac)
      {
        check_init_setup();
        return datasdyn_->get_random_time_step_factor() = rand_tsfac;
      }

      //! Get current refinement level for time step adaption
      [[nodiscard]] virtual int get_div_con_refine_level() const
      {
        check_init_setup();
        return datasdyn_->get_div_con_refine_level();
      }

      //! Set refinement level for time step adaption
      virtual int set_div_con_refine_level(int divconrefinementlevel)
      {
        check_init_setup();
        return datasdyn_->get_div_con_refine_level() = divconrefinementlevel;
      }

      //! Get step of current refinement level for time step adaption
      [[nodiscard]] virtual int get_div_con_num_fine_step() const
      {
        check_init_setup();
        return datasdyn_->get_div_con_num_fine_step();
      }

      //! Set step of current refinement level for time step adaption
      virtual int set_div_con_num_fine_step(int divconnumfinestep)
      {
        check_init_setup();
        return datasdyn_->get_div_con_num_fine_step() = divconnumfinestep;
      }

      /// set evaluation action
      void SetActionType(const Core::Elements::ActionType& action) override;

      // group id in nested parallelity
      [[nodiscard]] int GroupId() const;
      ///@}

      /// @name Structure with ale specific methods
      ///@{
      /// FixMe set/apply material displacements to structure field (structure with ale)
      void SetDispMatNp(Teuchos::RCP<Epetra_Vector> dispmatnp) override
      {
        FOUR_C_THROW("Not supported at the moment!");
      }

      /// FixMe write access to material displacements (strutcure with ale) at \f$t^{n+1}\f$
      Teuchos::RCP<Epetra_Vector> write_access_disp_mat_np() override
      {
        check_init_setup();
        FOUR_C_THROW("Not yet supported!");
        return Teuchos::null;
      }
      ///@}


      /// Time adaptivity (derived pure virtual functionality)
      /// @{
      /// Resize MStep Object due to time adaptivity in FSI (derived)
      void ResizeMStepTimAda() override;

      /// @}

      /// Output writer related routines (file and screen output)
      /// @{
      /// Access output object
      Teuchos::RCP<Core::IO::DiscretizationWriter> DiscWriter() override
      {
        return data_io().GetOutputPtr();
      }

      /// Calculate all output quantities depending on the constitutive model
      /// (and, hence, on a potential material history)
      void prepare_output(bool force_prepare_timestep) override;

      /// output results (implicit and explicit)
      virtual void Output() { Output(false); }
      void Output(bool forced_writerestart) override;

      /// Write Gmsh output for structural field
      void write_gmsh_struc_output_step() override;

      /// FixMe Check if there are any elements with the micro material definition.
      /// Maybe the detection can be moved to the element loop in the ad_str_structure_new.cpp.
      /// There is already one.
      bool HaveMicroMat() override
      {
        FOUR_C_THROW("Not yet considered!");
        return false;
      }

      /// create result test for encapsulated structure algorithm
      Teuchos::RCP<Core::UTILS::ResultTest> CreateFieldTest() override;

      /** \brief Get data that is written during restart
       *
       *  This routine is only for simple structure problems!
       *  \date 06/13
       *  \author biehler */
      void GetRestartData(Teuchos::RCP<int> step, Teuchos::RCP<double> time,
          Teuchos::RCP<Epetra_Vector> disnp, Teuchos::RCP<Epetra_Vector> velnp,
          Teuchos::RCP<Epetra_Vector> accnp, Teuchos::RCP<std::vector<char>> elementdata,
          Teuchos::RCP<std::vector<char>> nodedata) override;

      /** Read restart values
       *
       * \param stepn (in): restart step at \f${n}\f$
       */
      void read_restart(const int stepn) override;

      /// Set restart values (deprecated)
      void SetRestart(int stepn,                        //!< restart step at \f${n}\f$
          double timen,                                 //!< restart time at \f$t_{n}\f$
          Teuchos::RCP<Epetra_Vector> disn,             //!< restart displacements at \f$t_{n}\f$
          Teuchos::RCP<Epetra_Vector> veln,             //!< restart velocities at \f$t_{n}\f$
          Teuchos::RCP<Epetra_Vector> accn,             //!< restart accelerations at \f$t_{n}\f$
          Teuchos::RCP<std::vector<char>> elementdata,  //!< restart element data
          Teuchos::RCP<std::vector<char>> nodedata      //!< restart element data
          ) override;
      /// @}

      /// Biofilm related stuff
      /// @{
      /// FixMe set structure displacement vector due to biofilm growth
      void SetStrGrDisp(Teuchos::RCP<Epetra_Vector> struct_growth_disp) override
      {
        FOUR_C_THROW("Currently unsupported!");
      }
      /// @}

      /// @name Pure virtual adapter functions (have to be implemented in the derived classes)
      /// @{
      /// integrate the current step (implicit and explicit)
      virtual int IntegrateStep() = 0;
      /// right-hand-side of Newton's method (implicit only)
      Teuchos::RCP<const Epetra_Vector> RHS() override { return GetF(); };
      [[nodiscard]] virtual Teuchos::RCP<const Epetra_Vector> GetF() const = 0;
      /// @}

     public:
      /// @name External accessors for the class variables
      ///@{
      /// Get the indicator if we are currently restarting the simulation
      [[nodiscard]] inline const bool& IsRestarting() const { return isrestarting_; }

      /// Get the indicator if we need to restart the initial state
      [[nodiscard]] inline bool is_restarting_initial_state() const
      {
        return datasdyn_->is_restarting_initial_state();
      }

      /// Get TimIntBase data for global state quantities (read access)
      [[nodiscard]] Teuchos::RCP<const BaseDataGlobalState> get_data_global_state_ptr() const
      {
        check_init();
        return dataglobalstate_;
      }

      /// Get TimIntBase data for global state quantities (read & write access)
      Teuchos::RCP<BaseDataGlobalState>& get_data_global_state_ptr()
      {
        check_init();
        return dataglobalstate_;
      }

      [[nodiscard]] const BaseDataGlobalState& GetDataGlobalState() const
      {
        check_init();
        return *dataglobalstate_;
      }

      /// Get TimIntBase data for io quantities (read access)
      [[nodiscard]] Teuchos::RCP<const BaseDataIO> GetDataIOPtr() const
      {
        check_init();
        return dataio_;
      }

      [[nodiscard]] const BaseDataIO& GetDataIO() const
      {
        check_init();
        return *dataio_;
      }

      /// Get TimIntBase data or struct dynamics quantitites (read access)
      [[nodiscard]] Teuchos::RCP<const BaseDataSDyn> GetDataSDynPtr() const
      {
        check_init();
        return datasdyn_;
      }

      [[nodiscard]] const BaseDataSDyn& GetDataSDyn() const
      {
        check_init();
        return *datasdyn_;
      }

      /// return a reference to the Dirichlet Boundary Condition handler (read access)
      [[nodiscard]] const STR::Dbc& GetDBC() const
      {
        check_init_setup();
        return *dbc_ptr_;
      }

      /// return a reference to the Dirichlet Boundary Condition handler (write access)
      STR::Dbc& GetDBC()
      {
        check_init_setup();
        return *dbc_ptr_;
      }

      /// return a pointer to the Dirichlet Boundary Condition handler (read access)
      [[nodiscard]] Teuchos::RCP<const STR::Dbc> GetDBCPtr() const
      {
        check_init_setup();
        return dbc_ptr_;
      }

      /// return the integrator (read-only)
      [[nodiscard]] const STR::Integrator& Integrator() const
      {
        check_init_setup();
        return *int_ptr_;
      }

      /// Get the global state
      const BaseDataGlobalState& data_global_state() const
      {
        check_init();
        return *dataglobalstate_;
      }

      /// Get internal TimIntBase data for structural dynamics quantities (read and write access)
      BaseDataSDyn& DataSDyn()
      {
        check_init();
        return *datasdyn_;
      }

      /// return a pointer to the Dirichlet Boundary Condition handler (read and write access)
      const Teuchos::RCP<STR::Dbc>& DBCPtr()
      {
        check_init_setup();
        return dbc_ptr_;
      }

      [[nodiscard]] bool has_final_state_been_written() const override;

      /// get the indicator state
      [[nodiscard]] inline const bool& is_init() const { return isinit_; }

      /// get the indicator state
      [[nodiscard]] inline const bool& is_setup() const { return issetup_; }

      //! @name Attribute access functions
      //@{

      //! Provide Name
      virtual enum Inpar::STR::DynamicType MethodName() const = 0;

      //! Provide title
      std::string MethodTitle() const { return Inpar::STR::DynamicTypeString(MethodName()); }

      //! Return true, if time integrator is implicit
      virtual bool IsImplicit() const = 0;

      //! Return true, if time integrator is explicit
      virtual bool IsExplicit() const = 0;

      //! Provide number of steps, e.g. a single-step method returns 1,
      //! a \f$m\f$-multistep method returns \f$m\f$
      virtual int MethodSteps() const = 0;

      //! Give order of accuracy
      int method_order_of_accuracy() const
      {
        return std::min(method_order_of_accuracy_dis(), method_order_of_accuracy_vel());
      }

      //! Give local order of accuracy of displacement part
      virtual int method_order_of_accuracy_dis() const = 0;

      //! Give local order of accuracy of velocity part
      virtual int method_order_of_accuracy_vel() const = 0;

      //! Return linear error coefficient of displacements
      virtual double method_lin_err_coeff_dis() const = 0;

      //! Return linear error coefficient of velocities
      virtual double method_lin_err_coeff_vel() const = 0;

      //@}

      ///@}
     protected:
      /// Check if Init() and Setup() have been called, yet.
      inline void check_init_setup() const
      {
        FOUR_C_ASSERT(is_init() and is_setup(), "Call Init() and Setup() first!");
      }

      /// Check if Init() has been called
      inline void check_init() const { FOUR_C_ASSERT(is_init(), "Call Init() first!"); }

      /// Get the global state
      BaseDataGlobalState& data_global_state()
      {
        check_init();
        return *dataglobalstate_;
      }

      /// Get the pointer to global state
      const Teuchos::RCP<BaseDataGlobalState>& data_global_state_ptr() const
      {
        check_init();
        return dataglobalstate_;
      }

      /// Get internal TimIntBase data for io quantities (read and write access)
      BaseDataIO& data_io()
      {
        check_init();
        return *dataio_;
      }

      /// return a pointer to the input/output data container (read and write access)
      const Teuchos::RCP<BaseDataIO>& data_io_ptr()
      {
        check_init();
        return dataio_;
      }

      /// return a pointer to the structural dynamic data container (read and write access)
      const Teuchos::RCP<BaseDataSDyn>& data_s_dyn_ptr()
      {
        check_init();
        return datasdyn_;
      }

      /// return a reference to the Dirichlet Boundary Condition handler (read and write access)
      STR::Dbc& dbc()
      {
        check_init_setup();
        return *dbc_ptr_;
      }

      /// return a reference to the integrator (read and write access)
      STR::Integrator& integrator()
      {
        check_init_setup();
        return *int_ptr_;
      }

      /// return a pointer to the integrator (read and write access)
      const Teuchos::RCP<STR::Integrator>& integrator_ptr()
      {
        check_init_setup();
        return int_ptr_;
      }

      /** \brief Output to file
       *
       *  This routine always prints the last converged state, i.e.
       *  \f$D_{n}, V_{n}, A_{n}\f$.
       *
       *  \date 03/07
       *  \author mwgee (originally) */
      void output_step(bool forced_writerestart);

     private:
      /*! \brief Create a new input/output step in the output writer
       *
       * New step is created only once per time step. This is controlled by \c datawritten.
       * Do nothing if data has already been written in this time step.
       *
       * \param[in,out] Indicator whether data has already been written in this time step (true) or
       *                not (false)
       */
      void new_io_step(bool& datawritten);

      /// output of the current state
      void output_state();

      /** \brief output of the current state */
      void output_state(Core::IO::DiscretizationWriter& iowriter, bool write_owner) const;

      /** \brief output of the debug state */
      void OutputDebugState(
          Core::IO::DiscretizationWriter& iowriter, bool write_owner) const override;

      /// output during runtime
      void runtime_output_state();

      /// output reaction forces
      void output_reaction_forces();

      /// output element volumes
      void output_element_volume(Core::IO::DiscretizationWriter& iowriter) const;

      /// output stress and/or strain state
      void output_stress_strain();

      /// output energy
      void output_energy() const;

      /// output optional quantity
      void output_optional_quantity();

      /// write restart information
      void output_restart(bool& datawritten);

      /// add restart information to output state
      void add_restart_to_output_state();

      /** \brief set the number of nonlinear iterations of the last time step
       *
       *  \pre UpdateStepTime() must be called beforehand, otherwise the wrong
       *  step-id is considered.
       *
       *  \author hiermeier \date 11/17 */
      void set_number_of_nonlinear_iterations();

      /** \brief decide which contributions to the total system energy shall be
       *         computed and written during simulation
       *
       *  \author grill */
      void select_energy_types_to_be_written();

      /** \brief initialize file stream for energy values and write all the
       *         column headers for the previously selected energy contributions
       *         to be written separately
       *
       *  \author grill */
      void initialize_energy_file_stream_and_write_headers();

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
      inline void set_state_in_sync_with_nox_group(const bool insync)
      {
        state_is_insync_with_noxgroup_ = insync;
      }

      inline void throw_if_state_not_in_sync_with_nox_group() const
      {
        if (!state_is_insync_with_noxgroup_)
        {
          FOUR_C_THROW(
              " state has been requested but the manipulated state has\n"
              "not been communicated to NOX.\n"
              "Manipulations made in the state vector will have no effect.\n"
              "Call set_state(x) to synchronize the states stored in the global\n"
              "state object and in the NOX group!");
        }
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
  }     // namespace TimeInt
}  // namespace STR


FOUR_C_NAMESPACE_CLOSE

#endif

/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the structural time integration

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_STR_WRAPPER_HPP
#define FOUR_C_ADAPTER_STR_WRAPPER_HPP

#include "4C_config.hpp"

#include "4C_adapter_str_structure.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  /// This class is a wrapper of Adapter::Structure. It follows the "decorator" design pattern. This
  /// approach allows to dynamically add functionalities to an instance of Adapter::Structure. For
  /// example, Adapter::StructureTimeLoop implements the Integrate function for sequential time
  /// marching.
  class StructureWrapper : public Structure
  {
   public:
    /// constructor
    explicit StructureWrapper(Teuchos::RCP<Structure> structure) : structure_(std::move(structure))
    {
    }

    //! @name Construction
    //@{

    /*! \brief Setup all class internal objects and members

     Setup() is not supposed to have any input arguments !

     Must only be called after Init().

     Construct all objects depending on the parallel distribution and
     relying on valid maps like, e.g. the state vectors, system matrices, etc.

     Call all Setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning none
    \return void
    \date 08/16
    \author rauch  */
    void Setup() override { structure_->Setup(); };

    //@}

    //! @name Vector access
    //@{

    /// initial guess of Newton's method
    Teuchos::RCP<const Epetra_Vector> initial_guess() override
    {
      return structure_->initial_guess();
    }

    /// right-hand-side of Newton's method
    Teuchos::RCP<const Epetra_Vector> RHS() override { return structure_->RHS(); }

    /// unknown displacements at \f$t_{n+1}\f$
    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> Dispnp() const override
    {
      return structure_->Dispnp();
    }

    /// known displacements at \f$t_{n}\f$
    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> Dispn() const override
    {
      return structure_->Dispn();
    }

    /// unknown velocity at \f$t_{n+1}\f$
    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> Velnp() const override
    {
      return structure_->Velnp();
    }

    /// known velocity at \f$t_{n}\f$
    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> Veln() const override
    {
      return structure_->Veln();
    }

    /// known velocity at \f$t_{n-1}\f$
    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> Velnm() const override
    {
      return structure_->Velnm();
    }

    /// unknown acceleration at \f$t_{n+1}\f$
    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> Accnp() const override
    {
      return structure_->Accnp();
    }

    /// known acceleration at \f$t_{n}\f$
    [[nodiscard]] Teuchos::RCP<const Epetra_Vector> Accn() const override
    {
      return structure_->Accn();
    }

    //@}


    //! @name Misc
    //@{

    /// dof map of vector of unknowns
    Teuchos::RCP<const Epetra_Map> dof_row_map() override { return structure_->dof_row_map(); }

    /// dof map of vector of unknowns for multiple dof sets
    Teuchos::RCP<const Epetra_Map> dof_row_map(unsigned nds) override
    {
      return structure_->dof_row_map(nds);
    }

    /// view of dof map of vector of vector of unknowns
    const Epetra_Map* dof_row_map_view() override { return structure_->dof_row_map_view(); }

    /// domain map of system matrix
    [[nodiscard]] const Epetra_Map& DomainMap() const override { return structure_->DomainMap(); }

    /// direct access to system matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> system_matrix() override
    {
      return structure_->system_matrix();
    }

    /// direct access to system matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> block_system_matrix() override
    {
      return structure_->block_system_matrix();
    }

    /// switch structure field to block matrix
    void use_block_matrix(Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> domainmaps,
        Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> rangemaps) override
    {
      structure_->use_block_matrix(domainmaps, rangemaps);
    }

    // access to contact/meshtying bridge
    Teuchos::RCP<CONTACT::MeshtyingContactBridge> meshtying_contact_bridge() override
    {
      return structure_->meshtying_contact_bridge();
    }

    // access to locsys manager
    Teuchos::RCP<Core::Conditions::LocsysManager> LocsysManager() override
    {
      return structure_->LocsysManager();
    }

    /// direct access to discretization
    Teuchos::RCP<Core::FE::Discretization> discretization() override
    {
      return structure_->discretization();
    }

    /// read only access to discretization
    [[nodiscard]] virtual Teuchos::RCP<const Core::FE::Discretization> get_discretization() const
    {
      return structure_->discretization();
    }

    /// are there any algebraic constraints?
    bool HaveConstraint() override { return structure_->HaveConstraint(); }

    /// are there any spring dashpot BCs?
    bool HaveSpringDashpot() override { return structure_->HaveSpringDashpot(); }

    /// get constraint manager defined in the structure
    Teuchos::RCP<CONSTRAINTS::ConstrManager> get_constraint_manager() override
    {
      return structure_->get_constraint_manager();
    }

    /// get constraint manager defined in the structure
    Teuchos::RCP<CONSTRAINTS::SpringDashpotManager> get_spring_dashpot_manager() override
    {
      return structure_->get_spring_dashpot_manager();
    }

    /// get type of thickness scaling for thin shell structures
    Inpar::STR::StcScale get_stc_algo() override { return structure_->get_stc_algo(); }

    /// access to scaling matrix for STC
    Teuchos::RCP<Core::LinAlg::SparseMatrix> get_stc_mat() override
    {
      return structure_->get_stc_mat();
    }

    /// Return MapExtractor for Dirichlet boundary conditions
    Teuchos::RCP<const Core::LinAlg::MapExtractor> GetDBCMapExtractor() override
    {
      return structure_->GetDBCMapExtractor();
    }

    /// expand dirichlet bc map
    void AddDirichDofs(const Teuchos::RCP<const Epetra_Map> maptoadd) override
    {
      structure_->AddDirichDofs(maptoadd);
    };

    /// contract dirichlet bc map
    void RemoveDirichDofs(const Teuchos::RCP<const Epetra_Map> maptoremove) override
    {
      structure_->RemoveDirichDofs(maptoremove);
    };

    /// reset step and state vectors
    void Reset() override { structure_->Reset(); }

    /// reset last time step, needed for time step size adaptivity in FSI
    void reset_step() override { structure_->reset_step(); }

    //@}


    /// @name Time step helpers
    //@{

    /// return time integration factor
    [[nodiscard]] double TimIntParam() const override { return structure_->TimIntParam(); }

    //! Sets the current time \f$t_{n}\f$
    void set_time(const double time) override { structure_->set_time(time); }

    //! Sets the target time \f$t_{n+1}\f$ of this time step
    void SetTimen(const double time) override { structure_->SetTimen(time); }

    //! Sets the target step \f$n\f$
    void SetStep(int step) override { structure_->SetStep(step); }

    //! Sets the target step \f$n+1\f$
    void SetStepn(int step) override { structure_->SetStepn(step); }

    //! Return current time \f$t_{n}\f$
    [[nodiscard]] double TimeOld() const override { return structure_->TimeOld(); }

    //! Return target time \f$t_{n+1}\f$
    [[nodiscard]] double Time() const override { return structure_->Time(); }

    /// get upper limit of time range of interest
    [[nodiscard]] double GetTimeEnd() const override { return structure_->GetTimeEnd(); }

    //! Set upper limit of time range of interest //HACK for parameter continuation
    void SetTimeEnd(double timemax) override { structure_->SetTimeEnd(timemax); }

    /// get time step size \f$\Delta t_n\f$
    [[nodiscard]] double Dt() const override { return structure_->Dt(); }

    /// Return current step number $n$
    [[nodiscard]] int StepOld() const override { return structure_->StepOld(); }

    /// Return current step number $n+1$
    [[nodiscard]] int Step() const override { return structure_->Step(); }

    /// get number of time steps
    [[nodiscard]] int NumStep() const override { return structure_->NumStep(); }

    /// integrate from t1 to t2
    int Integrate() override { return structure_->Integrate(); }

    //! do something in case nonlinear solution does not converge for some reason
    Inpar::STR::ConvergenceStatus PerformErrorAction(
        Inpar::STR::ConvergenceStatus nonlinsoldiv) override
    {
      return structure_->PerformErrorAction(nonlinsoldiv);
    }

    /// tests if there are more time steps to do
    [[nodiscard]] bool not_finished() const override { return structure_->not_finished(); }

    /// set time step size
    void set_dt(const double dtnew) override { structure_->set_dt(dtnew); }

    /// start new time step
    void prepare_time_step() override { structure_->prepare_time_step(); }

    /// update displacment
    void update_state_incrementally(
        Teuchos::RCP<const Epetra_Vector> disi  ///< iterative solution increment
        ) override
    {
      structure_->update_state_incrementally(disi);
    }

    void determine_stress_strain() override { structure_->determine_stress_strain(); }

    /// update displacement and evaluate elements (implicit only)
    void evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc) override
    {
      structure_->evaluate(disiterinc);
    }

    /// don't update displacement but evaluate elements (implicit only)
    void evaluate() override { structure_->evaluate(); }

    /// update at time step end
    void Update() override { structure_->Update(); }

    /// update at time step end
    void Update(const double endtime) override { structure_->Update(endtime); }

    /// resize MStep objects for AB2
    void resize_m_step_tim_ada() override { structure_->resize_m_step_tim_ada(); }

    /// update iteration; add residual increment to Lagrange multipliers stored in Constraint
    /// manager
    void update_iter_incr_constr(Teuchos::RCP<Epetra_Vector> lagrincr) override
    {
      structure_->update_iter_incr_constr(lagrincr);
    }

    /// update iteration; add residual increment to pressures stored in 0D cardiovascular manager
    void update_iter_incr_cardiovascular0_d(Teuchos::RCP<Epetra_Vector> presincr) override
    {
      structure_->update_iter_incr_cardiovascular0_d(presincr);
    }

    /// access to output object
    Teuchos::RCP<Core::IO::DiscretizationWriter> disc_writer() override
    {
      return structure_->disc_writer();
    }

    /// prepare output (i.e. calculate stresses, strains, energies)
    void prepare_output(bool force_prepare) override { structure_->prepare_output(force_prepare); }

    /// Get restart data
    void get_restart_data(Teuchos::RCP<int> step, Teuchos::RCP<double> time,
        Teuchos::RCP<Epetra_Vector> disn, Teuchos::RCP<Epetra_Vector> veln,
        Teuchos::RCP<Epetra_Vector> accn, Teuchos::RCP<std::vector<char>> elementdata,
        Teuchos::RCP<std::vector<char>> nodedata) override
    {
      structure_->get_restart_data(step, time, disn, veln, accn, elementdata, nodedata);
    }

    /// output results
    void Output(bool forced_writerestart = false) override
    {
      structure_->Output(forced_writerestart);
    }

    /// Write Gmsh output for structural field
    void write_gmsh_struc_output_step() override { structure_->write_gmsh_struc_output_step(); }

    /// output results to screen
    void print_step() override { structure_->print_step(); }

    /// read restart information for given time step
    void read_restart(const int step) override { structure_->read_restart(step); }

    /// set restart information for parameter continuation
    void set_restart(int step, double time, Teuchos::RCP<Epetra_Vector> disn,
        Teuchos::RCP<Epetra_Vector> veln, Teuchos::RCP<Epetra_Vector> accn,
        Teuchos::RCP<std::vector<char>> elementdata,
        Teuchos::RCP<std::vector<char>> nodedata) override
    {
      structure_->set_restart(step, time, disn, veln, accn, elementdata, nodedata);
    }

    /// set the state of the nox group and the global state data container (implicit only)
    void set_state(const Teuchos::RCP<Epetra_Vector>& x) override { structure_->set_state(x); }

    /// set evaluation action
    void SetActionType(const Core::Elements::ActionType& action) override
    {
      structure_->SetActionType(action);
    }

    /// wrapper for things that should be done before prepare_time_step is called
    void PrePredict() override { structure_->PrePredict(); }

    /// wrapper for things that should be done before solving the nonlinear iterations
    void PreSolve() override { structure_->PreSolve(); }

    /// wrapper for things that should be done before updating
    void PreUpdate() override { structure_->PreUpdate(); }

    /// wrapper for things that should be done after solving the update
    void post_update() override { structure_->post_update(); }

    /// wrapper for things that should be done after the output
    void PostOutput() override { structure_->PostOutput(); }

    /// wrapper for things that should be done after the actual time loop is finished
    void PostTimeLoop() override { structure_->PostTimeLoop(); }

    //@}


    //! @name Solver calls
    //@{

    /// nonlinear solve
    Inpar::STR::ConvergenceStatus Solve() override { return structure_->Solve(); }

    //! linear structure solve with just an interface load
    Teuchos::RCP<Epetra_Vector> solve_relaxation_linear() override
    {
      return structure_->solve_relaxation_linear();
    }

    /// get the linear solver object used for this field
    Teuchos::RCP<Core::LinAlg::Solver> linear_solver() override
    {
      return structure_->linear_solver();
    }

    //@}


    //! @name volume coupled specific methods
    //@{

    /// set forces due to interface with fluid, the force is expected external-force-like
    void SetForceInterface(Teuchos::RCP<Epetra_MultiVector> iforce) override
    {
      structure_->SetForceInterface(iforce);
    }

    //! specific method for iterative staggered partitioned TSI
    //! will be obsolete after switch to new structural timint.
    void prepare_partition_step() override { structure_->prepare_partition_step(); }

    //@}


    //! @name Write access to field solution variables at \f$t^{n+1}\f$
    //@{

    /// write access to extract displacements at \f$t^{n+1}\f$
    Teuchos::RCP<Epetra_Vector> WriteAccessDispnp() override
    {
      return structure_->WriteAccessDispnp();
    }

    /// write access to extract velocities at \f$t^{n+1}\f$
    Teuchos::RCP<Epetra_Vector> WriteAccessVelnp() override
    {
      return structure_->WriteAccessVelnp();
    }

    /// write access to extract displacements at \f$t^{n}\f$
    Teuchos::RCP<Epetra_Vector> WriteAccessDispn() override
    {
      return structure_->WriteAccessDispn();
    }

    /// write access to extract velocities at \f$t^{n}\f$
    Teuchos::RCP<Epetra_Vector> WriteAccessVeln() override
    {
      return structure_->WriteAccessVelnp();
    }

    //@}

    /// extract rhs (used to calculate reaction force for post-processing)
    Teuchos::RCP<Epetra_Vector> freact() override { return structure_->freact(); }

    //! @name Structure with ale specific methods
    //@{

    /// material displacements (structure with ale)
    Teuchos::RCP<Epetra_Vector> DispMat() override { return structure_->DispMat(); }

    /// apply material displacements to structure field (structure with ale)
    void ApplyDisMat(Teuchos::RCP<Epetra_Vector> dismat) override
    {
      structure_->ApplyDisMat(dismat);
    }

    //@}


    /// create result test for encapsulated structure algorithm
    Teuchos::RCP<Core::UTILS::ResultTest> CreateFieldTest() override
    {
      return structure_->CreateFieldTest();
    }

    //! @name Biofilm specific methods
    //@{

    void SetStrGrDisp(Teuchos::RCP<Epetra_Vector> struct_growth_disp) override
    {
      structure_->SetStrGrDisp(struct_growth_disp);
    }

    //@}

    /// bool indicating if micro material is used
    bool HaveMicroMat() override { return structure_->HaveMicroMat(); }

    /// do we have this model
    bool HaveModel(Inpar::STR::ModelType model) override { return structure_->HaveModel(model); }

    /// return model evaluator
    STR::MODELEVALUATOR::Generic& ModelEvaluator(Inpar::STR::ModelType mtype) override
    {
      return structure_->ModelEvaluator(mtype);
    }

    [[nodiscard]] bool has_final_state_been_written() const override
    {
      return structure_->has_final_state_been_written();
    }

   protected:
    Teuchos::RCP<Structure> structure_;  ///< underlying structural time integration
  };


  /// Calculate increments from absolute values
  class StructureNOXCorrectionWrapper : public StructureWrapper
  {
   public:
    explicit StructureNOXCorrectionWrapper(Teuchos::RCP<Structure> structure)
        : StructureWrapper(structure)
    {
    }

    void prepare_time_step() override;

    //! evaluate() routine that can handle NOX step increments by computing the
    //! last iteration increment needed for structural evaluate() call
    void evaluate(Teuchos::RCP<const Epetra_Vector> disstepinc) override;

   private:
    /// sum of displacement increments already applied,
    ///
    /// there are two increments around
    ///
    /// x^n+1_i+1 = x^n+1_i + disiterinc  (also referred to as residual increment)
    ///
    /// x^n+1_i+1 = x^n     + disstepinc
    Teuchos::RCP<Epetra_Vector> disstepinc_;
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif

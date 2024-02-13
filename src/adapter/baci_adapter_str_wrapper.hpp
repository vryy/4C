/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the structural time integration

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_ADAPTER_STR_WRAPPER_HPP
#define BACI_ADAPTER_STR_WRAPPER_HPP

#include "baci_config.hpp"

#include "baci_adapter_str_structure.hpp"

BACI_NAMESPACE_OPEN

namespace ADAPTER
{
  /// This class is a wrapper of ADAPTER::Structure. It follows the "decorator" design pattern. This
  /// approach allows to dynamically add functionalities to an instance of ADAPTER::Structure. For
  /// example, ADAPTER::StructureTimeLoop implements the Integrate function for sequential time
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
    Teuchos::RCP<const Epetra_Vector> InitialGuess() override { return structure_->InitialGuess(); }

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
    Teuchos::RCP<const Epetra_Map> DofRowMap() override { return structure_->DofRowMap(); }

    /// dof map of vector of unknowns for multiple dof sets
    Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds) override
    {
      return structure_->DofRowMap(nds);
    }

    /// view of dof map of vector of vector of unknowns
    const Epetra_Map* DofRowMapView() override { return structure_->DofRowMapView(); }

    /// domain map of system matrix
    [[nodiscard]] const Epetra_Map& DomainMap() const override { return structure_->DomainMap(); }

    /// direct access to system matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> SystemMatrix() override
    {
      return structure_->SystemMatrix();
    }

    /// direct access to system matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> BlockSystemMatrix() override
    {
      return structure_->BlockSystemMatrix();
    }

    /// switch structure field to block matrix
    void UseBlockMatrix(Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> domainmaps,
        Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> rangemaps) override
    {
      structure_->UseBlockMatrix(domainmaps, rangemaps);
    }

    // access to contact/meshtying bridge
    Teuchos::RCP<CONTACT::MeshtyingContactBridge> MeshtyingContactBridge() override
    {
      return structure_->MeshtyingContactBridge();
    }

    // access to locsys manager
    Teuchos::RCP<DRT::UTILS::LocsysManager> LocsysManager() override
    {
      return structure_->LocsysManager();
    }

    /// direct access to discretization
    Teuchos::RCP<DRT::Discretization> Discretization() override
    {
      return structure_->Discretization();
    }

    /// read only access to discretization
    [[nodiscard]] virtual Teuchos::RCP<const DRT::Discretization> GetDiscretization() const
    {
      return structure_->Discretization();
    }

    /// are there any algebraic constraints?
    bool HaveConstraint() override { return structure_->HaveConstraint(); }

    /// are there any spring dashpot BCs?
    bool HaveSpringDashpot() override { return structure_->HaveSpringDashpot(); }

    /// get constraint manager defined in the structure
    Teuchos::RCP<CONSTRAINTS::ConstrManager> GetConstraintManager() override
    {
      return structure_->GetConstraintManager();
    }

    /// get constraint manager defined in the structure
    Teuchos::RCP<CONSTRAINTS::SpringDashpotManager> GetSpringDashpotManager() override
    {
      return structure_->GetSpringDashpotManager();
    }

    /// get type of thickness scaling for thin shell structures
    INPAR::STR::STC_Scale GetSTCAlgo() override { return structure_->GetSTCAlgo(); }

    /// access to scaling matrix for STC
    Teuchos::RCP<CORE::LINALG::SparseMatrix> GetSTCMat() override
    {
      return structure_->GetSTCMat();
    }

    /// Return MapExtractor for Dirichlet boundary conditions
    Teuchos::RCP<const CORE::LINALG::MapExtractor> GetDBCMapExtractor() override
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
    void ResetStep() override { structure_->ResetStep(); }

    //@}


    /// @name Time step helpers
    //@{

    /// return time integration factor
    [[nodiscard]] double TimIntParam() const override { return structure_->TimIntParam(); }

    //! Sets the current time \f$t_{n}\f$
    void SetTime(const double time) override { structure_->SetTime(time); }

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
    INPAR::STR::ConvergenceStatus PerformErrorAction(
        INPAR::STR::ConvergenceStatus nonlinsoldiv) override
    {
      return structure_->PerformErrorAction(nonlinsoldiv);
    }

    /// tests if there are more time steps to do
    [[nodiscard]] bool NotFinished() const override { return structure_->NotFinished(); }

    /// set time step size
    void SetDt(const double dtnew) override { structure_->SetDt(dtnew); }

    /// start new time step
    void PrepareTimeStep() override { structure_->PrepareTimeStep(); }

    /// update displacment
    void UpdateStateIncrementally(
        Teuchos::RCP<const Epetra_Vector> disi  ///< iterative solution increment
        ) override
    {
      structure_->UpdateStateIncrementally(disi);
    }

    void DetermineStressStrain() override { structure_->DetermineStressStrain(); }

    /// update displacement and evaluate elements (implicit only)
    void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc) override
    {
      structure_->Evaluate(disiterinc);
    }

    /// don't update displacement but evaluate elements (implicit only)
    void Evaluate() override { structure_->Evaluate(); }

    /// update at time step end
    void Update() override { structure_->Update(); }

    /// update at time step end
    void Update(const double endtime) override { structure_->Update(endtime); }

    /// resize MStep objects for AB2
    void ResizeMStepTimAda() override { structure_->ResizeMStepTimAda(); }

    /// update iteration; add residual increment to Lagrange multipliers stored in Constraint
    /// manager
    void UpdateIterIncrConstr(Teuchos::RCP<Epetra_Vector> lagrincr) override
    {
      structure_->UpdateIterIncrConstr(lagrincr);
    }

    /// update iteration; add residual increment to pressures stored in 0D cardiovascular manager
    void UpdateIterIncrCardiovascular0D(Teuchos::RCP<Epetra_Vector> presincr) override
    {
      structure_->UpdateIterIncrCardiovascular0D(presincr);
    }

    /// access to output object
    Teuchos::RCP<IO::DiscretizationWriter> DiscWriter() override
    {
      return structure_->DiscWriter();
    }

    /// prepare output (i.e. calculate stresses, strains, energies)
    void PrepareOutput(bool force_prepare) override { structure_->PrepareOutput(force_prepare); }

    /// Get restart data
    void GetRestartData(Teuchos::RCP<int> step, Teuchos::RCP<double> time,
        Teuchos::RCP<Epetra_Vector> disn, Teuchos::RCP<Epetra_Vector> veln,
        Teuchos::RCP<Epetra_Vector> accn, Teuchos::RCP<std::vector<char>> elementdata,
        Teuchos::RCP<std::vector<char>> nodedata) override
    {
      structure_->GetRestartData(step, time, disn, veln, accn, elementdata, nodedata);
    }

    /// output results
    void Output(bool forced_writerestart = false) override
    {
      structure_->Output(forced_writerestart);
    }

    /// Write Gmsh output for structural field
    void WriteGmshStrucOutputStep() override { structure_->WriteGmshStrucOutputStep(); }

    /// output results to screen
    void PrintStep() override { structure_->PrintStep(); }

    /// read restart information for given time step
    void ReadRestart(const int step) override { structure_->ReadRestart(step); }

    /// set restart information for parameter continuation
    void SetRestart(int step, double time, Teuchos::RCP<Epetra_Vector> disn,
        Teuchos::RCP<Epetra_Vector> veln, Teuchos::RCP<Epetra_Vector> accn,
        Teuchos::RCP<std::vector<char>> elementdata,
        Teuchos::RCP<std::vector<char>> nodedata) override
    {
      structure_->SetRestart(step, time, disn, veln, accn, elementdata, nodedata);
    }

    /// set the state of the nox group and the global state data container (implicit only)
    void SetState(const Teuchos::RCP<Epetra_Vector>& x) override { structure_->SetState(x); }

    /// set evaluation action
    void SetActionType(const DRT::ELEMENTS::ActionType& action) override
    {
      structure_->SetActionType(action);
    }

    /// wrapper for things that should be done before PrepareTimeStep is called
    void PrePredict() override { structure_->PrePredict(); }

    /// wrapper for things that should be done before solving the nonlinear iterations
    void PreSolve() override { structure_->PreSolve(); }

    /// wrapper for things that should be done before updating
    void PreUpdate() override { structure_->PreUpdate(); }

    /// wrapper for things that should be done after solving the update
    void PostUpdate() override { structure_->PostUpdate(); }

    /// wrapper for things that should be done after the output
    void PostOutput() override { structure_->PostOutput(); }

    /// wrapper for things that should be done after the actual time loop is finished
    void PostTimeLoop() override { structure_->PostTimeLoop(); }

    //@}


    //! @name Solver calls
    //@{

    /// nonlinear solve
    INPAR::STR::ConvergenceStatus Solve() override { return structure_->Solve(); }

    //! linear structure solve with just an interface load
    Teuchos::RCP<Epetra_Vector> SolveRelaxationLinear() override
    {
      return structure_->SolveRelaxationLinear();
    }

    /// get the linear solver object used for this field
    Teuchos::RCP<CORE::LINALG::Solver> LinearSolver() override
    {
      return structure_->LinearSolver();
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
    void PreparePartitionStep() override { structure_->PreparePartitionStep(); }

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
    Teuchos::RCP<Epetra_Vector> Freact() override { return structure_->Freact(); }

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
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest() override
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
    bool HaveModel(INPAR::STR::ModelType model) override { return structure_->HaveModel(model); }

    /// return model evaluator
    STR::MODELEVALUATOR::Generic& ModelEvaluator(INPAR::STR::ModelType mtype) override
    {
      return structure_->ModelEvaluator(mtype);
    }

    [[nodiscard]] bool HasFinalStateBeenWritten() const override
    {
      return structure_->HasFinalStateBeenWritten();
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

    void PrepareTimeStep() override;

    //! Evaluate() routine that can handle NOX step increments by computing the
    //! last iteration increment needed for structural Evaluate() call
    void Evaluate(Teuchos::RCP<const Epetra_Vector> disstepinc) override;

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
}  // namespace ADAPTER

BACI_NAMESPACE_CLOSE

#endif

/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the structural time integration

\level 1

\maintainer Anh-Tu Vuong

*/
/*----------------------------------------------------------------------*/

#ifndef ADAPTER_STRUCTURE_WRAPPER_H
#define ADAPTER_STRUCTURE_WRAPPER_H

#include "ad_str_structure.H"

namespace ADAPTER
{
  /// Just wrap, do nothing new, meant to be derived from
  class StructureWrapper : public Structure
  {
   public:
    /// constructor
    explicit StructureWrapper(Teuchos::RCP<Structure> structure) : structure_(structure) {}


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
    virtual void Setup() { structure_->Setup(); };

    //@}

    //! @name Vector access
    //@{

    /// initial guess of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() { return structure_->InitialGuess(); }

    /// right-hand-side of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> RHS() { return structure_->RHS(); }

    /// unknown displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp() const { return structure_->Dispnp(); }

    /// known displacements at \f$t_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispn() const { return structure_->Dispn(); }

    /// unknown velocity at \f$t_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Velnp() const { return structure_->Velnp(); }

    /// known velocity at \f$t_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Veln() const { return structure_->Veln(); }

    /// known velocity at \f$t_{n-1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Velnm() const { return structure_->Velnm(); }

    /// unknown acceleration at \f$t_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Accnp() const { return structure_->Accnp(); }

    /// known acceleration at \f$t_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Accn() const { return structure_->Accn(); }

    //@}


    //! @name Misc
    //@{

    /// dof map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() { return structure_->DofRowMap(); }

    /// dof map of vector of unknowns for multiple dof sets
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds)
    {
      return structure_->DofRowMap(nds);
    }

    /// view of dof map of vector of vector of unknowns
    virtual const Epetra_Map* DofRowMapView() { return structure_->DofRowMapView(); }

    /// domain map of system matrix
    virtual const Epetra_Map& DomainMap() const { return structure_->DomainMap(); }

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return structure_->SystemMatrix(); }

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix()
    {
      return structure_->BlockSystemMatrix();
    }

    /// switch structure field to block matrix
    virtual void UseBlockMatrix(Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps,
        Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps)
    {
      structure_->UseBlockMatrix(domainmaps, rangemaps);
      return;
    }

    // access to contact/meshtying bridge
    virtual Teuchos::RCP<CONTACT::MeshtyingContactBridge> MeshtyingContactBridge()
    {
      return structure_->MeshtyingContactBridge();
    }

    // access to locsys manager
    virtual Teuchos::RCP<DRT::UTILS::LocsysManager> LocsysManager()
    {
      return structure_->LocsysManager();
    }

    /// direct access to discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization()
    {
      return structure_->Discretization();
    }

    /// read only access to discretization
    virtual Teuchos::RCP<const DRT::Discretization> GetDiscretization() const
    {
      return structure_->Discretization();
    }

    /// are there any algebraic constraints?
    virtual bool HaveConstraint() { return structure_->HaveConstraint(); }

    /// are there any spring dashpot BCs?
    virtual bool HaveSpringDashpot() { return structure_->HaveSpringDashpot(); }

    /// get constraint manager defined in the structure
    virtual const Teuchos::RCP<UTILS::ConstrManager> GetConstraintManager()
    {
      return structure_->GetConstraintManager();
    }

    /// get constraint manager defined in the structure
    virtual const Teuchos::RCP<UTILS::SpringDashpotManager> GetSpringDashpotManager()
    {
      return structure_->GetSpringDashpotManager();
    }

    /// get type of thickness scaling for thin shell structures
    virtual INPAR::STR::STC_Scale GetSTCAlgo() { return structure_->GetSTCAlgo(); }

    /// access to scaling matrix for STC
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetSTCMat() { return structure_->GetSTCMat(); }

    /// Return MapExtractor for Dirichlet boundary conditions
    virtual const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
    {
      return structure_->GetDBCMapExtractor();
    }

    /// expand dirichlet bc map
    virtual void AddDirichDofs(const Teuchos::RCP<const Epetra_Map> maptoadd)
    {
      structure_->AddDirichDofs(maptoadd);
    };

    /// contract dirichlet bc map
    virtual void RemoveDirichDofs(const Teuchos::RCP<const Epetra_Map> maptoremove)
    {
      structure_->RemoveDirichDofs(maptoremove);
    };

    /// reset step and state vectors
    virtual void Reset() { structure_->Reset(); }

    /// reset last time step, needed for time step size adaptivity in FSI
    virtual void ResetStep() { structure_->ResetStep(); }

    //@}


    /// @name Time step helpers
    //@{

    /// return time integration factor
    virtual double TimIntParam() const { return structure_->TimIntParam(); }

    //! Sets the current time \f$t_{n}\f$
    virtual void SetTime(const double time)
    {
      structure_->SetTime(time);
      return;
    }

    //! Sets the target time \f$t_{n+1}\f$ of this time step
    virtual void SetTimen(const double time) { structure_->SetTimen(time); }

    //! Sets the target step \f$n\f$
    virtual void SetStep(int step)
    {
      structure_->SetStep(step);
      return;
    }

    //! Sets the target step \f$n+1\f$
    virtual void SetStepn(int step)
    {
      structure_->SetStepn(step);
      return;
    }

    //! Return current time \f$t_{n}\f$
    virtual double TimeOld() const { return structure_->TimeOld(); }

    //! Return target time \f$t_{n+1}\f$
    virtual double Time() const { return structure_->Time(); }

    /// get upper limit of time range of interest
    virtual double GetTimeEnd() const { return structure_->GetTimeEnd(); }

    //! Set upper limit of time range of interest //HACK for parameter continuation
    virtual void SetTimeEnd(double timemax) { structure_->SetTimeEnd(timemax); }

    /// get time step size \f$\Delta t_n\f$
    virtual double Dt() const { return structure_->Dt(); }

    /// Return current step number $n$
    virtual int StepOld() const { return structure_->StepOld(); }

    /// Return current step number $n+1$
    virtual int Step() const { return structure_->Step(); }

    /// get number of time steps
    virtual int NumStep() const { return structure_->NumStep(); }

    /// integrate from t1 to t2
    virtual int Integrate() { return structure_->Integrate(); }

    //! do something in case nonlinear solution does not converge for some reason
    virtual INPAR::STR::ConvergenceStatus PerformErrorAction(
        INPAR::STR::ConvergenceStatus nonlinsoldiv)
    {
      return structure_->PerformErrorAction(nonlinsoldiv);
    }

    /// tests if there are more time steps to do
    virtual bool NotFinished() const { return structure_->NotFinished(); }

    /// set time step size
    virtual void SetDt(const double dtnew)
    {
      structure_->SetDt(dtnew);
      return;
    }

    /// start new time step
    virtual void PrepareTimeStep() { structure_->PrepareTimeStep(); }

    /// update displacment
    void UpdateStateIncrementally(
        Teuchos::RCP<const Epetra_Vector> disi  ///< iterative solution increment
        ) override
    {
      structure_->UpdateStateIncrementally(disi);
    }

    /// update displacement and evaluate elements (implicit only)
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc)
    {
      structure_->Evaluate(disiterinc);
    }

    /// don't update displacement but evaluate elements (implicit only)
    virtual void Evaluate() { structure_->Evaluate(); }

    /// update at time step end
    virtual void Update() { structure_->Update(); }

    /// update at time step end
    virtual void Update(const double endtime) { structure_->Update(endtime); }

    /// resize MStep objects for AB2
    virtual void ResizeMStepTimAda() { structure_->ResizeMStepTimAda(); }

    /// update iteration; add residual increment to Lagrange multipliers stored in Constraint
    /// manager
    virtual void UpdateIterIncrConstr(Teuchos::RCP<Epetra_Vector> lagrincr)
    {
      structure_->UpdateIterIncrConstr(lagrincr);
    }

    /// update iteration; add residual increment to pressures stored in 0D cardiovascular manager
    virtual void UpdateIterIncrCardiovascular0D(Teuchos::RCP<Epetra_Vector> presincr)
    {
      structure_->UpdateIterIncrCardiovascular0D(presincr);
    }

    /// access to output object
    virtual Teuchos::RCP<IO::DiscretizationWriter> DiscWriter() { return structure_->DiscWriter(); }

    /// prepare output (i.e. calculate stresses, strains, energies)
    virtual void PrepareOutput()
    {
      structure_->PrepareOutput();
      return;
    }

    // Get restart data
    virtual void GetRestartData(Teuchos::RCP<int> step, Teuchos::RCP<double> time,
        Teuchos::RCP<Epetra_Vector> disn, Teuchos::RCP<Epetra_Vector> veln,
        Teuchos::RCP<Epetra_Vector> accn, Teuchos::RCP<std::vector<char>> elementdata,
        Teuchos::RCP<std::vector<char>> nodedata)
    {
      structure_->GetRestartData(step, time, disn, veln, accn, elementdata, nodedata);
      return;
    }

    /// output results
    virtual void Output(bool forced_writerestart = false)
    {
      structure_->Output(forced_writerestart);
      return;
    }

    /// Write Gmsh output for structural field
    virtual void writeGmshStrucOutputStep()
    {
      structure_->writeGmshStrucOutputStep();
      return;
    }

    /// output results to screen
    virtual void PrintStep()
    {
      structure_->PrintStep();
      return;
    }

    /// read restart information for given time step
    virtual void ReadRestart(const int step)
    {
      structure_->ReadRestart(step);
      return;
    }

    /// set restart information for parameter continuation
    virtual void SetRestart(int step, double time, Teuchos::RCP<Epetra_Vector> disn,
        Teuchos::RCP<Epetra_Vector> veln, Teuchos::RCP<Epetra_Vector> accn,
        Teuchos::RCP<std::vector<char>> elementdata, Teuchos::RCP<std::vector<char>> nodedata)
    {
      structure_->SetRestart(step, time, disn, veln, accn, elementdata, nodedata);
      return;
    }

    /// set the state of the nox group and the global state data container (implicit only)
    virtual void SetState(const Teuchos::RCP<Epetra_Vector>& x) { structure_->SetState(x); }

    /// set evaluation action
    virtual void SetActionType(const DRT::ELEMENTS::ActionType& action)
    {
      structure_->SetActionType(action);
    }

    /// wrapper for things that should be done before PrepareTimeStep is called
    virtual void PrePredict() { structure_->PrePredict(); }

    /// wrapper for things that should be done after PrepareTimeStep is called
    virtual void PostPredict() { structure_->PostPredict(); }

    /// wrapper for things that should be done before solving the nonlinear iterations
    virtual void PreSolve() { structure_->PreSolve(); }

    /// wrapper for things that should be done after solving the nonlinear iterations
    virtual void PostSolve() { structure_->PostSolve(); }

    /// wrapper for things that should be done before updating
    virtual void PreUpdate() { structure_->PreUpdate(); }

    /// wrapper for things that should be done after solving the update
    virtual void PostUpdate() { structure_->PostUpdate(); }

    /// wrapper for things that should be done before the output
    virtual void PreOutput() { structure_->PreOutput(); }

    /// wrapper for things that should be done after the output
    virtual void PostOutput() { structure_->PostOutput(); }

    /// wrapper for things that should be done after the actual time loop is finished
    virtual void PostTimeLoop() { structure_->PostTimeLoop(); }

    //@}


    //! @name Solver calls
    //@{

    /// nonlinear solve
    virtual INPAR::STR::ConvergenceStatus Solve() { return structure_->Solve(); }

    //! linear structure solve with just an interface load
    virtual Teuchos::RCP<Epetra_Vector> SolveRelaxationLinear()
    {
      return structure_->SolveRelaxationLinear();
    }

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() { return structure_->LinearSolver(); }

    //@}


    //! @name volume coupled specific methods
    //@{

    /// set forces due to interface with fluid, the force is expected external-force-like
    virtual void SetForceInterface(Teuchos::RCP<Epetra_MultiVector> iforce)
    {
      structure_->SetForceInterface(iforce);
    }

    //! specific method for iterative staggered partitioned TSI
    //! will be obsolete after switch to new structural timint.
    virtual void PreparePartitionStep() { structure_->PreparePartitionStep(); }

    //@}


    //! @name Write access to field solution variables at \f$t^{n+1}\f$
    //@{

    /// write access to extract displacements at \f$t^{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispnp()
    {
      return structure_->WriteAccessDispnp();
    }

    /// write access to extract velocities at \f$t^{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessVelnp()
    {
      return structure_->WriteAccessVelnp();
    }

    /// write access to extract displacements at \f$t^{n}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispn()
    {
      return structure_->WriteAccessDispn();
    }

    /// write access to extract velocities at \f$t^{n}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessVeln() { return structure_->WriteAccessVelnp(); }

    //@}

    /// extract rhs (used to calculate reaction force for post-processing)
    virtual Teuchos::RCP<Epetra_Vector> Freact() { return structure_->Freact(); }

    //! @name Structure with ale specific methods
    //@{

    /// material displacements (structure with ale)
    virtual Teuchos::RCP<Epetra_Vector> DispMat() { return structure_->DispMat(); }

    /// apply material displacements to structure field (structure with ale)
    virtual void ApplyDisMat(Teuchos::RCP<Epetra_Vector> dismat)
    {
      structure_->ApplyDisMat(dismat);
    }

    //@}


    /// create result test for encapsulated structure algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest()
    {
      return structure_->CreateFieldTest();
    }

    //! @name Biofilm specific methods
    //@{

    virtual void SetStrGrDisp(Teuchos::RCP<Epetra_Vector> struct_growth_disp)
    {
      structure_->SetStrGrDisp(struct_growth_disp);
      return;
    }

    //@}

    /// bool indicating if micro material is used
    virtual bool HaveMicroMat() { return structure_->HaveMicroMat(); }

    /// do we have this model
    virtual bool HaveModel(INPAR::STR::ModelType model) { return structure_->HaveModel(model); }

    /// return model evaluator
    virtual ::STR::MODELEVALUATOR::Generic& ModelEvaluator(INPAR::STR::ModelType mtype)
    {
      return structure_->ModelEvaluator(mtype);
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

    virtual void PrepareTimeStep();

    //! Evaluate() routine that can handle NOX step increments by computing the
    //! last iteration increment needed for structural Evaluate() call
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disstepinc);

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

#endif

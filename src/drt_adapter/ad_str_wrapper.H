/*----------------------------------------------------------------------*/
/*!
\file ad_str_wrapper.H

\brief Wrapper for the structural time integration

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
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


    //! @name Vector access
    //@{

    /// initial guess of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess()
    { return structure_->InitialGuess(); }

    /// right-hand-side of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> RHS()
    { return structure_->RHS(); }

    /// unknown displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp()
    { return structure_->Dispnp(); }

    /// known displacements at \f$t_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispn()
    { return structure_->Dispn(); }

    /// return displacements at all the time steps
    virtual Teuchos::RCP<STR::TimIntMStep<Epetra_Vector> > DispMStep()
    { return structure_->DispMStep(); }

    /// return times at all the time steps
    virtual Teuchos::RCP<STR::TimIntMStep<double> > TimeMStep()
    { return structure_->TimeMStep(); }

    //@}


    //! @name Misc
    //@{

    /// dof map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap()
    { return structure_->DofRowMap(); }

    /// dof map of vector of unknowns for multiple dof sets
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds)
    { return structure_->DofRowMap(nds); }

    /// domain map of system matrix
    virtual const Epetra_Map& DomainMap()
    { return structure_->DomainMap(); }

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix()
    { return structure_->SystemMatrix(); }

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix()
    { return structure_->BlockSystemMatrix(); }

    /// switch structure field to block matrix
    virtual void UseBlockMatrix(Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps, Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps)
    { return structure_->UseBlockMatrix(domainmaps, rangemaps); }

    /// recalculate structural matrices for tsi simulations
    virtual void TSIMatrix()
    { return structure_->TSIMatrix(); }

    // access to contact/meshtying manager
    virtual Teuchos::RCP<MORTAR::ManagerBase> ContactManager()
    { return structure_->ContactManager(); }

    // access to locsys manager
    virtual Teuchos::RCP<DRT::UTILS::LocsysManager> LocsysManager()
    { return structure_->LocsysManager(); }

    /// direct access to discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization()
    { return structure_->Discretization(); }

    /// are there any algebraic constraints?
    virtual bool HaveConstraint()
    { return structure_->HaveConstraint(); }

    /// get constraint manager defined in the structure
    virtual const Teuchos::RCP< UTILS::ConstrManager> GetConstraintManager()
    { return structure_->GetConstraintManager(); }

    /// get type of thickness scaling for thin shell structures
    virtual INPAR::STR::STC_Scale GetSTCAlgo()
    { return structure_->GetSTCAlgo(); }

    /// access to scaling matrix for STC
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetSTCMat()
    { return structure_->GetSTCMat(); }

    /// Return MapExtractor for Dirichlet boundary conditions
    virtual const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
    { return structure_->GetDBCMapExtractor(); }

    /// reset step and state vectors
    virtual void Reset()
    { structure_->Reset(); }

    /// integrate a single time step
    virtual int IntegrateStep()
    { return structure_->IntegrateStep(); }

    /*
    \brief Shorten the Dirichlet DOF set

    The method shortens the DOF set (map) which contains the DOFs
    subjected to Dirichlet boundary conditions. For instance, the method is
    called by the monolithic FSI where fluid boundary conditions on the
    FSI interface have to be removed if it is a fluid split.

    */
    virtual void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove)
    { return structure_->RemoveDirichCond(maptoremove); }

    //@}


    /// @name Time step helpers
    //@{

    /// return time integration factor
    virtual double TimIntParam()
    { return structure_->TimIntParam(); }

    //! Sets the current time needed for the inverse analysis
    virtual void SetTimen(double time) { structure_->SetTimen(time); }

    /// get the current time
    virtual double GetTime() const
    { return structure_->GetTime(); }

    /// get upper limit of time range of interest
    virtual double GetTimeEnd() const
    { return structure_->GetTimeEnd(); }

    //! Set upper limit of time range of interest //HACK for parameter continuation
    virtual void SetTimeEnd(double timemax)
    { structure_->SetTimeEnd(timemax); }

    /// get time step size \f$\Delta t_n\f$
    virtual double GetTimeStepSize() const
    { return structure_->GetTimeStepSize(); }

    /// get current time step \f$n\f$ (sometimes called \f$k\f$ [but not in Baci])
    virtual int GetTimeStep() const
    { return structure_->GetTimeStep(); }

    /// get number of time steps
    virtual int GetTimeNumStep() const
    { return structure_->GetTimeNumStep(); }

    /// integrate from t1 to t2
    virtual int Integrate() { return structure_->Integrate(); }

    /// set time step size
    virtual void SetTimeStepSize(double timestepsize)
    { return structure_->SetTimeStepSize(timestepsize); }

    /// set time step size and time
    virtual void SetInitialTimeStepAndTime(double timestepsize)
    { return structure_->SetInitialTimeStepAndTime(timestepsize); }

    /// start new time step
    virtual void PrepareTimeStep()
    { structure_->PrepareTimeStep(); }

    /// update displacement and evaluate elements
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc)
    { structure_->Evaluate(disiterinc); }

    /// update at time step end
    virtual void Update()
    { structure_->Update(); }

    /// update iteration; add residual increment to Lagrange multipliers stored in Constraint manager
    virtual void UpdateIterIncrConstr(Teuchos::RCP<Epetra_Vector> lagrincr)
    { structure_->UpdateIterIncrConstr(lagrincr); }

    /// access to output object
    virtual Teuchos::RCP<IO::DiscretizationWriter> DiscWriter()
    { return structure_->DiscWriter();}

    /// prepare output (i.e. calculate stresses, strains, energies)
    virtual void PrepareOutput()
    { return structure_->PrepareOutput(); }

    // Get restart data
    virtual void GetRestartData
    (
      Teuchos::RCP<int> step,
      Teuchos::RCP<double> time,
      Teuchos::RCP<Epetra_Vector> disn,
      Teuchos::RCP<Epetra_Vector> veln,
      Teuchos::RCP<Epetra_Vector> accn,
      Teuchos::RCP<std::vector<char> > elementdata
    )
    {return structure_->GetRestartData(step,time, disn,veln,accn,elementdata );}

    /// output results
    virtual void Output(bool forced_writerestart = false)
    { return structure_->Output(forced_writerestart); }

    /// read restart information for given time step
    virtual void ReadRestart(int step)
    { return structure_->ReadRestart(step); }

     /// set restart information for parameter continuation
    virtual void SetRestart
    (
      int step,
      double time,
      Teuchos::RCP<Epetra_Vector> disn,
      Teuchos::RCP<Epetra_Vector> veln,
      Teuchos::RCP<Epetra_Vector> accn,
      Teuchos::RCP<std::vector<char> > elementdata
    )
    {return structure_->SetRestart(step,time,disn,veln,accn,elementdata);}

    /// evaluate reference state (for frictional contact)
    virtual void EvaluateReferenceState()
    { structure_->EvaluateReferenceState(); }

    //@}


    //! @name Solver calls
    //@{

    /// nonlinear solve
    virtual int Solve()
    { return structure_->Solve(); }

    //! linear structure solve with just an interface load
    virtual Teuchos::RCP<Epetra_Vector> SolveRelaxationLinear()
    { return structure_->SolveRelaxationLinear(); }

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver()
    { return structure_->LinearSolver(); }

    //@}


    //! @name volume coupled specific methods
    //@{

    //! set volume coupling state (e.g. temperatures for tsi, fluid velocities for poro) in structural discretization
    virtual void SetCouplingState()
    { structure_->SetCouplingState(); }

    //! apply volume coupling state (e.g. temperatures for tsi, fluid velocities for poro) in structural time integrator
    virtual void ApplyCouplingState(Teuchos::RCP<const Epetra_Vector> velnp,   const std::string& name, unsigned dofset=1)
    { structure_->ApplyCouplingState(velnp,name,dofset); }

    /// set forces due to interface with fluid, the force is expected external-force-like
    virtual void SetForceInterface(Teuchos::RCP<Epetra_MultiVector> iforce)
    { structure_->SetForceInterface(iforce); }

    //! specific method for iterative staggered partitioned TSI
    virtual void PreparePartitionStep()
    { structure_->PreparePartitionStep(); }

    //@}


    //! @name Extract methods
    //@{

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> ExtractDispnp()
    { return structure_->ExtractDispnp(); }

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> ExtractDispn()
    { return structure_->ExtractDispn(); }

    /// extract velocities
    virtual Teuchos::RCP<Epetra_Vector> ExtractVeln()
    { return structure_->ExtractVeln(); }

    /// extract velocities
    virtual Teuchos::RCP<Epetra_Vector> ExtractVelnp()
    { return structure_->ExtractVelnp(); }

    /// extract velocities
    virtual Teuchos::RCP<Epetra_Vector> ExtractVelaf()
    { return structure_->ExtractVelaf(); }

    /// extract accelerations
    virtual Teuchos::RCP<Epetra_Vector> ExtractAccn()
    { return structure_->ExtractAccn(); }

    /// extract rhs (used to calculate reaction force for post-processing)
    virtual Teuchos::RCP<Epetra_Vector> Freact()
    { return structure_->Freact(); }

    //@}

    //! @name Structure with ale specific methods
    //@{

    /// material displacements (structure with ale)
    virtual Teuchos::RCP<Epetra_Vector> DispMat()
    { return structure_->DispMat(); }

    /// apply material displacements to structure field (structure with ale)
    virtual void ApplyDisMat(Teuchos::RCP<Epetra_Vector> dismat)
    { structure_->ApplyDisMat(dismat); }

    //@}


    /// create result test for encapsulated structure algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest()
    { return structure_->CreateFieldTest(); }

    //! @name Biofilm specific methods
    //@{

    virtual void SetStrGrDisp(Teuchos::RCP<Epetra_Vector> struct_growth_disp)
    { return structure_->SetStrGrDisp(struct_growth_disp); }

    //@}

  protected:

    Teuchos::RCP<Structure> structure_;   ///< underlying structural time integration
  };


  /// Calculate increments from absolute values
  class StructureNOXCorrectionWrapper : public StructureWrapper
  {
  public:

    explicit StructureNOXCorrectionWrapper(Teuchos::RCP<Structure> structure)
      : StructureWrapper(structure) {}

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
}

#endif

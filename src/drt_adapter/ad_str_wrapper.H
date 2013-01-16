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

    /// unknown displacements at \f$t_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispn()
    { return structure_->Dispn(); }

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
    virtual void IntegrateStep()
    { structure_->IntegrateStep(); }

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

    /// get the current time
    virtual double GetTime() const
    { return structure_->GetTime(); }

    /// get upper limit of time range of interest
    virtual double GetTimeEnd() const
    { return structure_->GetTimeEnd(); }

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
    virtual void Integrate() { structure_->Integrate(); }

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

    /// output results
    virtual void Output()
    { return structure_->Output(); }

    /// read restart information for given time step
    virtual void ReadRestart(int step)
    { return structure_->ReadRestart(step); }

    /// evaluate reference state (for frictional contact)
    virtual void EvaluateReferenceState()
    { structure_->EvaluateReferenceState(); }

    //@}


    //! @name Solver calls
    //@{

    /// nonlinear solve
    virtual void Solve()
    { return structure_->Solve(); }

    //! linear structure solve with just an interface load
    virtual Teuchos::RCP<Epetra_Vector> SolveRelaxationLinear()
    { return structure_->SolveRelaxationLinear(); }

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver()
    { return structure_->LinearSolver(); }

    //@}


    //! @name TSI specific methods
    //@{

    //! apply temperatures  for TSI
    virtual void ApplyTemperatures(Teuchos::RCP<const Epetra_Vector> temp)
    { structure_->ApplyTemperatures(temp); }

    /// set forces due to interface with fluid, the force is expected external-force-like
    virtual void SetForceInterface(Teuchos::RCP<Epetra_Vector> iforce)
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


    //! @name Poroelasticity specific methods
    //@{

    virtual void ApplyVelAndPress(Teuchos::RCP<const Epetra_Vector> velnp)
    { structure_->ApplyVelAndPress(velnp); }

    //@}

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

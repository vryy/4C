/*----------------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the ALE time integration

\level 2

\maintainer Matthias Mayr
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifndef ADAPTER_ALE_WRAPPER_H
#define ADAPTER_ALE_WRAPPER_H

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale.H"

#include "../drt_ale/ale_utils_mapextractor.H"

/*----------------------------------------------------------------------------*/
/* class definitions */
namespace ADAPTER
{
  /*! \brief Just a wrapper that does nothing, meant to be derived from
   *
   *  This wrapper just encapsulated the ADAPTER::Ale and implements all
   *  routines that are pure virtual in ADAPTER::Ale. For a specific ALE adapter
   *  just derive from this one and overload those routines you need with your
   *  problem specific routine.
   *
   *  \author mayr.mt \date 10/2014
   */
  class AleWrapper : public Ale
  {
   public:
    //! constructor
    explicit AleWrapper(Teuchos::RCP<Ale> ale) : ale_(ale) {}

    //! @name Vector access
    //@{

    //! initial guess of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() const { return ale_->InitialGuess(); }

    //! right-hand-side of Newton's method
    virtual Teuchos::RCP<const Epetra_Vector> RHS() const { return ale_->RHS(); }

    //! unknown displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp() const { return ale_->Dispnp(); }

    //! known displacements at \f$t_{n}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispn() const { return ale_->Dispn(); }

    //@}

    //! @name Misc
    //@{

    //! dof map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() const { return ale_->DofRowMap(); }

    //! direct access to system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return ale_->SystemMatrix(); }

    //! direct access to system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix()
    {
      return ale_->BlockSystemMatrix();
    }

    //! access to locsys manager
    virtual Teuchos::RCP<DRT::UTILS::LocsysManager> LocsysManager()
    {
      return ale_->LocsysManager();
    }

    //! direct access to discretization
    virtual Teuchos::RCP<const DRT::Discretization> Discretization() const
    {
      return ale_->Discretization();
    }

    /// writing access to discretization
    virtual Teuchos::RCP<DRT::Discretization> WriteAccessDiscretization()
    {
      return ale_->WriteAccessDiscretization();
    }

    //! Return MapExtractor for Dirichlet boundary conditions
    virtual const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
    {
      return ale_->GetDBCMapExtractor(ALE::UTILS::MapExtractor::dbc_set_std);
    }

    //! Return MapExtractor for Dirichlet boundary conditions in case of non-standard Dirichlet sets
    virtual const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor(
        ALE::UTILS::MapExtractor::AleDBCSetType dbc_type  ///< type of dbc set
    )
    {
      return ale_->GetDBCMapExtractor(dbc_type);
    }

    //! reset state vectors to zero
    virtual void Reset() { ale_->Reset(); }

    //! reset last time step, needed for time step size adaptivity of FSI
    virtual void ResetStep() { ale_->ResetStep(); }

    //@}

    //! @name Time step helpers
    //@{

    virtual void ResetTime(const double dtold) { ale_->ResetTime(dtold); }
    //! Return target time \f$t_{n+1}\f$
    virtual double Time() const { return ale_->Time(); }

    //! Return target step counter \f$step_{n+1}\f$
    virtual double Step() const { return ale_->Step(); }

    //! get time step size \f$\Delta t_n\f$
    virtual double Dt() const { return ale_->Dt(); }

    //! integrate from t1 to t2
    virtual int Integrate() { return ale_->Integrate(); }

    virtual void TimeStep(
        ALE::UTILS::MapExtractor::AleDBCSetType dbc_type = ALE::UTILS::MapExtractor::dbc_set_std)
    {
      ale_->TimeStep(dbc_type);
      return;
    }

    //! set time step size
    virtual void SetDt(const double dtnew  ///< new time step size (to be set)
    )
    {
      ale_->SetDt(dtnew);
    }

    //! Set time and step
    void SetTimeStep(const double time,  ///< simulation time (to be set)
        const int step                   ///< step number (to be set)
    )
    {
      ale_->SetTimeStep(time, step);
    }

    //! start new time step
    virtual void PrepareTimeStep() { ale_->PrepareTimeStep(); }

    //! update displacement and evaluate elements
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> stepinc =
                              Teuchos::null  ///< step increment such that \f$ x_{n+1}^{k+1} =
                                             ///< x_{n}^{converged}+ stepinc \f$
    )
    {
      Evaluate(stepinc, ALE::UTILS::MapExtractor::dbc_set_std);
    }

    //! update displacement and evaluate elements
    void Evaluate(
        Teuchos::RCP<const Epetra_Vector> stepinc,  ///< step increment such that \f$ x_{n+1}^{k+1}
                                                    ///< = x_{n}^{converged}+ stepinc \f$
        ALE::UTILS::MapExtractor::AleDBCSetType
            dbc_type  ///< application-specific type of Dirichlet set
    )
    {
      ale_->Evaluate(stepinc, dbc_type);
    }

    //! update at time step end
    virtual void Update() { ale_->Update(); }

    //! update at time step end
    virtual void UpdateIter() { ale_->UpdateIter(); }

    //! output results
    virtual void Output() { return ale_->Output(); }


    //! read restart information for given time step \p step
    virtual void ReadRestart(const int step  ///< step number to read restart from
    )
    {
      return ale_->ReadRestart(step);
    }

    //@}

    /// setup Dirichlet boundary condition map extractor
    virtual void SetupDBCMapEx(
        ALE::UTILS::MapExtractor::AleDBCSetType dbc_type =
            ALE::UTILS::MapExtractor::dbc_set_std,  //!< application-specific type of Dirichlet set
        Teuchos::RCP<const ALE::UTILS::MapExtractor> interface =
            Teuchos::null,  //!< interface for creation of additional, application-specific
                            //!< Dirichlet map extractors
        Teuchos::RCP<const ALE::UTILS::XFluidFluidMapExtractor> xff_interface =
            Teuchos::null  //!< interface for creation of a Dirichlet map extractor, taylored to
                           //!< XFFSI
    )
    {
      ale_->SetupDBCMapEx(dbc_type, interface, xff_interface);
    }

    //! @name Solver calls
    //@{

    //! nonlinear solve
    virtual int Solve() { return ale_->Solve(); }

    //! Access to linear solver of ALE field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() { return ale_->LinearSolver(); }

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Preconditioner> ConstPreconditioner()
    {
      return ale_->ConstPreconditioner();
    }

    //@}

    //! @name Write access to field solution variables at \f$t^{n+1}\f$
    //@{

    //! write access to extract displacements at \f$t^{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispnp() const
    {
      return ale_->WriteAccessDispnp();
    }

    //@}

    //! create result test for encapsulated structure algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest() { return ale_->CreateFieldTest(); }

    /*! \brief Create Systemmatrix
     *
     * We allocate the LINALG object just once, the result is an empty LINALG object.
     * Evaluate has to be called separately.
     *
     */
    void CreateSystemMatrix(
        Teuchos::RCP<const ALE::UTILS::MapExtractor> interface =
            Teuchos::null  ///< Blocksparsematrix if an interface is passed, SparseMatrix otherwise
    )
    {
      ale_->CreateSystemMatrix(interface);
    }

    //! update slave dofs for fsi simulations with ale mesh tying
    void UpdateSlaveDOF(Teuchos::RCP<Epetra_Vector>& a) { ale_->UpdateSlaveDOF(a); }

   private:
    Teuchos::RCP<Ale> ale_;  //!< underlying ALE time integration
  };

  //! Calculate increments from absolute values
  class AleNOXCorrectionWrapper : public AleWrapper  // ToDo (mayr) Do we really need this?
  {
   public:
    explicit AleNOXCorrectionWrapper(Teuchos::RCP<Ale> ale) : AleWrapper(ale) {}

    //! Prepare time step
    virtual void PrepareTimeStep();

    /*! \brief Evaluate() routine that can handle NOX step increments
     *
     *  We deal with NOX step increments by computing the last iteration increment
     *  needed for the ALE Evaluate() call.
     *
     *  The field solver always expects an iteration increment only. And
     *  there are Dirichlet conditions that need to be preserved. So take
     *  the sum of increments we get from NOX and apply the latest iteration
     *  increment only.
     *  Naming:
     *
     *  \f$x^{n+1}_{i+1} = x^{n+1}_i + iterinc\f$  (sometimes referred to as residual increment),
     * and
     *
     *  \f$x^{n+1}_{i+1} = x^n + stepinc\f$
     *
     *  \author mayr.mt \date 10/2014
     */
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> stepinc  ///< step increment
    );

   private:
    //! sum of displacement increments already applied,
    //!
    //! there are two increments around
    //!
    //! x^n+1_i+1 = x^n+1_i + stepinc  (also referred to as residual increment)
    //!
    //! x^n+1_i+1 = x^n     + disstepinc
    Teuchos::RCP<Epetra_Vector> stepinc_;
  };
}  // namespace ADAPTER

#endif  // #ifndef ADAPTER_ALE_WRAPPER_H

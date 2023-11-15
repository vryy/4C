/*----------------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the ALE time integration

\level 2

 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifndef BACI_ADAPTER_ALE_WRAPPER_H
#define BACI_ADAPTER_ALE_WRAPPER_H

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "baci_adapter_ale.H"
#include "baci_ale_utils_mapextractor.H"

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
    Teuchos::RCP<const Epetra_Vector> InitialGuess() const override { return ale_->InitialGuess(); }

    //! right-hand-side of Newton's method
    Teuchos::RCP<const Epetra_Vector> RHS() const override { return ale_->RHS(); }

    //! unknown displacements at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> Dispnp() const override { return ale_->Dispnp(); }

    //! known displacements at \f$t_{n}\f$
    Teuchos::RCP<const Epetra_Vector> Dispn() const override { return ale_->Dispn(); }

    //@}

    //! @name Misc
    //@{

    //! dof map of vector of unknowns
    Teuchos::RCP<const Epetra_Map> DofRowMap() const override { return ale_->DofRowMap(); }

    //! direct access to system matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> SystemMatrix() override
    {
      return ale_->SystemMatrix();
    }

    //! direct access to system matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> BlockSystemMatrix() override
    {
      return ale_->BlockSystemMatrix();
    }

    //! access to locsys manager
    Teuchos::RCP<DRT::UTILS::LocsysManager> LocsysManager() override
    {
      return ale_->LocsysManager();
    }

    //! direct access to discretization
    Teuchos::RCP<const DRT::Discretization> Discretization() const override
    {
      return ale_->Discretization();
    }

    /// writing access to discretization
    Teuchos::RCP<DRT::Discretization> WriteAccessDiscretization() override
    {
      return ale_->WriteAccessDiscretization();
    }

    //! Return MapExtractor for Dirichlet boundary conditions
    virtual Teuchos::RCP<const CORE::LINALG::MapExtractor> GetDBCMapExtractor()
    {
      return ale_->GetDBCMapExtractor(ALE::UTILS::MapExtractor::dbc_set_std);
    }

    //! Return MapExtractor for Dirichlet boundary conditions in case of non-standard Dirichlet sets
    Teuchos::RCP<const CORE::LINALG::MapExtractor> GetDBCMapExtractor(
        ALE::UTILS::MapExtractor::AleDBCSetType dbc_type  ///< type of dbc set
        ) override
    {
      return ale_->GetDBCMapExtractor(dbc_type);
    }

    //! reset state vectors to zero
    void Reset() override { ale_->Reset(); }

    //! reset last time step, needed for time step size adaptivity of FSI
    void ResetStep() override { ale_->ResetStep(); }

    //@}

    //! @name Time step helpers
    //@{

    void ResetTime(const double dtold) override { ale_->ResetTime(dtold); }
    //! Return target time \f$t_{n+1}\f$
    double Time() const override { return ale_->Time(); }

    //! Return target step counter \f$step_{n+1}\f$
    double Step() const override { return ale_->Step(); }

    //! get time step size \f$\Delta t_n\f$
    double Dt() const override { return ale_->Dt(); }

    //! integrate from t1 to t2
    int Integrate() override { return ale_->Integrate(); }

    void TimeStep(ALE::UTILS::MapExtractor::AleDBCSetType dbc_type =
                      ALE::UTILS::MapExtractor::dbc_set_std) override
    {
      ale_->TimeStep(dbc_type);
      return;
    }

    //! set time step size
    void SetDt(const double dtnew  ///< new time step size (to be set)
        ) override
    {
      ale_->SetDt(dtnew);
    }

    //! Set time and step
    void SetTimeStep(const double time,  ///< simulation time (to be set)
        const int step                   ///< step number (to be set)
        ) override
    {
      ale_->SetTimeStep(time, step);
    }

    //! start new time step
    void PrepareTimeStep() override { ale_->PrepareTimeStep(); }

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
        ) override
    {
      ale_->Evaluate(stepinc, dbc_type);
    }

    //! update at time step end
    void Update() override { ale_->Update(); }

    //! update at time step end
    void UpdateIter() override { ale_->UpdateIter(); }

    //! output results
    void Output() override { return ale_->Output(); }


    //! read restart information for given time step \p step
    void ReadRestart(const int step  ///< step number to read restart from
        ) override
    {
      return ale_->ReadRestart(step);
    }

    //@}

    /// setup Dirichlet boundary condition map extractor
    void SetupDBCMapEx(
        ALE::UTILS::MapExtractor::AleDBCSetType dbc_type =
            ALE::UTILS::MapExtractor::dbc_set_std,  //!< application-specific type of Dirichlet set
        Teuchos::RCP<const ALE::UTILS::MapExtractor> interface =
            Teuchos::null,  //!< interface for creation of additional, application-specific
                            //!< Dirichlet map extractors
        Teuchos::RCP<const ALE::UTILS::XFluidFluidMapExtractor> xff_interface =
            Teuchos::null  //!< interface for creation of a Dirichlet map extractor, taylored to
                           //!< XFFSI
        ) override
    {
      ale_->SetupDBCMapEx(dbc_type, interface, xff_interface);
    }

    //! @name Solver calls
    //@{

    //! nonlinear solve
    int Solve() override { return ale_->Solve(); }

    //! Access to linear solver of ALE field
    Teuchos::RCP<CORE::LINALG::Solver> LinearSolver() override { return ale_->LinearSolver(); }

    /// get the linear solver object used for this field
    Teuchos::RCP<CORE::LINALG::Preconditioner> ConstPreconditioner() override
    {
      return ale_->ConstPreconditioner();
    }

    //@}

    //! @name Write access to field solution variables at \f$t^{n+1}\f$
    //@{

    //! write access to extract displacements at \f$t^{n+1}\f$
    Teuchos::RCP<Epetra_Vector> WriteAccessDispnp() const override
    {
      return ale_->WriteAccessDispnp();
    }

    //@}

    //! create result test for encapsulated structure algorithm
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest() override { return ale_->CreateFieldTest(); }

    /*! \brief Create Systemmatrix
     *
     * We allocate the CORE::LINALG object just once, the result is an empty CORE::LINALG object.
     * Evaluate has to be called separately.
     *
     */
    void CreateSystemMatrix(
        Teuchos::RCP<const ALE::UTILS::MapExtractor> interface =
            Teuchos::null  ///< Blocksparsematrix if an interface is passed, SparseMatrix otherwise
        ) override
    {
      ale_->CreateSystemMatrix(interface);
    }

    //! update slave dofs for fsi simulations with ale mesh tying
    void UpdateSlaveDOF(Teuchos::RCP<Epetra_Vector>& a) override { ale_->UpdateSlaveDOF(a); }

   private:
    Teuchos::RCP<Ale> ale_;  //!< underlying ALE time integration
  };

  //! Calculate increments from absolute values
  class AleNOXCorrectionWrapper : public AleWrapper  // ToDo (mayr) Do we really need this?
  {
   public:
    explicit AleNOXCorrectionWrapper(Teuchos::RCP<Ale> ale) : AleWrapper(ale) {}

    //! Prepare time step
    void PrepareTimeStep() override;

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
    void Evaluate(Teuchos::RCP<const Epetra_Vector> stepinc  ///< step increment
        ) override;

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

#endif  // ADAPTER_ALE_WRAPPER_H

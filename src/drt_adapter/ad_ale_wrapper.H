/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_wrapper.H

 \brief Wrapper for the ALE time integration

 <pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifndef ADAPTER_ALE_WRAPPER_H
#define ADAPTER_ALE_WRAPPER_H

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale.H"

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
  explicit AleWrapper(Teuchos::RCP<Ale> ale)
  : ale_(ale)
  {
  }

  //! @name Vector access
  //@{

  //! initial guess of Newton's method
  virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() const
  {
    return ale_->InitialGuess();
  }

  //! right-hand-side of Newton's method
  virtual Teuchos::RCP<const Epetra_Vector> RHS() const
  {
    return ale_->RHS();
  }

  //! unknown displacements at \f$t_{n+1}\f$
  virtual Teuchos::RCP<const Epetra_Vector> Dispnp() const
  {
    return ale_->Dispnp();
  }

  //! known displacements at \f$t_{n}\f$
  virtual Teuchos::RCP<const Epetra_Vector> Dispn() const
  {
    return ale_->Dispn();
  }

  //@}

  //! @name Misc
  //@{

  //! dof map of vector of unknowns
  virtual Teuchos::RCP<const Epetra_Map> DofRowMap() const
  {
    return ale_->DofRowMap();
  }

  //! direct access to system matrix
  virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix()
  {
    return ale_->SystemMatrix();
  }

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
    return ale_->GetDBCMapExtractor();
  }

  //! reset state vectors to zero
  virtual void Reset()
  {
    ale_->Reset();
  }

  //! reset last time step, needed for time step size adaptivity of FSI
  virtual void ResetStep()
  {
    ale_->ResetStep();
  }

  //@}

  //! @name Time step helpers
  //@{

  virtual void ResetTime(const double dtold)
  {
    ale_->ResetTime(dtold);
  }
  //! Return target time \f$t_{n+1}\f$
  virtual double Time() const
  {
    return ale_->Time();
  }

  //! get time step size \f$\Delta t_n\f$
  virtual double Dt() const
  {
    return ale_->Dt();
  }

  //! integrate from t1 to t2
  virtual int Integrate()
  {
    return ale_->Integrate();
  }

  //! set time step size
  virtual void SetDt(const double dtnew)
  {
    ale_->SetDt(dtnew);
  }

  //! start new time step
  virtual void PrepareTimeStep()
  {
    ale_->PrepareTimeStep();
  }

  //! update displacement and evaluate elements
  virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc)
  {
    ale_->Evaluate(disiterinc);
  }

  //! update at time step end
  virtual void Update()
  {
    ale_->Update();
  }

  //! output results
  virtual void Output()
  {
    return ale_->Output();
  }


  //! read restart information for given time step #step
  virtual void ReadRestart(const int step)
  {
    return ale_->ReadRestart(step);
  }

  //@}

  /// displacements at the interface are prescribed
  virtual void ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp)
  {
    ale_->ApplyInterfaceDisplacements(idisp);
  }

  /// setup Dirichlet boundary condition map extractor
  virtual void SetupDBCMapEx(bool dirichletcond)
  {
    ale_->SetupDBCMapEx(dirichletcond);
  }

  //! @name Solver calls
  //@{

  //! nonlinear solve
  virtual int Solve()
  {
    return ale_->Solve();
  }

  //! Access to linear solver of ALE field
  virtual Teuchos::RCP<LINALG::Solver> LinearSolver()
  {
    return ale_->LinearSolver();
  }

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
  virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest()
  {
    return ale_->CreateFieldTest();
  }

  /*! \brief Create Systemmatrix
   *
   * We allocate the LINALG object just once, the result is an empty LINALG object.
   * Evaluate has to be called separately.
   *
   */
  virtual void CreateSystemMatrix(bool full = true ///< Sparsematrix (true) or Blocksparsematrix (false)
      )
  {
    ale_->CreateSystemMatrix(full);
  }

private:

  Teuchos::RCP<Ale> ale_; //!< underlying ALE time integration
};

//! Calculate increments from absolute values
class AleNOXCorrectionWrapper: public AleWrapper // ToDo (mayr) Do we really need this?
{
public:

  explicit AleNOXCorrectionWrapper(Teuchos::RCP<Ale> ale) :
      AleWrapper(ale)
  {
  }

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
   *  \f$x^{n+1}_{i+1} = x^{n+1}_i + iterinc\f$  (sometimes referred to as residual increment), and
   *
   *  \f$x^{n+1}_{i+1} = x^n + stepinc\f$
   *
   *  \author mayr.mt \date 10/2014
   */
  virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> stepinc);

private:

  //! sum of displacement increments already applied,
  //!
  //! there are two increments around
  //!
  //! x^n+1_i+1 = x^n+1_i + disiterinc  (also referred to as residual increment)
  //!
  //! x^n+1_i+1 = x^n     + disstepinc
  Teuchos::RCP<Epetra_Vector> stepinc_;
};
}

#endif // #ifndef ADAPTER_ALE_WRAPPER_H

/*----------------------------------------------------------------------*/
/*!
\file adapter_ale_springs_fixed_ref.H

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ADAPTER_ALE_SPRINGS_FIXED_REF_H
#define ADAPTER_ALE_SPRINGS_FIXED_REF_H

#include "adapter_ale.H"

namespace ADAPTER {


/// springs ALE algorithm
/*!
  Incremental springs algorithm with fixed initial reference and constant stiffness.

  When compared to the original springs, this algorithm shines, when periodic
  mesh deformation is expected. It does not "forget" its initial configuration,
  hence it can not degrade over time as the original springs algorithm.
  On the downside, the maximal mesh deformation that can be achieved is less than
  the original springs algorithm. If one time, large mesh deformation is expected,
  use the original springs.

  The implementation is copied from the classic_lin ALE implementation.
  Only the element call was modified to the springs_fixed_ref calculation.

  \author axel gerstenberger
 */
class AleSpringsFixedRef : public Ale
{
  friend class AleResultTest;

public:

  AleSpringsFixedRef(Teuchos::RCP<DRT::Discretization> actdis,
             Teuchos::RCP<LINALG::Solver> solver,
             Teuchos::RCP<ParameterList> params,
             Teuchos::RCP<IO::DiscretizationWriter> output,
             bool incremental,
             bool dirichletcond=true);

  /// build system matrix
  virtual void BuildSystemMatrix(bool full=true);

  /// a very simple time loop to be used for standalone ALE problems
  void Integrate();

  //! @name Algorithm core methods
  virtual void PrepareTimeStep();
  virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> ddisp);
  virtual void Solve();
  virtual void Update();
  virtual void Output();
  //@}

  /// get the linear solver object used for this field
  Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

  Teuchos::RCP<LINALG::Preconditioner> ConstPreconditioner() { return precond_; }

  void SetInterfaceMap(Teuchos::RCP<Epetra_Map> im);

  void ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp);
  void ApplyFreeSurfaceDisplacements(Teuchos::RCP<Epetra_Vector> fsdisp);

  virtual Teuchos::RCP<Epetra_Vector> ExtractDisplacement() const;

  //! @name access functions for monolithic FSI

  virtual Teuchos::RCP<Epetra_Vector> InitialGuess()
    { return residual_; }
  virtual Teuchos::RCP<Epetra_Vector> RHS() const
    { return residual_; }
  virtual Teuchos::RCP<Epetra_Vector> Disp()
    { return dispnp_; }
  virtual Teuchos::RCP<const Epetra_Map> DofRowMap() const
    { return Teuchos::rcp(discret_->DofRowMap(),false); }
  virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix()
  { return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_); }
  virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() const
  { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_); }

  virtual Teuchos::RCP<DRT::Discretization> Discretization()
    { return discret_; }

  //@}

  /// read restart for given step
  virtual void ReadRestart(int step);

  /// communication object at the interface
  const LINALG::MapExtractor& Interface() const { return interface_; }

  const LINALG::MapExtractor& FreeSurface() const { return freesurface_; }

private:

  /// the interface map setup for interface <-> full translation
  LINALG::MapExtractor interface_;

  LINALG::MapExtractor freesurface_;

  //! @name local copies of input parameters
  Teuchos::RCP<DRT::Discretization> discret_;
  Teuchos::RCP<LINALG::Solver> solver_;
  Teuchos::RCP<ParameterList> params_;
  Teuchos::RCP<IO::DiscretizationWriter> output_;
  //@}

  /// cached preconditioner on constant system matrix
  Teuchos::RCP<LINALG::Preconditioner> precond_;

  //! @name Algorithm core variables
  int step_;
  int numstep_;
  double time_;
  double maxtime_;
  double dt_;

  bool incremental_;

  Teuchos::RCP<LINALG::SparseOperator> sysmat_;

  //! maps for extracting Dirichlet and free DOF sets
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

  //Teuchos::RCP<Epetra_Vector>    zeros_;

  Teuchos::RCP<Epetra_Vector>    dispnp_;
  Teuchos::RCP<Epetra_Vector>    dispn_;
  //Teuchos::RCP<Epetra_Vector>    dispnm_;

  Teuchos::RCP<Epetra_Vector>    residual_;
  //@}

private:

  /// element call
  void EvaluateElements();

  //! @name restart
  int uprestart_;
  //@}

};

}

#endif
#endif

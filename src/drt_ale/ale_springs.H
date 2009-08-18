/*----------------------------------------------------------------------*/
/*!
\file ale_springs.H

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

#ifndef ALE_SPRINGS_H
#define ALE_SPRINGS_H

#include "ale.H"

namespace ALE {


/// springs ALE algorithm
/*!
  Incremental springs based implementation.
 */
class AleSprings : public Ale
{
  friend class AleResultTest;

public:

  AleSprings(Teuchos::RCP<DRT::Discretization> actdis,
             Teuchos::RCP<LINALG::Solver> solver,
             Teuchos::RCP<ParameterList> params,
             Teuchos::RCP<IO::DiscretizationWriter> output,
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

  Teuchos::RCP<LINALG::Preconditioner> ConstPreconditioner() { return Teuchos::null; }

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

  //! @name Algorithm core variables
  int step_;
  int numstep_;
  double time_;
  double maxtime_;
  double dt_;

  Teuchos::RCP<LINALG::SparseOperator> sysmat_;

  //! maps for extracting Dirichlet and free DOF sets
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

  //Teuchos::RCP<Epetra_Vector>    zeros_;

  Teuchos::RCP<Epetra_Vector>    dispnp_;
  Teuchos::RCP<Epetra_Vector>    dispn_;
  //Teuchos::RCP<Epetra_Vector>    dispnm_;
  Teuchos::RCP<Epetra_Vector>    incr_;

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

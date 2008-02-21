/*----------------------------------------------------------------------*/
/*!
\file

\brief Solve FSI problems using a Dirichlet-Neumann partitioning approach

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef FSI_ALE_H
#define FSI_ALE_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_Map.h>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_systemmatrix.H"
#include "../io/io_drt.H"

#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Required.H>

#include <set>

#include "fsi_utils.H"


namespace FSI {


/// General ALE interface
/*!
  An abstract base class meant to be derived for each ALE
  implementation we want to support.

  The purpose of the ALE field is to take an interface displacement and to
  return the deformed interface. All in all pretty simple.

  \warning Further cleanup is still needed.

  \sa Structure, Fluid
  \author u.kue
  \date 11/07
 */
class Ale
{
public:
  virtual ~Ale() {}

  //! @name Apply interface displacements

  /// displacements at the interface are prescribed
  virtual void ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp) = 0;

  //@}

  //! @name Time step helpers

  /// start a new time step
  virtual void PrepareTimeStep() = 0;

  /// calculate linear system
  virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> ddisp) const = 0;

  /// take the current solution to be the final one for this time step
  virtual void Update() = 0;

  /// write output
  virtual void Output() = 0;

  /// read restart for given step
  virtual void ReadRestart(int step) = 0;

  //@}

  //! @name Solver calls

  /// solve this time step with currently applied interface displacements
  virtual void Solve() = 0;

  //@}

  //! @name Extract displacement

  /// the whole displacement field is extracted here
  virtual Teuchos::RCP<Epetra_Vector> ExtractDisplacement() const = 0;

  //@}

  //! @name Vector access

  /// initial guess of Newton's method
  virtual Teuchos::RCP<Epetra_Vector> InitialGuess() = 0;

  /// rhs of Newton's method
  virtual Teuchos::RCP<Epetra_Vector> RHS() const = 0;

  /// unknown velocities (and pressures) at t(n+1)
  virtual Teuchos::RCP<Epetra_Vector> Disp() = 0;

  //@}

  //! @name Misc

  /// dof map of vector of unknowns
  virtual Teuchos::RCP<const Epetra_Map> DofRowMap() const = 0;

  /// direct access to system matrix
  virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() const = 0;

  /// Aii matrix block
  virtual const LINALG::SparseMatrix* InteriorMatrixBlock() const = 0;

  /// Aig matrix block
  virtual const LINALG::SparseMatrix* InterfaceMatrixBlock() const = 0;

  /// direct access to discretization
  virtual Teuchos::RCP<DRT::Discretization> Discretization() = 0;

  //virtual Teuchos::RCP<Epetra_Vector> StructCondRHS() const = 0;

  /// communication object at the interface
  virtual const LINALG::MapExtractor& Interface() const = 0;

  /// build system matrix
  virtual void BuildSystemMatrix(bool full=true) = 0;

  //@}
};


/// linear ALE algorithm
/*!
  The most basic, most simple implementation.
 */
class AleLinear : public Ale
{
  friend class AleResultTest;

public:

  AleLinear(Teuchos::RCP<DRT::Discretization> actdis,
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
  virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> ddisp) const;
  virtual void Solve();
  virtual void Update();
  virtual void Output();
  //@}

  void SetInterfaceMap(Teuchos::RCP<Epetra_Map> im);

  void ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp);

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
  virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() const
  { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_); }

  /// Aii matrix block
  virtual const LINALG::SparseMatrix* InteriorMatrixBlock() const;

  /// Aig matrix block
  virtual const LINALG::SparseMatrix* InterfaceMatrixBlock() const;

  virtual Teuchos::RCP<DRT::Discretization> Discretization()
    { return discret_; }
  //virtual Teuchos::RCP<Epetra_Vector> StructCondRHS() const;

  //@}

  /// read restart for given step
  virtual void ReadRestart(int step);

  /// communication object at the interface
  const LINALG::MapExtractor& Interface() const { return interface_; }

private:

  /// the interface map setup for interface <-> full translation
  LINALG::MapExtractor interface_;

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

  Teuchos::RCP<Epetra_Vector>    dirichtoggle_;
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
  int restartstep_;
  int uprestart_;
  //@}

};
}

#endif
#endif

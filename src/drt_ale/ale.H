/*----------------------------------------------------------------------*/
/*!
\file ale.H

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

#ifndef ALE_H
#define ALE_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_Map.h>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_precond.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/linalg_blocksparsematrix.H"
#include "../drt_io/io.H"

#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Required.H>

#include <set>



namespace ALE {


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

  /// displacements at the free surface are prescribed
  virtual void ApplyFreeSurfaceDisplacements(Teuchos::RCP<Epetra_Vector> fsdisp) = 0;

  //@}

  //! @name Time step helpers

  /// start a new time step
  virtual void PrepareTimeStep() = 0;

  /// calculate linear system
  virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> ddisp) = 0;

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

  /// get the linear solver object used for this field
  virtual Teuchos::RCP<LINALG::Solver> LinearSolver() = 0;

  /// get the linear solver object used for this field
  virtual Teuchos::RCP<LINALG::Preconditioner> ConstPreconditioner() = 0;

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
  virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() = 0;

  /// direct access to system matrix
  virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() const = 0;

  /// direct access to discretization
  virtual Teuchos::RCP<DRT::Discretization> Discretization() = 0;

  //virtual Teuchos::RCP<Epetra_Vector> StructCondRHS() const = 0;

  /// communication object at the interface
  virtual const LINALG::MapExtractor& Interface() const = 0;

  /// communication object at the free surface
  virtual const LINALG::MapExtractor& FreeSurface() const = 0;

  /// build system matrix
  virtual void BuildSystemMatrix(bool full=true) = 0;

  //@}
};


/// Base class of algorithms that use an ale field
class AleBaseAlgorithm
{
public:

  /// constructor
  explicit AleBaseAlgorithm(const Teuchos::ParameterList& prbdyn);

  /// virtual destructor to support polymorph destruction
  virtual ~AleBaseAlgorithm();

  /// ale field solver
  Ale& AleField() { return *ale_; }
  const Ale& AleField() const { return *ale_; }

private:

  /// setup ale algorithm (overriding some parameters with
  /// values specified in given problem-dependent ParameterList)
  void SetupAle(const Teuchos::ParameterList& prbdyn);

  /// ale field solver
  Teuchos::RCP<Ale> ale_;
};

}

#endif
#endif

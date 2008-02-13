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

  \note The fluid and structure implementations do not use such a
  class (yet). This is because we reuse the normal implementations
  here and do not support the exchange for another field
  implementation. Maybe we need this in the future, however.

  The purpose of the ALE field is to take an interface isplacement and
  to return the deformed interface. All in all pretty simple.
 */
class Ale
{
public:
  virtual ~Ale() {}

  /// start a new time step
  virtual void PrepareTimeStep() = 0;

  /// calculate linear system
  virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> ddisp) = 0;

  /// solve this time step with currently applied interface displacements
  virtual void Solve() = 0;

  /// take the current solution to be the final one for this time step
  virtual void Update() = 0;

  /// write output
  virtual void Output() = 0;

  //! @name Interface Map

  /// the map of all interface dofs will be set from the outside
  virtual void SetInterfaceMap(Teuchos::RCP<Epetra_Map> im) = 0;

  //@}

  //! @name Apply interface displacements

  /// displacements at the interface are prescribed
  virtual void ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp) = 0;

  //@}

  //! @name Extract displacement

  /// the whole displacement field is extracted here
  virtual Teuchos::RCP<Epetra_Vector> ExtractDisplacement() = 0;

  //@}

  //! @name access functions for monolithic FSI

  virtual Teuchos::RCP<Epetra_Vector> InitialGuess() = 0;
  virtual Teuchos::RCP<Epetra_Vector> RHS() = 0;
  virtual Teuchos::RCP<Epetra_Vector> Disp() = 0;
  virtual Teuchos::RCP<const Epetra_Map> DofRowMap() = 0;
  virtual Teuchos::RCP<Epetra_Operator> SysMat() = 0;
  virtual Teuchos::RCP<DRT::Discretization> Discretization() = 0;
  virtual Teuchos::RCP<Epetra_Vector> StructCondRHS() = 0;

  //@}

  /// read restart for given step
  virtual void ReadRestart(int step) = 0;

  /// communication object at the interface
  virtual const LINALG::MapExtractor& Interface() const = 0;
};


/// linear ALE algorithm
/*!
  The most basic, most simple implementation.
 */
class AleLinear : public Ale
{
#ifdef RESULTTEST
  friend class AleResultTest;
#endif

public:

  AleLinear(Teuchos::RCP<DRT::Discretization> actdis,
            Teuchos::RCP<LINALG::Solver> solver,
            Teuchos::RCP<ParameterList> params,
            Teuchos::RCP<IO::DiscretizationWriter> output,
            bool dirichletcond=true);

  /// a very simple time loop to be used for standalone ALE problems
  void Integrate();

  //! @name Algorithm core methods
  virtual void PrepareTimeStep();
  virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> ddisp);
  virtual void Solve();
  virtual void Update();
  virtual void Output();
  //@}

  void SetInterfaceMap(Teuchos::RCP<Epetra_Map> im);

  void ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp);

  virtual Teuchos::RCP<Epetra_Vector> ExtractDisplacement();

  //! @name access functions for monolithic FSI

  virtual Teuchos::RCP<Epetra_Vector> InitialGuess()
    { return residual_; }
  virtual Teuchos::RCP<Epetra_Vector> RHS()
    { return residual_; }
  virtual Teuchos::RCP<Epetra_Vector> Disp()
    { return dispnp_; }
  virtual Teuchos::RCP<const Epetra_Map> DofRowMap()
    { return Teuchos::rcp(discret_->DofRowMap(),false); }
  virtual Teuchos::RCP<Epetra_Operator> SysMat()
    { return sysmat_->EpetraOperator(); }
  virtual Teuchos::RCP<DRT::Discretization> Discretization()
    { return discret_; }
  virtual Teuchos::RCP<Epetra_Vector> StructCondRHS();

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

  Teuchos::RCP<LINALG::SparseMatrix> sysmat_;

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

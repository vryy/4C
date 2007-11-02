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
#include "../io/io_drt.H"

#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Required.H>

#include <set>


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

  /// solve this time step with currently applied interface displacements
  virtual void Solve() = 0;

  /// take the current solution to be the final one for this time step
  virtual void Update() = 0;

  /// write output
  virtual void Output() = 0;

  //! @name Interface Map

  /// the map of all interface dofs will be set from the outside
  virtual void SetInterfaceMap(Teuchos::RefCountPtr<Epetra_Map> im) = 0;

  //@}

  //! @name Apply interface displacements

  /// displacements at the interface are prescribed
  virtual void ApplyInterfaceDisplacements(Teuchos::RefCountPtr<Epetra_Vector> idisp) = 0;

  //@}

  //! @name Extract displacement

  /// the whole displacement field is extracted here
  virtual Teuchos::RefCountPtr<Epetra_Vector> ExtractDisplacement() = 0;

  //@}

  /// discretization access
  virtual const DRT::Discretization& Discretization() const = 0;
};


/// linear ALE algorithm
/*!
  The most basic, most simple implementation.
 */
class AleLinear : public Ale,
                  public NOX::Epetra::Interface::Required,
                  public NOX::Epetra::Interface::Jacobian
{
#ifdef RESULTTEST
  friend class AleResultTest;
#endif

public:

  AleLinear(Teuchos::RefCountPtr<DRT::Discretization> actdis,
            Teuchos::RefCountPtr<LINALG::Solver> solver,
            Teuchos::RefCountPtr<ParameterList> params,
            Teuchos::RefCountPtr<IO::DiscretizationWriter> output);

  /// a very simple time loop to be used for standalone ALE problems
  void Integrate();

  //! @name Algorithm core methods
  virtual void PrepareTimeStep();
  virtual void Solve();
  virtual void Update();
  virtual void Output();
  //@}

  //! @name Unused and nonworking
  virtual bool computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);
  virtual bool computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac);
  //@}

  void SetInterfaceMap(Teuchos::RefCountPtr<Epetra_Map> im);

  void ApplyInterfaceDisplacements(Teuchos::RefCountPtr<Epetra_Vector> idisp);

  const DRT::Discretization& Discretization() const { return *discret_; };

  virtual Teuchos::RefCountPtr<Epetra_Vector> ExtractDisplacement();

private:

  //! @name local copies of input parameters
  Teuchos::RefCountPtr<DRT::Discretization> discret_;
  Teuchos::RefCountPtr<LINALG::Solver> solver_;
  Teuchos::RefCountPtr<ParameterList> params_;
  Teuchos::RefCountPtr<IO::DiscretizationWriter> output_;
  //@}

  //! @name Algorithm core variables
  int step_;
  int nstep_;
  double time_;
  double maxtime_;
  double dt_;

  int maxentriesperrow_;

  Teuchos::RefCountPtr<Epetra_CrsMatrix> sysmat_;

  Teuchos::RefCountPtr<Epetra_Vector>    dirichtoggle_;
  //Teuchos::RefCountPtr<Epetra_Vector>    zeros_;

  Teuchos::RefCountPtr<Epetra_Vector>    dispnp_;
  Teuchos::RefCountPtr<Epetra_Vector>    dispn_;
  Teuchos::RefCountPtr<Epetra_Vector>    dispnm_;

  Teuchos::RefCountPtr<Epetra_Vector>    residual_;
  //@}

  /// interface dof map
  Teuchos::RefCountPtr<Epetra_Map> imeshmap_;

  /// communicator between interface dof map and field dof map
  Teuchos::RefCountPtr<Epetra_Import> extractor_;

private:

  /// element call
  void EvaluateElements();

  //! @name Unused and nonworking
  bool haveF_;
  bool haveJacobian_;
  //@}
};
}

#endif
#endif

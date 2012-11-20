/*----------------------------------------------------------------------*/
/*!
\file ale_lin.H

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef ALE_LIN_H
#define ALE_LIN_H

#include "ale.H"

namespace ALE {

/// forward declaration
namespace UTILS
{
  class XFluidFluidMapExtractor;
}


/// linear ALE algorithm
/*!
  The most basic, most simple implementation.
 */
class AleLinear : public Ale
{

public:

  AleLinear(Teuchos::RCP<DRT::Discretization> actdis,
            Teuchos::RCP<LINALG::Solver> solver,
            Teuchos::RCP<Teuchos::ParameterList> params,
            Teuchos::RCP<IO::DiscretizationWriter> output,
            bool incremental,
            bool dirichletcond);

  /// build system matrix
  virtual void BuildSystemMatrix(bool full=true);

  //! @name Algorithm core methods
  virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> ddisp);
  virtual void Solve();
  virtual void Update();
  virtual void Output();
  //@}

  //! @name Access methods

  /// the whole displacement field at time step \f$t^{n+1}\f$ is extracted here
  virtual Teuchos::RCP<Epetra_Vector> ExtractDispnp() const { return dispnp_; }

  /// the whole displacement field at time step \f$t^{n}\f$ is extracted here
  virtual Teuchos::RCP<Epetra_Vector> ExtractDispn() const { return dispn_; }

  //@}

  Teuchos::RCP<LINALG::Preconditioner> ConstPreconditioner() { return precond_; }

  virtual void SolveAleXFluidFluidFSI();

private:

  /// element call
  void EvaluateElements();

  bool incremental_;

  /// cached preconditioner on constant system matrix
  Teuchos::RCP<LINALG::Preconditioner> precond_;

  // toggle vector of xfluidfluid
  Teuchos::RCP<Epetra_Vector>    xfftoggle_;

  /// the interface map setup for xfluid-fluid interface
  Teuchos::RCP<ALE::UTILS::XFluidFluidMapExtractor> xffinterface_;

};

}

#endif

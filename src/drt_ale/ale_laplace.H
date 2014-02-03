/*----------------------------------------------------------------------*/
/*!
\file ale_laplace.H

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef ALE_LAPLACE_H
#define ALE_LAPLACE_H

#include "ale.H"

namespace ALE {


/// linear ALE algorithm
/*!
  The most basic, most simple implementation.
 */
class AleLaplace : public Ale
{

public:

  AleLaplace(Teuchos::RCP<DRT::Discretization> actdis,
            Teuchos::RCP<LINALG::Solver> solver,
            Teuchos::RCP<Teuchos::ParameterList> params,
            Teuchos::RCP<IO::DiscretizationWriter> output,
            bool incremental,
            bool dirichletcond);

  /// build system matrix
  virtual void BuildSystemMatrix(bool full=true);

  //! @name Algorithm core methods
  virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> ddisp, std::string incrementtype = std::string("step"));
  virtual void Solve();
  virtual void Update();
  //@}

  Teuchos::RCP<LINALG::Preconditioner> ConstPreconditioner() { return precond_; }

  //! @name Reading access to displacement
  //@{

  /// get the whole displacement field at time step \f$t^{n+1}\f$
  virtual Teuchos::RCP<const Epetra_Vector> Dispnp() const { return dispnp_; }

  /// get the whole displacement field at time step \f$t^{n}\f$
  virtual Teuchos::RCP<const Epetra_Vector> Dispn() const { return dispn_; }
  //@}

  //! @name Writing access to displacement

  /// write access to whole displacement field at time step \f$t^{n+1}\f$
  virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispnp() const { return dispnp_; }

  //@}

  // Solve call for biofilm growth
  virtual void SolveBioGr();

private:

  /// element call
  void EvaluateElements();

  /// cached preconditioner on constant system matrix
  Teuchos::RCP<LINALG::Preconditioner> precond_;

  bool incremental_;
};

}

#endif

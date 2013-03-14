/*----------------------------------------------------------------------*/
/*!
\file ale_springs_fixed_ref.H

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef ALE_SPRINGS_FIXED_REF_H
#define ALE_SPRINGS_FIXED_REF_H

#include "ale.H"

namespace ALE {


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

public:

  AleSpringsFixedRef(Teuchos::RCP<DRT::Discretization> actdis,
             Teuchos::RCP<LINALG::Solver> solver,
             Teuchos::RCP<Teuchos::ParameterList> params,
             Teuchos::RCP<IO::DiscretizationWriter> output,
             bool incremental,
             bool dirichletcond=true);

  /// build system matrix
  virtual void BuildSystemMatrix(bool full=true);

  //! @name Algorithm core methods
  virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> ddisp);
  virtual void Solve();
  virtual void Update();
  virtual void Output();
  //@}

  Teuchos::RCP<LINALG::Preconditioner> ConstPreconditioner() { return precond_; }

  //! @name Access methods

  /// the whole displacement field at time step \f$t^{n+1}\f$ is extracted here
  virtual Teuchos::RCP<Epetra_Vector> ExtractDispnp() const { return dispnp_; }

  /// the whole displacement field at time step \f$t^{n}\f$ is extracted here
  virtual Teuchos::RCP<Epetra_Vector> ExtractDispn() const { return dispn_; }

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

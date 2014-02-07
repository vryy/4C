/*----------------------------------------------------------------------*/
/*!
\file acou_impl_noli.H

\brief Acoustics Algorithm for implicit time integration for nonlinear wave equation

       Derived classes are implicit Euler, trapezoidal rule, BDF2,
       BDF3, BDF4 and several diagonally implicit RK-schemes

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef ACOU_IMPL_NOLI_H
#define ACOU_IMPL_NOLI_H


#include "acou_impl.H"

namespace ACOU
{

class TimIntImplNoli : public AcouImplicitTimeInt
{
public:
  /// constructor
  TimIntImplNoli(
      const Teuchos::RCP<DRT::DiscretizationHDG>&   actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output
      );

  /// override base class to do the nonlinear iteration loop
  virtual void Solve();

  /// override base class to include recalculation of system matrix
  virtual void UpdateInteriorVariablesAndAssemebleRHS();

  /// update increment
  void IterTraceUpdate();

  /// check convergence and number of iterations
  bool ConvergenceCheck(int itnum);

  // TODO: later on, this should be a base class again, from which all time integration schemes for nonlinear problems are derived
  /// return the name of the time integration scheme
  virtual std::string Name() { return "EulerNoli" ; }

private:

  /// tolerance for nonlinear iteration loop
  double ittol_;

  /// maximum number of nonlinear iterations
  int itmax_;

  /// flag for nonlinear loop, holds if converged or not
  bool stopnonliniter_;

  //! @name norms for convergence check
  double resnorm_;
  double incnorm_;
  //@}

  /// Nonlinear iteration increment vector for trace dofs
  Teuchos::RCP<Epetra_Vector>    veli_;
  /// Nonlinear iteration increment vector for interior dofs
  Teuchos::RCP<Epetra_Vector>    intveli_;
};

} // namespace ACOU

#endif /* ACOU_IMPL_NOLI_H */

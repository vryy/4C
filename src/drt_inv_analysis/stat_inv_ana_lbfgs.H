/*----------------------------------------------------------------------*/
/*!
 * \file stat_inv_ana_lbfgs.H
 *
<pre>
Maintainer: Sebastian Kehl
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef STAT_INV_ANA_LBFGS_H
#define STAT_INV_ANA_LBFGS_H

#include "stat_inv_analysis.H"

namespace DRT
{
  namespace UTILS
  {
    template <typename> class TimIntMStep;
  }
}

namespace STR
{
namespace INVANA
{

class StatInvAnaLBFGS : public StatInvAnalysis
{

private:

  //! Armijo rule based line search with polynomial models:
  int EvaluateArmijoRule(double* tauopt, int* numsteps);

  //! quadratic model:
  int polymod(double e_o, double dfp, double tau_n, double e_n, double blow, double bhigh, double* tauopt);

  //! cubic model:
  int polymod(double e_o, double dfp, double tau_n, double e_n, double blow, double bhigh, double tau_l, double e_l, double* tauopt);

  //! final print out of the results
  void Summarize();

  //! Print step
  void PrintOptStep(double stepsize, int numsteps);

  //! store new vector quantities
  void StoreVectors();

  //! compute new direction based on lbfgs formula
  void ComputeDirection();

  //! step size for the gradient scheme
  double stepsize_;

  //! maximum number of iterations
  int maxiter_;

  //! iteration counter
  int runc_;

  //! tolerance for convergence criterion
  double convtol_;

  //! current convergence criterion value
  double convcritc_;

  //! maximum size of storage (number of multivectors)
  int ssize_;

  //! current size of the storage
  int actsize_;

  //! storage of solution delta
  Teuchos::RCP<DRT::UTILS::TimIntMStep<Epetra_Vector> > sstore_;

  //! storage of gradient delta
  Teuchos::RCP<DRT::UTILS::TimIntMStep<Epetra_Vector> > ystore_;

  //! current direction of steepest descent maybe properly scaled by hessian approximation
  Teuchos::RCP<Epetra_MultiVector> p_;

  //! current step made based on direction p_ using sacled with a proper stepsize
  Teuchos::RCP<Epetra_MultiVector> step_;

public:

  //! standard constructor
  StatInvAnaLBFGS(Teuchos::RCP<DRT::Discretization> dis);

  ~StatInvAnaLBFGS() {};

  //! start the optimization routine
  void Optimize();

}; // class StatInvAnaLBFGS
}  // namespace INVANA
}  // namespace STR

#endif /*INV_ANALYSIS_LBFGSH*/

/*----------------------------------------------------------------------*/
/*!
 * \file stat_inv_ana_graddesc.H
 *
<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef STAT_INV_ANA_GRADDESC_H
#define STAT_INV_ANA_GRADDESC_H

#include "stat_inv_analysis.H"

namespace STR
{
namespace INVANA
{


  /*!
  \brief Steepest descent with armijo rule line search

  This is a standard steepest gradient descent implementation derived from StatInvAnalysis.

  Stepsize is computed such that the sufficient decrease condition is fulfilled.

  The LineSearch is based on poynomial approximations based on the available information
  in each line search step, i.e a quadratic model is used in the first step, while for
  steps >1 cubic models are applied.

  */

class StatInvAnaGradDesc : public StatInvAnalysis
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

  //! current direction of steepest descent
  Teuchos::RCP<Epetra_MultiVector> p_;

  //! current step made based on direction p_ sacled with a proper stepsize
  Teuchos::RCP<Epetra_MultiVector> step_;

public:

  //! standard constructor
  StatInvAnaGradDesc(Teuchos::RCP<DRT::Discretization> dis);

  ~StatInvAnaGradDesc() {};

  //! start the optimization routine
  void Optimize();

}; // class StatInvAnaGradDesc
}  // namespace INVANA
}  // namespace STR

#endif /*INV_ANALYSIS_GRADDESC_H*/

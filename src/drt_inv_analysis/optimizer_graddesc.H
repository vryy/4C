/*----------------------------------------------------------------------*/
/*!
 * \file optimizer_graddesc.H
 * \brief gradient descent optimization algorithm
 *
<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef OPTIMIZER_GRADDESC_H
#define OPTIMIZER_GRADDESC_H

#include "optimizer_base.H"

namespace STR
{
namespace INVANA
{


  /*!
  \brief Steepest gradient descent with armijo rule line search

  This is a standard steepest gradient descent implementation derived from OptimizerBase.

  Stepsize is computed such that the sufficient decrease condition is fulfilled.

  The LineSearch is based on poynomial approximations based on the available information
  in each line search step, i.e a quadratic model is used in the first step, while for
  steps >1 cubic models are applied.

  \author kehl
  */

class OptimizerGradDesc : public OptimizerBase
{
public:

  //! constructor
  OptimizerGradDesc(const Teuchos::ParameterList& invp);

  //! destructor
  virtual ~OptimizerGradDesc() {};

  //! use some specific optimization routine
  virtual void Integrate();

  //! read restart information
  virtual void ReadRestart(int run);

  //! write restart information
  virtual void WriteRestart();

  //! setup algorithm specific stuff
  virtual void Setup();

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

  //! Print line search step
  void PrintLSStep(double stepsize, int numsteps);

  //! current direction of steepest descent
  Teuchos::RCP<Epetra_MultiVector> p_;

  //! current step based on direction p_ scaled with a proper stepsize
  Teuchos::RCP<Epetra_MultiVector> step_;

};
}  // namespace INVANA
}  // namespace STR

#endif /*OPTIMIZER_GRADDESC_H*/

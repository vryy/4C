/*----------------------------------------------------------------------*/
/*!
 * \file optimizer_lbfgs.H
 * \brief LBFGS optimization algorithm
 *
<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef OPTIMIZER_LBFGS_H
#define OPTIMIZER_LBFGS_H

#include "optimizer_base.H"

namespace TIMINT
{
  template <typename> class TimIntMStep;
}

namespace INVANA
{

/*! \class OptimizerLBFGS
 * \brief Limited memory bfgs algorithm with armijo rule line search
 *
 * Implementation of the Limited Memory BFGS algorithm as e.g. given in:\n
 * Nocedal, Jorge. "Updating quasi-Newton matrices with limited storage." Mathematics of computation 35.151 (1980): 773-782.
 *
 * With a linesearch implemented from:\n
 * Kelley, Carl T. Iterative methods for optimization. Vol. 18. Siam, 1999.
 *
 * \author kehl
 */
class OptimizerLBFGS : public OptimizerBase
{
public:

  //! constructor
  OptimizerLBFGS(const Teuchos::ParameterList& invp);

  //! destructor
  virtual ~OptimizerLBFGS() {};

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

  //! Print optimization step
  void PrintOptStep(double stepsize, int numsteps);

  //! Print line search step
  void PrintLSStep(double stepsize, int numsteps);

  //! store new vector quantities
  void StoreVectors();

  //! compute new direction based on lbfgs formula
  void ComputeDirection();

  //! @name Ouput routines
  //@{

  //! write binary output at the end
  void WriteOutput();

  //@}

  //! initial scaling for the LBFGS two loop recursion algorithm
  bool initscal_;

  //! maximum size of storage (number of multivectors)
  int ssize_;

  //! current size of the storage
  int actsize_;

  //! storage of solution delta
  Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > sstore_;

  //! storage of gradient delta
  Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > ystore_;

  //! current direction of steepest descent maybe properly scaled by hessian approximation
  Teuchos::RCP<Epetra_MultiVector> p_;

  //! current step based on direction p_ scaled with a proper stepsize
  Teuchos::RCP<Epetra_MultiVector> step_;

};
}  // namespace INVANA

#endif /*OPTIMIZER_LBFGS_H*/

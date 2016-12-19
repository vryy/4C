/*----------------------------------------------------------------------*/
/*!
 * \file likelihood_evaluation.H
 * \brief Evaluation of mixture likelihoods
 *
<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/

#if __cplusplus >= 201103L

#ifndef LIKELIHOOD_EVALUATION_H_
#define LIKELIHOOD_EVALUATION_H_

#include <random>

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

namespace Teuchos
{
class ParameterList;
}

namespace INVANA
{

class InvanaBase;
class CholFactorBase;


/*! \class LogLikePrior
 * \brief LogLikelihood prior evaluator
 *
 * This class performs the evaluation of a gaussian prior
 * given in terms of mean and covariance.
 *
 * The covariance is given by the INVANA::Ifpack_ICT preconditioner
 * which keeps a (possibly incomplete) factorization of the covariance
 * matrix. It can thus be efficiently used to evaluate the log-
 * likelihood as well as to draw samples from the distribution.
 *
 * \author kehl \date 06/2016
 */
class LogLikePrior
{
public:

  //! @name Construction, destruction, setup
  //@{

  //! constructor
  LogLikePrior();

  //! destructor
  virtual ~LogLikePrior() {}

  //! initialize variables
  void Init(Teuchos::RCP<CholFactorBase> cov_factor, ///[in] factor of covariance matrix
      Teuchos::RCP<Epetra_Vector> pmean, ///[in] prior mean
      double covscale ///[in] scale facrtor for the covariance
      );
  //@}

  //! @name Evaluation
  //@{

  //! Evaluate the likelihood function
  void Evaluate(const Epetra_Vector& state, double* val);

  //! Draw a sample from this prior
  void DrawfromPrior(Epetra_Vector& draw);

  //! Draw a sample with this covariance (but scaled) and at a different mean
  void DrawProposal(const Epetra_Vector& mean, /// [in] mean
      const double& fac, /// [in] scaling factor for covariance
      Epetra_Vector& draw /// [out] proposal
      );
  //@}

  //! @name Access
  //@{

  //! the covariance factorization
  const CholFactorBase& Covariance() const {return *cov_factor_;}

  //! the mean
  const Epetra_Vector& Mean() const {return *pmean_;}

  //@}

private:

  //! @name prior evaluation helpers
  //@{

  /*! \brief cholesky factorization of the prior covariance
   *
   * \author kehl \date 06/2016
   */
  Teuchos::RCP<CholFactorBase> cov_factor_;

  //! prior mean
  Teuchos::RCP<Epetra_Vector> pmean_;

  //! covariance scale factor
  double cov_scale_;
  //@}

  //! @name random number generation
  //@{

  // random number generator
  std::mt19937 generator_;

  //! zero mean, variance 1, normal distribution
  std::normal_distribution<double> distribution_;
  //@}
};


/*! \class LogLikeMixture
 * \brief Mixture of Log-Likelihoods
 *
 * This class allows the evaluation of the mixture likelihood in the
 * SMC setting given by a posterior and a prior probablilty. How these
 * probabilities mix is not ofinterest for this class. It just evaluates
 * both of them by using INVANA::LogLikePrior for the prior and
 * INVANA::InvanaBase for the posterior.
 *
 * Since it is the "interface" to the prior and posterior probabilities of
 * particles in the INVANA::ParticleGroup it additionally gives access to
 * draws from the INVANA::LogLikePrior which is needed in the SMC-context.
 *
 * \author kehl \date 06/2016
 */
class LogLikeMixture
{
public:

  //! @name Construction, destruction, setup
  //@{

  //! constructor
  LogLikeMixture(const Teuchos::ParameterList& invp);

  //! destructor
  virtual ~LogLikeMixture() {}

  //! initialize variables
  void Init(Teuchos::RCP<INVANA::InvanaBase> posterior, ///[in] to evalaute the posterior
      Teuchos::RCP<LogLikePrior> prior ///[in] to evaluate the prior
      );

  //@}

  //! @name Evaluation methods
  //@{

  //! compute results due to state
  int EvaluateMixture(const Epetra_Vector& state, ///[in] new state to be processed
      double& val_posterior, /// [out] posterior value
      double& val_prior /// [out] prior value
      );

  //! compute results due to state and return also the primal solution
  int EvaluateMixture(const Epetra_Vector& state, ///[in] new state to be processed
      double& val_posterior, /// [out] posterior value
      double& val_prior, /// [out] prior value
      Epetra_Vector& solution /// [out] pde solution
      );

  //! draw from the prior distribution
  void DrawfromPrior(Epetra_Vector& draw);

  //! Draw a sample with this covariance (but scaled) and at a different mean
  void DrawProposal(const Epetra_Vector& mean, /// [in] mean
      const double& fac, /// [in] scaling factor for covariance
      Epetra_Vector& draw /// [out] proposal
      );
  //@}

  //! @name Query
  //@{

  //! Vector layout
  const Epetra_Map& StateMap();

  //! the posterior evaluator
  const INVANA::InvanaBase& EvalPost() const {return *eval_post_;}

  //! the prior evaluator
  const INVANA::LogLikePrior& EvalPrior() const {return *eval_prior_;}

  //@}

private:

  //! @name Likelihood evaluators
  //@{

  /*! \brief Posterior likelihood
   *
   * \author kehl \date 06/2016
   */
  Teuchos::RCP<INVANA::InvanaBase> eval_post_;

  /*! \brief Prior likelihood
   *
   * \author kehl \date 06/2016
   */
  Teuchos::RCP<LogLikePrior> eval_prior_;

  //@}

  //! @name posterior scaling
  //@{

  //! is the posterior scaled somehow
  bool prob_scaling_;

  //! scaling factor of the posterior in case
  double prob_scale_fac_;

  //@}
};
}



#endif /* LIKELIHOOD_EVALUATION_H_ */
#endif

/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of mixture likelihoods

\level 3

*/
/*----------------------------------------------------------------------*/

#if __cplusplus >= 201103L

#include "likelihood_evaluation.H"

#include "invana_base.H"  // for the posterior evaluation
#include "objective_funct.H"
#include "chol_factor_base.H"  // for the prior evaluation
#include "Epetra_CrsMatrix.h"
#include "../drt_inpar/inpar_invanalysis.H"

#include "../drt_lib/drt_dserror.H"
#include <random>
#include <chrono>

/*----------------------------------------------------------------------*/
INVANA::LogLikeMixture::LogLikeMixture(const Teuchos::ParameterList& invp)
    : eval_post_(Teuchos::null),
      eval_prior_(Teuchos::null),
      prob_scaling_(DRT::INPUT::IntegralValue<bool>(invp, "OBJECTIVEFUNCTSCAL")),
      prob_scale_fac_(1.0)
{
}

/*----------------------------------------------------------------------*/
void INVANA::LogLikeMixture::Init(
    Teuchos::RCP<INVANA::InvanaBase> posterior, Teuchos::RCP<INVANA::LogLikePrior> prior)
{
  eval_post_ = posterior;
  eval_prior_ = prior;

  if (prob_scaling_) prob_scale_fac_ = eval_post_->ObjectiveFunct()->GetScaleFac();

  return;
}

/*----------------------------------------------------------------------*/
int INVANA::LogLikeMixture::EvaluateMixture(
    const Epetra_Vector& state, double& val_posterior, double& val_prior, Epetra_Vector& solution)
{
  int err = EvaluateMixture(state, val_posterior, val_prior);

  int numvectors = eval_post_->GetPrimalVariables().NumVectors();
  solution.Update(1.0, *eval_post_->GetPrimalVariables()(numvectors - 1), 0.0);

  return err;
}

/*----------------------------------------------------------------------*/
int INVANA::LogLikeMixture::EvaluateMixture(
    const Epetra_Vector& state, double& val_posterior, double& val_prior)
{
  // evaluate the posterior
  double posterior = 0.0;
  int err = eval_post_->Evaluate(state, &posterior, Teuchos::null);

  // since InvanaBase was made for minimization ->
  val_posterior = posterior * (-1.0);

  // and since InvanaBase might have scaled the problem
  if (prob_scaling_) val_posterior /= prob_scale_fac_;

  // evaluation of the prior
  double prior = 0.0;
  eval_prior_->Evaluate(state, &prior);

  val_prior = prior;

  return err;
}

/*----------------------------------------------------------------------*/
const Epetra_Map& INVANA::LogLikeMixture::StateMap() { return *(eval_post_->VectorRowLayout()); }

/*----------------------------------------------------------------------*/
void INVANA::LogLikeMixture::DrawfromPrior(Epetra_Vector& draw)
{
  eval_prior_->DrawfromPrior(draw);

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::LogLikeMixture::DrawProposal(
    const Epetra_Vector& mean, const double& fac, Epetra_Vector& draw)
{
  eval_prior_->DrawProposal(mean, fac, draw);
  return;
}

/*----------------------------------------------------------------------*/
INVANA::LogLikePrior::LogLikePrior() : cov_scale_(1.0), generator_(), distribution_(0.0, 1.0) {}

/*----------------------------------------------------------------------*/
void INVANA::LogLikePrior::Init(
    Teuchos::RCP<CholFactorBase> cov_factor, Teuchos::RCP<Epetra_Vector> pmean, double covscale)
{
  cov_factor_ = cov_factor;
  pmean_ = pmean;
  cov_scale_ = covscale;

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  generator_.seed(seed);

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::LogLikePrior::DrawfromPrior(Epetra_Vector& draw)
{
  // the lower triangular factor
  const Epetra_CrsMatrix& L = cov_factor_->H();

  // populate a vector with uncorrelated normally distributed noise
  Epetra_Vector sample(L.RowMap());

  if (not L.RowMap().SameAs(draw.Map()))
    dserror("Pass in a vector matching the row layout of the covariance matrix!");

  double* values;
  sample.ExtractView(&values);
  for (int i = 0; i < sample.MyLength(); i++) values[i] = distribution_(generator_);

  // Correlate
  L.Multiply(false, sample, draw);

  // mean shift and scale
  double sqrtscale = sqrt(cov_scale_);
  draw.Update(1.0, *pmean_, sqrtscale);
  return;
}

/*----------------------------------------------------------------------*/
void INVANA::LogLikePrior::DrawProposal(
    const Epetra_Vector& mean, const double& fac, Epetra_Vector& draw)
{
  // fac should be a scale factor for the covariance
  // so for the factor it must be squarerooted
  double sqrtfac = sqrt(fac);

  // the lower triangular factor
  const Epetra_CrsMatrix& L = cov_factor_->H();

  // populate a vector with uncorrelated normally distributed noise
  Epetra_Vector sample(L.RowMap());

  if (not L.RowMap().SameAs(draw.Map()))
    dserror("Pass in a vector matching the row layout of the covariance matrix!");

  double* values;
  sample.ExtractView(&values);
  for (int i = 0; i < sample.MyLength(); i++) values[i] = distribution_(generator_);

  // Correlate
  L.Multiply(false, sample, draw);

  // mean shift
  draw.Update(1.0, mean, sqrtfac);
  return;
}

/*----------------------------------------------------------------------*/
void INVANA::LogLikePrior::Evaluate(const Epetra_Vector& state, double* val)
{
  Epetra_Vector tmp(state);
  tmp.Update(-1.0, *pmean_, 1.0);

  // Cov^-1*(state-mean)
  Epetra_Vector solv(tmp.Map(), false);
  cov_factor_->ApplyInverse(tmp, solv);

  // (state-mean)^T * Cov^-1 * (state-mean)
  double prior = 0.0;
  tmp.Dot(solv, &prior);

  *val = prior * (-1.0) / cov_scale_;

  return;
}

#endif

/*----------------------------------------------------------------------*/
/*!
\file regularization_tikhonov.cpp

\brief Tikhonov type regularization

<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>

!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "regularization_tikhonov.H"

#include "invana_utils.H"
#include "initial_guess.H"
#include "Teuchos_ParameterList.hpp"

/*----------------------------------------------------------------------*/
INVANA::RegularizationTikhonov::RegularizationTikhonov() :
  RegularizationBase()
{}

/*----------------------------------------------------------------------*/
void INVANA::RegularizationTikhonov::Setup(const Teuchos::ParameterList& invp)
{
  params_ = invp;
  cov_factor_ = INVANA::CreateICT_lowmem(initguess_->Covariance(),params_);

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::RegularizationTikhonov::Evaluate(const Epetra_MultiVector& theta, double* value)
{
  double val = 0.0;

  Epetra_MultiVector factorr(*initguess_->Mean());
  Epetra_MultiVector factorl(factorr.Map(),1,true);

  factorr.Update(1.0,theta,-1.0);
  cov_factor_->ApplyInverse(factorr,factorl);

  factorr.Dot(factorl,&val);
  *value += 0.5 * weight_ * val;

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::RegularizationTikhonov::EvaluateGradient(const Epetra_MultiVector& theta,
    Teuchos::RCP<Epetra_MultiVector> gradient)
{
  // the gradient of the regularization functional only
  Epetra_MultiVector tmp(*gradient);

  Epetra_MultiVector factorr(*initguess_->Mean());
  factorr.Update(1.0,theta,-1.0);

  cov_factor_->ApplyInverse(factorr,tmp);

  gradient->Update(weight_,tmp,1.0);

  return;
}


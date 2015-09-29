/*----------------------------------------------------------------------*/
/*!
\file regularization_tvdtikh.cpp

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            089 - 289-15271
</pre>

!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "regularization_tvdtikh.H"

/*----------------------------------------------------------------------*/
/* standard constructor                                     keh 10/14   */
/*----------------------------------------------------------------------*/
INVANA::RegularizationTotalVariationTikhonov::RegularizationTotalVariationTikhonov() :
  RegularizationBase()
{
  regularization_tvd_  = Teuchos::rcp(new INVANA::RegularizationTotalVariation());
  regularization_tikh_ = Teuchos::rcp(new INVANA::RegularizationTikhonov());
}

void INVANA::RegularizationTotalVariationTikhonov::Init(Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<ConnectivityData> connectivity)
{
  regularization_tvd_->Init(discret,connectivity);
  regularization_tikh_->Init(discret,connectivity);
  return;
}

/*----------------------------------------------------------------------*/
/* Setup                                                    keh 10/14   */
/*----------------------------------------------------------------------*/
void INVANA::RegularizationTotalVariationTikhonov::Setup(const Teuchos::ParameterList& invp)
{
  regularization_tvd_->Setup(invp);
  regularization_tikh_->Setup(invp);
  return;
}

/*----------------------------------------------------------------------*/
/* Evaluate                                                 keh 10/14   */
/*----------------------------------------------------------------------*/
void INVANA::RegularizationTotalVariationTikhonov::Evaluate(const Epetra_MultiVector& theta, double* value)
{
  regularization_tvd_->Evaluate(theta,value);
  regularization_tikh_->Evaluate(theta,value);
  return;
}

void INVANA::RegularizationTotalVariationTikhonov::EvaluateGradient(const Epetra_MultiVector& theta, Teuchos::RCP<Epetra_MultiVector> gradient)
{
  regularization_tvd_->EvaluateGradient(theta,gradient);
  regularization_tikh_->EvaluateGradient(theta,gradient);
  return;
}


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
#include "../linalg/linalg_mapextractor.H"
#include "matpar_manager.H" //for ConnectivityData
#include "../drt_lib/drt_dserror.H"

#include "Teuchos_ParameterList.hpp"

/*----------------------------------------------------------------------*/
INVANA::RegularizationTikhonov::RegularizationTikhonov() :
  RegularizationBase(),
mean_(0.0)
{}

/*----------------------------------------------------------------------*/
void INVANA::RegularizationTikhonov::Setup(const Teuchos::ParameterList& invp)
{
  weight_ = invp.get<double>("REG_WEIGHT");
  mean_ = invp.get<double>("REG_MEAN");
  return;
}

/*----------------------------------------------------------------------*/
void INVANA::RegularizationTikhonov::Evaluate(const Epetra_MultiVector& theta, double* value)
{
  double val = 0.0;

  Epetra_MultiVector meanvaluevector(*connectivity_->MapExtractor()->FullMap(),theta.NumVectors());
  meanvaluevector.PutScalar(mean_);
  meanvaluevector.Update(1.0,theta,-1.0);

  INVANA::MVNorm(meanvaluevector, *(connectivity_->MapExtractor()->FullMap()), 2, &val);
  *value += 0.5 * weight_ * val * val;

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::RegularizationTikhonov::EvaluateGradient(const Epetra_MultiVector& theta,
    Teuchos::RCP<Epetra_MultiVector> gradient)
{
  Epetra_MultiVector meanvaluevector(*connectivity_->MapExtractor()->FullMap(),theta.NumVectors());
  meanvaluevector.PutScalar(mean_);
  meanvaluevector.Update(1.0,theta,-1.0);

  gradient->Update(weight_, meanvaluevector, 1.0);

  return;
}


/*----------------------------------------------------------------------*/
/*!
 \file regularization_tikhonov.cpp

 <pre>
 Maintainer: Sebastian Kehl
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
/* constructor */
INVANA::RegularizationTikhonov::RegularizationTikhonov() :
  RegularizationBase()
{}

void INVANA::RegularizationTikhonov::Setup(const Teuchos::ParameterList& invp)
{
  weight_ = invp.get<double>("REG_WEIGHT_TIKH");
  meanvalue_ = invp.get<double>("MEAN_VAL_TIKH");
  return;
}

void INVANA::RegularizationTikhonov::Evaluate(const Epetra_MultiVector& theta, double* value)
{
  double val = 0.0;

  Epetra_MultiVector meanvaluevector(*connectivity_->MapExtractor()->FullMap(),theta.NumVectors());
  meanvaluevector.PutScalar(meanvalue_);
  meanvaluevector.Update(1.0,theta,-1.0);

  INVANA::MVNorm(meanvaluevector, *(connectivity_->MapExtractor()->FullMap()), 2, &val);
  *value += 0.5 * weight_ * val * val;

  return;
}

void INVANA::RegularizationTikhonov::EvaluateGradient(const Epetra_MultiVector& theta, Teuchos::RCP<Epetra_MultiVector> gradient)
{
  Epetra_MultiVector meanvaluevector(*connectivity_->MapExtractor()->FullMap(),theta.NumVectors());
  meanvaluevector.PutScalar(meanvalue_);
  meanvaluevector.Update(1.0,theta,-1.0);

  gradient->Update(weight_, meanvaluevector, 1.0);

  return;
}


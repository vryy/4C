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
STR::INVANA::RegularizationTikhonov::RegularizationTikhonov( const Teuchos::ParameterList& invp) :
  RegularizationBase(invp)
{
  return;
}

void STR::INVANA::RegularizationTikhonov::Evaluate(const Epetra_MultiVector& theta, double* value)
{
  double val = 0.0;
  STR::INVANA::MVNorm(theta, *(connectivity_->MapExtractor()->FullMap()), 2, &val);
  *value += 0.5 * weight_ * val * val;

  return;
}

void STR::INVANA::RegularizationTikhonov::EvaluateGradient(const Epetra_MultiVector& theta, Teuchos::RCP<Epetra_MultiVector> gradient)
{
  gradient->Update(weight_, theta, 1.0);

}


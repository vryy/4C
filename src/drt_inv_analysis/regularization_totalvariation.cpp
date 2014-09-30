/*----------------------------------------------------------------------*/
/*!
\file regularization_totalvariation.cpp

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>

!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "regularization_totalvariation.H"

#include "invana_utils.H"
#include "../drt_lib/drt_dserror.H"

#include "Teuchos_ParameterList.hpp"


/*----------------------------------------------------------------------*/
/* constructor */
STR::INVANA::RegularizationTotalVariation::RegularizationTotalVariation(const Teuchos::ParameterList& invp) :
  RegularizationBase(invp)
{
  return;
}

void STR::INVANA::RegularizationTotalVariation::Setup()
{

  dserror("must be implemented first");
  return;
}

void STR::INVANA::RegularizationTotalVariation::Evaluate(const Epetra_MultiVector& theta, double* value)
{
  return;
}

void STR::INVANA::RegularizationTotalVariation::EvaluateGradient(const Epetra_MultiVector& theta, Teuchos::RCP<Epetra_MultiVector> gradient)
{
  return;
}


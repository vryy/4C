/*----------------------------------------------------------------------*/
/*!
\file regularization_base.H

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>

!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "regularization_base.H"

#include "../linalg/linalg_mapextractor.H"
#include "Teuchos_ParameterList.hpp"


/*----------------------------------------------------------------------*/
/* constructor */
INVANA::RegularizationBase::RegularizationBase(const Teuchos::ParameterList& invp) :
discret_(Teuchos::null),
connectivity_(Teuchos::null),
weight_(0.0)
{
  weight_ = invp.get<double>("REG_WEIGHT");

  return;
}

void INVANA::RegularizationBase::Init(Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<ConnectivityData> connectivity)
{
  discret_=discret;
  connectivity_=connectivity;
}

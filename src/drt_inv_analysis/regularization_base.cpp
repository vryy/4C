/*----------------------------------------------------------------------*/
/*!
\file regularization_base.cpp

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
INVANA::RegularizationBase::RegularizationBase() :
discret_(Teuchos::null),
connectivity_(Teuchos::null),
weight_(0.0)
{}

void INVANA::RegularizationBase::Init(Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<ConnectivityData> connectivity)
{
  discret_=discret;
  connectivity_=connectivity;
  return;
}

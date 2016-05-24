/*----------------------------------------------------------------------*/
/*!
\file regularization_base.cpp

\brief Base class for regularization of optimization problems

<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "regularization_base.H"

#include "../linalg/linalg_mapextractor.H"
#include "Teuchos_ParameterList.hpp"


/*----------------------------------------------------------------------*/
/* constructor */
INVANA::RegularizationBase::RegularizationBase() :
connectivity_(Teuchos::null),
weight_(0.0)
{}

void INVANA::RegularizationBase::Init(Teuchos::RCP<ConnectivityData> connectivity)
{
  connectivity_=connectivity;
  return;
}

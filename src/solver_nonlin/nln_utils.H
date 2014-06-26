/*----------------------------------------------------------------------*/
/*!
\file nln_utils.H

\brief Collection of utilities for nonlinear solver objects

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef NLN_UTILS_H_
#define NLN_UTILS_H_

/*----------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------*/
/* forward declarations */

namespace Teuchos
{
  class ParameterList;
}

/*----------------------------------------------------------------------*/
/* definition of methods */
namespace NLNSOL
{
  namespace UTILS
  {
    Teuchos::RCP<Teuchos::ParameterList> CreateParamListFromXML();
  }
}

#endif /* NLN_UTILS_H_ */

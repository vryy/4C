/*!------------------------------------------------------------------------------------------------*
\file topopt_utils.H

\brief collection of functions in namespace topopt

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef TOPOPT_UTILS_H_
#define TOPOPT_UTILS_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

using namespace std;
using namespace Teuchos;

namespace DRT
{
  class Element;
}

namespace TOPOPT
{
  // prints the logo of topology optimization
  void printTopOptLogo();
}

#endif /* TOPOPT_UTILS_H_ */
#endif  // #ifdef CCADISCRET

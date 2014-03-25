/*----------------------------------------------------------------------*/
/*!
\file inpar_fpsi.H

<pre>
Maintainer: Andreas Rauch
            rauch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef INPAR_IMMERSED_H_
#define INPAR_IMMERSED_H_

#include "inpar_parameterlist_utils.H"

namespace INPAR
{

namespace IMMERSED
{

/*----------------------------------------------------------------------*
 | Coupling Methods                                                     |
 *----------------------------------------------------------------------*/
enum _IMMERSED_COUPLING
{
  partitioned,
  monolithic
};

enum _IMMERSED_COUPLING_SCHEME
{
  neumannneumann
};

enum _IMMERSED_PROJECTION
{
  shapefunctions,
  mortar
};

}// namespace IMMERSED
}// namespace INPAR
/*----------------------------------------------------------------------*/
#endif /* INPAR_IMMERSED_H_ */

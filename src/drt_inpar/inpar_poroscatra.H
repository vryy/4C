/*!------------------------------------------------------------------------------------------------*
\file inpar_poroscatra.H

\brief

<pre>
Maintainer: Anh-Tu Vuong
            vuong@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef INPAR_POROSCATRA_H_
#define INPAR_POROSCATRA_H_

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace INPAR
{

/// Input parameters for poroelasticity
namespace PORO_SCATRA
{

  /// Type of coupling strategy for poro scatra problems
  enum SolutionSchemeOverFields
  {
    Part_ScatraToPoro,
    Part_PoroToScatra,
    Part_TwoWay
 //   Monolithic
  };


}  // namespace PORO_SCATRA

}  // namespace INPAR

/*----------------------------------------------------------------------*/



#endif /* INPAR_POROSCATRA_H_ */

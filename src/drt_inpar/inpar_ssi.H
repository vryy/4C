/*!------------------------------------------------------------------------------------------------*
\file inpar_ssi.H

\brief

<pre>
Maintainer: Anh-Tu Vuong
            vuong@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef INPAR_SSI_H_
#define INPAR_SSI_H_

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace INPAR
{

/// Input parameters for poroelasticity
namespace SSI
{

  /// Type of coupling strategy for SSI problems
  enum SolutionSchemeOverFields
  {
    Part_ScatraToSolid,
    Part_SolidToScatra,
    Part_TwoWay
 //   SequStagg,
 //   IterStagg,
 //   Monolithic,
  };


}  // namespace SSI

}  // namespace INPAR

#endif /* INPAR_SSI_H_ */

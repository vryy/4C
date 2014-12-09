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
    ssi_OneWay_ScatraToSolid,
    ssi_OneWay_SolidToScatra,
    //ssi_SequStagg_ScatraToSolid,
    //ssi_SequStagg_SolidToScatra,
    ssi_IterStagg,
    ssi_IterStaggFixedRel_ScatraToSolid,
    ssi_IterStaggFixedRel_SolidToScatra,
    ssi_IterStaggAitken_ScatraToSolid,
    ssi_IterStaggAitken_SolidToScatra
    //IterStaggAitkenIrons,
    //Monolithic
  };


}  // namespace SSI

}  // namespace INPAR

#endif /* INPAR_SSI_H_ */

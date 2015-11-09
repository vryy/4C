/*!------------------------------------------------------------------------------------------------*
\file inpar_ssi.H

\brief

<pre>
Maintainer: Anh-Tu Vuong
            vuong@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef INPAR_SSI_H_
#define INPAR_SSI_H_

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT { namespace INPUT { class ConditionDefinition; } }

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace INPAR
{

/// Input parameters for solid-scatra-interaction
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

  /// Type of coupling strategy between the two fields of the SSI problems
  enum FieldCoupling
  {
    coupling_match,
    coupling_meshtying,
    coupling_volmortar
  };


  /// set the ssi parameters
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  /// set specific ssi conditions
  void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

}  // namespace SSI

}  // namespace INPAR

#endif /* INPAR_SSI_H_ */

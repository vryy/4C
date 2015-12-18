/*--------------------------------------------------------------------------*/
/*!
\file inpar_ehl.H

\brief Elastohydrodynamic lubrication (lubrication structure interaction) parameters

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_INPAR_INPAR_EHL_H_
#define SRC_DRT_INPAR_INPAR_EHL_H_


#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT { namespace INPUT { class ConditionDefinition; } }

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace INPAR
{

/// Input parameters for lubrication structure interaction
namespace EHL
{

  /// Type of coupling strategy for EHL problems
  enum SolutionSchemeOverFields
  {
    ehl_IterStagg,
    ehl_IterStaggFixedRel_LubToStr,
    ehl_IterStaggFixedRel_StrToLub,
    ehl_IterStaggAitken_LubToStr,
    ehl_IterStaggAitken_StrToLub
  };

  /// Type of coupling strategy between the two fields of the EHL problems
  enum FieldCoupling
  {
    coupling_none,
    coupling_matching
  };


  /// set the ehl parameters
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);


}  // namespace EHL

}  // namespace INPAR


#endif /* SRC_DRT_INPAR_INPAR_EHL_H_ */

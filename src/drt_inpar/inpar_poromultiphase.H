/*----------------------------------------------------------------------*/
/*!
 \file inpar_poromultiphase.H

 \brief input paramters for porous multiphase problem

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_INPAR_INPAR_POROMULTIPHASE_H_
#define SRC_DRT_INPAR_INPAR_POROMULTIPHASE_H_


#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT { namespace INPUT { class ConditionDefinition; } }

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace INPAR
{

/// Input parameters for solid-scatra-interaction
namespace POROMULTIPHASE
{

  /// Type of coupling strategy for POROMULTIPHASE problems
  enum SolutionSchemeOverFields
  {
    solscheme_undefined,
    solscheme_twoway
  };


  /// set the POROMULTIPHASE parameters
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

}  // namespace POROMULTIPHASE

}  // namespace INPAR



#endif /* SRC_DRT_INPAR_INPAR_POROMULTIPHASE_H_ */

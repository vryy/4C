/*----------------------------------------------------------------------*/
/*!
 \file inpar_poromultiphase_scatra.H

 \brief input parameters for porous multiphase problem with scalar transport

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_INPAR_INPAR_POROMULTIPHASE_SCATRA_H_
#define SRC_DRT_INPAR_INPAR_POROMULTIPHASE_SCATRA_H_

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace INPAR
{

/// Input parameters for poroelasticity with scalar transport
namespace POROMULTIPHASESCATRA
{

  /// Type of coupling strategy for poro scatra problems
  enum SolutionSchemeOverFields
  {
    solscheme_twoway
  };

  /// set the poroscatra parameters
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

}  // namespace PORO_SCATRA

}  // namespace INPAR



#endif /* SRC_DRT_INPAR_INPAR_POROMULTIPHASE_SCATRA_H_ */

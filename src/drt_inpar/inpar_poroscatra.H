/*!------------------------------------------------------------------------------------------------*
\file inpar_poroscatra.H

\brief

   \level 2

   \maintainer  Andreas Rauch
                rauch@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *------------------------------------------------------------------------------------------------*/

#ifndef INPAR_POROSCATRA_H_
#define INPAR_POROSCATRA_H_

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace INPAR
{
  /// Input parameters for poroelasticity with scalar transport
  namespace PORO_SCATRA
  {
    /// Type of coupling strategy for poro scatra problems
    enum SolutionSchemeOverFields
    {
      Monolithic,
      Part_ScatraToPoro,
      Part_PoroToScatra,
      Part_TwoWay
      //   Monolithic
    };

    /// set the poroscatra parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace PORO_SCATRA

}  // namespace INPAR

/*----------------------------------------------------------------------*/



#endif /* INPAR_POROSCATRA_H_ */

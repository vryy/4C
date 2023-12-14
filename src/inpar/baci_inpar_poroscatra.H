/*----------------------------------------------------------------------*/
/*! \file
\brief

\level 2


 *------------------------------------------------------------------------------------------------*/

#ifndef BACI_INPAR_POROSCATRA_H
#define BACI_INPAR_POROSCATRA_H

#include "baci_inpar_parameterlist_utils.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace INPAR
{
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



BACI_NAMESPACE_CLOSE

#endif  // INPAR_POROSCATRA_H

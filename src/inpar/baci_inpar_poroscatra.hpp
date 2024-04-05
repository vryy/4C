/*----------------------------------------------------------------------*/
/*! \file
\brief

\level 2


 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_INPAR_POROSCATRA_HPP
#define FOUR_C_INPAR_POROSCATRA_HPP

#include "baci_config.hpp"

#include "baci_utils_parameter_list.hpp"

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

#endif

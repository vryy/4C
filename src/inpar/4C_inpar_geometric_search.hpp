/*-----------------------------------------------------------*/
/*! \file

\brief input parameter for geometric search strategy

\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_INPAR_GEOMETRIC_SEARCH_HPP
#define FOUR_C_INPAR_GEOMETRIC_SEARCH_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Inpar::GeometricSearch
{
  //! set the parameters for the geometric search strategy
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

}  // namespace Inpar::GeometricSearch

FOUR_C_NAMESPACE_CLOSE

#endif

/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameters for model order reduction

\level 2

*/

/*----------------------------------------------------------------------*/
#ifndef FOUR_C_INPAR_MOR_HPP
#define FOUR_C_INPAR_MOR_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace Inpar
{
  namespace ModelOrderRed
  {
    /// Defines all valid parameters for model order reduction
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace ModelOrderRed
}  // namespace Inpar

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif

/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameters for model order reduction

\level 2

*/

/*----------------------------------------------------------------------*/
#ifndef FOUR_C_INPAR_MOR_HPP
#define FOUR_C_INPAR_MOR_HPP

#include "baci_config.hpp"

#include "baci_utils_parameter_list.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace MOR
  {
    /// Defines all valid parameters for model order reduction
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace MOR
}  // namespace INPAR

/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif

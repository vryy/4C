/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for cut library

\level 2


*/

/*----------------------------------------------------------------------*/
#ifndef FOUR_C_INPAR_CUT_HPP
#define FOUR_C_INPAR_CUT_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace Inpar
{
  namespace Cut
  {
    /// set the cut parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace Cut

}  // namespace Inpar

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif

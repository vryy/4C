#ifndef FOUR_C_INPAR_CUT_HPP
#define FOUR_C_INPAR_CUT_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace Inpar
{
  namespace Cut
  {
    /// set the cut parameters
    void set_valid_parameters(Teuchos::ParameterList& list);

  }  // namespace Cut

}  // namespace Inpar

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif

#ifndef FOUR_C_INPAR_IO_RUNTIME_OUTPUT_HPP
#define FOUR_C_INPAR_IO_RUNTIME_OUTPUT_HPP


/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace Inpar::IORuntimeOutput
{

  /// set the valid parameters related to writing of output at runtime
  void set_valid_parameters(Teuchos::ParameterList& list);

}  // namespace Inpar::IORuntimeOutput

FOUR_C_NAMESPACE_CLOSE

#endif

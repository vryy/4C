/*----------------------------------------------------------------------*/
/*! \file

\brief input parameters for output of a fluid field at runtime

\level 2

*/
/*----------------------------------------------------------------------*/
/* definitions */
#ifndef FOUR_C_INPAR_IO_RUNTIME_OUTPUT_FLUID_HPP
#define FOUR_C_INPAR_IO_RUNTIME_OUTPUT_FLUID_HPP


/*----------------------------------------------------------------------*/
/* headers */
#include "baci_config.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace IO_RUNTIME_OUTPUT
  {
    namespace FLUID
    {
      /// set the valid parameters related to writing of output at runtime
      void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    }  // namespace FLUID
  }    // namespace IO_RUNTIME_OUTPUT
}  // namespace INPAR

BACI_NAMESPACE_CLOSE

#endif

/*----------------------------------------------------------------------*/
/*! \file

\brief input parameters for output of beam discretization at runtime

\level 3

*/
/*----------------------------------------------------------------------*/
/* definitions */
#ifndef FOUR_C_INPAR_IO_RUNTIME_OUTPUT_STRUCTURE_BEAMS_HPP
#define FOUR_C_INPAR_IO_RUNTIME_OUTPUT_STRUCTURE_BEAMS_HPP


/*----------------------------------------------------------------------*/
/* headers */
#include "baci_config.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace IO_RUNTIME_OUTPUT
  {
    namespace BEAMS
    {
      /// set the valid parameters related to writing of output at runtime
      void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    }  // namespace BEAMS
  }    // namespace IO_RUNTIME_OUTPUT
}  // namespace INPAR

FOUR_C_NAMESPACE_CLOSE

#endif

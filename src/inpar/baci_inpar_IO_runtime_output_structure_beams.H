/*----------------------------------------------------------------------*/
/*! \file

\brief input parameters for output of beam discretization at runtime

\level 3

*/
/*----------------------------------------------------------------------*/
/* definitions */
#ifndef BACI_INPAR_IO_RUNTIME_OUTPUT_STRUCTURE_BEAMS_H
#define BACI_INPAR_IO_RUNTIME_OUTPUT_STRUCTURE_BEAMS_H


/*----------------------------------------------------------------------*/
/* headers */
#include "baci_config.H"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN


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

BACI_NAMESPACE_CLOSE

#endif

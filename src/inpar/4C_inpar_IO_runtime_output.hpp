/*----------------------------------------------------------------------*/
/*! \file

\brief input parameters for output of structural problem at runtime

\level 2

*/
/*----------------------------------------------------------------------*/
/* definitions */
#ifndef FOUR_C_INPAR_IO_RUNTIME_OUTPUT_HPP
#define FOUR_C_INPAR_IO_RUNTIME_OUTPUT_HPP


/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace INPAR::IO_RUNTIME_OUTPUT
{
  /// data format for written numeric data
  enum class OutputDataFormat
  {
    binary,
    ascii,
    vague
  };

  // Specify the output writer that shall be used
  enum class OutputWriter
  {
    none,
    vtu_per_rank  // Write one file per time step per rank in the vtu format
  };

  /// set the valid parameters related to writing of output at runtime
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

}  // namespace INPAR::IO_RUNTIME_OUTPUT

FOUR_C_NAMESPACE_CLOSE

#endif

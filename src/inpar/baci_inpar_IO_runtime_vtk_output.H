/*----------------------------------------------------------------------*/
/*! \file

\brief input parameters for VTK output of structural problem at runtime

\level 2

*/
/*----------------------------------------------------------------------*/
/* definitions */
#ifndef BACI_INPAR_IO_RUNTIME_VTK_OUTPUT_H
#define BACI_INPAR_IO_RUNTIME_VTK_OUTPUT_H


/*----------------------------------------------------------------------*/
/* headers */
#include "baci_config.H"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace INPAR::IO_RUNTIME_VTK
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

  /// set the valid parameters related to writing of VTK output at runtime
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

}  // namespace INPAR::IO_RUNTIME_VTK

BACI_NAMESPACE_CLOSE

#endif

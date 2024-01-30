/*----------------------------------------------------------------------*/
/*! \file
\brief input parameters monitoring dirichlet boundary conditions

\level 2

*/
/*----------------------------------------------------------------------*/
/* definitions */
#ifndef BACI_INPAR_IO_MONITOR_STRUCTURE_DBC_H
#define BACI_INPAR_IO_MONITOR_STRUCTURE_DBC_H


/*----------------------------------------------------------------------*/
/* headers */
#include "baci_config.H"

#include "baci_utils_exceptions.H"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace IO_MONITOR_STRUCTURE_DBC
  {
    /// data format for written numeric data
    enum FileType
    {
      csv,
      data
    };

    /// set the valid parameters related to writing of output at runtime
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace IO_MONITOR_STRUCTURE_DBC
}  // namespace INPAR

BACI_NAMESPACE_CLOSE

#endif

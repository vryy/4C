/*----------------------------------------------------------------------*/
/*! \file
\brief input parameters monitoring dirichlet boundary conditions

\level 2

*/
/*----------------------------------------------------------------------*/
/* definitions */
#ifndef INPAR_IO_MONITOR_STRUCTURE_DBC
#define INPAR_IO_MONITOR_STRUCTURE_DBC


/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../drt_lib/drt_dserror.H"

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

    /// set the valid parameters related to writing of VTK output at runtime
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace IO_MONITOR_STRUCTURE_DBC
}  // namespace INPAR

#endif

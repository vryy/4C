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
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace IO_RUNTIME_VTK
  {
    /// data format for written numeric data
    enum OutputDataFormat
    {
      binary,
      ascii,
      vague
    };

    /// set the valid parameters related to writing of VTK output at runtime
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace IO_RUNTIME_VTK
}  // namespace INPAR

#endif

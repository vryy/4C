/*----------------------------------------------------------------------*/
/*! \file
\brief input parameters for VTP output of structural problem at runtime

\level 2

*/
/*----------------------------------------------------------------------*/
/* definitions */
#ifndef BACI_INPAR_IO_RUNTIME_VTP_OUTPUT_STRUCTURE_H
#define BACI_INPAR_IO_RUNTIME_VTP_OUTPUT_STRUCTURE_H


/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace IO_RUNTIME_VTP_STRUCTURE
  {
    /// data format for written numeric data
    enum OutputDataFormat
    {
      binary,
      ascii,
      vague
    };

    /// set the valid parameters related to writing of VTP output at runtime
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace IO_RUNTIME_VTP_STRUCTURE
}  // namespace INPAR

#endif

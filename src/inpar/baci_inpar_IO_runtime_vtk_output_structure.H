/*----------------------------------------------------------------------*/
/*! \file

\brief input parameters for VTK output of structural problem at runtime

\level 2

*/
/*----------------------------------------------------------------------*/
/* definitions */
#ifndef BACI_INPAR_IO_RUNTIME_VTK_OUTPUT_STRUCTURE_H
#define BACI_INPAR_IO_RUNTIME_VTK_OUTPUT_STRUCTURE_H


/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace IO_RUNTIME_VTK
  {
    namespace STRUCTURE
    {
      /// set the valid parameters related to writing of VTK output at runtime
      void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    }  // namespace STRUCTURE
  }    // namespace IO_RUNTIME_VTK
}  // namespace INPAR

#endif

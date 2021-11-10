/*----------------------------------------------------------------------*/
/*! \file

\brief input parameters for VTK output of a fluid field at runtime

\level 2

*/
/*----------------------------------------------------------------------*/
/* definitions */
#ifndef INPAR_IO_RUNTIME_VTK_OUTPUT_FLUID_H
#define INPAR_IO_RUNTIME_VTK_OUTPUT_FLUID_H


/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace IO_RUNTIME_VTK
  {
    namespace FLUID
    {
      /// set the valid parameters related to writing of VTK output at runtime
      void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    }  // namespace FLUID
  }    // namespace IO_RUNTIME_VTK
}  // namespace INPAR

#endif

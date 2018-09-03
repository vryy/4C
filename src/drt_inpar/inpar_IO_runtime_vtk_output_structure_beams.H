/*----------------------------------------------------------------------*/
/*!
\file inpar_IO_runtime_vtk_output_structure_beams.H

\brief input parameters for VTK output of beam discretization at runtime

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/
/* definitions */
#ifndef INPAR_IO_RUNTIME_VTK_OUTPUT_STRUCTURE_BEAMS_H
#define INPAR_IO_RUNTIME_VTK_OUTPUT_STRUCTURE_BEAMS_H


/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>


/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace IO_RUNTIME_VTK
  {
    namespace BEAMS
    {
      /// set the valid parameters related to writing of VTK output at runtime
      void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    }  // namespace BEAMS
  }    // namespace IO_RUNTIME_VTK
}  // namespace INPAR

#endif

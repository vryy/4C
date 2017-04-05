/*----------------------------------------------------------------------*/
/*!
\file inpar_IO_runtime_vtk_output_structure.H

\brief input parameters for VTK output of structural problem at runtime

\level 2

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/
/* definitions */
#ifndef INPAR_IO_RUNTIME_VTK_OUTPUT_STRUCTURE_H
#define INPAR_IO_RUNTIME_VTK_OUTPUT_STRUCTURE_H


/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------*/
namespace INPAR
{
namespace IO_RUNTIME_VTK_STRUCTURE
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

}
}

#endif

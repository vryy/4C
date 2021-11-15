/*----------------------------------------------------------------------*/
/*! \file

\brief input parameters for VTK output of a fluid field at runtime

\level 2

*/
/*----------------------------------------------------------------------*/

#include "inpar_IO_runtime_vtk_output_fluid.H"

#include "drt_validparameters.H"
#include "inpar.H"
#include "inpar_parameterlist_utils.H"
#include "inpar_fluid.H"

#include <Teuchos_ParameterList.hpp>

namespace INPAR
{
  namespace IO_RUNTIME_VTK
  {
    namespace FLUID
    {
      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
      {
        using namespace DRT::INPUT;
        using Teuchos::setStringToIntegralParameter;
        using Teuchos::tuple;

        // related sublist
        Teuchos::ParameterList& sublist_IO = list->sublist("IO", false, "");
        Teuchos::ParameterList& sublist_IO_VTK =
            sublist_IO.sublist("RUNTIME VTK OUTPUT", false, "");
        Teuchos::ParameterList& sublist_IO_VTK_fluid = sublist_IO_VTK.sublist("FLUID", false, "");

        // whether to write output for fluid
        BoolParameter("OUTPUT_FLUID", "No", "write fluid output", &sublist_IO_VTK_fluid);

        // whether to write velocity state
        BoolParameter("VELOCITY", "No", "write velocity output", &sublist_IO_VTK_fluid);

        // whether to write pressure state
        BoolParameter("PRESSURE", "No", "write pressure output", &sublist_IO_VTK_fluid);

        // whether to write acceleration state
        BoolParameter("ACCELERATION", "No", "write acceleration output", &sublist_IO_VTK_fluid);

        // whether to write displacement state
        BoolParameter("DISPLACEMENT", "No", "write displacement output", &sublist_IO_VTK_fluid);

        // whether to write displacement state
        BoolParameter("GRIDVELOCITY", "No", "write grid velocity output", &sublist_IO_VTK_fluid);

        // whether to write element owner
        BoolParameter("ELEMENT_OWNER", "No", "write element owner", &sublist_IO_VTK_fluid);

        // whether to write element GIDs
        BoolParameter(
            "ELEMENT_GID", "No", "write baci internal element GIDs", &sublist_IO_VTK_fluid);

        // whether to write node GIDs
        BoolParameter("NODE_GID", "No", "write baci internal node GIDs", &sublist_IO_VTK_fluid);
      }
    }  // namespace FLUID
  }    // namespace IO_RUNTIME_VTK
}  // namespace INPAR

/*----------------------------------------------------------------------*/
/*! \file

\brief input parameters for output of a fluid field at runtime

\level 2

*/
/*----------------------------------------------------------------------*/

#include "baci_inpar_IO_runtime_output_fluid.hpp"

#include "baci_inpar.hpp"
#include "baci_inpar_fluid.hpp"
#include "baci_inpar_parameterlist_utils.hpp"
#include "baci_inpar_validparameters.hpp"

#include <Teuchos_ParameterList.hpp>

BACI_NAMESPACE_OPEN

namespace INPAR
{
  namespace IO_RUNTIME_OUTPUT
  {
    namespace FLUID
    {
      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
      {
        using namespace INPUT;
        using Teuchos::setStringToIntegralParameter;
        using Teuchos::tuple;

        // related sublist
        Teuchos::ParameterList& sublist_IO = list->sublist("IO", false, "");
        Teuchos::ParameterList& sublist_IO_output =
            sublist_IO.sublist("RUNTIME VTK OUTPUT", false, "");
        Teuchos::ParameterList& sublist_IO_output_fluid =
            sublist_IO_output.sublist("FLUID", false, "");

        // whether to write output for fluid
        BoolParameter("OUTPUT_FLUID", "No", "write fluid output", &sublist_IO_output_fluid);

        // whether to write velocity state
        BoolParameter("VELOCITY", "No", "write velocity output", &sublist_IO_output_fluid);

        // whether to write pressure state
        BoolParameter("PRESSURE", "No", "write pressure output", &sublist_IO_output_fluid);

        // whether to write acceleration state
        BoolParameter("ACCELERATION", "No", "write acceleration output", &sublist_IO_output_fluid);

        // whether to write displacement state
        BoolParameter("DISPLACEMENT", "No", "write displacement output", &sublist_IO_output_fluid);

        // whether to write displacement state
        BoolParameter("GRIDVELOCITY", "No", "write grid velocity output", &sublist_IO_output_fluid);

        // whether to write element owner
        BoolParameter("ELEMENT_OWNER", "No", "write element owner", &sublist_IO_output_fluid);

        // whether to write element GIDs
        BoolParameter(
            "ELEMENT_GID", "No", "write baci internal element GIDs", &sublist_IO_output_fluid);

        // whether to write node GIDs
        BoolParameter("NODE_GID", "No", "write baci internal node GIDs", &sublist_IO_output_fluid);
      }
    }  // namespace FLUID
  }    // namespace IO_RUNTIME_OUTPUT
}  // namespace INPAR

BACI_NAMESPACE_CLOSE

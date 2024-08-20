/*----------------------------------------------------------------------*/
/*! \file

\brief input parameters for output of a fluid field at runtime

\level 2

*/
/*----------------------------------------------------------------------*/

#include "4C_inpar_IO_runtime_output_fluid.hpp"

#include "4C_utils_parameter_list.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Inpar
{
  namespace IORuntimeOutput
  {
    namespace FLUID
    {
      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      void set_valid_parameters(Teuchos::RCP<Teuchos::ParameterList> list)
      {
        using Teuchos::setStringToIntegralParameter;
        using Teuchos::tuple;

        // related sublist
        Teuchos::ParameterList& sublist_IO = list->sublist("IO", false, "");
        Teuchos::ParameterList& sublist_IO_output =
            sublist_IO.sublist("RUNTIME VTK OUTPUT", false, "");
        Teuchos::ParameterList& sublist_IO_output_fluid =
            sublist_IO_output.sublist("FLUID", false, "");

        // whether to write output for fluid
        Core::UTILS::bool_parameter(
            "OUTPUT_FLUID", "No", "write fluid output", &sublist_IO_output_fluid);

        // whether to write velocity state
        Core::UTILS::bool_parameter(
            "VELOCITY", "No", "write velocity output", &sublist_IO_output_fluid);

        // whether to write pressure state
        Core::UTILS::bool_parameter(
            "PRESSURE", "No", "write pressure output", &sublist_IO_output_fluid);

        // whether to write acceleration state
        Core::UTILS::bool_parameter(
            "ACCELERATION", "No", "write acceleration output", &sublist_IO_output_fluid);

        // whether to write displacement state
        Core::UTILS::bool_parameter(
            "DISPLACEMENT", "No", "write displacement output", &sublist_IO_output_fluid);

        // whether to write displacement state
        Core::UTILS::bool_parameter(
            "GRIDVELOCITY", "No", "write grid velocity output", &sublist_IO_output_fluid);

        // whether to write element owner
        Core::UTILS::bool_parameter(
            "ELEMENT_OWNER", "No", "write element owner", &sublist_IO_output_fluid);

        // whether to write element GIDs
        Core::UTILS::bool_parameter(
            "ELEMENT_GID", "No", "write 4C internal element GIDs", &sublist_IO_output_fluid);

        // whether to write node GIDs
        Core::UTILS::bool_parameter(
            "NODE_GID", "No", "write 4C internal node GIDs", &sublist_IO_output_fluid);
      }
    }  // namespace FLUID
  }    // namespace IORuntimeOutput
}  // namespace Inpar

FOUR_C_NAMESPACE_CLOSE

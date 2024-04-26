/*----------------------------------------------------------------------*/
/*! \file

\brief input parameters for output of a fluid field at runtime

\level 2

*/
/*----------------------------------------------------------------------*/

#include "4C_inpar_IO_runtime_output_fluid.hpp"

#include "4C_inpar.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_utils_parameter_list.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

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
        CORE::UTILS::BoolParameter(
            "OUTPUT_FLUID", "No", "write fluid output", &sublist_IO_output_fluid);

        // whether to write velocity state
        CORE::UTILS::BoolParameter(
            "VELOCITY", "No", "write velocity output", &sublist_IO_output_fluid);

        // whether to write pressure state
        CORE::UTILS::BoolParameter(
            "PRESSURE", "No", "write pressure output", &sublist_IO_output_fluid);

        // whether to write acceleration state
        CORE::UTILS::BoolParameter(
            "ACCELERATION", "No", "write acceleration output", &sublist_IO_output_fluid);

        // whether to write displacement state
        CORE::UTILS::BoolParameter(
            "DISPLACEMENT", "No", "write displacement output", &sublist_IO_output_fluid);

        // whether to write displacement state
        CORE::UTILS::BoolParameter(
            "GRIDVELOCITY", "No", "write grid velocity output", &sublist_IO_output_fluid);

        // whether to write element owner
        CORE::UTILS::BoolParameter(
            "ELEMENT_OWNER", "No", "write element owner", &sublist_IO_output_fluid);

        // whether to write element GIDs
        CORE::UTILS::BoolParameter(
            "ELEMENT_GID", "No", "write 4C internal element GIDs", &sublist_IO_output_fluid);

        // whether to write node GIDs
        CORE::UTILS::BoolParameter(
            "NODE_GID", "No", "write 4C internal node GIDs", &sublist_IO_output_fluid);
      }
    }  // namespace FLUID
  }    // namespace IO_RUNTIME_OUTPUT
}  // namespace INPAR

FOUR_C_NAMESPACE_CLOSE

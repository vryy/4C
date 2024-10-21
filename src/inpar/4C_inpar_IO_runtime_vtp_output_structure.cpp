#include "4C_inpar_IO_runtime_vtp_output_structure.hpp"

#include "4C_io_geometry_type.hpp"
#include "4C_utils_parameter_list.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Inpar
{
  namespace IORuntimeVTPStructure
  {
    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    void set_valid_parameters(Teuchos::ParameterList& list)
    {
      using Teuchos::setStringToIntegralParameter;
      using Teuchos::tuple;

      // related sublist
      Teuchos::ParameterList& sublist_IO = list.sublist("IO", false, "");
      Teuchos::ParameterList& sublist_IO_VTP_structure =
          sublist_IO.sublist("RUNTIME VTP OUTPUT STRUCTURE", false, "");


      // output interval regarding steps: write output every INTERVAL_STEPS steps
      Core::Utils::int_parameter("INTERVAL_STEPS", -1,
          "write VTP output at runtime every INTERVAL_STEPS steps", &sublist_IO_VTP_structure);

      Core::Utils::int_parameter("STEP_OFFSET", 0,
          "An offset added to the current step to shift the steps to be written.",
          &sublist_IO_VTP_structure);

      // whether to write output in every iteration of the nonlinear solver
      Core::Utils::bool_parameter("EVERY_ITERATION", "No",
          "write output in every iteration of the nonlinear solver", &sublist_IO_VTP_structure);

      // write owner at every visualization point
      Core::Utils::bool_parameter(
          "OWNER", "No", "write owner of every point", &sublist_IO_VTP_structure);

      // write orientation at every visualization point
      Core::Utils::bool_parameter("ORIENTATIONANDLENGTH", "No", "write orientation at every point",
          &sublist_IO_VTP_structure);

      // write number of bonds at every visualization point
      Core::Utils::bool_parameter(
          "NUMBEROFBONDS", "No", "write number of bonds of every point", &sublist_IO_VTP_structure);

      // write force actin in linker
      Core::Utils::bool_parameter(
          "LINKINGFORCE", "No", "write force acting in linker", &sublist_IO_VTP_structure);
    }


  }  // namespace IORuntimeVTPStructure
}  // namespace Inpar

FOUR_C_NAMESPACE_CLOSE

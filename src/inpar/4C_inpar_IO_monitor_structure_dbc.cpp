#include "4C_inpar_IO_monitor_structure_dbc.hpp"

#include "4C_io_geometry_type.hpp"
#include "4C_utils_parameter_list.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Inpar
{
  namespace IOMonitorStructureDBC
  {
    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    void set_valid_parameters(Teuchos::ParameterList& list)
    {
      using Teuchos::setStringToIntegralParameter;
      using Teuchos::tuple;

      // related sublist
      Teuchos::ParameterList& sublist_IO = list.sublist("IO", false, "");
      Teuchos::ParameterList& sublist_IO_monitor_structure_dbc =
          sublist_IO.sublist("MONITOR STRUCTURE DBC", false, "");

      // output interval regarding steps: write output every INTERVAL_STEPS steps
      Core::Utils::int_parameter("INTERVAL_STEPS", -1,
          "write reaction force output every INTERVAL_STEPS steps",
          &sublist_IO_monitor_structure_dbc);

      // precision for file
      Core::Utils::int_parameter(
          "PRECISION_FILE", 16, "precision for written file", &sublist_IO_monitor_structure_dbc);

      // precision for screen
      Core::Utils::int_parameter("PRECISION_SCREEN", 5, "precision for written screen output",
          &sublist_IO_monitor_structure_dbc);

      // type of written output file
      setStringToIntegralParameter<Inpar::IOMonitorStructureDBC::FileType>("FILE_TYPE", "csv",
          "type of written output file",
          tuple<std::string>("csv", "CSV", "Csv", "data", "Data", "DATA"),
          tuple<Inpar::IOMonitorStructureDBC::FileType>(Inpar::IOMonitorStructureDBC::csv,
              Inpar::IOMonitorStructureDBC::csv, Inpar::IOMonitorStructureDBC::csv,
              Inpar::IOMonitorStructureDBC::data, Inpar::IOMonitorStructureDBC::data,
              Inpar::IOMonitorStructureDBC::data),
          &sublist_IO_monitor_structure_dbc);

      // whether to write output in every iteration of the nonlinear solver
      Core::Utils::bool_parameter("WRITE_HEADER", "No",
          "write information about monitored boundary condition to output file",
          &sublist_IO_monitor_structure_dbc);
    }
  }  // namespace IOMonitorStructureDBC
}  // namespace Inpar

FOUR_C_NAMESPACE_CLOSE

/*----------------------------------------------------------------------*/
/*! \file
\brief input parameters monitoring dirichlet boundary conditions

\level 2

*/
/*----------------------------------------------------------------------*/

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
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
    {
      using Teuchos::setStringToIntegralParameter;
      using Teuchos::tuple;

      // related sublist
      Teuchos::ParameterList& sublist_IO = list->sublist("IO", false, "");
      Teuchos::ParameterList& sublist_IO_monitor_structure_dbc =
          sublist_IO.sublist("MONITOR STRUCTURE DBC", false, "");

      // output interval regarding steps: write output every INTERVAL_STEPS steps
      Core::UTILS::IntParameter("INTERVAL_STEPS", -1,
          "write reaction force output every INTERVAL_STEPS steps",
          &sublist_IO_monitor_structure_dbc);

      // precision for file
      Core::UTILS::IntParameter(
          "PRECISION_FILE", 16, "precision for written file", &sublist_IO_monitor_structure_dbc);

      // precision for screen
      Core::UTILS::IntParameter("PRECISION_SCREEN", 5, "precision for written screen output",
          &sublist_IO_monitor_structure_dbc);

      // type of written output file
      setStringToIntegralParameter<int>("FILE_TYPE", "csv", "type of written output file",
          tuple<std::string>("csv", "CSV", "Csv", "data", "Data", "DATA"),
          tuple<int>(Inpar::IOMonitorStructureDBC::csv, Inpar::IOMonitorStructureDBC::csv,
              Inpar::IOMonitorStructureDBC::csv, Inpar::IOMonitorStructureDBC::data,
              Inpar::IOMonitorStructureDBC::data, Inpar::IOMonitorStructureDBC::data),
          &sublist_IO_monitor_structure_dbc);

      // whether to write output in every iteration of the nonlinear solver
      Core::UTILS::BoolParameter("WRITE_HEADER", "No",
          "write information about monitored boundary condition to output file",
          &sublist_IO_monitor_structure_dbc);
    }
  }  // namespace IOMonitorStructureDBC
}  // namespace Inpar

FOUR_C_NAMESPACE_CLOSE

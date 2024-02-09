/*----------------------------------------------------------------------*/
/*! \file
\brief input parameters monitoring dirichlet boundary conditions

\level 2

*/
/*----------------------------------------------------------------------*/

#include "baci_inpar_IO_monitor_structure_dbc.hpp"

#include "baci_inpar.hpp"
#include "baci_inpar_parameterlist_utils.hpp"
#include "baci_inpar_validparameters.hpp"

#include <Teuchos_ParameterList.hpp>

BACI_NAMESPACE_OPEN

namespace INPAR
{
  namespace IO_MONITOR_STRUCTURE_DBC
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
      Teuchos::ParameterList& sublist_IO_monitor_structure_dbc =
          sublist_IO.sublist("MONITOR STRUCTURE DBC", false, "");

      // output interval regarding steps: write output every INTERVAL_STEPS steps
      IntParameter("INTERVAL_STEPS", -1, "write reaction force output every INTERVAL_STEPS steps",
          &sublist_IO_monitor_structure_dbc);

      // precision for file
      IntParameter(
          "PRECISION_FILE", 16, "precision for written file", &sublist_IO_monitor_structure_dbc);

      // precision for screen
      IntParameter("PRECISION_SCREEN", 5, "precision for written screen output",
          &sublist_IO_monitor_structure_dbc);

      // type of written output file
      setStringToIntegralParameter<int>("FILE_TYPE", "csv", "type of written output file",
          tuple<std::string>("csv", "CSV", "Csv", "data", "Data", "DATA"),
          tuple<int>(INPAR::IO_MONITOR_STRUCTURE_DBC::csv, INPAR::IO_MONITOR_STRUCTURE_DBC::csv,
              INPAR::IO_MONITOR_STRUCTURE_DBC::csv, INPAR::IO_MONITOR_STRUCTURE_DBC::data,
              INPAR::IO_MONITOR_STRUCTURE_DBC::data, INPAR::IO_MONITOR_STRUCTURE_DBC::data),
          &sublist_IO_monitor_structure_dbc);

      // whether to write output in every iteration of the nonlinear solver
      BoolParameter("WRITE_HEADER", "No",
          "write information about monitored boundary condition to output file",
          &sublist_IO_monitor_structure_dbc);
    }
  }  // namespace IO_MONITOR_STRUCTURE_DBC
}  // namespace INPAR

BACI_NAMESPACE_CLOSE

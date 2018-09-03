/*----------------------------------------------------------------------*/
/*!
\file inpar_IO_monitor_structure_dbc.cpp

\brief input parameters monitoring dirichlet boundary conditions

\level 2

\maintainer Michael Hiermeier
*/
/*----------------------------------------------------------------------*/

#include "inpar_IO_monitor_structure_dbc.H"

#include "drt_validparameters.H"
#include "inpar.H"
#include "inpar_parameterlist_utils.H"

#include <Teuchos_ParameterList.hpp>

namespace INPAR
{
namespace IO_MONITOR_STRUCTURE_DBC
{

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
  {
    using namespace DRT::INPUT;
    using Teuchos::tuple;
    using Teuchos::setStringToIntegralParameter;

    Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
    Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

    // related sublist
    Teuchos::ParameterList& sublist_IO = list->sublist("IO",false,"");
    Teuchos::ParameterList& sublist_IO_monitor_structure_dbc =
        sublist_IO.sublist("MONITOR STRUCTURE DBC",false,"");


    // precision for file
    IntParameter( "PRECISION_FILE", 16,
        "precision for written file", &sublist_IO_monitor_structure_dbc );

    // precision for screen
    IntParameter( "PRECISION_SCREEN", 5,
        "precision for written screen output", &sublist_IO_monitor_structure_dbc );

    // type of written output file
    setStringToIntegralParameter<int>(
      "FILE_TYPE", "csv", "type of written output file",
      tuple<std::string>(
        "csv",
        "CSV",
        "Csv",
        "data",
        "Data",
        "DATA"),
      tuple<int>(
        INPAR::IO_MONITOR_STRUCTURE_DBC::csv,
        INPAR::IO_MONITOR_STRUCTURE_DBC::csv,
        INPAR::IO_MONITOR_STRUCTURE_DBC::csv,
        INPAR::IO_MONITOR_STRUCTURE_DBC::data,
        INPAR::IO_MONITOR_STRUCTURE_DBC::data,
        INPAR::IO_MONITOR_STRUCTURE_DBC::data),
      &sublist_IO_monitor_structure_dbc);


    // whether to write output in every iteration of the nonlinear solver
    setStringToIntegralParameter<int>("WRITE_HEADER","No",
                                 "write information about monitored boundary condition to output file",
                                 yesnotuple, yesnovalue, &sublist_IO_monitor_structure_dbc);

  }


}
}


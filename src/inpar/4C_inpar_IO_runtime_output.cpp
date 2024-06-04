/*----------------------------------------------------------------------*/
/*! \file

\brief input parameters for output of structural problem at runtime

\level 2

*/
/*----------------------------------------------------------------------*/

#include "4C_inpar_IO_runtime_output.hpp"

#include "4C_io_visualization_parameters.hpp"
#include "4C_utils_parameter_list.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace INPAR
{
  namespace IO_RUNTIME_OUTPUT
  {
    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
    {
      using Teuchos::setStringToIntegralParameter;
      using Teuchos::tuple;

      // related sublist
      Teuchos::ParameterList& sublist_IO = list->sublist("IO", false, "");
      Teuchos::ParameterList& sublist_IO_VTK_structure =
          sublist_IO.sublist("RUNTIME VTK OUTPUT", false, "");


      // output interval regarding steps: write output every INTERVAL_STEPS steps
      CORE::UTILS::IntParameter("INTERVAL_STEPS", -1,
          "write visualization output at runtime every INTERVAL_STEPS steps",
          &sublist_IO_VTK_structure);


      CORE::UTILS::IntParameter("STEP_OFFSET", 0,
          "An offset added to the current step to shift the steps to be written.",
          &sublist_IO_VTK_structure);


      // data format for written numeric data
      setStringToIntegralParameter<CORE::IO::OutputDataFormat>("OUTPUT_DATA_FORMAT", "binary",
          "data format for written numeric data", tuple<std::string>("binary", "ascii"),
          tuple<CORE::IO::OutputDataFormat>(
              CORE::IO::OutputDataFormat::binary, CORE::IO::OutputDataFormat::ascii),
          &sublist_IO_VTK_structure);

      // specify the maximum digits in the number of time steps that shall be written
      CORE::UTILS::IntParameter("TIMESTEP_RESERVE_DIGITS", 5,
          "Specify the maximum digits in the number of time steps that shall be written. This only "
          "affects the number of leading zeros in the output file names.",
          &sublist_IO_VTK_structure);

      // whether to write output in every iteration of the nonlinear solver
      CORE::UTILS::BoolParameter("EVERY_ITERATION", "No",
          "write output in every iteration of the nonlinear solver", &sublist_IO_VTK_structure);

      // virtual time increment that is added for each nonlinear output state
      CORE::UTILS::DoubleParameter("EVERY_ITERATION_VIRTUAL_TIME_INCREMENT", 1e-8,
          "Specify the virtual time increment that is added for each nonlinear output state",
          &sublist_IO_VTK_structure);

      // specify the maximum digits in the number of iterations that shall be written
      CORE::UTILS::IntParameter("EVERY_ITERATION_RESERVE_DIGITS", 4,
          "Specify the maximum digits in the number of iterations that shall be written. This only "
          "affects the number of leading zeros in the output file names.",
          &sublist_IO_VTK_structure);

      // specify the actual visualization writer
      setStringToIntegralParameter<CORE::IO::OutputWriter>("OUTPUT_WRITER", "vtu_per_rank",
          "Specify which output writer shall be used to write the visualization data to disk",
          tuple<std::string>("vtu_per_rank"),
          tuple<CORE::IO::OutputWriter>(CORE::IO::OutputWriter::vtu_per_rank),
          &sublist_IO_VTK_structure);
    }


  }  // namespace IO_RUNTIME_OUTPUT
}  // namespace INPAR

FOUR_C_NAMESPACE_CLOSE

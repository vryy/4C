/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Define a data container that stores visualization data in the VTU data format

\level 0

*/
/*-----------------------------------------------------------------------------------------------*/


#include "baci_io_visualization_parameters.H"

#include "baci_inpar_IO_runtime_vtk_output.H"
#include "baci_inpar_parameterlist_utils.H"
#include "baci_utils_exceptions.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

BACI_NAMESPACE_OPEN


/**
 *
 */
IO::VisualizationParameters IO::VisualizationParametersFactory(
    const Teuchos::ParameterList& visualization_ouput_parameter_list)
{
  IO::VisualizationParameters parameters;

  // Data format
  parameters.data_format_ = Teuchos::getIntegralValue<INPAR::IO_RUNTIME_VTK::OutputDataFormat>(
      visualization_ouput_parameter_list, "OUTPUT_DATA_FORMAT");

  // Number of digits to reserve for time step count
  parameters.digits_for_time_step_ =
      visualization_ouput_parameter_list.get<int>("TIMESTEP_RESERVE_DIGITS");

  // Parameters for output during the nonlinear solver
  parameters.every_iteration_virtual_time_increment_ =
      visualization_ouput_parameter_list.get<double>("EVERY_ITERATION_VIRTUAL_TIME_INCREMENT");
  parameters.digits_for_iteration_ =
      visualization_ouput_parameter_list.get<int>("EVERY_ITERATION_RESERVE_DIGITS");

  // This value can be overwritten from the physical field
  parameters.every_iteration_ =
      INPUT::IntegralValue<bool>(visualization_ouput_parameter_list, "EVERY_ITERATION");

  // Type of output writer
  const auto output_writer = Teuchos::getIntegralValue<INPAR::IO_RUNTIME_VTK::OutputWriter>(
      visualization_ouput_parameter_list, "OUTPUT_WRITER");
  if (output_writer == INPAR::IO_RUNTIME_VTK::OutputWriter::none)
  {
    dserror(
        "The visualization writer has to be set in the input file under IO/RUNTIME VTK "
        "OUTPUT/OUTPUT_WRITER");
  }
  parameters.writer_ = output_writer;

  return parameters;
}

/*
 *
 */
[[nodiscard]] int IO::GetTotalDigitsToReserveInTimeStep(
    const VisualizationParameters& visualization_parameters)
{
  if (visualization_parameters.every_iteration_)
  {
    return visualization_parameters.digits_for_time_step_ +
           visualization_parameters.digits_for_iteration_;
  }
  else
  {
    return visualization_parameters.digits_for_time_step_;
  }
}

/*
 *
 */
[[nodiscard]] std::pair<double, int> IO::GetTimeAndTimeStepIndexForOutput(
    const VisualizationParameters& visualization_parameters, const double time, const int step,
    const int iteration_number)
{
  if (iteration_number != 0 && !visualization_parameters.every_iteration_)
  {
    dserror(
        "An iteration index was passed to GetTimeAndTimeStepIndexForOutput, but output for "
        "each iteration is deactivated, this is a contradiction");
  }

  if (visualization_parameters.every_iteration_)
  {
    return {
        time + visualization_parameters.every_iteration_virtual_time_increment_ * iteration_number,
        step * std::pow(10, visualization_parameters.digits_for_iteration_) + iteration_number};
  }
  else
  {
    return {time, step};
  }
}
BACI_NAMESPACE_CLOSE

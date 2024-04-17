/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Define a data container that stores parameters for visualization output

\level 0

*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_IO_VISUALIZATION_PARAMETERS_HPP
#define FOUR_C_IO_VISUALIZATION_PARAMETERS_HPP


#include "baci_config.hpp"

#include "baci_inpar_IO_runtime_output.hpp"
#include "baci_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace IO
{
  class OutputControl;

  /**
   * @brief This struct holds parameters for visualization output
   */
  struct VisualizationParameters
  {
    //! Enum containing the type of output data format, i.e., binary or ascii.
    INPAR::IO_RUNTIME_OUTPUT::OutputDataFormat data_format_;

    //! Base output directory
    std::string directory_name_;

    //! Flag if output should be written for each nonlinear iteration
    bool every_iteration_;

    //! We need to add a small time increment for each iteration to avoid having "double" time steps
    double every_iteration_virtual_time_increment_;

    //! Only the prefix of the file name, without the directory
    std::string file_name_prefix_;

    //! Number of digits in the final "time_step_index" reserved for the nonlinear iteration
    int digits_for_iteration_;

    //! Number of digits in the final "time_step_index" reserved for the nonlinear iteration
    int digits_for_time_step_;

    //! Time the simulation is restarted from
    double restart_time_;

    //! In case of restart this prefix specifies the control file we read which might contain a path
    std::string restart_from_name_;

    //! Enum containing the output writer that shall be used
    INPAR::IO_RUNTIME_OUTPUT::OutputWriter writer_;
  };

  /**
   * @brief Create a container containing all visualization output parameters
   */
  [[nodiscard]] VisualizationParameters VisualizationParametersFactory(
      const Teuchos::ParameterList& visualization_ouput_parameter_list,
      const IO::OutputControl& output_control, double restart_time);

  /**
   * @brief Return the total number of digits to reserve in the time step numbering
   * @param visualization_parameters (in) Reference to a parameters container
   */
  [[nodiscard]] int GetTotalDigitsToReserveInTimeStep(
      const VisualizationParameters& visualization_parameters);

  /**
   * @brief Get the time step value and index that will be stored in the output file
   *
   * Since we have the option to output states during the nonlinear solution process, we use the
   * following convention for time step values and indices:
   * - For standard output the BACI internal time step value and index are used
   * - For output during each iteration, we increment the BACI internal time step value with a
   *   small increment for each nonlinear iteration, and we "add" the iteration index to the time
   *   step index, e.g., the time step index 0000270005 corresponds to the 5th nonlinear iteration
   *   in the 28th step -> 0000280000 is the converged 28th step.
   *
   * @param visualization_parameters (in) Reference to a parameters container
   * @param time (in) BACI internal time
   * @param step (in) BACI internal time step index of the last converged state
   * @param iteration_number (in) Number of nonlinear iteration
   */
  [[nodiscard]] std::pair<double, int> GetTimeAndTimeStepIndexForOutput(
      const VisualizationParameters& visualization_parameters, const double time, const int step,
      const int iteration_number = 0);
}  // namespace IO


FOUR_C_NAMESPACE_CLOSE

#endif

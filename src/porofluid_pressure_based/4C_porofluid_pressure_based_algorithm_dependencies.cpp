// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_porofluid_pressure_based_algorithm_dependencies.hpp"

#include "4C_global_data.hpp"
#include "4C_utils_function.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

PoroPressureBased::PorofluidAlgorithmDeps PoroPressureBased::make_algorithm_deps_from_problem(
    Global::Problem& problem)
{
  const bool has_artery_discretization = problem.does_exist_dis("artery");
  const std::shared_ptr<Core::FE::Discretization> artery_discretization =
      has_artery_discretization ? problem.get_dis("artery") : nullptr;

  return PorofluidAlgorithmDeps{
      .spatial_dimension = problem.n_dim(),
      .restart_step = problem.restart(),
      .verbosity =
          Teuchos::getIntegralValue<Core::IO::Verbositylevel>(problem.io_params(), "VERBOSITY"),
      .input_control_file = problem.input_control_file(),
      .output_control_file = problem.output_control_file(),
      .artery_discretization = artery_discretization,
      .runtime_vtk_output_parameters = &problem.io_params().sublist("RUNTIME VTK OUTPUT"),
      .artery_dynamic_parameters =
          has_artery_discretization ? &problem.arterial_dynamic_params() : nullptr,
      .function_manager = &problem.function_manager(),
      .function_of_space_time_by_id =
          [&problem](const int function_id) -> const Core::Utils::FunctionOfSpaceTime&
      { return problem.function_by_id<Core::Utils::FunctionOfSpaceTime>(function_id); },
      .solver_params_by_id = problem.solver_params_callback(),
      .add_field_test = [&problem](std::shared_ptr<Core::Utils::ResultTest> result_test)
      { problem.add_field_test(result_test); },
  };
}

FOUR_C_NAMESPACE_CLOSE

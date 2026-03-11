// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROFLUID_PRESSURE_BASED_ALGORITHM_DEPENDENCIES_HPP
#define FOUR_C_POROFLUID_PRESSURE_BASED_ALGORITHM_DEPENDENCIES_HPP

#include "4C_config.hpp"

#include "4C_io_pstream.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <functional>
#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  class InputControl;
  class OutputControl;
}  // namespace Core::IO

namespace Core::FE
{
  class Discretization;
}

namespace Core::Utils
{
  class FunctionManager;
  class FunctionOfSpaceTime;
  class ResultTest;
}  // namespace Core::Utils

namespace Global
{
  class Problem;
}  // namespace Global

namespace PoroPressureBased
{
  struct PorofluidAlgorithmDeps
  {
    int spatial_dimension = 0;
    int restart_step = 0;
    Core::IO::Verbositylevel verbosity = Core::IO::minimal;
    std::shared_ptr<Core::IO::InputControl> input_control_file;
    std::shared_ptr<Core::IO::OutputControl> output_control_file;
    std::shared_ptr<Core::FE::Discretization> artery_discretization;
    const Teuchos::ParameterList* runtime_vtk_output_parameters = nullptr;
    const Teuchos::ParameterList* artery_dynamic_parameters = nullptr;
    const Core::Utils::FunctionManager* function_manager = nullptr;
    std::function<const Core::Utils::FunctionOfSpaceTime&(int)> function_of_space_time_by_id;
    std::function<const Teuchos::ParameterList&(int)> solver_params_by_id;
    std::function<void(std::shared_ptr<Core::Utils::ResultTest>)> add_field_test;
  };

  PorofluidAlgorithmDeps make_algorithm_deps_from_problem(Global::Problem& problem);
}  // namespace PoroPressureBased

FOUR_C_NAMESPACE_CLOSE

#endif

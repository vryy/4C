// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_utils_function_library.hpp"

#include "4C_io_control.hpp"
#include "4C_io_file_reader.hpp"
#include "4C_io_input_spec_builders.hpp"
#include "4C_utils_cubic_spline_interpolation.hpp"
#include "4C_utils_function_manager.hpp"

#include <filesystem>
#include <memory>
#include <string>
#include <utility>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace
{

  std::shared_ptr<Core::Utils::FunctionOfScalar> create_library_function_scalar(
      const std::vector<Core::IO::InputParameterContainer>& parameters)
  {
    if (parameters.size() != 1) return nullptr;

    const auto& function_lin_def = parameters.front();

    if (function_lin_def.has_group("FASTPOLYNOMIAL"))
    {
      const auto& group = function_lin_def.group("FASTPOLYNOMIAL");
      std::vector<double> coefficients = group.get<std::vector<double>>("COEFF");

      return std::make_shared<Core::Utils::FastPolynomialFunction>(std::move(coefficients));
    }
    else if (function_lin_def.has_group("CUBIC_SPLINE_FROM_CSV"))
    {
      const auto& group = function_lin_def.group("CUBIC_SPLINE_FROM_CSV");
      const auto csv_file = group.get<std::filesystem::path>("CSV");

      // safety check
      if (csv_file.empty())
        FOUR_C_THROW("You forgot to specify the *.csv file for cubic spline interpolation!");

      return std::make_shared<Core::Utils::CubicSplineFromCSV>(csv_file.string());
    }
    else
    {
      return {nullptr};
    }
  }
}  // namespace


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::Utils::add_valid_library_functions(Core::Utils::FunctionManager& function_manager)
{
  using namespace IO::InputSpecBuilders;

  auto spec = one_of({
      group("FASTPOLYNOMIAL",
          {
              parameter<int>("NUMCOEFF"),
              parameter<std::vector<double>>("COEFF", {.size = from_parameter<int>("NUMCOEFF")}),
          }),
      group("CUBIC_SPLINE_FROM_CSV",
          {
              parameter<std::filesystem::path>("CSV"),
          }),
  });

  function_manager.add_function_definition(spec, create_library_function_scalar);
}


Core::Utils::FastPolynomialFunction::FastPolynomialFunction(std::vector<double> coefficients)
    : mypoly_(std::move(coefficients))
{
}

double Core::Utils::FastPolynomialFunction::evaluate(const double argument) const
{
  return mypoly_.evaluate(argument);
}

double Core::Utils::FastPolynomialFunction::evaluate_derivative(
    const double argument, const int deriv_order) const
{
  return mypoly_.evaluate_derivative(argument, deriv_order);
}


Core::Utils::CubicSplineFromCSV::~CubicSplineFromCSV() = default;

Core::Utils::CubicSplineFromCSV::CubicSplineFromCSV(const std::string& csv_file)
{
  auto vector_of_csv_columns = Core::IO::read_csv_as_columns(2, csv_file);

  cubic_spline_ = std::make_unique<Core::Utils::CubicSplineInterpolation>(
      Core::Utils::CubicSplineInterpolation(vector_of_csv_columns[0], vector_of_csv_columns[1]));
}

double Core::Utils::CubicSplineFromCSV::evaluate(const double scalar) const
{
  return cubic_spline_->evaluate(scalar);
}

double Core::Utils::CubicSplineFromCSV::evaluate_derivative(
    const double scalar, const int deriv_order) const
{
  return cubic_spline_->evaluate_derivative(scalar, deriv_order);
}

FOUR_C_NAMESPACE_CLOSE

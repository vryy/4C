/*----------------------------------------------------------------------*/
/*! \file

\brief Collection of functions that are not problem-specific

The functions in this file are not problem-specific and may be useful for a number of applications.


\level 3

*/
/*----------------------------------------------------------------------*/
#include "4C_utils_function_library.hpp"

#include "4C_io_control.hpp"
#include "4C_io_file_reader.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_utils_cubic_spline_interpolation.hpp"
#include "4C_utils_function.hpp"
#include "4C_utils_function_manager.hpp"

#include <Teuchos_RCP.hpp>

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace
{

  Teuchos::RCP<Core::UTILS::FunctionOfScalar> create_library_function_scalar(
      const std::vector<Input::LineDefinition>& function_line_defs)
  {
    if (function_line_defs.size() != 1) return Teuchos::null;

    const auto& function_lin_def = function_line_defs.front();

    if (function_lin_def.container().get_or<bool>("FASTPOLYNOMIAL", false))
    {
      std::vector<double> coefficients =
          function_lin_def.container().get<std::vector<double>>("COEFF");

      return Teuchos::make_rcp<Core::UTILS::FastPolynomialFunction>(std::move(coefficients));
    }
    else if (function_lin_def.container().get_or<bool>("CUBIC_SPLINE_FROM_CSV", false))
    {
      const auto csv_file = function_lin_def.container().get<std::filesystem::path>("CSV");

      // safety check
      if (csv_file.empty())
        FOUR_C_THROW("You forgot to specify the *.csv file for cubic spline interpolation!");

      return Teuchos::make_rcp<Core::UTILS::CubicSplineFromCSV>(csv_file.string());
    }
    else
      return {Teuchos::null};
  }
}  // namespace


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::UTILS::add_valid_library_functions(Core::UTILS::FunctionManager& function_manager)
{
  using namespace Input;

  LineDefinition fast_polynomial_funct =
      LineDefinition::Builder()
          .add_tag("FASTPOLYNOMIAL")
          .add_named_int("NUMCOEFF")
          .add_named_double_vector("COEFF", LengthFromIntNamed("NUMCOEFF"))
          .build();

  LineDefinition cubic_spline_from_csv_funct =
      LineDefinition::Builder().add_tag("CUBIC_SPLINE_FROM_CSV").add_named_string("CSV").build();

  function_manager.add_function_definition(
      {std::move(fast_polynomial_funct), std::move(cubic_spline_from_csv_funct)},
      create_library_function_scalar);
}


Core::UTILS::FastPolynomialFunction::FastPolynomialFunction(std::vector<double> coefficients)
    : mypoly_(std::move(coefficients))
{
}

double Core::UTILS::FastPolynomialFunction::evaluate(const double argument) const
{
  return mypoly_.evaluate(argument);
}

double Core::UTILS::FastPolynomialFunction::evaluate_derivative(
    const double argument, const int deriv_order) const
{
  return mypoly_.evaluate_derivative(argument, deriv_order);
}


Core::UTILS::CubicSplineFromCSV::CubicSplineFromCSV(const std::string& csv_file)
{
  auto vector_of_csv_columns = Core::IO::read_csv_as_columns(2, csv_file);

  cubic_spline_ = std::make_unique<Core::UTILS::CubicSplineInterpolation>(
      Core::UTILS::CubicSplineInterpolation(vector_of_csv_columns[0], vector_of_csv_columns[1]));
}


double Core::UTILS::CubicSplineFromCSV::evaluate(const double scalar) const
{
  return cubic_spline_->evaluate(scalar);
}


double Core::UTILS::CubicSplineFromCSV::evaluate_derivative(
    const double scalar, const int deriv_order) const
{
  return cubic_spline_->evaluate_derivative(scalar, deriv_order);
}

FOUR_C_NAMESPACE_CLOSE

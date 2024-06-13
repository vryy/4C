/*----------------------------------------------------------------------*/
/*! \file

\brief Collection of functions that are not problem-specific

The functions in this file are not problem-specific and may be useful for a number of applications.


\level 3

*/
/*----------------------------------------------------------------------*/
#include "4C_utils_function_library.hpp"

#include "4C_global_data.hpp"
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

  Teuchos::RCP<Core::UTILS::FunctionOfScalar> CreateLibraryFunctionScalar(
      const std::vector<Input::LineDefinition>& function_line_defs)
  {
    if (function_line_defs.size() != 1) return Teuchos::null;

    const auto& function_lin_def = function_line_defs.front();

    if (function_lin_def.has_named("FASTPOLYNOMIAL"))
    {
      std::vector<double> coefficients;
      function_lin_def.extract_double_vector("COEFF", coefficients);

      return Teuchos::rcp(new Core::UTILS::FastPolynomialFunction(std::move(coefficients)));
    }
    else if (function_lin_def.has_named("CUBIC_SPLINE_FROM_CSV"))
    {
      std::string csv_file;
      function_lin_def.extract_string("CSV", csv_file);

      // safety check
      if (csv_file.empty())
        FOUR_C_THROW("You forgot to specify the *.csv file for cubic spline interpolation!");

      // csv file needs to be placed in same folder as input file
      std::filesystem::path input_file_path =
          Global::Problem::Instance()->OutputControlFile()->InputFileName();
      const auto csv_file_path = input_file_path.replace_filename(csv_file);

      return Teuchos::rcp(new Core::UTILS::CubicSplineFromCSV(csv_file_path.string()));
    }
    else
      return {Teuchos::null};
  }
}  // namespace


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::UTILS::AddValidLibraryFunctions(Core::UTILS::FunctionManager& function_manager)
{
  using namespace Input;

  LineDefinition fast_polynomial_funct =
      LineDefinition::Builder()
          .add_tag("FASTPOLYNOMIAL")
          .add_named_int("NUMCOEFF")
          .add_named_double_vector("COEFF", LengthFromIntNamed("NUMCOEFF"))
          .Build();

  LineDefinition cubic_spline_from_csv_funct =
      LineDefinition::Builder().add_tag("CUBIC_SPLINE_FROM_CSV").add_named_string("CSV").Build();

  function_manager.add_function_definition(
      {std::move(fast_polynomial_funct), std::move(cubic_spline_from_csv_funct)},
      CreateLibraryFunctionScalar);
}


Core::UTILS::FastPolynomialFunction::FastPolynomialFunction(std::vector<double> coefficients)
    : mypoly_(std::move(coefficients))
{
}

double Core::UTILS::FastPolynomialFunction::Evaluate(const double argument) const
{
  return mypoly_.Evaluate(argument);
}

double Core::UTILS::FastPolynomialFunction::EvaluateDerivative(
    const double argument, const int deriv_order) const
{
  return mypoly_.EvaluateDerivative(argument, deriv_order);
}


Core::UTILS::CubicSplineFromCSV::CubicSplineFromCSV(const std::string& csv_file)
{
  auto vector_of_csv_columns = Core::IO::ReadCsvAsColumns(2, csv_file);

  cubic_spline_ = std::make_unique<Core::UTILS::CubicSplineInterpolation>(
      Core::UTILS::CubicSplineInterpolation(vector_of_csv_columns[0], vector_of_csv_columns[1]));
}


double Core::UTILS::CubicSplineFromCSV::Evaluate(const double scalar) const
{
  return cubic_spline_->Evaluate(scalar);
}


double Core::UTILS::CubicSplineFromCSV::EvaluateDerivative(
    const double scalar, const int deriv_order) const
{
  return cubic_spline_->EvaluateDerivative(scalar, deriv_order);
}

FOUR_C_NAMESPACE_CLOSE

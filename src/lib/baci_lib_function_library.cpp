/*----------------------------------------------------------------------*/
/*! \file

\brief Collection of functions that are not problem-specific

The functions in this file are not problem-specific and may be useful for a number of applications.


\level 3

*/
/*----------------------------------------------------------------------*/
#include "baci_lib_function_library.H"

#include "baci_io_control.H"
#include "baci_io_csv_reader.H"
#include "baci_lib_cubic_spline_interpolation.H"
#include "baci_lib_function.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_linedefinition.H"

#include <Teuchos_RCP.hpp>

#include <filesystem>
#include <utility>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::AddValidLibraryFunctionLines(DRT::INPUT::Lines& lines)
{
  using namespace DRT::INPUT;

  LineDefinition fastpolynomial_funct =
      LineDefinition::Builder()
          .AddTag("FASTPOLYNOMIAL")
          .AddNamedInt("NUMCOEFF")
          .AddNamedDoubleVector("COEFF", LengthFromIntNamed("NUMCOEFF"))
          .Build();

  LineDefinition cubicsplinefromcsv_funct =
      LineDefinition::Builder().AddTag("CUBIC_SPLINE_FROM_CSV").AddNamedString("CSV").Build();

  lines.Add(fastpolynomial_funct);
  lines.Add(cubicsplinefromcsv_funct);
}

Teuchos::RCP<DRT::UTILS::FunctionOfScalar> DRT::UTILS::TryCreateLibraryFunctionScalar(
    const std::vector<DRT::INPUT::LineDefinition>& function_line_defs)
{
  if (function_line_defs.size() != 1) return Teuchos::null;

  const auto& function_lin_def = function_line_defs.front();

  if (function_lin_def.HaveNamed("FASTPOLYNOMIAL"))
  {
    std::vector<double> coefficients;
    function_lin_def.ExtractDoubleVector("COEFF", coefficients);

    return Teuchos::rcp(new FastPolynomialFunction(std::move(coefficients)));
  }
  else if (function_lin_def.HaveNamed("CUBIC_SPLINE_FROM_CSV"))
  {
    std::string csv_file;
    function_lin_def.ExtractString("CSV", csv_file);

    // safety check
    if (csv_file.empty())
      dserror("You forgot to specify the *.csv file for cubic spline interpolation!");

    const std::string input_file = DRT::Problem::Instance()->OutputControlFile()->InputFileName();
    std::filesystem::path input_file_path =
        DRT::Problem::Instance()->OutputControlFile()->InputFileName();
    const auto csv_file_path = input_file_path.replace_filename(csv_file);

    return Teuchos::rcp(new DRT::UTILS::CubicSplineFromCSV(csv_file_path.string()));
  }
  else
    return {Teuchos::null};
}


DRT::UTILS::FastPolynomialFunction::FastPolynomialFunction(std::vector<double> coefficients)
    : mypoly_(std::move(coefficients))
{
}

double DRT::UTILS::FastPolynomialFunction::Evaluate(const double argument) const
{
  return mypoly_.Evaluate(argument);
}

double DRT::UTILS::FastPolynomialFunction::EvaluateDerivative(const double argument) const
{
  return mypoly_.EvaluateDerivative(argument, 1);
}


DRT::UTILS::CubicSplineFromCSV::CubicSplineFromCSV(const std::string& csv_file)
{
  auto vector_of_csv_columns = IO::ReadCsv(2, csv_file);

  cubic_spline_ = std::make_unique<CubicSplineInterpolation>(
      CubicSplineInterpolation(vector_of_csv_columns[0], vector_of_csv_columns[1]));
}


double DRT::UTILS::CubicSplineFromCSV::Evaluate(const double scalar) const
{
  return cubic_spline_->EvaluateScalar(scalar);
}


double DRT::UTILS::CubicSplineFromCSV::EvaluateDerivative(const double scalar) const
{
  return cubic_spline_->EvaluateScalarFirstDerivative(scalar);
}

/*----------------------------------------------------------------------*/
/*! \file

\brief Collection of functions that are not problem-specific

The functions in this file are not problem-specific and may be useful for a number of applications.


\level 3

*/
/*----------------------------------------------------------------------*/
#include <filesystem>
#include <utility>
#include "cubic_spline_interpolation.H"
#include "function.H"
#include "function_library.H"
#include "globalproblem.H"
#include "linedefinition.H"
#include "Teuchos_RCP.hpp"

#include "csv_reader.H"
#include "io_control.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::AddValidLibraryFunctionLines(Teuchos::RCP<DRT::INPUT::Lines> lines)
{
  DRT::INPUT::LineDefinition fastpolynomial_funct;
  fastpolynomial_funct.AddTag("FASTPOLYNOMIAL")
      .AddNamedInt("NUMCOEFF")
      .AddNamedDoubleVector("COEFF", "NUMCOEFF");

  DRT::INPUT::LineDefinition translatedfunction_funct;
  translatedfunction_funct.AddTag("TRANSLATEDFUNCTION").AddNamedInt("ORIGIN").AddNamedInt("LOCAL");

  DRT::INPUT::LineDefinition cubicsplinefromcsv_funct;
  cubicsplinefromcsv_funct.AddTag("CUBIC_SPLINE_FROM_CSV").AddNamedString("CSV");

  lines->Add(translatedfunction_funct);
  lines->Add(fastpolynomial_funct);
  lines->Add(cubicsplinefromcsv_funct);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime> DRT::UTILS::TryCreateLibraryFunction(
    Teuchos::RCP<DRT::INPUT::LineDefinition> function_lin_def, DRT::UTILS::FunctionManager& manager,
    const int index_current_funct_in_manager)
{
  if (function_lin_def->HaveNamed("TRANSLATEDFUNCTION"))
  {
    int origin, local;
    function_lin_def->ExtractInt("ORIGIN", origin);
    function_lin_def->ExtractInt("LOCAL", local);

    if (origin <= 0 or origin >= index_current_funct_in_manager)
    {
      dserror(
          "ORIGIN function ID (currently %d) must be positive and smaller than "
          "TRANSLATEDFUNCTION (currently %d).",
          origin, index_current_funct_in_manager);
    }
    if (local <= 0 or local >= index_current_funct_in_manager)
    {
      dserror(
          "LOCAL function ID (currently %d) must be positive and smaller than "
          "TRANSLATEDFUNCTION (currently %d).",
          local, index_current_funct_in_manager);
    }

    Teuchos::RCP<FunctionOfSpaceTime> origin_funct =
        Teuchos::rcpFromRef(manager.FunctionById<FunctionOfSpaceTime>(origin - 1));
    Teuchos::RCP<FunctionOfSpaceTime> local_funct =
        Teuchos::rcpFromRef(manager.FunctionById<FunctionOfSpaceTime>(local - 1));

    return Teuchos::rcp(new TranslatedFunction(origin_funct, local_funct));
  }
  else
  {
    return {Teuchos::null};
  }
}

Teuchos::RCP<DRT::UTILS::FunctionOfScalar> DRT::UTILS::TryCreateLibraryFunctionScalar(
    Teuchos::RCP<DRT::INPUT::LineDefinition> function_lin_def, DRT::UTILS::FunctionManager& manager,
    const int /*index_current_funct_in_manager*/)
{
  if (function_lin_def->HaveNamed("FASTPOLYNOMIAL"))
  {
    std::vector<double> coefficients;
    function_lin_def->ExtractDoubleVector("COEFF", coefficients);

    return Teuchos::rcp(new FastPolynomialFunction(std::move(coefficients)));
  }
  else if (function_lin_def->HaveNamed("CUBIC_SPLINE_FROM_CSV"))
  {
    std::string csv_file;
    function_lin_def->ExtractString("CSV", csv_file);

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


DRT::UTILS::TranslatedFunction::TranslatedFunction(
    Teuchos::RCP<FunctionOfSpaceTime> origin, Teuchos::RCP<FunctionOfSpaceTime> local)
    : originFunction_(std::move(origin)), localFunction_(std::move(local))
{
  if (originFunction_->NumberComponents() != nsd_originTranslation)
    dserror("Origin function needs to have exactly %d components but %d were given.",
        nsd_originTranslation, originFunction_->NumberComponents());
}

double DRT::UTILS::TranslatedFunction::Evaluate(
    const double* x, const double t, const std::size_t component)
{
  if (component < 0 or component >= localFunction_->NumberComponents())
    dserror("Component must be between 0 and %d but is %d.", localFunction_->NumberComponents(),
        component);

  std::array<double, 3> new_coord{};
  for (int i = 0; i < nsd_originTranslation; i++)
  {
    new_coord[i] = x[i] - originFunction_->Evaluate(x, t, i);
  }
  return localFunction_->Evaluate(new_coord.data(), t, component);
}

std::vector<double> DRT::UTILS::TranslatedFunction::EvaluateTimeDerivative(
    const double* x, const double t, const unsigned deg, const std::size_t component)
{
  if (deg == 0) return std::vector<double>{Evaluate(x, t, component)};

  if (deg != 1) dserror("Time derivative only implemented for degree <= 1");

  if (component < 0 or component >= localFunction_->NumberComponents())
    dserror("Component must be between 0 and %d but is %d.", localFunction_->NumberComponents(),
        component);

  std::vector<double> result(deg + 1);
  result[0] = Evaluate(x, t, component);

  std::array<double, 3> translatedCoord{};
  std::array<double, 3> translatedDeriv{};
  for (int i = 0; i < nsd_originTranslation; i++)
  {
    auto evalResult = originFunction_->EvaluateTimeDerivative(x, t, 1, i);
    translatedCoord[i] = x[i] - evalResult[0];
    translatedDeriv[i] = -evalResult[1];
  }

  auto localSpatialDeriv =
      localFunction_->EvaluateSpatialDerivative(translatedCoord.data(), t, component);
  auto localValues =
      localFunction_->EvaluateTimeDerivative(translatedCoord.data(), t, 1, component);
  // total time derivative according to chain rule
  // dh/dt = -df/dx*dg/dt + df/dt
  result[1] = localSpatialDeriv[0] * translatedDeriv[0] +
              localSpatialDeriv[1] * translatedDeriv[1] +
              localSpatialDeriv[2] * translatedDeriv[2] + localValues[1];
  return result;
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
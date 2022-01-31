/*----------------------------------------------------------------------*/
/*! \file

\brief Collection of functions that are not problem-specific

The functions in this file are not problem-specific and may be useful for a number of applications.


\level 3

*/
/*----------------------------------------------------------------------*/
#include <fstream>
#include <Sacado.hpp>
#include <utility>
#include "drt_function.H"
#include "drt_function_library.H"
#include "drt_globalproblem.H"


DRT::UTILS::FastPolynomialFunction::FastPolynomialFunction(std::vector<double>* coefficients)
    : Function(), mypoly_(new Polynomial(*coefficients))
{
}

/*---------------------------------------------------------------------*
 *---------------------------------------------------------------------*/
double DRT::UTILS::FastPolynomialFunction::Evaluate(const double argument) const
{
  return mypoly_->Evaluate(argument);
}

/*---------------------------------------------------------------------*
 *---------------------------------------------------------------------*/
double DRT::UTILS::FastPolynomialFunction::EvaluateDerivative(const double argument) const
{
  LINALG::Matrix<2, 1> derivs(false);
  mypoly_->Evaluate(argument, derivs);
  return derivs(1);
}

DRT::UTILS::TranslatedFunction::TranslatedFunction(
    Teuchos::RCP<Function> origin, Teuchos::RCP<Function> local)
    : originFunction_(std::move(origin)), localFunction_(std::move(local))
{
  if (originFunction_->NumberComponents() != nsd_originTranslation)
    dserror("Origin function needs to have exactly %d components but %d were given.",
        nsd_originTranslation, originFunction_->NumberComponents());
}

double DRT::UTILS::TranslatedFunction::Evaluate(const int index, const double* x, double t)
{
  if (index < 0 or index >= static_cast<int>(localFunction_->NumberComponents()))
    dserror("Index must be between 0 and %d but is %d.", localFunction_->NumberComponents(), index);

  std::array<double, 3> new_coord{};
  for (int i = 0; i < nsd_originTranslation; i++)
  {
    new_coord[i] = x[i] - originFunction_->Evaluate(i, x, t);
  }
  return localFunction_->Evaluate(index, new_coord.data(), t);
}

std::vector<double> DRT::UTILS::TranslatedFunction::EvaluateTimeDerivative(
    const int index, const double* x, const double t, const unsigned deg)
{
  if (deg == 0) return std::vector<double>{Evaluate(index, x, t)};

  if (deg != 1) dserror("Time derivative only implemented for degree <= 1");

  if (index < 0 or index >= static_cast<int>(localFunction_->NumberComponents()))
    dserror("Index must be between 0 and %d but is %d.", localFunction_->NumberComponents(), index);

  std::vector<double> result(deg + 1);
  result[0] = Evaluate(index, x, t);

  std::array<double, 3> translatedCoord{};
  std::array<double, 3> translatedDeriv{};
  for (int i = 0; i < nsd_originTranslation; i++)
  {
    auto evalResult = originFunction_->EvaluateTimeDerivative(i, x, t, 1);
    translatedCoord[i] = x[i] - evalResult[0];
    translatedDeriv[i] = -evalResult[1];
  }

  auto localSpatialDeriv =
      localFunction_->EvaluateSpatialDerivative(index, translatedCoord.data(), t);
  auto localValues = localFunction_->EvaluateTimeDerivative(index, translatedCoord.data(), t, 1);
  // total time derivative according to chain rule
  // dh/dt = -df/dx*dg/dt + df/dt
  result[1] = localSpatialDeriv[0] * translatedDeriv[0] +
              localSpatialDeriv[1] * translatedDeriv[1] +
              localSpatialDeriv[2] * translatedDeriv[2] + localValues[1];
  return result;
}

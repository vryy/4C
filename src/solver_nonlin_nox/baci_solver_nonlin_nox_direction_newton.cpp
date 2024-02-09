/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#include "baci_solver_nonlin_nox_direction_newton.hpp"

#include "baci_solver_nonlin_nox_group.hpp"

#include <NOX_GlobalData.H>
#include <NOX_Utils.H>

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Direction::Newton::Newton(
    const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& p)
    : ::NOX::Direction::Newton(gd, p), utils_(gd->getUtils())
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Direction::Newton::compute(::NOX::Abstract::Vector& dir,
    ::NOX::Abstract::Group& group, const ::NOX::Solver::Generic& solver)
{
  ::NOX::Abstract::Group::ReturnType status;

  // dynamic cast of the nox_abstract_group
  NOX::NLN::Group* nlnSoln = dynamic_cast<NOX::NLN::Group*>(&group);

  if (nlnSoln == nullptr)
  {
    throwError("compute", "dynamic_cast to nox_nln_group failed!");
  }

  // Compute F and Jacobian at current solution at once.
  status = nlnSoln->computeFandJacobian();
  if (status != ::NOX::Abstract::Group::Ok)
    throwError("compute", "Unable to compute F and/or Jacobian at once");

  // ------------------------------------------------
  // call base class version
  // ------------------------------------------------
  // NOTE: The base class calls first computeF() and in a second attempt
  // computeJacobian(). In many cases this is uneconomical and instead we do it
  // here at once. In this way the base class compute function calls perform a
  // direct return, because the isValid flags are already set to true!
  try
  {
    return ::NOX::Direction::Newton::compute(dir, group, solver);
  }
  catch (const char* e)
  {
    if (utils_->isPrintType(::NOX::Utils::Warning))
    {
      utils_->out() << e;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Direction::Newton::throwError(
    const std::string& functionName, const std::string& errorMsg)
{
  if (utils_->isPrintType(::NOX::Utils::Error))
    utils_->err() << "NOX::NLN::Direction::Newton::" << functionName << " - " << errorMsg
                  << std::endl;
  throw "NOX Error";
}

BACI_NAMESPACE_CLOSE

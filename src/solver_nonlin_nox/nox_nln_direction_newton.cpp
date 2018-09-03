/*-----------------------------------------------------------*/
/*!
\file nox_nln_direction_newton.cpp

\maintainer Michael Hiermeier

\date Aug 6, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_direction_newton.H"
#include "nox_nln_group.H"

#include <NOX_GlobalData.H>
#include <NOX_Utils.H>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Direction::Newton::Newton(
    const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& p)
    : NOX::Direction::Newton(gd, p), utils_(gd->getUtils())
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::Direction::Newton::compute(
    NOX::Abstract::Vector& dir, NOX::Abstract::Group& soln, const NOX::Solver::Generic& solver)
{
  NOX::Abstract::Group::ReturnType status;

  // dynamic cast of the nox_abstract_group
  NOX::NLN::Group* nlnSoln = dynamic_cast<NOX::NLN::Group*>(&soln);

  if (nlnSoln == NULL)
  {
    throwError("compute", "dynamic_cast to nox_nln_group failed!");
  }

  // Compute F and Jacobian at current solution at once.
  status = nlnSoln->computeFandJacobian();
  if (status != NOX::Abstract::Group::Ok)
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
    return NOX::Direction::Newton::compute(dir, soln, solver);
  }
  catch (const char* e)
  {
    if (utils_->isPrintType(NOX::Utils::Warning))
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
  if (utils_->isPrintType(NOX::Utils::Error))
    utils_->err() << "NOX::NLN::Direction::Newton::" << functionName << " - " << errorMsg
                  << std::endl;
  throw "NOX Error";
}

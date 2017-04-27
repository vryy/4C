/*!----------------------------------------------------------------------
\file fsi_nox_fixpoint.cpp

\brief Calculates the fix point direction.

\maintainer Andreas Rauch

\level 1

*----------------------------------------------------------------------*/

#include "fsi_nox_fixpoint.H"
#include <NOX_GlobalData.H>
#include <NOX_Abstract_Group.H>

NOX::FSI::FixPoint::FixPoint(const Teuchos::RCP<NOX::Utils>& utils,
                             Teuchos::ParameterList& params)
  : utils_(utils)
{
}


NOX::FSI::FixPoint::~FixPoint()
{
}


bool NOX::FSI::FixPoint::reset(const Teuchos::RCP<NOX::GlobalData>& gd,
                               Teuchos::ParameterList& params)
{
  utils_ = gd->getUtils();
  return true;
}


bool NOX::FSI::FixPoint::compute(NOX::Abstract::Vector& dir,
                                 NOX::Abstract::Group& soln,
                                 const NOX::Solver::Generic& solver)
{
  NOX::Abstract::Group::ReturnType status;

  // Compute F at current solution
  status = soln.computeF();
  if (status != NOX::Abstract::Group::Ok)
    throwError("compute", "Unable to compute F");

  // The residual is the direction.
  dir.update(1.0, soln.getF());

  return true;
}


bool NOX::FSI::FixPoint::compute(NOX::Abstract::Vector& dir,
                                 NOX::Abstract::Group& soln,
                                 const NOX::Solver::LineSearchBased& solver)
{
  return NOX::Direction::Generic::compute( dir, soln, solver );
}


void NOX::FSI::FixPoint::throwError(const std::string& functionName,
                                    const std::string& errorMsg)
{
    if (utils_->isPrintType(NOX::Utils::Error))
      utils_->err() << "FixPoint::" << functionName
                    << " - " << errorMsg << std::endl;
    throw "NOX Error";
}

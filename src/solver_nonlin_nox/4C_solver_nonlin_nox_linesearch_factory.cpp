/*-----------------------------------------------------------*/
/*! \file

\brief Factory to create the desired Line Search object.



\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_solver_nonlin_nox_linesearch_factory.hpp"  // class definition

#include "4C_solver_nonlin_nox_linesearch_backtrack.hpp"
#include "4C_utils_exceptions.hpp"

#include <NOX_Common.H>
#include <NOX_LineSearch_FullStep.H>
#include <NOX_StatusTest_Generic.H>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::LineSearch::Factory::Factory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::LineSearch::Generic> NOX::Nln::LineSearch::Factory::build_line_search(
    const Teuchos::RCP<::NOX::GlobalData>& gd,
    const Teuchos::RCP<::NOX::StatusTest::Generic> outerTests,
    const Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic> innerTests,
    Teuchos::ParameterList& params) const
{
  Teuchos::RCP<::NOX::LineSearch::Generic> line_search;

  std::string method = params.get("Method", "Full Step");

  // If we use not the full step method, a inner status test has to be provided!
  if (method != "Full Step") inner_status_test_is_required(innerTests);

  if (method == "Full Step")
    line_search = Teuchos::rcp(new ::NOX::LineSearch::FullStep(gd, params));
  else if (method == "Backtrack")
  {
    line_search =
        Teuchos::rcp(new NOX::Nln::LineSearch::Backtrack(gd, outerTests, innerTests, params));
  }
  else
  {
    std::ostringstream msg;
    msg << "Error - NOX::Nln::LineSearch::Factory::build_line_search() - The \"Method\" parameter "
           "\""
        << method << "\" is not a valid linesearch option. " << std::endl
        << "Please fix your parameter list!" << std::endl;
    FOUR_C_THROW(msg.str());
  }

  return line_search;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::LineSearch::Factory::inner_status_test_is_required(
    const Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>& innerTests) const
{
  if (innerTests.is_null())
    FOUR_C_THROW(
        "ERROR - NOX::Nln::LineSearch::Factory::inner_status_test_is_required -"
        " The inner status test pointer should be initialized at this point!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::LineSearch::Generic> NOX::Nln::LineSearch::BuildLineSearch(
    const Teuchos::RCP<::NOX::GlobalData>& gd,
    const Teuchos::RCP<::NOX::StatusTest::Generic> outerTests,
    const Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic> innerTests,
    Teuchos::ParameterList& lsparams)
{
  Factory factory;
  return factory.build_line_search(gd, outerTests, innerTests, lsparams);
}

FOUR_C_NAMESPACE_CLOSE

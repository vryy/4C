/*-----------------------------------------------------------*/
/*! \file

\brief Factory to create the desired Line Search object.



\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_linesearch_factory.H"  // class definition

#include "../drt_lib/drt_dserror.H"

#include <NOX_Common.H>
#include <NOX_StatusTest_Generic.H>

#include <Teuchos_ParameterList.hpp>

// All the different supported line searches
#include "nox_nln_linesearch_backtrack.H"
#include <NOX_LineSearch_FullStep.H>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::LineSearch::Factory::Factory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::LineSearch::Generic> NOX::NLN::LineSearch::Factory::BuildLineSearch(
    const Teuchos::RCP<NOX::GlobalData>& gd,
    const Teuchos::RCP<NOX::StatusTest::Generic> outerTests,
    const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> innerTests,
    Teuchos::ParameterList& lsparams)
{
  Teuchos::RCP<NOX::LineSearch::Generic> line_search;

  std::string method = lsparams.get("Method", "Full Step");

  // If we use not the full step method, a inner status test has to be provided!
  if (method != "Full Step") InnerStatusTestIsRequired(innerTests);

  if (method == "Full Step")
    line_search = Teuchos::rcp(new NOX::LineSearch::FullStep(gd, lsparams));
  else if (method == "Backtrack")
  {
    line_search =
        Teuchos::rcp(new NOX::NLN::LineSearch::Backtrack(gd, outerTests, innerTests, lsparams));
  }
  else
  {
    std::ostringstream msg;
    msg << "Error - NOX::NLN::LineSearch::Factory::BuildLineSearch() - The \"Method\" parameter \""
        << method << "\" is not a valid linesearch option. " << std::endl
        << "Please fix your parameter list!" << std::endl;
    dserror(msg.str());
  }

  return line_search;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::LineSearch::Factory::InnerStatusTestIsRequired(
    const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests) const
{
  if (innerTests.is_null())
    dserror(
        "ERROR - NOX::NLN::LineSearch::Factory::InnerStatusTestIsRequired -"
        " The inner status test pointer should be initialized at this point!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::LineSearch::Generic> NOX::NLN::LineSearch::BuildLineSearch(
    const Teuchos::RCP<NOX::GlobalData>& gd,
    const Teuchos::RCP<NOX::StatusTest::Generic> outerTests,
    const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> innerTests,
    Teuchos::ParameterList& lsparams)
{
  Factory factory;
  return factory.BuildLineSearch(gd, outerTests, innerTests, lsparams);
}

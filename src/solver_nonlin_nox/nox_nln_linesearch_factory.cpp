/*-----------------------------------------------------------*/
/*!
\file nox_nln_linesearch_factory.cpp

\maintainer Michael Hiermeier

\date Jun 11, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_linesearch_factory.H" // class definition
#include "nox_nln_globaldata.H"

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
    const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> innerTests,
    Teuchos::ParameterList& lsparams)
{
  Teuchos::RCP<NOX::LineSearch::Generic> line_search;

  std::string method = lsparams.get("Method", "Full Step");

  if (method == "Full Step")
    line_search = Teuchos::rcp(new NOX::LineSearch::FullStep(gd, lsparams));
  else if (method == "Backtrack")
    line_search = Teuchos::rcp(new NOX::NLN::LineSearch::Backtrack(gd, innerTests, lsparams));
  else
  {
    std::ostringstream msg;
        msg << "Error - NOX::NLN::LineSearch::Factory::BuildLineSearch() - The \"Method\" parameter \""
            << method << "\" is not a valid linesearch option. " << std::endl
            << "Please fix your parameter list!" << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
  }

  return line_search;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::LineSearch::Generic> NOX::NLN::LineSearch::BuildLineSearch(
    const Teuchos::RCP<NOX::GlobalData>& gd,
    const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> innerTests,
    Teuchos::ParameterList& lsparams)
{
  Factory factory;
  return factory.BuildLineSearch(gd, innerTests, lsparams);
}

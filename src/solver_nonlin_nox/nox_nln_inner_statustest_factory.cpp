/*-----------------------------------------------------------*/
/*!
\file nox_nln_inner_statustest_factory.cpp

\maintainer Michael Hiermeier

\date Aug 14, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_inner_statustest_factory.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_ParameterList.hpp>

#include <NOX_Utils.H>

// supported inner status tests
#include "nox_nln_inner_statustest_armijo.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::Factory::Factory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>
NOX::NLN::INNER::StatusTest::Factory::BuildInnerStatusTests(
    Teuchos::ParameterList& p,
    const NOX::Utils& u,
    std::map<std::string, Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> >* tagged_tests) const
{
  Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> status_test;

  std::string test_type = "???";

  if (Teuchos::isParameterType<std::string>(p, "Test Type"))
    test_type = Teuchos::get<std::string>(p, "Test Type");
  else
  {
    dserror("Error - The \"Test Type\" is a required parameter in the NOX::NLN::StatusTest::Factory!");
  }

  if (test_type == "Combo")
    status_test =  this->BuildComboTest(p, u, tagged_tests);
  else if (test_type == "Armijo")
    status_test = this->BuildArmijoTest(p, u);
  // supported StatusTests of the NOX::StatusTest classes for the inner check
  else if (test_type == "Stagnation" or
           test_type == "Divergence" or
           test_type == "MaxIters"   or
           test_type == "FiniteValue")
  {
    throwError("BuildInnerStatusTests()","Not yet supported");
  }
  else
  {
    std::ostringstream msg;
    msg << "The test type \"" << test_type << "\" is invalid!";
    throwError("BuildInnerStatusTests()",msg.str());
  }

  this->CheckAndTagTest(p, status_test, tagged_tests);

  return status_test;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>
NOX::NLN::INNER::StatusTest::Factory::BuildArmijoTest(
    Teuchos::ParameterList& p,
    const NOX::Utils& u) const
{
  // ------------------------------------------
  // Get c_1
  // Slope scaling parameter. Typical values
  // are 1.0e-4 for a Newton search direction
  // and 0.9 for a steepest descent or
  // BFGS-direction.
  // Anyway, for the latter cases I recommend
  // to use the (strong) Wolfe conditions,
  // or the Goldstein rule.
  // ------------------------------------------
  double c_1 = p.get("c_1",1.0e-4);

  // ------------------------------------------
  // Switch monotone behavior on and off
  // ------------------------------------------
  bool isMonotone = p.get("Monotone",true);

  // ------------------------------------------
  // Get the maximal length of the history
  // vector for a non-monotone Armijo rule
  // ------------------------------------------
  std::size_t maxHistSize = static_cast<std::size_t>(p.get<int>("Maximal History Length",1));

  Teuchos::RCP<NOX::NLN::INNER::StatusTest::Armijo> status_test = Teuchos::null;
  if (isMonotone)
    status_test =
        Teuchos::rcp(new NOX::NLN::INNER::StatusTest::Armijo(c_1));
  else
    status_test =
        Teuchos::rcp(new NOX::NLN::INNER::StatusTest::Armijo(c_1,isMonotone,maxHistSize));

  return status_test;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>
NOX::NLN::INNER::StatusTest::Factory::BuildComboTest(
    Teuchos::ParameterList& p,
    const NOX::Utils& u,
    std::map<std::string, Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> >* tagged_tests) const
{
  throwError("BuildComboTest()","Not yet implemented!");

  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::INNER::StatusTest::Factory::CheckAndTagTest(
    const Teuchos::ParameterList& p,
    const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& test,
    std::map<std::string, Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> >*
    tagged_tests) const
{
  if ( (Teuchos::isParameterType<std::string>(p, "Tag")) && (tagged_tests != NULL) ) {
    (*tagged_tests)[Teuchos::getParameter<std::string>(p, "Tag")] = test;
    return true;
  }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Factory::throwError(
    const std::string& functionName,
    const std::string& errorMsg) const
{
  std::ostringstream msg;
  msg << "ERROR - NOX::NLN::INNER::StatusTest::Factory::" << functionName
      << " - " << errorMsg << std::endl;
  dserror(msg.str());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>
NOX::NLN::INNER::StatusTest::BuildInnerStatusTests(
    Teuchos::ParameterList& p,
    const NOX::Utils& u,
    std::map<std::string, Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> >* tagged_tests)
{
  Factory factory;
  return factory.BuildInnerStatusTests(p,u,tagged_tests);
}

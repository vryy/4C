#include "4C_solver_nonlin_nox_inner_statustest_factory.hpp"

#include "4C_solver_nonlin_nox_enum_lists.hpp"
#include "4C_solver_nonlin_nox_inner_statustest_armijo.hpp"
#include "4C_solver_nonlin_nox_inner_statustest_upperbound.hpp"
#include "4C_solver_nonlin_nox_params_utils.hpp"
#include "4C_utils_exceptions.hpp"

#include <NOX_Abstract_Vector.H>
#include <NOX_Utils.H>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Inner::StatusTest::Factory::Factory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>
NOX::Nln::Inner::StatusTest::Factory::build_inner_status_tests(Teuchos::ParameterList& p,
    const ::NOX::Utils& u,
    std::map<std::string, Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>>* tagged_tests) const
{
  Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic> status_test;

  std::string test_type = "???";

  if (Teuchos::isParameterType<std::string>(p, "Test Type"))
    test_type = Teuchos::get<std::string>(p, "Test Type");
  else
  {
    FOUR_C_THROW(
        "Error - The \"Test Type\" is a required parameter in the NOX::Nln::StatusTest::Factory!");
  }

  if (test_type == "Armijo")
    status_test = this->build_armijo_test(p, u);
  else if (test_type == "UpperBound")
    status_test = this->build_upper_bound_test(p, u);
  // supported StatusTests of the NOX::StatusTest classes for the inner check
  else if (test_type == "Stagnation" or test_type == "Divergence" or test_type == "MaxIters" or
           test_type == "FiniteValue")
  {
    throw_error("build_inner_status_tests()", "Not yet supported");
  }
  else
  {
    std::ostringstream msg;
    msg << "The test type \"" << test_type << "\" is invalid!";
    throw_error("build_inner_status_tests()", msg.str());
  }

  this->check_and_tag_test(p, status_test, tagged_tests);

  return status_test;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>
NOX::Nln::Inner::StatusTest::Factory::build_armijo_test(
    Teuchos::ParameterList& p, const ::NOX::Utils& u) const
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
  const double c_1 = p.get<double>("c_1", 1.0e-4);

  // ------------------------------------------
  // Switch monotone behavior on and off
  // ------------------------------------------
  const bool isMonotone = p.get<bool>("Monotone", true);

  // ------------------------------------------
  // Get the maximal length of the history
  // vector for a non-monotone Armijo rule
  // ------------------------------------------
  std::size_t maxHistSize = static_cast<unsigned>(p.get<int>("Maximal History Length", 1));

  Teuchos::RCP<NOX::Nln::Inner::StatusTest::Armijo> status_test = Teuchos::null;
  if (isMonotone)
    status_test = Teuchos::make_rcp<NOX::Nln::Inner::StatusTest::Armijo>(c_1);
  else
    status_test =
        Teuchos::make_rcp<NOX::Nln::Inner::StatusTest::Armijo>(c_1, isMonotone, maxHistSize);

  return status_test;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>
NOX::Nln::Inner::StatusTest::Factory::build_upper_bound_test(
    Teuchos::ParameterList& p, const ::NOX::Utils& u) const
{
  // Get upper bound as specified in xml file
  double upperboundval = p.get("Value", 1.0e10);
  if (upperboundval < 1e-14)
  {
    std::string msg =
        "\"Value\" for Inner Status Test \"UpperBound\" must be "
        "greater than zero!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }

  const std::string& quantity_type_string = p.get("Quantity Type", "???");
  const NOX::Nln::StatusTest::QuantityType qtype =
      NOX::Nln::StatusTest::string_to_quantity_type(quantity_type_string);
  if (qtype == NOX::Nln::StatusTest::quantity_unknown)
  {
    std::ostringstream msg;
    FOUR_C_THROW(
        "The \"Quantity Type\" is a required parameter "
        "and the chosen option \"%s\" is invalid!",
        quantity_type_string.c_str());
  }

  // Norm Type
  std::string norm_type_string = p.get("Norm Type", "Two Norm");
  ::NOX::Abstract::Vector::NormType norm_type = ::NOX::Abstract::Vector::TwoNorm;
  if (norm_type_string == "Two Norm")
    norm_type = ::NOX::Abstract::Vector::TwoNorm;
  else if (norm_type_string == "One Norm")
    norm_type = ::NOX::Abstract::Vector::OneNorm;
  else if (norm_type_string == "Max Norm")
    norm_type = ::NOX::Abstract::Vector::MaxNorm;
  else
  {
    std::string msg = "\"Norm Type\" must be either \"Two Norm\", \"One Norm\", or \"Max Norm\"!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }

  Teuchos::RCP<NOX::Nln::Inner::StatusTest::UpperBound> status_test =
      Teuchos::make_rcp<NOX::Nln::Inner::StatusTest::UpperBound>(upperboundval, norm_type, qtype);

  return status_test;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::Nln::Inner::StatusTest::Factory::check_and_tag_test(const Teuchos::ParameterList& p,
    const Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>& test,
    std::map<std::string, Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>>* tagged_tests) const
{
  if ((Teuchos::isParameterType<std::string>(p, "Tag")) && (tagged_tests != nullptr))
  {
    (*tagged_tests)[Teuchos::getParameter<std::string>(p, "Tag")] = test;
    return true;
  }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Inner::StatusTest::Factory::throw_error(
    const std::string& functionName, const std::string& errorMsg) const
{
  std::ostringstream msg;
  msg << "ERROR - NOX::Nln::Inner::StatusTest::Factory::" << functionName << " - " << errorMsg
      << std::endl;
  FOUR_C_THROW(msg.str());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>
NOX::Nln::Inner::StatusTest::build_inner_status_tests(Teuchos::ParameterList& p,
    const ::NOX::Utils& u,
    std::map<std::string, Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>>* tagged_tests)
{
  Factory factory;
  return factory.build_inner_status_tests(p, u, tagged_tests);
}

FOUR_C_NAMESPACE_CLOSE

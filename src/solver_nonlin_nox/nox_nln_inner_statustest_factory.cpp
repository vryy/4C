/*-----------------------------------------------------------*/
/*!
\file nox_nln_inner_statustest_factory.cpp

\brief factory for user defined NOX inner status tests

\maintainer Michael Hiermeier

\date Aug 14, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_inner_statustest_factory.H"
#include "nox_nln_params_utils.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_ParameterList.hpp>

#include <NOX_Utils.H>
#include <NOX_Abstract_Vector.H>
#include <Teuchos_StandardParameterEntryValidators.hpp>

// supported inner status tests
#include "nox_nln_inner_statustest_armijo.H"
#include "nox_nln_inner_statustest_filter.H"
#include "nox_nln_inner_statustest_upperbound.H"
#include "nox_nln_inner_statustest_combo.H"
#include "nox_nln_inner_statustest_volume_change.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::Factory::Factory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>
NOX::NLN::INNER::StatusTest::Factory::BuildInnerStatusTests(Teuchos::ParameterList& p,
    const NOX::Utils& u,
    std::map<std::string, Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>>* tagged_tests) const
{
  Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> status_test;

  std::string test_type = "???";

  if (Teuchos::isParameterType<std::string>(p, "Test Type"))
    test_type = Teuchos::get<std::string>(p, "Test Type");
  else
  {
    dserror(
        "Error - The \"Test Type\" is a required parameter in the NOX::NLN::StatusTest::Factory!");
  }

  if (test_type == "Combo")
    status_test = this->BuildComboTest(p, u, tagged_tests);
  else if (test_type == "Armijo")
    status_test = this->BuildArmijoTest(p, u);
  else if (test_type == "Filter")
    status_test = this->BuildFilterTest(p, u);
  else if (test_type == "UpperBound")
    status_test = this->BuildUpperBoundTest(p, u);
  else if (test_type == "VolumeChange")
    status_test = this->BuildVolumeChangeTest(p, u);
  // supported StatusTests of the NOX::StatusTest classes for the inner check
  else if (test_type == "Stagnation" or test_type == "Divergence" or test_type == "MaxIters" or
           test_type == "FiniteValue")
  {
    throwError("BuildInnerStatusTests()", "Not yet supported");
  }
  else
  {
    std::ostringstream msg;
    msg << "The test type \"" << test_type << "\" is invalid!";
    throwError("BuildInnerStatusTests()", msg.str());
  }

  this->CheckAndTagTest(p, status_test, tagged_tests);

  return status_test;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>
NOX::NLN::INNER::StatusTest::Factory::BuildArmijoTest(
    Teuchos::ParameterList& p, const NOX::Utils& u) const
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

  Teuchos::RCP<NOX::NLN::INNER::StatusTest::Armijo> status_test = Teuchos::null;
  if (isMonotone)
    status_test = Teuchos::rcp(new NOX::NLN::INNER::StatusTest::Armijo(c_1));
  else
    status_test =
        Teuchos::rcp(new NOX::NLN::INNER::StatusTest::Armijo(c_1, isMonotone, maxHistSize));

  return status_test;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>
NOX::NLN::INNER::StatusTest::Factory::BuildFilterTest(
    Teuchos::ParameterList& p, const NOX::Utils& u) const
{
  FilterParams fparams;

  std::vector<Teuchos::RCP<NOX::MeritFunction::Generic>>& infeasibility_vec =
      fparams.infeasibility_vec_;

  unsigned count = 0;
  std::ostringstream infeasibility_count;

  infeasibility_count << "Infeasibility Function " << count;
  while (p.isSublist(infeasibility_count.str()))
  {
    const Teuchos::ParameterList& p_infeasibility = p.sublist(infeasibility_count.str());

    // build infeasibility function object
    infeasibility_vec.push_back(
        Teuchos::rcp(new NOX::NLN::MeritFunction::Infeasibility(p_infeasibility, u)));

    /// clear and increase infeasibility function count string
    infeasibility_count.str("");
    infeasibility_count << "Infeasibility Function " << ++count;
  }

  fparams.weight_objective_func_ = p.get<double>("Objective Function Weight", 1.0);
  fparams.weight_infeasibility_func_ = p.get<double>("Infeasibility Function Weight", 1.0);

  fparams.sf_ = p.get<double>("Ftype Condition Exponent s_f", 2.3);
  fparams.st_ = p.get<double>("Ftype Condition Exponent s_t", 1.1);

  /* Safety factor gamma_alpha to compensate for the neglected higher order
   * terms in the chosen model equation during the minimal step length
   * approximation. */
  fparams.gamma_alpha_ = p.get<double>("Gamma Alpha", 0.05);

  // read second order correction parameters
  fparams.use_soc_ = p.get<bool>("Second Order Correction", true);
  fparams.soc_type_ = NOX::NLN::CorrectionType::vague;
  if (fparams.use_soc_)
  {
    // setup validator
    fparams.soc_type_ =
        NOX::NLN::PUTILS::SetAndValidate<NOX::NLN::CorrectionType>(p, "SOC Type", "auto",
            "Second order correction type. Per default the "
            "SOC type is set to \"auto\".",
            Teuchos::tuple<std::string>("auto", "cheap", "full"),
            Teuchos::tuple<NOX::NLN::CorrectionType>(NOX::NLN::CorrectionType::soc_automatic,
                NOX::NLN::CorrectionType::soc_cheap, NOX::NLN::CorrectionType::soc_full));
  }

  // build internal Armijo test
  fparams.armijo_ = BuildArmijoTest(p.sublist("Filter-Armijo"), u);

  // blocking parameters
  {
    int tmp = p.get<int>("Consecutive Blocking Iterates", 3);
    if (tmp < 1) dserror("The Consecutive Blocking Iterates must be greater or equal to 1.");
    fparams.consecutive_blocking_iterates_ = tmp;

    tmp = p.get<int>("Consecutive Blocking LS Steps", 7);
    if (tmp < 1) dserror("The Consecutive Blocking LS Steps must be greater or equal to 1.");
    fparams.consecutive_blocking_ls_steps_ = tmp;

    fparams.max_theta_blocking_red_ = p.get<double>("Max Theta Blocking Reduction", 0.25);
    if (fparams.max_theta_blocking_red_ > 1.0 or fparams.max_theta_blocking_red_ < 0.0)
      dserror("The Max Theta Blocking Reduction must be between 0.0 and 1.0.");

    fparams.init_max_theta_blocking_scaling_ = p.get<double>("Init Max Theta Blocking Scale", 2.0);
    if (fparams.init_max_theta_blocking_scaling_ < 0.0)
      dserror("The Initial Max Theta Blocking Scaling must be greater than 0.0.");
  }

  Teuchos::RCP<NOX::NLN::INNER::StatusTest::Filter> status_test_ptr =
      Teuchos::rcp(new NOX::NLN::INNER::StatusTest::Filter(fparams, u));

  return status_test_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>
NOX::NLN::INNER::StatusTest::Factory::BuildUpperBoundTest(
    Teuchos::ParameterList& p, const NOX::Utils& u) const
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
  const NOX::NLN::StatusTest::QuantityType qtype =
      NOX::NLN::StatusTest::String2QuantityType(quantity_type_string);
  if (qtype == NOX::NLN::StatusTest::quantity_unknown)
  {
    std::ostringstream msg;
    dserror(
        "The \"Quantity Type\" is a required parameter "
        "and the chosen option \"%s\" is invalid!",
        quantity_type_string.c_str());
  }

  // Norm Type
  std::string norm_type_string = p.get("Norm Type", "Two Norm");
  NOX::Abstract::Vector::NormType norm_type = NOX::Abstract::Vector::TwoNorm;
  if (norm_type_string == "Two Norm")
    norm_type = NOX::Abstract::Vector::TwoNorm;
  else if (norm_type_string == "One Norm")
    norm_type = NOX::Abstract::Vector::OneNorm;
  else if (norm_type_string == "Max Norm")
    norm_type = NOX::Abstract::Vector::MaxNorm;
  else
  {
    std::string msg = "\"Norm Type\" must be either \"Two Norm\", \"One Norm\", or \"Max Norm\"!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }

  Teuchos::RCP<NOX::NLN::INNER::StatusTest::UpperBound> status_test =
      Teuchos::rcp(new NOX::NLN::INNER::StatusTest::UpperBound(upperboundval, norm_type, qtype));

  return status_test;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>
NOX::NLN::INNER::StatusTest::Factory::BuildComboTest(Teuchos::ParameterList& p, const NOX::Utils& u,
    std::map<std::string, Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>>* tagged_tests) const
{
  if (not p.isParameter("Combo Type")) dserror("You have to specify a \"Combo Type\".");


  NOX::StatusTest::Combo::ComboType combo_type =
      NOX::NLN::PUTILS::SetAndValidate<NOX::StatusTest::Combo::ComboType>(p, "Combo Type", "AND",
          "Combination type to combine multiple inner "
          "status tests. Possible choices are AND and OR.",
          Teuchos::tuple<std::string>("AND", "OR"),
          Teuchos::tuple<NOX::StatusTest::Combo::ComboType>(
              NOX::StatusTest::Combo::AND, NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::NLN::INNER::StatusTest::Combo> combo_test =
      Teuchos::rcp(new NOX::NLN::INNER::StatusTest::Combo(combo_type, &u));

  int i = 0;
  std::ostringstream subtest_name;
  subtest_name << "Test " << i;
  while (p.isSublist(subtest_name.str()))
  {
    Teuchos::ParameterList& subtest_list = p.sublist(subtest_name.str(), true);

    Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> subtest =
        this->BuildInnerStatusTests(subtest_list, u, tagged_tests);

    combo_test->addStatusTest(subtest);

    // increase iterator
    ++i;
    subtest_name.str("");
    subtest_name << "Test " << i;
  }

  return combo_test;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::INNER::StatusTest::Factory::CheckAndTagTest(const Teuchos::ParameterList& p,
    const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& test,
    std::map<std::string, Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>>* tagged_tests) const
{
  if ((Teuchos::isParameterType<std::string>(p, "Tag")) && (tagged_tests != NULL))
  {
    (*tagged_tests)[Teuchos::getParameter<std::string>(p, "Tag")] = test;
    return true;
  }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Factory::throwError(
    const std::string& functionName, const std::string& errorMsg) const
{
  std::ostringstream msg;
  msg << "ERROR - NOX::NLN::INNER::StatusTest::Factory::" << functionName << " - " << errorMsg
      << std::endl;
  dserror(msg.str());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>
NOX::NLN::INNER::StatusTest::BuildInnerStatusTests(Teuchos::ParameterList& p, const NOX::Utils& u,
    std::map<std::string, Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>>* tagged_tests)
{
  Factory factory;
  return factory.BuildInnerStatusTests(p, u, tagged_tests);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>
NOX::NLN::INNER::StatusTest::Factory::BuildVolumeChangeTest(
    Teuchos::ParameterList& p, const NOX::Utils& u) const
{
  VolumeChangeParams vcparams;
  vcparams.upper_bound_ = p.get<double>("Upper Bound", 2.0);
  vcparams.lower_bound_ = p.get<double>("Lower Bound", 0.5);

  return Teuchos::rcp(new VolumeChange(vcparams, u));
}

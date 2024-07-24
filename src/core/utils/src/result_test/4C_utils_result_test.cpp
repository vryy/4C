/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of general result test framework

\level 1


*/
/*----------------------------------------------------------------------*/

#include "4C_utils_result_test.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_io_pstream.hpp"
#include "4C_utils_exceptions.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN


Core::UTILS::ResultTest::ResultTest(std::string name) : myname_(std::move(name)) {}

void Core::UTILS::ResultTest::test_element(
    const Core::IO::InputParameterContainer& container, int& nerr, int& test_count)
{
  FOUR_C_THROW("no element test available");
}

void Core::UTILS::ResultTest::test_node(
    const Core::IO::InputParameterContainer& container, int& nerr, int& test_count)
{
  FOUR_C_THROW("no node test available");
}

void Core::UTILS::ResultTest::test_node_on_geometry(
    const Core::IO::InputParameterContainer& container, int& nerr, int& test_count,
    const std::vector<std::vector<std::vector<int>>>& nodeset)
{
  FOUR_C_THROW("no geometry test available");
}

void Core::UTILS::ResultTest::test_special(
    const Core::IO::InputParameterContainer& container, int& nerr, int& test_count)
{
  FOUR_C_THROW("no special case test available");
}

void Core::UTILS::ResultTest::test_special(const Core::IO::InputParameterContainer& container,
    int& nerr, int& test_count, int& unevaluated_test_count)
{
  test_special(container, nerr, test_count);
}


int Core::UTILS::ResultTest::compare_values(
    double actresult, std::string type, const Core::IO::InputParameterContainer& container)
{
  int gid;

  if (type != "SPECIAL")
  {
    gid = container.get<int>(type);
  }
  std::string quantity = container.get<std::string>("QUANTITY");
  double givenresult = container.get<double>("VALUE");
  double tolerance = container.get<double>("TOLERANCE");
  // safety check
  if (tolerance <= 0.) FOUR_C_THROW("Tolerance for result test must be strictly positive!");
  // name is an optional input argument!
  auto name = container.get_or<std::string>("NAME", "");

  // return value (0 if results are correct, 1 if results are not correct)
  int ret = 0;

  // prepare std::string stream 'msghead' containing general information on the current test
  std::stringstream msghead;
  msghead << std::left << std::setw(9) << myname_ << ": " << std::left << std::setw(8)
          << quantity.c_str();

  if (name != "") msghead << "(" << name << ")";

  if (type != "SPECIAL")
  {
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    msghead << " at " << type << " " << std::right << std::setw(3) << gid;
  }
  else
  {
    msghead << "\t";
  }

  // write something to screen depending if the result check was ok or not
  if (std::isnan(actresult))
  {
    // Result is 'not a number'
    std::cout << msghead.str() << "\t is NAN!\n";
    ret = 1;
  }
  else if (std::abs(actresult - givenresult) > tolerance)
  {
    // Result is wrong
    std::cout << msghead.str() << "\t is WRONG --> actresult=" << std::setw(24)
              << std::setprecision(17) << std::scientific << actresult
              << ", givenresult=" << std::setw(24) << givenresult << ", abs(diff)=" << std::setw(24)
              << std::abs(actresult - givenresult) << " >" << std::setw(24) << tolerance << "\n";

    ret = 1;
  }
  else
  {
    // Result is correct
    std::cout << msghead.str() << "\t is CORRECT"
              << ", abs(diff)=" << std::setw(24) << std::setprecision(17) << std::scientific
              << std::abs(actresult - givenresult) << " <" << std::setw(24) << tolerance << "\n";
  }

  return ret;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Core::UTILS::ResultTest::match(const Core::IO::InputParameterContainer& container)
{
  return container.get_or<bool>(myname_, false);
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::UTILS::ResultTestManager::add_field_test(Teuchos::RCP<ResultTest> test)
{
  fieldtest_.push_back(test);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::UTILS::ResultTestManager::test_all(const Epetra_Comm& comm)
{
  int nerr = 0;                      // number of tests with errors
  int test_count = 0;                // number of tests performed
  int uneval_test_count = 0;         // number of unevaluated tests
  const int size = results_.size();  // total number of tests

  if (comm.MyPID() == 0) Core::IO::cout << "\nChecking results of " << size << " tests:\n";

  for (auto& result : results_)
  {
    for (const auto& fieldtest : fieldtest_)
    {
      if (fieldtest->match(result.container()))
      {
        if (result.container().get_if<int>("ELEMENT") != nullptr)
          fieldtest->test_element(result.container(), nerr, test_count);
        else if (result.container().get_if<int>("NODE") != nullptr)
          fieldtest->test_node(result.container(), nerr, test_count);
        else if (result.container().get_if<int>("LINE") != nullptr ||
                 result.container().get_if<int>("SURFACE") != nullptr ||
                 result.container().get_if<int>("VOLUME") != nullptr)
          fieldtest->test_node_on_geometry(result.container(), nerr, test_count, get_node_set());
        else
          fieldtest->test_special(result.container(), nerr, test_count, uneval_test_count);
      }
    }
  }

  // print number of unevaluated tests to screen
  int guneval_test_count = 0;
  comm.SumAll(&uneval_test_count, &guneval_test_count, 1);
  if (guneval_test_count > 0 and comm.MyPID() == 0)
    Core::IO::cout << guneval_test_count << " tests stay unevaluated" << Core::IO::endl;

  // determine the total number of errors
  int numerr;
  comm.SumAll(&nerr, &numerr, 1);

  if (numerr > 0)
  {
    FOUR_C_THROW("Result check failed with %d errors out of %d tests", numerr, size);
  }

  /* test_count == -1 means we had a special test routine. It's thus
   * illegal to use both a special routine and single tests. But who
   * wants that? */
  int count;
  if (test_count > -1)
  {
    int lcount = test_count + uneval_test_count;
    comm.SumAll(&lcount, &count, 1);

    /* It's indeed possible to count more tests than expected if you
     * happen to test values of a boundary element. We don't want this
     * FOUR_C_THROW to go off in that case. */
    if (count < size)
    {
      FOUR_C_THROW("expected %d tests but performed %d", size, count);
    }
  }

  if (comm.MyPID() == 0) Core::IO::cout << "\nOK (" << count << ")\n";
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::UTILS::ResultTestManager::set_parsed_lines(std::vector<Input::LineDefinition> results)
{
  results_ = std::move(results);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::UTILS::ResultTestManager::set_node_set(
    const std::vector<std::vector<std::vector<int>>>& nodeset)
{
  nodeset_ = std::move(nodeset);
}

FOUR_C_NAMESPACE_CLOSE

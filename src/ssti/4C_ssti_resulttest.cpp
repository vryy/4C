/*----------------------------------------------------------------------*/
/*! \file
\brief result testing functionality for scalar-structure-thermo interaction problems

\level 2


*/
/*----------------------------------------------------------------------*/
#include "4C_ssti_resulttest.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_ssti_algorithm.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSTI::SSTIResultTest::SSTIResultTest(const SSTI::SSTIAlgorithm& ssti_algorithm)
    : Core::UTILS::ResultTest("SSTI"), ssti_algorithm_(ssti_algorithm)
{
}

/*-------------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------------*/
void SSTI::SSTIResultTest::test_special(Input::LineDefinition& res, int& nerr, int& test_count)
{
  // make sure that quantity is tested only by one processor
  if (ssti_algorithm_.get_comm().MyPID() == 0)
  {
    // extract name of quantity to be tested
    std::string quantity;
    res.extract_string("QUANTITY", quantity);

    // get result to be tested
    const double result = result_special(quantity);

    // compare values
    const int err = compare_values(result, "SPECIAL", res);
    nerr += err;
    ++test_count;
  }
}

/*-------------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------------*/
double SSTI::SSTIResultTest::result_special(const std::string& quantity) const
{
  double result(0.0);

  // number of Newton-Raphson iterations (monolithic SSTI) in last time step
  if (quantity == "numiterlastnonlinearsolve")
  {
    result = static_cast<double>(ssti_algorithm_.iter());
  }

  // test total number of time steps
  else if (!quantity.compare(0, 7, "numstep"))
  {
    result = static_cast<double>(ssti_algorithm_.step());
  }
  // catch unknown quantity strings
  else
  {
    FOUR_C_THROW(
        "Quantity '%s' not supported by result testing functionality for scalar-structure-thermo "
        "interaction!",
        quantity.c_str());
  }

  return result;
}
FOUR_C_NAMESPACE_CLOSE

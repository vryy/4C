/*----------------------------------------------------------------------*/
/*! \file
\brief result testing functionality for scalar-structure-thermo interaction problems

\level 2


*/
/*----------------------------------------------------------------------*/
#include "ssti_resulttest.H"

#include "ssti_algorithm.H"

#include "drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSTI::SSTIResultTest::SSTIResultTest(const SSTI::SSTIAlgorithm& ssti_algorithm)
    : DRT::ResultTest("SSTI"), ssti_algorithm_(ssti_algorithm)
{
}

/*-------------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------------*/
void SSTI::SSTIResultTest::TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // make sure that quantity is tested only by one processor
  if (ssti_algorithm_.Comm().MyPID() == 0)
  {
    // extract name of quantity to be tested
    std::string quantity;
    res.ExtractString("QUANTITY", quantity);

    // get result to be tested
    const double result = ResultSpecial(quantity);

    // compare values
    const int err = CompareValues(result, "SPECIAL", res);
    nerr += err;
    ++test_count;
  }
}

/*-------------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------------*/
double SSTI::SSTIResultTest::ResultSpecial(const std::string& quantity) const
{
  double result(0.0);

  // number of Newton-Raphson iterations (monolithic SSTI) in last time step
  if (quantity == "numiterlastnonlinearsolve")
  {
    result = static_cast<double>(ssti_algorithm_.Iter());
  }

  // test total number of time steps
  else if (!quantity.compare(0, 7, "numstep"))
  {
    result = static_cast<double>(ssti_algorithm_.Step());
  }
  // catch unknown quantity strings
  else
  {
    dserror(
        "Quantity '%s' not supported by result testing functionality for scalar-structure-thermo "
        "interaction!",
        quantity.c_str());
  }

  return result;
}
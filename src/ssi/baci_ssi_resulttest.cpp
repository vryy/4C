/*----------------------------------------------------------------------*/
/*! \file
\brief result testing functionality for scalar-structure interaction problems

\level 2


*/
/*----------------------------------------------------------------------*/
#include "baci_ssi_resulttest.H"

#include "baci_adapter_scatra_base_algorithm.H"
#include "baci_adapter_str_ssiwrapper.H"
#include "baci_lib_linedefinition.H"
#include "baci_linear_solver_method_linalg.H"
#include "baci_scatra_timint_implicit.H"
#include "baci_ssi_monolithic.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::SSIResultTest::SSIResultTest(const Teuchos::RCP<const SSI::SSIBase> ssi_base)
    : DRT::ResultTest("SSI"), ssi_base_(ssi_base)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double SSI::SSIResultTest::ResultSpecial(const std::string& quantity) const
{
  // initialize variable for result
  double result(0.0);

  // number of outer coupling iterations (partitioned SSI) or Newton-Raphson iterations (monolithic
  // SSI) in last time step
  if (quantity == "numiterlastnonlinearsolve")
  {
    result = static_cast<double>(ssi_base_->IterationCount());
  }

  // number of iterations performed by linear solver during last Newton-Raphson iteration
  // (monolithic SSI only)
  else if (quantity == "numiterlastlinearsolve")
  {
    result = static_cast<double>(SSIMono().Solver().getNumIters());
  }

  // test total number of time steps
  else if (!quantity.compare(0, 7, "numstep"))
  {
    result = static_cast<double>(ssi_base_->Step());
  }
  // catch unknown quantity strings
  else
  {
    dserror(
        "Quantity '%s' not supported by result testing functionality for scalar-structure "
        "interaction!",
        quantity.c_str());
  }

  return result;
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
const SSI::SSIMono& SSI::SSIResultTest::SSIMono() const
{
  const auto* const ssi_mono = dynamic_cast<const SSI::SSIMono* const>(ssi_base_.get());
  if (ssi_mono == nullptr)
    dserror("Couldn't access time integrator for monolithic scalar-structure interaction!");
  return *ssi_mono;
}

/*-------------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------------*/
void SSI::SSIResultTest::TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // make sure that quantity is tested only by one processor
  if (ssi_base_->Comm().MyPID() == 0)
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

/*----------------------------------------------------------------------*/
/*!
\file sti_resulttest.cpp

\brief result testing functionality for scatra-thermo interaction problems

\level 2

\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
*/
/*----------------------------------------------------------------------*/
#include "sti_resulttest.H"

#include "sti_monolithic.H"

#include "../drt_lib/drt_linedefinition.H"

#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 01/17 |
 *----------------------------------------------------------------------*/
STI::STIResultTest::STIResultTest(const Teuchos::RCP<STI::Algorithm>&
        sti_algorithm  //!< time integrator for scatra-thermo interaction
    )
    // call base class constructor
    : DRT::ResultTest("STI"),

      // store pointer to time integrator for scatra-thermo interaction
      sti_algorithm_(sti_algorithm)
{
  return;
}


/*-------------------------------------------------------------------------------------*
 | test special quantity not associated with a particular element or node   fang 01/17 |
 *-------------------------------------------------------------------------------------*/
void STI::STIResultTest::TestSpecial(
    DRT::INPUT::LineDefinition& res,  //!< input file line containing result test specification
    int& nerr,                        //!< number of failed result tests
    int& test_count                   ///< number of result tests
)
{
  // make sure that quantity is tested only by one processor
  if (sti_algorithm_->Comm().MyPID() == 0)
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

  return;
}  // STI::STIResultTest::TestSpecial


/*----------------------------------------------------------------------*
 | get special result to be tested                           fang 01/17 |
 *----------------------------------------------------------------------*/
double STI::STIResultTest::ResultSpecial(
    const std::string& quantity  //! name of quantity to be tested
    ) const
{
  // initialize variable for result
  double result(0.);

  // number of Newton-Raphson iterations in last time step
  if (quantity == "numiterlastnonlinearsolve") result = (double)sti_algorithm_->Iter();

  // number of iterations performed by linear solver during last Newton-Raphson iteration
  else if (quantity == "numiterlastlinearsolve")
  {
    // safety check
    if (STIMonolithic().Solver().Params().get("solver", "none") != "aztec")
      dserror(
          "Must have Aztec solver for result test involving number of solver iterations during "
          "last Newton-Raphson iteration!");
    result = (double)STIMonolithic().Solver().getNumIters();
  }

  // catch unknown quantity strings
  else
    dserror(
        "Quantity '%s' not supported by result testing functionality for scatra-thermo "
        "interaction!",
        quantity.c_str());

  return result;
}  // STI::STIResultTest::ResultSpecial


/*------------------------------------------------------------------------------*
 | return time integrator for monolithic scatra-thermo interaction   fang 09/17 |
 *------------------------------------------------------------------------------*/
const STI::Monolithic& STI::STIResultTest::STIMonolithic() const
{
  const STI::Monolithic* const sti_monolithic =
      dynamic_cast<const STI::Monolithic* const>(sti_algorithm_.get());
  if (sti_monolithic == NULL)
    dserror("Couldn't access time integrator for monolithic scatra-thermo interaction!");
  return *sti_monolithic;
}

/*----------------------------------------------------------------------*/
/*!
\file ssi_monolithic_resulttest.cpp

\brief result testing functionality for monolithic scalar-structure interaction problems

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#include "ssi_monolithic_resulttest.H"

#include "ssi_monolithic.H"

#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 11/17 |
 *----------------------------------------------------------------------*/
SSI::SSI_Mono_ResultTest::SSI_Mono_ResultTest(
    const Teuchos::RCP<const SSI::SSI_Mono>   ssi_mono   //!< time integrator for monolithic scalar-structure interaction
    )
    // call base class constructor
  : DRT::ResultTest("SSI"),

    // store pointer to time integrator for monolithic scalar-structure interaction
    ssi_mono_(ssi_mono)
{
  return;
}


/*-------------------------------------------------------------------------------------*
 | test special quantity not associated with a particular element or node   fang 11/17 |
 *-------------------------------------------------------------------------------------*/
void SSI::SSI_Mono_ResultTest::TestSpecial(
    DRT::INPUT::LineDefinition&   res,         //!< input file line containing result test specification
    int&                          nerr,        //!< number of failed result tests
    int&                          test_count   ///< number of result tests
    )
{
  // make sure that quantity is tested only by one processor
  if(ssi_mono_->Comm().MyPID() == 0)
  {
    // extract name of quantity to be tested
    std::string quantity;
    res.ExtractString("QUANTITY",quantity);

    // get result to be tested
    const double result = ResultSpecial(quantity);

    // compare values
    const int err = CompareValues(result,"SPECIAL",res);
    nerr += err;
    ++test_count;
  }

  return;
} // SSI::SSI_Mono_ResultTest::TestSpecial


/*----------------------------------------------------------------------*
 | get special result to be tested                           fang 11/17 |
 *----------------------------------------------------------------------*/
double SSI::SSI_Mono_ResultTest::ResultSpecial(
    const std::string&   quantity   //! name of quantity to be tested
    ) const
{
  // initialize variable for result
  double result(0.);

  // number of Newton-Raphson iterations in last time step
  if(quantity == "numiterlastnonlinearsolve")
    result = (double) ssi_mono_->Iter();

  // catch unknown quantity strings
  else
    dserror("Quantity '%s' not supported by result testing functionality for monolithic scalar-structure interaction!",quantity.c_str());

  return result;
} // SSI::SSI_Mono_ResultTest::ResultSpecial

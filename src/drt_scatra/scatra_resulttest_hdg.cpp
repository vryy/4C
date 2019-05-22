/*----------------------------------------------------------------------*/
/*!

\brief result tests for HDG problems

\level 3

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/
#include "scatra_timint_hdg.H"
#include "scatra_resulttest_hdg.H"

/*----------------------------------------------------------------------*
 | constructor                                           hoermann 09/15 |
 *----------------------------------------------------------------------*/
SCATRA::HDGResultTest::HDGResultTest(const Teuchos::RCP<ScaTraTimIntImpl> timint)
    : ScaTraResultTest::ScaTraResultTest(timint),
      scatratiminthdg_(Teuchos::rcp_dynamic_cast<const TimIntHDG>(timint))

{
  return;
}


/*----------------------------------------------------------------------*
 | get nodal result to be tested                         hoermann 09/15 |
 *----------------------------------------------------------------------*/
double SCATRA::HDGResultTest::ResultNode(
    const std::string quantity,  //! name of quantity to be tested
    DRT::Node* node              //! node carrying the result to be tested
    ) const
{
  // initialize variable for result
  double result(0.);

  // extract row map from solution vector
  const Epetra_BlockMap& phinpmap = scatratiminthdg_->InterpolatedPhinp()->Map();

  // test result value of single scalar field (averaged value on element node is tested)
  if (quantity == "phi")
    result = (*scatratiminthdg_->InterpolatedPhinp())[phinpmap.LID(node->Id())];
  // catch unknown quantity strings
  else
    dserror("Quantity '%s' not supported in result test!", quantity.c_str());

  return result;
}  // SCATRA::HDGResultTest::ResultNode

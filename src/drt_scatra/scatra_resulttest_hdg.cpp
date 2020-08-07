/*----------------------------------------------------------------------*/
/*! \file

\brief result tests for HDG problems

\level 3


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
  errors_ = scatratiminthdg_->ComputeError();
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
  else if (quantity == "abs_L2error_phi")
    result = std::sqrt((*errors_)[0]);
  else if (quantity == "rel_L2error_phi")
    result = std::sqrt((*errors_)[0] / (*errors_)[1]);
  else if (quantity == "abs_L2error_gradPhi")
    result = std::sqrt((*errors_)[2]);
  else if (quantity == "rel_L2error_gradPhi")
    result = std::sqrt((*errors_)[2] / (*errors_)[3]);
  // catch unknown quantity strings
  else
    dserror("Quantity '%s' not supported in result test!", quantity.c_str());

  return result;
}  // SCATRA::HDGResultTest::ResultNode

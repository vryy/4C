/*--------------------------------------------------------------------------*/
/*!
\file lubrication_resulttest.cpp

\brief testing of lubrication calculation results

\level 3

\maintainer Mostafa Faraji

*/
/*--------------------------------------------------------------------------*/

#include "lubrication_timint_implicit.H"
#include "lubrication_resulttest.H"

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"


/*----------------------------------------------------------------------*
 | ctor                                                     wirtz 11/15 |
 *----------------------------------------------------------------------*/
LUBRICATION::ResultTest::ResultTest(Teuchos::RCP<TimIntImpl> lubrication)
    : DRT::ResultTest("LUBRICATION"),
      dis_(lubrication->Discretization()),
      mysol_(lubrication->Prenp()),
      mynumiter_(lubrication->IterNum())
{
  return;
}


/*----------------------------------------------------------------------*
 | test node                                                wirtz 11/15 |
 *----------------------------------------------------------------------*/
void LUBRICATION::ResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS", dis);
  if (dis != dis_->Name()) return;

  int node;
  res.ExtractInt("NODE", node);
  node -= 1;

  int havenode(dis_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  dis_->Comm().SumAll(&havenode, &isnodeofanybody, 1);

  if (isnodeofanybody == 0)
  {
    dserror("Node %d does not belong to discretization %s", node + 1, dis_->Name().c_str());
  }
  else
  {
    if (dis_->HaveGlobalNode(node))
    {
      DRT::Node* actnode = dis_->gNode(node);

      // Here we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->Owner() != dis_->Comm().MyPID()) return;

      // extract name of quantity to be tested
      std::string quantity;
      res.ExtractString("QUANTITY", quantity);

      // get result to be tested
      const double result = ResultNode(quantity, actnode);

      nerr += CompareValues(result, "NODE", res);
      test_count++;
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | get nodal result to be tested                            wirtz 11/15 |
 *----------------------------------------------------------------------*/
double LUBRICATION::ResultTest::ResultNode(
    const std::string quantity,  //! name of quantity to be tested
    DRT::Node* node              //! node carrying the result to be tested
    ) const
{
  // initialize variable for result
  double result(0.);

  // extract row map from solution vector
  const Epetra_BlockMap& prenpmap = mysol_->Map();

  // test result value of pressure field
  if (quantity == "pre") result = (*mysol_)[prenpmap.LID(dis_->Dof(0, node, 0))];

  // catch unknown quantity strings
  else
    dserror("Quantity '%s' not supported in result test!", quantity.c_str());

  return result;
}  // LUBRICATION::ResultTest::ResultNode


/*-------------------------------------------------------------------------------------*
 | test special quantity not associated with a particular element or node  wirtz 11/15 |
 *-------------------------------------------------------------------------------------*/
void LUBRICATION::ResultTest::TestSpecial(
    DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // make sure that quantity is tested only once
  if (dis_->Comm().MyPID() == 0)
  {
    // extract name of quantity to be tested
    std::string quantity;
    res.ExtractString("QUANTITY", quantity);

    // get result to be tested
    const double result = ResultSpecial(quantity);

    // compare values
    const int err = CompareValues(result, "SPECIAL", res);
    nerr += err;
    test_count++;
  }

  return;
}


/*----------------------------------------------------------------------*
 | get special result to be tested                          wirtz 11/15 |
 *----------------------------------------------------------------------*/
double LUBRICATION::ResultTest::ResultSpecial(
    const std::string quantity  //! name of quantity to be tested
    ) const
{
  // initialize variable for result
  double result(0.);

  if (quantity == "numiterlastnewton") result = (double)mynumiter_;
  // catch unknown quantity strings
  else
    dserror("Quantity '%s' not supported in result test!", quantity.c_str());

  return result;
}  // LUBRICATION::ResultTest::ResultSpecial

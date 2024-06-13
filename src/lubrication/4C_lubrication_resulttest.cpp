/*--------------------------------------------------------------------------*/
/*! \file

\brief testing of lubrication calculation results

\level 3


*/
/*--------------------------------------------------------------------------*/

#include "4C_lubrication_resulttest.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_lubrication_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | ctor                                                     wirtz 11/15 |
 *----------------------------------------------------------------------*/
LUBRICATION::ResultTest::ResultTest(Teuchos::RCP<TimIntImpl> lubrication)
    : Core::UTILS::ResultTest("LUBRICATION"),
      dis_(lubrication->discretization()),
      mysol_(lubrication->Prenp()),
      mynumiter_(lubrication->IterNum())
{
  return;
}


/*----------------------------------------------------------------------*
 | test node                                                wirtz 11/15 |
 *----------------------------------------------------------------------*/
void LUBRICATION::ResultTest::test_node(Input::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.extract_string("DIS", dis);
  if (dis != dis_->Name()) return;

  int node;
  res.extract_int("NODE", node);
  node -= 1;

  int havenode(dis_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  dis_->Comm().SumAll(&havenode, &isnodeofanybody, 1);

  if (isnodeofanybody == 0)
  {
    FOUR_C_THROW("Node %d does not belong to discretization %s", node + 1, dis_->Name().c_str());
  }
  else
  {
    if (dis_->HaveGlobalNode(node))
    {
      Core::Nodes::Node* actnode = dis_->gNode(node);

      // Here we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->Owner() != dis_->Comm().MyPID()) return;

      // extract name of quantity to be tested
      std::string quantity;
      res.extract_string("QUANTITY", quantity);

      // get result to be tested
      const double result = result_node(quantity, actnode);

      nerr += compare_values(result, "NODE", res);
      test_count++;
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | get nodal result to be tested                            wirtz 11/15 |
 *----------------------------------------------------------------------*/
double LUBRICATION::ResultTest::result_node(
    const std::string quantity,  //! name of quantity to be tested
    Core::Nodes::Node* node      //! node carrying the result to be tested
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
    FOUR_C_THROW("Quantity '%s' not supported in result test!", quantity.c_str());

  return result;
}  // LUBRICATION::ResultTest::ResultNode


/*-------------------------------------------------------------------------------------*
 | test special quantity not associated with a particular element or node  wirtz 11/15 |
 *-------------------------------------------------------------------------------------*/
void LUBRICATION::ResultTest::TestSpecial(Input::LineDefinition& res, int& nerr, int& test_count)
{
  // make sure that quantity is tested only once
  if (dis_->Comm().MyPID() == 0)
  {
    // extract name of quantity to be tested
    std::string quantity;
    res.extract_string("QUANTITY", quantity);

    // get result to be tested
    const double result = result_special(quantity);

    // compare values
    const int err = compare_values(result, "SPECIAL", res);
    nerr += err;
    test_count++;
  }

  return;
}


/*----------------------------------------------------------------------*
 | get special result to be tested                          wirtz 11/15 |
 *----------------------------------------------------------------------*/
double LUBRICATION::ResultTest::result_special(
    const std::string quantity  //! name of quantity to be tested
) const
{
  // initialize variable for result
  double result(0.);

  if (quantity == "numiterlastnewton") result = (double)mynumiter_;
  // catch unknown quantity strings
  else
    FOUR_C_THROW("Quantity '%s' not supported in result test!", quantity.c_str());

  return result;
}  // LUBRICATION::ResultTest::result_special

FOUR_C_NAMESPACE_CLOSE

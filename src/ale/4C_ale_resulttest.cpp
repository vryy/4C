/*----------------------------------------------------------------------------*/
/*! \file

\brief Result tests for pure ALE problems

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "4C_ale_resulttest.hpp"

#include "4C_ale.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ALE::AleResultTest::AleResultTest(ALE::Ale& ale)
    : CORE::UTILS::ResultTest("ALE"), aledis_(ale.discretization()), dispnp_(ale.Dispnp())
{
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::AleResultTest::test_node(INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS", dis);
  if (dis != aledis_->Name()) return;

  int node;
  res.ExtractInt("NODE", node);
  node -= 1;

  int havenode(aledis_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  aledis_->Comm().SumAll(&havenode, &isnodeofanybody, 1);

  if (isnodeofanybody == 0)
  {
    FOUR_C_THROW("Node %d does not belong to discretization %s", node + 1, aledis_->Name().c_str());
  }
  else
  {
    if (aledis_->HaveGlobalNode(node))
    {
      CORE::Nodes::Node* actnode = aledis_->gNode(node);

      // Here we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->Owner() != aledis_->Comm().MyPID()) return;

      double result = 0.;

      const Epetra_BlockMap& dispnpmap = dispnp_->Map();

      std::string position;
      res.ExtractString("QUANTITY", position);
      if (position == "dispx")
      {
        result = (*dispnp_)[dispnpmap.LID(aledis_->Dof(actnode, 0))];
      }
      else if (position == "dispy")
      {
        result = (*dispnp_)[dispnpmap.LID(aledis_->Dof(actnode, 1))];
      }
      else if (position == "dispz")
      {
        result = (*dispnp_)[dispnpmap.LID(aledis_->Dof(actnode, 2))];
      }
      else
      {
        FOUR_C_THROW("Quantity '%s' not supported in ALE testing", position.c_str());
      }

      nerr += compare_values(result, "NODE", res);
      test_count++;
    }
  }
}

FOUR_C_NAMESPACE_CLOSE

/*----------------------------------------------------------------------------*/
/*!
\file ale_resulttest.cpp

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "ale_resulttest.H"

#include "ale.H"

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ALE::AleResultTest::AleResultTest(ALE::Ale& ale)
    : DRT::ResultTest("ALE"), aledis_(ale.Discretization()), dispnp_(ale.Dispnp())
{
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ALE::AleResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
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
    dserror("Node %d does not belong to discretization %s", node + 1, aledis_->Name().c_str());
  }
  else
  {
    if (aledis_->HaveGlobalNode(node))
    {
      DRT::Node* actnode = aledis_->gNode(node);

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
        dserror("Quantity '%s' not supported in ALE testing", position.c_str());
      }

      nerr += CompareValues(result, "NODE", res);
      test_count++;
    }
  }
}

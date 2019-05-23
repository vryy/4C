/*----------------------------------------------------------------------*/
/*!
\file acou_resulttest.cpp

\brief testing of acoustical calculation results

\level 2

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
*/
/*----------------------------------------------------------------------*/


#include "acou_resulttest.H"
#include "acou_impl.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*
 |                                                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
ACOU::AcouResultTest::AcouResultTest(AcouTimeInt& acoualgo) : DRT::ResultTest("ACOUSTIC")
{
  dis_ = acoualgo.Discretization();
  mysol_ = LINALG::CreateVector(*(dis_->NodeRowMap()), true);
  acoualgo.NodalPressureField(mysol_);
}

/*----------------------------------------------------------------------*
 |                                                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
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

      // Here, we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->Owner() != dis_->Comm().MyPID()) return;

      double result = 0.;
      const Epetra_BlockMap& map = mysol_->Map();
      std::string position;
      res.ExtractString("QUANTITY", position);

      if (position == "pressure")
      {
        result = (*mysol_)[map.LID(
            actnode->Id())];  // Here, we got a node map, not a dof map with dofs based on nodes!
      }
      else
      {
        dserror("Quantity '%s' not supported in result-test of acoustic transport problems",
            position.c_str());
      }

      nerr += CompareValues(result, "NODE", res);
      test_count++;
    }
  }
  return;
}

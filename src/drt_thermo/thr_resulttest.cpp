/*----------------------------------------------------------------------*/
/*!
\file thr_resulttest.cpp

\brief tesing of structure calculation results

\level 1
<pre>
\maintainer Christoph Meier
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 08/09 |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                   dano 08/09 |
 *----------------------------------------------------------------------*/
#include <string>

#include "thr_resulttest.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |                                                           dano 08/09 |
 *----------------------------------------------------------------------*/
THR::ResultTest::ResultTest(TimInt& tintegrator) : DRT::ResultTest("THERMAL")
{
  temp_ = tintegrator.Temp();
  rate_ = tintegrator.Rate();
  thrdisc_ = tintegrator.Discretization();
}

/*----------------------------------------------------------------------*
 |                                                           dano 08/09 |
 *----------------------------------------------------------------------*/
void THR::ResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS", dis);
  if (dis != thrdisc_->Name()) return;

  int node;
  res.ExtractInt("NODE", node);
  node -= 1;

  int havenode(thrdisc_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  thrdisc_->Comm().SumAll(&havenode, &isnodeofanybody, 1);

  if (isnodeofanybody == 0)
  {
    dserror("Node %d does not belong to discretization %s", node + 1, thrdisc_->Name().c_str());
  }
  else
  {
    // this implementation does not allow testing of heatfluxes
    if (thrdisc_->HaveGlobalNode(node))
    {
      const DRT::Node* actnode = thrdisc_->gNode(node);

      // Here we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->Owner() != thrdisc_->Comm().MyPID()) return;

      std::string position;
      res.ExtractString("QUANTITY", position);
      bool unknownpos = true;  // make sure the result value std::string can be handled
      double result = 0.0;     // will hold the actual result of run

      // test temperature
      if (temp_ != Teuchos::null)
      {
        const Epetra_BlockMap& tempmap = temp_->Map();

        if (position == "temp")
        {
          unknownpos = false;
          result = (*temp_)[tempmap.LID(thrdisc_->Dof(0, actnode, 0))];
        }
      }

      // test temperature rates
      if (rate_ != Teuchos::null)
      {
        const Epetra_BlockMap& ratemap = rate_->Map();

        if (position == "rate")
        {
          unknownpos = false;
          result = (*rate_)[ratemap.LID(thrdisc_->Dof(0, actnode, 0))];
        }
      }

      // catch position strings, which are not handled by thermo result test
      if (unknownpos) dserror("Quantity '%s' not supported in thermo testing", position.c_str());

      // compare values
      const int err = CompareValues(result, "NODE", res);
      nerr += err;
      test_count++;
    }
  }
}  // TestNode

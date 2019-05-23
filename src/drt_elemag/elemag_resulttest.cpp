/*----------------------------------------------------------------------*/
/*!
\file elemag_resulttest.cpp

\brief testing of electromagnetic calculation results

\level 2

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            089 - 289-15244
*/
/*----------------------------------------------------------------------*/


#include "elemag_resulttest.H"
#include "elemag_timeint.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*
 |                                                     berardocco 10/18 |
 *----------------------------------------------------------------------*/
ELEMAG::ElemagResultTest::ElemagResultTest(ElemagTimeInt& elemagalgo)
    : DRT::ResultTest("ELECTROMAGNETIC")
{
  dis_ = elemagalgo.Discretization();
  // mysol_ = LINALG::CreateVector(*(dis_->NodeRowMap()), true);
  error_ = Teuchos::rcp(new Epetra_SerialDenseVector(4));
  error_ = elemagalgo.ComputeError();
  // elemagalgo.NodalPressureField(mysol_);
}

/*----------------------------------------------------------------------*
 |                                                     berardocco 10/18 |
 *----------------------------------------------------------------------*/
void ELEMAG::ElemagResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
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
      // const Epetra_BlockMap& map = mysol_->Map();
      std::string variable;
      res.ExtractString("QUANTITY", variable);

      if (variable == "electric")
      {
        std::cout << "This error function has not yet been implemented." << std::endl;
        // result = (*mysol_)[map.LID(actnode->Id())];  // Here, we got a node map, not a dof map
        // with dofs based on nodes!
      }
      else if (variable == "L2electric")
      {
        result = sqrt((*error_)[0]);
      }
      else if (variable == "L2magnetic")
      {
        result = sqrt((*error_)[2]);
      }

      else
      {
        dserror("Quantity '%s' not supported in result-test of elemagstic transport problems",
            variable.c_str());
      }

      nerr += CompareValues(result, "NODE", res);
      test_count++;
    }
  }
  return;
}
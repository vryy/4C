/*!------------------------------------------------------------------------------------------------*
\file opti_resulttest.cpp

\brief Result test for optimization algorithms

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "opti_resulttest.H"
#include "opti_GCMMA.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
OPTI::OptiResultTest::OptiResultTest(
    const GCMMA& optimizer
)
: DRT::ResultTest("OPTI"),
  optidis_(optimizer.Discretization()),
  sol_(optimizer.X())
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void OPTI::OptiResultTest::TestNode(
    DRT::INPUT::LineDefinition& res,
    int& nerr,
    int& test_count
)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS",dis);
  if (dis != optidis_->Name())
    return;

  int nodeGid;
  res.ExtractInt("NODE",nodeGid);
  nodeGid -= 1;

  int havenode(optidis_->HaveGlobalNode(nodeGid));
  int isnodeofanybody(0);
  optidis_->Comm().SumAll(&havenode,&isnodeofanybody,1);

  if (isnodeofanybody==0)
  {
    dserror("Node %d does not belong to discretization %s",nodeGid+1,optidis_->Name().c_str());
  }
  else
  {
    if (optidis_->HaveGlobalNode(nodeGid))
    {
      const DRT::Node* node = optidis_->gNode(nodeGid);

      // Test only, if actnode is a row node
      if (node->Owner() != optidis_->Comm().MyPID())
        return;

      double result = 0.;

      const Epetra_BlockMap& optimap = sol_->Map();

      std::string position;
      res.ExtractString("QUANTITY",position);
      if (position=="x")
        result = (*sol_)[optimap.LID(optidis_->Dof(0,node,0))];
      else
        dserror("Quantity '%s' not supported in fluid testing", position.c_str());

      nerr += CompareValues(result, res);
      test_count++;
    }
  }
}

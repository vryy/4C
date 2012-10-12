/*!------------------------------------------------------------------------------------------------*
\file opti_resulttest.cpp

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


/*!------------------------------------------------------------------------------------------------*
\file topopt_fluidAdjointResulttest.cpp

\brief

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


//#include <string>
//
#include "opti_resulttest.H"
#include "opti_GCMMA.H"
#include "../drt_lib/drt_discret.H"
//#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
OPTI::OptiResultTest::OptiResultTest(
    const GCMMA& optimizer
)
: optidis_(optimizer.Discretization()),
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
  int dis;
  res.ExtractInt("DIS",dis);
  if (dis != 1)
    dserror("fix me: only one optimization discretization supported for testing");

  int nodeGid;
  res.ExtractInt("NODE",nodeGid);
  nodeGid -= 1;

  if (optidis_->HaveGlobalNode(nodeGid))
  {
    const DRT::Node* node = optidis_->gNode(nodeGid);

    // Test only, if actnode is a row node
    if (node->Owner() != optidis_->Comm().MyPID())
      return;

    double result = 0.;

    const Epetra_BlockMap& optimap = sol_->Map();

    std::string position;
    res.ExtractString("POSITION",position);
    if (position=="x")
      result = (*sol_)[optimap.LID(optidis_->Dof(0,node,0))];
    else
      dserror("position '%s' not supported in fluid testing", position.c_str());

    nerr += CompareValues(result, res);
    test_count++;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool OPTI::OptiResultTest::Match(DRT::INPUT::LineDefinition& res)
{
  return res.HaveNamed("OPTI");
}

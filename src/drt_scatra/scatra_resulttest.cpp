/*----------------------------------------------------------------------*/
/*!
\file scatra_resulttest.cpp

\brief testing of scalar transport calculation results

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_timint_implicit.H"
#include "scatra_resulttest.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SCATRA::ScaTraResultTest::ScaTraResultTest(ScaTraTimIntImpl& scatra)
{
  dis_    = scatra.Discretization();
  mysol_  = scatra.Phinp();
  myflux_ = scatra.Flux();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  int dis;
  res.ExtractInt("DIS",dis);
  if (dis != 1)
    dserror("fix me: only one scalar transport discretization supported for testing");

  int node;
  res.ExtractInt("NODE",node);
  node -= 1;

  if (dis_->HaveGlobalNode(node))
  {
    DRT::Node* actnode = dis_->gNode(node);

    // Strange! It seems we might actually have a global node around
    // even if it does not belong to us. But here we are just
    // interested in our nodes!
    if (actnode->Owner() != dis_->Comm().MyPID())
      return;

    double result = 0.;
    const Epetra_BlockMap& phinpmap = mysol_->Map();
    std::string position;
    res.ExtractString("POSITION",position);

    // test result value of single scalar field
    if (position=="phi")
      result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,0))];
    // test result values for a system of scalars
    else if (position=="phi1")
      result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,0))];
    else if (position=="phi2")
      result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,1))];
    else if (position=="phi3")
      result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,2))];
    else if (position=="phi4")
      result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,3))];
    else if (position=="phi5")
      result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,4))];
    else if (position=="phi6")
      result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,5))];
    else if (position=="phi7")
      result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,6))];
    // we support only testing of fluxes for the first scalar
    else if (position=="fluxx")
      result = ((*myflux_)[0])[phinpmap.LID(dis_->Dof(0,actnode,0))];
    else if (position=="fluxy")
      result = ((*myflux_)[1])[phinpmap.LID(dis_->Dof(0,actnode,0))];
    else if (position=="fluxz")
      result = ((*myflux_)[2])[phinpmap.LID(dis_->Dof(0,actnode,0))];
    else
    {
      dserror("position '%s' not supported in result-test of scalar transport problems", position.c_str());
    }

    nerr += CompareValues(result, res);
    test_count++;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool SCATRA::ScaTraResultTest::Match(DRT::INPUT::LineDefinition& res)
{
  return res.HaveNamed("SCATRA");
}


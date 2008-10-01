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

#ifdef CCADISCRET

#include "scatra_resulttest.H"

#ifdef PARALLEL
#include <mpi.h>
#endif


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SCATRA::ScaTraResultTest::ScaTraResultTest(ScaTraTimIntImpl& scatra)
{
  dis_    = scatra.discret_;
  mysol_  = scatra.phinp_ ;
  myflux_ = scatra.CalcFlux();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraResultTest::TestNode(const RESULTDESCR* res, int& nerr, int& test_count)
{
  if (res->dis != 0)
    dserror("fix me: only one scalar transport discretization supported for testing");

  if (dis_->HaveGlobalNode(res->node))
  {
    DRT::Node* actnode = dis_->gNode(res->node);

    // Strange! It seems we might actually have a global node around
    // even if it does not belong to us. But here we are just
    // interested in our nodes!
    if (actnode->Owner() != dis_->Comm().MyPID())
      return;

    double result = 0.;
    const Epetra_BlockMap& phinpmap = mysol_->Map();
    string position = res->position;

    // we only can test the first(!) scalar of a system (e.g. first dof at every node)
    if (position=="phi")
      result = (*mysol_)[phinpmap.LID(dis_->Dof(actnode,0))];
    else if (position=="phi1")
      result = (*mysol_)[phinpmap.LID(dis_->Dof(actnode,0))];
    else if (position=="phi2")
      result = (*mysol_)[phinpmap.LID(dis_->Dof(actnode,1))];
    else if (position=="phi3")
      result = (*mysol_)[phinpmap.LID(dis_->Dof(actnode,2))];
    else if (position=="phi4")
      result = (*mysol_)[phinpmap.LID(dis_->Dof(actnode,3))];
    else if (position=="phi5")
      result = (*mysol_)[phinpmap.LID(dis_->Dof(actnode,4))];
    else if (position=="phi6")
      result = (*mysol_)[phinpmap.LID(dis_->Dof(actnode,5))];
    // we rely on the fact, that we have 
    else if (position=="fluxx")
      result = (*myflux_)[0][phinpmap.LID(dis_->Dof(actnode,0))];
    else if (position=="fluxy")
      result = (*myflux_)[1][phinpmap.LID(dis_->Dof(actnode,0))];
    else if (position=="fluxz")
      result = (*myflux_)[2][phinpmap.LID(dis_->Dof(actnode,0))];
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
bool SCATRA::ScaTraResultTest::Match(const RESULTDESCR* res)
{
  return res->field==scatra;
}


#endif /* CCADISCRET       */

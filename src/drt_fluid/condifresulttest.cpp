/*----------------------------------------------------------------------*/
/*!
\file condifresulttest.cpp

\brief tesing of convection-diffusion calculation results

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <string>

#include "condifresulttest.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C" /* stuff which is c and is accessed from c++ */
{
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ConDifResultTest::ConDifResultTest(CondifImplicitTimeInt& condif)
{
  dis_=condif.discret_;
  mysol_   =condif.phinp_ ;
  myflux_ = condif.CalcFlux();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ConDifResultTest::ConDifResultTest(CondifGenAlphaIntegration& condif)
{
  dis_=condif.discret_;
  mysol_   =condif.phinp_ ;
  dserror("flux calculation method missing for GenAlpha-Result test");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ConDifResultTest::TestNode(const RESULTDESCR* res, int& nerr, int& test_count)
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

    if (position=="phi")
    {
      result = (*mysol_)[phinpmap.LID(dis_->Dof(actnode,0))];
    }
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
bool ConDifResultTest::Match(const RESULTDESCR* res)
{
  return res->field==scatra;
}


#endif /* CCADISCRET       */

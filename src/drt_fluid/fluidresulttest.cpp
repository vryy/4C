/*----------------------------------------------------------------------*/
/*!
\file fluidresulttest.cpp

\brief tesing of fluid calculation results

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <string>

#include "fluidresulttest.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C" /* stuff which is c and is accessed from c++ */
{
//#include "../headers/standardtypes.h"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FluidResultTest::FluidResultTest(FluidImplicitTimeInt& fluid)
{
  fluiddis_=fluid.discret_;
  mysol_   =fluid.velnp_ ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FluidResultTest::FluidResultTest(FluidGenAlphaIntegration& fluid)
{
  fluiddis_=fluid.discret_;
  mysol_   =fluid.velnp_ ;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FluidResultTest::TestNode(RESULTDESCR* res, int& nerr, int& test_count)
{
  if (res->dis != 0)
    dserror("fix me: only one fluid discretization supported for testing");

  if (fluiddis_->HaveGlobalNode(res->node))
  {
    DRT::Node* actnode = fluiddis_->gNode(res->node);

    // Strange! It seems we might actually have a global node around
    // even if it does not belong to us. But here we are just
    // interested in our nodes!
    if (actnode->Owner() != fluiddis_->Comm().MyPID())
      return;

    double result = 0.;

    const Epetra_BlockMap& velnpmap = mysol_->Map();

    string position = res->position;
    if (position=="velx")
    {
      result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(actnode,0))];
    }
    else if (position=="vely")
    {
      result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(actnode,1))];
    }
    else if (position=="velz")
    {
      result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(actnode,2))];
    }
    else if (position=="pressure")
    {
      if (genprob.ndim==2)
      {
        if (fluiddis_->NumDof(actnode)<3)
          dserror("too few dofs at node %d for pressure testing",actnode->Id());
        result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(actnode,2))];
      }
      else
      {
        if (fluiddis_->NumDof(actnode)<4)
          dserror("too few dofs at node %d for pressure testing",actnode->Id());
        result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(actnode,3))];
      }
    }
    else
    {
      dserror("position '%s' not supported in fluid testing", position.c_str());
    }

    nerr += CompareValues(result, res);
    test_count++;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FluidResultTest::Match(RESULTDESCR* res)
{
  return res->field==fluid;
}


#endif /* CCADISCRET       */

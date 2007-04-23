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
#ifdef TRILINOS_PACKAGE
#ifdef D_FLUID
#ifdef RESULTTEST

#include <string>

#include "fluidresulttest.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"

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
  : fluid_(fluid)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FluidResultTest::TestNode(RESULTDESCR* res, int& nerr, int& test_count)
{
  if (res->dis != 0)
    dserror("fix me: only one fluid discretization supported for testing");

  if (fluid_.discret_->HaveGlobalNode(res->node))
  {
    DRT::Node* actnode = fluid_.discret_->gNode(res->node);

    // Strange! It seems we might actually have a global node around
    // even if it does not belong to us. But here we are just
    // interested in our nodes!
    if (actnode->Owner() != fluid_.discret_->Comm().MyPID())
      return;

    double result = 0.;

    const Epetra_BlockMap& velnpmap = fluid_.velnp_->Map();

    string position = res->position;
    if (position=="velx")
    {
      result = (*fluid_.velnp_)[velnpmap.LID(actnode->Dof()[0])];
    }
    else if (position=="vely")
    {
      result = (*fluid_.velnp_)[velnpmap.LID(actnode->Dof()[1])];
    }
    else if (position=="velz")
    {
      result = (*fluid_.velnp_)[velnpmap.LID(actnode->Dof()[2])];
    }
    else if (position=="pressure")
    {
      if (genprob.ndim==2)
      {
        if (actnode->Dof().NumDof()<3)
          dserror("too few dofs at node %d for pressure testing",actnode->Id());
        result = (*fluid_.velnp_)[velnpmap.LID(actnode->Dof()[2])];
      }
      else
      {
        if (actnode->Dof().NumDof()<4)
          dserror("too few dofs at node %d for pressure testing",actnode->Id());
        result = (*fluid_.velnp_)[velnpmap.LID(actnode->Dof()[3])];
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


#endif
#endif /* D_FLUID          */
#endif /* TRILINOS_PACKAGE */
#endif /* CCADISCRET       */

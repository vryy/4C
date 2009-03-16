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
#include "../drt_lib/drt_globalproblem.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::FluidResultTest::FluidResultTest(FluidImplicitTimeInt& fluid)
{
  fluiddis_= fluid.discret_;
  mysol_   = fluid.velnp_ ;
  mytraction_ = fluid.CalcStresses();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::FluidResultTest::FluidResultTest(FluidGenAlphaIntegration& fluid)
{
  fluiddis_= fluid.discret_;
  mysol_   = fluid.velnp_;
  mytraction_ = fluid.CalcStresses();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::FluidResultTest::TestNode(const RESULTDESCR* res, int& nerr, int& test_count)
{
  if (res->dis != 0)
    dserror("fix me: only one fluid discretization supported for testing");

  if (fluiddis_->HaveGlobalNode(res->node))
  {
    const DRT::Node* actnode = fluiddis_->gNode(res->node);

    // Test only, if actnode is a row node
    if (actnode->Owner() != fluiddis_->Comm().MyPID())
      return;

    double result = 0.;

    const Epetra_BlockMap& velnpmap = mysol_->Map();

    const int numdim = DRT::Problem::Instance()->ProblemSizeParams().get<int>("DIM");

    const string position = res->position;
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
      if (numdim==2)
        dserror("Cannot test result for velz in 2D case.");
      result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(actnode,2))];
    }
    else if (position=="pressure")
    {
      if (numdim==2)
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
    else if (position=="tractionx")
      result = (*mytraction_)[(mytraction_->Map()).LID(fluiddis_->Dof(actnode,0))];
    else if (position=="tractiony")
      result = (*mytraction_)[(mytraction_->Map()).LID(fluiddis_->Dof(actnode,1))];
    else if (position=="tractionz")
    {
      if (numdim==2)
        dserror("Cannot test result for tractionz in 2D case.");
      result = (*mytraction_)[(mytraction_->Map()).LID(fluiddis_->Dof(actnode,2))];
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
bool FLD::FluidResultTest::Match(const RESULTDESCR* res)
{
  return res->field==fluid;
}


#endif /* CCADISCRET       */

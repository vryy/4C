/*----------------------------------------------------------------------*/
/*!
\file fluidxfluidresulttest.cpp

\brief tesing of fluid calculation results

<pre>
Maintainer: Shadan Shahmiri
            shahmiri@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <string>

#include "fluidxfluidresulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_linedefinition.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::FluidXFluidResultTest::FluidXFluidResultTest(FluidXFluidImplicitTimeInt& fluid)
{
    fluiddis_= fluid.fluiddis_;
    xfluiddis_= fluid.xfluiddis_;
    fluidsol_   = fluid.fluidstate_.velnp_ ;
    xfluidsol_   = fluid.xfluidstate_.velnp_ ;
    mytraction_ = fluid.CalcStresses();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::FluidXFluidResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  int dis;
  res.ExtractInt("DIS",dis);
//  if (dis != 1)
//        dserror("fix me: only one fluid discretization supported for testing");

  int node;
  res.ExtractInt("NODE",node);
  //node -= 1;

  if (dis == 1)
  {  
    if (fluiddis_->HaveGlobalNode(node))
    {
        const DRT::Node* actnode = fluiddis_->gNode(node);

        // Test only, if actnode is a row node
        if (actnode->Owner() != fluiddis_->Comm().MyPID())
            return;

        double result = 0.;

        const Epetra_BlockMap& velnpmap = fluidsol_->Map();

        const int numdim = DRT::Problem::Instance()->ProblemSizeParams().get<int>("DIM");

    std::string position;
    res.ExtractString("POSITION",position);
        if (position=="velx")
        {
            result = (*fluidsol_)[velnpmap.LID(fluiddis_->Dof(actnode,0))];
        }
        else if (position=="vely")
        {
            result = (*fluidsol_)[velnpmap.LID(fluiddis_->Dof(actnode,1))];
        }
        else if (position=="velz")
        {
            if (numdim==2)
                dserror("Cannot test result for velz in 2D case.");
            result = (*fluidsol_)[velnpmap.LID(fluiddis_->Dof(actnode,2))];
        }
        else if (position=="pressure")
        {
            if (numdim==2)
            {
                if (fluiddis_->NumDof(actnode)<3)
                    dserror("too few dofs at node %d for pressure testing",actnode->Id());
                result = (*fluidsol_)[velnpmap.LID(fluiddis_->Dof(actnode,2))];
            }
            else
            {
                if (fluiddis_->NumDof(actnode)<4)
                    dserror("too few dofs at node %d for pressure testing",actnode->Id());
                result = (*fluidsol_)[velnpmap.LID(fluiddis_->Dof(actnode,3))];
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
  else if(dis == 2)
  {
    if (xfluiddis_->HaveGlobalNode(node))
    {
        const DRT::Node* actnode = xfluiddis_->gNode(node);

        // Test only, if actnode is a row node
        if (actnode->Owner() != xfluiddis_->Comm().MyPID())
            return;

        double result = 0.;

        const Epetra_BlockMap& velnpmap = xfluidsol_->Map();

        const int numdim = DRT::Problem::Instance()->ProblemSizeParams().get<int>("DIM");

    std::string position;
    res.ExtractString("POSITION",position);
        if (position=="velx")
        {
            result = (*xfluidsol_)[velnpmap.LID(xfluiddis_->Dof(actnode,0))];
        }
        else if (position=="vely")
        {
            result = (*xfluidsol_)[velnpmap.LID(xfluiddis_->Dof(actnode,1))];
        }
        else if (position=="velz")
        {
            if (numdim==2)
                dserror("Cannot test result for velz in 2D case.");
            result = (*xfluidsol_)[velnpmap.LID(xfluiddis_->Dof(actnode,2))];
        }
        else if (position=="pressure")
        {
            if (numdim==2)
            {
                if (xfluiddis_->NumDof(actnode)<3)
                    dserror("too few dofs at node %d for pressure testing",actnode->Id());
                result = (*xfluidsol_)[velnpmap.LID(xfluiddis_->Dof(actnode,2))];
            }
            else
            {
                if (xfluiddis_->NumDof(actnode)<4)
                    dserror("too few dofs at node %d for pressure testing",actnode->Id());
                result = (*xfluidsol_)[velnpmap.LID(xfluiddis_->Dof(actnode,3))];
            }
        }
        else if (position=="tractionx")
            result = (*mytraction_)[(mytraction_->Map()).LID(xfluiddis_->Dof(actnode,0))];
        else if (position=="tractiony")
            result = (*mytraction_)[(mytraction_->Map()).LID(xfluiddis_->Dof(actnode,1))];
        else if (position=="tractionz")
        {
            if (numdim==2)
                dserror("Cannot test result for tractionz in 2D case.");
            result = (*mytraction_)[(mytraction_->Map()).LID(xfluiddis_->Dof(actnode,2))];
        }
        else
        {
            dserror("position '%s' not supported in fluid testing", position.c_str());
        }

        nerr += CompareValues(result, res);
        test_count++;
    } 
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FLD::FluidXFluidResultTest::Match(DRT::INPUT::LineDefinition& res)
{
  return res.HaveNamed("FLUID");
}


#endif /* CCADISCRET       */

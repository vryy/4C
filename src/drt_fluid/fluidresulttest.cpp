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


#include <string>

#include "fluidresulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_linedefinition.H"
#include "fluidimplicitintegration.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::FluidResultTest::FluidResultTest(FluidImplicitTimeInt& fluid)
  : DRT::ResultTest("FLUID")
{
    fluiddis_= fluid.discret_;
    mysol_   = fluid.velnp_ ;
    mytraction_ = fluid.CalcStresses();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::FluidResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS",dis);
  if (dis != fluiddis_->Name())
    return;

  int node;
  res.ExtractInt("NODE",node);
  node -= 1;

  int havenode(fluiddis_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  fluiddis_->Comm().SumAll(&havenode,&isnodeofanybody,1);

  if (isnodeofanybody==0)
  {
    dserror("Node %d does not belong to discretization %s",node+1,fluiddis_->Name().c_str());
  }
  else
  {
    if (fluiddis_->HaveGlobalNode(node))
    {
      const DRT::Node* actnode = fluiddis_->gNode(node);

      // Here we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->Owner() != fluiddis_->Comm().MyPID())
        return;

      double result = 0.;

      const Epetra_BlockMap& velnpmap = mysol_->Map();

      const int numdim = DRT::Problem::Instance()->NDim();

      std::string position;
      res.ExtractString("QUANTITY",position);
      if (position=="velx")
      {
        result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(0,actnode,0))];
      }
      else if (position=="vely")
      {
        result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(0,actnode,1))];
      }
      else if (position=="velz")
      {
        if (numdim==2)
          dserror("Cannot test result for velz in 2D case.");
        result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(0,actnode,2))];
      }
      else if (position=="pressure")
      {
        if (numdim==2)
        {
          if (fluiddis_->NumDof(0,actnode)<3)
            dserror("too few dofs at node %d for pressure testing",actnode->Id());
          result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(0,actnode,2))];
        }
        else
        {
          if (fluiddis_->NumDof(0,actnode)<4)
            dserror("too few dofs at node %d for pressure testing",actnode->Id());
          result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(0,actnode,3))];
        }
      }
      else if (position=="tractionx")
        result = (*mytraction_)[(mytraction_->Map()).LID(fluiddis_->Dof(0,actnode,0))];
      else if (position=="tractiony")
        result = (*mytraction_)[(mytraction_->Map()).LID(fluiddis_->Dof(0,actnode,1))];
      else if (position=="tractionz")
      {
        if (numdim==2)
          dserror("Cannot test result for tractionz in 2D case.");
        result = (*mytraction_)[(mytraction_->Map()).LID(fluiddis_->Dof(0,actnode,2))];
      }
      else
      {
        dserror("Quantity '%s' not supported in fluid testing", position.c_str());
      }

      nerr += CompareValues(result, res);
      test_count++;
    }
  }
}

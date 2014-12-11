/*----------------------------------------------------------------------*/
/*!
\file xfluidresulttest.cpp

\brief tesing of fluid calculation results

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/


#include <string>

#include "xfluid.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_linedefinition.H"

#include "xfluidresulttest.H"


FLD::XFluidResultTest::XFluidResultTest( XFluid * xfluid )
  : DRT::ResultTest("FLUID"),
    discret_( *xfluid->discret_ ),
    velnp_( xfluid->state_->velnp_ )
{
}

void FLD::XFluidResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS",dis);
  if (dis != discret_.Name())
    return;

  int node;
  res.ExtractInt("NODE",node);
  node -= 1;

  int havenode(discret_.HaveGlobalNode(node));
  int isnodeofanybody(0);
  discret_.Comm().SumAll(&havenode,&isnodeofanybody,1);

  if (isnodeofanybody==0)
  {
    dserror("Node %d does not belong to discretization %s",node+1,discret_.Name().c_str());
  }
  else
  {
    if (discret_.HaveGlobalNode(node))
    {
      DRT::Node* actnode = discret_.gNode(node);

      if (actnode->Owner() != discret_.Comm().MyPID())
        return;

      double result = 0.;

      const Epetra_BlockMap& velnpmap = velnp_->Map();

      std::string position;
      res.ExtractString("QUANTITY",position);
      if (position=="velx")
      {
        result = (*velnp_)[velnpmap.LID(discret_.Dof(actnode,0))];
      }
      else if (position=="vely")
      {
        result = (*velnp_)[velnpmap.LID(discret_.Dof(actnode,1))];
      }
      else if (position=="velz")
      {
        result = (*velnp_)[velnpmap.LID(discret_.Dof(actnode,2))];
      }
      else if (position=="pressure")
      {
        result = (*velnp_)[velnpmap.LID(discret_.Dof(actnode,3))];
      }
      else
      {
        dserror("Quantity '%s' not supported in ale testing", position.c_str());
      }

      nerr += CompareValues(result, "NODE", res);
      test_count++;
    }
  }
}

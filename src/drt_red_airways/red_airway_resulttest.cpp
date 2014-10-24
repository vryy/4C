/*----------------------------------------------------------------------*/
/*!
\file red_airway_resulttest.cpp

\brief testing of Red_Airway calculation results

<pre>
Maintainer: Christian Roth
            roth@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/
/*----------------------------------------------------------------------*/


#include "airwayimplicitintegration.H"
#include "red_airway_resulttest.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
AIRWAY::RedAirwayResultTest::RedAirwayResultTest(RedAirwayImplicitTimeInt& airways)
  : DRT::ResultTest("RED_AIRWAY")
{
  dis_    = airways.Discretization();
  mysol_  = airways.Pnp();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void AIRWAY::RedAirwayResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS",dis);
  if (dis != dis_->Name())
    return;

  int node;
  res.ExtractInt("NODE",node);
  node -= 1;

  int havenode(dis_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  dis_->Comm().SumAll(&havenode,&isnodeofanybody,1);

  if (isnodeofanybody==0)
  {
    dserror("Node %d does not belong to discretization %s",node+1,dis_->Name().c_str());
  }
  else
  {
    if (dis_->HaveGlobalNode(node))
    {
      DRT::Node* actnode = dis_->gNode(node);

      // Here we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->Owner() != dis_->Comm().MyPID())
        return;

      double result = 0.;
      const Epetra_BlockMap& pnpmap = mysol_->Map();
      std::string position;
      res.ExtractString("QUANTITY",position);

      // test result value of single scalar field
      if (position=="pressure")
        result = (*mysol_)[pnpmap.LID(dis_->Dof(actnode,0))];
      // test result values for a system of scalars
      else
      {
        dserror("Quantity '%s' not supported in result-test of red_airway transport problems", position.c_str());
      }

      nerr += CompareValues(result, "NODE", res);
      test_count++;
    }
  }
}

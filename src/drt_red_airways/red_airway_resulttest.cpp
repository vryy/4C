/*----------------------------------------------------------------------*/
/*!
\file red_airway_resulttest.cpp

\brief testing of Red_Ariway calculation results

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
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
{
  dis_    = airways.Discretization();
  mysol_  = airways.Pnp();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void AIRWAY::RedAirwayResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{

  int dis;
  res.ExtractInt("DIS",dis);
  if (dis != 1)
    dserror("fix me: only one red_airway discretization supported for testing");

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
    const Epetra_BlockMap& pnpmap = mysol_->Map();
    std::string position;
    res.ExtractString("POSITION",position);

    // test result value of single scalar field
    if (position=="pressure")
      result = (*mysol_)[pnpmap.LID(dis_->Dof(actnode,0))];
    // test result values for a system of scalars
    else 
    {
      dserror("position '%s' not supported in result-test of red_airway transport problems", position.c_str());
    }

    nerr += CompareValues(result, res);
    test_count++;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool AIRWAY::RedAirwayResultTest::Match(DRT::INPUT::LineDefinition& res)
{
  return res.HaveNamed("RED_AIRWAY");
}



/*----------------------------------------------------------------------*/
/*!
\file scatra_resulttest.cpp

\brief testing of scalar transport calculation results

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_timint_implicit.H"
#include "scatra_timint_elch.H"
#include "scatra_resulttest.H"

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SCATRA::ScaTraResultTest::ScaTraResultTest(Teuchos::RCP<ScaTraTimIntImpl> scatra)
  : DRT::ResultTest("SCATRA")
{
  dis_    = scatra->Discretization();
  mysol_  = scatra->Phinp();
  myvan_  = scatra->Phiatmeshfreenodes();
  myflux_ = scatra->Flux();
  if(Teuchos::rcp_dynamic_cast<ScaTraTimIntElch>(scatra) == Teuchos::null)
    myelectkin_ = Teuchos::null;
  else
    myelectkin_ = Teuchos::rcp_dynamic_cast<ScaTraTimIntElch>(scatra)->OutputElectrodeInfoBoundary(false,false);
  mystrgrowth_ = scatra->StrGrowth();
  myfldgrowth_ = scatra->FldGrowth();
  mynumiter_ =scatra->IterNum();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScaTraResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
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
      const Epetra_BlockMap& phinpmap = mysol_->Map();
      std::string position;
      res.ExtractString("QUANTITY",position);

      // test result value of single scalar field
      if (position=="phi")
        result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,0))];
      // test result values for a system of scalars
      else if (position=="phi1")
        result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,0))];
      else if (position=="phi2")
        result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,1))];
      else if (position=="phi3")
        result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,2))];
      else if (position=="phi4")
        result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,3))];
      else if (position=="phi5")
        result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,4))];
      else if (position=="phi6")
        result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,5))];
      else if (position=="phi7")
        result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,6))];
      else if (position=="phi8")
        result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,7))];
      else if (position=="phi9")
        result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,8))];
      else if (position=="phi10")
        result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,9))];
      else if (position=="phi11")
        result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,10))];
      else if (position=="phi12")
        result = (*mysol_)[phinpmap.LID(dis_->Dof(0,actnode,11))];
      else if (position=="van_phi")
        result = (*myvan_)[phinpmap.LID(dis_->Dof(0,actnode,0))];
      // we support only testing of fluxes for the first scalar
      else if (position=="fluxx")
        result = ((*myflux_)[0])[phinpmap.LID(dis_->Dof(0,actnode,0))];
      else if (position=="fluxy")
        result = ((*myflux_)[1])[phinpmap.LID(dis_->Dof(0,actnode,0))];
      else if (position=="fluxz")
        result = ((*myflux_)[2])[phinpmap.LID(dis_->Dof(0,actnode,0))];
      else if (position=="meanc")
        result = (*myelectkin_)[0];
      else if (position=="meaneta")
        result = (*myelectkin_)[1];
      else if (position=="meancur")
        result = (*myelectkin_)[2];
      // test result values for biofilm growth (scatra structure and scatra fluid)
      else if (position=="scstr_growth_displx")
        result = ((*mystrgrowth_)[0])[phinpmap.LID(dis_->Dof(0,actnode,0))];
      else if (position=="scstr_growth_disply")
        result = ((*mystrgrowth_)[1])[phinpmap.LID(dis_->Dof(0,actnode,0))];
      else if (position=="scstr_growth_displz")
        result = ((*mystrgrowth_)[2])[phinpmap.LID(dis_->Dof(0,actnode,0))];
      else if (position=="scfld_growth_displx")
        result = ((*myfldgrowth_)[0])[phinpmap.LID(dis_->Dof(0,actnode,0))];
      else if (position=="scfld_growth_disply")
        result = ((*myfldgrowth_)[1])[phinpmap.LID(dis_->Dof(0,actnode,0))];
      else if (position=="scfld_growth_displz")
        result = ((*myfldgrowth_)[2])[phinpmap.LID(dis_->Dof(0,actnode,0))];
      else
      {
        dserror("Quantity '%s' not supported in result-test of scalar transport problems", position.c_str());
      }

      nerr += CompareValues(result, "NODE", res);
      test_count++;
    }
  }
}

void SCATRA::ScaTraResultTest::TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  std::string quantity;
  res.ExtractString("QUANTITY",quantity);
  bool unknownquantity = true; // make sure the result value std::string can be handled
  double result = 0.0;    // will hold the actual result of run

  // test for time step size
  if ( quantity == "numiterlastnewton" )
  {
    unknownquantity = false;
    result = (double)mynumiter_;
  }

  // catch quantity strings, which are not handled by scatra result test
  if ( unknownquantity )
    dserror("Quantity '%s' not supported in scatra testing", quantity.c_str());

  // compare values
  const int err = CompareValues(result, "SPECIAL", res);
  nerr += err;
  test_count++;

  return;
}

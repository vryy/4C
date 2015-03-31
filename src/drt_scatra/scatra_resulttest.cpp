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
SCATRA::ScaTraResultTest::ScaTraResultTest(Teuchos::RCP<ScaTraTimIntImpl> scatra) :
DRT::ResultTest("SCATRA"),
dis_(scatra->Discretization()),
mysol_(scatra->Phinp()),
myvan_(scatra->Phiatmeshfreenodes()),
myflux_(scatra->Flux()),
mystrgrowth_(scatra->StrGrowth()),
myfldgrowth_(scatra->FldGrowth()),
mynumiter_(scatra->IterNum())
{
  return;
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

      // extract name of quantity to be tested
      std::string quantity;
      res.ExtractString("QUANTITY",quantity);

      // get result to be tested
      const double result = ResultNode(quantity,actnode);

      nerr += CompareValues(result, "NODE", res);
      test_count++;
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | get nodal result to be tested                             fang 03/15 |
 *----------------------------------------------------------------------*/
const double SCATRA::ScaTraResultTest::ResultNode(
    const std::string   quantity,   //! name of quantity to be tested
    DRT::Node*          node        //! node carrying the result to be tested
    ) const
{
  // initialize variable for result
  double result(0.);

  // extract row map from solution vector
  const Epetra_BlockMap& phinpmap = mysol_->Map();

  // test result value of single scalar field
  if(quantity == "phi")
    result = (*mysol_)[phinpmap.LID(dis_->Dof(0,node,0))];

  // test result values for a system of scalars
  else if(quantity == "phi1")
    result = (*mysol_)[phinpmap.LID(dis_->Dof(0,node,0))];
  else if(quantity == "phi2")
    result = (*mysol_)[phinpmap.LID(dis_->Dof(0,node,1))];
  else if(quantity == "phi3")
    result = (*mysol_)[phinpmap.LID(dis_->Dof(0,node,2))];
  else if(quantity == "phi4")
    result = (*mysol_)[phinpmap.LID(dis_->Dof(0,node,3))];
  else if(quantity == "phi5")
    result = (*mysol_)[phinpmap.LID(dis_->Dof(0,node,4))];
  else if(quantity == "phi6")
    result = (*mysol_)[phinpmap.LID(dis_->Dof(0,node,5))];
  else if(quantity == "phi7")
    result = (*mysol_)[phinpmap.LID(dis_->Dof(0,node,6))];
  else if(quantity == "phi8")
    result = (*mysol_)[phinpmap.LID(dis_->Dof(0,node,7))];
  else if(quantity == "phi9")
    result = (*mysol_)[phinpmap.LID(dis_->Dof(0,node,8))];
  else if(quantity == "phi10")
    result = (*mysol_)[phinpmap.LID(dis_->Dof(0,node,9))];
  else if(quantity == "phi11")
    result = (*mysol_)[phinpmap.LID(dis_->Dof(0,node,10))];
  else if(quantity == "phi12")
    result = (*mysol_)[phinpmap.LID(dis_->Dof(0,node,11))];
  else if(quantity == "van_phi")
    result = (*myvan_)[phinpmap.LID(dis_->Dof(0,node,0))];

  // we support only testing of fluxes for the first scalar
  else if(quantity == "fluxx")
    result = ((*myflux_)[0])[phinpmap.LID(dis_->Dof(0,node,0))];
  else if(quantity == "fluxy")
    result = ((*myflux_)[1])[phinpmap.LID(dis_->Dof(0,node,0))];
  else if(quantity == "fluxz")
    result = ((*myflux_)[2])[phinpmap.LID(dis_->Dof(0,node,0))];

  // test result values for biofilm growth (scatra structure and scatra fluid)
  else if(quantity == "scstr_growth_displx")
    result = ((*mystrgrowth_)[0])[phinpmap.LID(dis_->Dof(0,node,0))];
  else if(quantity == "scstr_growth_disply")
    result = ((*mystrgrowth_)[1])[phinpmap.LID(dis_->Dof(0,node,0))];
  else if(quantity == "scstr_growth_displz")
    result = ((*mystrgrowth_)[2])[phinpmap.LID(dis_->Dof(0,node,0))];
  else if(quantity == "scfld_growth_displx")
    result = ((*myfldgrowth_)[0])[phinpmap.LID(dis_->Dof(0,node,0))];
  else if(quantity == "scfld_growth_disply")
    result = ((*myfldgrowth_)[1])[phinpmap.LID(dis_->Dof(0,node,0))];
  else if(quantity == "scfld_growth_displz")
    result = ((*myfldgrowth_)[2])[phinpmap.LID(dis_->Dof(0,node,0))];
  // catch unknown quantity strings
  else
    dserror("Quantity '%s' not supported in result test!", quantity.c_str());

  return result;
} // SCATRA::ScaTraResultTest::ResultNode


/*-------------------------------------------------------------------------------------*
 | test special quantity not associated with a particular element or node   fang 03/15 |
 *-------------------------------------------------------------------------------------*/
void SCATRA::ScaTraResultTest::TestSpecial(
    DRT::INPUT::LineDefinition&   res,
    int&                          nerr,
    int&                          test_count
    )
{
  // make sure that quantity is tested only once
  if(dis_->Comm().MyPID() == 0)
  {
    // extract name of quantity to be tested
    std::string quantity;
    res.ExtractString("QUANTITY",quantity);

    // get result to be tested
    const double result = ResultSpecial(quantity);

    // compare values
    const int err = CompareValues(result,"SPECIAL",res);
    nerr += err;
    test_count++;
  }

  return;
}


/*----------------------------------------------------------------------*
 | get special result to be tested                           fang 03/15 |
 *----------------------------------------------------------------------*/
const double SCATRA::ScaTraResultTest::ResultSpecial(
    const std::string   quantity   //! name of quantity to be tested
    ) const
{
  // initialize variable for result
  double result(0.);

  if(quantity == "numiterlastnewton")
    result = (double) mynumiter_;
  // catch unknown quantity strings
  else
    dserror("Quantity '%s' not supported in result test!", quantity.c_str());

  return result;
} // SCATRA::ScaTraResultTest::ResultSpecial

/*----------------------------------------------------------------------*/
/*!
\file scatra_resulttest.cpp

\brief testing of scalar transport calculation results

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
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
SCATRA::ScaTraResultTest::ScaTraResultTest(Teuchos::RCP<ScaTraTimIntImpl> scatratimint) :
DRT::ResultTest("SCATRA"),
scatratimint_(scatratimint)
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
  if (dis != scatratimint_->Discretization()->Name())
    return;

  int node;
  res.ExtractInt("NODE",node);
  node -= 1;

  int havenode(scatratimint_->Discretization()->HaveGlobalNode(node));
  int isnodeofanybody(0);
  scatratimint_->Discretization()->Comm().SumAll(&havenode,&isnodeofanybody,1);

  if (isnodeofanybody==0)
  {
    dserror("Node %d does not belong to discretization %s",node+1,scatratimint_->Discretization()->Name().c_str());
  }
  else
  {
    if (scatratimint_->Discretization()->HaveGlobalNode(node))
    {
      DRT::Node* actnode = scatratimint_->Discretization()->gNode(node);

      // Here we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->Owner() != scatratimint_->Discretization()->Comm().MyPID())
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
double SCATRA::ScaTraResultTest::ResultNode(
    const std::string   quantity,   //! name of quantity to be tested
    DRT::Node*          node        //! node carrying the result to be tested
    ) const
{
  // initialize variable for result
  double result(0.);

  // extract row map from solution vector
  const Epetra_BlockMap& phinpmap = scatratimint_->Phinp()->Map();

  // test result value of single scalar field
  if(quantity == "phi")
    result = (*scatratimint_->Phinp())[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,0))];

  // test result values for a system of scalars
  else if(quantity == "phi1")
    result = (*scatratimint_->Phinp())[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,0))];
  else if(quantity == "phi2")
    result = (*scatratimint_->Phinp())[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,1))];
  else if(quantity == "phi3")
    result = (*scatratimint_->Phinp())[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,2))];
  else if(quantity == "phi4")
    result = (*scatratimint_->Phinp())[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,3))];
  else if(quantity == "phi5")
    result = (*scatratimint_->Phinp())[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,4))];
  else if(quantity == "phi6")
    result = (*scatratimint_->Phinp())[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,5))];
  else if(quantity == "phi7")
    result = (*scatratimint_->Phinp())[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,6))];
  else if(quantity == "phi8")
    result = (*scatratimint_->Phinp())[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,7))];
  else if(quantity == "phi9")
    result = (*scatratimint_->Phinp())[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,8))];
  else if(quantity == "phi10")
    result = (*scatratimint_->Phinp())[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,9))];
  else if(quantity == "phi11")
    result = (*scatratimint_->Phinp())[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,10))];
  else if(quantity == "phi12")
    result = (*scatratimint_->Phinp())[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,11))];
  else if(quantity == "van_phi")
    result = (*scatratimint_->Phiatmeshfreenodes())[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,0))];

  // we support only testing of fluxes for the first scalar
  else if(quantity == "fluxx")
    result = ((*scatratimint_->Flux())[0])[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,0))];
  else if(quantity == "fluxy")
    result = ((*scatratimint_->Flux())[1])[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,0))];
  else if(quantity == "fluxz")
    result = ((*scatratimint_->Flux())[2])[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,0))];

  // test result values for biofilm growth (scatra structure and scatra fluid)
  else if(quantity == "scstr_growth_displx")
    result = ((*scatratimint_->StrGrowth())[0])[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,0))];
  else if(quantity == "scstr_growth_disply")
    result = ((*scatratimint_->StrGrowth())[1])[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,0))];
  else if(quantity == "scstr_growth_displz")
    result = ((*scatratimint_->StrGrowth())[2])[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,0))];
  else if(quantity == "scfld_growth_displx")
    result = ((*scatratimint_->FldGrowth())[0])[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,0))];
  else if(quantity == "scfld_growth_disply")
    result = ((*scatratimint_->FldGrowth())[1])[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,0))];
  else if(quantity == "scfld_growth_displz")
    result = ((*scatratimint_->FldGrowth())[2])[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,0))];

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
  // make sure that quantity is tested only on specified discretization and only by one processor
  std::string disname;
  res.ExtractString("DIS",disname);
  if(disname == scatratimint_->Discretization()->Name() and scatratimint_->Discretization()->Comm().MyPID() == 0)
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
double SCATRA::ScaTraResultTest::ResultSpecial(
    const std::string   quantity   //! name of quantity to be tested
    ) const
{
  // initialize variable for result
  double result(0.);

  // number of Newton-Raphson iterations in last time step
  if(quantity == "numiterlastnewton")
    result = (double) scatratimint_->IterNum();

  // total or mean value of transported scalar in subdomain or entire domain
  else if(!quantity.compare(0,8,"totalphi") or !quantity.compare(0,7,"meanphi"))
  {
    // examples of possible specifications:
    // 1.) quantity == "{total,mean}phi", suffix == "": total or mean value of transported scalar with ID 0 in entire domain (with ID -1)
    // 2.) quantity == "{total,mean}phia", suffix == "a": total or mean value of transported scalar with ID a in entire domain (with ID -1)
    // 3.) quantity == "{total,mean}phi_b", suffix == "_b": total or mean value of transported scalar with ID 0 in subdomain with ID b
    // 4.) quantity == "{total,mean}phia_b", suffix == "a_b": total or mean value of transported scalar with ID a in subdomain with ID b
    // other specifications are invalid and result in an dserror
    std::string suffix;
    if(!quantity.compare(0,8,"totalphi"))
      suffix = quantity.substr(8);
    else
      suffix = quantity.substr(7);

    const char* index(suffix.c_str());
    char* locator(NULL);
    unsigned species(0), domain(-1);

    if(*index == '\0')
    {
      // case 1 from list above: do nothing
    }

    else
    {
      if(*index != '_')
      {
        // cases 2 and 4 from list above
        // extract species ID
        species = strtol(index,&locator,10)-1;
        if(locator == index)
          dserror("Couldn't read species ID!");

        // move locator to position behind species ID
        index = locator;
      }

      if(*index == '\0')
      {
        // case 2 from list above: do nothing
      }

      else if(*index == '_')
      {
        // cases 3 and 4 from list above
        // move locator to domain ID, i.e., to position behind underscore
        ++index;

        // extract domain ID
        domain = strtol(index,&locator,10);
        if(locator == index)
          dserror("Couldn't read domain ID!");
      }

      // catch invalid syntax
      else
        dserror("Wrong syntax!");
    }

    // extract map with relevant result from scalar transport time integrator
    const std::map<const int,std::vector<double> >* map(NULL);
    if(!quantity.compare(0,8,"totalphi"))
      map = &scatratimint_->TotalScalars();
    else
      map = &scatratimint_->MeanScalars();

    // extract relevant result from map
    std::map<const int,std::vector<double> >::const_iterator iterator = map->find(domain);
    if(iterator == map->end() or species < 0 or species >= iterator->second.size())
      dserror("Couldn't extract total or mean value of transported scalar with ID %d inside domain with ID %d from time integrator!",species,domain);
    result = iterator->second[species];
  }

  // catch unknown quantity strings
  else
    dserror("Quantity '%s' not supported in result test!", quantity.c_str());

  return result;
} // SCATRA::ScaTraResultTest::ResultSpecial

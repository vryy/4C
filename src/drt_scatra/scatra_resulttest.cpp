/*----------------------------------------------------------------------*/
/*!
\file scatra_resulttest.cpp

\brief testing of scalar transport calculation results

\level 1

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#include "scatra_resulttest.H"
#include "scatra_timint_implicit.H"
#include "scatra_timint_meshtying_strategy_s2i.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_linedefinition.H"

#include "../linalg/linalg_solver.H"

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

  // test result value for a system of scalars
  else if(!quantity.compare(0,3,"phi"))
  {
    // read species ID
    std::string k_string = quantity.substr(3);
    char* locator(NULL);
    int k = strtol(k_string.c_str(),&locator,10) - 1;

    // safety checks
    if(locator == k_string.c_str())
      dserror("Couldn't read species ID!");
    if(scatratimint_->Discretization()->NumDof(0,node)<=k)
      dserror("Species ID is larger than number of DOFs of node!");

    // extract result
    result = (*scatratimint_->Phinp())[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,k))];
  }

  // test domain or boundary flux
  else if(!quantity.compare(0,12,"flux_domain_") or !quantity.compare(0,14,"flux_boundary_"))
  {
    // read species ID
    std::string suffix;
    if(!quantity.compare(0,12,"flux_domain_"))
      suffix = quantity.substr(12);
    else
      suffix = quantity.substr(14);
    char* locator(NULL);
    int k = strtol(suffix.c_str(),&locator,10) - 1;

    // safety checks
    if(locator == suffix.c_str())
      dserror("Couldn't read species ID!");
    if(scatratimint_->Discretization()->NumDof(0,node)<=k)
      dserror("Species ID is larger than number of DOFs of node!");

    // read spatial dimension
    ++locator;
    unsigned dim(-1);
    if(*locator == 'x')
      dim = 0;
    else if(*locator == 'y')
      dim = 1;
    else if(*locator == 'z')
      dim = 2;
    else
      dserror("Invalid syntax!");

    // safety check
    if(*(++locator) != '\0')
      dserror("Invalid syntax!");

    // extract result
    if(!quantity.compare(0,12,"flux_domain_"))
      result = ((*scatratimint_->FluxDomain())[dim])[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,k))];
    else
      result = ((*scatratimint_->FluxBoundary())[dim])[phinpmap.LID(scatratimint_->Discretization()->Dof(0,node,k))];
  }

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

  // test scatra-scatra interface layer thickness
  else if(quantity == "s2ilayerthickness")
  {
    // extract scatra-scatra interface meshtying strategy class
    const Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> strategy = Teuchos::rcp_dynamic_cast<const SCATRA::MeshtyingStrategyS2I>(scatratimint_->Strategy());
    if(strategy == Teuchos::null)
      dserror("Couldn't extract scatra-scatra interface meshtying strategy class!");

    // extract state vector of discrete scatra-scatra interface layer thicknesses
    // depending on whether monolithic or semi-implicit solution approach is used
    Teuchos::RCP<const Epetra_Vector> s2igrowthvec(Teuchos::null);
    switch(strategy->IntLayerGrowthEvaluation())
    {
      case INPAR::S2I::growth_evaluation_monolithic:
      {
        s2igrowthvec = strategy->GrowthVarNp();
        break;
      }

      case INPAR::S2I::growth_evaluation_semi_implicit:
      {
        s2igrowthvec = strategy->GrowthVarN();
        break;
      }

      default:
      {
        dserror("Can't test scatra-scatra interface layer thickness!");
        break;
      }
    }

    // safety check
    if(s2igrowthvec == Teuchos::null)
      dserror("Couldn't extract state vector of discrete scatra-scatra interface layer thicknesses!");

    // extract result
    result = (*s2igrowthvec)[scatratimint_->Discretization()->DofRowMap(2)->LID(scatratimint_->Discretization()->Dof(2,node,0))];
  }

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
  // make sure that quantity is tested only on specified discretization
  std::string disname;
  res.ExtractString("DIS",disname);
  if(disname == scatratimint_->Discretization()->Name())
  {
    // extract name of quantity to be tested
    std::string quantity;
    res.ExtractString("QUANTITY",quantity);

    // get result to be tested on all processors
    const double result = ResultSpecial(quantity);

    // compare values on first processor
    if(scatratimint_->Discretization()->Comm().MyPID() == 0)
    {
      const int err = CompareValues(result,"SPECIAL",res);
      nerr += err;
      test_count++;
    }
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

  // total values, mean values, relative L2 errors, or relative H1 errors of scalar fields in subdomain or entire domain
  else if(!quantity.compare(0,8,"totalphi") or !quantity.compare(0,7,"meanphi") or !quantity.compare(0,7,"L2error") or !quantity.compare(0,7,"H1error"))
  {
    // examples of possible specifications:
    // 1.) quantity == "{{total,mean}phi,{L2,H1}error}", suffix == "": total value, mean value, relative L2 error, or relative H1 error of scalar with ID 0 in entire domain (with ID -1)
    // 2.) quantity == "{{total,mean}phi,{L2,H1}error}a", suffix == "a": total value, mean value, relative L2 error, or relative H1 error of scalar with ID a in entire domain (with ID -1)
    // 3.) quantity == "{{total,mean}phi,{L2,H1}error}_b", suffix == "_b": total value, mean value, relative L2 error, or relative H1 error of scalar with ID 0 in subdomain with ID b
    // 4.) quantity == "{{total,mean}phi,{L2,H1}error}a_b", suffix == "a_b": total value, mean value, relative L2 error, or relative H1 error of scalar with ID a in subdomain with ID b
    // other specifications are invalid and result in a dserror
    std::string suffix;
    if(!quantity.compare(0,8,"totalphi"))
      suffix = quantity.substr(8);
    else
      suffix = quantity.substr(7);

    const char* index(suffix.c_str());
    char* locator(NULL);
    int species(0), domain(-1);

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

        // move index to position behind species ID
        index = locator;
      }

      if(*index == '\0')
      {
        // case 2 from list above: do nothing
      }

      else if(*index == '_')
      {
        // cases 3 and 4 from list above
        // move index to domain ID, i.e., to position behind underscore
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

    // total or mean value
    if(!quantity.compare(0,8,"totalphi") or !quantity.compare(0,7,"meanphi"))
    {
      // extract map with relevant result from scalar transport time integrator
      const std::map<const int,std::vector<double> >* map(NULL);
      if(!quantity.compare(0,8,"totalphi"))
        map = &scatratimint_->TotalScalars();
      else
        map = &scatratimint_->MeanScalars();

      // extract relevant result from map
      std::map<const int,std::vector<double> >::const_iterator iterator = map->find(domain);
      if(iterator == map->end() or species < 0 or species >= (int) iterator->second.size())
        dserror("Couldn't extract total or mean value of transported scalar with ID %d inside domain with ID %d from time integrator!",species,domain);
      result = iterator->second[species];
    }

    // relative L2 or H1 error
    else
    {
      // global error
      if(domain == -1)
      {
        if(!quantity.compare(0,7,"L2error"))
          result = (*scatratimint_->RelErrors())[species*2];
        else
          result = (*scatratimint_->RelErrors())[species*2+1];
      }

      // error inside subdomain
      else
      {
        if(!quantity.compare(0,7,"L2error"))
          result = (*scatratimint_->RelErrors())[domain*scatratimint_->NumDofPerNode()*2+species*2];
        else
          result = (*scatratimint_->RelErrors())[domain*scatratimint_->NumDofPerNode()*2+species*2+1];
      }
    }
  }

  // number of iterations performed by linear solver during last Newton-Raphson iteration
  else if(quantity == "numiterlastsolve")
  {
    // safety check
    if(scatratimint_->Solver()->Params().get("solver","none") != "aztec")
      dserror("Must have Aztec solver for result test involving number of solver iterations during last Newton-Raphson iteration!");
    result = (double) scatratimint_->Solver()->getNumIters();
  }

  // test parallel distribution of scatra-scatra coupling interface
  else if(!quantity.compare(0,9,"s2inumdof"))
  {
    // extract string part behind "s2inumdof"
    std::string suffix = quantity.substr(9);

    // check syntax
    if(suffix.compare(0,4,"_int"))
      dserror("Wrong syntax!");

    // initialize auxiliary variables
    const char* index(suffix.substr(4).c_str());
    char* locator(NULL);

    // extract interface ID
    const int interface = strtol(index,&locator,10);
    if(locator == index)
      dserror("Couldn't read interface ID!");

    // extract string part behind interface ID
    suffix.assign(locator);

    // check syntax
    if(suffix.compare(0,5,"_proc"))
      dserror("Wrong syntax!");

    // move index to processor ID
    index = suffix.substr(5).c_str();

    // extract processor ID
    const int processor = strtol(index,&locator,10);
    if(locator == index)
      dserror("Couldn't read processor ID!");
    if(processor >= scatratimint_->Discretization()->Comm().NumProc())
      dserror("Invalid processor ID!");

    // check syntax
    if(*locator != '\0')
      dserror("Wrong syntax!");

    // extract scatra-scatra interface meshtying strategy class
    const Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> strategy = Teuchos::rcp_dynamic_cast<const SCATRA::MeshtyingStrategyS2I>(scatratimint_->Strategy());
    if(strategy == Teuchos::null)
      dserror("Couldn't extract scatra-scatra interface meshtying strategy class!");

    // extract number of degrees of freedom owned by specified processor at specified scatra-scatra coupling interface
    result = strategy->MortarDiscretization(interface).DofRowMap()->NumMyElements();
    scatratimint_->Discretization()->Comm().Broadcast(&result,1,processor);
  }

  // catch unknown quantity strings
  else
    dserror("Quantity '%s' not supported in result test!", quantity.c_str());

  return result;
} // SCATRA::ScaTraResultTest::ResultSpecial

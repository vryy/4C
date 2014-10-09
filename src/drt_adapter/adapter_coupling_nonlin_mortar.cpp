/*!----------------------------------------------------------------------
\file adapter_coupling_nonlin_mortar.cpp

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  headers                                                  farah 10/14|
 *----------------------------------------------------------------------*/
#include "adapter_coupling_nonlin_mortar.H"

#include "../drt_contact/contact_interface.H"
#include "../drt_contact/contact_node.H"
#include "../drt_contact/contact_element.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 |  ctor                                                     farah 10/14|
 *----------------------------------------------------------------------*/
ADAPTER::CouplingNonLinMortar::CouplingNonLinMortar()
{
  //empty...
}


/*----------------------------------------------------------------------*
 |  setup for nonlinear mortar framework                     farah 10/14|
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingNonLinMortar::Setup(
    Teuchos::RCP<DRT::Discretization>   masterdis,
    Teuchos::RCP<DRT::Discretization>   slavedis,
    std::vector<int>                    coupleddof,
    const std::string&                  couplingcond)
{
  //TODO: extend this to sliding ale + ALE-dis
  // vector coupleddof defines degree of freedom which are coupled (1: coupled; 0: not coupled), e.g.:
  // - fluid 3D meshtying: coupleddof = [1, 1, 1, 1] -> all degrees of freedom (velocity and pressure) are coupled
  // - fluid 3D meshtying: coupleddof = [1, 1, 1, 0] -> only velocity degrees of freedom are coupled
  // - fsi 3D: coupleddof = [1, 1, 1] -> at the interface only displacements are coupled
  // - ....

  myrank_= masterdis->Comm().MyPID();
  comm_  = Teuchos::rcp(masterdis->Comm().Clone());

  // initialize maps for row nodes
  std::map<int, DRT::Node*> masternodes;
  std::map<int, DRT::Node*> slavenodes;

  // initialize maps for column nodes
  std::map<int, DRT::Node*> mastergnodes;
  std::map<int, DRT::Node*> slavegnodes;

  //initialize maps for elements
  std::map<int, Teuchos::RCP<DRT::Element> > masterelements;
  std::map<int, Teuchos::RCP<DRT::Element> > slaveelements;

  // Coupling condition is defined by "MORTAR COUPLING CONDITIONS"
  // There is only one discretization (masterdis == slavedis). Therefore, the node set have to be
  // separated beforehand.
  if(couplingcond=="Mortar")
  {
    std::vector<DRT::Condition*> conds;
    std::vector<DRT::Condition*> conds_master(0);
    std::vector<DRT::Condition*> conds_slave(0);
    masterdis->GetCondition(couplingcond, conds);

    for (unsigned i=0; i<conds.size(); i++)
    {
      const std::string* side = conds[i]->Get<std::string>("Side");

      if (*side == "Master")
        conds_master.push_back(conds[i]);

      if(*side == "Slave")
        conds_slave.push_back(conds[i]);
    }

    // Fill maps based on condition for master side (masterdis == slavedis)
    DRT::UTILS::FindConditionObjects(*masterdis, masternodes, mastergnodes, masterelements, conds_master);

    // Fill maps based on condition for slave side (masterdis == slavedis)
    DRT::UTILS::FindConditionObjects(*slavedis, slavenodes, slavegnodes, slaveelements, conds_slave);
  }
  // Coupling condition is defined by "FSI COUPLING CONDITIONS"
  // There are two discretizations for the master and slave side. Therefore, the master/slave nodes
  // are chosen based on the discretization.
  else
  {
    // Fill maps based on condition for master side (masterdis != slavedis)
    DRT::UTILS::FindConditionObjects(*masterdis, masternodes, mastergnodes, masterelements, couplingcond);

    // Fill maps based on condition for slave side (masterdis != slavedis)
    DRT::UTILS::FindConditionObjects(*slavedis, slavenodes, slavegnodes, slaveelements, couplingcond);
  }

  // get mortar coupling parameters
  const Teuchos::ParameterList& inputmortar = DRT::Problem::Instance()->MortarCouplingParams();
  Teuchos::ParameterList input;
  input.setParameters(inputmortar);

  // is this a nurbs problem?
  std::string distype = DRT::Problem::Instance()->SpatialApproximation();
  if(distype=="Nurbs")
  {
    // ***
    dserror("nurbs for fsi mortar not supported!");
    input.set<bool>("NURBS",true);
  }
  else
    input.set<bool>("NURBS",false);

  // check for invalid parameter values
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(input,"SHAPEFCN") != INPAR::MORTAR::shape_dual)
    if(myrank_== 0) dserror("Mortar coupling adapter only works for dual shape functions");
  if (DRT::INPUT::IntegralValue<int>(input,"LM_NODAL_SCALE")==true)
    if(myrank_== 0) dserror("Mortar coupling adapter does not work with LM_NODAL_SCALE");

  // check for parallel redistribution (only if more than 1 proc)
  bool parredist = false;
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(input,"PARALLEL_REDIST")
      != INPAR::MORTAR::parredist_none)
    if (comm_->NumProc()>1) parredist = true;

  // get problem dimension (2D or 3D) and create (MORTAR::MortarInterface)
  const int dim = DRT::Problem::Instance()->NDim();

  // create an empty mortar interface
  // (To be on the safe side we still store all interface nodes and elements
  // fully redundant here in the mortar ADAPTER. This makes applications such
  // as SlidingALE much easier, whereas it would not be needed for others.)
  INPAR::MORTAR::RedundantStorage redundant =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::RedundantStorage>(input,"REDUNDANT_STORAGE");

  Teuchos::RCP<CONTACT::CoInterface> interface =
      Teuchos::rcp(new CONTACT::CoInterface(0, *comm_, dim, input,false, redundant));

  // number of dofs per node based on the coupling vector coupleddof
  int dof = coupleddof.size();
//  if((masterdis->NumDof(masterdis->lRowNode(0))!=dof and slavewithale==true and slidingale==false) or
//      (slavedis->NumDof(slavedis->lRowNode(0))!=dof and slavewithale==false and slidingale==false))
//  {
//    dserror("The size of the coupling vector coupleddof and dof defined in the discretization does not fit!! \n"
//            "dof defined in the discretization: %i \n"
//            "length of coupleddof: %i",masterdis->NumDof(masterdis->lRowNode(0)), dof);
//  }

  // special case: sliding ale
  // In the sliding ale framework two mortar discretizations are generated from identical
  // masterelement and slaveelement sets. Since node-, dof- and element ids of the original elements are
  // the same, an offset have to be defined
  int nodeoffset=0;
  int dofoffset=0;
//  if(slidingale==true)
//  {
//    nodeoffset = masterdis->NodeRowMap()->MaxAllGID()+1;
//    dofoffset = masterdis->DofRowMap()->MaxAllGID()+1;
//  }

  // number of coupled dofs (defined in coupleddof by a 1)
  int numcoupleddof = 0;
  for(int ii=0; ii<dof; ++ii)
    if(coupleddof[ii]==1) numcoupleddof+=1;

  // feeding master nodes to the interface including ghosted nodes
  std::map<int, DRT::Node*>::const_iterator nodeiter;
  for (nodeiter = mastergnodes.begin(); nodeiter != mastergnodes.end(); ++nodeiter)
  {
    DRT::Node* node = nodeiter->second;
    // vector containing only the gids of the coupled dofs (size numcoupleddof)
    std::vector<int> dofids(numcoupleddof);
    int ii=0;
    for (int k=0;k<dof;++k)
    {
      // Should this dof be coupled? (==1),
      if (coupleddof[k]==1)
      {
        // get the gid of the coupled dof (size dof)
        // and store it in the vector dofids containing only coupled dofs (size numcoupleddof)
        dofids[ii] = masterdis->Dof(node)[k];
        ii +=1;
      }
    }
    Teuchos::RCP<CONTACT::CoNode> cnode = Teuchos::rcp(
                new CONTACT::CoNode(node->Id(), node->X(), node->Owner(),
                    numcoupleddof, dofids, false,false));

    interface->AddCoNode(cnode);
  }

  // feeding slave nodes to the interface including ghosted nodes
  for (nodeiter = slavegnodes.begin(); nodeiter != slavegnodes.end(); ++nodeiter)
  {
    DRT::Node* node = nodeiter->second;
    // vector containing only the gids of the coupled dofs (size numcoupleddof)
    std::vector<int> dofids(numcoupleddof);
    int ii=0;
    for (int k=0;k<dof;++k)
    {
      // Should this dof be coupled? (==1)
      if (coupleddof[k]==1)
      {
        // get the gid of the coupled dof (size dof)
        // and store it in the vector dofids containing only coupled dofs (size numcoupleddof)
        dofids[ii] = slavedis->Dof(node)[k]+dofoffset;
        ii += 1;
      }
    }
    Teuchos::RCP<CONTACT::CoNode> cnode = Teuchos::rcp(
                new CONTACT::CoNode(node->Id(), node->X(), node->Owner(),
                    numcoupleddof, dofids, true,true));

    interface->AddCoNode(cnode);
  }

  // We need to determine an element offset to start the numbering of the slave
  // mortar elements AFTER the master mortar elements in order to ensure unique
  // eleIDs in the interface discretization. The element offset equals the
  // overall number of master mortar elements (which is not equal to the number
  // of elements in the field that is chosen as master side).
  //
  // If masterdis==slavedis, the element numbering is right without offset
  int eleoffset = 0;
  if(masterdis.get()!=slavedis.get())
  {
    int nummastermtreles = masterelements.size();
    comm_->SumAll(&nummastermtreles,&eleoffset,1);
  }

//  if(slidingale==true)
//    eleoffset = masterdis->ElementRowMap()->MaxAllGID()+1;

  // feeding master elements to the interface
  std::map<int, Teuchos::RCP<DRT::Element> >::const_iterator elemiter;
  for (elemiter = masterelements.begin(); elemiter != masterelements.end(); ++elemiter)
  {
    Teuchos::RCP<DRT::Element> ele = elemiter->second;
    Teuchos::RCP<CONTACT::CoElement> cele = Teuchos::rcp(
                new CONTACT::CoElement(ele->Id(), ele->Owner(), ele->Shape(),
                    ele->NumNode(), ele->NodeIds(), false));

    interface->AddCoElement(cele);
  }

  // feeding slave elements to the interface
  for (elemiter = slaveelements.begin(); elemiter != slaveelements.end(); ++elemiter)
  {
    Teuchos::RCP<DRT::Element> ele = elemiter->second;

    // Here, we have to distinguish between standard and sliding ale since mortar elements are generated
    // from the identical element sets in the case of sliding ale
    // Therefore, we introduce an element offset AND a node offset for the the slave mortar elements
    if(true)//(slidingale==false)
    {
      Teuchos::RCP<CONTACT::CoElement> cele = Teuchos::rcp(
                  new CONTACT::CoElement(ele->Id(), ele->Owner(), ele->Shape(),
                      ele->NumNode(), ele->NodeIds(), true));

      interface->AddCoElement(cele);
    }
    else
    {
      std::vector<int> nidsoff;
      for(int i=0; i<ele->NumNode(); i++)
      {
        nidsoff.push_back(ele->NodeIds()[ele->NumNode()-1-i]+nodeoffset);
      }

      Teuchos::RCP<CONTACT::CoElement> cele = Teuchos::rcp(
                  new CONTACT::CoElement(ele->Id() + eleoffset, ele->Owner(), ele->Shape(),
                      ele->NumNode(), &(nidsoff[0]), true));

      interface->AddCoElement(cele);
    }
  }

  // finalize the contact interface construction
  interface->FillComplete();

  // store old row maps (before parallel redistribution)
  slavedofrowmap_  = Teuchos::rcp(new Epetra_Map(*interface->SlaveRowDofs()));
  masterdofrowmap_ = Teuchos::rcp(new Epetra_Map(*interface->MasterRowDofs()));

  // print parallel distribution
  interface->PrintParallelDistribution(1);

  //**********************************************************************
  // PARALLEL REDISTRIBUTION OF INTERFACE
  //**********************************************************************
  if (parredist && comm_->NumProc()>1)
  {
    // redistribute optimally among all procs
    interface->Redistribute(1);

    // call fill complete again
    interface->FillComplete();

    // print parallel distribution again
    interface->PrintParallelDistribution(1);
  }

  // store interface
  interface_ = interface;

  return;
}





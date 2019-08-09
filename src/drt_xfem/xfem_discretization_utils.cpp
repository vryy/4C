/*----------------------------------------------------------------------*/
/*!

\brief Basic discretization-related tools used in XFEM routines

\level 1

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
*/
/*----------------------------------------------------------------------*/

#include "xfem_discretization_utils.H"

#include "../drt_lib/drt_discret_faces.H"
#include "../drt_lib/drt_discret_xfem.H"

#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_utils_rebalancing.H"
#include "../drt_lib/drt_dofset_fixed_size.H"

#include "../linalg/linalg_utils.H"

#include "../drt_io/io_gmsh.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::UTILS::PrintDiscretizationToStream(Teuchos::RCP<DRT::Discretization> dis,
    const std::string& disname, bool elements, bool elecol, bool nodes, bool nodecol, bool faces,
    bool facecol, std::ostream& s, std::map<int, LINALG::Matrix<3, 1>>* curr_pos)
{
  if (elements)
  {
    // draw bg elements with associated gid
    s << "View \" " << disname;
    if (elecol)
    {
      s << " col e->Id() \" {\n";
      for (int i = 0; i < dis->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = dis->lColElement(i);
        if (curr_pos == NULL)
          IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
        else
          IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
      };
    }
    else
    {
      s << " row e->Id() \" {\n";
      for (int i = 0; i < dis->NumMyRowElements(); ++i)
      {
        const DRT::Element* actele = dis->lRowElement(i);
        if (curr_pos == NULL)
          IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
        else
          IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
      };
    }
    s << "};\n";
  }

  if (nodes)
  {
    s << "View \" " << disname;
    if (nodecol)
    {
      s << " col n->Id() \" {\n";
      for (int i = 0; i < dis->NumMyColNodes(); ++i)
      {
        const DRT::Node* actnode = dis->lColNode(i);
        LINALG::Matrix<3, 1> pos(true);

        if (curr_pos != NULL)
        {
          const LINALG::Matrix<3, 1>& curr_x = curr_pos->find(actnode->Id())->second;
          pos(0) = curr_x(0);
          pos(1) = curr_x(1);
          pos(2) = curr_x(2);
        }
        else
        {
          const LINALG::Matrix<3, 1> x(actnode->X());
          pos(0) = x(0);
          pos(1) = x(1);
          pos(2) = x(2);
        }
        IO::GMSH::cellWithScalarToStream(DRT::Element::point1, actnode->Id(), pos, s);
      }
    }
    else
    {
      s << " row n->Id() \" {\n";
      for (int i = 0; i < dis->NumMyRowNodes(); ++i)
      {
        const DRT::Node* actnode = dis->lRowNode(i);
        LINALG::Matrix<3, 1> pos(true);

        if (curr_pos != NULL)
        {
          const LINALG::Matrix<3, 1>& curr_x = curr_pos->find(actnode->Id())->second;
          pos(0) = curr_x(0);
          pos(1) = curr_x(1);
          pos(2) = curr_x(2);
        }
        else
        {
          const LINALG::Matrix<3, 1> x(actnode->X());
          pos(0) = x(0);
          pos(1) = x(1);
          pos(2) = x(2);
        }
        IO::GMSH::cellWithScalarToStream(DRT::Element::point1, actnode->Id(), pos, s);
      }
    }
    s << "};\n";
  }

  if (faces)
  {
    // cast to DiscretizationXFEM
    Teuchos::RCP<DRT::DiscretizationFaces> xdis =
        Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(dis, true);
    if (xdis == Teuchos::null)
      dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");

    s << "View \" " << disname;

    if (xdis->FilledExtension() == true)  // faces output
    {
      if (facecol)
      {
        s << " col f->Id() \" {\n";
        for (int i = 0; i < xdis->NumMyColFaces(); ++i)
        {
          const DRT::Element* actele = xdis->lColFace(i);
          if (curr_pos == NULL)
            IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
          else
            IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
        };
      }
      else
      {
        s << " row f->Id() \" {\n";
        for (int i = 0; i < xdis->NumMyRowFaces(); ++i)
        {
          const DRT::Element* actele = xdis->lRowFace(i);
          if (curr_pos == NULL)
            IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
          else
            IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
        };
      }
      s << "};\n";
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::UTILS::XFEMDiscretizationBuilder::SetupXFEMDiscretization(
    const Teuchos::ParameterList& xgen_params, Teuchos::RCP<DRT::Discretization> dis,
    int numdof) const
{
  Teuchos::RCP<DRT::DiscretizationXFEM> xdis =
      Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(dis, false);
  //
  if (xdis == Teuchos::null)
  {
    dserror("No XFEM discretization for XFEM problem available!");

    // REMARK: standard fluid could also step into this routine, as a special case! (remove dserror)
    if (!dis->Filled()) dis->FillComplete();

    return;
  }

  if (!xdis->Filled()) xdis->FillComplete();

  const Epetra_Map* noderowmap = xdis->NodeRowMap();
  if (noderowmap == NULL) dserror("we expect a fill-complete call before!");

  // now we can reserve dofs for xfem discretization
  int nodeindexrange =
      noderowmap->MaxAllGID() - noderowmap->MinAllGID() + 1;  // if id's are not continuous numbered
  int maxNumMyReservedDofsperNode = (xgen_params.get<int>("MAX_NUM_DOFSETS")) * numdof;
  Teuchos::RCP<DRT::FixedSizeDofSet> maxdofset =
      Teuchos::rcp(new DRT::FixedSizeDofSet(maxNumMyReservedDofsperNode, nodeindexrange));

  const int fluid_nds = 0;
  xdis->ReplaceDofSet(fluid_nds, maxdofset, true);  // fluid dofset has nds = 0
  std::vector<int> nds;
  nds.push_back(fluid_nds);
  xdis->InitialFillComplete(nds);

  // print all dofsets
  xdis->GetDofSetProxy()->PrintAllDofsets(xdis->Comm());

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::UTILS::XFEMDiscretizationBuilder::SetupXFEMDiscretization(
    const Teuchos::ParameterList& xgen_params, Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<DRT::Discretization> embedded_dis, const std::string& embedded_cond_name,
    int numdof) const
{
  if (!embedded_dis->Filled()) embedded_dis->FillComplete();

  Teuchos::RCP<DRT::DiscretizationXFEM> xdis =
      Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(dis, true);
  if (!xdis->Filled()) xdis->FillComplete();

  // get fluid mesh conditions: hereby we specify standalone embedded discretizations
  std::vector<DRT::Condition*> conditions;
  xdis->GetCondition(embedded_cond_name, conditions);

  std::vector<std::string> conditions_to_copy;
  xdis->GetConditionNames(conditions_to_copy);

  SplitDiscretizationByCondition(xdis, embedded_dis, conditions, conditions_to_copy);

  SetupXFEMDiscretization(xgen_params, xdis, numdof);

  DRT::UTILS::PrintParallelDistribution(*dis);
  DRT::UTILS::PrintParallelDistribution(*embedded_dis);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XFEM::UTILS::XFEMDiscretizationBuilder::SetupXFEMDiscretization(
    const Teuchos::ParameterList& xgen_params, Teuchos::RCP<DRT::DiscretizationInterface> src_dis,
    Teuchos::RCP<DRT::DiscretizationInterface> target_dis,
    const std::vector<DRT::Condition*>& boundary_conds) const
{
  Teuchos::RCP<DRT::Discretization> src_dis_ptr =
      Teuchos::rcp_dynamic_cast<DRT::Discretization>(src_dis, true);
  Teuchos::RCP<DRT::Discretization> target_dis_ptr =
      Teuchos::rcp_dynamic_cast<DRT::Discretization>(target_dis, true);
  return SetupXFEMDiscretization(xgen_params, src_dis_ptr, target_dis_ptr, boundary_conds);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XFEM::UTILS::XFEMDiscretizationBuilder::SetupXFEMDiscretization(
    const Teuchos::ParameterList& xgen_params, Teuchos::RCP<DRT::Discretization> src_dis,
    Teuchos::RCP<DRT::Discretization> target_dis,
    const std::vector<DRT::Condition*>& boundary_conds) const
{
  if (!target_dis->Filled()) target_dis->FillComplete();

  if (!src_dis->Filled()) src_dis->FillComplete();

  // get the number of DoF's per node
  int gid_node = src_dis->NodeRowMap()->MinMyGID();
  DRT::Node* node_ptr = src_dis->gNode(gid_node);
  int num_dof_per_node = src_dis->NumDof(node_ptr);

  std::vector<std::string> conditions_to_copy;
  src_dis->GetConditionNames(conditions_to_copy);

  SplitDiscretizationByBoundaryCondition(src_dis, target_dis, boundary_conds, conditions_to_copy);

  if (!Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(src_dis).is_null())
    SetupXFEMDiscretization(xgen_params, src_dis, num_dof_per_node);
  if (!Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(target_dis).is_null())
    SetupXFEMDiscretization(xgen_params, target_dis, num_dof_per_node);

  DRT::UTILS::PrintParallelDistribution(*src_dis);
  DRT::UTILS::PrintParallelDistribution(*target_dis);

  return num_dof_per_node;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::UTILS::XFEMDiscretizationBuilder::SplitDiscretizationByCondition(
    Teuchos::RCP<DRT::Discretization> sourcedis, Teuchos::RCP<DRT::Discretization> targetdis,
    std::vector<DRT::Condition*>& conditions,
    const std::vector<std::string>& conditions_to_copy) const
{
  // row node map (id -> pointer)
  std::map<int, DRT::Node*> sourcenodes;

  // column node map
  std::map<int, DRT::Node*> sourcegnodes;

  // element map
  std::map<int, Teuchos::RCP<DRT::Element>> sourceelements;

  // find conditioned nodes (owned and ghosted) and elements
  DRT::UTILS::FindConditionObjects(
      *sourcedis, sourcenodes, sourcegnodes, sourceelements, conditions);

  SplitDiscretization(
      sourcedis, targetdis, sourcenodes, sourcegnodes, sourceelements, conditions_to_copy);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::UTILS::XFEMDiscretizationBuilder::SplitDiscretization(
    Teuchos::RCP<DRT::Discretization> sourcedis, Teuchos::RCP<DRT::Discretization> targetdis,
    const std::map<int, DRT::Node*>& sourcenodes, const std::map<int, DRT::Node*>& sourcegnodes,
    const std::map<int, Teuchos::RCP<DRT::Element>>& sourceelements,
    const std::vector<std::string>& conditions_to_copy) const
{
  if (!sourcedis->Filled()) dserror("sourcedis is not filled");
  const int myrank = targetdis->Comm().MyPID();

  const int numothernoderow = sourcedis->NumMyRowNodes();
  const int numothernodecol = sourcedis->NumMyColNodes();

  // add the conditioned elements
  for (std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator sourceele_iter =
           sourceelements.begin();
       sourceele_iter != sourceelements.end(); ++sourceele_iter)
  {
    if (sourceele_iter->second->Owner() == myrank)
    {
      targetdis->AddElement(Teuchos::rcp(sourceele_iter->second->Clone(), false));
    }
  }

  // row/col sets of conditioned node ids
  std::set<int> condnoderowset;
  std::set<int> condnodecolset;
  // row/col vectors of target node ids
  std::vector<int> targetnoderowvec;
  targetnoderowvec.reserve(sourcenodes.size());
  std::vector<int> targetnodecolvec;
  targetnodecolvec.reserve(sourcegnodes.size());

  // ------------------------------------------------------------------------
  // add conditioned nodes and fill the id vectors
  // ------------------------------------------------------------------------
  for (std::map<int, DRT::Node*>::const_iterator sourcegnode_iter = sourcegnodes.begin();
       sourcegnode_iter != sourcegnodes.end(); ++sourcegnode_iter)
  {
    const int nid = sourcegnode_iter->first;
    if (sourcegnode_iter->second->Owner() == myrank)
    {
      Teuchos::RCP<DRT::Node> sourcegnode =
          Teuchos::rcp(new DRT::Node(nid, sourcegnode_iter->second->X(), myrank));
      targetdis->AddNode(sourcegnode);
      condnoderowset.insert(nid);
      targetnoderowvec.push_back(nid);
    }
    condnodecolset.insert(nid);
    targetnodecolvec.push_back(nid);
  }

  // ------------------------------------------------------------------------
  // copy selected conditions to the new discretization
  // ------------------------------------------------------------------------
  for (std::vector<std::string>::const_iterator conditername = conditions_to_copy.begin();
       conditername != conditions_to_copy.end(); ++conditername)
  {
    std::vector<DRT::Condition*> conds;
    sourcedis->GetCondition(*conditername, conds);
    for (unsigned i = 0; i < conds.size(); ++i)
    {
      Teuchos::RCP<DRT::Condition> cond_to_copy =
          SplitCondition(conds[i], targetnodecolvec, targetdis->Comm());
      if (not cond_to_copy.is_null()) targetdis->SetCondition(*conditername, cond_to_copy);
    }
  }

  Redistribute(targetdis, targetnoderowvec, targetnodecolvec);

  // ------------------------------------------------------------------------
  // remove all nodes from the condnodecol and condnoderow sets, which also
  // belong to a not deleted source element
  // ------------------------------------------------------------------------
  for (unsigned j = 0; j < static_cast<unsigned>(sourcedis->NumMyColElements()); ++j)
  {
    int source_ele_gid = sourcedis->ElementColMap()->GID(j);
    // continue, if we are going to delete this element
    if (sourceelements.find(source_ele_gid) != sourceelements.end()) continue;
    DRT::Element* source_ele = sourcedis->gElement(source_ele_gid);
    const int* nid = source_ele->NodeIds();
    for (unsigned i = 0; i < static_cast<unsigned>(source_ele->NumNode()); ++i)
    {
      // Remove all nodes from the condition sets, which should stay in
      // the source discretization, since they belong to elements
      // which are not going to be deleted!
      std::set<int>::iterator pos = condnodecolset.find(nid[i]);
      if (pos != condnodecolset.end()) condnodecolset.erase(pos);
      pos = condnoderowset.find(nid[i]);
      if (pos != condnoderowset.end()) condnoderowset.erase(pos);
    }
  }

  // row/col vectors of non-conditioned node ids
  std::vector<int> othernoderowvec;
  othernoderowvec.reserve(numothernoderow - condnoderowset.size());
  std::vector<int> othernodecolvec;
  othernodecolvec.reserve(numothernodecol - condnodecolset.size());

  // determine non-conditioned nodes
  for (int lid = 0; lid < sourcedis->NodeColMap()->NumMyElements(); ++lid)
  {
    const int nid = sourcedis->NodeColMap()->GID(lid);

    // if we erase this node, we do not add it and just go on
    if (condnodecolset.find(nid) != condnodecolset.end()) continue;

    othernodecolvec.push_back(nid);

    if (sourcedis->NodeRowMap()->LID(nid) > -1) othernoderowvec.push_back(nid);
  }
  // delete conditioned nodes, which are not connected to any unconditioned elements
  for (std::set<int>::iterator it = condnodecolset.begin(); it != condnodecolset.end(); ++it)
    if (not sourcedis->DeleteNode(*it)) dserror("Node %d could not be deleted!", *it);

  // delete conditioned elements from source discretization
  for (std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator sourceele_iter =
           sourceelements.begin();
       sourceele_iter != sourceelements.end(); ++sourceele_iter)
  {
    sourcedis->DeleteElement(sourceele_iter->first);
  }

  // ------------------------------------------------------------------------
  // validate the source conditions
  // ------------------------------------------------------------------------
  std::vector<std::string> src_conditions;
  sourcedis->GetConditionNames(src_conditions);
  for (std::vector<std::string>::const_iterator conditername = src_conditions.begin();
       conditername != src_conditions.end(); ++conditername)
  {
    std::vector<DRT::Condition*> conds;
    sourcedis->GetCondition(*conditername, conds);
    std::vector<Teuchos::RCP<DRT::Condition>> src_conds(conds.size(), Teuchos::null);
    for (unsigned i = 0; i < conds.size(); ++i)
      src_conds[i] = SplitCondition(conds[i], othernodecolvec, sourcedis->Comm());
    sourcedis->ReplaceConditions(*conditername, src_conds);
  }
  // re-partioning
  Redistribute(sourcedis, othernoderowvec, othernodecolvec);


  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::UTILS::XFEMDiscretizationBuilder::Redistribute(Teuchos::RCP<DRT::Discretization> dis,
    std::vector<int>& noderowvec, std::vector<int>& nodecolvec) const
{
  dis->CheckFilledGlobally();

  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(dis->Comm().Clone());

  Teuchos::RCP<Epetra_Map> noderowmap =
      Teuchos::rcp(new Epetra_Map(-1, noderowvec.size(), &noderowvec[0], 0, *comm));

  Teuchos::RCP<Epetra_Map> nodecolmap =
      Teuchos::rcp(new Epetra_Map(-1, nodecolvec.size(), &nodecolvec[0], 0, *comm));
  if (!dis->Filled()) dis->Redistribute(*noderowmap, *nodecolmap);

  Teuchos::RCP<Epetra_Map> elerowmap = Teuchos::rcp(new Epetra_Map(*dis->ElementRowMap()));
  DRT::UTILS::REBALANCING::ComputeRebalancedNodeMaps(
      dis, elerowmap, noderowmap, nodecolmap, comm, false, comm->NumProc());

  Teuchos::RCP<Epetra_Map> roweles = Teuchos::null;
  Teuchos::RCP<Epetra_Map> coleles = Teuchos::null;
  dis->BuildElementRowColumn(*noderowmap, *nodecolmap, roweles, coleles);

  dis->ExportRowNodes(*noderowmap);
  dis->ExportRowElements(*roweles);

  dis->ExportColumnNodes(*nodecolmap);
  dis->ExportColumnElements(*coleles);

  dis->FillComplete();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::UTILS::XFEMDiscretizationBuilder::SplitDiscretizationByBoundaryCondition(
    const Teuchos::RCP<DRT::Discretization>& sourcedis,
    const Teuchos::RCP<DRT::Discretization>& targetdis,
    const std::vector<DRT::Condition*>& boundary_conds,
    const std::vector<std::string>& conditions_to_copy) const
{
  if (not sourcedis->Filled()) dserror("sourcedis is not filled");
  const int myrank = targetdis->Comm().MyPID();

  // element map
  std::map<int, Teuchos::RCP<DRT::Element>> src_cond_elements;

  // find conditioned nodes (owned and ghosted) and elements
  DRT::UTILS::FindConditionObjects(src_cond_elements, boundary_conds);

  std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator cit;
  std::map<int, Teuchos::RCP<DRT::Element>> src_elements;
  // row node map (id -> pointer)
  std::map<int, DRT::Node*> src_my_gnodes;
  std::vector<int> condnoderowvec;
  // column node map
  std::map<int, DRT::Node*> src_gnodes;
  std::vector<int> condnodecolvec;
  // find all parent elements
  for (cit = src_cond_elements.begin(); cit != src_cond_elements.end(); ++cit)
  {
    DRT::FaceElement* src_face_element = dynamic_cast<DRT::FaceElement*>(cit->second.get());
    if (src_face_element == NULL)
      dserror("Dynamic cast failed! The src element %d is no DRT::FaceElement!", cit->second->Id());
    // get the parent element
    DRT::Element* src_ele = src_face_element->ParentElement();
    int src_ele_gid = src_face_element->ParentElementId();
    src_elements[src_ele_gid] = Teuchos::rcp<DRT::Element>(src_ele, false);
    const int* n = src_ele->NodeIds();
    for (unsigned i = 0; i < static_cast<unsigned>(src_ele->NumNode()); ++i)
    {
      const int gid = n[i];
      if (sourcedis->HaveGlobalNode(gid))
      {
        DRT::Node* node = sourcedis->gNode(gid);
        src_gnodes[gid] = node;

        if (node->Owner() == myrank) src_my_gnodes[gid] = node;
      }
      else
        dserror("All nodes of known elements must be known!");
    }
  }

  SplitDiscretization(
      sourcedis, targetdis, src_my_gnodes, src_gnodes, src_elements, conditions_to_copy);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Condition> XFEM::UTILS::XFEMDiscretizationBuilder::SplitCondition(
    const DRT::Condition* src_cond, const std::vector<int>& nodecolvec,
    const Epetra_Comm& comm) const
{
  const std::vector<int>* cond_node_gids = src_cond->Nodes();
  std::set<int> nodecolset;
  nodecolset.insert(nodecolvec.begin(), nodecolvec.end());

  int lcount = 0;
  int gcount = 0;
  for (unsigned i = 0; i < cond_node_gids->size(); ++i)
  {
    int ngid = cond_node_gids->at(i);
    // add the node GID, if it is also a part of the new discretization
    if (nodecolset.find(ngid) != nodecolset.end()) lcount++;
  }

  comm.SumAll(&lcount, &gcount, 1);
  // return a Teuchos::null pointer, if there is nothing to copy
  if (gcount == 0) return Teuchos::null;

  // copy and keep this src condition
  return Teuchos::rcp(new DRT::Condition(*src_cond));
}

///*----------------------------------------------------------------------------*
// *----------------------------------------------------------------------------*/
// Teuchos::RCP<DRT::Condition> XFEM::UTILS::XFEMDiscretizationBuilder::
//    SplitCondition(
//    const DRT::Condition& src_cond,
//    const std::vector<int>& nodecolvec) const
//{
//  const std::vector<int>* cond_node_gids = src_cond.Nodes();
//  std::set<int> nodecolset;
//  nodecolset.insert(nodecolvec.begin(),nodecolvec.end());
//
//  Teuchos::RCP<std::vector<int> > keep_node_col_gids =
//      Teuchos::rcp(new std::vector<int>(0));
//  keep_node_col_gids->reserve(cond_node_gids->size());
//
//  for (unsigned i=0;i<cond_node_gids->size();++i)
//  {
//    int ngid = cond_node_gids->at(i);
//    // add the node GID, if it is also a part of the new discretization
//    if (nodecolset.find(ngid)!=nodecolset.end())
//    {
//      keep_node_col_gids->push_back(ngid);
//    }
//  }
//  // return a Teuchos::null pointer, if there is nothing to copy
//  if (keep_node_col_gids->size()==0)
//    return Teuchos::null;
//
//  // create a new target condition from source condition
//  Teuchos::RCP<DRT::Condition> target_cond =
//      Teuchos::rcp(new DRT::Condition(src_cond));
//  // overwrite node ids
//  target_cond->Add("Node Ids",keep_node_col_gids);
//
//  return target_cond;
//}

/*----------------------------------------------------------------------*/
/*!
\file drt_condition_utils.cpp

\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_condition_utils.H"
#include "standardtypes_cpp.H"
//#include "adapter_coupling_mortar.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"

#include <map>
#include <set>
#include <string>
#include <vector>
#include <algorithm>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::FindInterfaceObjects(const DRT::Discretization& dis,
                                      map<int, DRT::Node*>& nodes,
                                      map<int, RCP<DRT::Element> >& elements,
                                      const string& condname)
{
  int myrank = dis.Comm().MyPID();
  vector<DRT::Condition*> conds;
  dis.GetCondition(condname, conds);
  for (unsigned i = 0; i < conds.size(); ++i)
  {
    // get this condition's nodes
    const vector<int>* n = conds[i]->Nodes();
    for (unsigned j=0; j < n->size(); ++j)
    {
      const int gid = (*n)[j];
      if (dis.HaveGlobalNode(gid) and dis.gNode(gid)->Owner() == myrank)
      {
        nodes[gid] = dis.gNode(gid);
      }
    }

    // get this condition's elements
    map< int, RCP< DRT::Element > >& geo = conds[i]->Geometry();
    map< int, RCP< DRT::Element > >::iterator iter, pos;
    pos = elements.begin();
    for (iter = geo.begin(); iter != geo.end(); ++iter)
    {
      if (iter->second->Owner() == myrank)
      {
        pos = elements.insert(pos, *iter);
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::FindInterfaceObjects(const DRT::Discretization& dis,
                                      map<int, RCP<DRT::Element> >& elements,
                                      const string& condname)
{
  vector<DRT::Condition*> conds;
  dis.GetCondition(condname, conds);
  for (unsigned i = 0; i < conds.size(); ++i)
  {
    // get this condition's elements
    map< int, RCP< DRT::Element > >& geo = conds[i]->Geometry();
    map< int, RCP< DRT::Element > >::iterator iter, pos;
    pos = elements.begin();
    for (iter = geo.begin(); iter != geo.end(); ++iter)
    {
      // get all elements locally known, including ghost elements
      pos = elements.insert(pos, *iter);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> DRT::UTILS::CreateDiscretizationFromCondition(
    Teuchos::RCP<DRT::Discretization>  sourcedis,
        const string&                  condname,
        const string&                  discret_name,
        const string&                  element_name,
        const vector<string>&          conditions_to_copy
        )
{
  RCP<Epetra_Comm> com = rcp(sourcedis->Comm().Clone());
  RCP<DRT::Discretization> conditiondis = rcp(new DRT::Discretization(discret_name,com));

  // make sure connectivity is all set
  // we don't care, whether dofs exist or not
  if (!sourcedis->Filled())
    dserror("sourcedis is not filled");

  const int myrank = conditiondis->Comm().MyPID();
  const Epetra_Map* sourcenoderowmap = sourcedis->NodeRowMap();

  // We need to test for all elements (including ghosted ones) to
  // catch all nodes
  map<int, RCP<DRT::Element> >  sourceelements;
  DRT::UTILS::FindInterfaceObjects(*sourcedis, sourceelements, condname);

  set<int> rownodeset;
  set<int> colnodeset;

  // construct new elements
  for (map<int, RCP<DRT::Element> >::const_iterator sourceele_iter = sourceelements.begin();
       sourceele_iter != sourceelements.end();
       ++sourceele_iter)
  {
    const RCP<DRT::Element> sourceele = sourceele_iter->second;

    // get global node ids
    vector<int> nids;
    nids.reserve(sourceele->NumNode());
    transform(sourceele->Nodes(), sourceele->Nodes()+sourceele->NumNode(),
              back_inserter(nids), mem_fun(&DRT::Node::Id));

    if (std::count_if(nids.begin(), nids.end(), DRT::UTILS::MyGID(sourcenoderowmap))==0)
    {
      dserror("no own node in element %d", sourceele->Id());
    }

    if (std::count_if(nids.begin(), nids.end(),
                      DRT::UTILS::MyGID(sourcedis->NodeColMap())) < static_cast<int>(nids.size()))
    {
      dserror("element %d has remote non-ghost nodes",sourceele->Id());
    }

    copy(nids.begin(), nids.end(),
         inserter(colnodeset, colnodeset.begin()));

    // copy node ids of condition ele to rownodeset but leave those that do
    // not belong to this processor
    remove_copy_if(nids.begin(), nids.end(),
                   inserter(rownodeset, rownodeset.begin()),
                   not1(DRT::UTILS::MyGID(sourcenoderowmap)));

    // Do not clone ghost elements here! Those will be handled by the
    // discretization itself.
    if (sourceele->Owner()==myrank)
    {
      // create an element with the same global element id
      RCP<DRT::Element> condele = DRT::UTILS::Factory(element_name, "Polynomial", sourceele->Id(), myrank);

      // set the same global node ids to the new element
      condele->SetNodeIds(nids.size(), &nids[0]);

      // add element
      conditiondis->AddElement(condele);
    }
  }

  // construct new nodes, which use the same global id as the source nodes
  for (int i=0; i<sourcenoderowmap->NumMyElements(); ++i)
  {
    const int gid = sourcenoderowmap->GID(i);
    if (rownodeset.find(gid)!=rownodeset.end())
    {
      const DRT::Node* sourcenode = sourcedis->lRowNode(i);
      conditiondis->AddNode(rcp(new DRT::Node(gid, sourcenode->X(), myrank)));
    }
  }

  // we get the node maps almost for free
  vector<int> condnoderowvec(rownodeset.begin(), rownodeset.end());
  rownodeset.clear();
  RCP<Epetra_Map> condnoderowmap = rcp(new Epetra_Map(-1,
                                                      condnoderowvec.size(),
                                                      &condnoderowvec[0],
                                                      0,
                                                      conditiondis->Comm()));
  condnoderowvec.clear();

  vector<int> condnodecolvec(colnodeset.begin(), colnodeset.end());
  colnodeset.clear();
  RCP<Epetra_Map> condnodecolmap = rcp(new Epetra_Map(-1,
                                                      condnodecolvec.size(),
                                                      &condnodecolvec[0],
                                                      0,
                                                      conditiondis->Comm()));
  condnodecolvec.clear();

  // copy selected conditions to the new discretization
  for (vector<string>::const_iterator conditername = conditions_to_copy.begin();
       conditername != conditions_to_copy.end();
       ++conditername)
  {
    vector<DRT::Condition*> conds;
    sourcedis->GetCondition(*conditername, conds);
    for (unsigned i=0; i<conds.size(); ++i)
    {
      // We use the same nodal ids and therefore we can just copy the conditions.
      conditiondis->SetCondition(*conditername, rcp(new DRT::Condition(*conds[i])));
    }
  }

  // redistribute nodes to column (ghost) map
  RedistributeWithNewNodalDistribution(*conditiondis, *condnoderowmap, *condnodecolmap);
  conditiondis->FillComplete();

  return conditiondis;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::RedistributeWithNewNodalDistribution(
    DRT::Discretization&     dis,
    const Epetra_Map&        noderowmap,
    const Epetra_Map&        nodecolmap
    )
{
  // redistribute nodes to column (ghost) map
  dis.ExportColumnNodes(nodecolmap);

  Teuchos::RCP< Epetra_Map > elerowmap;
  Teuchos::RCP< Epetra_Map > elecolmap;

  // now we have all elements in a linear map roweles
  // build resonable maps for elements from the
  // already valid and final node maps
  dis.BuildElementRowColumn(noderowmap, nodecolmap, elerowmap, elecolmap);

  // we can now export elements to resonable row element distribution
  dis.ExportRowElements(*elerowmap);

  // export to the column map / create ghosting of elements
  dis.ExportColumnElements(*elecolmap);

}


/*----------------------------------------------------------------------
 * collects elements by labels (have to be implemented in the           *
 * corresponding condition)                                             *
 *----------------------------------------------------------------------*/
void DRT::UTILS::CollectElementsByConditionLabel(
    const DRT::Discretization&           discret,
    std::map<int,std::set<int> >&        elementsByLabel,
    const string&                        name)
{
  // Reset
  elementsByLabel.clear();
  // get condition
  vector< DRT::Condition * >  conditions;
  discret.GetCondition (name, conditions);
  
  // collect elements by xfem coupling label
  for(vector<DRT::Condition*>::const_iterator conditer = conditions.begin(); conditer!= conditions.end(); ++conditer)
  {
    DRT::Condition* condition = *conditer;
    const int label = condition->GetInt("label");
    const vector<int> geometryMap = *condition->Nodes();
    vector<int>::const_iterator iterNode;
    for(iterNode = geometryMap.begin(); iterNode != geometryMap.end(); ++iterNode )
    {
      const int nodegid = *iterNode;
      const DRT::Node* node = discret.gNode(nodegid);
      const DRT::Element*const* elements = node->Elements();
      for (int iele=0;iele < node->NumElement(); ++iele)
      {
        const DRT::Element* element = elements[iele];
        elementsByLabel[label].insert(element->Id());
      }
    }
  }
  int numOfCollectedIds = 0;
  for(unsigned int i = 0; i < elementsByLabel.size(); i++)
    numOfCollectedIds += elementsByLabel[i].size();
  
  if(discret.NumMyColElements() != numOfCollectedIds)
    dserror("not all elements collected.");
}

#endif

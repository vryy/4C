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
//#include "adapter_coupling_mortar.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"

// we need to know all element types for the ale mesh creation
#include "../drt_f2/fluid2.H"
#include "../drt_f3/fluid3.H"

#include "../drt_ale2/ale2.H"
#include "../drt_ale3/ale3.H"

#include <map>
#include <set>
#include <string>
#include <vector>
#include <algorithm>

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;



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
// create ale discretization parallel to the fluid one
/*----------------------------------------------------------------------*/
void DRT::UTILS::CreateAleDiscretization()
{
  RefCountPtr<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(genprob.numff,0);
  RefCountPtr<DRT::Discretization> aledis   = DRT::Problem::Instance()->Dis(genprob.numaf,0);

  if (!fluiddis->Filled()) fluiddis->FillComplete();

  if (aledis->NumGlobalElements() or aledis->NumGlobalNodes())
  {
    dserror("There are %d elements and %d nodes in empty ale discretization. Panic.",
            aledis->NumGlobalElements(), aledis->NumGlobalNodes());
  }

  int myrank = aledis->Comm().MyPID();

  vector<int> egid;
  egid.reserve(fluiddis->NumMyRowElements());

  vector<string> aletype;
  aletype.reserve(fluiddis->NumMyRowElements());

  set<int> rownodeset;
  set<int> colnodeset;
  const Epetra_Map* noderowmap = fluiddis->NodeRowMap();

  // Loop all fluid elements and find the ones that live on an ale
  // mesh. Here we do ugly castings. Please do not take this for an
  // example of how to code in baci!

  // We need to test for all elements (including ghosted ones) to
  // catch all nodes attached to ale elements
  int numelements = fluiddis->NumMyColElements();

  for (int i=0; i<numelements; ++i)
  {
    DRT::Element* actele = fluiddis->lColElement(i);
    bool isale = false;
    bool found = false;
    bool myele = fluiddis->ElementRowMap()->MyGID(actele->Id());

#ifdef D_FLUID2
    DRT::ELEMENTS::Fluid2* f2 = dynamic_cast<DRT::ELEMENTS::Fluid2*>(actele);
    if (not found and f2!=NULL)
    {
      found = true;
      isale = f2->IsAle();
      if (isale and myele)
        aletype.push_back("ALE2");
    }
#endif

#ifdef D_FLUID3
    DRT::ELEMENTS::Fluid3* f3 = dynamic_cast<DRT::ELEMENTS::Fluid3*>(actele);
    if (not found and f3!=NULL)
    {
      found = true;
      isale = f3->IsAle();
      if (isale and myele)
        aletype.push_back("ALE3");
    }
#endif

    if (not found)
      dserror("unsupported fluid element type '%s'", typeid(*actele).name());

    if (isale)
    {
      if (myele)
        egid.push_back(actele->Id());

      // copy node ids of actele to rownodeset but leave those that do
      // not belong to this processor
      remove_copy_if(actele->NodeIds(), actele->NodeIds()+actele->NumNode(),
                     inserter(rownodeset, rownodeset.begin()),
                     not1(DRT::UTILS::MyGID(noderowmap)));

      copy(actele->NodeIds(), actele->NodeIds()+actele->NumNode(),
           inserter(colnodeset, colnodeset.begin()));
    }
  }

  // construct ale nodes
  for (int i=0; i<noderowmap->NumMyElements(); ++i)
  {
    int gid = noderowmap->GID(i);
    if (rownodeset.find(gid)!=rownodeset.end())
    {
      DRT::Node* fluidnode = fluiddis->lRowNode(i);
      aledis->AddNode(rcp(new DRT::Node(gid, fluidnode->X(), myrank)));
    }
  }

  // we get the node maps almost for free

  vector<int> alenoderowvec(rownodeset.begin(), rownodeset.end());
  rownodeset.clear();
  RefCountPtr<Epetra_Map> alenoderowmap = rcp(new Epetra_Map(-1,
                                                             alenoderowvec.size(),
                                                             &alenoderowvec[0],
                                                             0,
                                                             aledis->Comm()));
  alenoderowvec.clear();

  vector<int> alenodecolvec(colnodeset.begin(), colnodeset.end());
  colnodeset.clear();
  RefCountPtr<Epetra_Map> alenodecolmap = rcp(new Epetra_Map(-1,
                                                             alenodecolvec.size(),
                                                             &alenodecolvec[0],
                                                             0,
                                                             aledis->Comm()));
  alenodecolvec.clear();

  // now do the elements

  // We must not add a new material type here because that might move
  // the internal material vector. And each element material might
  // have a pointer to that vector. Too bad.
  // So we search for a StVenantKirchhoff material and take the first
  // one we find.
  int nummat = DRT::Problem::Instance()->NumMaterials();
  int matnr = -1;
  for (int i=0; i<nummat; ++i)
  {
    if (DRT::Problem::Instance()->Material(i).mattyp==m_stvenant)
    {
      // For historical reasons material numbers are given in FORTRAN
      // style.
      matnr = i+1;
      break;
    }
  }
  if (matnr==-1)
    dserror("No StVenantKirchhoff material defined. Cannot generate ale mesh.");

  // construct ale elements
  // The order of the ale elements might be different from that of the
  // fluid elements. We don't care. There are not dofs to these
  // elements.
  for (unsigned i=0; i<egid.size(); ++i)
  {
    DRT::Element* fluidele = fluiddis->gElement(egid[i]);

    // create the ale element with the same global element id
    RefCountPtr<DRT::Element> aleele = DRT::UTILS::Factory(aletype[i],"Polynomial",egid[i], myrank);

    // get global node ids of fluid element
    vector<int> nids;
    nids.reserve(fluidele->NumNode());
    transform(fluidele->Nodes(), fluidele->Nodes()+fluidele->NumNode(),
              back_inserter(nids), mem_fun(&DRT::Node::Id));

    // set the same global node ids to the ale element
    aleele->SetNodeIds(nids.size(), &nids[0]);

    // We need to set material and gauss points to complete element setup.
    // This is again really ugly as we have to extract the actual
    // element type in order to access the material property
#ifdef D_ALE
    DRT::ELEMENTS::Ale2* ale2 = dynamic_cast<DRT::ELEMENTS::Ale2*>(aleele.get());
    if (ale2!=NULL)
    {
      ale2->SetMaterial(matnr);
    }
    else
#endif
    {
#ifdef D_ALE
      DRT::ELEMENTS::Ale3* ale3 = dynamic_cast<DRT::ELEMENTS::Ale3*>(aleele.get());
      if (ale3!=NULL)
      {
        ale3->SetMaterial(matnr);
      }
      else
#endif
      {
        dserror("unsupported ale element type '%s'", typeid(*aleele).name());
      }
    }

    // add ale element
    aledis->AddElement(aleele);
  }

  // conditions

  // copy the conditions to the ale discretization
  vector<DRT::Condition*> cond;
  fluiddis->GetCondition("FSICoupling", cond);
  for (unsigned i=0; i<cond.size(); ++i)
  {
    // We use the same nodal ids and therefore we can just copy the
    // conditions.
    aledis->SetCondition("FSICoupling", rcp(new DRT::Condition(*cond[i])));
  }

  cond.clear();
  fluiddis->GetCondition("FREESURFCoupling", cond);
  for (unsigned i=0; i<cond.size(); ++i)
  {
    // We use the same nodal ids and therefore we can just copy the
    // conditions.
    aledis->SetCondition("FREESURFCoupling", rcp(new DRT::Condition(*cond[i])));
  }

  cond.clear();
  fluiddis->GetCondition("ALEDirichlet", cond);
  for (unsigned i=0; i<cond.size(); ++i)
  {
    // We use the same nodal ids and therefore we can just copy the
    // conditions. But here we rename it. So we have nice dirichlet
    // conditions at the ale.
    aledis->SetCondition("Dirichlet", rcp(new DRT::Condition(*cond[i])));
  }

  // now care about the parallel distribution
  //

  // Right now all fluid elements must be ale enabled, otherwise we
  // get a very nasty parallel bug!

#if 0
  // At first make sure we have the same starting point on all
  // processors! This is cruical as the ALE field might be smaller
  // than the fluid field and there might be processors that do not
  // have ALE nodes and elements. These are not reset yet!

  aledis->Reset();
#endif

  // redistribute nodes to column (ghost) map
  RedistributeWithNewNodalDistribution(*aledis, *alenoderowmap, *alenodecolmap);
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

  if (!sourcedis->Filled())
    sourcedis->FillComplete();

  const int myrank = conditiondis->Comm().MyPID();
  const Epetra_Map* sourcenoderowmap = sourcedis->NodeRowMap();

  // We need to test for all elements (including ghosted ones) to
  // catch all nodes
  map<int, RCP<DRT::Element> >  sourceelements;
  DRT::UTILS::FindInterfaceObjects(*sourcedis, sourceelements, condname);

  set<int> rownodeset;
  set<int> colnodeset;

  // construct new elements
  for (map<int, RCP<DRT::Element> >::const_iterator cuttereleiter = sourceelements.begin();
       cuttereleiter != sourceelements.end();
       ++cuttereleiter)
  {
    const RCP<DRT::Element> sourceele = cuttereleiter->second;

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

    // copy node ids of cutterele to rownodeset but leave those that do
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

      // set the same global node ids to the ale element
      condele->SetNodeIds(nids.size(), &nids[0]);

      // add boundary element
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

  // Now we are done. :)
  dis.FillComplete();
}

#endif

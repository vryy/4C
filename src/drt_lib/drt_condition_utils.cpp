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




void DRT::UTILS::FindInterfaceObjects(
    const DRT::Discretization& dis,
    map<int, DRT::Node*>& nodes,
    map<int, RefCountPtr<DRT::Element> >& elements,
    const string&              condname
    )
{
    int myrank = dis.Comm().MyPID();
    vector<DRT::Condition*> conds;
    dis.GetCondition( condname, conds );
    for ( unsigned i = 0; i < conds.size(); ++i )
    {
        // get this condition's nodes
        const vector<int>* n = conds[i]->Nodes();
        for ( unsigned j = 0; j < n->size(); ++j )
        {
            const int gid = (*n)[j];
            if ( dis.HaveGlobalNode( gid ) and dis.gNode( gid )->Owner() == myrank )
            {
                nodes[gid] = dis.gNode( gid );
            }
        }

        // get this condition's elements
        map< int, RefCountPtr< DRT::Element > > geo = conds[i]->Geometry();
        map< int, RefCountPtr< DRT::Element > >::iterator iter, pos;
        pos = elements.begin();
        for ( iter = geo.begin(); iter != geo.end(); ++iter )
        {
            if ( iter->second->Owner() == myrank )
            {
                pos = elements.insert( pos, *iter );
            }
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
    RefCountPtr<DRT::Element> aleele = DRT::UTILS::Factory(aletype[i], egid[i], myrank);

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

  aledis->ExportColumnNodes(*alenodecolmap);

  RefCountPtr< Epetra_Map > elerowmap;
  RefCountPtr< Epetra_Map > elecolmap;

  // now we have all elements in a linear map roweles
  // build resonable maps for elements from the
  // already valid and final node maps
  // note that nothing is actually redistributed in here
  aledis->BuildElementRowColumn(*alenoderowmap, *alenodecolmap, elerowmap, elecolmap);

  // we can now export elements to resonable row element distribution
  aledis->ExportRowElements(*elerowmap);

  // export to the column map / create ghosting of elements
  aledis->ExportColumnElements(*elecolmap);

  // Now we are done. :)
  aledis->FillComplete();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
RCP<DRT::Discretization> DRT::UTILS::CreateDiscretizationFromCondition(
        RCP<DRT::Discretization>  cutterdis,
        const string&             condname,
        const string&             discret_name,
        const string&             element_name
        )
{
  RCP<Epetra_Comm> com = rcp(cutterdis->Comm().Clone());

  RCP<DRT::Discretization> boundarydis = rcp(new DRT::Discretization(discret_name,com));

  if (!cutterdis->Filled()) cutterdis->FillComplete();

  const int myrank = boundarydis->Comm().MyPID();

  if (myrank == 0)
  {
      cout << "creating discretization <"<< discret_name <<"> from condition <" << condname <<">" << endl;
  }
//  vector< DRT::Condition* >      xfemConditions;
//  cutterdis->GetCondition ("XFEMCoupling", xfemConditions);
//
//  if(xfemConditions.size()==0)
//      cout << "number of fsi xfem conditions = 0 --> empty boundary discretization will be created" << endl;

  // vector with boundary ele id's
  vector<int> egid;
  //egid.reserve(cutterdis->NumMyRowElements());

  set<int> rownodeset;
  set<int> colnodeset;
  const Epetra_Map* cutternoderowmap = cutterdis->NodeRowMap();

  // Loop all cutter elements and find the ones that live on an ale
  // mesh.
  // We need to test for all elements (including ghosted ones) to
  // catch all nodes attached to cutter elements
  map<int, DRT::Node*>          cutternodes;
  map<int, RCP<DRT::Element> >  cutterelements;
  DRT::UTILS::FindInterfaceObjects(*cutterdis, cutternodes, cutterelements, condname);

  // Loop all cutter elements

  // We need to test for all elements (including ghosted ones) to
  // catch all nodes attached to ale elements
  //const int numelements = cutterdis->NumMyColElements();

  map<int, RCP<DRT::Element> >::const_iterator cuttereleiter;
  for ( cuttereleiter = cutterelements.begin(); cuttereleiter != cutterelements.end(); ++cuttereleiter)
  {
    const RCP<DRT::Element> cutterele = cuttereleiter->second;
//    cout << "cutterele" << endl;
//    cout << (*cutterele) << endl;

    egid.push_back(cutterele->Id());

    // copy node ids of cutterele to rownodeset but leave those that do
    // not belong to this processor
    remove_copy_if(cutterele->NodeIds(), cutterele->NodeIds()+cutterele->NumNode(),
                   inserter(rownodeset, rownodeset.begin()),
                   not1(DRT::UTILS::MyGID(cutternoderowmap)));

    copy(cutterele->NodeIds(), cutterele->NodeIds()+cutterele->NumNode(),
        inserter(colnodeset, colnodeset.begin()));
  }

  // construct boundary nodes, which use the same global id as the cutter nodes
  for (int i=0; i<cutternoderowmap->NumMyElements(); ++i)
  {
    const int gid = cutternoderowmap->GID(i);
    if (rownodeset.find(gid)!=rownodeset.end())
    {
      const DRT::Node* cutternode = cutterdis->lRowNode(i);
      boundarydis->AddNode(rcp(new DRT::Node(gid, cutternode->X(), myrank)));
    }
  }

  // we get the node maps almost for free
  vector<int> boundarynoderowvec(rownodeset.begin(), rownodeset.end());
  rownodeset.clear();
  RCP<Epetra_Map> boundarynoderowmap = rcp(new Epetra_Map(-1,
                                                             boundarynoderowvec.size(),
                                                             &boundarynoderowvec[0],
                                                             0,
                                                             boundarydis->Comm()));
  boundarynoderowvec.clear();

  vector<int> boundarynodecolvec(colnodeset.begin(), colnodeset.end());
  colnodeset.clear();
  RCP<Epetra_Map> boundarynodecolmap = rcp(new Epetra_Map(-1,
                                                             boundarynodecolvec.size(),
                                                             &boundarynodecolvec[0],
                                                             0,
                                                             boundarydis->Comm()));
  boundarynodecolvec.clear();

  // now do the elements

  // construct boundary elements
  // The order of the boundary elements might be different from that of the
  // cutter elements. We don't care. There are not dofs to these
  // elements.
  for (unsigned i=0; i<egid.size(); ++i)
  {
    RCP<DRT::Element> cutterele = cutterelements[i];

    // create an element with the same global element id
    RCP<DRT::Element> boundaryele = DRT::UTILS::Factory(element_name, egid[i], myrank);

    // get global node ids of fluid element
    vector<int> nids;
    nids.reserve(cutterele->NumNode());
    transform(cutterele->Nodes(), cutterele->Nodes()+cutterele->NumNode(),
              back_inserter(nids), mem_fun(&DRT::Node::Id));

    // set the same global node ids to the ale element
    boundaryele->SetNodeIds(nids.size(), &nids[0]);

    // add boundary element
    boundarydis->AddElement(boundaryele);
//    cout << "boundary element:" << endl;
//    cout << (*boundaryele) << endl;
  }

  // conditions

  // copy the conditions to the boundary discretization
  // note, the condition is still named after the structure,
  // but that does not seem to matter in the subsequent computations
  vector<DRT::Condition*> conds;
  cutterdis->GetCondition(condname, conds);
  for (unsigned i=0; i<conds.size(); ++i)
  {
    // We use the same nodal ids and therefore we can just copy the
    // conditions.
    boundarydis->SetCondition(condname, rcp(new DRT::Condition(*conds[i])));
  }
  conds.clear();

  cutterdis->GetCondition("XFEMCoupling", conds);
  cout <<  "\n number of XFEMCoupling conditions (cutterdis)   " << conds.size() << endl;
  for (unsigned i=0; i<conds.size(); ++i)
  {
    // We use the same nodal ids and therefore we can just copy the
    // conditions.
    boundarydis->SetCondition("XFEMCoupling", rcp(new DRT::Condition(*conds[i])));
  }
  conds.clear();

  boundarydis->GetCondition("XFEMCoupling", conds);
  cout <<  "\n number of XFEMCoupling conditions (boundarydis) " << conds.size() << endl;
  
  // now care about the parallel distribution
  //

  // Right now all fluid elements must be ale enabled, otherwise we
  // get a very nasty parallel bug!

#if 0
  // At first make sure we have the same starting point on all
  // processors! This is cruical as the ALE field might be smaller
  // than the fluid field and there might be processors that do not
  // have ALE nodes and elements. These are not reset yet!

  boundarydis->Reset();
#endif

  // redistribute nodes to column (ghost) map

  boundarydis->ExportColumnNodes(*boundarynodecolmap);

  RefCountPtr< Epetra_Map > boundaryelerowmap;
  RefCountPtr< Epetra_Map > boundaryelecolmap;

  // now we have all elements in a linear map roweles
  // build resonable maps for elements from the
  // already valid and final node maps
  // note that nothing is actually redistributed in here
  boundarydis->BuildElementRowColumn(*boundarynoderowmap, *boundarynodecolmap, boundaryelerowmap, boundaryelecolmap);

  // we can now export elements to resonable row element distribution
  boundarydis->ExportRowElements(*boundaryelerowmap);

  // export to the column map / create ghosting of elements
  boundarydis->ExportColumnElements(*boundaryelecolmap);

  // Now we are done. :)
  boundarydis->FillComplete();
  //cout << (*boundarydis) << endl;

  return boundarydis;
}

#endif

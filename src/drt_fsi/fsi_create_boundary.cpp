/*----------------------------------------------------------------------*/
/*!
\file fsi_create_boundary.cpp

\brief Basis of all FSI algorithms

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <string>
#include <vector>
#include <set>
#include <functional>

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "fsi_dyn.H"
#include "fsi_dirichletneumann.H"
#include "fsi_monolithicoverlap.H"
#include "fsi_structureale.H"
#include "fsi_fluid_ale.H"
#include "fsi_utils.H"
#include "fsi_coupling_mortar.H"

#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_fluid/fluidresulttest.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "../drt_lib/drt_globalproblem.H"

// we need to know all element types for the boundary mesh creation
#include "../drt_xfem/bele3.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;



/*----------------------------------------------------------------------*/
// create boundary discretization from structure surface
/*----------------------------------------------------------------------*/
void CreateBoundaryDiscretization(RCP<DRT::Discretization> cutterdis)
{
  
  RCP<DRT::Discretization> boundarydis = DRT::Problem::Instance()->Dis(genprob.numbf,0);

  if (!cutterdis->Filled()) cutterdis->FillComplete(true,true,true,false);

  const int myrank = boundarydis->Comm().MyPID();
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
  FSI::FindInterfaceObjects(*cutterdis, cutternodes, cutterelements);
  
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
                   not1(FSI::UTILS::MyGID(cutternoderowmap)));

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

    // create the ale element with the same global element id
    RCP<DRT::Element> boundaryele = DRT::UTILS::Factory("BELE3", egid[i], myrank);

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
  cutterdis->GetCondition("FSICoupling", conds);
  for (unsigned i=0; i<conds.size(); ++i)
  {
    // We use the same nodal ids and therefore we can just copy the
    // conditions.
    boundarydis->SetCondition("FSICoupling", rcp(new DRT::Condition(*conds[i])));
  }
  conds.clear();

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
  boundarydis->FillComplete(true,true,true,false);
  cout << (*boundarydis) << endl;
}

#endif

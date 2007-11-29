
#ifdef CCADISCRET

#include <string>
#include <vector>
#include <set>
#include <functional>

#include "fsi_dyn.H"
#include "fsi_dirichletneumann.H"
#include "fsi_utils.H"

#include "../drt_mfsi/mfsi_algorithm.H"

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

// we need to know all element types for the ale mesh creation
#include "../drt_f2/fluid2.H"
#include "../drt_f3/fluid3.H"

#include "../drt_ale2/ale2.H"
#include "../drt_ale3/ale3.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;



/*----------------------------------------------------------------------*/
// create ale discretization parallel to the fluid one
/*----------------------------------------------------------------------*/
void CreateAleDiscretization()
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
    DRT::Elements::Fluid2* f2 = dynamic_cast<DRT::Elements::Fluid2*>(actele);
    if (not found and f2!=NULL)
    {
      found = true;
      isale = f2->IsAle();
      if (isale and myele)
        aletype.push_back("ALE2");
    }
#endif

#ifdef D_FLUID3
    DRT::Elements::Fluid3* f3 = dynamic_cast<DRT::Elements::Fluid3*>(actele);
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
                     not1(FSI::Utils::MyGID(noderowmap)));

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
    RefCountPtr<DRT::Element> aleele = DRT::Utils::Factory(aletype[i], egid[i], myrank);

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
    DRT::Elements::Ale2* ale2 = dynamic_cast<DRT::Elements::Ale2*>(aleele.get());
    if (ale2!=NULL)
    {
      ale2->SetMaterial(matnr);
    }
    else
#endif
    {
#ifdef D_ALE
      DRT::Elements::Ale3* ale3 = dynamic_cast<DRT::Elements::Ale3*>(aleele.get());
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
// entry point for FSI in DRT
/*----------------------------------------------------------------------*/
void fsi_ale_drt()
{
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  RefCountPtr<DRT::Discretization> aledis = DRT::Problem::Instance()->Dis(genprob.numaf,0);
  if (!aledis->Filled()) aledis->FillComplete();

  // create ale elements if the ale discretization is empty
  if (aledis->NumGlobalNodes()==0)
    CreateAleDiscretization();

  FSI_DYNAMIC *fsidyn = alldyn[3].fsidyn;

  if (fsidyn->ifsi != fsi_iter_monolithic)
  {
    Teuchos::RefCountPtr<FSI::DirichletNeumannCoupling> fsi = rcp(new FSI::DirichletNeumannCoupling(comm));

    if (genprob.restart)
    {
      // read the restart information, set vectors and variables
      fsi->ReadRestart(genprob.restart);
    }

    fsi->Timeloop(fsi);

#ifdef RESULTTEST
    DRT::ResultTestManager testmanager(comm);
    testmanager.AddFieldTest(rcp(new FluidResultTest(fsi->FluidField())));
    testmanager.TestAll();
#endif
  }
  else
  {
    Teuchos::RCP<MFSI::Algorithm> mfsi = rcp(new MFSI::Algorithm(comm));
    mfsi->Timeloop();
  }
}

#endif

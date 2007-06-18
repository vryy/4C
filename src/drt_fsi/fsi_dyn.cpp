
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include <string>
#include <vector>
#include <set>
#include <functional>

#include "fsi_dyn.H"
#include "fsi_dirichletneumann.H"
#include "fsi_utils.H"

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

  set<int> ngidset;
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

    DRT::Elements::Fluid2* f2 = dynamic_cast<DRT::Elements::Fluid2*>(actele);
    if (not found and f2!=NULL)
    {
      found = true;
      isale = f2->IsAle();
      if (isale and myele)
        aletype.push_back("ALE2");
    }

    DRT::Elements::Fluid3* f3 = dynamic_cast<DRT::Elements::Fluid3*>(actele);
    if (not found and f3!=NULL)
    {
      found = true;
      isale = f3->IsAle();
      if (isale and myele)
        aletype.push_back("ALE3");
    }

    if (not found)
      dserror("unsupported fluid element type '%s'", typeid(*actele).name());

    if (isale)
    {
      if (myele)
        egid.push_back(actele->Id());

      // copy node ids of actele to ngidset but leave those that do
      // not belong to this processor
      remove_copy_if(actele->NodeIds(), actele->NodeIds()+actele->NumNode(),
                     inserter(ngidset, ngidset.begin()),
                     not1(FSI::Utils::MyGID(noderowmap)));
    }
  }

  // construct ale nodes
  // Make sure the order of the ale nodes will be just the same as
  // those of the fluid nodes. This is important to preserve the
  // system matrix sparseness.
  for (int i=0; i<noderowmap->NumMyElements(); ++i)
  {
    int gid = noderowmap->GID(i);
    if (ngidset.find(gid)!=ngidset.end())
    {
      DRT::Node* fluidnode = fluiddis->lRowNode(i);
      aledis->AddNode(rcp(new DRT::Node(gid, fluidnode->X(), myrank)));
    }
  }

  // add a new material to use with the ale elements
  _MATERIAL localmat;
  localmat.mattyp = m_stvenant;
  localmat.m.stvenant = new _STVENANT();
  localmat.m.stvenant->youngs = 1.;
  localmat.m.stvenant->possionratio = 0.;
  localmat.m.stvenant->density = 0.;
  localmat.m.stvenant->thermexpans = 0.;
  DRT::Problem::Instance()->AddMaterial(localmat);

  // Material numbers are stored in fortran style, so we take the
  // number of materials after we added our new one.
  int matnr = DRT::Problem::Instance()->NumMaterials();

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
    aleele->SetNodeIds(fluidele->NumNode(), &nids[0]);

    // We need to set material and gauss points to complete element setup.
    // This is again really ugly as we have to extract the actual
    // element type in order to access the material property
    DRT::Elements::Ale2* ale2 = dynamic_cast<DRT::Elements::Ale2*>(aleele.get());
    if (ale2!=NULL)
    {
      ale2->SetMaterial(matnr);
      ale2->SetGaussPoints();
    }
    else
    {
      DRT::Elements::Ale3* ale3 = dynamic_cast<DRT::Elements::Ale3*>(aleele.get());
      if (ale3!=NULL)
      {
        ale3->SetMaterial(matnr);
        ale3->SetGaussPoints();
      }
      else
      {
        dserror("unsupported ale element type '%s'", typeid(*aleele).name());
      }
    }

    // add ale element
    aledis->AddElement(aleele);
  }

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

  // create ale elements if the ale discretization is empty
  if (DRT::Problem::Instance()->Dis(genprob.numaf,0)->NumGlobalNodes()==0)
    CreateAleDiscretization();

  Teuchos::RefCountPtr<FSI::DirichletNeumannCoupling> fsi = rcp(new FSI::DirichletNeumannCoupling(comm));

  fsi->Timeloop(fsi);

#ifdef RESULTTEST
  DRT::ResultTestManager testmanager(comm);
  testmanager.AddFieldTest(rcp(new FluidResultTest(fsi->FluidField())));
  testmanager.TestAll();
#endif
}

#endif
#endif

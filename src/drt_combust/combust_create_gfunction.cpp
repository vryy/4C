/*!----------------------------------------------------------------------
\file combust_create_gfunction.H

\brief create G-function discretization based on fluid discretization

	The G-function field is a convection-diffusion field. Based on 
	the fluid discretization, the G-function discretization is 
	created. BACI knows two discretizations.

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifdef PARALLEL
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "combust_create_gfunction.H"
#include <string>
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_lib/drt_condition_utils.H"

/*----------------------------------------------------------------------*
 | element types for fluid and G-function discretization    henke 06/08 |
 *----------------------------------------------------------------------*/
#include "../drt_f2/fluid2.H"
#include "../drt_f3/fluid3.H"
#include "../drt_condif2/condif2.H"
#include "../drt_condif3/condif3.H"

namespace COMBUST
{

/*----------------------------------------------------------------------*
 | create a G-function discretization from the fluid        henke 06/08 |
 *----------------------------------------------------------------------*/
void CreateGfuncDiscretization(int disnumff,int disnumcdf)
{
  if (disnumff == disnumcdf) dserror("Got identical discretization ids");

  RefCountPtr<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(disnumff,0);
  RefCountPtr<DRT::Discretization> condifdis   = DRT::Problem::Instance()->Dis(disnumcdf,0);

  if (!fluiddis->Filled()) fluiddis->FillComplete();

  // is the second discretization really empty?
  if (condifdis->NumGlobalElements() or condifdis->NumGlobalNodes())
  {
    dserror("There are %d elements and %d nodes in empty discretization. Panic.",
            condifdis->NumGlobalElements(), condifdis->NumGlobalNodes());
  }

  int myrank = condifdis->Comm().MyPID();

  vector<int> egid;
  egid.reserve(fluiddis->NumMyRowElements());

  vector<string> condiftype;
  condiftype.reserve(fluiddis->NumMyRowElements());

  set<int> rownodeset;
  set<int> colnodeset;
  const Epetra_Map* noderowmap = fluiddis->NodeRowMap();

  // Loop all fluid elements and find the ones that live on an condif
  // mesh. Here we do ugly castings. Please do not take this for an
  // example of how to code in baci!

  // We need to test for all elements (including ghosted ones) to
  // catch all nodes attached to condif elements
  int numelements = fluiddis->NumMyColElements();

  for (int i=0; i<numelements; ++i)
  {
    DRT::Element* actele = fluiddis->lColElement(i);
    bool found = false;
    bool myele = fluiddis->ElementRowMap()->MyGID(actele->Id());

#ifdef D_FLUID2
    DRT::ELEMENTS::Fluid2* f2 = dynamic_cast<DRT::ELEMENTS::Fluid2*>(actele);
    if (not found and f2!=NULL)
    {
      found = true;
      condiftype.push_back("CONDIF2");
    }
#endif

#ifdef D_FLUID3
    DRT::ELEMENTS::Fluid3* f3 = dynamic_cast<DRT::ELEMENTS::Fluid3*>(actele);
    if (not found and f3!=NULL)
    {
      found = true;
      condiftype.push_back("CONDIF3");
    }
#endif

    if (not found)
      dserror("unsupported fluid element type '%s'", typeid(*actele).name());

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

  // construct condif nodes
  for (int i=0; i<noderowmap->NumMyElements(); ++i)
  {
    int gid = noderowmap->GID(i);
    if (rownodeset.find(gid)!=rownodeset.end())
    {
      DRT::Node* fluidnode = fluiddis->lRowNode(i);
      condifdis->AddNode(rcp(new DRT::Node(gid, fluidnode->X(), myrank)));
    }
  }

  // we get the node maps almost for free

  vector<int> condifnoderowvec(rownodeset.begin(), rownodeset.end());
  rownodeset.clear();
  RefCountPtr<Epetra_Map> condifnoderowmap = rcp(new Epetra_Map(-1,
                                                             condifnoderowvec.size(),
                                                             &condifnoderowvec[0],
                                                             0,
                                                             condifdis->Comm()));
  condifnoderowvec.clear();

  vector<int> condifnodecolvec(colnodeset.begin(), colnodeset.end());
  colnodeset.clear();
  RefCountPtr<Epetra_Map> condifnodecolmap = rcp(new Epetra_Map(-1,
                                                             condifnodecolvec.size(),
                                                             &condifnodecolvec[0],
                                                             0,
                                                             condifdis->Comm()));
  condifnodecolvec.clear();

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
    if (DRT::Problem::Instance()->Material(i).mattyp==m_condif)
    {
      // For historical reasons material numbers are given in FORTRAN
      // style.
      matnr = i+1;
      break;
    }
  }
  if (matnr==-1)
    dserror("No ConDif material defined. Cannot generate convection-diffusion mesh.");

  // construct condif elements
  // The order of the condif elements might be different from that of the
  // fluid elements. We don't care. There are not dofs to these
  // elements.
  for (unsigned i=0; i<egid.size(); ++i)
  {
    DRT::Element* fluidele = fluiddis->gElement(egid[i]);

    // create the condif element with the same global element id
    RefCountPtr<DRT::Element> condifele = DRT::UTILS::Factory(condiftype[i],
                                                              "Polynomial" ,
                                                              egid[i]      ,
                                                              myrank       );

    // get global node ids of fluid element
    vector<int> nids;
    nids.reserve(fluidele->NumNode());
    transform(fluidele->Nodes(), fluidele->Nodes()+fluidele->NumNode(),
              back_inserter(nids), mem_fun(&DRT::Node::Id));

    // set the same global node ids to the condif element
    condifele->SetNodeIds(nids.size(), &nids[0]);

    // We need to set material and gauss points to complete element setup.
    // This is again really ugly as we have to extract the actual
    // element type in order to access the material property

#ifdef D_FLUID2
    DRT::ELEMENTS::Condif2* condif2 = dynamic_cast<DRT::ELEMENTS::Condif2*>(condifele.get());
    if (condif2!=NULL)
    {
      condif2->SetMaterial(matnr);
    }
    else
#endif
    {
#ifdef D_FLUID3
      DRT::ELEMENTS::Condif3* condif3 = dynamic_cast<DRT::ELEMENTS::Condif3*>(condifele.get());
      if (condif3!=NULL)
      {
        condif3->SetMaterial(matnr);
      }
      else
#endif
      {
        dserror("unsupported element type '%s'", typeid(*condifele).name());
      }
    }

    // add new condif element to discretization
    condifdis->AddElement(condifele);

  }

  // conditions

  // copy the conditions to the condif discretization
  vector<DRT::Condition*> cond;
  fluiddis->GetCondition("TransportDirichlet", cond);
  for (unsigned i=0; i<cond.size(); ++i)
  {
    // We use the same nodal ids and therefore we can just copy the
    // conditions. But here we rename it. So we have nice dirichlet
    // conditions at the condif discretization.
    condifdis->SetCondition("Dirichlet", rcp(new DRT::Condition(*cond[i])));
    cout<<"...transferred ConDif Dirichlet condition no. "<<i+1<<endl;
  }

#if 0
  // a small hack !!
  fluiddis->GetCondition("FluidStressCalc", cond);
  for (unsigned i=0; i<cond.size(); ++i)
  {
    // We use the same nodal ids and therefore we can just copy the
    // conditions. But here we rename it. So we have nice dirichlet
    // conditions at the condif.
    condifdis->SetCondition("InitialField", rcp(new DRT::Condition(*cond[i])));
  }

#endif

  // now care about the parallel distribution

#if 0
  // At first make sure we have the same starting point on all
  // processors! This is cruical as the ALE field might be smaller
  // than the fluid field and there might be processors that do not
  // have ALE nodes and elements. These are not reset yet!

  condifdis->Reset();
#endif

  // redistribute nodes to column (ghost) map
  DRT::UTILS::RedistributeWithNewNodalDistribution(*condifdis, *condifnoderowmap, *condifnodecolmap);
}
} // namespace COMBUST

#endif  // #ifdef CCADISCRET

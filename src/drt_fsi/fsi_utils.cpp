
#ifdef CCADISCRET

#include <map>
#include <set>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>

#include "fsi_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/standardtypes_cpp.H"

#include <Epetra_CrsMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <NOX.H>
#include <NOX_Epetra.H>
#include <Epetra_SerialDenseMatrix.h>

// we need to know all element types for the ale mesh creation
#include "../drt_f2/fluid2.H"
#include "../drt_f2/fluid2_nurbs.H"
#include "../drt_f3/fluid3.H"

#include "../drt_ale2/ale2.H"

#include "../drt_ale2/ale2_nurbs.H"
#include "../drt_ale3/ale3.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::DumpJacobian(NOX::Epetra::Interface::Required& interface,
                              double alpha,
                              double beta,
                              Teuchos::RefCountPtr<Epetra_Vector> soln,
                              std::string filename)
{
  // that's really stupid again
  const Epetra_BlockMap& bmap = soln->Map();
  Epetra_Map map(bmap.NumGlobalElements(),
                 bmap.NumMyElements(),
                 bmap.MyGlobalElements(),
                 0,
                 bmap.Comm());

  RefCountPtr<Epetra_CrsMatrix> jacobian = rcp(new Epetra_CrsMatrix(Copy, map, map.NumGlobalElements()));

  int nummyelements = map.NumMyElements();
  int mypos = LINALG::FindMyPos(nummyelements, map.Comm());
  double eta = 0.0;

  Epetra_Vector fo(*soln);
  Epetra_Vector fp(*soln);
  Epetra_Vector Jc(*soln);

  // Compute the RHS at the initial solution
  interface.computeF(*soln, fo, NOX::Epetra::Interface::Required::FD_Res);

  Epetra_Vector x_perturb = *soln;

  for (int i=0; i<map.NumGlobalElements(); ++i)
  {
    if (map.Comm().MyPID()==0)
      cout << "calculate column " << i << "\n";

    int proc = 0;
    int idx = 0;
    if (i>=mypos and i<mypos+nummyelements)
    {
      eta = alpha*(*soln)[i-mypos] + beta;
      x_perturb[i-mypos] += eta;
      idx = map.GID(i-mypos);
      proc = map.Comm().MyPID();
    }

    // Find what proc eta is on
    int broadcastProc = 0;
    map.Comm().SumAll(&proc, &broadcastProc, 1);

    // Send the perturbation variable, eta, to all processors
    map.Comm().Broadcast(&eta, 1, broadcastProc);

    map.Comm().Broadcast(&idx, 1, broadcastProc);

    // Compute the perturbed RHS
    interface.computeF(x_perturb, fp, NOX::Epetra::Interface::Required::FD_Res);

    // Compute the column k of the Jacobian
    Jc.Update(1.0, fp, -1.0, fo, 0.0);
    Jc.Scale( 1.0/eta );

    // Insert nonzero column entries into the jacobian
    for (int j = 0; j < map.NumMyElements(); ++j)
    {
      int gid = map.GID(j);
      if (Jc[j] != 0.0)
      {
        int err = jacobian->SumIntoGlobalValues(gid,1,&Jc[j],&idx);
        if (err>0)
        {
          err = jacobian->InsertGlobalValues(gid,1,&Jc[j],&idx);
        }
        if (err != 0)
          dserror("Assembly failed");
      }
    }

    // Unperturb the solution vector
    x_perturb = *soln;
  }

  jacobian->FillComplete();

  EpetraExt::RowMatrixToMatlabFile(filename.c_str(),*jacobian);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>
FSI::UTILS::ShiftMap(Teuchos::RCP<const Epetra_Map> emap,
                     const std::vector<Teuchos::RCP<const Epetra_Map> >& vecSpaces)
{
  int maxgid = 0;
  for (unsigned i=0; i<vecSpaces.size(); ++i)
  {
    maxgid = max(maxgid,vecSpaces[i]->MaxAllGID());
  }

  std::vector<int> gids;
  gids.reserve(emap->NumMyElements());
  std::transform(emap->MyGlobalElements(),
                 emap->MyGlobalElements()+emap->NumMyElements(),
                 std::back_inserter(gids),
                 std::bind2nd(std::plus<int>(),maxgid+1-emap->MinAllGID()));

  return Teuchos::rcp(new Epetra_Map(-1,gids.size(),&gids[0],0,emap->Comm()));
}


/*----------------------------------------------------------------------*/
// create ale discretization parallel to the fluid one
/*----------------------------------------------------------------------*/
void FSI::UTILS::CreateAleDiscretization()
{
  RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(genprob.numff,0);
  RCP<DRT::Discretization> aledis   = DRT::Problem::Instance()->Dis(genprob.numaf,0);

  // try to cast fluiddis to NurbsDiscretisation
  DRT::NURBS::NurbsDiscretization* nurbsdis
    =
    dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*(fluiddis)));

  if (!fluiddis->Filled()) fluiddis->FillComplete(false,true,true);

  if (aledis->NumGlobalElements() or aledis->NumGlobalNodes())
  {
    dserror("There are %d elements and %d nodes in empty ale discretization. Panic.",
            aledis->NumGlobalElements(), aledis->NumGlobalNodes());
  }

  const int myrank = aledis->Comm().MyPID();

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
  const int numelements = fluiddis->NumMyColElements();

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
  if(nurbsdis==NULL)
  {
    for (int i=0; i<noderowmap->NumMyElements(); ++i)
    {
      const int gid = noderowmap->GID(i);
      if (rownodeset.find(gid)!=rownodeset.end())
      {
        DRT::Node* fluidnode = fluiddis->lRowNode(i);
        aledis->AddNode(rcp(new DRT::Node(gid, fluidnode->X(), myrank)));
      }
    }
  }
  else
  {
    for (int i=0; i<noderowmap->NumMyElements(); ++i)
    {
      const int gid = noderowmap->GID(i);
      if (rownodeset.find(gid)!=rownodeset.end())
      {
        DRT::NURBS::ControlPoint* fluidnode 
        =
          dynamic_cast<DRT::NURBS::ControlPoint* >(fluiddis->lRowNode(i));
          aledis->AddNode(rcp(new DRT::NURBS::ControlPoint(gid, fluidnode->X(),fluidnode->W(),myrank)));
      }
    }
  }
  

  // we get the node maps almost for free

  vector<int> alenoderowvec(rownodeset.begin(), rownodeset.end());
  rownodeset.clear();
  RCP<Epetra_Map> alenoderowmap = rcp(new Epetra_Map(-1,
						     alenoderowvec.size(),
						     &alenoderowvec[0],
						     0,
						     aledis->Comm()));
  alenoderowvec.clear();

  vector<int> alenodecolvec(colnodeset.begin(), colnodeset.end());
  colnodeset.clear();
  RCP<Epetra_Map> alenodecolmap = rcp(new Epetra_Map(-1,
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
  const int matnr = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_stvenant);
  if (matnr==-1)
    dserror("No StVenantKirchhoff material defined. Cannot generate ale mesh.");

  // construct ale elements
  // The order of the ale elements might be different from that of the
  // fluid elements. We don't care. There are not dofs to these
  // elements.

  for (std::size_t i=0; i<egid.size(); ++i)
  {
    const DRT::Element* fluidele = fluiddis->gElement(egid[i]);


    string eletype = "Polynomial";
    if(nurbsdis!=NULL)
    {
      if(fluidele->NumNode()==9)
      {
        eletype="NURBS9";
      }
      else if(fluidele->NumNode()==4)
      {
        eletype="NURBS4";
      }
      else if(fluidele->NumNode()==27)
      {
        eletype="NURBS27";
      }
      else
      {
        dserror("unknown type of nurbs element\n");
      }
    }


    // create the ale element with the same global element id
    RCP<DRT::Element> aleele = DRT::UTILS::Factory(aletype[i],eletype,egid[i], myrank);

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
    if(nurbsdis==NULL)
    {
      DRT::ELEMENTS::Ale2* ale2 = dynamic_cast<DRT::ELEMENTS::Ale2*>(aleele.get());
      if (ale2!=NULL)
      {
        ale2->SetMaterial(matnr);
      }
      else
      {
        DRT::ELEMENTS::Ale3* ale3 = dynamic_cast<DRT::ELEMENTS::Ale3*>(aleele.get());
        if (ale3!=NULL)
        {
          ale3->SetMaterial(matnr);
        }
        else
        {
          dserror("unsupported ale element type '%s'", typeid(*aleele).name());
        }
      }
    }
    else
    {
      DRT::ELEMENTS::NURBS::Ale2Nurbs* ale2 = dynamic_cast<DRT::ELEMENTS::NURBS::Ale2Nurbs*>(aleele.get());
      if (ale2!=NULL)
      {
        ale2->SetMaterial(matnr);
      }
      else
      {
        dserror("unsupported ale element type '%s'", typeid(*aleele).name());
      }
    }
#endif


    // add ale element
    aledis->AddElement(aleele);
  }

  // conditions

  // copy the conditions to the ale discretization
  vector<DRT::Condition*> cond;
  fluiddis->GetCondition("FSICoupling", cond);
  for (std::size_t i=0; i<cond.size(); ++i)
  {
    // We use the same nodal ids and therefore we can just copy the
    // conditions.
    aledis->SetCondition("FSICoupling", rcp(new DRT::Condition(*cond[i])));
  }

  cond.clear();
  fluiddis->GetCondition("FREESURFCoupling", cond);
  for (std::size_t i=0; i<cond.size(); ++i)
  {
    // We use the same nodal ids and therefore we can just copy the
    // conditions.
    aledis->SetCondition("FREESURFCoupling", rcp(new DRT::Condition(*cond[i])));
  }

  cond.clear();
  fluiddis->GetCondition("ALEDirichlet", cond);
  for (std::size_t i=0; i<cond.size(); ++i)
  {
    // We use the same nodal ids and therefore we can just copy the
    // conditions. But here we rename it. So we have nice dirichlet
    // conditions at the ale.
    aledis->SetCondition("Dirichlet", rcp(new DRT::Condition(*cond[i])));
  }

  cond.clear();
  fluiddis->GetCondition("SurfacePeriodic", cond);
  for (std::size_t i=0; i<cond.size(); ++i)
  {
    // We use the same nodal ids and therefore we can just copy the
    // conditions.
    aledis->SetCondition("SurfacePeriodic", rcp(new DRT::Condition(*cond[i])));
  }

  cond.clear();
  fluiddis->GetCondition("LinePeriodic", cond);
  for (std::size_t i=0; i<cond.size(); ++i)
  {
    // We use the same nodal ids and therefore we can just copy the
    // conditions.
    aledis->SetCondition("LinePeriodic", rcp(new DRT::Condition(*cond[i])));
  }

  cond.clear();
  fluiddis->GetCondition("ElectrodeKinetics", cond);
  for (std::size_t i=0; i<cond.size(); ++i)
  {
    // We use the same nodal ids and therefore we can just copy the
    // conditions. 
    // Needed for Electrochemistry simulations with moving boundaries
    aledis->SetCondition("ElectrodeKinetics", rcp(new DRT::Condition(*cond[i])));
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
  DRT::UTILS::RedistributeWithNewNodalDistribution(*aledis, *alenoderowmap, *alenodecolmap);
  aledis->FillComplete();

  // finally do knot vectors in the nurbs case

  if(nurbsdis!=NULL)
  {
    DRT::NURBS::NurbsDiscretization* nurbsaledis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*(aledis)));

    if(nurbsaledis==NULL)
    {
      dserror("Nurbs fluid discretisation but no nurbs ALE discretisation\n");
    }
    
    Teuchos::RCP<DRT::NURBS::Knotvector> knots
      =
      Teuchos::rcp(new DRT::NURBS::Knotvector(*(nurbsdis->GetKnotVector())));

    // reset offsets
    int smallest_gid_in_dis=nurbsaledis->ElementRowMap()->MinAllGID();
    knots->FinishKnots(smallest_gid_in_dis);


    nurbsaledis->SetKnotVector(knots);
  }


}

#endif

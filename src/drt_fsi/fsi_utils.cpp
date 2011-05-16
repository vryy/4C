
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
#include "../linalg/linalg_utils.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_mat/matpar_bundle.H"

#include <Epetra_CrsMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <NOX.H>
#include <NOX_Epetra.H>
#include <Epetra_SerialDenseMatrix.h>

// we need to know all element types for the ale mesh creation
#include "../drt_f3/fluid3.H"
#include "../drt_f3/xfluid3.H"

#include "../drt_ale2/ale2.H"
#include "../drt_ale2/ale2_nurbs.H"
#include "../drt_ale3/ale3.H"
#include "../drt_ale3/ale3_nurbs.H"

#include "../drt_so3/so_surface.H"
#include "../drt_so3/so_line.H"

#include "fsi_debugwriter.H"
#include "../drt_geometry/searchtree.H"
#include "../drt_adapter/adapter_fluid_ale.H"
#include "../drt_adapter/adapter_coupling_mortar.H"
#include "../drt_mortar/mortar_interface.H"
#include "../drt_adapter/adapter_fluid.H"
#include "../drt_adapter/adapter_structure.H"
#include "../drt_ale/ale.H"

#include "../drt_io/io_control.H"

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
/*
void FSI::UTILS::CreateAleDiscretization()
{
  Teuchos::RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(genprob.numff,0);
  Teuchos::RCP<DRT::Discretization> aledis   = DRT::Problem::Instance()->Dis(genprob.numaf,0);

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
  Teuchos::RCP<Epetra_Map> alenoderowmap = rcp(new Epetra_Map(-1,
						     alenoderowvec.size(),
						     &alenoderowvec[0],
						     0,
						     aledis->Comm()));
  alenoderowvec.clear();

  vector<int> alenodecolvec(colnodeset.begin(), colnodeset.end());
  colnodeset.clear();
  Teuchos::RCP<Epetra_Map> alenodecolmap = rcp(new Epetra_Map(-1,
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
    Teuchos::RCP<DRT::Element> aleele = DRT::UTILS::Factory(aletype[i],eletype,egid[i], myrank);

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
        DRT::ELEMENTS::NURBS::Ale3Nurbs* ale3 = dynamic_cast<DRT::ELEMENTS::NURBS::Ale3Nurbs*>(aleele.get());
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
  fluiddis->GetCondition("StructAleCoupling", cond);
  for (unsigned i=0; i<cond.size(); ++i)
  {
    // We use the same nodal ids and therefore we can just copy the
    // conditions.
    aledis->SetCondition("StructAleCoupling", rcp(new DRT::Condition(*cond[i])));
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
*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<string,string> FSI::UTILS::AleFluidCloneStrategy::ConditionsToCopy()
{
  std::map<string,string> conditions_to_copy;

  conditions_to_copy.insert(pair<string,string>("ALEDirichlet","Dirichlet"));
  conditions_to_copy.insert(pair<string,string>("FSICoupling","FSICoupling"));
  conditions_to_copy.insert(pair<string,string>("FREESURFCoupling","FREESURFCoupling"));
  conditions_to_copy.insert(pair<string,string>("StructAleCoupling","StructAleCoupling"));
  conditions_to_copy.insert(pair<string,string>("LinePeriodic","LinePeriodic"));
  conditions_to_copy.insert(pair<string,string>("SurfacePeriodic","SurfacePeriodic"));
  conditions_to_copy.insert(pair<string,string>("ElectrodeKinetics","ElectrodeKinetics"));
  conditions_to_copy.insert(pair<string,string>("XFEMCoupling","XFEMCoupling"));

  return conditions_to_copy;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::AleFluidCloneStrategy::CheckMaterialType(const int matid)
{
  // no check implemented for ALE materials!
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::AleFluidCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele,
    DRT::Element* oldele,
    const int matid,
    const bool nurbsdis)
{
  // We must not add a new material type here because that might move
  // the internal material vector. And each element material might
  // have a pointer to that vector. Too bad.
  // So we search for a StVenantKirchhoff material and take the first
  // one we find.
  // => matid from outside remains unused!
  const int matnr = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_stvenant);
  if (matnr==-1)
    dserror("No StVenantKirchhoff material defined. Cannot generate ale mesh.");

#ifdef D_ALE
    if(nurbsdis==false)
    {
      DRT::ELEMENTS::Ale2* ale2 = dynamic_cast<DRT::ELEMENTS::Ale2*>(newele.get());
      if (ale2!=NULL)
      {
        ale2->SetMaterial(matnr);
      }
      else
      {
        DRT::ELEMENTS::Ale3* ale3 = dynamic_cast<DRT::ELEMENTS::Ale3*>(newele.get());
        if (ale3!=NULL)
        {
          ale3->SetMaterial(matnr);
        }
        else
        {
          dserror("unsupported ale element type '%s'", typeid(*newele).name());
        }
      }
    }
    else
    {
      DRT::ELEMENTS::NURBS::Ale2Nurbs* ale2 = dynamic_cast<DRT::ELEMENTS::NURBS::Ale2Nurbs*>(newele.get());
      if (ale2!=NULL)
      {
        ale2->SetMaterial(matnr);
      }
      else
      {
        DRT::ELEMENTS::NURBS::Ale3Nurbs* ale3 = dynamic_cast<DRT::ELEMENTS::NURBS::Ale3Nurbs*>(newele.get());

        if(ale3!=NULL)
        {
          ale3->SetMaterial(matnr);
        }
        else
        {
          dserror("unsupported ale element type '%s'", typeid(*newele).name());
        }
      }
    }
#endif

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::UTILS::AleFluidCloneStrategy::DetermineEleType(
    DRT::Element* actele,
    const bool ismyele,
    vector<string>& eletype)
{
  bool isale = false;
  bool found = false;

#ifdef D_FLUID3
    DRT::ELEMENTS::Fluid3* f3 = dynamic_cast<DRT::ELEMENTS::Fluid3*>(actele);
    if (not found and f3!=NULL)
    {
      const int  nsd = DRT::UTILS::getDimension(f3->Shape());
      found = true;
      isale = f3->IsAle();

      if (isale and ismyele)
      {
        if (nsd == 3)
          eletype.push_back("ALE3");
        else if (nsd == 2)
          eletype.push_back("ALE2");
        else
          dserror("%i D Dimension not supported", nsd);
      }
    }
    DRT::ELEMENTS::XFluid3* xf3 = dynamic_cast<DRT::ELEMENTS::XFluid3*>(actele);
    if (not found and xf3!=NULL)
    {
      const int  nsd = DRT::UTILS::getDimension(xf3->Shape());
      found = true;
      isale = xf3->IsAle();

      if (isale and ismyele)
      {
        if (nsd == 3)
          eletype.push_back("ALE3");
        else if (nsd == 2)
          eletype.push_back("ALE2");
        else
          dserror("%i D Dimension not supported", nsd);
      }
    }
#endif

    if (not found)
      dserror("unsupported fluid element type '%s'", typeid(*actele).name());

  return isale;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// class SlideAleUtils
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::UTILS::SlideAleUtils::SlideAleUtils
(
    Teuchos::RCP<DRT::Discretization> structdis,
    Teuchos::RCP<DRT::Discretization> fluiddis,
    ADAPTER::CouplingMortar& coupsf,
    bool structcoupmaster,
    INPAR::FSI::SlideALEProj aleproj
)
:
aletype_(aleproj),
slideeleredmap_(null)
{
  structcoupmaster_ =  structcoupmaster;

  // declare struct objects in interface
  map<int, RefCountPtr<DRT::Element> > structelements;
  map<int, RefCountPtr<DRT::Element> > structmelements;
  map<int, RefCountPtr<DRT::Element> > structdelements;
  map<int, DRT::Node*> structnodes; // dummy map
  map<int, DRT::Node*> structmnodes; // partial map of sticking structure nodes
  map<int, DRT::Node*> structdnodes; // partial map of centerdisp structure nodes
  map<int, DRT::Node*> structgnodes; // complete map of strucutre nodes

  //initialize struct objects in interface
  DRT::UTILS::FindConditionObjects(*structdis, structnodes, structgnodes, structelements,"FSICoupling");
  DRT::UTILS::FindConditionObjects(*structdis, structnodes, structmnodes, structmelements,"FSICouplingNoSlide");
  DRT::UTILS::FindConditionObjects(*structdis, structnodes, structdnodes, structdelements,"FSICouplingCenterDisp");
  istructdispnodes_ = structdnodes;
  istructdispeles_ = structdelements;
  istructslidnodes_ = structgnodes;
  istructslideles_ = structelements;

  vector<int> slideeleidvector;

  map<int, DRT::Node*>::iterator nit;
  for ( nit=structmnodes.begin() ; nit != structmnodes.end(); nit++ )
  {
    int err = istructslidnodes_.erase((*nit).first);
    if (!err)
      dserror("Non sliding interface has to be a subset of FSI-interface or empty");
  }

  map<int, RefCountPtr<DRT::Element> >::iterator eit;
  for ( eit=structmelements.begin() ; eit != structmelements.end(); eit++ )
  {
    int err = istructslideles_.erase((*eit).first);
    if (!err)
      dserror("Non sliding interface has to be a subset of FSI-interface or empty");
  }

  //build a redundant map of ids of all sliding elements
  for ( eit=istructslideles_.begin() ; eit != istructslideles_.end(); eit++ )
  {
    //build slideeleidvector with unique distribution. Otherwise, AllreduceEMap() will complain in DEBUG
    if (structdis->Comm().MyPID()==(*eit).second->Owner())
      slideeleidvector.push_back((*eit).first);
  }
  const Epetra_Map slideelemap (-1, slideeleidvector.size(), &slideeleidvector[0], 0, structdis->Comm());
  slideeleredmap_ = LINALG::AllreduceEMap(slideelemap);

  // declare struct objects in interface
  map<int, RefCountPtr<DRT::Element> > fluidelements;
  map<int, RefCountPtr<DRT::Element> > fluidmelements;
  // declare fluid objects in interface
  map<int, DRT::Node*> fluidnodes;  // complete map of fluid nodes
  map<int, DRT::Node*> fluidmnodes; // partial map of sticking fluid nodes
  map<int, DRT::Node*> fluidgnodes; // projection nodes

  //initialize struct objects in interface
  DRT::UTILS::FindConditionObjects(*fluiddis, fluidnodes, fluidgnodes, fluidelements,"FSICoupling");
  DRT::UTILS::FindConditionObjects(*fluiddis, fluidnodes, fluidmnodes, fluidmelements,"FSICouplingNoSlide");
  ifluidconfnodes_ = fluidmnodes;
  ifluidslidnodes_ = fluidnodes;
  ifluidslideles_ = fluidelements;

  for ( eit=fluidmelements.begin() ; eit != fluidmelements.end(); eit++ )
  {
    int err = ifluidslideles_.erase((*eit).first);
    if (!err)
      dserror("Non sliding interface has to be a subset of FSI-interface or empty");
  }

  for ( nit=ifluidconfnodes_.begin() ; nit != ifluidconfnodes_.end(); nit++ )
  {
    int err = ifluidslidnodes_.erase((*nit).first);
    if (!err)
      dserror("Non sliding interface has to be a subset of FSI-interface or empty");
  }

  Teuchos::RCP<Epetra_Map> structdofrowmap;
  Teuchos::RCP<Epetra_Map> fluiddofrowmap;


  // useful displacement vectors
  if (structcoupmaster_)
  {
    structdofrowmap_ = coupsf.MasterDofRowMap();
    fluiddofrowmap_ = coupsf.SlaveDofRowMap();
  }
  else
  {
    structdofrowmap_ = coupsf.SlaveDofRowMap();
    fluiddofrowmap_ = coupsf.MasterDofRowMap();
  }

  Teuchos::RCP<Epetra_Map> dofrowmap = LINALG::MergeMap(*structdofrowmap_,*fluiddofrowmap_, true);
  idispms_ = LINALG::CreateVector(*dofrowmap, true);

  iprojhist_=Teuchos::rcp(new Epetra_Vector(*fluiddofrowmap_,true));

  centerdisptotal_.resize(genprob.ndim);

  RedundantElements(coupsf,structdis->Comm());

  maxmindist_ = 1.0e-1;

  coupff_.Setup(*fluiddis);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::SlideAleUtils::Remeshing
(
    ADAPTER::Structure& structure,
    Teuchos::RCP<DRT::Discretization> fluiddis,
    Teuchos::RCP<Epetra_Vector> idispale,
    Teuchos::RCP<Epetra_Vector> iprojdispale,
    ADAPTER::CouplingMortar& coupsf,
    const Epetra_Comm& comm
)
{

#ifdef D_SOLID3
  Teuchos::RCP<Epetra_Vector> idisptotal = structure.ExtractInterfaceDispnp();
  const int dim = genprob.ndim;

  //project sliding fluid nodes onto struct interface surface
  SlideProjection(structure,
                  fluiddis,
                  idispale,
                  iprojdispale,
                  coupsf,
                  comm);

  //For the NON sliding ALE Nodes, use standard ALE displacements

  map<int, DRT::Node*>::const_iterator nodeiter;
  for (nodeiter = ifluidconfnodes_.begin(); nodeiter != ifluidconfnodes_.end(); ++nodeiter)
  {
    DRT::Node* node = nodeiter->second;
    vector<int> lids(dim);
    for(int p=0; p<dim; p++)
      //lids of gids of node
      lids[p] = fluiddofrowmap_->LID((fluiddis->Dof(node))[p]);

    // current coord of ale node = ref coord + ifluid_
    vector <double> finaldxyz(dim);

    for(int p=0; p<dim; p++)
      finaldxyz[p] =   (*idispale)[(lids[p])];

    int err = iprojdispale->ReplaceMyValues(dim, &finaldxyz[0], &lids[0]);
    if (err == 1) dserror("error while replacing values");

  }

  //merge displacement values of interface nodes (struct+fluid) into idispms_ for mortar
  idispms_->Scale(0.0);

  Teuchos::RCP<Epetra_Map> dofrowmap = LINALG::MergeMap(*structdofrowmap_,*fluiddofrowmap_, true);
  Teuchos::RCP<Epetra_Import> msimpo = rcp (new Epetra_Import(*dofrowmap,*structdofrowmap_));
  Teuchos::RCP<Epetra_Import> slimpo = rcp (new Epetra_Import(*dofrowmap,*fluiddofrowmap_));

  idispms_ -> Import(*idisptotal,*msimpo,Add);
  idispms_ -> Import(*iprojdispale,*slimpo,Add);

  iprojhist_->Update(1.0,*iprojdispale,0.0);
#endif
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::SlideAleUtils::EvaluateMortar
(
    Teuchos::RCP<Epetra_Vector> idisptotal,
    Teuchos::RCP<Epetra_Vector> ifluid,
    ADAPTER::CouplingMortar& coupsf
)
{

  //merge displacement values of interface nodes (struct+fluid) into idispms_ for mortar
  idispms_->Scale(0.0);

  Teuchos::RCP<Epetra_Map> dofrowmap = LINALG::MergeMap(*structdofrowmap_,*fluiddofrowmap_, true);
  Teuchos::RCP<Epetra_Import> msimpo = rcp (new Epetra_Import(*dofrowmap,*structdofrowmap_));
  Teuchos::RCP<Epetra_Import> slimpo = rcp (new Epetra_Import(*dofrowmap,*fluiddofrowmap_));

  idispms_ -> Import(*idisptotal,*msimpo,Add);
  idispms_ -> Import(*ifluid,*slimpo,Add);

  //new D,M,Dinv out of disp of struct and fluid side
  //coupsf.Evaluate(idispms_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::SlideAleUtils::EvaluateFluidMortar
(
    Teuchos::RCP<Epetra_Vector> ima,
    Teuchos::RCP<Epetra_Vector> isl
)
{
  //new D,M,Dinv out of fluid disp before and after sliding
  coupff_.Evaluate(ima,isl);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::UTILS::SlideAleUtils::InterpolateFluid
(
    Teuchos::RCP<const Epetra_Vector> uold
)
{
  Teuchos::RCP<Epetra_Vector> unew = coupff_.MasterToSlave(uold);
  unew->ReplaceMap(uold->Map());

  return unew;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
vector<double> FSI::UTILS::SlideAleUtils::Centerdisp
(
    ADAPTER::Structure& structure,
    const Epetra_Comm& comm
)
{
  Teuchos::RCP<DRT::Discretization> structdis = structure.Discretization();

  Teuchos::RCP<Epetra_Vector> idispn = structure.ExtractInterfaceDispn();
  Teuchos::RCP<Epetra_Vector> idisptotal = structure.ExtractInterfaceDispnp();
  Teuchos::RCP<Epetra_Vector> idispstep = structure.ExtractInterfaceDispnp();

  idispstep->Update(-1.0, *idispn, 1.0);

  const int dim = genprob.ndim;
  // get structure and fluid discretizations  and set stated for element evaluation

  const Teuchos::RCP<Epetra_Vector> idisptotalcol = LINALG::CreateVector(*structdis->DofColMap(),false);
  LINALG::Export(*idisptotal,*idisptotalcol);
  structdis->SetState("displacementtotal",idisptotalcol);
  const Teuchos::RCP<Epetra_Vector> idispstepcol = LINALG::CreateVector(*structdis->DofColMap(),false);
  LINALG::Export(*idispstep,*idispstepcol);
  structdis->SetState("displacementincr",idispstepcol);

  //define stuff needed by the elements
  Teuchos::ParameterList params;
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  //prepare variables for length (2D) or area (3D) of the interface
  vector<double> mycenterdisp(dim);
  vector<double> centerdisp(dim);
  double mylengthcirc = 0.0;
  double lengthcirc = 0.0;

  //calculating the center displacement by evaluating structure interface elements
  map<int, RefCountPtr<DRT::Element> >::const_iterator elemiter;
  for (elemiter = istructdispeles_.begin(); elemiter != istructdispeles_.end(); ++elemiter)
  {

    RefCountPtr<DRT::Element> iele = elemiter->second;
    vector<int> lm;
    vector<int> lmowner;
    vector<int> lmstride;
    iele->LocationVector(*structdis,lm,lmowner,lmstride);
    elevector2.Size(1);   //length of circ with gaussinteg
    elevector3.Size(dim);   //centerdisp part of ele

    params.set<string>("action","calc_struct_centerdisp");
    int err = iele->Evaluate(params,*structdis,lm,elematrix1,elematrix2,elevector1,elevector2,elevector3);
    if (err) dserror("error while evaluating elements");
    mylengthcirc += elevector2[0];

    //disp of the interface
    for (int i=0; i<dim ;i++)
    {
      mycenterdisp[i] += elevector3[i];
    }
  } //end of ele loop
  structdis->ClearState();

  //Communicate to 'assemble' length and center displacements
  comm.SumAll(&mylengthcirc, &lengthcirc, 1);
  comm.SumAll(&mycenterdisp[0], &centerdisp[0], dim);

  if (lengthcirc <= 1.0E-6)
    dserror("Zero interface length!");

  //calculating the final disp of the interface and summation over all time steps
  for (int i=0; i<dim ;i++)
  {
    centerdisp[i] = centerdisp[i] / lengthcirc;
    centerdisptotal_[i] +=centerdisp[i];
  }

  return centerdisp;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<int,LINALG::Matrix<3,1> > FSI::UTILS::SlideAleUtils::CurrentStructPos
(
    Teuchos::RCP<Epetra_Vector> reddisp,
    DRT::Discretization& interfacedis,
    double& maxcoord
)
{
  std::map<int,LINALG::Matrix<3,1> > currentpositions;
  map<int, Teuchos::RCP<DRT::Element> >::const_iterator eleiter;

  // map with fully reduced struct element distribution
  for (eleiter = structreduelements_.begin(); eleiter!=structreduelements_.end(); eleiter++)
  {
    Teuchos::RCP<DRT::Element> tmpele = eleiter->second;

    const int* n = tmpele->NodeIds();

    // fill currentpositions
    for (int j=0; j < tmpele->NumNode(); j++)
    {
      const int gid = n[j];
      const DRT::Node* node = interfacedis.gNode(gid);
      vector<int> lm;
      lm.reserve(3);
      // extract global dof ids
      interfacedis.Dof(node, lm);
      vector<double> mydisp(3);
      LINALG::Matrix<3,1> currpos;

      DRT::UTILS::ExtractMyValues(*reddisp,mydisp,lm);

      for (int a=0; a<3; a++)
      {
        currpos(a,0) = node->X()[a] + mydisp[a];
      }
      if (abs(currpos(2,0))>maxcoord)
        maxcoord=abs(currpos(2,0));
      currentpositions[node->Id()] = currpos;
    }
  }

  return currentpositions;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::UTILS::SlideAleUtils::SlideProjection
(
    ADAPTER::Structure& structure,
    Teuchos::RCP<DRT::Discretization> fluiddis,
    Teuchos::RCP<Epetra_Vector> idispale,
    Teuchos::RCP<Epetra_Vector> iprojdispale,
    ADAPTER::CouplingMortar& coupsf,
    const Epetra_Comm& comm

)
{
  const int dim = genprob.ndim;

  Teuchos::RCP<Epetra_Vector> idispnp = structure.ExtractInterfaceDispnp();

  // Redistribute displacement of structnodes on the interface to all processors.
  Teuchos::RCP<Epetra_Import> interimpo = rcp (new Epetra_Import(*structfullnodemap_,*structdofrowmap_));
  Teuchos::RCP<Epetra_Vector> reddisp = LINALG::CreateVector(*structfullnodemap_,true);
  reddisp -> Import(*idispnp,*interimpo,Add);

  DRT::Discretization& interfacedis = coupsf.Interface()->Discret();
  double rotrat=0.0;
  //currentpositions of struct nodes for the search tree (always 3 coordinates)
  std::map<int,LINALG::Matrix<3,1> > currentpositions =
      CurrentStructPos(reddisp, interfacedis, rotrat);

  //calculate structural interface center of gravity
  vector<double> centerdisp = Centerdisp(structure, comm);

  Teuchos::RCP<Epetra_Vector> frotfull = LINALG::CreateVector(*fluiddofrowmap_,true);
  if (aletype_==INPAR::FSI::ALEprojection_rot_z || aletype_==INPAR::FSI::ALEprojection_rot_zsphere)
  {
    Rotation(coupsf.Interface()->Discret(), idispale, comm, rotrat, frotfull);
  }

  // Project fluid nodes onto the struct interface
  //init of search tree
  Teuchos::RCP<GEO::SearchTree> searchTree = rcp(new GEO::SearchTree(0));
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofEles(structreduelements_, currentpositions);

  if(dim==2)
    searchTree->initializeTreeSlideALE(rootBox, structreduelements_, GEO::TreeType(GEO::QUADTREE));
  else if(dim==3)
    searchTree->initializeTreeSlideALE(rootBox, structreduelements_, GEO::TreeType(GEO::OCTTREE));
  else dserror("wrong dimension");

  // translation + projection
  map<int, DRT::Node*>::const_iterator nodeiter;
  for (nodeiter = ifluidslidnodes_.begin(); nodeiter != ifluidslidnodes_.end(); ++nodeiter)
  {

    DRT::Node* node = nodeiter->second;
    vector<int> lids(dim);
    for(int p=0; p<dim; p++)
    //lids of gids of node
    lids[p] = (fluiddofrowmap_)->LID((fluiddis->Dof(node))[p]);

    // current coord of ale node.
    // Initialize as coordinates of current node, which is extremely important for 2D!
    LINALG::Matrix<3,1> alenodecurr (node->X());

    // compute ALE position to project from
    if  (aletype_==INPAR::FSI::ALEprojection_curr)
    {
      // current coord of ale node = ref + centerdispincr + history
      for(int p=0; p<dim; p++)
        alenodecurr(p,0) =  (node->X()[p]) + centerdisp[p] + 1.0* (*iprojhist_)[(lids[p])];
    }
    else if (aletype_==INPAR::FSI::ALEprojection_ref)
    {
      // current coord of ale node = ref + centerdisp
      for(int p=0; p<dim; p++)
        alenodecurr(p,0) =   node->X()[p] + centerdisptotal_[p];
    }
    else if (aletype_==INPAR::FSI::ALEprojection_rot_z || aletype_==INPAR::FSI::ALEprojection_rot_zsphere)
    {
      // current coord of ale node = ref + centerdisp
      for(int p=0; p<dim; p++)
      {
        alenodecurr(p,0) =
          node->X()[p] + (*idispale)[(lids[p])] - 1.0*rotrat*(*frotfull)[(lids[p])];
      }
    }
    else
      dserror("you should not turn up here!");


    // final displacement of projection
    vector <double> finaldxyz(dim);

    //search for near elements next to the query point (ie within a radius of 2x maxmindist)
    std::map<int,std::set<int> >  closeeles =
        searchTree->searchElementsInRadius(interfacedis,currentpositions,alenodecurr,maxmindist_,0);
    //if no close elements could be found, try with a much larger radius and print a warning
    if (closeeles.empty())
    {
      cout<<"WARNING: no elements found in radius r="<<maxmindist_<<". Will try once with a bigger radius!"<<endl;
      closeeles = searchTree->searchElementsInRadius(interfacedis,currentpositions,alenodecurr,100.0*maxmindist_,0);
      maxmindist_ *= 10.0;

      // if still no element is found, complain about it!
      if (closeeles.empty())
        dserror("No elements in a large radius! Should not happen!");
    }
    //search for the nearest point to project on
    LINALG::Matrix<3,1> minDistCoords;
    if(dim == 2)
    {
      GEO::nearest2DObjectInNode(rcp(&interfacedis,false), structreduelements_, currentpositions,
          closeeles, alenodecurr, minDistCoords);
      finaldxyz[0] = minDistCoords(0,0) - node->X()[0];
      finaldxyz[1] = minDistCoords(1,0) - node->X()[1];
    }
    else
    {
      GEO::nearest3DObjectInNode(rcp(&interfacedis,false), structreduelements_, currentpositions,
          closeeles, alenodecurr, minDistCoords);
      finaldxyz[0] = minDistCoords(0,0) - node->X()[0];
      finaldxyz[1] = minDistCoords(1,0) - node->X()[1];
      finaldxyz[2] = minDistCoords(2,0) - node->X()[2];

    }

    //store displacement into parallel vector
    int err = iprojdispale->ReplaceMyValues(dim, &finaldxyz[0], &lids[0]);
    if (err == 1) dserror("error while replacing values");
  }
}

void FSI::UTILS::SlideAleUtils::RedundantElements
(
    ADAPTER::CouplingMortar& coupsf,
    const Epetra_Comm& comm
)
{
  // We need the structure elements (NOT THE MORTAR-ELEMENTS!) on every processor for the projection of the fluid nodes.
  // Furthermore we need the current position of the structnodes on every processor.
  // Elements provided by interface discretization, necessary maps provided by interface.

  int soffset = 0;
  int foffset = 0;
  if (structcoupmaster_)
  {
    structfullnodemap_ =  LINALG::AllreduceEMap(*(coupsf.Interface()->MasterRowDofs()));
    structfullelemap_ =  LINALG::AllreduceEMap(*(coupsf.Interface()->MasterRowElements()));
    fluidfullnodemap_ =  LINALG::AllreduceEMap(*(coupsf.Interface()->SlaveRowDofs()));
    fluidfullelemap_ =  LINALG::AllreduceEMap(*(coupsf.Interface()->SlaveRowElements()));
    soffset = 0;
    foffset = fluidfullelemap_->MinMyGID();
  }
  else
  {
    fluidfullnodemap_ =  LINALG::AllreduceEMap(*(coupsf.Interface()->MasterRowDofs()));
    fluidfullelemap_ =  LINALG::AllreduceEMap(*(coupsf.Interface()->MasterRowElements()));
    structfullnodemap_ =  LINALG::AllreduceEMap(*(coupsf.Interface()->SlaveRowDofs()));
    structfullelemap_ =  LINALG::AllreduceEMap(*(coupsf.Interface()->SlaveRowElements()));
    soffset = structfullelemap_->MinMyGID();
    foffset = 0;
  }

  DRT::Discretization& interfacedis = coupsf.Interface()->Discret();


  // build redundant version istructslideles_;
  map<int, Teuchos::RCP<DRT::Element> >::iterator mapit;
  vector<int> vstruslideleids; // vector for ele ids
  for (mapit = istructslideles_.begin();mapit != istructslideles_.end();mapit++)
  {
    vstruslideleids.push_back(mapit->first);
  }
  int globsum=0;
  int partsum=(vstruslideleids.size());
  comm.SumAll(&partsum,&globsum,1);
  //map with ele ids
  Epetra_Map mstruslideleids(globsum,vstruslideleids.size(),&(vstruslideleids[0]),0,comm);
  //redundant version of it
  Epetra_Map redmstruslideleids (*LINALG::AllreduceEMap(mstruslideleids));

  int dim = genprob.ndim;
  for (int eleind = 0; eleind<redmstruslideleids.NumMyElements(); eleind++)
  {
    {
      DRT::Element* tmpele = interfacedis.gElement(redmstruslideleids.GID(eleind)+soffset);
      if (dim == 3)
      {
        structreduelements_[tmpele->Id()]= rcp(new DRT::ELEMENTS::StructuralSurface(
          tmpele->Id(),tmpele->Owner(),tmpele->NumNode(),tmpele->NodeIds(),tmpele->Nodes(),&(*tmpele),0));
      }
      else if (dim == 2)
      {
        structreduelements_[tmpele->Id()]= rcp(new DRT::ELEMENTS::StructuralLine(
          tmpele->Id(),tmpele->Owner(),tmpele->NumNode(),tmpele->NodeIds(),tmpele->Nodes(),&(*tmpele),0));
      }
    }
  }


  //rebuild fluidslideeles_ with StructuralSurface Elements
  for (mapit = ifluidslideles_.begin();mapit != ifluidslideles_.end();mapit++)
  {
    DRT::Element* tmpele = interfacedis.gElement(mapit->first+foffset);
    if (dim == 3)
    {
      ifluidslidstructeles_[tmpele->Id()]= rcp(new DRT::ELEMENTS::StructuralSurface(
        tmpele->Id(),tmpele->Owner(),tmpele->NumNode(),tmpele->NodeIds(),tmpele->Nodes(),&(*tmpele),0));
    }
    else if (dim == 2)
    {
      ifluidslidstructeles_[tmpele->Id()]= rcp(new DRT::ELEMENTS::StructuralLine(
        tmpele->Id(),tmpele->Owner(),tmpele->NumNode(),tmpele->NodeIds(),tmpele->Nodes(),&(*tmpele),0));
    }
  }
}


void FSI::UTILS::SlideAleUtils::Rotation
(
    DRT::Discretization& mtrdis,      ///< fluid discretization
    Teuchos::RCP<Epetra_Vector> idispale,            ///< vector of ALE displacements
    const Epetra_Comm& comm,                         ///< communicator
    double& rotrat,                                  ///< rotation ratio of tangential displacements
    Teuchos::RCP<Epetra_Vector> rotfull              ///< vector of full displacements in tangential directions
)
{
  double maxcoord=rotrat;

  Teuchos::RCP<Epetra_Vector> idispstep = LINALG::CreateVector(*fluiddofrowmap_,false);
  idispstep->Update(1.0, *idispale, -1.0, *iprojhist_, 0.0);

  // get structure and fluid discretizations  and set state for element evaluation
  const Teuchos::RCP<Epetra_Vector> idispstepcol = LINALG::CreateVector(*mtrdis.DofColMap(),false);
  LINALG::Export(*idispstep,*idispstepcol);
  const Teuchos::RCP<Epetra_Vector> idispnpcol = LINALG::CreateVector(*mtrdis.DofColMap(),false);
  LINALG::Export(*idispale,*idispnpcol);

  mtrdis.SetState("displacementnp",idispnpcol);
  mtrdis.SetState("displacementincr",idispstepcol);

  //prepare variables for length (2D) or area (3D) of the interface
  double myrotation = 0.0;
  double rotation = 0.0;
  double mylengthcirc = 0.0;
  double lengthcirc = 0.0;

  //calculating the center displacement by evaluating structure interface elements
  map<int, RefCountPtr<DRT::Element> >::const_iterator elemiter;
  for (elemiter = ifluidslidstructeles_.begin(); elemiter != ifluidslidstructeles_.end(); ++elemiter)
  {
    //define stuff needed by the elements
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2;
    Epetra_SerialDenseVector elevector3;
    Teuchos::ParameterList params;

    RefCountPtr<DRT::Element> iele = elemiter->second;
    vector<int> lm;
    vector<int> lmowner;
    vector<int> lmstride;
    iele->LocationVector(mtrdis,lm,lmowner,lmstride);
    elevector2.Size(1);   //circumference (2D) or surface area (3D) of the considered elements
    elevector3.Size(1);   //normalized displacement in tangential direction ('rotation')

    params.set<string>("action","calc_struct_rotation");
    params.set<double>("maxcoord",maxcoord);
    params.set<INPAR::FSI::SlideALEProj>("aletype",aletype_);
    int err = iele->Evaluate(params,mtrdis,lm,elematrix1,elematrix2,elevector1,elevector2,elevector3);
    if (err)
      dserror("error while evaluating elements");

    mylengthcirc += elevector2[0];
    //disp of the interface
    myrotation += elevector3[0];
  } //end of ele loop

  //Communicate to 'assemble' length and center displacements
  comm.SumAll(&mylengthcirc, &lengthcirc, 1);
  comm.SumAll(&myrotation, &rotation, 1);

  if (lengthcirc <= 1.0E-6)
    dserror("Zero interface length!");

  //calculating the final disp of the interface and summation over all time steps
  rotrat = rotation / lengthcirc;

  // second round!
  // compute correction displacement to account for rotation
  for (elemiter = ifluidslidstructeles_.begin(); elemiter != ifluidslidstructeles_.end(); ++elemiter)
  {
    //define stuff needed by the elements
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2;
    Epetra_SerialDenseVector elevector3;
    Teuchos::ParameterList params;

    RefCountPtr<DRT::Element> iele = elemiter->second;
    vector<int> lm;
    vector<int> lmowner;
    vector<int> lmstride;
    iele->LocationVector(mtrdis,lm,lmowner,lmstride);
    elevector1.Size(lm.size());

    params.set<string>("action","calc_undo_struct_rotation");
    params.set<double>("maxcoord",maxcoord);
    params.set<INPAR::FSI::SlideALEProj>("aletype",aletype_);
    int err = iele->Evaluate(params,mtrdis,lm,elematrix1,elematrix2,elevector1,elevector2,elevector3);
    if (err)
      dserror("error while evaluating elements");

    LINALG::Assemble(*rotfull,elevector1,lm,lmowner);
  }

  mtrdis.ClearState();

  return;
}

void FSI::UTILS::SlideAleUtils::OutputRestart
(
  IO::DiscretizationWriter& output
)
{
  output.WriteVector("projhist", iprojhist_);

  return;
}

void FSI::UTILS::SlideAleUtils::ReadRestart
(
    IO::DiscretizationReader& reader
)
{
  reader.ReadVector(iprojhist_, "projhist");
}

#endif

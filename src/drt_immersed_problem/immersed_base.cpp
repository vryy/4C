/*!----------------------------------------------------------------------
\file immersed_base.cpp

\brief base class for all immersed algorithms

<pre>
Maintainers: Andreas Rauch & Anh-Tu Vuong
             {rauch,vuong}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289--15240 / 15264
</pre>
*----------------------------------------------------------------------*/
#include "immersed_base.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_calc_utils.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_element.H"
#include "../drt_fluid_ele/fluid_ele_immersed.H"
#include "../drt_fluid_ele/fluid_ele_action.H"

IMMERSED::ImmersedBase::ImmersedBase()
{

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<int> IMMERSED::ImmersedBase::DetermineImmersionDomain(Teuchos::RCP<DRT::Discretization> backgrounddis, Teuchos::RCP<DRT::Discretization> immerseddis, bool firstcall)
{
# ifdef DEBUG
  int notconvergedcounter = 0;
# endif

  if(backgrounddis->Comm().MyPID() == 0)
  {
    std::cout<<"################################################################################################"<<std::endl;
    std::cout<<"###   Determine " << backgrounddis->Name() <<" elements in which the "<<immerseddis->Name()<<" is immersed ..."<<std::endl;
    std::cout<<"################################################################################################"<<std::endl;
  }

  std::set<int> nodeset;

  // get gids of column elements of backgrounddis
  const Epetra_Map* backgroundelecolmap = backgrounddis->ElementColMap();
  int myglobalcolelementsize = backgroundelecolmap->NumMyElements();
  std::vector<int> myglobalcolelements(myglobalcolelementsize);
  backgroundelecolmap->MyGlobalElements(&myglobalcolelements[0]);

  // get possible elements being intersected by immersed structure
  DRT::Condition* searchbox = backgrounddis->GetCondition("ImmersedSearchbox");
  std::map<int,Teuchos::RCP<DRT::Element> >& searchboxgeom = searchbox->Geometry();

  // get node ids of immersed discretization
  const Epetra_Map* noderowmap = immerseddis->NodeColMap();
  int mynoderowmapsize = noderowmap ->NumMyElements();
  std::vector<int> myglobalelements(mynoderowmapsize);
  noderowmap->MyGlobalElements(&myglobalelements[0]);

  std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;
  std::vector<double> my_displacements_np;
  double xi[DRT::Problem::Instance()->NDim()];
  double x [DRT::Problem::Instance()->NDim()];
  int inode = 0;

  Teuchos::RCP<const Epetra_Vector> displacements_np;
  if(firstcall)
    displacements_np = Teuchos::rcp(new const Epetra_Vector(*immerseddis->DofColMap(),true));
  else
    displacements_np = immerseddis->GetState("dispnp");

  //std::cout<<*immerseddis->DofRowMap()<<std::endl;
  /////////////////////////////////////////////////////
  // loop over all immersed nodes
  /////////////////////////////////////////////////////
  for(int i=0;i<mynoderowmapsize;i++)
  {//std::cout<<"i="<<i<<std::endl;
    DRT::Node* immersednode = immerseddis->gNode(myglobalelements[i]);
    if(immersednode == NULL)
      dserror("Could not get node with GID %d",immersednode->Id());
    const double* X = immersednode -> X();
    DRT::Element* adjacentelement = immersednode->Elements()[0];
    if(adjacentelement == NULL)
      dserror("Could not get adjacent element to node with GID %d",immersednode->Id());

    adjacentelement->LocationVector(*immerseddis,lm,lmowner,lmstride);
    my_displacements_np.resize(lm.size());
    if((int)my_displacements_np.size() != adjacentelement->NumNode()*immerseddis->NumDof(immersednode))
      dserror("my_displacements_np has less capacity than the numnode*numdofpernode");
    DRT::UTILS::ExtractMyValues(*displacements_np,my_displacements_np,lm);

    // get node id on adjacentelement of immersednode
    for (int j=0;j<adjacentelement->NumNode();++j)
    {
      if (adjacentelement->NodeIds()[j]==immersednode->Id())
        inode = j;
    }
    // update node position X -> x
    {
      for (int idof=0;idof<DRT::Problem::Instance()->NDim();++idof)
      {
        x[idof] = X[idof] + my_displacements_np[inode*immerseddis->NumDof(immersednode)+idof];
      }
    }
    lmowner.clear();
    lmstride.clear();
    lm.clear();
    ////////////////////////////////////////////
    // loop over all background elements
    ////////////////////////////////////////////
    for (curr=searchboxgeom.begin(); curr!=searchboxgeom.end(); ++curr)
    {
      bool converged = false;

      DRT::Element::DiscretizationType distype = immerseddis->gElement(0)->Shape();
      switch(distype)
      {
      case DRT::Element::hex8 :
      {
        MORTAR::UTILS::GlobalToLocal<DRT::Element::hex8>(*(curr->second),&x[0],&xi[0],converged);
        break;
      }
      default:
      {
        dserror("DISTYPE NOT SUPPORTED YET. PLEASE CREATE ENTRY IN THIS SWITCH-CASE STATEMENT");
        break;
      }
      }

# ifdef DEBUG
      if(!converged)
      {
        notconvergedcounter ++;
//        std::cout<<" Map immersed node with GID "<<immerseddis->gNode(myglobalelements[i])->Id()<<" to element with GID "<<curr->second->Id()<<std::endl;
//        dserror("MAPPING FROM IMMERSED NODE TO BACKGROUNDELEMENT DID NOT CONVERGE");
      }
# endif

      if ((abs(xi[0])-1.0)<1e-12 and (abs(xi[1])-1.0)<1e-12 and (abs(xi[2])-1.0)<1e-12)
      {// node i lies in element curr
        Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidImmersed>(curr->second)->SetIsImmersed(1);
        const int* nodes;
        nodes = curr->second->NodeIds();
        for(int k=0;k<curr->second->NumNode();++k)
        {
          nodeset.insert(nodes[k]);
        }
      }

    }// loop over background elements in searchbox
  }// loop over immersed nodes

  std::vector<int> nodevector(nodeset.size()); // variable to return
  std::copy(nodeset.begin(), nodeset.end(), nodevector.begin());


# ifdef DEBUG
    std::cout<<"PROC "<<backgrounddis->Comm().MyPID()<<" "<<notconvergedcounter<<" mappings did not converge"<<std::endl;
    std::cout<<"PROC "<<backgrounddis->Comm().MyPID()<<" : searchboxgeom.size() = "<<searchboxgeom.size()<<std::endl;
    std::cout<<"PROC "<<backgrounddis->Comm().MyPID()<<" : identified "<<nodeset.size()<<" nodes in immersion domain."<<std::endl;
    if(nodeset.size() != nodevector.size())
      dserror("nodeset and nodevector must have same size");
#endif

  return nodevector;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<int> IMMERSED::ImmersedBase::DetermineImmersionBoundaryDomain(Teuchos::RCP<DRT::Discretization> backgrounddis, Teuchos::RCP<DRT::Discretization> immerseddis,  const std::string& condname, bool firstcall)
{
# ifdef DEBUG
  int notconvergedcounter = 0;
# endif

  if(backgrounddis->Comm().MyPID() == 0)
  {
    std::cout<<"################################################################################################"<<std::endl;
    std::cout<<"###   Determine " << backgrounddis->Name() <<" elements in which the "<<immerseddis->Name()<<" boundary is immersed ..."<<std::endl;
    std::cout<<"################################################################################################"<<std::endl;
  }

  std::set<int> nodeset;

  const Epetra_Comm& comm = backgrounddis->Comm();

  DRT::Condition* immersedcond = immerseddis->GetCondition(condname);
  std::map<int,Teuchos::RCP<DRT::Element> >& immersedgeom = immersedcond->Geometry();

  // get gids of column elements of backgrounddis
  const Epetra_Map* backgroundelecolmap = backgrounddis->ElementColMap();
  int myglobalcolelementsize = backgroundelecolmap->NumMyElements();
  std::vector<int> myglobalcolelements(myglobalcolelementsize);
  backgroundelecolmap->MyGlobalElements(&myglobalcolelements[0]);

  // get possible elements being intersected by immersed structure
  DRT::Condition* searchbox = backgrounddis->GetCondition("ImmersedSearchbox");
  std::map<int,Teuchos::RCP<DRT::Element> >& searchboxgeom = searchbox->Geometry();

  std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
  std::map<int,Teuchos::RCP<DRT::Element> >::iterator geomcurr;
  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;
  std::vector<double> my_displacements_np;
  double xi[DRT::Problem::Instance()->NDim()];
  double x [DRT::Problem::Instance()->NDim()];
  int inode = 0;

  Teuchos::RCP<const Epetra_Vector> displacements_np;
  if(firstcall)
    displacements_np = Teuchos::rcp(new const Epetra_Vector(*immerseddis->DofColMap(),true));
  else
    displacements_np = immerseddis->GetState("dispnp");

  //std::cout<<*immerseddis->DofRowMap()<<std::endl;
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // loop over all immersed boundary geometry and fill vector with unique node ids
  /////////////////////////////////////////////////////////////////////////////////////////////////
  std::set<int> boundarynodeids;
  for (geomcurr=immersedgeom.begin(); geomcurr!=immersedgeom.end(); ++geomcurr)
  {
    for (int i=0;i<geomcurr->second->NumNode();++i)
    {
      DRT::Node* immersednode = geomcurr->second->Nodes()[i];
      if(immersednode == NULL)
        dserror("Could not get node with GID %d",immersednode->Id());
      boundarynodeids.insert(immersednode->Id());
    }
  }
  // gather global vector from proc local vectors
  LINALG::GatherAll(boundarynodeids,comm);
  std::vector<int> boundarynodevector(boundarynodeids.size());
  std::copy(boundarynodeids.begin(), boundarynodeids.end(), boundarynodevector.begin());

#ifdef DEBUG
  int mysize = boundarynodevector.size();
  int globalsize = 0;
  comm.SumAll(&mysize,&globalsize,1);
  std::cout<<"PROC "<<comm.MyPID()<<" : "<< "local number of bundarynodes = "<<mysize<<std::endl;
  std::cout<<"PROC "<<comm.MyPID()<<" : "<< "global number of bundarynodes = "<<globalsize<<std::endl;
  std::cout<<"PROC "<<comm.MyPID()<<" : "<< "global number of gathererd bundarynodes = "<<boundarynodeids.size()<<std::endl;
#endif


  ///////////////////////////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////
  // loop over all immersed boundary nodes
  /////////////////////////////////////////////////////
  for(int i=0;i< (int)boundarynodeids.size();i++)
  {//std::cout<<"i="<<i<<std::endl;
    DRT::Node* immersednode = immerseddis->gNode(boundarynodevector[i]);
    if(immersednode == NULL)
      dserror("Could not get node with GID %d",immersednode->Id());
    const double* X = immersednode -> X();
    DRT::Element* adjacentelement = immersednode->Elements()[0];
    if(adjacentelement == NULL)
      dserror("Could not get adjacent element to node with GID %d",immersednode->Id());

    adjacentelement->LocationVector(*immerseddis,lm,lmowner,lmstride);
    my_displacements_np.resize(lm.size());
    if((int)my_displacements_np.size() != adjacentelement->NumNode()*immerseddis->NumDof(immersednode))
      dserror("my_displacements_np has less capacity than the numnode*numdofpernode");
    DRT::UTILS::ExtractMyValues(*displacements_np,my_displacements_np,lm);

    // get node id on adjacentelement of immersednode
    for (int j=0;j<adjacentelement->NumNode();++j)
    {
      if (adjacentelement->NodeIds()[j]==immersednode->Id())
        inode = j;
    }
    // update node position X -> x
    {
      for (int idof=0;idof<DRT::Problem::Instance()->NDim();++idof)
      {
        x[idof] = X[idof] + my_displacements_np[inode*immerseddis->NumDof(immersednode)+idof];
      }
    }
    lmowner.clear();
    lmstride.clear();
    lm.clear();
    ////////////////////////////////////////////
    // loop over all background elements
    ////////////////////////////////////////////
    for (curr=searchboxgeom.begin(); curr!=searchboxgeom.end(); ++curr)
    {
      bool converged = false;

      DRT::Element::DiscretizationType distype = immerseddis->gElement(0)->Shape();
      switch(distype)
      {
      case DRT::Element::hex8 :
      {
        MORTAR::UTILS::GlobalToLocal<DRT::Element::hex8>(*(curr->second),&x[0],&xi[0],converged);
        break;
      }
      default:
      {
        dserror("DISTYPE NOT SUPPORTED YET. PLEASE CREATE ENTRY IN THIS SWITCH-CASE STATEMENT");
        break;
      }
      }

# ifdef DEBUG
      if(!converged)
      {
        notconvergedcounter ++;
//        std::cout<<" Map immersed node with GID "<<immerseddis->gNode(myglobalelements[i])->Id()<<" to element with GID "<<curr->second->Id()<<std::endl;
//        dserror("MAPPING FROM IMMERSED NODE TO BACKGROUNDELEMENT DID NOT CONVERGE");
      }
# endif

      if ((abs(xi[0])-1.0)<1e-12 and (abs(xi[1])-1.0)<1e-12 and (abs(xi[2])-1.0)<1e-12)
      {// node i lies in element curr
        Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidImmersed>(curr->second)->SetIsImmersedBoundary(1);
        const int* nodes;
        nodes = curr->second->NodeIds();
        for(int k=0;k<curr->second->NumNode();++k)
        {
          nodeset.insert(nodes[k]);
        }
      }

    }// loop over background elements in searchbox
  }// loop over immersed nodes

  std::vector<int> nodevector(nodeset.size()); // variable to return
  std::copy(nodeset.begin(), nodeset.end(), nodevector.begin());


# ifdef DEBUG
    std::cout<<"PROC "<<comm.MyPID()<<" "<<notconvergedcounter<<" mappings did not converge"<<std::endl;
    std::cout<<"PROC "<<comm.MyPID()<<" : searchboxgeom.size() = "<<searchboxgeom.size()<<std::endl;
    std::cout<<"PROC "<<comm.MyPID()<<" : identified "<<nodeset.size()<<" nodes in immersion domain."<<std::endl;
    if(nodeset.size() != nodevector.size())
      dserror("nodeset and nodevector must have same size");
#endif

  return nodevector;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::CreateGhosting(Teuchos::RCP<DRT::Discretization> distobeghosted)
{
  if(distobeghosted->Comm().MyPID() == 0)
  {
    std::cout<<"################################################################################################"<<std::endl;
    std::cout<<"###   Ghost discretization "<<distobeghosted->Name()<<" redundantly on all procs ... "<<std::endl;
    std::cout<<"################################################################################################"<<std::endl;
  }

  std::vector<int> allproc(distobeghosted->Comm().NumProc());
  for (int i=0; i<distobeghosted->Comm().NumProc(); ++i) allproc[i] = i;

  // fill my own row node ids
  const Epetra_Map* noderowmap = distobeghosted->NodeRowMap();
  std::vector<int> sdata;
  for (int i=0; i<noderowmap->NumMyElements(); ++i)
  {
    int gid = noderowmap->GID(i);
    DRT::Node* node = distobeghosted->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    sdata.push_back(gid);
  }

  // gather all master row node gids redundantly
  std::vector<int> rdata;
  LINALG::Gather<int>(sdata,rdata,(int)allproc.size(),&allproc[0],distobeghosted->Comm());

  // build new node column map (on ALL processors)
  Teuchos::RCP<Epetra_Map> newnodecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)rdata.size(),&rdata[0],0,distobeghosted->Comm()));
  sdata.clear();
  rdata.clear();

  // fill my own row element ids
  const Epetra_Map* elerowmap  = distobeghosted->ElementRowMap();
  sdata.resize(0);
  for (int i=0; i<elerowmap->NumMyElements(); ++i)
  {
    int gid = elerowmap->GID(i);
    DRT::Element* ele = distobeghosted->gElement(gid);
    if (!ele) dserror("ERROR: Cannot find element with gid %",gid);
    sdata.push_back(gid);
  }

  // gather all gids of elements redundantly
  rdata.resize(0);
  LINALG::Gather<int>(sdata,rdata,(int)allproc.size(),&allproc[0],distobeghosted->Comm());

  // build new element column map (on ALL processors)
  Teuchos::RCP<Epetra_Map> newelecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)rdata.size(),&rdata[0],0,distobeghosted->Comm()));
  sdata.clear();
  rdata.clear();
  allproc.clear();

  // redistribute the discretization of the interface according to the
  // new node / element column layout (i.e. master = full overlap)
  distobeghosted->ExportColumnNodes(*newnodecolmap);
  distobeghosted->ExportColumnElements(*newelecolmap);

  distobeghosted->FillComplete();

#ifdef DEBUG
  int nummycolnodes = newnodecolmap->NumMyElements();
  int sizelist[distobeghosted->Comm().NumProc()];
  distobeghosted->Comm().GatherAll(&nummycolnodes,&sizelist[0],1);
  std::cout<<"PROC "<<distobeghosted->Comm().MyPID()<<" : "<<nummycolnodes<<" colnodes"<<std::endl;
  distobeghosted->Comm().Barrier(); // wait for procs
  for(int k=1;k<distobeghosted->Comm().NumProc();++k)
  {
    if(sizelist[k-1]!=nummycolnodes)
      dserror("Since whole dis is ghosted every processor should have the same number of colnodes. This is not the case! Fix this!");
  }
#endif

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::CreateVolumeCondition(Teuchos::RCP<DRT::Discretization> dis, std::vector<int> dvol_fenode, DRT::Condition::ConditionType condtype, std::string condname)
{
  // determine id of condition
  std::multimap<std::string,Teuchos::RCP<DRT::Condition> > allconditions;
  allconditions = dis->GetAllConditions();
  int id = (int)allconditions.size();
  id += 1;

  // build condition
  bool buildgeometry = true; // needed for now to check number of elements in neumannnaumann.cpp
  Teuchos::RCP<DRT::Condition> condition =
          Teuchos::rcp(new DRT::Condition(id,condtype,buildgeometry,DRT::Condition::Volume));

  // add nodes to conditions
   condition->Add("Node Ids",dvol_fenode);

   // add condition to discretization
   dis->SetCondition(condname,condition);

   // fill complete if necessary
   if (!dis->Filled())
     dis -> FillComplete();

   //debug
#ifdef DEBUG
   std::cout<<"PROC "<<dis->Comm().MyPID()<<" : Number of conditioned elements: "<<dis->GetCondition(condname)->Geometry().size()<<" ("<<condname<<")"<<std::endl;
#endif

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::UpdateVolumeCondition(Teuchos::RCP<DRT::Discretization> dis, std::vector<int> dvol_fenode, DRT::Condition::ConditionType condtype, std::string condname)
{

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedBase::InterpolateToImmersedNodes(Teuchos::RCP<DRT::Discretization> backgrounddis, std::vector<int> dvol_fenode, std::string condname)
{
  Teuchos::ParameterList params;
  params.set<int>("action", FLD::interpolatetoimmersednodes);
  DRT::AssembleStrategy strategy(
      0,              // dofset for row
      0,              // dofset for column
      Teuchos::null,  // matrix 1
      Teuchos::null,  //
      Teuchos::null,  // vector 1 (to be filled)
      Teuchos::null,  //
      Teuchos::null   //
  );

  backgrounddis->EvaluateCondition(params, strategy, condname);

  return Teuchos::null;
}

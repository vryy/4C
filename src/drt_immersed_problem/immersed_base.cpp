/*!----------------------------------------------------------------------
\file immersed_base.cpp

\brief base class for all immersed algorithms

\level 2

\maintainer Andreas Rauch
            http://www.lnm.mw.tum.de
            089 - 289 -15240
*----------------------------------------------------------------------*/
#include "immersed_base.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_element.H"
#include "../drt_fluid_ele/fluid_ele_immersed.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_geometry/position_array.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_lib/drt_element.H"
#include "../linalg/linalg_utils.H"


IMMERSED::ImmersedBase::ImmersedBase()
{

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
//std::vector<int> IMMERSED::ImmersedBase::DetermineImmersionDomain(Teuchos::RCP<DRT::Discretization> backgrounddis, Teuchos::RCP<DRT::Discretization> immerseddis,DRT::AssembleStrategy* strategy, bool firstcall)
//{
//# ifdef DEBUG
//  int notconvergedcounter = 0;
//# endif
//
//  if(backgrounddis->Comm().MyPID() == 0)
//  {
//    std::cout<<"################################################################################################"<<std::endl;
//    std::cout<<"###   Determine " << backgrounddis->Name() <<" elements in which the "<<immerseddis->Name()<<" is immersed ..."<<std::endl;
//    std::cout<<"################################################################################################"<<std::endl;
//  }
//
//  std::set<int> nodeset;
//
//  // get gids of column elements of backgrounddis
//  const Epetra_Map* backgroundelecolmap = backgrounddis->ElementColMap();
//  int myglobalcolelementsize = backgroundelecolmap->NumMyElements();
//  std::vector<int> myglobalcolelements(myglobalcolelementsize);
//  backgroundelecolmap->MyGlobalElements(&myglobalcolelements[0]);
//
//  // pointer to background element
//  Teuchos::RCP<DRT::Element> ele;
//
//  // get possible elements being intersected by immersed structure
//  DRT::Condition* searchbox = backgrounddis->GetCondition("ImmersedSearchbox");
//  std::map<int,Teuchos::RCP<DRT::Element> >& searchboxgeom = searchbox->Geometry();
//
//  // get node ids of immersed discretization
//  const Epetra_Map* nodecolmap = immerseddis->NodeColMap();
//  int mynoderowmapsize = nodecolmap ->NumMyElements();
//  std::vector<int> myglobalelements(mynoderowmapsize);
//  nodecolmap->MyGlobalElements(&myglobalelements[0]);
//
//  std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
//  std::vector<int> lm;
//  std::vector<int> lmowner;
//  std::vector<int> lmstride;
//  std::vector<double> my_displacements_np;
//  double xi[DRT::Problem::Instance()->NDim()];
//  double x [DRT::Problem::Instance()->NDim()];
//  int inode = 0;
//
//  Teuchos::RCP<const Epetra_Vector> displacements_np;
//  if(firstcall)
//    displacements_np = Teuchos::rcp(new const Epetra_Vector(*immerseddis->DofColMap(),true));
//  else
//    displacements_np = immerseddis->GetState("dispnp");
//
//  //std::cout<<*immerseddis->DofRowMap()<<std::endl;
//  /////////////////////////////////////////////////////
//  // loop over all immersed nodes
//  /////////////////////////////////////////////////////
//  for(int i=0;i<mynoderowmapsize;i++)
//  {//std::cout<<"i="<<i<<std::endl;
//    DRT::Node* immersednode = immerseddis->gNode(myglobalelements[i]);
//
//#ifdef DEBUG
//    if(immersednode == NULL)
//      dserror("Could not get node with GID %d",immersednode->Id());
//#endif
//    // get initial coordinates and arbitrary [0]-th adjacent element to this node
//    const double* X = immersednode -> X();
//    DRT::Element* adjacentelement = immersednode->Elements()[0];
//
//#ifdef DEBUG
//    if(adjacentelement == NULL)
//      dserror("Could not get adjacent element to node with GID %d",immersednode->Id());
//#endif
//
//    // get data from adjacent element
//    adjacentelement->LocationVector(*immerseddis,lm,lmowner,lmstride);
//    my_displacements_np.resize(lm.size());
//
//#ifdef DEBUG
//    if((int)my_displacements_np.size() != adjacentelement->NumNode()*immerseddis->NumDof(immersednode))
//      dserror("my_displacements_np has less capacity than the numnode*numdofpernode");
//#endif
//
//    // get displacements from adjacent element
//    DRT::UTILS::ExtractMyValues(*displacements_np,my_displacements_np,lm);
//
//    // get node id on adjacentelement of immersednode
//    for (int j=0;j<adjacentelement->NumNode();++j)
//    {
//      if (adjacentelement->NodeIds()[j]==immersednode->Id())
//        inode = j;
//    }
//
//    // update node position X -> x
//    {
//      for (int idof=0;idof<DRT::Problem::Instance()->NDim();++idof)
//      {
//        x[idof] = X[idof] + my_displacements_np[inode*immerseddis->NumDof(immersednode)+idof];
//      }
//    }
//    lmowner.clear();
//    lmstride.clear();
//    lm.clear();
//
//    ////////////////////////////////////////////
//    // loop over all background elements
//    //
//    //
//    ////////////////////////////////////////////
//    for (curr=searchboxgeom.begin(); curr!=searchboxgeom.end(); ++curr)
//    {
//      bool converged = false;
//      //std::cout<<"PROC "<<backgrounddis->Comm().MyPID()<<" : "<<colele<<std::endl;
//      ele=curr->second;
//
//      // get shape
//      DRT::Element::DiscretizationType distype = immerseddis->gElement(0)->Shape();
//
//      switch(distype)
//      {
//      case DRT::Element::hex8 :
//      {
//        MORTAR::UTILS::GlobalToLocal<DRT::Element::hex8>(*ele,&x[0],&xi[0],converged);
//        break;
//      }
//      default:
//      {
//        dserror("DISTYPE NOT SUPPORTED YET. PLEASE CREATE ENTRY IN THIS SWITCH-CASE STATEMENT");
//        break;
//      }
//      }
//
//# ifdef DEBUG
//      if(!converged)
//      {
//        notconvergedcounter ++;
////        std::cout<<" Map immersed node with GID "<<immerseddis->gNode(myglobalelements[i])->Id()<<" to element with GID "<<curr->second->Id()<<std::endl;
////        dserror("MAPPING FROM IMMERSED NODE TO BACKGROUNDELEMENT DID NOT CONVERGE");
//      }
//# endif
//
//      if ((abs(xi[0])-1.0)<1e-12 and (abs(xi[1])-1.0)<1e-12 and (abs(xi[2])-1.0)<1e-12)
//      {
//        // -> node i lies in element background element
//        Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidImmersedBase>(ele)->SetIsImmersed(1);
//
//        // fill nodeset with node ids of background eles
//        const int* nodes;
//        nodes = ele->NodeIds();
//        for(int k=0;k<ele->NumNode();++k)
//        {
//          nodeset.insert(nodes[k]);
//        }
//      }
//
//    }// loop over background elements in searchbox
//  }// loop over immersed nodes
//
//  std::vector<int> nodevector(nodeset.size()); // variable to return
//  std::copy(nodeset.begin(), nodeset.end(), nodevector.begin());
//
//
//# ifdef DEBUG
//    std::cout<<"PROC "<<backgrounddis->Comm().MyPID()<<" "<<notconvergedcounter<<" mappings did not converge"<<std::endl;
//    std::cout<<"PROC "<<backgrounddis->Comm().MyPID()<<" : searchboxgeom.size() = "<<searchboxgeom.size()<<std::endl;
//    std::cout<<"PROC "<<backgrounddis->Comm().MyPID()<<" : identified "<<nodeset.size()<<" nodes in immersion domain."<<std::endl;
//    if(nodeset.size() != nodevector.size())
//      dserror("nodeset and nodevector must have same size");
//#endif
//
//  return nodevector;
//}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
//std::vector<int> IMMERSED::ImmersedBase::DetermineImmersionBoundaryDomain(Teuchos::RCP<DRT::Discretization> backgrounddis,
//                                                                          Teuchos::RCP<DRT::Discretization> immerseddis,
//                                                                          const std::string& condname,
//                                                                          bool gpversion,
//                                                                          bool firstcall,
//                                                                          const Teuchos::RCP<std::map<int,Teuchos::RCP<std::map<int,Teuchos::RCP<std::vector<double> > > > > >& gpmap)
//{
//# ifdef DEBUG
//  int notconvergedcounter = 0;
//# endif
//
//  if(backgrounddis->Comm().MyPID() == 0)
//  {
//    std::cout<<"################################################################################################"<<std::endl;
//    std::cout<<"###   Determine " << backgrounddis->Name() <<" elements in which the "<<immerseddis->Name()<<" boundary is immersed ..."<<std::endl;
//    std::cout<<"################################################################################################"<<std::endl;
//  }
//
//  std::set<int> nodeset;
//  double xigp = 0.577350269; // 1/sqrt(3)
//
//  const Epetra_Comm& comm = backgrounddis->Comm();
//
//  DRT::Condition* immersedcond = immerseddis->GetCondition(condname);
//  std::map<int,Teuchos::RCP<DRT::Element> >& immersedgeom = immersedcond->Geometry();
//
//  // get gids of column elements of backgrounddis
//  const Epetra_Map* backgroundelecolmap = backgrounddis->ElementColMap();
//  int myglobalcolelementsize = backgroundelecolmap->NumMyElements();
//  std::vector<int> myglobalcolelements(myglobalcolelementsize);
//  backgroundelecolmap->MyGlobalElements(&myglobalcolelements[0]);
//
//  // get possible elements being intersected by immersed structure
//  DRT::Condition* searchbox = backgrounddis->GetCondition("ImmersedSearchbox");
//  std::map<int,Teuchos::RCP<DRT::Element> >& searchboxgeom = searchbox->Geometry();
//
//  std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
//  std::map<int,Teuchos::RCP<DRT::Element> >::iterator geomcurr;
//  std::vector<int> lm;
//  std::vector<int> lmowner;
//  std::vector<int> lmstride;
//  std::vector<double> my_displacements_np;
//  double xi[DRT::Problem::Instance()->NDim()];
//  double x [DRT::Problem::Instance()->NDim()];
//  int inode = 0;
//
//  Teuchos::RCP<const Epetra_Vector> displacements_np;
//  if(firstcall)
//    displacements_np = Teuchos::rcp(new const Epetra_Vector(*immerseddis->DofColMap(),true));
//  else
//    displacements_np = immerseddis->GetState("dispnp");
//
//
//  if (gpversion == false)
//  {
//    std::set<int> boundarynodeids;
//
//    //////////////////////////////////////////////////////////////////////////////////////////////////
//    // loop over all immersed boundary geometry and fill vector with unique node ids
//    // and fill vector with global coordinates of element integration points
//    /////////////////////////////////////////////////////////////////////////////////////////////////
//    for (geomcurr=immersedgeom.begin(); geomcurr!=immersedgeom.end(); ++geomcurr)
//    {
//      for (int i=0;i<geomcurr->second->NumNode();++i)
//      {
//        DRT::Node* immersednode = geomcurr->second->Nodes()[i];
//        if(immersednode == NULL)
//          dserror("Could not get node with GID %d",immersednode->Id());
//        boundarynodeids.insert(immersednode->Id());
//      }
//    }
//    // gather global vector from proc local vectors
//    LINALG::GatherAll(boundarynodeids,comm);
//
//    std::vector<int> boundarynodevector(boundarynodeids.size());
//    std::copy(boundarynodeids.begin(), boundarynodeids.end(), boundarynodevector.begin());
//
//#ifdef DEBUG
//    int mysize = boundarynodevector.size();
//    int globalsize = 0;
//    comm.SumAll(&mysize,&globalsize,1);
//    std::cout<<"PROC "<<comm.MyPID()<<" : "<< "local number of bundarynodes = "<<mysize<<std::endl;
//    std::cout<<"PROC "<<comm.MyPID()<<" : "<< "global number of bundarynodes = "<<globalsize<<std::endl;
//    std::cout<<"PROC "<<comm.MyPID()<<" : "<< "global number of gathererd bundarynodes = "<<boundarynodeids.size()<<std::endl;
//#endif
//
//    ///////////////////////////////////////////////////////////////////////////////////////////////////
//
//    std::vector<int> boundaryelementids;
//    std::vector<int> boundaryelementinpointglobalcoords;
//
//    /////////////////////////////////////////////////////
//    // loop over all immersed boundary nodes
//    /////////////////////////////////////////////////////
//    for(int i=0;i< (int)boundarynodeids.size();i++)
//    {
//      DRT::Node* immersednode = immerseddis->gNode(boundarynodevector[i]);
//      if(immersednode == NULL)
//        dserror("Could not get node with GID %d",immersednode->Id());
//      const double* X = immersednode -> X();
//      DRT::Element* adjacentelement = immersednode->Elements()[0];
//      if(adjacentelement == NULL)
//        dserror("Could not get adjacent element to node with GID %d",immersednode->Id());
//
//      adjacentelement->LocationVector(*immerseddis,lm,lmowner,lmstride);
//      my_displacements_np.resize(lm.size());
//      if((int)my_displacements_np.size() != adjacentelement->NumNode()*immerseddis->NumDof(immersednode))
//        dserror("my_displacements_np has less capacity than the numnode*numdofpernode");
//      DRT::UTILS::ExtractMyValues(*displacements_np,my_displacements_np,lm);
//
//      // get node id on adjacentelement of immersednode
//      for (int j=0;j<adjacentelement->NumNode();++j)
//      {
//        if (adjacentelement->NodeIds()[j]==immersednode->Id())
//          inode = j;
//      }
//      // update node position X -> x of boundary element
//      {
//        for (int idof=0;idof<DRT::Problem::Instance()->NDim();++idof)
//        {
//          x[idof] = X[idof] + my_displacements_np[inode*immerseddis->NumDof(immersednode)+idof];
//        }
//      }
//      lmowner.clear();
//      lmstride.clear();
//      lm.clear();
//
//      ////////////////////////////////////////////
//      // loop over all background elements
//      ////////////////////////////////////////////
//      for (curr=searchboxgeom.begin(); curr!=searchboxgeom.end(); ++curr)
//      {
//        bool converged = false;
//
//        DRT::Element::DiscretizationType distype = immerseddis->gElement(0)->Shape();
//        switch(distype)
//        {
//        case DRT::Element::hex8 :
//        {
//          MORTAR::UTILS::GlobalToLocal<DRT::Element::hex8>(*(curr->second),&x[0],&xi[0],converged);
//          break;
//        }
//        default:
//        {
//          dserror("DISTYPE NOT SUPPORTED YET. PLEASE CREATE ENTRY IN THIS SWITCH-CASE STATEMENT");
//          break;
//        }
//        }
//
//# ifdef DEBUG
//        if(!converged)
//        {
//          notconvergedcounter ++;
//          //        std::cout<<" Map immersed node with GID "<<immerseddis->gNode(myglobalelements[i])->Id()<<" to element with GID "<<curr->second->Id()<<std::endl;
//          //        dserror("MAPPING FROM IMMERSED NODE TO BACKGROUNDELEMENT DID NOT CONVERGE");
//        }
//# endif
//
//        if ((abs(xi[0])-1.0)<1e-12 and (abs(xi[1])-1.0)<1e-12 and (abs(xi[2])-1.0)<1e-12)
//        {// node i lies in element curr
//          Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidImmersedBase>(curr->second)->SetBoundaryIsImmersed(1);
//          const int* nodes;
//          nodes = curr->second->NodeIds();
//          for(int k=0;k<curr->second->NumNode();++k)
//          {
//            nodeset.insert(nodes[k]);
//          }
//        }
//      }// loop over background elements in searchbox
//    }// loop over immersed nodes
//  }// gpversion = false
//
//  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  // alternative loop
//  // gpversion : checks not for nodes in background elements , but for integration points of boundary elements in background elements
//  //
//  //                           | xi_2
//  //                           |
//  //              _____________|______________
//  //             |             |              |
//  //             |             |              |
//  //             |     1       |         2    |
//  //             |     X       |         X    |
//  //             |             |              |
//  //             |             |              |
//  //      _____________________|______________|_______ xi_1
//  //             |             |              |
//  //             |             |              |
//  //             |     4       |         3    |
//  //             |     X       |         X    |
//  //             |             |              |
//  //             |             |              |
//  //             |_____________|______________|
//  //                           |
//  //
//  else if (gpversion == true)
//  {
//    // key is element id, entries are vectors with current positions of bdry int points
//    std::multimap<int,std::vector<double> > mygeometrygpcurrentposition;
//
//    const int numproc  = comm.NumProc();
//    const int myrank   = comm.MyPID();                     // me
//    const int torank   = (myrank + 1) % numproc;           // to
//    const int fromrank = (myrank + numproc - 1) % numproc; // from
//
//    DRT::Exporter exporter(comm);
//
//    int mygeomsize = immersedgeom.size();
//    int maxsize = -1;
//    comm.MaxAll(&mygeomsize,&maxsize,1);
//    geomcurr=immersedgeom.begin();
//    /////////////////////////////////////////////////////
//    // loop over all immersed boundary elements
//    /////////////////////////////////////////////////////
//    //for (geomcurr=immersedgeom.begin(); geomcurr!=immersedgeom.end(); ++geomcurr)
//    for(int geomcount=0;geomcount<maxsize;++geomcount)
//    {
//      if(geomcurr != immersedgeom.end())
//      {
//        Teuchos::RCP<DRT::Element> immersedelement = geomcurr->second;
//        if(immersedelement == Teuchos::null)
//          dserror("Could not get element with GID %d",immersedelement->Id());
//
//        immersedelement->LocationVector(*immerseddis,lm,lmowner,lmstride);
//        my_displacements_np.resize(lm.size());
//        DRT::UTILS::ExtractMyValues(*displacements_np,my_displacements_np,lm);
//
//        ////////////////////////////////////////////////////////
//        // loop over all gps of immersed element on every proc
//        ////////////////////////////////////////////////////////
//        for(int gp = 0; gp < 4; gp++)
//        {
//          double xibdry[2];
//
//          if(gp==0)       {xibdry[0]=( xigp); xibdry[1]=( xigp);}
//          else if (gp==1) {xibdry[0]=( xigp); xibdry[1]=(-xigp);}
//          else if (gp==2) {xibdry[0]=(-xigp); xibdry[1]=(-xigp);}
//          else if (gp==3) {xibdry[0]=(-xigp); xibdry[1]=( xigp);}
//
//          DRT::Element::DiscretizationType distype = immersedelement->Shape();
//          switch(distype)
//          {
//          case DRT::Element::quad4 :
//          {
//            MORTAR::UTILS::LocalToCurrentGlobal<DRT::Element::quad4>(*immersedelement,3,&xibdry[0],my_displacements_np,&x[0]);
//            std::vector<double> xvec(4);
//            xvec[0]=x[0];
//            xvec[1]=x[1];
//            xvec[2]=x[2];
//            xvec[3]=(double)gp; // gp id
//            mygeometrygpcurrentposition.insert(std::pair<int,std::vector<double> >(geomcurr->second->Id(),xvec));
//            break;
//          }
//          default:
//          {
//            dserror("DISTYPE NOT SUPPORTED YET. PLEASE CREATE ENTRY IN THIS SWITCH-CASE STATEMENT");
//            break;
//          }
//          }
//        }
//      }// if mypid has still an element "geomcount" to send around
//
//      ///////////////////////////////////////////////////////////////////////////////////////////
//      /////////
//      /////////     round robin loop
//      /////////
//      ///////// in the end each rank should store the same mygeometrygpcurrentposition multimap
//      /////////
//      ///////////////////////////////////////////////////////////////////////////////////////////
//      //////////////////////////
//      // preparation
//      //////////////////////////
//      int idtosend = geomcurr->first;
//      std::vector<std::vector<double> >gpstosend;
//
//      std::multimap<int,std::vector<double> >::iterator fit = mygeometrygpcurrentposition.find(idtosend);
//      std::pair <std::multimap<int,std::vector<double> >::iterator, std::multimap<int,std::vector<double> >::iterator> range;
//      range = mygeometrygpcurrentposition.equal_range(fit->first);
//      for (std::multimap<int,std::vector<double> >::iterator it=range.first; it!=range.second; ++it)
//      {
//        std::vector<double> gpstosendcoords;
//        for(int dim=0;dim<3;++dim)
//          gpstosendcoords.push_back(it->second[dim]);
//        gpstosendcoords.push_back(it->second[3]); // gp id
//        gpstosend.push_back(gpstosendcoords);
//      }
////      std::cout<<"PROC "<<myrank<<" idtosend "<<idtosend<<std::endl;
////      for (int gp=0;gp<(int)gpstosend.size();++gp)
////      {
////        std::cout<<"PROC "<<myrank<<" gpstosend : ";
////        for(int dim=0;dim<(int)gpstosend[gp].size();++dim)
////        {
////          std::cout<<gpstosend[gp][dim];
////        }
////        std::cout<<" "<<std::endl;
////      }
//
//      //////////////////////////////////
//      // actual loop
//      /////////////////////////////////
//      for (int irobin = 0; irobin < numproc-1; ++irobin)
//      {
//        std::vector<char> sdata;
//        std::vector<char> rdata;
//        int recievedid;
//        std::vector<std::vector<double> > recievedgps;
//
//        // ---- pack data for sending -----
//        {
//          DRT::PackBuffer data;
//          data.StartPacking();
//          data.AddtoPack(idtosend);
//          {
//            for(int count=0;count<(int)gpstosend.size();++count)
//            {
//              for(int dim=0;dim<3;++dim)
//              {
//                data.AddtoPack(gpstosend[count][dim]);
//              }
//              data.AddtoPack((int)gpstosend[count][3]); // gp id
//            }
//          }
//          std::swap(sdata, data());
//        }
//
//        // ---- send ----
//        MPI_Request request;
//        exporter.ISend(myrank, torank, &(sdata[0]), (int)sdata.size(), 1234, request);
//
//        // ---- receive ----
//        int length = rdata.size();
//        int tag = -1;
//        int from = -1;
//        exporter.ReceiveAny(from,tag,rdata,length);
//        if (tag != 1234 or from != fromrank)
//          dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank, myrank, from, myrank);
//
//        // ---- unpack data -----
//        std::vector<char>::size_type position = 0;
//        recievedid = DRT::ParObject::ExtractInt(position,rdata);
//        while (position < rdata.size())
//        {
//          std::vector<double> gpcoord;
//          gpcoord.push_back(DRT::ParObject::ExtractDouble(position,rdata));
//          gpcoord.push_back(DRT::ParObject::ExtractDouble(position,rdata));
//          gpcoord.push_back(DRT::ParObject::ExtractDouble(position,rdata));
//          gpcoord.push_back((double)DRT::ParObject::ExtractInt(position,rdata));
//          recievedgps.push_back(gpcoord);
//        }
////        std::cout<<"PROC: "<<comm.MyPID()<<"  Recieved Id -> "<<recievedid<<std::endl;
//
//        idtosend = recievedid;
//        gpstosend.swap(recievedgps);
////
////        std::cout<<"PROC "<<myrank<<" idtosend "<<idtosend<<std::endl;
////        for (int gp=0;gp<(int)gpstosend.size();++gp)
////        {
////          std::cout<<"PROC "<<myrank<<" gpstosend : ";
////          for(int dim=0;dim<(int)gpstosend[gp].size();++dim)
////          {
////            std::cout<<gpstosend[gp][dim];
////          }
////          std::cout<<" "<<std::endl;
////        }
//
//        // wait for all communication to finish
//        exporter.Wait(request);
//        comm.Barrier();
//
//        if(!immersedgeom.count(recievedid) and !mygeometrygpcurrentposition.count(recievedid))
//        {
//          mygeometrygpcurrentposition.insert(std::pair<int,std::vector<double> >(idtosend,gpstosend[0]));
//          mygeometrygpcurrentposition.insert(std::pair<int,std::vector<double> >(idtosend,gpstosend[1]));
//          mygeometrygpcurrentposition.insert(std::pair<int,std::vector<double> >(idtosend,gpstosend[2]));
//          mygeometrygpcurrentposition.insert(std::pair<int,std::vector<double> >(idtosend,gpstosend[3]));
//        }
//
//        } // end for irobin
//
//      if(geomcurr != immersedgeom.end())
//      {
//        geomcurr++;
//      }
//    } // end loop over all immersed boundary elements
//
////    ////////////////////////////////////////////////////////
////    // OUTPUT geometrycurentpositions to all procs
////    ////////////////////////////////////////////////////////
////    if(comm.MyPID()==-1)
////    {
////      std::multimap<int,std::vector<double> >::iterator fit = mygeometrygpcurrentposition.begin(); // iterator to my first element of local multimap
////      int mysize = mygeometrygpcurrentposition.size();
////      for (int i=0; i<mysize; ++i)
////      {
////        std::pair <std::multimap<int,std::vector<double> >::iterator, std::multimap<int,std::vector<double> >::iterator> range;
////        range = mygeometrygpcurrentposition.equal_range(fit->first);
////
////        std::cout<<"PROC: "<<comm.MyPID()<<" ID: "<<fit->first<<" IntPoints: ";
////        for (std::multimap<int,std::vector<double> >::iterator it=range.first; it!=range.second; ++it)
////        {
////          //[coord_1, coord_2, coord_3, gp_id]
////          std::cout<<" "<<std::setprecision(6)<<it->second[0]<<" "<<std::setprecision(6)<<it->second[1]<<" "<<std::setprecision(6)<<it->second[2]<<" ID "<<(int)it->second[3]<<"|";
////        }
////        std::cout<<""<<std::endl;
////        fit++;
////      }
////    }
////    comm.Barrier();
////    dserror("");
//
//# ifdef DEBUG
//    std::cout<<"PROC "<<comm.MyPID()<<" : size of mygeometrygpcurrentposition : "<<mygeometrygpcurrentposition.size()<<std::endl;
//# endif
//
//    //////////////////////////////////////////////////////////////////////
//    // loop over all background elements
//    //
//    // determine background elements in which there are int points of
//    // the boundary of the immersed dis
//    //
//    /////////////////////////////////////////////////////////////////////
//    for (curr=searchboxgeom.begin(); curr!=searchboxgeom.end(); ++curr)
//    {
//      int immersedeleid=-1;
//      int recentimmersedeleid=-1;
//
//      for(std::multimap<int,std::vector<double> >::iterator fit=mygeometrygpcurrentposition.begin();fit!=mygeometrygpcurrentposition.end();++fit)
//      {
//        immersedeleid = fit->first;
//
//        // immerse ele id already stored in gpmap
//        if(immersedeleid == recentimmersedeleid)
////=======
////      if ((abs(xi[0])-1.0)<1e-12 and (abs(xi[1])-1.0)<1e-12 and (abs(xi[2])-1.0)<1e-12)
////      {// node i lies in element curr
////        Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidImmersedBase>(curr->second)->SetBoundaryIsImmersed(1);
////        const int* nodes;
////        nodes = curr->second->NodeIds();
////        for(int k=0;k<curr->second->NumNode();++k)
////>>>>>>> .r19849
//        {
//          // do nothing because background ele already in map and immersedele already in map (makes no sense) multimap issue wich is to be surveyed in the future
//        }
//        else
//        {
//          double gpid = -111.;
//          std::pair <std::multimap<int,std::vector<double> >::iterator, std::multimap<int,std::vector<double> >::iterator> range;
//          range = mygeometrygpcurrentposition.equal_range(fit->first);
//          for (std::multimap<int,std::vector<double> >::iterator it=range.first; it!=range.second; ++it)
//          {
//            bool converged = false;
//
//            x[0]=it->second[0];
//            x[1]=it->second[1];
//            x[2]=it->second[2];
//
//            DRT::Element::DiscretizationType distype = immerseddis->gElement(0)->Shape();
//            switch(distype)
//            {
//            case DRT::Element::hex8 :
//            {
//              MORTAR::UTILS::GlobalToLocal<DRT::Element::hex8>(*(curr->second),&x[0],&xi[0],converged);
//              break;
//            }
//            default:
//            {
//              dserror("DISTYPE NOT SUPPORTED YET. PLEASE CREATE ENTRY IN THIS SWITCH-CASE STATEMENT");
//              break;
//            }
//            }
//
//# ifdef DEBUG
//            if(!converged)
//            {
//              notconvergedcounter ++;
//              //        std::cout<<" Map immersed node with GID "<<immerseddis->gNode(myglobalelements[i])->Id()<<" to element with GID "<<curr->second->Id()<<std::endl;
//              //        dserror("MAPPING FROM IMMERSED NODE TO BACKGROUNDELEMENT DID NOT CONVERGE");
//            }
//# endif
//
//            if ((abs(xi[0])-1.0)<1e-12 and (abs(xi[1])-1.0)<1e-12 and (abs(xi[2])-1.0)<1e-12)
//            {// gp lies in element curr
//              if (Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidImmersedBase>(curr->second)->IsBoundaryImmersed()==0)
//              {
//                Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidImmersedBase>(curr->second)->SetBoundaryIsImmersed(1);
//                const int* nodes;
//                nodes = curr->second->NodeIds();
//                for(int k=0;k<curr->second->NumNode();++k)
//                {
//                  nodeset.insert(nodes[k]);
//                }
//              }
//
//              // here the gpmap is filled. it is needed outside this function to interpolate
//              // the background quantities to the gauss point coordinates in the background
//              // parameter space.
//              if (gpmap != Teuchos::null)
//              {
//                if(gpmap->find(curr->second->Id())==gpmap->end()) // id not yet in map
//                {
//                  Teuchos::RCP<std::vector<double> > tempvec = Teuchos::rcp(new std::vector<double> );
//                  Teuchos::RCP<std::map<int, Teuchos::RCP<std::vector<double> > > > tempgpmap = Teuchos::rcp(new std::map<int, Teuchos::RCP<std::vector<double> > >);
//                  tempgpmap->insert(std::pair<int,Teuchos::RCP<std::vector<double> > >(immersedeleid,tempvec));
//                  gpmap->insert(std::pair<int,Teuchos::RCP<std::map<int,Teuchos::RCP<std::vector<double> > > > >(curr->second->Id(),tempgpmap));
//                }
//                else if(gpmap->at(curr->second->Id())->find(immersedeleid)==gpmap->at(curr->second->Id())->end())
//                {
//                  Teuchos::RCP<std::vector<double> > tempvec = Teuchos::rcp(new std::vector<double>);
//                  gpmap->at(curr->second->Id())->insert(std::pair<int,Teuchos::RCP<std::vector<double> > >(immersedeleid,tempvec));;
//                }
//
//                gpid=it->second[3];
//                gpmap->at(curr->second->Id())->at(immersedeleid)->push_back(xi[0]);
//                gpmap->at(curr->second->Id())->at(immersedeleid)->push_back(xi[1]);
//                gpmap->at(curr->second->Id())->at(immersedeleid)->push_back(xi[2]);
//                gpmap->at(curr->second->Id())->at(immersedeleid)->push_back(gpid);
//
//                recentimmersedeleid = immersedeleid;
//              }
//            }// if not yet in map
//          }// lies in element?
//        }// loop over range of key
//      }// loop over all map keys
//    }// loop over background elements in searchbox
//  } // gpversion = true
//
//  std::vector<int> nodevector(nodeset.size()); // variable to return
//  std::copy(nodeset.begin(), nodeset.end(), nodevector.begin());
//
//
//# ifdef DEBUG
//    std::cout<<"PROC "<<comm.MyPID()<<" "<<notconvergedcounter<<" mappings did not converge"<<std::endl;
//    std::cout<<"PROC "<<comm.MyPID()<<" : searchboxgeom.size() = "<<searchboxgeom.size()<<std::endl;
//    std::cout<<"PROC "<<comm.MyPID()<<" : immersedgeom.size() = "<<immersedgeom.size()<<std::endl;
//    std::cout<<"PROC "<<comm.MyPID()<<" : identified "<<nodeset.size()<<" nodes in immersion domain."<<std::endl;
//    if(nodeset.size() != nodevector.size())
//      dserror("nodeset and nodevector must have same size");
//    if (gpmap != Teuchos::null)
//    {
//      std::cout<<"PROC "<<comm.MyPID()<<" : size of gpmap = "<<gpmap->size()<<std::endl;
////      /////////////////////////////////////////////////
////      ////
////      ////  display gpmap
////      ////
////      /////////////////////////////////////////////////
////      int i = 0;
////      int gpcounter = 0;
////      int immersedcounter = 0;
////      std::map<int,Teuchos::RCP<std::vector<double> > >::iterator curr;
////      std::map<int,Teuchos::RCP<std::map<int,Teuchos::RCP<std::vector<double> > > > >::iterator backgrdcurr;
////      for(backgrdcurr=gpmap->begin();backgrdcurr!=gpmap->end();++backgrdcurr)
////      { ++i;
////        std::cout<<"PROC "<<comm.MyPID()<<": gpmap->at("<<backgrdcurr->first<<") (entry " <<i<<") :"<<std::endl;
////        for(curr=backgrdcurr->second->begin();curr!=backgrdcurr->second->end();++curr)
////        {
////          ++immersedcounter;
////          std::cout<<"    immersedeleid: "<<curr->first<<" with gp coords : [";
////          for(int k=0;k<(int)curr->second->size();++k)
////          {
////            ++gpcounter;
////            if(k>0 and k%4==3)
////              std::cout<<" gpID: "<<(int)curr->second->at(k)<<"] [";
////            else
////              std::cout<<curr->second->at(k)<<" ";
////          }
////          std::cout<<"\n "<<std::endl;
////        }
////      }
////      std::cout<<"DISPLAYED "<<i<<" BACKGROUND ELEMENT IDs"<<std::endl;
////      std::cout<<"DISPLAYED "<<gpcounter/4<<" GAUSS POINT COORDS"<<std::endl;
////      std::cout<<"DISPLAYED "<<immersedcounter<<" IMMERSED ELEMENTS"<<std::endl;
//
//    }
//
//    comm.Barrier();
//
//#endif
//
//  return nodevector;
//}


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
//void IMMERSED::ImmersedBase::UpdateVolumeCondition(Teuchos::RCP<DRT::Discretization> dis, std::vector<int> dvol_fenode, DRT::Condition::ConditionType condtype, std::string condname)
//{
//
//}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::EvaluateImmersed(Teuchos::ParameterList& params,
                                              Teuchos::RCP<DRT::Discretization> dis,
                                              DRT::AssembleStrategy* strategy,
                                              std::map<int,std::set<int> >* elementstoeval,
                                              Teuchos::RCP<GEO::SearchTree> structsearchtree,
                                              std::map<int,LINALG::Matrix<3,1> >* currpositions_struct,
                                              int action,
                                              bool evaluateonlyboundary)
{
  // pointer to element
  DRT::Element* ele;

  for(std::map<int, std::set<int> >::const_iterator closele = elementstoeval->begin(); closele != elementstoeval->end(); closele++)
  {
    for(std::set<int>::const_iterator eleIter = (closele->second).begin(); eleIter != (closele->second).end(); eleIter++)
    {
      ele=dis->gElement(*eleIter);

      DRT::ELEMENTS::FluidImmersedBase* immersedelebase = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(ele);
      if(immersedelebase==NULL)
        dserror("dynamic cast from DRT::Element* to DRT::ELEMENTS::FluidImmersedBase* failed");

      // evaluate this element and fill vector with immersed dirichlets
      int row = strategy->FirstDofSet();
      int col = strategy->SecondDofSet();

      params.set<int>("action",action);
      params.set<Teuchos::RCP<GEO::SearchTree> >("structsearchtree_rcp",structsearchtree);
      params.set<std::map<int,LINALG::Matrix<3,1> >* >("currpositions_struct",currpositions_struct);
      params.set<int>("Physical Type",INPAR::FLUID::poro_p1);
      if(dis->Name()=="fluid")
        params.set<std::string>("immerseddisname","structure");
      else if (dis->Name()=="porofluid")
        params.set<std::string>("immerseddisname","cell");
      else
        dserror("no corresponding immerseddisname set for this type of backgrounddis!");

      DRT::Element::LocationArray la(1);
      immersedelebase->LocationVector(*dis,la,false);
      strategy->ClearElementStorage( la[row].Size(), la[col].Size() );

      if(!evaluateonlyboundary)
        immersedelebase->Evaluate(params,*dis,la[0].lm_,
            strategy->Elematrix1(),
            strategy->Elematrix2(),
            strategy->Elevector1(),
            strategy->Elevector2(),
            strategy->Elevector3());
      else
      {
        if(immersedelebase->IsBoundaryImmersed())
          immersedelebase->Evaluate(params,*dis,la[0].lm_,
              strategy->Elematrix1(),
              strategy->Elematrix2(),
              strategy->Elevector1(),
              strategy->Elevector2(),
              strategy->Elevector3());
      }

      strategy->AssembleVector1( la[row].lm_, la[row].lmowner_ );
    }
  }
  return;
} // EvaluateImmersed

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::EvaluateImmersedNoAssembly(Teuchos::ParameterList& params,
                                              Teuchos::RCP<DRT::Discretization> dis,
                                              std::map<int,std::set<int> >* elementstoeval,
                                              Teuchos::RCP<GEO::SearchTree> structsearchtree,
                                              std::map<int,LINALG::Matrix<3,1> >* currpositions_struct,
                                              int action
                                                        )
{
  // pointer to element
  DRT::Element* ele;

  for(std::map<int, std::set<int> >::const_iterator closele = elementstoeval->begin(); closele != elementstoeval->end(); closele++)
  {
    for(std::set<int>::const_iterator eleIter = (closele->second).begin(); eleIter != (closele->second).end(); eleIter++)
    {
      ele=dis->gElement(*eleIter);

      DRT::ELEMENTS::FluidImmersedBase* immersedelebase = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(ele);
      if(immersedelebase==NULL)
        dserror("dynamic cast from DRT::Element* to DRT::ELEMENTS::FluidImmersedBase* failed");

      // provide important objects to ParameterList
      params.set<int>("action",action);
      params.set<Teuchos::RCP<GEO::SearchTree> >("structsearchtree_rcp",structsearchtree);
      params.set<std::map<int,LINALG::Matrix<3,1> >* >("currpositions_struct",currpositions_struct);
      params.set<int>("Physical Type",INPAR::FLUID::poro_p1);
      if(dis->Name()=="fluid")
        params.set<std::string>("immerseddisname","structure");
      else if (dis->Name()=="porofluid")
        params.set<std::string>("immerseddisname","cell");
      else
        dserror("no corresponding immerseddisname set for this type of backgrounddis!");

      // evaluate the element
      Epetra_SerialDenseMatrix dummymat;
      Epetra_SerialDenseVector dummyvec;

      DRT::Element::LocationArray la(1);
      immersedelebase->LocationVector(*dis,la,false);

      immersedelebase->Evaluate(params,*dis,la[0].lm_,dummymat,dummymat,dummyvec,dummyvec,dummyvec);
    }
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::EvaluateScaTraWithInternalCommunication(Teuchos::RCP<DRT::Discretization> dis,
                                                                     const Teuchos::RCP<const DRT::Discretization> idis,
                                                                     DRT::AssembleStrategy* strategy,
                                                                     std::map<int,std::set<int> >* elementstoeval,
                                                                     Teuchos::RCP<GEO::SearchTree> structsearchtree,
                                                                     std::map<int,LINALG::Matrix<3,1> >* currpositions_struct,
                                                                     Teuchos::ParameterList& params,
                                                                     bool evaluateonlyboundary)
{
  // pointer to element
  DRT::Element* ele;
  DRT::Element* iele;

  for(std::map<int, std::set<int> >::const_iterator closele = elementstoeval->begin(); closele != elementstoeval->end(); closele++)
  {
    for(std::set<int>::const_iterator eleIter = (closele->second).begin(); eleIter != (closele->second).end(); eleIter++)
    {
      ele=dis->gElement(*eleIter);
      iele=idis->gElement(*eleIter);

      DRT::ELEMENTS::FluidImmersedBase* immersedelebase = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(iele);
      if(immersedelebase==NULL)
        dserror("dynamic cast from DRT::Element* to DRT::ELEMENTS::FluidImmersedBase* failed");

      // evaluate this element and fill vector with immersed dirichlets
      int row = strategy->FirstDofSet();
      int col = strategy->SecondDofSet();

      params.set<Teuchos::RCP<GEO::SearchTree> >("structsearchtree_rcp",structsearchtree);
      params.set<std::map<int,LINALG::Matrix<3,1> >* >("currpositions_struct",currpositions_struct);
      params.set<int>("Physical Type",INPAR::FLUID::poro_p1);

      DRT::Element::LocationArray la(dis->NumDofSets());
      ele->LocationVector(*dis,la,false);
      strategy->ClearElementStorage( la[row].Size(), la[col].Size() );

      if(!evaluateonlyboundary)
        ele->Evaluate(params,*dis,la,
            strategy->Elematrix1(),
            strategy->Elematrix2(),
            strategy->Elevector1(),
            strategy->Elevector2(),
            strategy->Elevector3());
      else
      {
        if(immersedelebase->IsBoundaryImmersed())
          ele->Evaluate(params,*dis,la,
              strategy->Elematrix1(),
              strategy->Elematrix2(),
              strategy->Elevector1(),
              strategy->Elevector2(),
              strategy->Elevector3());
      }

      strategy->AssembleVector1( la[row].lm_, la[row].lmowner_ );
    }
  }
} // EvaluateWithInternalCommunication

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/// Reduces to standard EvaluateCondition on one proc.
/// Evaluate a specific condition using assemble strategy allowing communication at element level
/// until every conditioned element is evaluated. Needed especially during interpolation from an
/// other discretization to the conditioned elements (e.g. in immersed method).
/// The integration point of a conditioned element requesting a quantity may be owned by another
/// proc as the interpolating element providing this quantity.  rauch 05/14
void IMMERSED::ImmersedBase::EvaluateInterpolationCondition
(
    Teuchos::RCP<DRT::Discretization> evaldis,
    Teuchos::ParameterList& params,
    DRT::AssembleStrategy & strategy,
    const std::string& condstring,
    const int condid
)
{
# ifdef DEBUG
  if (!(evaldis->Filled()) ) dserror("FillComplete() was not called");
  if (!(evaldis->HaveDofs())) dserror("AssignDegreesOfFreedom() was not called");
# endif

  int row = strategy.FirstDofSet();
  int col = strategy.SecondDofSet();

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  params.set<int>("dummy_call",0);

  DRT::Element::LocationArray la(evaldis->NumDofSets());

  std::multimap<std::string,Teuchos::RCP<DRT::Condition> >::iterator fool;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (fool=evaldis->GetAllConditions().begin(); fool!=evaldis->GetAllConditions().end(); ++fool)
  {
    if (fool->first == condstring)
    {
      DRT::Condition& cond = *(fool->second);
      if (condid == -1 || condid ==cond.GetInt("ConditionID"))
      {
        std::map<int,Teuchos::RCP<DRT::Element> >& geom = cond.Geometry();
        if (geom.empty()) dserror("evaluation of condition with empty geometry on proc %d",evaldis->Comm().MyPID());

        std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;

        // Evaluate Loadcurve if defined. Put current load factor in parameterlist
        const std::vector<int>* curve  = cond.Get<std::vector<int> >("curve");
        int curvenum = -1;
        if (curve) curvenum = (*curve)[0];
        double curvefac = 1.0;
        if (curvenum>=0 && usetime)
          curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

        // Get ConditionID of current condition if defined and write value in parameterlist
        const std::vector<int>*    CondIDVec  = cond.Get<std::vector<int> >("ConditionID");
        if (CondIDVec)
        {
          params.set("ConditionID",(*CondIDVec)[0]);
          char factorname[30];
          sprintf(factorname,"LoadCurveFactor %d",(*CondIDVec)[0]);
          params.set(factorname,curvefac);
        }
        else
        {
          params.set("LoadCurveFactor",curvefac);
        }
        params.set<Teuchos::RCP<DRT::Condition> >("condition", fool->second);

        int mygeometrysize=-1234;
        if(geom.empty()==true)
          mygeometrysize=0;
        else
          mygeometrysize=geom.size();
        int maxgeometrysize=-1234;
        evaldis->Comm().MaxAll(&mygeometrysize,&maxgeometrysize,1);
        curr=geom.begin();

#ifdef DEBUG
        std::cout<<"PROC "<<evaldis->Comm().MyPID()<<": mygeometrysize = "<<mygeometrysize<<" maxgeometrysize = "<<maxgeometrysize<<std::endl;
#endif


        // enter loop on every proc until the last proc evaluated his last geometry element
        // because there is communication happening inside
        for (int i=0;i<maxgeometrysize;++i)
        {
          if(i>=mygeometrysize)
            params.set<int>("dummy_call",1);

          // get element location vector and ownerships
          // the LocationVector method will return the the location vector
          // of the dofs this condition is meant to assemble into.
          // These dofs do not need to be the same as the dofs of the element
          // (this is the standard case, though). Special boundary conditions,
          // like weak dirichlet conditions, assemble into the dofs of the parent element.
          curr->second->LocationVector(*evaldis,la,false,condstring,params);

          // get dimension of element matrices and vectors
          // Reshape element matrices and vectors and init to zero

          strategy.ClearElementStorage( la[row].Size(), la[col].Size() );

          // call the element specific evaluate method
          int err = curr->second->Evaluate(params,*evaldis,la,
              strategy.Elematrix1(),
              strategy.Elematrix2(),
              strategy.Elevector1(),
              strategy.Elevector2(),
              strategy.Elevector3());
          if (err) dserror("error while evaluating elements");

          // assemble every element contribution only once
          // do not assemble after dummy call for internal communication
          if(i<mygeometrysize)
          {
            // assembly
            int eid = curr->second->Id();
            strategy.AssembleMatrix1( eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_ );
            strategy.AssembleMatrix2( eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_ );
            strategy.AssembleVector1( la[row].lm_, la[row].lmowner_ );
            strategy.AssembleVector2( la[row].lm_, la[row].lmowner_ );
            strategy.AssembleVector3( la[row].lm_, la[row].lmowner_ );
          }

          // go to next element
          if (i<(mygeometrysize-1))
            ++curr;

        } // for 0 to max. geometrysize over all procs
      } // if check of condid successful
    } // if condstring found
  } //for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::SearchPotentiallyCoveredBackgrdElements(
        std::map<int,std::set<int> >* current_subset_tofill,
        Teuchos::RCP<GEO::SearchTree> backgrd_SearchTree,
        const DRT::Discretization& dis,
        const std::map<int, LINALG::Matrix<3, 1> >& currentpositions,
        const LINALG::Matrix<3, 1>& point,
        const double radius,
        const int label)
{

  *current_subset_tofill = backgrd_SearchTree->searchElementsInRadius(dis,currentpositions,point,radius,label);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::EvaluateSubsetElements(Teuchos::ParameterList& params,
                                                    Teuchos::RCP<DRT::Discretization> dis,
                                                    std::map<int,std::set<int> >& elementstoeval,
                                                    int action)
{
  // pointer to element
  DRT::Element* ele;

  // initialize location array
  DRT::Element::LocationArray la(1);

  for(std::map<int, std::set<int> >::const_iterator closele = elementstoeval.begin(); closele != elementstoeval.end(); closele++)
  {
    for(std::set<int>::const_iterator eleIter = (closele->second).begin(); eleIter != (closele->second).end(); eleIter++)
    {
        ele = dis->gElement(*eleIter);

        Epetra_SerialDenseMatrix dummymatrix;
        Epetra_SerialDenseVector dummyvector;
        ele->Evaluate(params,*dis,la,dummymatrix,dummymatrix,dummyvector,dummyvector,dummyvector);

    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::CreateGhosting(const Teuchos::RCP<DRT::Discretization> distobeghosted)
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
  int nummycolelements = newelecolmap->NumMyElements();
  int sizelist[distobeghosted->Comm().NumProc()];
  distobeghosted->Comm().GatherAll(&nummycolnodes,&sizelist[0],1);
  std::cout<<"PROC "<<distobeghosted->Comm().MyPID()<<" : "<<nummycolnodes<<" colnodes"<<std::endl;
  distobeghosted->Comm().Barrier();
  std::cout<<"PROC "<<distobeghosted->Comm().MyPID()<<" : "<<nummycolelements<<" colelements"<<std::endl;
  distobeghosted->Comm().Barrier();
  std::cout<<"PROC "<<distobeghosted->Comm().MyPID()<<" : "<<distobeghosted->lColElement(0)->Nodes()[0]->Id()<<" first ID of first node of first colele"<<std::endl;
  distobeghosted->Comm().Barrier(); // wait for procs
  for(int k=1;k<distobeghosted->Comm().NumProc();++k)
  {
    if(sizelist[k-1]!=nummycolnodes)
      dserror("Since whole dis is ghosted every processor should have the same number of colnodes. This is not the case! Fix this!");
  }
#endif

}

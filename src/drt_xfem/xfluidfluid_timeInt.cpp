/*!
\file xfluidfluid_timeInt.classes

\brief provides the xfluidfluid timeIntegration

<pre>
Maintainer: Shadan Shahmiri
            shahmiri@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include "xfluidfluid_timeInt.H"
#include "xfem_fluidwizard.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_cut/cut_boundingbox.H"
#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_position.H"
#include "../drt_cut/cut_point.H"
#include "../drt_cut/cut_element.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_interface.H"
#include "../drt_fluid_ele/fluid_ele_factory.H"

#include <iostream>

XFEM::XFluidFluidTimeIntegration::XFluidFluidTimeIntegration(
  const RCP<DRT::Discretization> bgdis,
  const RCP<DRT::Discretization> embdis,
  RCP<XFEM::FluidWizard>         wizard,
  int                            step,
  enum INPAR::XFEM::XFluidFluidTimeInt xfem_timeintapproach,
  const ParameterList&              params
  ) :
  embdis_(embdis),
  step_(step)
  {
    myrank_ = bgdis->Comm().MyPID();
    numproc_ = bgdis->Comm().NumProc();
    CreateBgNodeMaps(bgdis,wizard);
    currentbgdofmap_ = bgdis->DofRowMap();
    timeintapproach_ = xfem_timeintapproach;
    params_ = params;

    return;

  } // end constructor

// -------------------------------------------------------------------
// map of standard node ids and their dof-gids in for this time step
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::CreateBgNodeMaps(const RCP<DRT::Discretization> bgdis,
                                                        RCP<XFEM::FluidWizard>         wizard)
{
  const Epetra_Map* noderowmap = bgdis->NodeRowMap();

  // map of standard nodes and their dof-ids
  for (int lid=0; lid<noderowmap->NumMyPoints(); lid++)
  {
    int gid;
    // get global id of a node
    gid = noderowmap->GID(lid);
    // get the node
    DRT::Node * node = bgdis->gNode(gid);
    GEO::CUT::Node * n = wizard->GetNode(node->Id());
    if (n!=NULL) // xfem nodes
    {
      // set of volumecells which belong to the node. The size of the
      // vector depends on the combined volumecells which are around
      // the node. For the std nodes the size is 1.
      std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> > vcs= n->DofCellSets();
      std::vector<int> parentelements;

      for (size_t i=0; i<vcs.size(); i++ )
      {
        // we just need the first vc to find out the parent element.
        // for hex20 we have sets of volumecells otherwise for linear
        // elements the size of the set is one

        // get the i volumecell set, which could be a set of volumecells
        std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> myvcsets = vcs.at(i);

        //get the first set
        std::set<GEO::CUT::plain_volumecell_set>::iterator s = myvcsets.begin();

        // get the first volume cell (vc)
        GEO::CUT::plain_volumecell_set k = *s;
        GEO::CUT::plain_volumecell_set::iterator it = k.begin();

        GEO::CUT::VolumeCell * vc = *it;

        // get the element handle and the nds vector
        std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
        std::vector< std::vector< int > >   nds_sets;
        DRT::Element * actele = bgdis->gElement(vc->ParentElement()->Id());
        parentelements.push_back(actele->Id());
        GEO::CUT::ElementHandle * e = wizard->GetElement( actele );
        e->GetVolumeCellsDofSets( cell_sets, nds_sets );

        parenteletondsset_[actele->Id()] = nds_sets;

      }

      nodetoparentele_[gid] = parentelements;


      GEO::CUT::Point * p = n->point();
      GEO::CUT::Point::PointPosition pos = p->Position();
      if (pos==GEO::CUT::Point::outside and bgdis->NumDof(node) != 0) //std
      {
        //cout << node->Id() << "std!! size()" << vcs.size() << endl;
        //cout << " outside " << pos <<  " "<< node->Id() << endl;
        vector<int> gdofs = bgdis->Dof(node);
        stdnodenp_[gid] = gdofs;
      }
      else if (pos==GEO::CUT::Point::inside and  bgdis->NumDof(node) == 0) //void
      {
        //cout << " inside " <<  pos << " " << node->Id() << endl;
      }
      else if (pos==GEO::CUT::Point::inside and  bgdis->NumDof(node) != 0) //enriched
      {
        //cout << " inside enriched" <<  pos << " " << node->Id() << endl;
        vector<int> gdofs = bgdis->Dof(node);
        enrichednodenp_[gid] = gdofs;
      }
      else if (pos==GEO::CUT::Point::oncutsurface and  bgdis->NumDof(node) != 0)
      {
        vector<int> gdofs = bgdis->Dof(node);
        stdnodenp_[gid] = gdofs;
      }
      // this case only happens if the two fluid domains are both the same
      // size (void)
      else if (pos==GEO::CUT::Point::oncutsurface and  bgdis->NumDof(node) == 0)
      {
        //cout << " on the surface " <<  pos << " " << node->Id() << endl;
      }
      else
      {
        cout << "  hier ?! " <<  pos << " " <<  node->Id() <<  " numdof: " << bgdis->NumDof(node) << endl;
      }
    }
    else if( bgdis->NumDof(node) != 0) // no xfem node
    {
      vector<int> gdofs = bgdis->Dof(node);
      stdnodenp_[gid] = gdofs;
    }
    else
      cout << " why here? " << "node "<<  node->Id()<< " " <<  bgdis->NumDof(node) <<  endl;
  }

#ifdef PARALLEL
  // build a reduced map of all noderowmaps of all processors
  RCP<Epetra_Map> allnoderowmap = LINALG::AllreduceEMap(*noderowmap);
  // gather the informations of all processors
  DRT::Exporter ex(*noderowmap,*allnoderowmap,bgdis->Comm());
  ex.Export(stdnodenp_);
  ex.Export(enrichednodenp_);
#endif

//  debug output
//   for(map<int, vector<int> >::iterator iter = nodetoparentele_.begin(); iter!= nodetoparentele_.end();
//       iter++)
//   {
//     cout << "ngid " <<  iter->first << endl;
//     vector<int> fff= iter->second;
//     for (int u=0;u<fff.size();++u)
//       cout << fff.at(u);
//     cout << endl;
//   }
//   for(std::map<int, std::vector< std::vector< int > > >::iterator iter = parenteletondsset_.begin(); iter!= parenteletondsset_.end();
//       iter++)
//   {
//     cout << "elegid " << iter->first << endl;
//     cout << "size nds " << iter->second.size() << endl;
//     for (int s=0; s < iter->second.size(); ++s)
//     {
//       const std::vector<int> & nds = iter->second[s];
//       for (int ttt=0; ttt<nds.size(); ++ttt)
//         cout << nds.at(ttt) << " " ;
//       cout << endl;
//     }
//   }
//     for (int i=0; i<bgdis->NumMyColNodes(); ++i)
//     {
//        const DRT::Node* actnode = bgdis->lColNode(i);
//        map<int, vector<int> >::const_iterator iter = stdnoden_.find(actnode->Id());
//        map<int, vector<int> >::const_iterator iter2 = enrichednoden_.find(actnode->Id());
//        map<int, vector<int> >::const_iterator iter3 = stdnodenp_.find(actnode->Id());
//        map<int, vector<int> >::const_iterator iter4 = enrichednodenp_.find(actnode->Id());
//        if (iter2 != enrichednoden_.end()) cout  << " enrichned n : " << actnode->Id() << " "  ;
//        if (iter2 == enrichednoden_.end() and iter == stdnoden_.end()) cout  << " void n :" <<  actnode->Id() << " "  ;
//        if (iter4 != enrichednodenp_.end()) cout  << " enrichned np : " << actnode->Id() << " "  ;
//        if (iter4 == enrichednodenp_.end() and iter3 == stdnodenp_.end()) cout  << " void np :" <<  actnode->Id() << " "  ;
//        if (iter  != stdnoden_.end()) cout  << " std n : " << actnode->Id() << " "  ;
//        if (iter3 != stdnodenp_.end()) cout  << " std np :" <<  actnode->Id() << " "  ;
//    }
}

// -------------------------------------------------------------------
// save the old maps
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SaveBgNodeMaps()
{
  // save the old maps and clear the maps for the new cut
  // (all maps are related to the background fluid)
  stdnoden_ = stdnodenp_;
  enrichednoden_ = enrichednodenp_;
  stdnodenp_.clear();
  enrichednodenp_.clear();
}

// -------------------------------------------------------------------
// - Save the old maps of bg nodes
// - Create new map of bg nodes
// - Gmsh-output
// -------------------------------------------------------------------
int XFEM::XFluidFluidTimeIntegration::SaveAndCreateNewBgNodeMaps(RCP<DRT::Discretization> bgdis,
                                                                 RCP<XFEM::FluidWizard>   wizard)
{

  // save the old maps and clear the maps for the new cut
  // (all maps are related to the background fluid)
  SaveBgNodeMaps();

  // Create new maps
  CreateBgNodeMaps(bgdis,wizard);

  oldbgdofmap_ = currentbgdofmap_;
  currentbgdofmap_ =  bgdis->DofRowMap();

  if ((oldbgdofmap_->SameAs(*currentbgdofmap_) and (stdnoden_ == stdnodenp_) and
       (enrichednoden_ == enrichednodenp_)))
    samemaps_ = true;
  else samemaps_ = false;

  GmshOutput(bgdis);

  return samemaps_;
}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void  XFEM::XFluidFluidTimeIntegration::CreateBgNodeMapsForRestart(RCP<DRT::Discretization> bgdis,
                                                                   RCP<XFEM::FluidWizard>   wizard)
{

  // Create new maps
  CreateBgNodeMaps(bgdis,wizard);

  currentbgdofmap_ =  bgdis->DofRowMap();

}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SetNewBgStatevectorAndProjectEmbToBg(const RCP<DRT::Discretization>        bgdis,
                                                                            Teuchos::RCP<Epetra_Vector>           bgstatevn,
                                                                            Teuchos::RCP<Epetra_Vector>           bgstatevnp,
                                                                            Teuchos::RCP<Epetra_Vector>           embstatevn,
                                                                            Teuchos::RCP<Epetra_Vector>           aledispn)
{
  if (timeintapproach_ == INPAR::XFEM::Xff_TimeInt_FullProj or timeintapproach_ == INPAR::XFEM::Xff_TimeInt_ProjIfMoved)
    SetNewBgStatevectorFullProjection(bgdis, bgstatevn, bgstatevnp, embstatevn, aledispn);
  else if(timeintapproach_ == INPAR::XFEM::Xff_TimeInt_KeepGhostValues)
    SetNewBgStatevectorKeepGhostValues(bgdis, bgstatevn, bgstatevnp, embstatevn, aledispn);
  else
    dserror("xfem time integration approach unknown!");
}

// -------------------------------------------------------------------
// Always do the projection from embedded fluid. Also for the enriched
// nodes.
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SetNewBgStatevectorFullProjection(const RCP<DRT::Discretization>        bgdis,
                                                                         Teuchos::RCP<Epetra_Vector>           bgstatevn,
                                                                         Teuchos::RCP<Epetra_Vector>           bgstatevnp,
                                                                         Teuchos::RCP<Epetra_Vector>           embstatevn,
                                                                         Teuchos::RCP<Epetra_Vector>           aledispn)
{

  // coordinates of bg-nodes which need the projection from embedded dis
  std::vector<LINALG::Matrix<3,1> > bgnodes_coords;
  // the vector containing interpolated values from embedded dis
  std::vector<LINALG::Matrix<4,1> > interpolated_vecs;
  // bg-node ids which have no history and need a projection
  std::vector<int> bgnodeidwithnohistory;

  projectednodeids_.clear();

  // loop over bg-row-nodes of each processor
  for (int lnid=0; lnid<bgdis->NumMyRowNodes(); lnid++)
  {
    DRT::Node* bgnode = bgdis->lRowNode(lnid);
    map<int, vector<int> >::const_iterator iterstn = stdnoden_.find(bgnode->Id());
    map<int, vector<int> >::const_iterator iterstnp = stdnodenp_.find(bgnode->Id());
    map<int, vector<int> >::const_iterator iteren = enrichednoden_.find(bgnode->Id());
    map<int, vector<int> >::const_iterator iterenp = enrichednodenp_.find(bgnode->Id());

    // Transfer the dofs:
    // n:std -> n+1:std, n:std -> n+1:enriched
    if ((iterstn != stdnoden_.end() and iterstnp != stdnodenp_.end()) or
        (iterstn != stdnoden_.end() and iterenp != enrichednodenp_.end()))
    {
      int numsets = bgdis->NumDof(bgnode)/4;
      vector<int> gdofsn = iterstn->second;

      //TODO!! die richtige dofs von bgstatevn rauspicke, wenn mehrere
      //dofsets vorhanden sind

      // Information
      if (iterstn->second.size()>4)
        cout << RED_LIGHT << "BUG: more standard sets!!!! "<< "Node GID " << bgnode->Id() << END_COLOR << endl;

      if (numsets > 1)
        cout << GREEN_LIGHT << "Info: more dofsets in transfer.. " <<  "Node GID " << bgnode->Id() << END_COLOR << endl;

      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp,bgstatevn);

    }
    // Project dofs from embdis to bgdis:
    // n:void -> n+1:std, n:enriched -> n+1:enriched,
    // n:enriched -> n+1: std, n:void ->  n+1:enriched
    else if (((iterstn == stdnoden_.end() and iteren == enrichednoden_.end()) and iterstnp != stdnodenp_.end()) or
             (iteren != enrichednoden_.end() and iterenp != enrichednodenp_.end()) or
             (iteren != enrichednoden_.end() and iterstnp != stdnodenp_.end()) or
             ((iterstn == stdnoden_.end() and iteren == enrichednoden_.end()) and iterenp != enrichednodenp_.end()))
    {
      // save the information needed from bg-nodes to call the method CommunicateNodes
      // set of are projected nodes
      projectednodeids_.insert(bgnode->Id());

      LINALG::Matrix<3,1> bgnodecords(true);
      bgnodecords(0,0) = bgnode->X()[0];
      bgnodecords(1,0) = bgnode->X()[1];
      bgnodecords(2,0) = bgnode->X()[2];

      // collect the coordinates of bg-nodes without history values
      bgnodes_coords.push_back(bgnodecords);

      // collect the bg-node ids which need a projection
      bgnodeidwithnohistory.push_back(bgnode->Id());

      // the vector of interpolated values
      LINALG::Matrix<4,1>    interpolatedvec(true);
      interpolated_vecs.push_back(interpolatedvec);

      int numsets = bgdis->NumDof(bgnode)/4;
      if( numsets > 1 )
        cout << GREEN_LIGHT << "Info: more dofsets in projection.. " <<  "Node GID " << bgnode->Id() << END_COLOR << endl;

    }
    //do nothing:
    //n: void->n+1: void, n:std->n+1:void, n:enriched->n+1:void
    else if ( (iterstn == stdnoden_.end() and iteren == enrichednoden_.end()
               and iterstnp == stdnodenp_.end() and iterenp == enrichednodenp_.end() ) or
              (iterstn != stdnoden_.end() and (iterstnp == stdnodenp_.end() and iterenp == enrichednodenp_.end())) or
              (iteren != enrichednoden_.end() and (iterenp == enrichednodenp_.end() and iterstnp == stdnodenp_.end())) )
    {

      if( bgdis->NumDof(bgnode) > 0 )
        cout << RED_LIGHT <<  "BUG:: in do nothig!!" << " Node GID " << bgnode->Id() <<  END_COLOR<< endl;
      //cout << "do nothing" << bgnode->Id() << endl; ;
    }
    else
      cout << "warum bin ich da?! " <<  bgdis->NumDof(bgnode)   << " " <<    bgnode->Id() <<  endl;
  }

  // call the Round Robin Communicator
  CommunicateNodes(bgdis,bgnodes_coords,interpolated_vecs,bgnodeidwithnohistory,embstatevn,aledispn,bgstatevnp,bgstatevn);

}//XFEM::XFluidFluidTimeIntegration::SetNewBgStatevectorFullProjection


//-----------------------------------------------------------------------------
// Transfer data in the round robin communication pattern:
//
// In this function the coordinates of the background-nodes without history
// are send in a round robin pattern to all processors to find the embedded
// -element where this nodes lay at the previous time step. Once this element
// is found we do the projection.
//
// bgnodes_coords(in):           coordinates of bg-nodes without history
// interpolated_vecs(in/out):    interpolated vector
// bgnodeidwithnohistory(in):    bg-node-ids with no history
// embstatevn(in):               vector from which we want to interpolate values
// aledispn(in):                 displacement vector of embedded-dis
// bgstatevnp(in/out):           epetra state vector which needs interpolation
// bgstatevn(in):                epetra state vector needed to transfer
//                               enriched values if no embedded element is found
//------------------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::CommunicateNodes(const RCP<DRT::Discretization>        bgdis,
                                                        std::vector<LINALG::Matrix<3,1> >   & bgnodes_coords,
                                                        std::vector<LINALG::Matrix<4,1> >   & interpolated_vecs,
                                                        std::vector<int>                    & bgnodeidwithnohistory,
                                                        Teuchos::RCP<Epetra_Vector>           embstatevn,
                                                        Teuchos::RCP<Epetra_Vector>           aledispn,
                                                        Teuchos::RCP<Epetra_Vector>           bgstatevnp,
                                                        Teuchos::RCP<Epetra_Vector>           bgstatevn)
 {

  // get number of processors and the current processors id
  int numproc=embdis_->Comm().NumProc();

  //information how many processors work at all
  vector<int> allproc(numproc);

  // create an exporter for point to point comunication
  DRT::Exporter exporter(embdis_->Comm());

  // necessary variables
  MPI_Request request;

  // define send and receive blocks
  vector<char> sblock;
  vector<char> rblock;

  // vector which identifies if a bg-node has already interpolated values
  vector<int> NodeDone;
  // initialize NodeDone with zeros at the beginning
  for(size_t i=0; i<bgnodeidwithnohistory.size(); ++i)
    NodeDone.push_back(0);

  //----------------------------------------------------------------------
  // communication is done in a round robin loop
  //----------------------------------------------------------------------
  for (int np=0;np<numproc+1;++np)
  {
    // in the first step, we cannot receive anything
    if(np >0)
    {
      ReceiveBlock(rblock,exporter,request);

      // Unpack info from the receive block from the last proc
      vector<LINALG::Matrix<3,1> > stuff_coord;
      vector<LINALG::Matrix<4,1> > stuff_interpolatedvecs;
      vector<int>  stuff_bgnodeidswithnohistory;
      vector<int>  stuff_nodedone;

      vector<char>::size_type position = 0;
      DRT::ParObject::ExtractfromPack(position,rblock,stuff_coord);
      DRT::ParObject::ExtractfromPack(position,rblock,stuff_interpolatedvecs);
      DRT::ParObject::ExtractfromPack(position,rblock,stuff_bgnodeidswithnohistory);
      DRT::ParObject::ExtractfromPack(position,rblock,stuff_nodedone);
      bgnodes_coords = stuff_coord;
      interpolated_vecs = stuff_interpolatedvecs;
      bgnodeidwithnohistory = stuff_bgnodeidswithnohistory;
      NodeDone = stuff_nodedone;
    }

    // in the last step, we keep everything on this proc
    if(np < numproc)
    {
      // -----------------------
      // do what we wanted to do
      FindEmbeleAndInterpolatevalues(bgnodes_coords,interpolated_vecs,NodeDone,embstatevn,aledispn);

      // Pack info into block to sendit
      PackValues(bgnodes_coords,interpolated_vecs,bgnodeidwithnohistory,NodeDone,sblock);

      // add size  to sendblock
      SendBlock(sblock,exporter,request);
    }
  } // end of loop over processors

  //set the interpolated values to bgstatevnp
  for (size_t i=0; i<bgnodeidwithnohistory.size(); ++i)
  {
    DRT::Node* bgnode = bgdis->gNode(bgnodeidwithnohistory.at(i));
    // number of dof-sets
    int numsets = bgdis->NumDof(bgnode)/4;

    if( numsets > 1 )
      cout << GREEN_LIGHT << "Info: more dofsets in projection.. " <<  "Node GID " << bgnode->Id() << END_COLOR << endl;
    int offset = 0;

    // if interpolated values are available
    if (NodeDone.at(i) == 1)
    {
      for (int set=0; set<numsets; set++)
      {
        (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[offset+0])] = interpolated_vecs.at(i)(0);
        (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[offset+1])] = interpolated_vecs.at(i)(1);
        (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[offset+2])] = interpolated_vecs.at(i)(2);
        (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[offset+3])] = interpolated_vecs.at(i)(3);
        offset += 4;
      }
    }
    // if no embedded element is found try to find an enriched value
    else
    {
      cout << RED_LIGHT << "NodeDone " << NodeDone.at(i) << END_COLOR << endl;

      map<int, vector<int> >::const_iterator iterstn = stdnoden_.find(bgnode->Id());
      map<int, vector<int> >::const_iterator iterstnp = stdnodenp_.find(bgnode->Id());
      map<int, vector<int> >::const_iterator iteren = enrichednoden_.find(bgnode->Id());
      map<int, vector<int> >::const_iterator iterenp = enrichednodenp_.find(bgnode->Id());

      if ((iteren != enrichednoden_.end() and iterenp != enrichednodenp_.end())
          or ((iteren != enrichednoden_.end() and iterstnp != stdnodenp_.end())))
      {
        cout << RED_LIGHT <<  "CHECK: Took enriched values !!" << " Node GID " << bgnode->Id() << END_COLOR<< endl;
        vector<int> gdofsn = iteren->second;

        WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp,bgstatevn);
      }
      else
      {
        cout << YELLOW_LIGHT << " Warning: No patch element found for the node " << bgnode->Id() ;
        if ((iterstn == stdnoden_.end() and iteren == enrichednoden_.end()) and iterstnp != stdnodenp_.end())
          cout << YELLOW_LIGHT << " n:void -> n+1:std  " << END_COLOR<< endl;
        else if (iteren != enrichednoden_.end() and iterenp != enrichednodenp_.end())
          cout << YELLOW_LIGHT << " n:enriched -> n+1:enriched " << END_COLOR<< endl;
        else if (iteren != enrichednoden_.end() and iterstnp != stdnodenp_.end())
          cout << YELLOW_LIGHT << " n:enriched -> n+1: std " << END_COLOR<< endl;
        else if ((iterstn == stdnoden_.end() and iteren == enrichednoden_.end()) and iterenp != enrichednodenp_.end())
          cout << YELLOW_LIGHT << " n:void ->  n+1:enriched " << END_COLOR<< endl;
      }
    }
  }// end of loop over bgnodes without history

 }//XFEM::XFluidFluidTimeIntegration::CommunicateNodes

//---------------------------------------------------------
// receive a block in the round robin communication pattern
//---------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::ReceiveBlock(vector<char>   & rblock,
                                                    DRT::Exporter  & exporter,
                                                    MPI_Request    & request)
{
  // get number of processors and the current processors id
  int numproc=embdis_->Comm().NumProc();
  int myrank =embdis_->Comm().MyPID();

  // necessary variables
  int         length =-1;
  int         frompid=(myrank+numproc-1)%numproc;
  int         tag    =frompid;


  // receive from predecessor
  exporter.ReceiveAny(frompid,tag,rblock,length);

// #ifdef DEBUG
//   cout << "----receiving " << rblock.size() <<  " bytes: to proc " << myrank << " from proc " << frompid << endl;
// #endif

  if(tag!=(myrank+numproc-1)%numproc)
  {
    dserror("received wrong message (ReceiveAny)");
  }

  exporter.Wait(request);

  // for safety
  exporter.Comm().Barrier();

  return;
} // XFEM::XFluidFluidTimeIntegration::ReceiveBlock

//---------------------------------------------------------
// send a block in the round robin communication pattern
//---------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SendBlock(vector<char>  & sblock  ,
                                                 DRT::Exporter & exporter,
                                                 MPI_Request   & request )
{
  // get number of processors and the current processors id
  int numproc=embdis_->Comm().NumProc();
  int myrank =embdis_->Comm().MyPID();

  // Send block to next proc.
  int         tag    =myrank;
  int         frompid=myrank;
  int         topid  =(myrank+1)%numproc;

// #ifdef DEBUG
//   cout << "----sending " << sblock.size() <<  " bytes: from proc " << myrank << " to proc " << topid << endl;
// #endif

  exporter.ISend(frompid,topid,
                 &(sblock[0]),sblock.size(),
                 tag,request);

  // for safety
  exporter.Comm().Barrier();

  return;
} // XFluidFluidTimeIntegration::SendBlock

//---------------------------------------------------------
// pack values in the round robin communication pattern
//---------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::PackValues(std::vector<LINALG::Matrix<3,1> > & bgnodes_coords,
                                                  std::vector<LINALG::Matrix<4,1> > & interpolatedvec,
                                                  std::vector<int>                  &  bgnodeidwithnohistory,
                                                  vector<int>                       & NodeDone,
                                                  vector<char>                      & sblock)
{
  // Pack info into block to send
  DRT::PackBuffer data;
  DRT::ParObject::AddtoPack(data,bgnodes_coords);
  DRT::ParObject::AddtoPack(data,interpolatedvec);
  DRT::ParObject::AddtoPack(data,bgnodeidwithnohistory);
  DRT::ParObject::AddtoPack(data,NodeDone);
  data.StartPacking();

  DRT::ParObject::AddtoPack(data,bgnodes_coords);
  DRT::ParObject::AddtoPack(data,interpolatedvec);
  DRT::ParObject::AddtoPack(data,bgnodeidwithnohistory);
  DRT::ParObject::AddtoPack(data,NodeDone);
  swap( sblock, data() );

} // XFluidFluidTimeIntegration::PackValue

//--------------------------------------------------------
//
//--------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::FindEmbeleAndInterpolatevalues(std::vector<LINALG::Matrix<3,1> > & bgnodes_coords,
                                                                      std::vector<LINALG::Matrix<4,1> > & interpolated_vecs,
                                                                      vector<int>                       & NodeDone,
                                                                      Teuchos::RCP<Epetra_Vector>        embstatevn,
                                                                      Teuchos::RCP<Epetra_Vector>        aledispn)

{
  // loop over bgnodes_coords
   for (size_t i=0; i<bgnodes_coords.size(); ++i)
  {
    // bgnode coordinate
    LINALG::Matrix<3,1> bgnodecords = bgnodes_coords.at(i);
    // interpolated vector which is zero at the beginning
    LINALG::Matrix<4,1> interpolatedvec = interpolated_vecs.at(i);

    bool insideelement = false;
    // check all embedded elements to find the right one
    for (int e=0; e<embdis_->NumMyColElements(); e++)
    {
      DRT::Element* pele = embdis_->lColElement(e);
      insideelement = ComputeSpacialToElementCoordAndProject(pele,bgnodecords,interpolatedvec,embstatevn,aledispn,embdis_);
      if (insideelement and NodeDone.at(i)==0)
      {
        // Set NodeDone to 1 if the embedded-element is found
        NodeDone.at(i) = 1;

        //set the interpolated values
        interpolated_vecs.at(i)(0) = interpolatedvec(0);
        interpolated_vecs.at(i)(1) = interpolatedvec(1);
        interpolated_vecs.at(i)(2) = interpolatedvec(2);
        interpolated_vecs.at(i)(3) = interpolatedvec(3);
      }
    }
  }
}

// -------------------------------------------------------------------
// If enriched values are available keep them (no projection from
// embedded fluid)
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SetNewBgStatevectorKeepGhostValues(const RCP<DRT::Discretization>        bgdis,
                                                                          Teuchos::RCP<Epetra_Vector>           bgstatevn,
                                                                          Teuchos::RCP<Epetra_Vector>           bgstatevnp,
                                                                          Teuchos::RCP<Epetra_Vector>           embstatevn,
                                                                          Teuchos::RCP<Epetra_Vector>           aledispn)
{

  // coordinates of bg-nodes which need the projection from embedded dis
  std::vector<LINALG::Matrix<3,1> > bgnodes_coords;
  // the vector containing interpolated values from embedded dis
  std::vector<LINALG::Matrix<4,1> > interpolated_vecs;
  // bg-node ids which have no history and need a projection
  std::vector<int> bgnodeidwithnohistory;

  // loop over bg-row-nodes of each processor
  for (int lnid=0; lnid<bgdis->NumMyRowNodes(); lnid++)
  {
    DRT::Node* bgnode = bgdis->lRowNode(lnid);
    map<int, vector<int> >::const_iterator iterstn = stdnoden_.find(bgnode->Id());
    map<int, vector<int> >::const_iterator iterstnp = stdnodenp_.find(bgnode->Id());
    map<int, vector<int> >::const_iterator iteren = enrichednoden_.find(bgnode->Id());
    map<int, vector<int> >::const_iterator iterenp = enrichednodenp_.find(bgnode->Id());

    // Transfer the dofs:
    // n:std -> n+1:std, n:std -> n+1:enriched
    if ((iterstn != stdnoden_.end() and iterstnp != stdnodenp_.end()) or
        (iterstn != stdnoden_.end() and iterenp != enrichednodenp_.end()))
    {
      int numsets = bgdis->NumDof(bgnode)/4;
      vector<int> gdofsn = iterstn->second;

      //TODO!! die richtige dofs von bgstatevn rauspicke, wenn mehrere
      //dofsets vorhanden sind

      // Information
      if (iterstn->second.size()>4)
        cout << RED_LIGHT << "BUG: more standard sets!!!! "<< "Node GID " << bgnode->Id() << END_COLOR << endl;

      if (numsets > 1)
        cout << GREEN_LIGHT << "Info: more dofsets in transfer.. " <<  "Node GID " << bgnode->Id() << END_COLOR << endl;

      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp,bgstatevn);

    }
    // Project dofs from embdis to bgdis:
    // n:void -> n+1:std, n:void ->  n+1:enriched
    else if (((iterstn == stdnoden_.end() and iteren == enrichednoden_.end()) and iterstnp != stdnodenp_.end()) or
             ((iterstn == stdnoden_.end() and iteren == enrichednoden_.end()) and iterenp != enrichednodenp_.end()))
    {
      // save the information needed from bg-nodes to call the method CommunicateNodes
      // set of are projected nodes
      projectednodeids_.insert(bgnode->Id());

      LINALG::Matrix<3,1> bgnodecords(true);
      bgnodecords(0,0) = bgnode->X()[0];
      bgnodecords(1,0) = bgnode->X()[1];
      bgnodecords(2,0) = bgnode->X()[2];

      // collect the coordinates of bg-nodes without history values
      bgnodes_coords.push_back(bgnodecords);

      // collect the bg-node ids which need a projection
      bgnodeidwithnohistory.push_back(bgnode->Id());

      // the vector of interpolated values
      LINALG::Matrix<4,1>    interpolatedvec(true);
      interpolated_vecs.push_back(interpolatedvec);

      int numsets = bgdis->NumDof(bgnode)/4;

      if( numsets > 1 )
        cout << GREEN_LIGHT << "Info: more dofsets in projection.. " <<  "Node GID " << bgnode->Id() << END_COLOR << endl;

    }
    //keep the ghost dofs:
    //n: enriched -> n+1:enriched, n: enriched -> n+1: std
    else if ((iteren != enrichednoden_.end() and iterenp != enrichednodenp_.end()) or
             (iteren != enrichednoden_.end() and iterstnp != stdnodenp_.end()))
    {
#ifdef DOFSETS_NEW
      int numsets = bgdis->NumDof(bgnode)/4;
      if (numsets > 1)
        cout << RED_LIGHT << "ghost-fluid-approach just available for one dofset!" << endl;

      vector<int> gdofsn = iteren->second;
      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp,bgstatevn);
#else
      vector<int> gdofsn = iteren->second;
      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp,bgstatevn);

#endif
    }
    //do nothing:
    //n: void->n+1: void, n:std->n+1:void, n:enriched->n+1:void
    else if ( (iterstn == stdnoden_.end() and iteren == enrichednoden_.end()
               and iterstnp == stdnodenp_.end() and iterenp == enrichednodenp_.end() ) or
              (iterstn != stdnoden_.end() and (iterstnp == stdnodenp_.end() and iterenp == enrichednodenp_.end())) or
              (iteren != enrichednoden_.end() and (iterenp == enrichednodenp_.end() and iterstnp == stdnodenp_.end())) )
    {

      if( bgdis->NumDof(bgnode) > 0 )
        cout << RED_LIGHT <<  "BUG:: in do nothig!!" << " Node GID " << bgnode->Id() <<  END_COLOR<< endl;
      //cout << "do nothing" << bgnode->Id() << endl; ;
    }
    else
      cout << "warum bin ich da?! " <<  bgdis->NumDof(bgnode)   << " " <<    bgnode->Id() <<  endl;
  }

  // call the Round Robin Communicator
  CommunicateNodes(bgdis,bgnodes_coords,interpolated_vecs,bgnodeidwithnohistory,embstatevn,aledispn,bgstatevnp,bgstatevn);

}//XFEM::XFluidFluidTimeIntegration::SetNewBgStatevectorKeepGhostValues

//-------------------------------------------------------------------
// Write the values of node from bgstatevn to bgstatevnp
//--------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::WriteValuestoBgStateVector(const RCP<DRT::Discretization>   bgdis,
                                                                  DRT::Node*                       bgnode,
                                                                  vector<int>                      gdofs_n,
                                                                  Teuchos::RCP<Epetra_Vector>      bgstatevnp,
                                                                  Teuchos::RCP<Epetra_Vector>      bgstatevn)
{
  int numsets = bgdis->NumDof(bgnode)/4;

  if (numsets > 1)
    cout << GREEN_LIGHT << "Info: more dofsets in transfer.. " <<  "Node GID " << bgnode->Id() << END_COLOR << endl;

  int offset = 0;
  for (int set=0; set<numsets; set++)
  {
    (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[offset+0])] =
      (*bgstatevn)[bgstatevn->Map().LID(gdofs_n[0])];
    (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[offset+1])] =
      (*bgstatevn)[bgstatevn->Map().LID(gdofs_n[1])];
    (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[offset+2])] =
      (*bgstatevn)[bgstatevn->Map().LID(gdofs_n[2])];
    (*bgstatevnp)[bgstatevnp->Map().LID(bgdis->Dof(bgnode)[offset+3])] =
      (*bgstatevn)[bgstatevn->Map().LID(gdofs_n[3])];
    offset += 4;
  }

}//XFEM::XFluidFluidTimeIntegration::WriteValuestoBgStateVector

// -------------------------------------------------------------------
// In this function:
//
// - pele (in): The element where we want to extract our values from
// - x    (in): The coordinates of the node without values (x,y,z)
// - interpolatedvec (out): The vector after interpolation
// - embstate_n (in): Source state vector
// - embeddedddisp (in): The displacement of the embedded fluid at the
//   time that it is considered as source of the interpolation
// - sourcedis (in): the source discretization where we get the values
//   from
//
//   Check whether the unknown node lies in the pele. If yes fill the
//   interpolatedvec and return true.
//
// -------------------------------------------------------------------
bool XFEM::XFluidFluidTimeIntegration::ComputeSpacialToElementCoordAndProject(DRT::Element*                       pele,
                                                                              LINALG::Matrix<3,1>&                x,
                                                                              LINALG::Matrix<4,1>&                interpolatedvec,
                                                                              Teuchos::RCP<Epetra_Vector>         embstate_n,
                                                                              Teuchos::RCP<Epetra_Vector>         embeddeddisp,
                                                                              Teuchos::RCP<DRT::Discretization>   sourcedis)
{
  // numnodes of the embedded element
  const std::size_t numnode = pele->NumNode();
  LINALG::SerialDenseMatrix pxyze(3,numnode);
  DRT::Node** pelenodes = pele->Nodes();

  std::vector<double> myval(4);
  std::vector<double> mydisp(4);
  std::vector<int> pgdofs(4);
  LINALG::Matrix<3,1> xsi(true);

  bool inside = false;

  switch ( pele->Shape() )
  {
    case DRT::Element::hex8:
    {
      const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
      //const int numnodes = pele->NumNode();
      LINALG::Matrix<4,numnodes> veln(true);
      LINALG::Matrix<4,numnodes> disp;

      // loop over the nodes of the patch element and save the veln and
      // coordinates of all nodes
      for (int inode = 0; inode < numnodes; ++inode)
      {
        sourcedis->Dof(pelenodes[inode],0,pgdofs);
        DRT::UTILS::ExtractMyValues(*embstate_n,myval,pgdofs);
        veln(0,inode) = myval[0];
        veln(1,inode) = myval[1];
        veln(2,inode) = myval[2];
        veln(3,inode) = myval[3];

        // get the coordinates of patch element
        pxyze(0,inode) = pelenodes[inode]->X()[0];
        pxyze(1,inode) = pelenodes[inode]->X()[1];
        pxyze(2,inode) = pelenodes[inode]->X()[2];

        if (sourcedis->Name() == "xfluid") // embedded fluid
        {
          DRT::UTILS::ExtractMyValues(*embeddeddisp,mydisp,pgdofs);
          disp(0,inode) = mydisp[0];
          disp(1,inode) = mydisp[1];
          disp(2,inode) = mydisp[2];

          // add the current displacement to the coordinates
          pxyze(0,inode) += disp(0,inode);
          pxyze(1,inode) += disp(1,inode);
          pxyze(2,inode) += disp(2,inode);
        }
      }

      // check whether the xfemnode (the node with no values) is in the element
      LINALG::Matrix< 3,numnodes > xyzem(pxyze,true);
      GEO::CUT::Position <DRT::Element::hex8> pos(xyzem,x);
      double tol = 1e-10;
      bool insideelement = pos.ComputeTol(tol);
      //bool insideelement = pos.Compute();

      // von ursula
      // bool in  = GEO::currentToVolumeElementCoordinates(DRT::Element::hex8, pxyze, x, xsi);
      // in  = GEO::checkPositionWithinElementParameterSpace(xsi, DRT::Element::hex8);

      if (insideelement)
      {
        // get the coordinates of x in element coordinates of patch element pele (xsi)
        xsi = pos.LocalCoordinates();
        // evaluate shape function
        LINALG::SerialDenseVector shp(numnodes);
        DRT::UTILS::shape_function_3D( shp, xsi(0,0), xsi(1,0), xsi(2,0), DRT::Element::hex8 );
        // Interpolate
        for (int inode = 0; inode < numnodes; ++inode){
          for (std::size_t isd = 0; isd < 4; ++isd){
            interpolatedvec(isd) += veln(isd,inode)*shp(inode);
          }
        }
        inside = true;
        break;
      }
      else
      {
        inside = false;
        break;
      }
    }
    case DRT::Element::hex20:
    case DRT::Element::hex27:
    {
      dserror("No support for hex20 and hex27!");
      break;
    }
    default:
    {
      dserror("Element-type not supported here!");
    }
  }
  return inside;
}//ComputeSpacialToElementCoordAndProject

// -------------------------------------------------------------------
// Needed for interpolation-based fluid-fluid-fsi
//
// In this function:
//
// - pele (in): The element where we want to extract our values from
// - x    (in): The coordinates of the node without values (x,y,z)
// - interpolatedvec (out): The vector after interpolation
// - embstate_n (in): Source state vector
// - embeddedddisp (in): The displacement of the embedded fluid at the
//   time that it is considered as source of the interpolation
// - sourcedis (in): the source discretization where we get the values
//   from
//
//   Check whether the unknown node lies in the pele. If yes fill the
//   interpolatedvec and return true.
//
// -------------------------------------------------------------------
bool XFEM::XFluidFluidTimeIntegration::ComputeSpacialToElementCoordAndProject2(DRT::Element*                       pele,
                                                                               LINALG::Matrix<3,1>&                x,
                                                                               LINALG::Matrix<4,1>&                interpolatedvec,
                                                                               Teuchos::RCP<Epetra_Vector>         embstate_n,
                                                                               Teuchos::RCP<Epetra_Vector>         embeddeddisp,
                                                                               Teuchos::RCP<DRT::Discretization>   sourcedis)
{
  const std::size_t numnode = pele->NumNode();
  LINALG::SerialDenseMatrix pxyze(3,numnode);
  DRT::Node** pelenodes = pele->Nodes();

  std::vector<double> myval(4);
  std::vector<double> mydisp(4);
  std::vector<int> pgdofs(4);
  LINALG::Matrix<3,1> xsi(true);

  bool inside = false;

  switch ( pele->Shape() )
  {
    case DRT::Element::hex8:
    {
//      cout <<" bg ele id" << pele->Id() << endl;
      const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
      LINALG::Matrix<4,numnodes> veln(true);
      LINALG::Matrix<4,numnodes> disp;

      for (int inode = 0; inode < numnodes; ++inode)
      {
        // the enriched nodes are also included
        if (sourcedis->NumDof(pelenodes[inode]) > 0)
        {
          sourcedis->Dof(pelenodes[inode],0,pgdofs);
          DRT::UTILS::ExtractMyValues(*embstate_n,myval,pgdofs);
          veln(0,inode) = myval[0];
          veln(1,inode) = myval[1];
          veln(2,inode) = myval[2];
          veln(3,inode) = myval[3];

          // get the coordinates of patch element
          pxyze(0,inode) = pelenodes[inode]->X()[0];
          pxyze(1,inode) = pelenodes[inode]->X()[1];
          pxyze(2,inode) = pelenodes[inode]->X()[2];
        }
      }


      // check whether the node with no values is in the element
      LINALG::Matrix< 3,numnodes > xyzem(pxyze,true);
      GEO::CUT::Position <DRT::Element::hex8> pos(xyzem,x);
      double tol = 1e-10;
      bool insideelement = pos.ComputeTol(tol);

      if (insideelement)
      {
        // get the coordinates of x in element coordinates of patch element pele (xsi)
        xsi = pos.LocalCoordinates();
        // evaluate shape function
        LINALG::SerialDenseVector shp(numnodes);
        DRT::UTILS::shape_function_3D( shp, xsi(0,0), xsi(1,0), xsi(2,0), DRT::Element::hex8 );
        // Interpolate
        for (int inode = 0; inode < numnodes; ++inode){
          for (std::size_t isd = 0; isd < 4; ++isd){
            interpolatedvec(isd) += veln(isd,inode)*shp(inode);
          }
        }
        inside = true;
        break;
      }
      else
      {
        inside = false;
        break;
      }
    }
    case DRT::Element::hex20:
    case DRT::Element::hex27:
    {
      dserror("No support for hex20 and hex27!");
      break;
    }
    default:
    {
      dserror("Element-type not supported here!");
    }
  }
  return inside;
}//ComputeSpacialToElementCoordAndProject2

// -------------------------------------------------------------------
// Needed for interpolation-based fluid-fluid-fsi
//
// In this function:
// loop over all embedded fluid nodes and find the embedded element
// where this node was located before we solved the ale again.(the
// displacement is saved at aledispnpoldstate_). If found interpolate
// the values on the new embedded fluid vector. If not found we have
// to look for a matching element in the background fluid. Right now
// we just take the old values of the embedded mesh: TODO
//
// statevbg_n (in): source state vector from gbackground fluid
// statevemb_n (in): source state vector from embedded fluid
// statevembnew_n (out)
// aledispnp (in): displacement of ale at the current time (time of
//               interpolation)
//
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SetNewEmbStatevector(const RCP<DRT::Discretization> bgdis,
                                                            Teuchos::RCP<Epetra_Vector>    statevbg_n,
                                                            Teuchos::RCP<Epetra_Vector>    statevemb_n,
                                                            Teuchos::RCP<Epetra_Vector>    statevembnew_n,
                                                            Teuchos::RCP<Epetra_Vector>    aledispnp,
                                                            Teuchos::RCP<Epetra_Vector>    aledispnpoldstate)
{
  std::vector<double> mydisp(4);
  std::vector<int> pgdofs(4);

  // Gmsh----------------------------------------------------------
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("emb_element_node_id", 0, 0, 0, bgdis->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  {
    {
      // draw embedded elements with associated gid at the old  position
      gmshfilecontent << "View \" " << "emb Element(old)->Id() \" {\n";
      for (int i=0; i<embdis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = embdis_->lColElement(i);
        const DRT::Node*const* pelenodes = actele->Nodes();
        std::map<int,LINALG::Matrix<3,1> > mapofnodepos; //node id-> position

        vector<int> lm;
        vector<int> lmowner;
        vector<int> lmstride;
        actele->LocationVector(*embdis_, lm, lmowner, lmstride);

        vector<double> myolddisp(lm.size());
        DRT::UTILS::ExtractMyValues(*aledispnpoldstate, myolddisp, lm);

        for (int inode = 0; inode < actele->NumNode(); ++inode)
        {
          // the coordinates of the actuall node
          LINALG::Matrix<3,1> inodepos(true);
          inodepos(0,0) = pelenodes[inode]->X()[0] + myolddisp[0+(inode*4)];
          inodepos(1,0) = pelenodes[inode]->X()[1] + myolddisp[1+(inode*4)];
          inodepos(2,0) = pelenodes[inode]->X()[2] + myolddisp[2+(inode*4)];

          mapofnodepos[pelenodes[inode]->Id()] = inodepos;
        }
        IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, mapofnodepos, gmshfilecontent);
      };
      gmshfilecontent << "};\n";
    }

    {
      // draw embedded elements with associated gid at the current position
      gmshfilecontent << "View \" " << "emb Element(new)->Id() \" {\n";
      for (int i=0; i<embdis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = embdis_->lColElement(i);
        const DRT::Node*const* pelenodes = actele->Nodes();
        std::map<int,LINALG::Matrix<3,1> > mapofnodepos; //node id-> position

        vector<int> lm;
        vector<int> lmowner;
        vector<int> lmstride;
        actele->LocationVector(*embdis_, lm, lmowner, lmstride);

        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*aledispnp, mydisp, lm);

        for (int inode = 0; inode < actele->NumNode(); ++inode)
        {
          // the coordinates of the actuall node
          LINALG::Matrix<3,1> inodepos(true);
          inodepos(0,0) = pelenodes[inode]->X()[0] + mydisp[0+(inode*4)];
          inodepos(1,0) = pelenodes[inode]->X()[1] + mydisp[1+(inode*4)];
          inodepos(2,0) = pelenodes[inode]->X()[2] + mydisp[2+(inode*4)];

          mapofnodepos[pelenodes[inode]->Id()] = inodepos;
        }
        IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, mapofnodepos, gmshfilecontent);
      };
      gmshfilecontent << "};\n";
    }

    {
      // draw embedded nodes with associated gid at the current position
      gmshfilecontent << "View \" " << "emb Node(new)->Id() \" {\n";
      for (int i=0; i<embdis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = embdis_->lColElement(i);
        const DRT::Node*const* pelenodes = actele->Nodes();
        std::map<int,LINALG::Matrix<3,1> > mapofnodepos; //node id-> position

        vector<int> lm;
        vector<int> lmowner;
        vector<int> lmstride;
        actele->LocationVector(*embdis_, lm, lmowner, lmstride);

        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*aledispnp, mydisp, lm);

        for (int inode = 0; inode < actele->NumNode(); ++inode)
        {
          // the coordinates of the actuall node
          LINALG::Matrix<3,1> inodepos(true);
          inodepos(0,0) = pelenodes[inode]->X()[0] + mydisp[0+(inode*4)];
          inodepos(1,0) = pelenodes[inode]->X()[1] + mydisp[1+(inode*4)];
          inodepos(2,0) = pelenodes[inode]->X()[2] + mydisp[2+(inode*4)];

          IO::GMSH::cellWithScalarToStream(DRT::Element::point1, pelenodes[inode]->Id(), inodepos, gmshfilecontent);
        }
      };
      gmshfilecontent << "};\n";
    }
  }
  gmshfilecontent.close();
  // --------------------------------------------------------------

  for (int lnid=0; lnid<embdis_->NumMyRowNodes(); lnid++)
  {
    DRT::Node* embnode = embdis_->lRowNode(lnid);

    DRT::UTILS::ConditionSelector conds(*embdis_, "FSICoupling");
    if (conds.ContainsNode(embnode->Id())==false)//if no fsi node
    {
      //cout << "unknown node id " << embnode->Id() << endl;
      embdis_->Dof(embnode,0,pgdofs);

      DRT::UTILS::ExtractMyValues(*aledispnp,mydisp,pgdofs);

      // the coordinates of the unknown node
      LINALG::Matrix<3,1> embnodecords(true);
      embnodecords(0,0) = embnode->X()[0] + mydisp[0];
      embnodecords(1,0) = embnode->X()[1] + mydisp[1];
      embnodecords(2,0) = embnode->X()[2] + mydisp[2];

      bool insideelement = false;
      int count = 0;
      // check all embedded elements to find the right one

      for (int e=0; e<embdis_->NumMyColElements(); e++)
      {
        DRT::Element* pele = embdis_->lColElement(e);
        LINALG::Matrix<4,1> interpolatedvec(true);
        //cout << "pele id " << pele->Id() << endl;

        insideelement = ComputeSpacialToElementCoordAndProject(pele,embnodecords,interpolatedvec,statevemb_n,aledispnpoldstate,embdis_);
        if (insideelement)
        {
          // here set state
          (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[0])] = interpolatedvec(0);
          (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[1])] = interpolatedvec(1);
          (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[2])] = interpolatedvec(2);
          (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[3])] = interpolatedvec(3);
          break;
        }
        count ++;
      }
      if (count == embdis_->NumMyColElements()) //if no embedded elements found
      {
        cout << "no embedded element found for: " <<  embnode->Id() << endl;
         count = 0;
//         // check all background elements to find the right one
         for (int e=0; e<bgdis->NumMyColElements(); e++)
         {
            DRT::Element* bgele = bgdis->lColElement(e);
            LINALG::Matrix<4,1> interpolatedvec(true);

            insideelement = ComputeSpacialToElementCoordAndProject2(bgele,embnodecords,interpolatedvec,statevbg_n,aledispnpoldstate,bgdis);
//            if (insideelement)
//            {
//              // here set state
//              (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[0])] = interpolatedvec(0);
//              (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[1])] = interpolatedvec(1);
//              (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[2])] = interpolatedvec(2);
//              (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[3])] = interpolatedvec(3);
//              break;
//            }
//            count ++;
        }
//          if (count == bgdis->NumMyColElements())
//            cout << "no bg elements found for: " << embnode->Id() << endl;
        }
    }
  }

}//SetNewEmbStatevector

// -------------------------------------------------------------------
// Gmsh Output
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::GmshOutput(const RCP<DRT::Discretization>        bgdis)
{
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("std_enriched_map", step_, 30, 0, bgdis->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  {
    gmshfilecontent << "View \" " << "std/enriched/void n\" {\n";
    for (int i=0; i<bgdis->NumMyColNodes(); ++i)
    {
      int kind = 0;
      const DRT::Node* actnode = bgdis->lColNode(i);
      const LINALG::Matrix<3,1> pos(actnode->X());
      map<int, vector<int> >::const_iterator iter = stdnoden_.find(actnode->Id());
      map<int, vector<int> >::const_iterator iteren = enrichednoden_.find(actnode->Id());
      if (iter != stdnoden_.end()) kind = 1;//std
      if (iteren != enrichednoden_.end()) kind = 2; // enriched
      IO::GMSH::cellWithScalarToStream(DRT::Element::point1, kind, pos, gmshfilecontent);
    }
    gmshfilecontent << "};\n";
  }
  {
    gmshfilecontent << "View \" " << "std/enriched/void n+1\" {\n";
    for (int i=0; i<bgdis->NumMyColNodes(); ++i)
    {
      int kind = 0;
      const DRT::Node* actnode = bgdis->lColNode(i);
      const LINALG::Matrix<3,1> pos(actnode->X());
      map<int, vector<int> >::const_iterator iter = stdnodenp_.find(actnode->Id());
      map<int, vector<int> >::const_iterator iteren = enrichednodenp_.find(actnode->Id());
      if (iter != stdnodenp_.end()) kind = 1;//std
      if (iteren != enrichednodenp_.end()) kind = 2; // enriched
      IO::GMSH::cellWithScalarToStream(DRT::Element::point1, kind, pos, gmshfilecontent);
    }
    gmshfilecontent << "};\n";
  }
}

// -----------------------------------------------------------------------
// all cut elements at tn with full dofs are included in incompressibility
// patch. If there  are nodes in projectednodeids which still are not in
// incompressibility patch we add them afterwards. These added Elements should
// have the full dofs again.
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::PatchelementForIncompressibility(const RCP<DRT::Discretization>     bgdis,
                                                                        XFEM::FluidWizard                  wizard_n,
                                                                        XFEM::FluidWizard                  wizard_np,
                                                                        Teuchos::RCP<LINALG::MapExtractor> dbcmaps)
{
  ///////////////////////////////////////////////////
  // finde den patch at time t_n

  // some temporary variables to build the patch
  std::vector<int> incompdofs;
  std::vector<int> incompdofs_vel;
  std::vector<int> incompdofs_vel_without_dirichlet;

  std::set<int> incompdofs_set;
  std::set<int> incompdofs_vel_set;
  std::set<int> incompdofs_vel_without_dirichlet_set;

  std::set<int> incompelementids_set;
  std::vector<int> incompelementids;

  std::set<int> incompnodeids_set;

  // call loop over elements
  const int numele = bgdis->NumMyColElements();
  for (int i=0; i<numele; ++i)
  {
    DRT::Element* actele = bgdis->lRowElement(i);
    int numnodes = actele->NumNode();

    // get the element handle of actele at time tn
    GEO::CUT::ElementHandle * e = wizard_n.GetElement( actele );

    const DRT::Node* const* nodesofele = actele->Nodes();

    // We are looking for elements whose all nodes have dofs at tn+1.
    // So check whether the nodes of this elements have dofs
    int numnodeswithdofs = 0;
    for(int inode = 0; inode < numnodes; inode++)
    {
      if(bgdis->NumDof(nodesofele[inode]) > 0)
        numnodeswithdofs++;
    }

    // xfem element
    if ( e!=NULL )
    {
      // select all cut elements at t_n, which has complete dofs at t_n+1
      if (e->IsCut() and (numnodeswithdofs == numnodes))
      {
        //insert the element
        incompelementids_set.insert(actele->Id());

        //get the nodes adjust to this element
        const int* nodeids = actele->NodeIds();

        for(int inode=0; inode<actele->NumNode(); ++inode)
        {
          const int nodegid = nodeids[inode];

          //insert the node
          incompnodeids_set.insert(nodegid);

          //get the dofs of the node
          const std::vector<int> mydofs(bgdis->Dof(bgdis->gNode(nodegid)));

          // insert the dofs
          for (size_t k=0; k<mydofs.size(); ++k)
          {
            incompdofs_set.insert(mydofs.at(k));
          }

          // insert just the velocity dofs
          int numsets = mydofs.size()/4;

          incompdofs_vel_set.insert(mydofs.at((numsets-1)+0));
          incompdofs_vel_set.insert(mydofs.at((numsets-1)+1));
          incompdofs_vel_set.insert(mydofs.at((numsets-1)+2));
        }
      }
    }
  }

  // debug output
//   for(set<int>::iterator iter = incompelementids_set.begin(); iter!= incompelementids_set.end();
//       iter++)
//   {
//     set<int>::const_iterator iterpatchele = incompelementids_set.find(*iter);
//     cout << "all eles vorher" << *iter << endl;
//     DRT::Node ** elenodes = bgdis->gElement(*iter)->Nodes();
//     for(int inode=0; inode<bgdis->gElement(*iter)->NumNode(); ++inode)
//     {
//       cout << bgdis->NumDof(elenodes[inode]) << endl;
//     }
//   }

  ////////////////////////////////////////////////////
  //  check if all projected nodes are included
  for(set<int>::iterator iter = projectednodeids_.begin(); iter!= projectednodeids_.end();
      iter++)
  {
    set<int>::const_iterator iterpatchnodes = incompnodeids_set.find(*iter);

    if (iterpatchnodes == incompnodeids_set.end())
    {
      cout << "STOP!! Nodes found which were in projected set but not in incompressibility patch!" << *iter << endl;

      //get all adjacent elements of this node
      int numberOfElements = bgdis->gNode(*iter)->NumElement();
      const DRT::Element* const* elements = bgdis->gNode(*iter)->Elements();

      // loop over adjacent elements of this node
      for(int ele_current=0; ele_current<numberOfElements; ele_current++)
      {
        // get the adjacent element
        const DRT::Element* ele_adj = elements[ele_current];

        // get vector of pointers of node (for this element)
        const DRT::Node* const* nodesofadjele = ele_adj->Nodes();
        const int numberOfNodes = ele_adj->NumNode();

        int numnodeswithdofs = 0;
        // loop nodes of this element to know whether it has complete dofs
        for(int vec_it = 0; vec_it < numberOfNodes; vec_it++)
        {
          // check if the nodes of this elements have dofs
          if(bgdis->NumDof(nodesofadjele[vec_it]) > 0)
            numnodeswithdofs++;
        }

        // if all nodes of this elements have dofs add it to incomp patch
        if (numnodeswithdofs == numberOfNodes )
        {
          incompelementids_set.insert(ele_adj->Id());
          cout << "element found " << ele_adj->Id() << endl;

          //get the nodes of this element
          const int* nodeids = ele_adj->NodeIds();

          for(int inode=0; inode<ele_adj->NumNode(); ++inode)
          {
            const int nodegid = nodeids[inode];

            //insert the node
            incompnodeids_set.insert(nodegid);

            //get the dofs of the node
            const std::vector<int> mydofs(bgdis->Dof(bgdis->gNode(nodegid)));

            // insert the dofs
            for (size_t k=0; k<mydofs.size(); ++k)
              incompdofs_set.insert(mydofs.at(k));

            // insert just the velocity dofs
            int numsets = mydofs.size()/4;
            incompdofs_vel_set.insert(mydofs.at((numsets-1)+0));
            incompdofs_vel_set.insert(mydofs.at((numsets-1)+1));
            incompdofs_vel_set.insert(mydofs.at((numsets-1)+2));

          }
          // break the element loop, we've already found one element for this node
          break;
        }
      }
    }
  }

  ////////////////////////////////////
  //check it again..
  for(set<int>::iterator iter = projectednodeids_.begin(); iter!= projectednodeids_.end();
      iter++)
  {
    set<int>::const_iterator iterpatchnodes = incompnodeids_set.find(*iter);

    if (iterpatchnodes == incompnodeids_set.end())
      dserror("BUG!! Nodes found which were in projected set but not in incompressibility patch!",*iter);
  }

// debug output
//   for(set<int>::iterator iter = incompelementids_set.begin(); iter!= incompelementids_set.end();
//       iter++)
//   {
//     set<int>::const_iterator iterpatchele = incompelementids_set.find(*iter);
//     cout << "all eles" << *iter << endl;
//     DRT::Node ** elenodes = bgdis->gElement(*iter)->Nodes();
//     for(int inode=0; inode<bgdis->gElement(*iter)->NumNode(); ++inode)
//     {
//       cout << bgdis->NumDof(elenodes[inode]) << endl;
//     }
//   }

  incompdofs_vel_without_dirichlet_set = incompdofs_vel_set;
  // throw the Dirichlet values out of incompdofs_vel map
  for(std::set<int>::iterator it=incompdofs_vel_set.begin(); it!=incompdofs_vel_set.end(); ++it)
  {
    if(dbcmaps->CondMap()->MyGID(*it))
      incompdofs_vel_without_dirichlet_set.erase(*it);
  }


  ///////////////////////////////////////
  // convert the sets to vectors
  std::copy(incompelementids_set.begin(), incompelementids_set.end(), std::back_inserter(incompelementids));
  std::copy(incompdofs_set.begin(), incompdofs_set.end(), std::back_inserter(incompdofs));
  std::copy(incompdofs_vel_without_dirichlet_set.begin(), incompdofs_vel_without_dirichlet_set.end(), std::back_inserter(incompdofs_vel));

  // Gather all informations from all processors
  vector<int> incompdofsAllproc;
  vector<int> incompveldofsAllproc;

  //information how many processors work at all
  vector<int> allproc(bgdis->Comm().NumProc());

  //in case of n processors allproc becomes a vector with entries (0,1,...,n-1)
  for (int i=0; i<bgdis->Comm().NumProc(); ++i) allproc[i] = i;

  //gathers information of all processors
  LINALG::Gather<int>(incompdofs,incompdofsAllproc,(int)bgdis->Comm().NumProc(),&allproc[0],bgdis->Comm());
  LINALG::Gather<int>(incompdofs_vel,incompveldofsAllproc,(int)bgdis->Comm().NumProc(),&allproc[0],bgdis->Comm());
  LINALG::Gather<int>(incompelementids,incompelementidsAllproc_,(int)bgdis->Comm().NumProc(),&allproc[0],bgdis->Comm());

  incompressibilitydofrowmap_ = rcp(new Epetra_Map(-1, incompdofsAllproc.size(), &incompdofsAllproc[0], 0, bgdis->Comm()));
  incompressibilityveldofrowmap_ = rcp(new Epetra_Map(-1, incompveldofsAllproc.size(), &incompveldofsAllproc[0], 0, bgdis->Comm()));

  // Gmsh debug output
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("incom_patch", step_, 5, 0, bgdis->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    {
      // draw bg elements with associated gid
      gmshfilecontent << "View \" " << "bg Element->Id() \" {\n";
      for (int i=0; i<bgdis->NumMyColElements(); ++i)
      {
        DRT::Element* actele = bgdis->lColElement(i);
//         GEO::CUT::ElementHandle * e = wizard_n.GetElement( actele );
        set<int>::const_iterator iter = incompelementids_set.find(actele->Id());

        if ( iter != incompelementids_set.end())
          IO::GMSH::elementAtInitialPositionToStream(1.0, actele, gmshfilecontent);
        else
          IO::GMSH::elementAtInitialPositionToStream(0.0, actele, gmshfilecontent);
      };
      gmshfilecontent << "};\n";
    }
  }

}
// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::EnforceIncompressibility(const RCP<DRT::Discretization>  bgdis,
                                                                XFEM::FluidWizard               wizard,
                                                                Teuchos::RCP<Epetra_Vector>     initialvel)
{
  // Map extractor of incompressibility patch
  LINALG::MapExtractor   incompmapextr = LINALG::MapExtractor(*bgdis->DofRowMap(),incompressibilityveldofrowmap_);


  //the original velocity vector which we want to improve
  Teuchos::RCP<Epetra_Vector> vel_org = LINALG::CreateVector(*incompressibilityveldofrowmap_,true);
  vel_org = incompmapextr.ExtractCondVector(initialvel);


  // Problem definition:
  //
  // Find u_p*, u_p: interpolated vector
  // min || u_p* - u_p ||^2
  //
  // s.t.  _
  //      |
  //      | div u_p* dOmega = 0
  //     _|
  //
  // equal to: cT.u_p* = 0


  ////////////////////////////////
  // Find the vector c:
  Teuchos::RCP<Epetra_Vector> C = LINALG::CreateVector(*incompressibilitydofrowmap_,true);
  Teuchos::RCP<Epetra_Vector> C_vel = LINALG::CreateVector(*incompressibilityveldofrowmap_,true);
  DRT::Element::LocationArray la( 1 );

  // loop over column elements of continuity patch
  const int numcolele = bgdis->NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* actele = bgdis->lColElement(i);

    // check if it is a conti patch element
    bool patchelement = true;
    size_t count = 0;
    for ( vector<int>::iterator it=incompelementidsAllproc_.begin(); it != incompelementidsAllproc_.end(); it++ )
    {
      count++;
      // if patch element
      if (*it == actele->Id())
        break;
      else if(count == incompelementidsAllproc_.size())
        patchelement = false;
    }


    if (patchelement == false) continue;

    Teuchos::RCP<MAT::Material> mat = actele->Material();
    DRT::ELEMENTS::Fluid * ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actele );

    GEO::CUT::ElementHandle * e = wizard.GetElement( actele );
    Epetra_SerialDenseVector C_elevec;
    // xfem element
    if ( e!=NULL )
    {
#ifdef DOFSETS_NEW
      std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
      std::vector< std::vector<int> > nds_sets;
      std::vector<std::vector< DRT::UTILS::GaussIntegration > >intpoints_sets;
      std::string VolumeCellGaussPointBy =  params_.sublist("XFEM").get<string>("VOLUME_GAUSS_POINTS_BY");

      e->GetCellSets_DofSets_GaussPoints( cell_sets, nds_sets, intpoints_sets, VolumeCellGaussPointBy );

      if(cell_sets.size() != intpoints_sets.size()) dserror("number of cell_sets and intpoints_sets not equal!");
      if(cell_sets.size() != nds_sets.size()) dserror("number of cell_sets and nds_sets not equal!");

      int set_counter = 0;

      for( std::vector< GEO::CUT::plain_volumecell_set>::iterator s=cell_sets.begin();
           s!=cell_sets.end();
           s++)
      {
        //       GEO::CUT::plain_volumecell_set & cells = *s;
        const std::vector<int> & nds = nds_sets[set_counter];

        // get element location vector, dirichlet flags and ownerships
        actele->LocationVector(*bgdis,nds,la,false);

        // number of dofs for background element
        // ndof contains all dofs of velocity and pressure. (we need just velocity)
        const size_t ndof   = la[0].lm_.size();
        C_elevec.Reshape(ndof,1);


        for( unsigned cellcount=0;cellcount!=cell_sets[set_counter].size();cellcount++ )
        {
          // call element method
          DRT::ELEMENTS::FluidFactory::ProvideImpl(actele->Shape(), "xfem")->CalculateContinuityXFEM(ele,
                                                                                           *bgdis,
                                                                                           la[0].lm_,
                                                                                           C_elevec,
                                                                                           intpoints_sets[set_counter][cellcount]);
        }
        set_counter++; //Shadan: make sure this is to be uncommented
      }
#else
      GEO::CUT::plain_volumecell_set cells;
      std::vector<DRT::UTILS::GaussIntegration> intpoints;
      std::vector<std::vector<double> > refEqns;
      std::string VolumeCellGaussPointBy =  params_.sublist("XFEM").get<string>("VOLUME_GAUSS_POINTS_BY");
      e->VolumeCellGaussPoints( cells, intpoints, refEqns, VolumeCellGaussPointBy);//modify gauss type

      int count = 0;
      for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        GEO::CUT::VolumeCell * vc = *i;
        if ( vc->Position()==GEO::CUT::Point::outside )
        {
          // one set of dofsets
          std::vector<int>  ndstest;
          for (int t=0;t<8; ++t)
            ndstest.push_back(0);

          actele->LocationVector(*bgdis,ndstest,la,false);

          const size_t ndof   = la[0].lm_.size();   // number of dofs for background element
          C_elevec.Reshape(ndof,1);


          // call element method
          DRT::ELEMENTS::FluidFactory::ProvideImpl(actele->Shape(), "xfem")->CalculateContinuityXFEM(ele,
                                                                                             *bgdis,
                                                                                             la[0].lm_,
                                                                                             C_elevec,
                                                                                             intpoints[count]);
        }
      }

      // if we have only one volume cell and it is inside (at void domain) we
      // need to handle it like a not-xfem element
      if ( cells.size()==1 )
      {
        GEO::CUT::VolumeCell * vc = *cells.begin();
        if ( vc->Position()==GEO::CUT::Point::inside )
        {
          // one set of dofsets
          std::vector<int>  ndstest;
          for (int t=0;t<8; ++t)
            ndstest.push_back(0);

          actele->LocationVector(*bgdis,ndstest,la,false);

          const size_t ndof   = la[0].lm_.size();   // number of dofs for background element
          C_elevec.Reshape(ndof,1);


          // call element method
          DRT::ELEMENTS::FluidFactory::ProvideImpl(actele->Shape(), "xfem")->CalculateContinuityXFEM(ele,
                                                                                             *bgdis,
                                                                                             la[0].lm_,
                                                                                             C_elevec);
        }
      }

#endif
    }
    // no xfem element
    else
    {
      // get element location vector, dirichlet flags and ownerships
      actele->LocationVector(*bgdis,la,false);

      const size_t ndof   = la[0].lm_.size();
      C_elevec.Reshape(ndof,1);

      DRT::ELEMENTS::FluidFactory::ProvideImpl(actele->Shape(), "xfem")->CalculateContinuityXFEM(ele,
                                                                                         *bgdis,
                                                                                         la[0].lm_,
                                                                                         C_elevec);


    }
    LINALG::Assemble(*C, C_elevec, la[0].lm_, la[0].lmowner_);
  }

  // The vector C still includes the pressure degrees of freedom which are
  // zero. So we need to eliminate them from C so that the vector C has just
  // the velocity degrees of freedom -> C_vel

  // C_vel should on all processors then this will work in parallel
  for (int iter=0; iter<C->MyLength();++iter)
  {
     int Cgid = C->Map().GID(iter);
     if (C_vel->Map().MyGID(Cgid))
       (*C_vel)[C_vel->Map().LID(Cgid)] = (*C)[C->Map().LID(Cgid)];
  }


  ////////////////////////////////////////////////////////////////
  // Transform the constraint cT.u* = 0 to (Qc)T.Qu* = 0
  // First with appropiate Givens Rotationen we annul all of the
  // entires of c but the first entry and ... Then we do the same
  // for all entries of c besides c_2 and so on.
  // next := entry to eliminate, next_us := entry to make next to zero
  // The vector we find at last is Qc
  //

  int next = 0;
  int next_us = 0;
  double TOL = 1.0e-16;

  // Initialize Q as identity matrix
  Teuchos::RCP<Epetra_CrsMatrix> Q = rcp(new Epetra_CrsMatrix(Copy,*incompressibilityveldofrowmap_,C_vel->MyLength()));

  double ones = 1.0;
  for(int i=0; i<C_vel->MyLength(); ++i)
  {
    int myGID = incompressibilityveldofrowmap_->GID(i);
    Q->InsertGlobalValues(myGID,1,&ones,&myGID);
  }

  Q->FillComplete();
  Teuchos::RCP<LINALG::SparseMatrix> Q_spr = rcp(new LINALG::SparseMatrix(Q));
  Q_spr->Complete();


  int mylengthminus = C_vel->MyLength()-1;
  while (next < mylengthminus)
  {
    // find the next entry which has to be eliminated and is not zero
    while (-TOL<(*C_vel)[next] && (*C_vel)[next]<TOL)
      next++;

    next_us = next+1;
    while (-TOL<(*C_vel)[next_us] && (*C_vel)[next_us]<TOL)
      next_us++;

    if(next_us < C_vel->MyLength())
    {
      double Nenner = sqrt((*C_vel)[next]*(*C_vel)[next]+(*C_vel)[next_us]*(*C_vel)[next_us]);
      Nenner = 1/Nenner;

      // c = x_i/(sqrt(x_i^2+x_j^2))
      // s = x_j/(sqrt(x_i^2+x_j^2))
      double cphi = Nenner*(*C_vel)[next];
      double sphi = Nenner*(*C_vel)[next_us];


      // Update of C_vel (which is after every update Qc)
      /* _     _
        | s  -c || x_i |    | s*x_i-c*x_j |   |      0      |
        |       ||     |  = |             | = |             |
        | c   s || x_j |    | c*x_i+s*x_j |   | c*x_i+s*x_j |
         -     -
        x_i: next
        x_j: next_us
      */

      (*C_vel)[next]  = 0.0;
      (*C_vel)[next_us] = cphi*(*C_vel)[next] + sphi*(*C_vel)[next_us];

      // Update Q = Qn..Q3*Q2*Q1

      // build the rotation matrix for current next and next_us
      Teuchos::RCP<Epetra_CrsMatrix> Q_i = rcp(new Epetra_CrsMatrix(Copy,*incompressibilityveldofrowmap_,2));

      // first build Q_i as an identity matrix
      double myval = 1.0;
      for(int j=0; j<C_vel->MyLength(); ++j)
      {
        int myGID = incompressibilityveldofrowmap_->GID(j);
        Q_i->InsertGlobalValues(myGID,1,&myval,&myGID);
      }

      int nextGID = incompressibilityveldofrowmap_->GID(next);
      int next_usGID = incompressibilityveldofrowmap_->GID(next_us);

      // set the rotation values of Q_i
      double cphi_min = -cphi;
      Q_i->ReplaceGlobalValues(nextGID,1,&sphi,&nextGID);
      Q_i->InsertGlobalValues(nextGID,1,&cphi_min,&next_usGID);
      Q_i->InsertGlobalValues(next_usGID,1,&cphi,&nextGID);
      Q_i->ReplaceGlobalValues(next_usGID,1,&sphi,&next_usGID);

      Q_i->FillComplete(*incompressibilityveldofrowmap_,*incompressibilityveldofrowmap_);

      // build a sparse matrix out of Q_i
      Teuchos::RCP<LINALG::SparseMatrix> Q_i_spr = rcp(new LINALG::SparseMatrix(Q_i));
      Q_i_spr->Complete(*incompressibilityveldofrowmap_,*incompressibilityveldofrowmap_);

      // build the final Q  = Qn..Q3*Q2*Q1
      Teuchos::RCP<LINALG::SparseMatrix> Q_final = rcp(new LINALG::SparseMatrix(*incompressibilityveldofrowmap_,C_vel->MyLength(),false,true));
      Q_final = MLMultiply(*Q_i_spr,false,*Q_spr,false,false,false);
      // save this old Q_final
      Q_spr = Q_final;
    }
    else
      break;
  }

  //string sysmat = "/home/shahmiri/work_tmp/sysmat";
  //Teuchos::RCP<LINALG::SparseMatrix> sysmatmatrixmatlab = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix);
  //LINALG::PrintMatrixInMatlabFormat(sysmat,*Q_spr->EpetraMatrix(),true);


//   int maxnumberofentries = 1000;
//   int numofnonzeros = 1000;
//   vector<double> Q_Zeile(maxnumberofentries, 0.0);
//   vector<int> indices(maxnumberofentries, 0);
//   (Q_spr->EpetraMatrix())->ExtractMyRowCopy(next,1000,numofnonzeros,&Q_Zeile[0]);


  ////////////////////////////////////////////////////////////////
  // find u* from Qu* != Qu
  // Note: From (Qc)T.Qu* = 0  we know that Qu*(next) is zero
  //
  // | Qu1 |   | Qu1    |       | Qu1 |   | Qu1    |   |    0   |
  // | Qu2 |   | Qu2    |       | Qu2 |   | Qu2    |   |    0   |
  // |  .  | = |  .     |   <=> |  .  | = |  .     | - |    0   |
  // |  .  |   |  .     |       |  .  |   |  .     |   |    0   |
  // |  0  |   | Qu_next|       |  0  |   | Qu_next|   | Qu_next|
  //
  // =>
  //               |  0     |       |x ...   next1||    0   |
  //               |  0     |       |  .     next2||    0   |
  // u* = Q'Qu - Q'|  0     | = u - |   .     .   ||    0   |
  //               |  0     |       |         .   ||    0   |
  //               | Qu_next|       |        nextn|| Qu_next|
  //

  // Qu
  Teuchos::RCP<Epetra_Vector> Qunext = LINALG::CreateVector(*incompressibilityveldofrowmap_,true);
  (Q_spr->EpetraMatrix())->Multiply(false,*vel_org,*Qunext);

  // we need the last entry of Qu (Qu_next)
  for (int i=0; i<Qunext->MyLength();++i)
  {
    int mygid = Qunext->Map().GID(i);
    int nextgid = Qunext->Map().GID(next);
    if ( mygid != nextgid)
      (*Qunext)[Qunext->Map().LID(mygid)] = 0.0;
  }

  Teuchos::RCP<Epetra_Vector> QTQu = LINALG::CreateVector(*incompressibilityveldofrowmap_,true);
  (Q_spr->EpetraMatrix())->Multiply(true,*Qunext,*QTQu);

  Teuchos::RCP<Epetra_Vector> u_incomp = LINALG::CreateVector(*incompressibilityveldofrowmap_,true);
  u_incomp->Update(1.0, *vel_org, -1.0, *QTQu, 0.0);
}

// // -------------------------------------------------------------------
// //
// // -------------------------------------------------------------------
// void FLD::XFluidFluid::CreatePatchBoxes(std::map<int, GEO::CUT::BoundingBox> & patchboxes)
// {
//   // get column version of the displacement vector
//   Teuchos::RCP<const Epetra_Vector> col_embfluiddisp =
//     DRT::UTILS::GetColVersionOfRowVector(embdis_, aledispnm_);

//   // Map of all boxes of embedded fluid discretization
//   for (int pele=0; pele<embdis_->NumMyColElements(); ++pele)
//   {
//     const DRT::Element* actpele = embdis_->lColElement(pele);
//     const DRT::Node*const* pelenodes = actpele->Nodes();

//     vector<int> lm;
//     vector<int> lmowner;
//     vector<int> lmstride;
//     actpele->LocationVector(*embdis_, lm, lmowner, lmstride);

//     vector<double> mydisp(lm.size());
//     DRT::UTILS::ExtractMyValues(*col_embfluiddisp, mydisp, lm);

//     GEO::CUT::BoundingBox patchbox;
//     //patchboxes
//     for (int pnode = 0; pnode < actpele->NumNode(); ++pnode)
//     {
//       // the coordinates of the actuall node
//       LINALG::Matrix<3,1> pnodepos(true);
//       pnodepos(0,0) = pelenodes[pnode]->X()[0] + mydisp[0+(pnode*4)];
//       pnodepos(1,0) = pelenodes[pnode]->X()[1] + mydisp[1+(pnode*4)];
//       pnodepos(2,0) = pelenodes[pnode]->X()[2] + mydisp[2+(pnode*4)];

//       // fill the patchbox
//       patchbox.AddPoint(pnodepos);
//     }
// //     cout << actpele->Id() << endl;
// //     patchbox.Print();
//     patchboxes[actpele->Id()] = patchbox;
//   }
// }


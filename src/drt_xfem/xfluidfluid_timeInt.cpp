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
#include <EpetraExt_MatrixMatrix.h>
#include "xfluidfluid_timeInt.H"
#include "xfem_fluidwizard.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dofset_transparent_independent.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_pstream.H"
#include "../drt_cut/cut_boundingbox.H"
#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_position.H"
#include "../drt_cut/cut_point.H"
#include "../drt_cut/cut_element.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_cut.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_interface.H"
#include "../drt_fluid_ele/fluid_ele_factory.H"
#include "../drt_fluid/fluid_utils.H"
#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"

#include <iostream>

XFEM::XFluidFluidTimeIntegration::XFluidFluidTimeIntegration(
  const Teuchos::RCP<DRT::Discretization> bgdis,
  const Teuchos::RCP<DRT::Discretization> embdis,
  Teuchos::RCP<XFEM::FluidWizardMesh>         wizard,
  int                            step,
  enum INPAR::XFEM::XFluidFluidTimeInt xfem_timeintapproach,
  const Teuchos::ParameterList&              params
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

  Teuchos::ParameterList&   params_xfem  = params_.sublist("XFEM");
  gmsh_debug_out_ = (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_DEBUG_OUT");

  Teuchos::ParameterList&   params_xf_gen = params_.sublist("XFLUID DYNAMIC/GENERAL");
  searchradius_fac_= params_xf_gen.get<double>("XFLUIDFLUID_SEARCHRADIUS");

  // find the radius of the search tree
  SearchRadius();

  return;

} // end constructor

// -------------------------------------------------------------------
// map of standard node ids and their dof-gids in for this time step
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::CreateBgNodeMaps(const Teuchos::RCP<DRT::Discretization> bgdis,
                                                        Teuchos::RCP<XFEM::FluidWizardMesh>         wizard)
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
        e->GetVolumeCellsDofSets( cell_sets, nds_sets, false); //(include_inner=false)

        parenteletondsset_[actele->Id()] = nds_sets;

      }

      nodetoparentele_[gid] = parentelements;


      GEO::CUT::Point * p = n->point();
      GEO::CUT::Point::PointPosition pos = p->Position();
      if (pos==GEO::CUT::Point::outside and bgdis->NumDof(node) != 0) //std
      {
        //std::cout << node->Id() << "std!! size()" << vcs.size() << std::endl;
        //std::cout << " outside " << pos <<  " "<< node->Id() << std::endl;
        std::vector<int> gdofs = bgdis->Dof(node);
        stdnodenp_[gid] = gdofs;
      }
      else if (pos==GEO::CUT::Point::inside and  bgdis->NumDof(node) == 0) //void
      {
        //std::cout << " inside " <<  pos << " " << node->Id() << std::endl;
      }
      else if (pos==GEO::CUT::Point::inside and  bgdis->NumDof(node) != 0) //enriched
      {
        //std::cout << " inside enriched" <<  pos << " " << node->Id() << std::endl;
        std::vector<int> gdofs = bgdis->Dof(node);
        enrichednodenp_[gid] = gdofs;
      }
      else if (pos==GEO::CUT::Point::oncutsurface and  bgdis->NumDof(node) != 0)
      {
        std::vector<int> gdofs = bgdis->Dof(node);
        stdnodenp_[gid] = gdofs;
      }
      // this case only happens if the two fluid domains are both the same
      // size (void)
      else if (pos==GEO::CUT::Point::oncutsurface and  bgdis->NumDof(node) == 0)
      {
        //std::cout << " on the surface " <<  pos << " " << node->Id() << std::endl;
      }
      else
      {
        IO::cout << "  here ?! " <<  pos << " " <<  node->Id() <<  " numdof: " << bgdis->NumDof(node) << IO::endl;
      }
    }
    else if( bgdis->NumDof(node) != 0) // no xfem node
    {
      std::vector<int> gdofs = bgdis->Dof(node);
      stdnodenp_[gid] = gdofs;
    }
    else
      IO::cout << " why here? " << "node "<<  node->Id()<< " " <<  bgdis->NumDof(node) <<  IO::endl;
  }

#ifdef PARALLEL
  // build a reduced map of all noderowmaps of all processors
  Teuchos::RCP<Epetra_Map> allnoderowmap = LINALG::AllreduceEMap(*noderowmap);
  // gather the informations of all processors
  DRT::Exporter ex(*noderowmap,*allnoderowmap,bgdis->Comm());
  ex.Export(stdnodenp_);
  ex.Export(enrichednodenp_);
#endif

//  debug output
//   for(std::map<int, std::vector<int> >::iterator iter = nodetoparentele_.begin(); iter!= nodetoparentele_.end();
//       iter++)
//   {
//     std::cout << "ngid " <<  iter->first << std::endl;
//     std::vector<int> fff= iter->second;
//     for (int u=0;u<fff.size();++u)
//       std::cout << fff.at(u);
//     std::cout << std::endl;
//   }
//   for(std::map<int, std::vector< std::vector< int > > >::iterator iter = parenteletondsset_.begin(); iter!= parenteletondsset_.end();
//       iter++)
//   {
//     std::cout << "elegid " << iter->first << std::endl;
//     std::cout << "size nds " << iter->second.size() << std::endl;
//     for (int s=0; s < iter->second.size(); ++s)
//     {
//       const std::vector<int> & nds = iter->second[s];
//       for (int ttt=0; ttt<nds.size(); ++ttt)
//         std::cout << nds.at(ttt) << " " ;
//       std::cout << std::endl;
//     }
//   }
//     for (int i=0; i<bgdis->NumMyColNodes(); ++i)
//     {
//        const DRT::Node* actnode = bgdis->lColNode(i);
//        std::map<int, std::vector<int> >::const_iterator iter = stdnoden_.find(actnode->Id());
//        std::map<int, std::vector<int> >::const_iterator iter2 = enrichednoden_.find(actnode->Id());
//        std::map<int, std::vector<int> >::const_iterator iter3 = stdnodenp_.find(actnode->Id());
//        std::map<int, std::vector<int> >::const_iterator iter4 = enrichednodenp_.find(actnode->Id());
//        if (iter2 != enrichednoden_.end()) std::cout  << " enrichned n : " << actnode->Id() << " "  ;
//        if (iter2 == enrichednoden_.end() and iter == stdnoden_.end()) std::cout  << " void n :" <<  actnode->Id() << " "  ;
//        if (iter4 != enrichednodenp_.end()) std::cout  << " enrichned np : " << actnode->Id() << " "  ;
//        if (iter4 == enrichednodenp_.end() and iter3 == stdnodenp_.end()) std::cout  << " void np :" <<  actnode->Id() << " "  ;
//        if (iter  != stdnoden_.end()) std::cout  << " std n : " << actnode->Id() << " "  ;
//        if (iter3 != stdnodenp_.end()) std::cout  << " std np :" <<  actnode->Id() << " "  ;
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
int XFEM::XFluidFluidTimeIntegration::SaveAndCreateNewBgNodeMaps(Teuchos::RCP<DRT::Discretization> bgdis,
                                                                 Teuchos::RCP<XFEM::FluidWizardMesh>   wizard)
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

  if (gmsh_debug_out_)
    GmshOutput(bgdis);

  return samemaps_;
}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void  XFEM::XFluidFluidTimeIntegration::CreateBgNodeMapsForRestart(Teuchos::RCP<DRT::Discretization> bgdis,
                                                                   Teuchos::RCP<XFEM::FluidWizardMesh>   wizard)
{

  // Create new maps
  CreateBgNodeMaps(bgdis,wizard);

  currentbgdofmap_ =  bgdis->DofRowMap();

}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SetNewBgStatevectorAndProjectEmbToBg(const Teuchos::RCP<DRT::Discretization>        bgdis,
                                                                            Teuchos::RCP<Epetra_Vector>           bgstatevn1,
                                                                            Teuchos::RCP<Epetra_Vector>           bgstatevnp1,
                                                                            Teuchos::RCP<Epetra_Vector>           embstatevn1,
                                                                            Teuchos::RCP<Epetra_Vector>           bgstatevn2,
                                                                            Teuchos::RCP<Epetra_Vector>           bgstatevnp2,
                                                                            Teuchos::RCP<Epetra_Vector>           embstatevn2,
                                                                            Teuchos::RCP<Epetra_Vector>           bgstatevn3,
                                                                            Teuchos::RCP<Epetra_Vector>           bgstatevnp3,
                                                                            Teuchos::RCP<Epetra_Vector>           embstatevn3,
                                                                            Teuchos::RCP<Epetra_Vector>           bgstatevn4,
                                                                            Teuchos::RCP<Epetra_Vector>           bgstatevnp4,
                                                                            Teuchos::RCP<Epetra_Vector>           embstatevn4,
                                                                            Teuchos::RCP<Epetra_Vector>           bgstatevn5,
                                                                            Teuchos::RCP<Epetra_Vector>           bgstatevnp5,
                                                                            Teuchos::RCP<Epetra_Vector>           embstatevn5,
                                                                            Teuchos::RCP<Epetra_Vector>           aledispn)
{
  if (timeintapproach_ == INPAR::XFEM::Xff_TimeInt_FullProj or
      timeintapproach_ == INPAR::XFEM::Xff_TimeInt_ProjIfMoved or
      timeintapproach_ == INPAR::XFEM::Xff_TimeInt_IncompProj)
    SetNewBgStatevectorFullProjection(bgdis,
                              bgstatevn1, bgstatevnp1, embstatevn1,
                              bgstatevn2, bgstatevnp2, embstatevn2,
                              bgstatevn3, bgstatevnp3, embstatevn3,
                              bgstatevn4, bgstatevnp4, embstatevn4,
                              bgstatevn5, bgstatevnp5, embstatevn5,
                              aledispn);
  else if(timeintapproach_ == INPAR::XFEM::Xff_TimeInt_KeepGhostValues)
    SetNewBgStatevectorKeepGhostValues(bgdis,
                               bgstatevn1, bgstatevnp1, embstatevn1,
                               bgstatevn2, bgstatevnp2, embstatevn2,
                               bgstatevn3, bgstatevnp3, embstatevn3,
                               bgstatevn4, bgstatevnp4, embstatevn4,
                               bgstatevn5, bgstatevnp5, embstatevn5,
                               aledispn);
  else
    dserror("xfem time integration approach unknown!");
}

// -------------------------------------------------------------------
// Always do the projection from embedded fluid. Also for the enriched
// nodes.
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SetNewBgStatevectorFullProjection(const Teuchos::RCP<DRT::Discretization>       bgdis,
                                                             Teuchos::RCP<Epetra_Vector>          bgstatevn1,
                                                             Teuchos::RCP<Epetra_Vector>          bgstatevnp1,
                                                             Teuchos::RCP<Epetra_Vector>          embstatevn1,
                                                             Teuchos::RCP<Epetra_Vector>          bgstatevn2,
                                                             Teuchos::RCP<Epetra_Vector>          bgstatevnp2,
                                                             Teuchos::RCP<Epetra_Vector>          embstatevn2,
                                                             Teuchos::RCP<Epetra_Vector>          bgstatevn3,
                                                             Teuchos::RCP<Epetra_Vector>          bgstatevnp3,
                                                             Teuchos::RCP<Epetra_Vector>          embstatevn3,
                                                             Teuchos::RCP<Epetra_Vector>          bgstatevn4,
                                                             Teuchos::RCP<Epetra_Vector>          bgstatevnp4,
                                                             Teuchos::RCP<Epetra_Vector>          embstatevn4,
                                                             Teuchos::RCP<Epetra_Vector>          bgstatevn5,
                                                             Teuchos::RCP<Epetra_Vector>          bgstatevnp5,
                                                             Teuchos::RCP<Epetra_Vector>          embstatevn5,
                                                                         Teuchos::RCP<Epetra_Vector>          aledispn)
{

  // coordinates of bg-nodes which need the projection from embedded dis
  std::vector<LINALG::Matrix<3,1> > bgnodes_coords;
  // the vector containing interpolated values for each state vector from embedded dis. For each node one interpolated value..
  // Entries 0 to 3: velnp_, Entries 4 to 7: veln_,  Entries 8 to 11: velnm_, Entries 12 to 15: accn_, Entries 16 to 19: accnp_
  std::vector<LINALG::Matrix<20,1> > interpolated_vecs;
  // bg-node ids which have no history and need a projection
  std::vector<int> bgnodeidwithnohistory;

  projectednodeids_.clear();

  // loop over bg-row-nodes of each processor
  for (int lnid=0; lnid<bgdis->NumMyRowNodes(); lnid++)
  {
    DRT::Node* bgnode = bgdis->lRowNode(lnid);
    std::map<int, std::vector<int> >::const_iterator iterstn = stdnoden_.find(bgnode->Id());
    std::map<int, std::vector<int> >::const_iterator iterstnp = stdnodenp_.find(bgnode->Id());
    std::map<int, std::vector<int> >::const_iterator iteren = enrichednoden_.find(bgnode->Id());
    std::map<int, std::vector<int> >::const_iterator iterenp = enrichednodenp_.find(bgnode->Id());

    // Transfer the dofs:
    // n:std -> n+1:std, n:std -> n+1:enriched
    if ((iterstn != stdnoden_.end() and iterstnp != stdnodenp_.end()) or
        (iterstn != stdnoden_.end() and iterenp != enrichednodenp_.end()))
    {
      std::vector<int> gdofsn = iterstn->second;

#ifdef DEBUG
      if (iterstn->second.size()>4)
        IO::cout << "INFO: more standard sets!!!!"<< "Node GID " << bgnode->Id() <<" size " <<
          iterstn->second.size()  << IO::endl;
#endif

//      int numsets = bgdis->NumDof(bgnode)/4;
//      if (numsets > 1)
//         std::cout << GREEN_LIGHT << "Info: more dofsets in transfer.. " <<  "Node GID " << bgnode->Id() << END_COLOR << std::endl;

      //  the first set of bgstatevn is transfered to bgstatevnp.
      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp1,bgstatevn1);
      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp2,bgstatevn2);
      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp3,bgstatevn3);
      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp4,bgstatevn4);
      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp5,bgstatevn5);

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

      // the vector of interpolated values of all state vector (max five state vectors) for every node
      LINALG::Matrix<20,1>    interpolatedvec(true);
      interpolated_vecs.push_back(interpolatedvec);

    }
    //do nothing:
    //n: void->n+1: void, n:std->n+1:void, n:enriched->n+1:void
    else if ( (iterstn == stdnoden_.end() and iteren == enrichednoden_.end()
               and iterstnp == stdnodenp_.end() and iterenp == enrichednodenp_.end() ) or
              (iterstn != stdnoden_.end() and (iterstnp == stdnodenp_.end() and iterenp == enrichednodenp_.end())) or
              (iteren != enrichednoden_.end() and (iterenp == enrichednodenp_.end() and iterstnp == stdnodenp_.end())) )
    {

      if( bgdis->NumDof(bgnode) > 0 )
        IO::cout << "BUG:: in do nothig!!" << " Node GID " << bgnode->Id() << IO::endl;
    }
    else
      IO::cout << "warum bin ich da?! " <<  bgdis->NumDof(bgnode)   << " " <<    bgnode->Id() <<  IO::endl;
  }


  // call the Round Robin Communicator
  CommunicateNodes(bgdis,bgnodes_coords,interpolated_vecs,
               bgnodeidwithnohistory,aledispn,
               embstatevn1,bgstatevnp1,bgstatevn1,
               embstatevn2,bgstatevnp2,bgstatevn2,
               embstatevn3,bgstatevnp3,bgstatevn3,
               embstatevn4,bgstatevnp4,bgstatevn4,
               embstatevn5,bgstatevnp5,bgstatevn5);

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
void XFEM::XFluidFluidTimeIntegration::CommunicateNodes(const Teuchos::RCP<DRT::Discretization>        bgdis,
                                                        std::vector<LINALG::Matrix<3,1> >   & bgnodes_coords,
                                                        std::vector<LINALG::Matrix<20,1> >  & interpolated_vecs,
                                                        std::vector<int>                    & bgnodeidwithnohistory,
                                                        Teuchos::RCP<Epetra_Vector>           aledispn,
                                                        Teuchos::RCP<Epetra_Vector>           embstatevn1,
                                                        Teuchos::RCP<Epetra_Vector>           bgstatevnp1,
                                                        Teuchos::RCP<Epetra_Vector>           bgstatevn1,
                                                        Teuchos::RCP<Epetra_Vector>           embstatevn2,
                                                        Teuchos::RCP<Epetra_Vector>           bgstatevnp2,
                                                        Teuchos::RCP<Epetra_Vector>           bgstatevn2,
                                                        Teuchos::RCP<Epetra_Vector>           embstatevn3,
                                                        Teuchos::RCP<Epetra_Vector>           bgstatevnp3,
                                                        Teuchos::RCP<Epetra_Vector>           bgstatevn3,
                                                        Teuchos::RCP<Epetra_Vector>           embstatevn4,
                                                        Teuchos::RCP<Epetra_Vector>           bgstatevnp4,
                                                        Teuchos::RCP<Epetra_Vector>           bgstatevn4,
                                                        Teuchos::RCP<Epetra_Vector>           embstatevn5,
                                                        Teuchos::RCP<Epetra_Vector>           bgstatevnp5,
                                                        Teuchos::RCP<Epetra_Vector>           bgstatevn5)
 {

  // get number of processors and the current processors id
  int numproc=embdis_->Comm().NumProc();

  //information how many processors work at all
  std::vector<int> allproc(numproc);

  // create an exporter for point to point comunication
  DRT::Exporter exporter(embdis_->Comm());

  // necessary variables
  MPI_Request request;

  // define send and receive blocks
  std::vector<char> sblock;
  std::vector<char> rblock;

  // vector which identifies if a bg-node has already interpolated values
  std::vector<int> NodeDone;
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
      std::vector<LINALG::Matrix<3,1> > stuff_coord;
      // we have state vectors, each has 4 dofs (20 interpolated values)
      std::vector<LINALG::Matrix<20,1> > stuff_interpolatedvecs;
      std::vector<int>  stuff_bgnodeidswithnohistory;
      std::vector<int>  stuff_nodedone;

      std::vector<char>::size_type position = 0;
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
      FindEmbeleAndInterpolatevalues(bgnodes_coords,interpolated_vecs,NodeDone,
                                 aledispn,embstatevn1,embstatevn2,embstatevn3,embstatevn4,embstatevn5);

      // Pack info into block to sendit
      PackValues(bgnodes_coords,interpolated_vecs,bgnodeidwithnohistory,NodeDone,sblock);

      // add size  to sendblock
      SendBlock(sblock,exporter,request);
    }
  } // end of loop over processors

  //----------------------------------------------------------------------------------------------
  //set the interpolated values to bgstatevnp
  //---------------------------------------------------------------------------------------------
  for (size_t i=0; i<bgnodeidwithnohistory.size(); ++i)
  {
    DRT::Node* bgnode = bgdis->gNode(bgnodeidwithnohistory.at(i));
    // number of dof-sets
    int numsets = bgdis->NumDof(bgnode)/4;

    int offset = 0;

    // if interpolated values are available
    if (NodeDone.at(i) == 1)
    {
      for (int set=0; set<numsets; set++)
      {
        // offset for different state vectors
        int countvecs = 0;
      for (std::size_t isd = 0; isd < 4; ++isd){
        (*bgstatevnp1)[bgstatevnp1->Map().LID(bgdis->Dof(bgnode)[offset+isd])] = interpolated_vecs.at(i)(isd+countvecs);
      }

      countvecs += 4;
      for (std::size_t isd = 0; isd < 4; ++isd){
        (*bgstatevnp2)[bgstatevnp2->Map().LID(bgdis->Dof(bgnode)[offset+isd])] = interpolated_vecs.at(i)(isd+countvecs);
      }

      countvecs += 4;
      for (std::size_t isd = 0; isd < 4; ++isd){
        (*bgstatevnp3)[bgstatevnp2->Map().LID(bgdis->Dof(bgnode)[offset+isd])] = interpolated_vecs.at(i)(isd+countvecs);
      }

      countvecs += 4;
      for (std::size_t isd = 0; isd < 4; ++isd){
        (*bgstatevnp4)[bgstatevnp2->Map().LID(bgdis->Dof(bgnode)[offset+isd])] = interpolated_vecs.at(i)(isd+countvecs);
      }

      countvecs += 4;
      for (std::size_t isd = 0; isd < 4; ++isd){
        (*bgstatevnp5)[bgstatevnp2->Map().LID(bgdis->Dof(bgnode)[offset+isd])] = interpolated_vecs.at(i)(isd+countvecs);
      }

        offset += 4;
      }
    }
    // if no embedded element is found try to find an enriched value
    else
    {
      std::map<int, std::vector<int> >::const_iterator iterstn = stdnoden_.find(bgnode->Id());
      std::map<int, std::vector<int> >::const_iterator iterstnp = stdnodenp_.find(bgnode->Id());
      std::map<int, std::vector<int> >::const_iterator iteren = enrichednoden_.find(bgnode->Id());
      std::map<int, std::vector<int> >::const_iterator iterenp = enrichednodenp_.find(bgnode->Id());

      if ((iteren != enrichednoden_.end() and iterenp != enrichednodenp_.end())
          or ((iteren != enrichednoden_.end() and iterstnp != stdnodenp_.end())))
      {
        IO::cout << "CHECK: Took enriched values !!" << " Node GID " << bgnode->Id() << IO::endl;
        IO::cout << " Warning: You may need to make your search radius bigger in the dat-file!" << IO::endl;
        std::vector<int> gdofsn = iteren->second;

        WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp1,bgstatevn1);
        WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp2,bgstatevn2);
        WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp3,bgstatevn3);
        WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp4,bgstatevn4);
        WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp5,bgstatevn5);
      }
      else
      {
        IO::cout << " Warning: No patch element found for the node " << bgnode->Id();
        if ((iterstn == stdnoden_.end() and iteren == enrichednoden_.end()) and iterstnp != stdnodenp_.end())
          IO::cout << " n:void -> n+1:std  " << IO::endl;
        else if (iteren != enrichednoden_.end() and iterenp != enrichednodenp_.end())
          IO::cout << " n:enriched -> n+1:enriched " << IO::endl;
        else if (iteren != enrichednoden_.end() and iterstnp != stdnodenp_.end())
          IO::cout << " n:enriched -> n+1: std " << IO::endl;
        else if ((iterstn == stdnoden_.end() and iteren == enrichednoden_.end()) and iterenp != enrichednodenp_.end())
          IO::cout << " n:void ->  n+1:enriched " << IO::endl;
        IO::cout << " Warning: You may need to make your search radius bigger in the dat-file!" << IO::endl;
      }
    }
  }// end of loop over bgnodes without history

 }//XFEM::XFluidFluidTimeIntegration::CommunicateNodes

//---------------------------------------------------------
// receive a block in the round robin communication pattern
//---------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::ReceiveBlock(std::vector<char>   & rblock,
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

#ifdef DEBUG
  IO::cout << "----receiving " << rblock.size() <<  " bytes: to proc " << myrank << " from proc " << frompid << IO::endl;
#endif

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
void XFEM::XFluidFluidTimeIntegration::SendBlock(std::vector<char>  & sblock  ,
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

#ifdef DEBUG
   IO::cout << "----sending " << sblock.size() <<  " bytes: from proc " << myrank << " to proc " << topid << IO::endl;
#endif

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
void XFEM::XFluidFluidTimeIntegration::PackValues(std::vector<LINALG::Matrix<3,1> >  & bgnodes_coords,
                                                  std::vector<LINALG::Matrix<20,1> > & interpolatedvec,
                                                  std::vector<int>                   & bgnodeidwithnohistory,
                                                  std::vector<int>                   & NodeDone,
                                                  std::vector<char>                  & sblock)
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
void XFEM::XFluidFluidTimeIntegration::FindEmbeleAndInterpolatevalues(std::vector<LINALG::Matrix<3,1> >  & bgnodes_coords,
                                                                      std::vector<LINALG::Matrix<20,1> > & interpolated_vecs,
                                                                      std::vector<int>                   & NodeDone,
                                                                      Teuchos::RCP<Epetra_Vector>          aledispn,
                                                                      Teuchos::RCP<Epetra_Vector>          embstatevn1,
                                                                      Teuchos::RCP<Epetra_Vector>          embstatevn2,
                                                                      Teuchos::RCP<Epetra_Vector>          embstatevn3,
                                                                      Teuchos::RCP<Epetra_Vector>          embstatevn4,
                                                                      Teuchos::RCP<Epetra_Vector>          embstatevn5)

{
  //init of 3D search tree
  Teuchos::RCP<GEO::SearchTree> searchTree = Teuchos::rcp(new GEO::SearchTree(5));

  // find current positions for emb fluid discretization
  std::map<int,LINALG::Matrix<3,1> > currentpositions;
  for (int lid = 0; lid < embdis_->NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = embdis_->lColNode(lid);
    LINALG::Matrix<3,1> currpos;
    std::vector<int> pgdofs(4);
    std::vector<double> mydisp(4);

    // get the current displacement
    embdis_->Dof(node,0,pgdofs);
    DRT::UTILS::ExtractMyValues(*aledispn,mydisp,pgdofs);

    currpos(0) = node->X()[0]+mydisp.at(0);
    currpos(1) = node->X()[1]+mydisp.at(1);
    currpos(2) = node->X()[2]+mydisp.at(2);

    currentpositions[node->Id()] = currpos;
  }

  // find the bounding box of the elements and initialize the search tree
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(*embdis_,currentpositions);
  searchTree->initializeTree(rootBox,*embdis_,GEO::TreeType(GEO::OCTTREE));

  // loop over bgnodes_coords
  for (size_t i=0; i<bgnodes_coords.size(); ++i)
  {
    // indicates that we found the embedded element, the background node is covered by
    bool insideelement = false;

    // bgnode coordinate
    LINALG::Matrix<3,1> bgnodecords = bgnodes_coords.at(i);
    // interpolated vector which is zero at the beginning
    LINALG::Matrix<20,1> interpolatedvec = interpolated_vecs.at(i);

    //search for near elements to the background node's coord
    std::map<int,std::set<int> >  closeeles =
        searchTree->searchElementsInRadius(*embdis_,currentpositions,bgnodecords,minradius_,0);


    // Remark: it could be that closeles is empty on one processor but still has elements on other processors.

    if(closeeles.empty() == false)
    {
      for(std::map<int, std::set<int> >::const_iterator closele = closeeles.begin(); closele != closeeles.end(); closele++)
      {
        if (insideelement) break;
        for(std::set<int>::const_iterator eleIter = (closele->second).begin(); eleIter != (closele->second).end(); eleIter++)
        {
          if (insideelement) break;
          DRT::Element* pele = embdis_->gElement(*eleIter); // eleIter is the gid of the pele
          insideelement = ComputeSpacialToElementCoordAndProject(pele,bgnodecords,interpolatedvec,
                                                             embstatevn1,embstatevn2,embstatevn3,embstatevn4,embstatevn5,
                                                               aledispn,embdis_);

          if (insideelement and NodeDone.at(i)==0)
            {
              // Set NodeDone to 1 if the embedded-element is found
              NodeDone.at(i) = 1;
              for (size_t j=0; j<20; ++j)
              {
                //set the interpolated values
                interpolated_vecs.at(i)(j) = interpolatedvec(j);
              }
            }
        }
      }
    }
    else IO::cout << "The search radius is empty on one processor! You may need to change the XFLUIDFLUID_SEARCHRADIUS is your dat-file."<< IO::endl;
  }
}

// -------------------------------------------------------------------
// If enriched values are available keep them (no projection from
// embedded fluid)
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SetNewBgStatevectorKeepGhostValues(const Teuchos::RCP<DRT::Discretization>        bgdis,
                                                                                Teuchos::RCP<Epetra_Vector>           bgstatevn1,
                                                                                Teuchos::RCP<Epetra_Vector>           bgstatevnp1,
                                                                                Teuchos::RCP<Epetra_Vector>           embstatevn1,
                                                                                Teuchos::RCP<Epetra_Vector>           bgstatevn2,
                                                                                Teuchos::RCP<Epetra_Vector>           bgstatevnp2,
                                                                                Teuchos::RCP<Epetra_Vector>           embstatevn2,
                                                                                Teuchos::RCP<Epetra_Vector>           bgstatevn3,
                                                                                Teuchos::RCP<Epetra_Vector>           bgstatevnp3,
                                                                                Teuchos::RCP<Epetra_Vector>           embstatevn3,
                                                                                Teuchos::RCP<Epetra_Vector>           bgstatevn4,
                                                                                Teuchos::RCP<Epetra_Vector>           bgstatevnp4,
                                                                                Teuchos::RCP<Epetra_Vector>           embstatevn4,
                                                                                Teuchos::RCP<Epetra_Vector>           bgstatevn5,
                                                                                Teuchos::RCP<Epetra_Vector>           bgstatevnp5,
                                                                                Teuchos::RCP<Epetra_Vector>           embstatevn5,
                                                                                Teuchos::RCP<Epetra_Vector>           aledispn)
{

  // coordinates of bg-nodes which need the projection from embedded dis
  std::vector<LINALG::Matrix<3,1> > bgnodes_coords;
  // the vector containing interpolated values for each state vector from embedded dis. For each node one interpolated value..
  // Entries 0 to 3: velnp_, Entries 4 to 7: veln_,  Entries 8 to 11: velnm_, Entries 12 to 15: accn_, Entries 16 to 19: accnp_
  std::vector<LINALG::Matrix<20,1> > interpolated_vecs;
  // bg-node ids which have no history and need a projection
  std::vector<int> bgnodeidwithnohistory;

  // loop over bg-row-nodes of each processor
  for (int lnid=0; lnid<bgdis->NumMyRowNodes(); lnid++)
  {
    DRT::Node* bgnode = bgdis->lRowNode(lnid);
    std::map<int, std::vector<int> >::const_iterator iterstn = stdnoden_.find(bgnode->Id());
    std::map<int, std::vector<int> >::const_iterator iterstnp = stdnodenp_.find(bgnode->Id());
    std::map<int, std::vector<int> >::const_iterator iteren = enrichednoden_.find(bgnode->Id());
    std::map<int, std::vector<int> >::const_iterator iterenp = enrichednodenp_.find(bgnode->Id());

    // Transfer the dofs:
    // n:std -> n+1:std, n:std -> n+1:enriched
    if ((iterstn != stdnoden_.end() and iterstnp != stdnodenp_.end()) or
        (iterstn != stdnoden_.end() and iterenp != enrichednodenp_.end()))
    {
      //int numsets = bgdis->NumDof(bgnode)/4;
      std::vector<int> gdofsn = iterstn->second;

      //TODO!! die richtige dofs von bgstatevn rauspicke, wenn mehrere
      //dofsets vorhanden sind

      // Information
#ifdef DEBUG
      if (iterstn->second.size()>4)
        IO::cout << " INFO: more standard sets!!!! "<< "Node GID " << bgnode->Id() << IO::endl;
#endif

      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp1,bgstatevn1);
      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp2,bgstatevn2);
      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp3,bgstatevn3);
      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp4,bgstatevn4);
      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp5,bgstatevn5);

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
      LINALG::Matrix<20,1>    interpolatedvec(true);
      interpolated_vecs.push_back(interpolatedvec);

      //int numsets = bgdis->NumDof(bgnode)/4;

//      if( numsets > 1 )
//        std::cout << GREEN_LIGHT << "Info: more dofsets in projection.. " <<  "Node GID " << bgnode->Id() << END_COLOR << std::endl;

    }
    //keep the ghost dofs:
    //n: enriched -> n+1:enriched, n: enriched -> n+1: std
    else if ((iteren != enrichednoden_.end() and iterenp != enrichednodenp_.end()) or
             (iteren != enrichednoden_.end() and iterstnp != stdnodenp_.end()))
    {
      int numsets = bgdis->NumDof(bgnode)/4;
      if (numsets > 1)
        IO::cout << "ghost-fluid-approach just available for one dofset!" << IO::endl;

      std::vector<int> gdofsn = iteren->second;
      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp1,bgstatevn1);
      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp2,bgstatevn2);
      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp3,bgstatevn3);
      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp4,bgstatevn4);
      WriteValuestoBgStateVector(bgdis,bgnode,gdofsn,bgstatevnp5,bgstatevn5);

    }
    //do nothing:
    //n: void->n+1: void, n:std->n+1:void, n:enriched->n+1:void
    else if ( (iterstn == stdnoden_.end() and iteren == enrichednoden_.end()
               and iterstnp == stdnodenp_.end() and iterenp == enrichednodenp_.end() ) or
              (iterstn != stdnoden_.end() and (iterstnp == stdnodenp_.end() and iterenp == enrichednodenp_.end())) or
              (iteren != enrichednoden_.end() and (iterenp == enrichednodenp_.end() and iterstnp == stdnodenp_.end())) )
    {

      if( bgdis->NumDof(bgnode) > 0 )
        IO::cout  <<  "BUG:: in do nothig!!" << " Node GID " << bgnode->Id() <<  IO::endl;
      //std::cout << "do nothing" << bgnode->Id() << std::endl; ;
    }
    else
      IO::cout << "warum bin ich da?! " <<  bgdis->NumDof(bgnode)   << " " <<    bgnode->Id() <<  IO::endl;
  }

  // call the Round Robin Communicator

  CommunicateNodes(bgdis,
               bgnodes_coords,interpolated_vecs,bgnodeidwithnohistory,
               aledispn,
               embstatevn1,bgstatevnp1,bgstatevn1,
               embstatevn2,bgstatevnp2,bgstatevn2,
               embstatevn3,bgstatevnp3,bgstatevn3,
               embstatevn4,bgstatevnp4,bgstatevn4,
               embstatevn5,bgstatevnp5,bgstatevn5);

}//XFEM::XFluidFluidTimeIntegration::SetNewBgStatevectorKeepGhostValues

//-------------------------------------------------------------------
// Write the values of node from bgstatevn to bgstatevnp
//--------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::WriteValuestoBgStateVector(const Teuchos::RCP<DRT::Discretization>   bgdis,
                                                                  DRT::Node*                       bgnode,
                                                                  std::vector<int>                 gdofs_n,
                                                                  Teuchos::RCP<Epetra_Vector>      bgstatevnp,
                                                                  Teuchos::RCP<Epetra_Vector>      bgstatevn)
{
  int numsets = bgdis->NumDof(bgnode)/4;

#ifdef DEBUG
  if (numsets > 1)
    IO::cout << "Info: more dofsets in transfer.. " <<  "Node GID " << bgnode->Id() << IO::endl;
#endif

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
                                                                              LINALG::Matrix<20,1>&               interpolatedvec,
                                                                              Teuchos::RCP<Epetra_Vector>         embstate_n1,
                                                                              Teuchos::RCP<Epetra_Vector>         embstate_n2,
                                                                              Teuchos::RCP<Epetra_Vector>         embstate_n3,
                                                                              Teuchos::RCP<Epetra_Vector>         embstate_n4,
                                                                              Teuchos::RCP<Epetra_Vector>         embstate_n5,
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
      LINALG::Matrix<4,numnodes> veln1(true);
      LINALG::Matrix<4,numnodes> veln2(true);
      LINALG::Matrix<4,numnodes> veln3(true);
      LINALG::Matrix<4,numnodes> veln4(true);
      LINALG::Matrix<4,numnodes> veln5(true);

      LINALG::Matrix<4,numnodes> disp;

      // loop over the nodes of the patch element and save the veln and
      // coordinates of all nodes
      for (int inode = 0; inode < numnodes; ++inode)
      {
        sourcedis->Dof(pelenodes[inode],0,pgdofs);

        if(embstate_n1 != Teuchos::null)
        {
          DRT::UTILS::ExtractMyValues(*embstate_n1,myval,pgdofs);
          for (std::size_t isd = 0; isd < 4; ++isd)
            veln1(isd,inode) = myval[isd];
        }

        if(embstate_n2 != Teuchos::null)
        {
          DRT::UTILS::ExtractMyValues(*embstate_n2,myval,pgdofs);
          for (std::size_t isd = 0; isd < 4; ++isd)
            veln2(isd,inode) = myval[isd];
        }

        if(embstate_n3 != Teuchos::null)
        {
          DRT::UTILS::ExtractMyValues(*embstate_n3,myval,pgdofs);
          for (std::size_t isd = 0; isd < 4; ++isd)
            veln3(isd,inode) = myval[isd];
        }

        if(embstate_n4 != Teuchos::null)
        {
          DRT::UTILS::ExtractMyValues(*embstate_n4,myval,pgdofs);
          for (std::size_t isd = 0; isd < 4; ++isd)
            veln4(isd,inode) = myval[isd];
        }

        if(embstate_n5 != Teuchos::null)
        {
          DRT::UTILS::ExtractMyValues(*embstate_n5,myval,pgdofs);
          for (std::size_t isd = 0; isd < 4; ++isd)
            veln5(isd,inode) = myval[isd];
        }

        // get the coordinates of patch element
        pxyze(0,inode) = pelenodes[inode]->X()[0];
        pxyze(1,inode) = pelenodes[inode]->X()[1];
        pxyze(2,inode) = pelenodes[inode]->X()[2];

        if (sourcedis->Name() == "fluid") // embedded fluid
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
          int offset = 0;
          if(embstate_n1 != Teuchos::null){
            for (std::size_t isd = 0; isd < 4; ++isd){
              interpolatedvec(isd+offset) += veln1(isd,inode)*shp(inode);
            }
          }

          offset += 4;
          if(embstate_n2 != Teuchos::null){
            for (std::size_t isd = 0; isd < 4; ++isd){
              interpolatedvec(isd+offset) += veln2(isd,inode)*shp(inode);
            }
          }

          offset += 4;
          if(embstate_n3 != Teuchos::null){
            for (std::size_t isd = 0; isd < 4; ++isd){
              interpolatedvec(isd+offset) += veln3(isd,inode)*shp(inode);
            }
          }

          offset += 4;
          if(embstate_n4 != Teuchos::null){
            for (std::size_t isd = 0; isd < 4; ++isd){
              interpolatedvec(isd+offset) += veln4(isd,inode)*shp(inode);
            }
          }

          offset += 4;
          if(embstate_n5 != Teuchos::null){
            for (std::size_t isd = 0; isd < 4; ++isd){
              interpolatedvec(isd+offset) += veln5(isd,inode)*shp(inode);
            }
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
//      std::cout <<" bg ele id" << pele->Id() << std::endl;
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
// we just take the old values of the embedded mesh
//
// statevbg_n (in): source state vector from gbackground fluid
// statevemb_n (in): source state vector from embedded fluid
// statevembnew_n (out)
// aledispnp (in): displacement of ale at the current time (time of
//               interpolation)
//
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SetNewEmbStatevector(const Teuchos::RCP<DRT::Discretization> bgdis,
                                                            Teuchos::RCP<Epetra_Vector>    statevbg_n,
                                                            Teuchos::RCP<Epetra_Vector>    statevemb_n,
                                                            Teuchos::RCP<Epetra_Vector>    statevembnew_n,
                                                            Teuchos::RCP<Epetra_Vector>    aledispnp,
                                                            Teuchos::RCP<Epetra_Vector>    aledispnpoldstate)
{
  std::vector<double> mydisp(4);
  std::vector<int> pgdofs(4);

  if (gmsh_debug_out_)
    GmshOutputForInterpolateFSI(bgdis,aledispnp,aledispnpoldstate);

  // dummy state vectors
  Teuchos::RCP<Epetra_Vector>    vec1;
  Teuchos::RCP<Epetra_Vector>    vec2;
  Teuchos::RCP<Epetra_Vector>    vec3;
  Teuchos::RCP<Epetra_Vector>    vec4;


  for (int lnid=0; lnid<embdis_->NumMyRowNodes(); lnid++)
  {
    DRT::Node* embnode = embdis_->lRowNode(lnid);

    DRT::UTILS::ConditionSelector conds(*embdis_, "FSICoupling");
    if (conds.ContainsNode(embnode->Id())==false)//if no fsi node
    {
      //std::cout << "unknown node id " << embnode->Id() << std::endl;
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
        LINALG::Matrix<20,1> interpolatedvec(true);
        //std::cout << "pele id " << pele->Id() << std::endl;

        insideelement = ComputeSpacialToElementCoordAndProject(
            pele,embnodecords,interpolatedvec,
            statevemb_n,vec1, vec2, vec3, vec4,
            aledispnpoldstate,embdis_);
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
        IO::cout << "no embedded element found for: " <<  embnode->Id() << IO::endl;
         count = 0;
         // check all background elements to find the right one
         for (int e=0; e<bgdis->NumMyColElements(); e++)
         {
            DRT::Element* bgele = bgdis->lColElement(e);
            LINALG::Matrix<4,1> interpolatedvec(true);

            insideelement = ComputeSpacialToElementCoordAndProject2(bgele,embnodecords,interpolatedvec,statevbg_n,aledispnpoldstate,bgdis);
            // this is commented out because this method is just called for time integration of interpolated-Ale approach which is just used for
            // 1D-examples and for every other complicated fsi simulation the Ale-partitioned is used
/*            if (insideelement)
            {
              // here set state
              (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[0])] = interpolatedvec(0);
              (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[1])] = interpolatedvec(1);
              (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[2])] = interpolatedvec(2);
             (*statevembnew_n)[statevembnew_n->Map().LID(embdis_->Dof(embnode)[3])] = interpolatedvec(3);
              break;
           }*/
           count ++;
       }
/*          if (count == bgdis->NumMyColElements())
            std::cout << "no bg elements found for: " << embnode->Id() << std::endl;*/
        }
    }
  }

}//SetNewEmbStatevector

// ------------------------------------------------------------------------
// find an appropriate radius for the search tree. The minimum radius is the
// max diameter of the surfaces of the first element of embedded discretization
// ------------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::SearchRadius()
{
    const DRT::Element* actele = embdis_->lColElement(0);
    const DRT::Node* const* nodes = actele->Nodes();

  LINALG::Matrix<3,8> xyze(true);

    switch ( actele->Shape() )
    {
      case DRT::Element::hex8:
      {
        for(int i=0; i<8; i++)
        {
          const double* x = nodes[i]->X();
          xyze(0,i)=x[0];
          xyze(1,i)=x[1];
          xyze(2,i)=x[2];
        }
        break;
      }
      default:
        dserror("Element-type not supported here!");
    }

    double line0 = 0.0;
    double line1 = 0.0;
    double line2 = 0.0;
    double line3 = 0.0;
    double line4 = 0.0;
    double line5 = 0.0;

    for(int i=0; i<3; i++)
    {
      line0 += (xyze(i,0) - xyze(i,2))*(xyze(i,0) - xyze(i,2));
      line1 += (xyze(i,1) - xyze(i,4))*(xyze(i,1) - xyze(i,4));
      line2 += (xyze(i,1) - xyze(i,6))*(xyze(i,1) - xyze(i,6));
      line3 += (xyze(i,3) - xyze(i,6))*(xyze(i,3) - xyze(i,6));
      line4 += (xyze(i,0) - xyze(i,7))*(xyze(i,0) - xyze(i,7));
      line5 += (xyze(i,4) - xyze(i,6))*(xyze(i,4) - xyze(i,6));
    }

    line0 = sqrt(line0);
    line1 = sqrt(line1);
    line2 = sqrt(line2);
    line3 = sqrt(line3);
    line4 = sqrt(line4);
    line5 = sqrt(line5);

    minradius_ =  searchradius_fac_*std::max(line0,std::max(line1,std::max(line2,std::max(line3,std::max(line4,line5)))));

}//SearchRadius

// -------------------------------------------------------------------
// Gmsh Output
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::GmshOutput(const Teuchos::RCP<DRT::Discretization>    bgdis)
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
      std::map<int, std::vector<int> >::const_iterator iter = stdnoden_.find(actnode->Id());
      std::map<int, std::vector<int> >::const_iterator iteren = enrichednoden_.find(actnode->Id());
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
      std::map<int, std::vector<int> >::const_iterator iter = stdnodenp_.find(actnode->Id());
      std::map<int, std::vector<int> >::const_iterator iteren = enrichednodenp_.find(actnode->Id());
      if (iter != stdnodenp_.end()) kind = 1;//std
      if (iteren != enrichednodenp_.end()) kind = 2; // enriched
      IO::GMSH::cellWithScalarToStream(DRT::Element::point1, kind, pos, gmshfilecontent);
    }
    gmshfilecontent << "};\n";
  }
}
// -------------------------------------------------------------------
// Gmsh Output for interpolated-Ale FSI-Approach
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::GmshOutputForInterpolateFSI(const Teuchos::RCP<DRT::Discretization>    bgdis,
                                                               Teuchos::RCP<Epetra_Vector>    aledispnp,
                                                                   Teuchos::RCP<Epetra_Vector>    aledispnpoldstate)
{
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

        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        actele->LocationVector(*embdis_, lm, lmowner, lmstride);

        std::vector<double> myolddisp(lm.size());
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

          std::vector<int> lm;
          std::vector<int> lmowner;
          std::vector<int> lmstride;
          actele->LocationVector(*embdis_, lm, lmowner, lmstride);

          std::vector<double> mydisp(lm.size());
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

          std::vector<int> lm;
          std::vector<int> lmowner;
          std::vector<int> lmstride;
          actele->LocationVector(*embdis_, lm, lmowner, lmstride);

          std::vector<double> mydisp(lm.size());
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
}
// -----------------------------------------------------------------------
// all cut elements at tn with full dofs are included in incompressibility
// patch. If there  are nodes in projectednodeids which still are not in
// incompressibility patch we add them afterwards. These added Elements should
// have the full dofs again.
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::PatchelementForIncompressibility(const Teuchos::RCP<DRT::Discretization>     bgdis,
                                                                        Teuchos::RCP<XFEM::FluidWizardMesh>             wizard_n,
                                                                        Teuchos::RCP<XFEM::FluidWizardMesh>             wizard_np,
                                                                        Teuchos::RCP<LINALG::MapExtractor> dbcmaps )
{
  //---------------------------------------------
  // find the patch at time t_n

  // delete the elements and nodes of the last time step
  incompnodeids_set_.clear();
  incompelementids_set_.clear();


  // call loop over elements
  const int numele = bgdis->NumMyRowElements();
  for (int i=0; i<numele; ++i)
  {
    DRT::Element* actele = bgdis->lRowElement(i);
    int numnodes = actele->NumNode();

    // get the element handle of actele at time tn
    GEO::CUT::ElementHandle * e = wizard_n->GetElement( actele );

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
        incompelementids_set_.insert(actele->Id());

        //get the nodes adjust to this element
        const int* nodeids = actele->NodeIds();

        for(int inode=0; inode<actele->NumNode(); ++inode)
        {
          const int nodegid = nodeids[inode];

          //insert the node
          incompnodeids_set_.insert(nodegid);
        }
      }
    }
  }

  //------------------------------------------------
  //  check if all projected nodes are included
  for(std::set<int>::iterator iter = projectednodeids_.begin(); iter!= projectednodeids_.end();
      iter++)
  {
    std::set<int>::const_iterator iterpatchnodes = incompnodeids_set_.find(*iter);

    if (iterpatchnodes == incompnodeids_set_.end())
    {
      IO::cout << "STOP!! Nodes found which were in projected set but not in incompressibility patch! " << *iter << IO::endl;

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
          incompelementids_set_.insert(ele_adj->Id());
          IO::cout << "element found " << ele_adj->Id() << IO::endl;

          //get the nodes of this element
          const int* nodeids = ele_adj->NodeIds();

          for (int inode=0; inode<ele_adj->NumNode(); ++inode)
          {
            const int nodegid = nodeids[inode];

            //insert the node
            incompnodeids_set_.insert(nodegid);
          }
          // break the element loop, we've already found one element for this node
          break;
        }
      }
    }
  }

  //------------------------------------------
  //check it again..
  for(std::set<int>::iterator iter = projectednodeids_.begin(); iter!= projectednodeids_.end();
      iter++)
  {
    std::set<int>::const_iterator iterpatchnodes = incompnodeids_set_.find(*iter);

    if (iterpatchnodes == incompnodeids_set_.end())
      dserror("BUG!! Nodes found which were in projected set but not in incompressibility patch!",*iter);
  }


//  debug output
//     for(std::set<int>::iterator iter = incompelementids_set_.begin(); iter!= incompelementids_set_.end();
//         iter++)
//     {
//       std::set<int>::const_iterator iterpatchele = incompelementids_set_.find(*iter);
//       std::cout << "all eles" << *iter << std::endl;
//       DRT::Node ** elenodes = bgdis->gElement(*iter)->Nodes();
//       for(int inode=0; inode<bgdis->gElement(*iter)->NumNode(); ++inode)
//       {
//         std::cout <<  "bgnode id " <<elenodes[inode]->Id() ;
//         std::cout << " dofs " << bgdis->NumDof(elenodes[inode]) << std::endl;
//       }
//     }


  //---------------------------------
  // Gmsh debug output
  if (gmsh_debug_out_)
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("incom_patch", step_, 5, 0, bgdis->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    {
      // draw bg elements with associated gid
      gmshfilecontent << "View \" " << "bg Element->Id() \" {\n";
      for (int i=0; i<bgdis->NumMyColElements(); ++i)
      {
        DRT::Element* actele = bgdis->lColElement(i);
//         GEO::CUT::ElementHandle * e = wizard_n->GetElement( actele );
        std::set<int>::const_iterator iter = incompelementids_set_.find(actele->Id());
        if ( iter != incompelementids_set_.end())
          IO::GMSH::elementAtInitialPositionToStream(1.0, actele, gmshfilecontent);
        else
          IO::GMSH::elementAtInitialPositionToStream(0.0, actele, gmshfilecontent);
      };
      gmshfilecontent << "};\n";
    }
  }
}

// -------------------------------------------------------------------
// build an incompressibility discretization
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::PrepareIncompDiscret(const Teuchos::RCP<DRT::Discretization>     bgdis,
                                                            Teuchos::RCP<XFEM::FluidWizardMesh>             wizard_np)
{

  // generate an empty boundary discretisation
  incompdis_ = Teuchos::rcp(new DRT::Discretization((std::string)"incompressibility discretisation",
                                           Teuchos::rcp(bgdis->Comm().Clone())));

  std::set<int> incompelementids_set_all;
  std::set<int> incompnodeids_set_all;

   // Gather all informations from all processors
  std::vector<int> incompdofsAllproc;
  std::vector<int> incompveldofsAllproc;

  // information how many processors work at all
  std::vector<int> allproc(bgdis->Comm().NumProc());

  // in case of n processors allproc becomes a vector with entries (0,1,...,n-1)
  for (int i=0; i<bgdis->Comm().NumProc(); ++i) allproc[i] = i;

  LINALG::Gather<int>(incompelementids_set_,incompelementids_set_all,(int)bgdis->Comm().NumProc(),&allproc[0],bgdis->Comm());
  LINALG::Gather<int>(incompnodeids_set_,incompnodeids_set_all,(int)bgdis->Comm().NumProc(),&allproc[0],bgdis->Comm());


  // determine sets of col und row nodes
  std::set<int> adjacent_row;
  std::set<int> adjacent_col;


  // loop all column elements and label all row nodes next to a MHD node
  for (int i=0; i<bgdis->NumMyColElements(); ++i)
  {
    DRT::Element* actele = bgdis->lColElement(i);

    // get the node ids of this elements
    const int  numnode = actele->NumNode();
    const int* nodeids = actele->NodeIds();

    bool found=false;

    std::set<int>::const_iterator iter = incompelementids_set_all.find(actele->Id());
    if ( iter != incompelementids_set_all.end()) found=true;

    if(found==true)
    {
      // loop nodeids
      for(int rr=0;rr<numnode;++rr)
      {
        int gid=nodeids[rr];

        if ((bgdis->NodeRowMap())->LID(gid)>-1)
        {
          adjacent_row.insert(gid);
        }
        adjacent_col.insert(gid);
      }
    }
  }

  // add nodes to incompressibility discretisation
  for(std::set<int>::iterator id = adjacent_row.begin();
      id!=adjacent_row.end(); ++id)
  {
    DRT::Node* actnode=bgdis->gNode(*id);

    Teuchos::RCP<DRT::Node> incompnode =Teuchos::rcp(actnode->Clone());

    incompdis_->AddNode(incompnode);
  }


  // loop all row elements and add all elements with a MHD node
  for (int i=0; i<bgdis->NumMyRowElements(); ++i)
  {
    DRT::Element* actele = bgdis->lRowElement(i);

    bool found=false;

    // check if incompressibility element
    std::set<int>::const_iterator iter = incompelementids_set_all.find(actele->Id());
    if ( iter != incompelementids_set_all.end()) found=true;

    // yes, we have a MHD condition
    if(found==true)
    {
      Teuchos::RCP<DRT::Element> incompele =Teuchos::rcp(actele->Clone());

      incompdis_->AddElement(incompele);
    }
  }

  //incompelementids_set_ needs a full NodeRowMap and a NodeColMap
  Teuchos::RCP<Epetra_Map> newrownodemap;
  Teuchos::RCP<Epetra_Map> newcolnodemap;

  std::vector<int> rownodes;

  // convert std::set to std::vector
  for(std::set<int>::iterator id = adjacent_row.begin();
      id!=adjacent_row.end();
      ++id)
  {
    rownodes.push_back(*id);
  }

  // build noderowmap for new distribution of nodes
  newrownodemap = Teuchos::rcp(new Epetra_Map(-1,
                                     rownodes.size(),
                                     &rownodes[0],
                                     0,
                                     incompdis_->Comm()));


  std::vector<int> colnodes;
  for(std::set<int>::iterator id = adjacent_col.begin();
      id!=adjacent_col.end();
      ++id)
  {
    colnodes.push_back(*id);
  }

  // build nodecolmap for new distribution of nodes
  newcolnodemap = Teuchos::rcp(new Epetra_Map(-1,
                                     colnodes.size(),
                                     &colnodes[0],
                                     0,
                                     incompdis_->Comm()));

  incompdis_->Redistribute(*newrownodemap,*newcolnodemap,false,false,false);
  Teuchos::RCP<DRT::DofSet> newdofset=Teuchos::rcp(new DRT::TransparentIndependentDofSet(bgdis,true,wizard_np));
  incompdis_->ReplaceDofSet(newdofset); // do not call this with true!!
  incompdis_->FillComplete();
}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void XFEM::XFluidFluidTimeIntegration::EvaluateIncompressibility(const Teuchos::RCP<DRT::Discretization>  bgdis,
                                                                 Teuchos::RCP<XFEM::FluidWizardMesh>          wizard)
{

  C_ = LINALG::CreateVector(*incompdis_->DofRowMap(),true);

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

  //---------------------------------------
  // Find the vector c:
  //
  //   __   _                  __   _                 __   _
  //   \   |  dN_A             \   |  dN_A            \   |  dN_A
  // ( /   |  ---- dOmega_e1,  /   |  ---- dOmega_e1, /   |  ---- dOmega_e1, ...)
  //   -- -   dx               -- -   dy              -- -   dz
  //   A                       A                      A

  DRT::Element::LocationArray la( 1 );

  // loop over column elements of bgdis
  const int numcolele = incompdis_->NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* actele = incompdis_->lColElement(i);

    Teuchos::RCP<MAT::Material> mat = actele->Material();
    DRT::ELEMENTS::Fluid * ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actele );

    GEO::CUT::ElementHandle * e = wizard->GetElement( actele );

    Epetra_SerialDenseVector C_elevec;
    // xfem element
    if ( e!=NULL )
    {

      std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
      std::vector< std::vector<int> > nds_sets;
      std::vector<std::vector< DRT::UTILS::GaussIntegration > >intpoints_sets;
      INPAR::CUT::VCellGaussPts VolumeCellGaussPointBy = DRT::INPUT::IntegralValue<INPAR::CUT::VCellGaussPts>(params_.sublist("XFEM"), "VOLUME_GAUSS_POINTS_BY");

      e->GetCellSets_DofSets_GaussPoints( cell_sets, nds_sets, intpoints_sets, VolumeCellGaussPointBy, false); //(include_inner=false)

      if(cell_sets.size() != intpoints_sets.size()) dserror("number of cell_sets and intpoints_sets not equal!");
      if(cell_sets.size() != nds_sets.size()) dserror("number of cell_sets and nds_sets not equal!");

      int set_counter = 0;

      if (cell_sets.size() == 0)
      {
        IO::cout << "Warning: Element " << actele->Id() << " has all it's volume-cells in void. Check if all nodes are "<<
          "included in other elements." << IO::endl;
        continue;
      }

      for( std::vector< GEO::CUT::plain_volumecell_set>::iterator s=cell_sets.begin();
           s!=cell_sets.end();
           s++)
      {
        //GEO::CUT::plain_volumecell_set & cells = *s;
        const std::vector<int> & nds = nds_sets[set_counter];

        // get element location vector, dirichlet flags and ownerships
        actele->LocationVector(*incompdis_,nds,la,false);

        // number of dofs for background element
        // ndof contains all dofs of velocity and pressure. (we need just the velocity)
        const size_t ndof  = la[0].lm_.size();
        C_elevec.Reshape(ndof,1);

        for( unsigned cellcount=0;cellcount!=cell_sets[set_counter].size();cellcount++ )
        {
          // call element method
          DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(actele->Shape(), "xfem")->CalculateContinuityXFEM(ele,
                                                                                                         *incompdis_,
                                                                                                         la[0].lm_,
                                                                                                         C_elevec,
                                                                                                         intpoints_sets[set_counter][cellcount]);


        }
        set_counter += 1;
      }

    }
    else
    {
      // get element location vector, dirichlet flags and ownerships
      actele->LocationVector(*incompdis_,la,false);

      const size_t ndof = la[0].lm_.size();

      C_elevec.Reshape(ndof,1);

      DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(actele->Shape(), "xfem")->CalculateContinuityXFEM(ele,
                                                                                                     *incompdis_,
                                                                                                     la[0].lm_,
                                                                                                     C_elevec);

    }
    LINALG::Assemble(*C_, C_elevec, la[0].lm_, la[0].lmowner_);
  }
}

// -------------------------------------------------------------------
// Solve the incompressibility Optimitzation problem
// -------------------------------------------------------------------

void  XFEM::XFluidFluidTimeIntegration::SolveIncompOptProb(Teuchos::RCP<Epetra_Vector>   initialvel)
{
  // ----------------------
  // Prepare the C vector:
  // The vector C still includes the pressure degrees of freedom which are
  // zero. So we need to eliminate them from C so that the vector C has just
  // the velocity degrees of freedom -> C_vel
  LINALG::MapExtractor      velpressplitter;
  int numdim = 3;
  FLD::UTILS::SetupFluidSplit(*incompdis_, numdim, 1,velpressplitter);
  Teuchos::RCP<Epetra_Vector> C_vel = velpressplitter.ExtractOtherVector(C_);

  // ----------------------
  // Create needed C_vel-dofmaps
  std::vector<int> C_vel_dofids;

  for(int i=0; i<C_vel->MyLength(); ++i)
  {
    C_vel_dofids.push_back(C_vel->Map().GID(i));
  }

  // build dofrowmap for velocity dofs
  Teuchos::RCP<Epetra_Map> veldofrowmap = Teuchos::rcp(new Epetra_Map(-1,
                                                            C_vel_dofids.size(),
                                                            &C_vel_dofids[0],
                                                            0,
                                                            incompdis_->Comm()));

  Teuchos::RCP<Epetra_Vector> C_vel_copy =  LINALG::CreateVector(*veldofrowmap,true);
  C_vel_copy->Update(1.0, *C_vel, 0.0);

  // -----------------------
  // Initialize Q as identity matrix
  int maxnumberofentries = C_vel->MyLength()*C_vel->MyLength();

  Teuchos::RCP<Epetra_CrsMatrix> Q = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*veldofrowmap,maxnumberofentries,false));
  double ones = 1.0;
  for(int i=0; i<C_vel->MyLength(); ++i)
  {
    int myGID = C_vel->Map().GID(i);
    int err = Q->InsertGlobalValues(myGID,1,&ones,&myGID);
    if (err<0) dserror("Epetra_CrsMatrix::InsertGlobalValues returned err=%d",err);
  }
  Q->FillComplete();

  Teuchos::RCP<LINALG::SparseMatrix> Q_spr = Teuchos::rcp(new LINALG::SparseMatrix(Q));
  Q_spr->Complete();

  // ------------------------
  // Transform the constraint cT.u* = 0 to (Qc)T.Qu* = 0
  //
  // With appropiate Givens Rotations we annul all of the entires of c
  // but the last entry. First we start with the first entry (next).
  // The entry "next" will be annulated through the second
  // entry (next_us). We do the same for all entries of c until all entries
  // of c are zero beside the last entry.
  //
  // next := entry to eliminate, next_us := entry to make next to zero
  //
  // The vector we find at last is Qc


  int pair = veldofrowmap->NumMyElements()-1;
  int maxpair;
  incompdis_->Comm().MaxAll(&pair, &maxpair, 1);

  int next_us = 0;
  double TOL = 1.0e-14;
  int maxgid = C_vel->Map().MaxMyGID();

  // the loop over all pairs of C_vel
  for (int next=0; next<maxpair; ++next)
  {
    // build the rotation matrix for current next and next_us
    Teuchos::RCP<Epetra_CrsMatrix> Q_i = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*veldofrowmap,2));
    // first build Q_i as an identity matrix
    double myval = 1.0;
    for(int j=0; j<C_vel->MyLength(); ++j)
    {
      int myGID = C_vel->Map().GID(j);
      Q_i->InsertGlobalValues(myGID,1,&myval,&myGID);
    }

    if (next < ((C_vel->MyLength())-1))
    {
      // if the next values is not zero find next+1 and do the givens,
      // otherwiese Q_i remains an idendity matrix
      if ( abs((*C_vel)[next]) > TOL )
      {
        next_us = next+1;
        while ( abs((*C_vel)[next_us]) < TOL )
        {
          next_us++;
        }

        if ( next_us <= maxgid )
        {
          // std::cout << "next yy "<< next << " " << (*C_vel)[next] << std::endl;
          // std::cout << "next_us yy "<< next_us << " " <<  (*C_vel)[next_us] << std::endl;

          double Nenner = sqrt((*C_vel)[next]*(*C_vel)[next]+
                               (*C_vel)[next_us]*(*C_vel)[next_us]);
          Nenner = 1/Nenner;

          double cphi = Nenner*(*C_vel)[next_us];
          double sphi = Nenner*(*C_vel)[next];

          int nextGID = veldofrowmap->GID(next);
          int next_usGID = veldofrowmap->GID(next_us);

          double sphi_min = -sphi;
          Q_i->ReplaceGlobalValues(nextGID,1,&cphi,&nextGID);
          Q_i->InsertGlobalValues(nextGID,1,&sphi_min,&next_usGID);
          Q_i->InsertGlobalValues(next_usGID,1,&sphi,&nextGID);
          Q_i->ReplaceGlobalValues(next_usGID,1,&cphi,&next_usGID);
        }
      }
    }
    // do nothing or wait until all processors are finished
    else { }

    incompdis_->Comm().Barrier();
    Q_i->FillComplete(*veldofrowmap,*veldofrowmap);

    Teuchos::RCP<LINALG::SparseMatrix> Q_i_spr = Teuchos::rcp(new LINALG::SparseMatrix(Q_i));
    Q_i_spr->Complete(*veldofrowmap,*veldofrowmap);

    // Update of C_vel (which is after every update Qc)
    /* _     _
       | s  -c || x_i |    | s*x_i-c*x_j |   |      0      |
       |       ||     |  = |             | = |             |
       | c   s || x_j |    | c*x_i+s*x_j |   | c*x_i+s*x_j |
       -     -
       x_i: next
       x_j: next_us
    */

    Teuchos::RCP<Epetra_Vector> C_vel_test =  LINALG::CreateVector(*veldofrowmap,true);
    (Q_i_spr->EpetraMatrix())->Multiply(false,*C_vel,*C_vel_test);
    C_vel->Update(1.0,*C_vel_test,0.0);

    // build the final Q  = Qn..Q3*Q2*Q1
    Teuchos::RCP<LINALG::SparseMatrix> Q_final = Teuchos::rcp(new LINALG::SparseMatrix(*veldofrowmap,C_vel->MyLength(),false,true));
    Q_final = Multiply(*Q_i_spr,false,*Q_spr,false,false,false);

    // save this old Q_final
    Q_spr = Q_final;
  }

  incompdis_->Comm().Barrier();

  //-----------------------------------------------------
  // Commuinate the last entry of C_vel of each processor

  // build an allreduced vector of all last entries of C_vel
  std::vector<int> C_vellast;
  std::vector<int> C_vellast_All;

  if ((maxgid>-1) and ((*C_vel)[veldofrowmap->LID(maxgid)]>TOL))
    C_vellast.push_back(maxgid);

  //information how many processors work at all
  std::vector<int> allproc(incompdis_->Comm().NumProc());

  LINALG::Gather<int>(C_vellast,C_vellast_All,(int)incompdis_->Comm().NumProc(),
                      &allproc[0],incompdis_->Comm());

  // build dofrowmap of last entries
  Teuchos::RCP<Epetra_Map> lastdofrowmap = Teuchos::rcp(new Epetra_Map(-1,
                                                             C_vellast_All.size(),
                                                             &C_vellast_All[0],
                                                             0,
                                                             incompdis_->Comm()));

  const Epetra_Map lastallreduced = *LINALG::AllreduceOverlappingEMap(*lastdofrowmap);
  Teuchos::RCP<Epetra_Vector> C_vel_last =  LINALG::CreateVector(lastallreduced,true);
  int mylastGID = lastallreduced.MaxMyGID();
  LINALG::Export(*C_vel,*C_vel_last);

  // std::cout << "C_vellast " << *C_vel_last << std::endl;
  // std::cout << "length " << C_vel_last->MyLength() << std::endl;

  // loop over all processors
  if (C_vel_last->MyLength() > 1)
  {
    for(int pr=0; pr<incompdis_->Comm().NumProc()-1; ++pr)
    {
      // build the rotation matrix for current next and next_us
      Teuchos::RCP<Epetra_CrsMatrix> Q_i = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*veldofrowmap,2));
      // first build Q_i as an identity matrix
      double myval = 1.0;
      for(int j=0; j<C_vel->MyLength(); ++j)
      {
        int myGID = C_vel->Map().GID(j);
        Q_i->InsertGlobalValues(myGID,1,&myval,&myGID);
      }

      if(C_vel_last->Map().MyLID(pr+1))
      {
        double next_val = (*C_vel_last)[pr];
        double nextus_val = (*C_vel_last)[pr+1];

        // std::cout << "val " << next_val << " " << nextus_val << std::endl;

        double Nenner = sqrt(next_val*next_val+nextus_val*nextus_val);
        Nenner = 1/Nenner;

        double cphi = Nenner*nextus_val;
        double sphi = Nenner*next_val;

        int nextGID = lastallreduced.GID(pr);
        int next_usGID = lastallreduced.GID(pr+1);
        mylastGID = lastallreduced.GID(pr+1);

        double sphi_min = -sphi;
        Q_i->ReplaceGlobalValues(nextGID,1,&cphi,&nextGID);
        Q_i->InsertGlobalValues(nextGID,1,&sphi_min,&next_usGID);
        Q_i->InsertGlobalValues(next_usGID,1,&sphi,&nextGID);
        Q_i->ReplaceGlobalValues(next_usGID,1,&cphi,&next_usGID);
        Q_i->FillComplete(*veldofrowmap,*veldofrowmap);


        Teuchos::RCP<LINALG::SparseMatrix> Q_i_spr = Teuchos::rcp(new LINALG::SparseMatrix(Q_i));
        Q_i_spr->Complete(*veldofrowmap,*veldofrowmap);

        Teuchos::RCP<Epetra_Vector> C_vel_test =  LINALG::CreateVector(*veldofrowmap,true);
        (Q_i_spr->EpetraMatrix())->Multiply(false,*C_vel,*C_vel_test);
        C_vel->Update(1.0,*C_vel_test,0.0);

        // build the final Q  = Qn..Q3*Q2*Q1
        Teuchos::RCP<LINALG::SparseMatrix> Q_final = Teuchos::rcp(new LINALG::SparseMatrix(*veldofrowmap,C_vel->MyLength(),false,true));
        Q_final = Multiply(*Q_i_spr,false,*Q_spr,false,false,false);

        // save this old Q_final
        Q_spr = Q_final;
      }
    }
  }

  //------------------------------------------------------------
  // Update the initial velocity vector:
  //
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

  // the original velocity vector which we want to improve
  Teuchos::RCP<Epetra_Vector> vel_org = LINALG::CreateVector(*veldofrowmap,true);
  LINALG::Export(*initialvel,*vel_org);

  // Qu
  Teuchos::RCP<Epetra_Vector> Qunext = LINALG::CreateVector(*veldofrowmap,true);
  (Q_spr->EpetraMatrix())->Multiply(false,*vel_org,*Qunext);

  // we need the last entry of Qu (Qu_next)
  for (int i=0; i<Qunext->MyLength(); ++i)
  {
     int mygid = Qunext->Map().GID(i);
     if ( mygid != mylastGID)
       (*Qunext)[Qunext->Map().LID(mygid)] = 0.0;
  }

  Teuchos::RCP<Epetra_Vector> QTQu = LINALG::CreateVector(*veldofrowmap,true);
  (Q_spr->EpetraMatrix())->Multiply(true,*Qunext,*QTQu);

  Teuchos::RCP<Epetra_Vector> u_incomp = LINALG::CreateVector(*veldofrowmap,true);
  u_incomp->Update(1.0, *vel_org, 0.0);
  u_incomp->Update(-1.0, *QTQu, 1.0);

  //incompressibility check before solving the optimization problem
  double sum = 0.0;
  vel_org->Dot(*C_vel_copy, &sum);

  double sum_opt = 0.0;
  u_incomp->Dot(*C_vel_copy, &sum_opt);

  if (myrank_ == 0)
  {
    IO::cout << " Incompressibility Check.. " << IO::endl;
    IO::cout << " Original:  "  << sum << ",  After solving optimization problem: "  << sum_opt << IO::endl;
  }

  LINALG::Export(*(u_incomp),*(initialvel));

}


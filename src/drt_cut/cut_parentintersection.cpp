/*-----------------------------------------------------------------------------------------------*/
/**
\file cut_parentintersection.cpp

\brief provides the basic functionality for cutting a mesh

<pre>
\maintainer Benedikt Schott and Magnus Winter
            schott@lnm.mw.tum.de, winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\level 2
*/
/*------------------------------------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_discret.H"

#include "cut_volumecell.H"

#include "cut_parentintersection.H"

/*------------------------------------------------------------------------------------------------*
 * Create nodal dofset sets within the parallel cut framework
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::ParentIntersection::CreateNodalDofSet( bool include_inner, const DRT::Discretization& dis)
{
  dis.Comm().Barrier();

  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- CreateNodalDofSet" );


  std::set<int> eids; // eids of elements that are involved in CUT and include ele_vc_set_inside/outside (no duplicates!)

  Mesh & m = NormalMesh();

  const INPAR::CUT::NodalDofSetStrategy strategy = options_.GetNodalDofSetStrategy();

  // nodes used for CUT std::map<node->ID, Node>, shadow nodes have ID<0
  std::map<int, Node* > nodes;
  m.GetNodeMap( nodes );

  //===============
  // STEP 1: create for each node involved in CUT nodal cell sets for all nodes of adjacent elements of this element
  //         nodal cell sets are sets of volumecells, that are connected via subelements in a real element
  //         ° for linear elements a set contains just one single volumecell in its plain_volumecell_set(with usually more than one integrationcells)
  //         ° for quadratic elements a set contains connected volumecell sets (sorted for inside and outside connections)
  //
  //         FindDOFSets() finds also the connections of volumecell sets between adjacent elements
  //         finally, each found DOFSet around a 1-ring of this node maintains its own set of DOFs
  //===============
  // for each node (which is not a shadow node) we have to built nodalDofSets
  // each connected set of volumecells (connected via adjacent elements - not only via subelements -) gets an own dofset
  for ( std::map<int, Node* >::iterator i=nodes.begin();
      i!=nodes.end();
      ++i )
  {
    Node * n = i->second;
    int n_gid = n->Id();

    std::vector<int> sourrounding_elements;

    // get all adjacent elements to this node if this is a real (- not a shadow -) node
    if( n_gid >= 0 )
    {
      DRT::Node * node = dis.gNode(n_gid);

      // get adjacent elements for this node
      const DRT::Element*const* adjelements = node->Elements();

      for (int iele=0; iele < node->NumElement(); iele++)
      {
        int adj_eid =  adjelements[iele]->Id();

        // get its elementhandle
        GEO::CUT::ElementHandle * e = GetElement( adj_eid );

        if( e!=NULL )
        {
          sourrounding_elements.push_back(adj_eid);
        }
      } // end loop over adjacent elements



      // each node stores all its sets of volumecells,
      // includes all volumecell_sets that are connected (within a whole adjacent element) via subelements, inside and outside sets
      // (appended for all adjacent elements of a node)
      std::map<Node*, std::vector<plain_volumecell_set> > nodal_cell_sets_inside;    // for each node a vector of volumecell_sets
      std::map<Node*, std::vector<plain_volumecell_set> > nodal_cell_sets_outside;   // for each node a vector of volumecell_sets


      // includes all volumecell_sets that are connected (within a whole adjacent element) via subelements, inside and outside sets
      std::vector<plain_volumecell_set> cell_sets;  // vector of all volumecell sets connected within elements adjacent to this node

      //split for inside and outside
      std::vector<plain_volumecell_set> cell_sets_inside; // sets of volumecells connected between subelements
      std::vector<plain_volumecell_set> cell_sets_outside;

      FindNodalCellSets(include_inner,
          eids,
          sourrounding_elements,
          nodal_cell_sets_inside,
          nodal_cell_sets_outside,
          cell_sets_inside,
          cell_sets_outside,
          cell_sets );

#if(0)
      // Gmsh output for the NodalCellSet of node with id n_gid
      if(n_gid == 6)
      {
        DumpGmshNodalCellSet(nodal_cell_sets_outside, dis);
        DumpGmshCellSets(cell_sets_outside, dis);
      }
#endif


      // finds also the connections of volumecell sets between adjacent elements
      // finally, each found DOFSet around a 1-ring of the node maintains its own set of DOFs
      if( include_inner )
      {
        //WARNING:
        //This is necessary to have to set the "standard values" at the first DOF-set, and the
        //ghost values on the following DOFS.
        // It is important for result check AND later for time-integration!
        if(n->Position()==Point::outside)
        {
          n->FindDOFSetsNEW( nodal_cell_sets_outside, cell_sets_outside);
          n->FindDOFSetsNEW( nodal_cell_sets_inside, cell_sets_inside);
        }
        else
        {
          n->FindDOFSetsNEW( nodal_cell_sets_inside, cell_sets_inside);
          n->FindDOFSetsNEW( nodal_cell_sets_outside, cell_sets_outside);
        }
      }
      else
      {
        n->FindDOFSetsNEW( nodal_cell_sets_outside, cell_sets_outside);
      }

      // sort the dofsets for this node after FindDOFSetsNEW
      n->SortNodalDofSets();


      if(strategy == INPAR::CUT::NDS_Strategy_OneDofset_PerNodeAndPosition)
      {

        /* combine the (ghost) dofsets for this node w.r.t each phase to avoid
         * multiple ghost nodal dofsets for a certain phase */
        n->CollectNodalDofSets();

      } // otherwise do nothing

    } // end if n_gid >= 0
  } // end loop over nodes


  //===============
  // STEP 2: for each element that contains volumecell_sets (connections via subelements...),
  // all nodes of this element have to know the dofset_number for each set of volumecells
  //===============
  for( std::set<int>::iterator i= eids.begin(); i!= eids.end(); i++)
  {
    TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- STEP 2" );


    int eid = *i;

    // get the nodes of this element
    // get the element via discret
    DRT::Element * e = dis.gElement(eid);

    if( e == NULL) dserror(" element not found, this should not be! ");

    // get the nodes of this element
    int numnode = e->NumNode();
    const int* nids = e->NodeIds();
    std::vector<Node*> nodes(numnode);

    for(int i=0; i< numnode; i++)
    {
      Node * node = GetNode(nids[i]);

      if(node == NULL) dserror("node not found!");

      nodes[i] = node;
    }

    ElementHandle* eh = GetElement(eid);

    // get inside and outside cell_sets connected within current element

    const std::vector<plain_volumecell_set> & ele_vc_sets_inside = eh->GetVcSetsInside();
    const std::vector<plain_volumecell_set> & ele_vc_sets_outside = eh->GetVcSetsOutside();

    std::vector<std::vector<int> > & nodaldofset_vc_sets_inside  = eh->GetNodalDofSet_VcSets_Inside();
    std::vector<std::vector<int> > & nodaldofset_vc_sets_outside = eh->GetNodalDofSet_VcSets_Outside();

    std::vector<std::map<int,int> > & vcsets_nid_dofsetnumber_map_toComm_inside  = eh->Get_NodeDofsetMap_VcSets_Inside_forCommunication();
    std::vector<std::map<int,int> > & vcsets_nid_dofsetnumber_map_toComm_outside = eh->Get_NodeDofsetMap_VcSets_Outside_forCommunication();

    if(include_inner )
    {
      ConnectNodalDOFSets(
          nodes,
          include_inner,
          dis,
          ele_vc_sets_inside,
          nodaldofset_vc_sets_inside,
          vcsets_nid_dofsetnumber_map_toComm_inside);
    }

    ConnectNodalDOFSets(
        nodes,
        include_inner,
        dis,
        ele_vc_sets_outside,
        nodaldofset_vc_sets_outside,
        vcsets_nid_dofsetnumber_map_toComm_outside);



  }


#if(0)
  DumpGmshNumDOFSets("numdofsets_debug", include_inner, dis);
#endif
}


/*--------------------------------------------------------------------------------------*
 | fill parallel DofSetData with information that has to be communicated   schott 03/12 |
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::ParentIntersection::FillParallelDofSetData(
    std::vector<Teuchos::RCP<DofSetData> >& parallel_dofSetData,
    const DRT::Discretization& dis,
    bool include_inner
)
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- FillParallelDofSetData" );

  // find volumecell sets and non-row nodes for that dofset numbers has to be communicated parallel
  // the communication is done element wise for all its sets of volumecells when there is a non-row node in this element
  for ( int k = 0; k < dis.NumMyColElements(); ++k )
  {
    DRT::Element* ele = dis.lColElement( k );
    int eid = ele->Id();
    GEO::CUT::ElementHandle * e = GetElement( eid );

    if(e!=NULL)
    {

      if( include_inner)
      {
        // get inside cell_sets connected within current element
        const std::vector<plain_volumecell_set> & ele_vc_sets_inside = e->GetVcSetsInside();
        std::vector<std::map<int,int> > & vcsets_nid_dofsetnumber_map_toComm_inside = e->Get_NodeDofsetMap_VcSets_Inside_forCommunication();

        int set_index=0;
        // decide for each set of connected volumecells, if communication is necessary
        for(std::vector<std::map<int,int> >::iterator set_it = vcsets_nid_dofsetnumber_map_toComm_inside.begin();
            set_it != vcsets_nid_dofsetnumber_map_toComm_inside.end();
            set_it++)
        {

          // does the current set contain dofset data to communicate
          if(set_it->size() > 0)
          {
            // communicate data for the first Volumecell in this set
            // REMARK: all cells contained in a set carry the same dofset information

            //first vc in set
            VolumeCell * cell = *(ele_vc_sets_inside[set_index].begin());

            if(cell == NULL) dserror("pointer to first Volumecell of set is NULL!");

            CreateParallelDofSetDataVC(parallel_dofSetData, eid, set_index, true, cell, *set_it);

          }

          set_index++;
        }
      }

      // standard case for outside elements
      {
        // get outside cell_sets connected within current element
        const std::vector<plain_volumecell_set> & ele_vc_sets_outside = e->GetVcSetsOutside();
        std::vector<std::map<int,int> > & vcsets_nid_dofsetnumber_map_toComm_outside = e->Get_NodeDofsetMap_VcSets_Outside_forCommunication();

        int set_index=0;
        // decide for each set of connected volumecells, if communication is necessary
        for(std::vector<std::map<int,int> >::iterator set_it = vcsets_nid_dofsetnumber_map_toComm_outside.begin();
            set_it != vcsets_nid_dofsetnumber_map_toComm_outside.end();
            set_it++)
        {

          // does the current set contain dofset data to communicate
          if(set_it->size() > 0)
          {
            // communicate data for the first Volumecell in this set
            // REMARK: all cells contained in a set carry the same dofset information

            //first vc in set
            VolumeCell * cell = *(ele_vc_sets_outside[set_index].begin());

            if(cell == NULL) dserror("pointer to first Volumecell of set is NULL!");

            CreateParallelDofSetDataVC(parallel_dofSetData, eid, set_index, false, cell, *set_it);

          }



          set_index++;
        }
      }

    }

  }// end col elements

}


/*--------------------------------------------------------------------------------------*
 | create parallel DofSetData for a volumecell that has to be communicated schott 03/12 |
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::ParentIntersection::CreateParallelDofSetDataVC(
    std::vector<Teuchos::RCP<DofSetData> >& parallel_dofSetData,
    int                           eid,
    int                           set_index,
    bool                          inside,
    VolumeCell *                  cell,
    std::map<int,int>&            node_dofset_map
    )
{
    if(node_dofset_map.size() > 0)
    {
      // get volumecell information
      // REMARK: identify volumecells using the volumecells (its facets) points
      plain_point_set cut_points; // use sets here such that common points of facets are not stored twice

      {
        // get all the facets points
        const plain_facet_set & facets = cell->Facets();

        for(plain_facet_set::const_iterator i = facets.begin(); i!=facets.end(); ++i)
        {
          Facet* f = *i;

          // decide which points has to be send!!
          // Points, CornerPoints, AllPoints
          const std::vector<Point*> & facetpoints = f->Points();

          std::copy( facetpoints.begin(), facetpoints.end(), std::inserter( cut_points, cut_points.begin() ) );
        }
      }

      LINALG::Matrix<3,1> coords(true);

      std::vector<LINALG::Matrix<3,1> > cut_points_coords(cut_points.size(), coords );

      for(plain_point_set::iterator p=cut_points.begin(); p!=cut_points.end(); ++p)
      {
        const int idx = std::distance(cut_points.begin(),p);
        LINALG::Matrix<3,1> & xyz = cut_points_coords[idx];

        std::copy( (*p)->X(), (*p)->X()+3, &xyz( 0, 0 ) );
      }


      // get the parent element Id
      // REMARK: for quadratic elements use the eid for the base element, not -1 for subelements

      // create dofset data for this volumecell for Communication
      parallel_dofSetData.push_back( Teuchos::rcp(new DofSetData( set_index, inside, cut_points_coords, eid, node_dofset_map) ) );

    }
    else dserror("communication for empty node-dofset map not necessary!");

}


/*--------------------------------------------------------------------------------------*
 | find cell sets around each node (especially for quadratic elements)     schott 03/12 |
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::ParentIntersection::FindNodalCellSets( bool include_inner,
                                                    std::set<int> & eids,
                                                    std::vector<int> & sourrounding_elements,
                                                    std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets_inside,
                                                    std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets_outside,
                                                    std::vector<plain_volumecell_set> & cell_sets_inside,
                                                    std::vector<plain_volumecell_set> & cell_sets_outside,
                                                    std::vector<plain_volumecell_set> & cell_sets )
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- FindNodalCellSets" );

  for ( std::vector<int>::iterator i=sourrounding_elements.begin(); i!=sourrounding_elements.end(); ++i )
  {
    int eid = *i;

    ElementHandle * e = GetElement(eid);

    const std::vector<plain_volumecell_set> & ele_vc_sets_inside   = e->GetVcSetsInside();
    const std::vector<plain_volumecell_set> & ele_vc_sets_outside  = e->GetVcSetsOutside();
//    e->VolumeCellSets( include_inner, ele_vc_sets_inside, ele_vc_sets_outside);

    // copy into cell_sets that collects all sets of adjacent elements
    if ( include_inner )
    {
      std::copy( ele_vc_sets_inside.begin(), ele_vc_sets_inside.end(), std::inserter( cell_sets, cell_sets.end() ) );
      std::copy( ele_vc_sets_inside.begin(), ele_vc_sets_inside.end(), std::inserter( cell_sets_inside, cell_sets_inside.end() ) );
    }

    std::copy( ele_vc_sets_outside.begin(), ele_vc_sets_outside.end(), std::inserter( cell_sets, cell_sets.end() ) );
    std::copy( ele_vc_sets_outside.begin(), ele_vc_sets_outside.end(), std::inserter( cell_sets_outside, cell_sets_outside.end() ) );


    if(   (ele_vc_sets_inside.size() >0 and include_inner)
        or ele_vc_sets_outside.size()>0)
    {
      eids.insert(eid); // no duplicates in std::set
    }

    const std::vector<Node*> & nodes = e->Nodes();


    for ( std::vector<Node*>::const_iterator n=nodes.begin();
        n!=nodes.end();
        ++n )
    {
      Node * node = *n;

      // call once for inside and once for outside
      {
        if(include_inner)
        {
          node->AssignNodalCellSet(ele_vc_sets_inside, nodal_cell_sets_inside);
        }

        node->AssignNodalCellSet(ele_vc_sets_outside, nodal_cell_sets_outside);

      }
    } // end loop over nodes of current sourrounding element


  }

}

/*--------------------------------------------------------------------------------------*
 | connect sets of volumecells for neighboring elements around a node      schott 03/12 |
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::ParentIntersection::ConnectNodalDOFSets(
    std::vector<Node *> &                     nodes,
    bool                                      include_inner,
    const DRT::Discretization&                dis,
    const std::vector<plain_volumecell_set> & connected_vc_sets, // connections of volumecells within one element connected via subelements
    std::vector<std::vector<int> > &          nodaldofset_vc_sets,
    std::vector<std::map<int,int> >&          vcsets_nid_dofsetnumber_map_toComm
)
{

  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- ConnectNodalDOFSets" );


  for(std::vector<plain_volumecell_set>::const_iterator s=connected_vc_sets.begin(); // connections within this element
      s!=connected_vc_sets.end();
      s++)
  {
    const plain_volumecell_set & cells = *s; // this is one connection of volumecells, connected via subelements, within one element

    std::vector<int > nds;

#ifdef PARALLEL
    // fill the map with nids, whose dofsets for the current set of volumecells has to filled by the nodes row proc
    // initialize the value (dofset_number with -1)
    std::map<int,int> nids_dofsetnumber_map_toComm;
#endif

    // find this plain_volumecell_set in dof_cellsets_ vector of each node
    {
      for ( std::vector<Node*>::iterator i=nodes.begin();
          i!=nodes.end();
          ++i )
      {
        //                Node * n = *i;
        //
        //                if( n->Id() >= 0) nds.push_back( n->DofSetNumberNEW( cells ) );
        //                else dserror("node with negative Id gets no dofnumber!");

        Node * n = *i;

        int nid = n->Id();

        if( nid >= 0)
        {
#ifdef PARALLEL
          DRT::Node* drt_node = dis.gNode(nid);

          // decide if the information for this cell has to be ordered from row-node or not
          //REMARK:
          if(drt_node->Owner() == dis.Comm().MyPID())
          {
            nds.push_back( n->DofSetNumberNEW( cells ) );
          }
          else
          {
            // insert the required pair of nid and unset dofsetnumber value (-1)
            nids_dofsetnumber_map_toComm[nid] = -1;

            // set dofset number to minus one, not a valid dofset number
            nds.push_back(-1);

          }
#else
          nds.push_back( n->DofSetNumberNEW( cells ) );
#endif
        }
        else dserror("node with negative Id gets no dofnumber!");

      }

    }

    vcsets_nid_dofsetnumber_map_toComm.push_back(nids_dofsetnumber_map_toComm);

    // set the nds vector for each volumecell of the current set
    for(plain_volumecell_set::const_iterator c=cells.begin(); c!=cells.end(); c++)
    {
      VolumeCell* cell = *c;
      cell->SetNodalDofSet(nds);
    }

    nodaldofset_vc_sets.push_back(nds);

  }
}


/*------------------------------------------------------------------------------------------------*
 * standard Cut routine for parallel XFSI, XFLUIDFLUID and Level set cut where dofsets and        *
 * node positions have to be parallelized                                            schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::ParentIntersection::Cut_Finalize( bool include_inner,
                                               INPAR::CUT::VCellGaussPts VCellgausstype,
                                               INPAR::CUT::BCellGaussPts BCellgausstype,
                                               bool tetcellsonly,
                                               bool screenoutput)
{

  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 6/6 --- Cut_Finalize" );

  if(myrank_==0 and screenoutput) IO::cout << "\t * 6/6 Cut_Finalize ...\t";

//  const double t_start = Teuchos::Time::wallTime();

  //----------------------------------------------------------

  Mesh & m = NormalMesh();

#ifdef DEBUGCUTLIBRARY
  include_inner = true; //Needed for debugging the volume
#endif

  if(VCellgausstype==INPAR::CUT::VCellGaussPts_Tessellation)
  {
    TEUCHOS_FUNC_TIME_MONITOR( "XFEM::FluidWizard::Cut::Tessellation" );
    m.CreateIntegrationCells( 0 ,tetcellsonly ); // boundary cells will be created within TetMesh.CreateElementTets
    //m.RemoveEmptyVolumeCells();

//    // Test and Debug:
    m.TestElementVolume( true, VCellgausstype );
    if(myrank_==0 and screenoutput) IO::cout << "\n\t *     TestElementVolume ...";
    m.TestFacetArea();
    if(myrank_==0 and screenoutput) IO::cout << "\n\t *     TestFacetArea ...";

#ifdef DEBUGCUTLIBRARY
    //m.TestVolumeSurface(); //Broken test, needs to be fixed for proper usage.
    m.TestFacetArea();
    m.TestElementVolume( true, VCellgausstype );
    DebugCut(m);
#endif
  }
  else if(VCellgausstype==INPAR::CUT::VCellGaussPts_MomentFitting)
  {
    TEUCHOS_FUNC_TIME_MONITOR( "XFEM::FluidWizard::Cut::MomentFitting" );
    m.MomentFitGaussWeights(include_inner, BCellgausstype);
#ifdef DEBUGCUTLIBRARY
    m.TestFacetArea();
#endif
  }
  else if(VCellgausstype==INPAR::CUT::VCellGaussPts_DirectDivergence)
  {
    TEUCHOS_FUNC_TIME_MONITOR( "XFEM::FluidWizard::Cut::DirectDivergence" );

#ifndef LOCAL
    if(myrank_==0)
      std::cout << "\t DirectDivergence Cut with Global configuration!" << std::endl;
#endif
    m.DirectDivergenceGaussRule(include_inner, BCellgausstype);

//    m.TestElementVolume( true, VCellgausstype );
//    m.TestFacetArea();

#ifdef DEBUGCUTLIBRARY
    m.TestElementVolume( true, VCellgausstype );
    m.TestFacetArea();
#endif
  }
  else
    dserror("Undefined option of volumecell gauss points generation");


  //----------------------------------------------------------
//
//  const double t_diff = Teuchos::Time::wallTime()-t_start;
//  if ( myrank_ == 0 and screenoutput)
//  {
//    IO::cout << " Success (" << t_diff  <<  " secs)" << IO::endl;
//  }

#if(0)
  std::cout << "\n XFEM::FluidWizard::Quadrature construction time = " << t_diff <<"\n";
#endif


#ifdef DEBUGCUTLIBRARY
  std::stringstream str;
  str << "fullfinalcutoutput." << myrank_ << ".pos";
  NormalMesh().DumpGmsh(str.str());
#endif
}

/*--------------------------------------------------------------------------------------*
 * get the node based on node id
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Node * GEO::CUT::ParentIntersection::GetNode( int nid ) const
{
  return mesh_.GetNode( nid );
}

/*--------------------------------------------------------------------------------------*
 * get the mesh's side based on node ids and return the side
 *-------------------------------------------------------------------------------------*/
GEO::CUT::SideHandle * GEO::CUT::ParentIntersection::GetSide( std::vector<int>& nodeids ) const
{
  return mesh_.GetSide( nodeids );
}

/*--------------------------------------------------------------------------------------*
 * get the mesh's side based on side id and return the sidehandle
 *-------------------------------------------------------------------------------------*/
GEO::CUT::SideHandle * GEO::CUT::ParentIntersection::GetSide( int sid ) const
{
  return mesh_.GetSide( sid );
}

/*--------------------------------------------------------------------------------------*
 * get the mesh's element based on element id
 *-------------------------------------------------------------------------------------*/
GEO::CUT::ElementHandle * GEO::CUT::ParentIntersection::GetElement( int eid ) const
{
  return mesh_.GetElement( eid );
}

/*--------------------------------------------------------------------------------------*
 * print cell statistics
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::ParentIntersection::PrintCellStats()
{
  NormalMesh().PrintCellStats();
}

/*--------------------------------------------------------------------------------------*
 * write gmsh debug output for nodal cell sets
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::ParentIntersection::DumpGmshNodalCellSet( std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets, const DRT::Discretization & dis )
{
    std::string filename = "cut_test"; //DRT::Problem::Instance()->OutputControlFile()->FileName();
    std::stringstream str;
    str << filename
        << "CUT_NodalCellSet."
        << dis.Comm().MyPID()
        << ".pos";


    std::string name = str.str();

    std::ofstream file( name.c_str() );



    // Gmsh output for the sets of volumecells (connected within a global element) that are assigned to a node
    // all the cells of a set get the node id of the node they are assigned to

    file << "View \"NodalCellSet\" {\n";

    for(std::map<Node*, std::vector<plain_volumecell_set> >::iterator i=nodal_cell_sets.begin(); i!=nodal_cell_sets.end(); i++)
    {
        Node* n= i->first;

        int nid = n->Id();

        std::vector<plain_volumecell_set> & sets = i->second;

        for(std::vector<plain_volumecell_set>::iterator s=sets.begin(); s!=sets.end(); s++)
        {
            const plain_volumecell_set & volumes = *s;

            for ( plain_volumecell_set::const_iterator i=volumes.begin(); i!=volumes.end(); ++i )
            {
                VolumeCell * vc = *i;

                const plain_integrationcell_set & integrationcells = vc->IntegrationCells();
                for ( plain_integrationcell_set::const_iterator i=integrationcells.begin();
                    i!=integrationcells.end();
                    ++i )
                {
                    IntegrationCell * ic = *i;
                    ic->DumpGmsh( file, &nid );
                }
            }

        }
    }

    file << "};\n";



    // Gmsh output, additional information (node Ids)

    file << "View \"NodeID\" {\n";

    for(std::map<Node*, std::vector<plain_volumecell_set> >::iterator i=nodal_cell_sets.begin(); i!=nodal_cell_sets.end(); i++)
    {
        Node* n= i->first;

        int nid = n->Id();


        Point * p = n->point();
        const double * x = p->X();

        // output just for real nodes of elements, not for shadow nodes
        file << "SP(" << x[0] << "," << x[1] << "," << x[2] << "){" << nid << "};\n";

    }

    file << "};\n";

}

/*--------------------------------------------------------------------------------------*
 * write gmsh debug output for CellSets
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::ParentIntersection::DumpGmshCellSets( std::vector<plain_volumecell_set> & cell_sets, const DRT::Discretization & dis )
{
    std::string filename = "cut_test"; //DRT::Problem::Instance()->OutputControlFile()->FileName();
    std::stringstream str;
    str << filename
        << "CUT_CellSets."
        << dis.Comm().MyPID()
        << ".pos";


    std::string name = str.str();

    std::ofstream file( name.c_str() );


    // Gmsh output for all sets of connected volumecells (connected within a global element) that are assigned to a node

    plain_volumecell_set cells;

    // get cell_sets as a plain_volume_set
    for(std::vector<plain_volumecell_set>::iterator i=cell_sets.begin(); i!=cell_sets.end(); i++)
    {
      std::copy((*i).begin(), (*i).end(), std::inserter( cells, cells.begin() ) );
    }

    file << "View \"CellSet\" {\n";
    int count=0;

    for ( plain_volumecell_set::const_iterator i=cells.begin(); i!=cells.end(); ++i )
    {
        count++;
        VolumeCell * vc = *i;

        const plain_integrationcell_set & integrationcells = vc->IntegrationCells();
        for ( plain_integrationcell_set::const_iterator i=integrationcells.begin();
              i!=integrationcells.end();
              ++i )
        {
            IntegrationCell * ic = *i;
            ic->DumpGmsh( file, &count );
        }
    }


    file << "};\n";
}


/*--------------------------------------------------------------------------------------*
 * write gmsh cut output for number of dofsets and the connected vc sets
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::ParentIntersection::DumpGmshNumDOFSets(std::string filename, bool include_inner, const DRT::Discretization & dis )
{
    std::stringstream str;
    str << filename
        << ".CUT_NumDOFSets."
        << dis.Comm().MyPID()
        << ".pos";


    Mesh & m = NormalMesh();



    std::string name = str.str();

    std::ofstream file( name.c_str() );


    // Gmsh output for all sets of connected volumecells (connected within a global element) separated for inside and outside
    // each set gets its own number ( inside (negative) , outside(positive) )

    file << "View \"ConnectedVcSets\" {\n";
    int count_inside=-1;
    int count_outside=0;

    int num_row_ele = dis.NumMyRowElements();

    for( int lid=0; lid< num_row_ele; lid++) //std::set<int>::iterator i= eids.begin(); i!= eids.end(); i++)
    {
      DRT::Element * e = dis.lRowElement(lid);
        int eid = e->Id();

        ElementHandle* eh = GetElement(eid);

        if(eh != NULL)
        {
        // get inside and outside cell_sets connected within current element
        const std::vector<plain_volumecell_set> & ele_vc_sets_inside  = eh->GetVcSetsInside();
        const std::vector<plain_volumecell_set> & ele_vc_sets_outside = eh->GetVcSetsOutside();


        for ( std::vector<plain_volumecell_set>::const_iterator i=ele_vc_sets_outside.begin();
              i!=ele_vc_sets_outside.end();
              ++i )
        {
            plain_volumecell_set volumes = *i;

            for ( plain_volumecell_set::const_iterator i=volumes.begin(); i!=volumes.end(); ++i )
            {
                VolumeCell * vc = *i;

                const plain_integrationcell_set & integrationcells = vc->IntegrationCells();
                for ( plain_integrationcell_set::const_iterator i=integrationcells.begin();
                      i!=integrationcells.end();
                      ++i )
                {
                    IntegrationCell * ic = *i;
                    ic->DumpGmsh( file, &count_outside );
                }
            }
            count_outside +=1;
        }

        // for inside cells
        if(include_inner)
        {
            for ( std::vector<plain_volumecell_set>::const_iterator i=ele_vc_sets_inside.begin();
                  i!=ele_vc_sets_inside.end();
                  ++i )
            {
                const plain_volumecell_set & volumes = *i;

                for ( plain_volumecell_set::const_iterator i=volumes.begin(); i!=volumes.end(); ++i )
                {
                    VolumeCell * vc = *i;

                    const plain_integrationcell_set & integrationcells = vc->IntegrationCells();
                    for ( plain_integrationcell_set::const_iterator i=integrationcells.begin();
                          i!=integrationcells.end();
                          ++i )
                    {
                        IntegrationCell * ic = *i;
                        ic->DumpGmsh( file, &count_inside );
                    }
                }
                count_inside -=1;
            }
        }
        }
    }

    file << "};\n";



    // Gmsh output for all dof sets of one node (connected via adjacent elements of this node)
    // each set gets its own number ( inside (negative) , outside(positive) )


//    file << "View \"DofSets for special node\" {\n";
//
//
//    Node* n = m.GetNode( nid );
//    Node* n = m.GetNode( 70 );
//    std::vector<std::set<plain_volumecell_set> > dof_cellsets = n->DofCellSets();
//
//    int count=0;
////    int count_vc = 0;
//    for(std::vector<std::set<plain_volumecell_set> >::iterator i=dof_cellsets.begin(); i!=dof_cellsets.end(); i++)
//    {
//        std::set<plain_volumecell_set> & cellset = *i;
//
//        for(std::set<plain_volumecell_set>::iterator vc_set=cellset.begin(); vc_set!=cellset.end(); vc_set++)
//        {
//            const plain_volumecell_set & volumes = *vc_set;
//
//            for ( plain_volumecell_set::const_iterator i=volumes.begin(); i!=volumes.end(); ++i )
//            { //count_vc++;
//                VolumeCell * vc = *i;
//
//                const plain_integrationcell_set & integrationcells = vc->IntegrationCells();
//                for ( plain_integrationcell_set::const_iterator i=integrationcells.begin();
//                      i!=integrationcells.end();
//                      ++i )
//                {
//                    IntegrationCell * ic = *i;
//                    ic->DumpGmsh( file, &count );
//
//                }
//
//            }
//        }
//        count++;
//    }
//    file << "};\n";

    // nodes used for CUT std::map<node->ID, Node>, shadow nodes have ID<0
    // print the dofsets just for the row nodes
    std::map<int, Node* > nodes;
    m.GetNodeMap( nodes );

    file << "View \"NumDofSets\" {\n";
    for ( std::map<int, Node* >::iterator i=nodes.begin(); i!=nodes.end(); ++i )
    {
      int nid = i->first;

      if(nid >= 0)
      {
        if(dis.NodeRowMap()->LID(nid) == -1) continue; // non-local row node

        Node * n = i->second;
        Point * p = n->point();
        const double * x = p->X();

        // output just for real nodes of elements, not for shadow nodes
        if( n->Id() >= 0) file << "SP(" << x[0] << "," << x[1] << "," << x[2] << "){" << n->NumDofSets() << "};\n";
      }

    }
    file << "};\n";
}

/*--------------------------------------------------------------------------------------*
 * write gmsh output for volumecells
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::ParentIntersection::DumpGmshVolumeCells( std::string name, bool include_inner )
{
  NormalMesh().DumpGmshVolumeCells( name, include_inner );
}

/*--------------------------------------------------------------------------------------*
 * write gmsh output for volumecells
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::ParentIntersection::DumpGmshIntegrationCells( std::string name )
{
  NormalMesh().DumpGmshIntegrationCells( name );
}

/*--------------------------------------------------------------------------------------*
 * write gmsh output for volumecells
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::ParentIntersection::DumpGmshVolumeCells( std::string name )
{
  NormalMesh().DumpGmshVolumeCells( name );
}


/*--------------------------------------------------------------------------------------*
 * status
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::ParentIntersection::Status(INPAR::CUT::VCellGaussPts gausstype)
{
#ifdef DEBUG
  NormalMesh().Status();

#ifdef DEBUGCUTLIBRARY
  NormalMesh().DumpGmsh( "mesh.pos" );

  if(gausstype==INPAR::CUT::VCellGaussPts_Tessellation)
  {
    DumpGmshIntegrationCells( "integrationcells.pos" );
    DumpGmshVolumeCells("volumecells.pos");
  }
  else
    DumpGmshVolumeCells("volumecells.pos");
#endif
#endif
}

/*
 * If the cut is done by Tesselation then the volume of the volumecells,
 *   can be compared to DirectDivergence.
 */
void GEO::CUT::ParentIntersection::DebugCut(Mesh &m)
{
  // TEST IF VOLUMES PREDICTED OF CELLS ARE SAME:
   //###########################################################
   std::vector<double> tessVol,dirDivVol;
   double diff_tol = 10000*BASICTOL;//1e-20; //How sharp to test for.

   const std::list<Teuchos::RCP<VolumeCell> >  & other_cells = m.VolumeCells();

   for ( std::list<Teuchos::RCP<VolumeCell> >::const_iterator i=other_cells.begin();
       i!=other_cells.end();
       ++i )
   {
     GEO::CUT::VolumeCell * vc = &**i;
     tessVol.push_back(vc->Volume());

     // TEST THAT TESSELATION DOES NOT HAVE EMPTY VolumeCells.
     //++++++++++++++++++++
//     const GEO::CUT::plain_integrationcell_set & integrationcells = vc->IntegrationCells();

//     if(integrationcells.size()==0)
//     {
//       throw std::runtime_error("VolumeCell contains 0 integration cells.");
//     }
     //++++++++++++++++++++
   }
   //------------------------------------------------------

   //Cut with DirectDivergence as well
   for ( std::list<Teuchos::RCP<VolumeCell> >::const_iterator i=other_cells.begin();
       i!=other_cells.end();
       ++i )
   {
     GEO::CUT::VolumeCell * vc = &**i;
     vc->DirectDivergenceGaussRule(vc->ParentElement(),m,true);
     dirDivVol.push_back(vc->Volume());
   }

   // Test Direct Divergence against Tesselation
   //------------------------------------------------------
   bool error = false;
   double sumtessvol = 1.0;
   for(unsigned i=0;i<tessVol.size();i++)
     sumtessvol+=tessVol[i];

   //Estimate the total volume of a element...
   sumtessvol = sumtessvol/(tessVol.size()/2.0);

   std::cout << "Tolerance: " << diff_tol*sumtessvol << std::endl;
   std::cout<<"the volumes predicted by\n Tessellation \t DirectDivergence\n";
   for(unsigned i=0;i<tessVol.size();i++)
   {
     std::cout<< std::setprecision(17) << std::scientific << tessVol[i]<<"\t"<<dirDivVol[i];//<<"\n";
     if( fabs(tessVol[i]-dirDivVol[i]) > diff_tol*(sumtessvol)) // (1e-8) diff_tol*(sumtessvol)) // Is the relative difference greater than the Tolerance?
     {
       std::cout << "\t difference: " << fabs(tessVol[i]-dirDivVol[i]);// << "\n";
       error = true;
     }
     else if(std::isnan(tessVol[i]) or std::isnan(dirDivVol[i]))
       error = true;
     std::cout << "\n";
   }
   if(error)
   {
     throw std::runtime_error("Volume predicted by one of the methods is wrong.");
     //std::cout << "Volume predicted by one of the methods is wrong." << std::endl;
   }
   //------------------------------------------------------
   //###########################################################
}

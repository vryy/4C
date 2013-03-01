/*!-----------------------------------------------------------------------------------------------*
\file cut_meshintersection.cpp

\brief provides the basic functionality for cutting a mesh

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "cut_integrationcell.H"
#include "cut_volumecell.H"

#include "cut_meshintersection.H"

#include "../drt_fluid/xfluid_defines.H"

#include <Teuchos_TimeMonitor.hpp>

/*-----------------------------------------------------------------------------------------*
 * add this background element if it falls within the bounding box of cut mesh
 * If it is not within BB, this element is never cut
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::ElementHandle * GEO::CUT::MeshIntersection::AddElement( int eid,
                                                                  const std::vector<int> & nids,
                                                                  const Epetra_SerialDenseMatrix & xyz,
                                                                  DRT::Element::DiscretizationType distype )
{
  for ( std::vector<Teuchos::RCP<MeshHandle> >::iterator i=cut_mesh_.begin();
        i!=cut_mesh_.end();
        ++i )
  {
    MeshHandle & cut_mesh_handle = **i;
    Mesh & cut_mesh = cut_mesh_handle.LinearMesh();

    // element is to be added only when it falls within bounding box
    // generated over cut mesh. otherwise this element is never cut
    if ( cut_mesh.WithinBB( xyz ) )
    {
      int numnode = nids.size();
      if ( numnode != xyz.N() )
      {
        throw std::runtime_error( "node coordiante number mismatch" );
      }

      // make sure all nodes are there
      for ( int i=0; i<numnode; ++i )
      {
        NormalMesh().GetNode( nids[i], &xyz( 0, i ) );
//         if ( n==NULL )
//         {
//           // if there is no node with that id but a node at the given
//           // location, the element is illegal and cannot be created
//           return;
//         }
      }

      // create element
      return mesh_.CreateElement( eid, nids, distype );
    }
  }
  return NULL;
}

/*-----------------------------------------------------------------------------------------*
 * add a side of the cut mesh and return the sidehandle (e.g. quadratic sidehandle for quadratic sides)
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::SideHandle * GEO::CUT::MeshIntersection::AddCutSide( int sid,
                                                               const std::vector<int> & nids,
                                                               DRT::Element::DiscretizationType distype,
                                                               int mi )
{
  // create side
  return cut_mesh_[mi]->CreateSide( sid, nids, distype );
}

/*-----------------------------------------------------------------------------------------*
 * add a side of the cut mesh and return the sidehandle (e.g. quadratic sidehandle for quadratic sides)
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::SideHandle * GEO::CUT::MeshIntersection::AddCutSide( int sid,
                                                               const std::vector<int> & nids,
                                                               const Epetra_SerialDenseMatrix & xyz,
                                                               DRT::Element::DiscretizationType distype,
                                                               int mi )
{
  Mesh & cut_mesh = CutMesh( mi );

  int numnode = nids.size();
  if ( numnode != xyz.N() )
  {
    throw std::runtime_error( "node coordiante number mismatch" );
  }

//   PointSet nodalpoints;

  // make sure all nodes are there
  for ( int i=0; i<numnode; ++i )
  {
    cut_mesh.GetNode( nids[i], &xyz( 0, i ) );
//     nodalpoints.insert( n->point() );
//     if ( n==NULL )
//     {
//       // if there is no node with that id but a node at the given location,
//       // the side is illegal and cannot be created
//       return;
//     }
  }

//   // do not create degenerated cut sides
//   if ( nodalpoints.size() < nids.size() )
//     return;

  // create side
  return cut_mesh_[mi]->CreateSide( sid, nids, distype );
}

/*------------------------------------------------------------------------------------------------*
 * standard Cut routine for two phase flow and combustion where dofsets and node positions        *
 * have not to be computed, standard cut for cut_est                                 schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::Cut( bool include_inner, std::string VCellgausstype, std::string BCellgausstype )
{
  Status();

  // cut the mesh and create cutlines, facets, volumecells
  Cut_Mesh( include_inner );

  // determine inside-outside position and dofset-data, parallel communication if required
  Cut_Positions_Dofsets( include_inner );

  // create integration points and/or subtetrahedralization
  Cut_Finalize( include_inner, VCellgausstype, BCellgausstype);

  // DumpGmshVolumeCells("CUT_vc", true);
  // DumpGmshIntegrationCells("CUT_intcells");


  Status(VCellgausstype);
}




/*------------------------------------------------------------------------------------------------*
 * standard Cut routine for parallel XFSI and XFLUIDFLUID where dofsets and node positions        *
 * have to be parallelized                                                           schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::Cut_Mesh( bool include_inner)
{

  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 1/3 --- Cut_Mesh" );


  if(myrank_==0) IO::cout << "\n\t ... 1/3 Cut_Mesh";

  const double t_start = Teuchos::Time::wallTime();

  //----------------------------------------------------------

  Mesh & m = NormalMesh();

  plain_element_set elements_done;

  // loop cut sides and cut against elements at the same position in space
  for ( std::vector<Teuchos::RCP<MeshHandle> >::iterator i=cut_mesh_.begin();
        i!=cut_mesh_.end();
        ++i )
  {
    MeshHandle & cut_mesh_handle = **i;
    Mesh & cut_mesh = cut_mesh_handle.LinearMesh();
    cut_mesh.Cut( m, elements_done, 0 );
  }

  m.MakeCutLines();
  m.MakeFacets();
  m.MakeVolumeCells();

  //----------------------------------------------------------

  const double t_diff = Teuchos::Time::wallTime()-t_start;
  if ( myrank_ == 0 )
  {
    IO::cout << " ... Success (" << t_diff  <<  " secs)\n";
  }

}

/*------------------------------------------------------------------------------------------------*
 * Routine for deciding the inside-outside position. This creates the dofset data,                *
 * also in parallel                                                                  schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::Cut_Positions_Dofsets( bool include_inner )
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 2/3 --- Cut_Positions_Dofsets (serial)" );

  if(myrank_==0) IO::cout << "\n\t ... 2/3 Cut_Positions_Dofsets";

  const double t_start = Teuchos::Time::wallTime();

  //----------------------------------------------------------


  Mesh & m = NormalMesh();

  if ( options_.FindPositions() )
  {

    // find inside and outside positions of nodes, corresponding facets and volumecells using a DFS-algorithm
    m.FindNodePositions();

    // find the positions for all remaining facets ( and points, volumecells)
    m.FindFacetPositions();

    // find number and connection of dofsets at nodes from cut volumes, also in parallel
    m.FindNodalDOFSets( include_inner );

  }


  //----------------------------------------------------------

   const double t_diff = Teuchos::Time::wallTime()-t_start;
   if ( myrank_ == 0 )
   {
     IO::cout << " ... Success (" << t_diff  <<  " secs)";
   }
}



/*------------------------------------------------------------------------------------------------*
 * standard Cut routine for parallel XFSI and XFLUIDFLUID where dofsets and node positions        *
 * have to be parallelized                                                           schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::Cut_Finalize( bool include_inner, std::string VCellgausstype, std::string BCellgausstype)
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 3/3 --- Cut_Finalize" );

  if(myrank_==0) IO::cout << "\n\t ... 3/3 Cut_Finalize";

  const double t_start = Teuchos::Time::wallTime();

  //----------------------------------------------------------

  Mesh & m = NormalMesh();

  if(VCellgausstype=="Tessellation")
  {
    TEUCHOS_FUNC_TIME_MONITOR( "XFEM::FluidWizard::Cut::Tessellation" );
    m.CreateIntegrationCells( 0, false ); // boundary cells will be created within TetMesh.CreateElementTets
    //m.RemoveEmptyVolumeCells();

#ifdef DEBUGCUTLIBRARY
    //m.TestVolumeSurface();
    m.TestFacetArea();
#endif
    m.SimplifyIntegrationCells();

#ifdef DEBUGCUTLIBRARY
    m.TestElementVolume( true );
#endif
  }
  else if(VCellgausstype=="MomentFitting")
  {
    TEUCHOS_FUNC_TIME_MONITOR( "XFEM::FluidWizard::Cut::MomentFitting" );
    m.MomentFitGaussWeights(include_inner, BCellgausstype);
  }
  else if(VCellgausstype=="DirectDivergence")
  {
    TEUCHOS_FUNC_TIME_MONITOR( "XFEM::FluidWizard::Cut::DirectDivergence" );
    m.DirectDivergenceGaussRule(include_inner, BCellgausstype);
  }
  else
    dserror("Undefined option of volumecell gauss points generation");


  //----------------------------------------------------------

  const double t_diff = Teuchos::Time::wallTime()-t_start;
  if ( myrank_ == 0 )
  {
    IO::cout << " ... Success (" << t_diff  <<  " secs)\n";
  }

#if(0)
  std::cout << "\n XFEM::FluidWizard::Quadrature construction time = " << t_diff <<"\n";
#endif
}


/*------------------------------------------------------------------------------------------------*
 * Create nodal dofset sets within the parallel cut framework
 *------------------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::CreateNodalDofSetNEW( bool include_inner, DRT::Discretization& dis)
{

    std::set<int> eids; // eids of elements that are involved in CUT and include ele_vc_set_inside/outside (no duplicates!)

    Mesh & m = NormalMesh();

    // nodes used for CUT map<node->ID, Node>, shadow nodes have ID<0
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

        std::map<int, ElementHandle*> sourrounding_elements;

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
                    sourrounding_elements.insert(std::pair<int, ElementHandle*>(adj_eid, e));
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
            if(include_inner)
            {
                n->FindDOFSetsNEW( nodal_cell_sets_inside, cell_sets_inside);
            }

            n->FindDOFSetsNEW( nodal_cell_sets_outside, cell_sets_outside);

            // sort the dofsets for this node after FindDOFSetsNEW
            n->SortDOFCellSets();

        } // end if n_gid >= 0

    } // end loop over nodes


    //===============
    // STEP 2: for each element that contains volumecell_sets (connections via subelements...),
    // all nodes of this element have to now the dofset_number for each set of volumecells
    //===============
    for( std::set<int>::iterator i= eids.begin(); i!= eids.end(); i++)
    {
        int eid = *i;

        // get the nodes of this element
        // get the element via discret
        DRT::Element * e = dis.gElement(eid);

        if( e == NULL) dserror(" element not found, this should not be! ");

        std::vector<Node * > nodes;

        // get the nodes of this element
        int numnode = e->NumNode();
        const int* nids = e->NodeIds();

        for(int i=0; i< numnode; i++)
        {
            Node * node = GetNode(nids[i]);

            if(node == NULL) dserror("node not found!");

            nodes.push_back( node );
        }

        if((int)nodes.size() != numnode) dserror("number of nodes not equal to the required number of elemental nodes");

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
          ConnectNodalDOFSets(nodes,
                              include_inner,
                              dis,
                              ele_vc_sets_inside,
                              nodaldofset_vc_sets_inside,
                              vcsets_nid_dofsetnumber_map_toComm_inside);
        }

        ConnectNodalDOFSets(nodes,
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
void GEO::CUT::MeshIntersection::FillParallelDofSetData(RCP<std::vector<DofSetData> > parallel_dofSetData,
                                                        DRT::Discretization& dis)
{

  //TODO: set bool include_inner
  bool include_inner = false;

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
void GEO::CUT::MeshIntersection::CreateParallelDofSetDataVC( RCP<std::vector<DofSetData> > parallel_dofSetData,
                                                             int                           eid,
                                                             int                           set_index,
                                                             bool                          inside,
                                                             VolumeCell *                  cell,
                                                             std::map<int,int>&            node_dofset_map)
{

    if(node_dofset_map.size() > 0)
    {
      // get volumcell information
      // REMARK: identify volumecells using the volumecells (its facets) points
      std::vector<Point*> cut_points;
      {
        // get all the facets points
        const plain_facet_set facets = cell->Facets();

        for(plain_facet_set::const_iterator i = facets.begin(); i!=facets.end(); ++i)
        {
          Facet* f = *i;

          // decide which points has to be send!!
          // Points, CornerPoints, AllPoints
          std::vector<Point*> facetpoints = f->Points();

          std::copy( facetpoints.begin(), facetpoints.end(), std::inserter( cut_points, cut_points.begin() ) );
        }
      }

      std::vector<LINALG::Matrix<3,1> > cut_points_coords_;

      for(std::vector<Point*>::iterator p=cut_points.begin(); p!=cut_points.end(); ++p)
      {
        LINALG::Matrix<3,1> coords(true);
        const double * x = (*p)->X();
        for(int i=0; i<3; i++)
        {
          coords(i) = x[i];
        }
        cut_points_coords_.push_back(coords);
      }



      // get the parent element Id
//      int peid = cell->ParentElement()->Id();
      // REMARK: for quadratic elements use the eid for the base element, not -1 for subelements
      int peid = eid;

      // create dofset data for this volumecell for Communication
      parallel_dofSetData->push_back( DofSetData( set_index, inside, cut_points_coords_, peid, node_dofset_map) );

    }
    else dserror("communication for empty node-dofset map not necessary!");

}


/*--------------------------------------------------------------------------------------*
 | find cell sets around each node (especially for quadratic elements)     schott 03/12 |
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::FindNodalCellSets( bool include_inner,
                                                    std::set<int> & eids,
                                                    std::map<int, ElementHandle*> & sourrounding_elements,
                                                    std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets_inside,
                                                    std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets_outside,
                                                    std::vector<plain_volumecell_set> & cell_sets_inside,
                                                    std::vector<plain_volumecell_set> & cell_sets_outside,
                                                    std::vector<plain_volumecell_set> & cell_sets )
{

    for ( std::map<int,ElementHandle*>::iterator i=sourrounding_elements.begin(); i!=sourrounding_elements.end(); ++i )
    {
        ElementHandle * e = i->second;

        std::vector<plain_volumecell_set> ele_vc_sets_inside;
        std::vector<plain_volumecell_set> ele_vc_sets_outside;
        e->VolumeCellSets( include_inner, ele_vc_sets_inside, ele_vc_sets_outside);

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
            eids.insert(i->first); // no duplicates in std::set
        }

        const std::vector<Node*> & nodes = e->Nodes();


        for ( std::vector<Node*>::const_iterator i=nodes.begin();
              i!=nodes.end();
              ++i )
        {
            Node * n = *i;

            // call once for inside and once for outside
            {
                if(include_inner)
                {
                    n->AssignNodalCellSet(ele_vc_sets_inside, nodal_cell_sets_inside);
                }

                n->AssignNodalCellSet(ele_vc_sets_outside, nodal_cell_sets_outside);

            }


        } // end loop over nodes of current sourrounding element


    }

}

/*--------------------------------------------------------------------------------------*
 | connect sets of volumecells for neighboring elements around a node      schott 03/12 |
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::ConnectNodalDOFSets( std::vector<Node *> &                     nodes,
                                                      bool                                      include_inner,
                                                      DRT::Discretization&                      dis,
                                                      const std::vector<plain_volumecell_set> & connected_vc_sets,
                                                      std::vector<std::vector<int> > &          nodaldofset_vc_sets,
                                                      std::vector<std::map<int,int> >&          vcsets_nid_dofsetnumber_map_toComm)
{

    for(std::vector<plain_volumecell_set>::const_iterator s=connected_vc_sets.begin();
        s!=connected_vc_sets.end();
        s++)
    {
        plain_volumecell_set cells = *s; // this is one connection of volumecells, connected via subelements, within one element

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

              DRT::Node* drt_node = dis.gNode(nid);

              if( nid >= 0)
              {
#ifdef PARALLEL
                // decide if the information for this cell has to be ordered from row-node or not
                //REMARK:
                if(drt_node->Owner() == dis.Comm().MyPID())
                {
                  nds.push_back( n->DofSetNumberNEW( cells ) );
                }
                else
                {
                  // insert the required pair of nid and unset dofsetnumber value (-1)
                  nids_dofsetnumber_map_toComm.insert(std::pair<int,int>(nid,-1));

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
        for(plain_volumecell_set::iterator c=cells.begin(); c!=cells.end(); c++)
        {
          VolumeCell* cell = *c;
          cell->SetNodalDofSet(nds);
        }

        nodaldofset_vc_sets.push_back(nds);

    }
}


/*--------------------------------------------------------------------------------------*
 * get the node based on node id
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Node * GEO::CUT::MeshIntersection::GetNode( int nid ) const
{
  return mesh_.GetNode( nid );
}

/*--------------------------------------------------------------------------------------*
 * get the mesh's side based on node ids and return the side
 *-------------------------------------------------------------------------------------*/
GEO::CUT::Side * GEO::CUT::MeshIntersection::GetSide( std::vector<int>& nodeids ) const
{
  return mesh_.GetSide( nodeids );
}

/*--------------------------------------------------------------------------------------*
 * get the mesh's side based on side id and return the sidehandle
 *-------------------------------------------------------------------------------------*/
GEO::CUT::SideHandle * GEO::CUT::MeshIntersection::GetSide( int sid ) const
{
  return mesh_.GetSide( sid );
}

/*--------------------------------------------------------------------------------------*
 * get the mesh's element based on element id
 *-------------------------------------------------------------------------------------*/
GEO::CUT::ElementHandle * GEO::CUT::MeshIntersection::GetElement( int eid ) const
{
  return mesh_.GetElement( eid );
}

/*--------------------------------------------------------------------------------------*
 * get the cut mesh's side based on side id
 *-------------------------------------------------------------------------------------*/
GEO::CUT::SideHandle * GEO::CUT::MeshIntersection::GetCutSide( int sid, int mi ) const
{
  return cut_mesh_[mi]->GetSide( sid );
}

/*--------------------------------------------------------------------------------------*
 * print cell statistics
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::PrintCellStats()
{
  NormalMesh().PrintCellStats();
}

void GEO::CUT::MeshIntersection::Status(std::string gausstype)
{
#ifdef DEBUG
  NormalMesh().Status();
  for ( std::vector<Teuchos::RCP<MeshHandle> >::iterator i=cut_mesh_.begin();
        i!=cut_mesh_.end();
        ++i )
  {
    MeshHandle & cut_mesh_handle = **i;
    Mesh & cut_mesh = cut_mesh_handle.LinearMesh();
    cut_mesh.Status();
  }

#ifdef DEBUGCUTLIBRARY
  NormalMesh().DumpGmsh( "mesh.pos" );
  int count = 0;
  for ( std::vector<Teuchos::RCP<MeshHandle> >::iterator i=cut_mesh_.begin();
        i!=cut_mesh_.end();
        ++i )
  {
    MeshHandle & cut_mesh_handle = **i;
    Mesh & cut_mesh = cut_mesh_handle.LinearMesh();
    std::stringstream str;
    str << "cut_mesh" << count << ".pos";
    cut_mesh.DumpGmsh( str.str().c_str() );
    count++;
  }

  //NormalMesh().DumpGmshVolumeCells( "volumecells" );
  if(gausstype=="Tessellation")
  {
    DumpGmshIntegrationCells( "integrationcells.pos" );
    DumpGmshVolumeCells("volumecells.pos");
  }
  else if(gausstype=="MomentFitting")
    DumpGmshVolumeCells("volumecells.pos");
#endif
#endif
}


/*--------------------------------------------------------------------------------------*
 * write gmsh debug output for nodal cell sets
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::DumpGmshNodalCellSet( std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets, DRT::Discretization & dis )
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
void GEO::CUT::MeshIntersection::DumpGmshCellSets( std::vector<plain_volumecell_set> & cell_sets,  DRT::Discretization & dis )
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
void GEO::CUT::MeshIntersection::DumpGmshNumDOFSets(std::string filename, bool include_inner, DRT::Discretization & dis )
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
        std::vector<plain_volumecell_set> ele_vc_sets_inside;
        std::vector<plain_volumecell_set> ele_vc_sets_outside;

        eh->VolumeCellSets( include_inner, ele_vc_sets_inside, ele_vc_sets_outside);


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
////	  int count_vc = 0;
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

    // nodes used for CUT map<node->ID, Node>, shadow nodes have ID<0
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
void GEO::CUT::MeshIntersection::DumpGmshVolumeCells( std::string name, bool include_inner )
{
  NormalMesh().DumpGmshVolumeCells( name, include_inner );
}

/*--------------------------------------------------------------------------------------*
 * write gmsh output for volumecells
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::DumpGmshIntegrationCells( std::string name )
{
  NormalMesh().DumpGmshIntegrationCells( name );
}

/*--------------------------------------------------------------------------------------*
 * write gmsh output for volumecells
 *-------------------------------------------------------------------------------------*/
void GEO::CUT::MeshIntersection::DumpGmshVolumeCells( std::string name )
{
  NormalMesh().DumpGmshVolumeCells( name );
}

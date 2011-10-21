
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"

#include "cut_integrationcell.H"
#include "cut_volumecell.H"

#include "cut_meshintersection.H"

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

GEO::CUT::SideHandle * GEO::CUT::MeshIntersection::AddCutSide( int sid,
                                                               const std::vector<int> & nids,
                                                               DRT::Element::DiscretizationType distype,
                                                               int mi )
{
  // create side
  return cut_mesh_[mi]->CreateSide( sid, nids, distype );
}

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

void GEO::CUT::MeshIntersection::Cut( bool include_inner )
{
  Status();

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

  if ( options_.FindPositions() )
  {
    // find inside and outside positions of nodes
    m.FindNodePositions();

    // find number and connection of dofsets at nodes from cut volumes
//    m.FindNodalDOFSets( include_inner );
  }

  m.CreateIntegrationCells( 0, false );
  //m.RemoveEmptyVolumeCells();
  //m.MomentFitGaussWeights();

#ifdef DEBUGCUTLIBRARY
  //m.TestVolumeSurface();
  m.TestFacetArea();
#endif
  m.SimplifyIntegrationCells();

  Status();

#ifdef DEBUGCUTLIBRARY
  m.TestElementVolume( true );
#endif
}



void GEO::CUT::MeshIntersection::CreateNodalDofSetNEW( bool include_inner, DRT::Discretization & dis )
{

    std::set<int> eids; // eids of elements that are involved in CUT and include ele_vc_set_inside/ouside (no duplicates!)

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
                    sourrounding_elements.insert(pair<int, ElementHandle*>(adj_eid, e));
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

        std::vector<std::vector<int> > & nodaldofset_vc_sets_inside = eh->GetNodalDofSet_VcSets_Inside();
        std::vector<std::vector<int> > & nodaldofset_vc_sets_outside = eh->GetNodalDofSet_VcSets_Outside();

        if(include_inner )
        {
        	ConnectNodalDOFSets(nodes, include_inner, ele_vc_sets_inside, nodaldofset_vc_sets_inside);
        }

        ConnectNodalDOFSets(nodes, include_inner, ele_vc_sets_outside, nodaldofset_vc_sets_outside);

    }

#if(0)
    DumpGmshNumDOFSets(include_inner, eids, dis);
#endif
}


// set nodal cell set inside and outside
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

//	             cout << "\t number of connected_sets in element "<< i->first << " inside\t" << cell_sets_inside.size() << endl;
//	             cout << "\t number of connected_sets in element "<< i->first << " outside\t" << cell_sets_outside.size() << endl;


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



void GEO::CUT::MeshIntersection::ConnectNodalDOFSets( std::vector<Node *> & nodes, bool include_inner,
		const std::vector<plain_volumecell_set> & connected_vc_sets,
		std::vector<std::vector<int> > &    nodaldofset_vc_sets)
{

    for(std::vector<plain_volumecell_set>::const_iterator s=connected_vc_sets.begin();
        s!=connected_vc_sets.end();
        s++)
    {
        plain_volumecell_set cells = *s; // this is one connection of volumecells, connected via subelements, within one element

        std::vector<int > nds;

        // find this plain_volumecell_set in dof_cellsets_ vector of each node
        {
            for ( std::vector<Node*>::iterator i=nodes.begin();
                  i!=nodes.end();
                  ++i )
            {
                Node * n = *i;

                if( n->Id() >= 0) nds.push_back( n->DofSetNumberNEW( cells ) );
                else dserror("node with negative Id gets no dofnumber!");
            }

        }

        nodaldofset_vc_sets.push_back(nds);

    }
}


GEO::CUT::Node * GEO::CUT::MeshIntersection::GetNode( int nid ) const
{
  return mesh_.GetNode( nid );
}

GEO::CUT::ElementHandle * GEO::CUT::MeshIntersection::GetElement( int eid ) const
{
  return mesh_.GetElement( eid );
}

GEO::CUT::SideHandle * GEO::CUT::MeshIntersection::GetCutSide( int sid, int mi ) const
{
  return cut_mesh_[mi]->GetSide( sid );
}

void GEO::CUT::MeshIntersection::PrintCellStats()
{
  NormalMesh().PrintCellStats();
}

void GEO::CUT::MeshIntersection::Status()
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
  DumpGmshIntegrationCells( "integrationcells.pos" );
#endif
#endif
}


void GEO::CUT::MeshIntersection::DumpGmshNodalCellSet( std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets, DRT::Discretization & dis )
{
    std::string filename = "hallo"; //DRT::Problem::Instance()->OutputControlFile()->FileName();
    std::stringstream str;
    str << filename
        << ".NodalCellSet."
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


void GEO::CUT::MeshIntersection::DumpGmshCellSets( std::vector<plain_volumecell_set> & cell_sets,  DRT::Discretization & dis )
{
    std::string filename = "hallo"; //DRT::Problem::Instance()->OutputControlFile()->FileName();
    std::stringstream str;
    str << filename
        << ".CellSets."
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



//void GEO::CUT::MeshIntersection::DumpGmshNumDOFSets( bool include_inner, std::set<int> & eids, DRT::Discretization & dis )
void GEO::CUT::MeshIntersection::DumpGmshNumDOFSets(std::string filename, bool include_inner, DRT::Discretization & dis )
{
//    std::string filename = "hallo"; //DRT::Problem::Instance()->OutputControlFile()->FileName();
    std::stringstream str;
    str << filename
        << ".NumDOFSets."
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
    std::map<int, Node* > nodes;
    m.GetNodeMap( nodes );

    file << "View \"NumDofSets\" {\n";
    for ( std::map<int, Node* >::iterator i=nodes.begin(); i!=nodes.end(); ++i )
    {
        Node * n = i->second;
        Point * p = n->point();
        const double * x = p->X();

        // output just for real nodes of elements, not for shadow nodes
        if( n->Id() >= 0) file << "SP(" << x[0] << "," << x[1] << "," << x[2] << "){" << n->NumDofSets() << "};\n";
    }
    file << "};\n";
}


void GEO::CUT::MeshIntersection::DumpGmshVolumeCells( std::string name, bool include_inner )
{
  NormalMesh().DumpGmshVolumeCells( name, include_inner );
}

void GEO::CUT::MeshIntersection::DumpGmshIntegrationCells( std::string name )
{
  NormalMesh().DumpGmshIntegrationCells( name );
}

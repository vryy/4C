
#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyTraits.hpp>

#include "cut_tetmeshintersection.H"
#include "cut_side.H"
#include "cut_node.H"
#include "cut_facet.H"
#include "cut_volumecell.H"
#include "cut_integrationcell.H"
#include "cut_boundarycell.H"
#include "cut_pointpool.H"
#include "cut_options.H"

GEO::CUT::TetMeshIntersection::TetMeshIntersection( const Options & options,
                                                    Element * element,
                                                    const std::vector<std::vector<int> > & tets,
                                                    const std::vector<int> & accept_tets,
                                                    const std::vector<Point*> & points,
                                                    const std::set<Side*> & cut_sides,
                                                    bool levelset )
  : pp_( Teuchos::rcp( new PointPool ) ),
    mesh_( options, 1, pp_ ),
    cut_mesh_( options, 1, pp_, true )
{

  for ( std::vector<Point*>::const_iterator i=points.begin(); i!=points.end(); ++i )
  {
    Point * p = *i;
    Node * n = mesh_.GetNode( i - points.begin(), p->X() );
    Point * np = n->point();
    //np->Position( p->Position() );
    Register( p, np );
  }

  for ( std::vector<std::vector<int> >::const_iterator i=tets.begin(); i!=tets.end(); ++i )
  {
    const std::vector<int> & tet = *i;
    unsigned id = i-tets.begin();
    if ( accept_tets[id] )
      mesh_.GetElement( id, tet, *shards::getCellTopologyData< shards::Tetrahedron<4> >(), accept_tets[id] );
  }

  const std::set<Facet*> & element_facets = element->Facets();

  // Triangulated cut facets need to be converted to tri cut sides. This is
  // done after the normal cut sides.

  std::vector<Facet*> triangulated;
  std::map<Point*, Node*> nodemap;

  for ( std::set<Side*>::const_iterator i=cut_sides.begin(); i!=cut_sides.end(); ++i )
  {
    Side * s = *i;

    std::set<Facet*> facets;
    const std::vector<Facet*> & side_facets = s->Facets();
    for ( std::vector<Facet*>::const_iterator i=side_facets.begin(); i!=side_facets.end(); ++i )
    {
      Facet * f = *i;
      if ( element_facets.count( f ) > 0 )
      {
        facets.insert( f );
      }
    }

    for ( std::set<Facet*>::iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;

      if ( levelset or f->IsTriangulated() )
      {
        triangulated.push_back( f );
        std::set<Point*> points;
        f->AllPoints( points );
        for ( std::set<Point*>::iterator i=points.begin(); i!=points.end(); ++i )
        {
          Point *  p = *i;
          nodemap[ToChild( p )];
        }
      }
      else
      {
        CopyCutSide( s, f );
      }
    }
  }

  // create nodes

  cut_mesh_.NewNodesFromPoints( nodemap );

  // do triangulated facets (create extra cut sides)

  for ( std::vector<Facet*>::iterator i=triangulated.begin(); i!=triangulated.end(); ++i )
  {
    Facet * f = *i;
    Side * s = f->ParentSide();

    if ( f->IsTriangulated() )
    {
      const std::vector<std::vector<Point*> > & triangulation = f->Triangulation();
      for ( std::vector<std::vector<Point*> >::const_iterator i=triangulation.begin();
            i!=triangulation.end();
            ++i )
      {
        const std::vector<Point*> & tri = *i;
        if ( tri.size()!=3 )
          throw std::runtime_error( "tri3 expected" );
        std::vector<Node*> nodes;
        nodes.reserve( 3 );
        for ( std::vector<Point*>::const_iterator i=tri.begin(); i!=tri.end(); ++i )
        {
          Point * p = *i;
          nodes.push_back( nodemap[ToChild( p )] );
        }
        Side * cs = cut_mesh_.GetSide( s->Id(), nodes, shards::getCellTopologyData< shards::Triangle<3> >() );
        side_parent_to_child_[s].push_back( cs );
        for ( std::vector<Node*>::iterator i=nodes.begin(); i!=nodes.end(); ++i )
        {
          Node * n = *i;
          Point::InsertCut( NULL, cs, n );
        }
      }
    }
    else
    {
      const std::vector<Point*> & points = f->CornerPoints();

      std::vector<Node*> nodes;
      nodes.reserve( points.size() );
      for ( std::vector<Point*>::const_iterator i=points.begin(); i!=points.end(); ++i )
      {
        Point * p = *i;
        nodes.push_back( nodemap[ToChild( p )] );
      }

      switch ( points.size() )
      {
      case 2:
        // Degenerated nonsense. Why does that happen?

        // make sure the mapping entry exists
        side_parent_to_child_[s];
        break;
      case 3:
      {
        Side * cs = cut_mesh_.GetSide( s->Id(), nodes, shards::getCellTopologyData< shards::Triangle<3> >() );
        side_parent_to_child_[s].push_back( cs );
        for ( std::vector<Node*>::iterator i=nodes.begin(); i!=nodes.end(); ++i )
        {
          Node * n = *i;
          Point::InsertCut( NULL, cs, n );
        }
        break;
      }
      default:
        throw std::runtime_error( "facet with more that three points" );
      }
    }
  }

  Status();
}

void GEO::CUT::TetMeshIntersection::Cut( Mesh & parent_mesh, Element * element, const std::set<VolumeCell*> & parent_cells, int count, bool levelset )
{
  mesh_.Status();

#ifdef DEBUGCUTLIBRARY
  {
    std::stringstream str;
    str << "tetmesh" << count << ".pos";
    mesh_.DumpGmsh( str.str().c_str() );
  }
  {
    std::stringstream str;
    str << "tetcutmesh" << count << ".pos";
    cut_mesh_.DumpGmsh( str.str().c_str() );
  }
#endif

  std::set<Element*> elements_done;
  cut_mesh_.Cut( mesh_, elements_done );

  mesh_.Status();

  mesh_.MakeCutLines();
  mesh_.MakeFacets();
  mesh_.MakeVolumeCells();

  std::map<VolumeCell*, ChildCell> cellmap;

  MapVolumeCells( parent_mesh, element, parent_cells, cellmap );

#ifdef DEBUGCUTLIBRARY
  mesh_.DumpGmsh( "mesh.pos" );
#endif

#if 0
  // Not needed on this level.
  if ( mesh_.CreateOptions().FindPositions() )
  {
    mesh_.FindNodePositions();
  }
#endif

#ifdef DEBUGCUTLIBRARY
  mesh_.DumpGmsh( "mesh.pos" );
#endif

  mesh_.CreateIntegrationCells( count, levelset );

  Fill( parent_mesh, element, parent_cells, cellmap );
}

void GEO::CUT::TetMeshIntersection::MapVolumeCells( Mesh & parent_mesh, Element * element, const std::set<VolumeCell*> & parent_cells, std::map<VolumeCell*, ChildCell> & cellmap )
{
  std::set<VolumeCell*> done_child_cells;

  SeedCells( parent_mesh, parent_cells, cellmap, done_child_cells );

  int nonnodecells = 0;

  for ( std::set<VolumeCell*>::const_iterator i=parent_cells.begin();
        i!=parent_cells.end();
        ++i )
  {
    VolumeCell * vc = *i;
    ChildCell & cc = cellmap[vc];
    std::set<VolumeCell*> & childset = cc.cells_;

    if ( childset.size() > 0 )
    {
      Fill( vc, cc );
      std::copy( childset.begin(), childset.end(),
                 std::inserter( done_child_cells, done_child_cells.begin() ) );
    }
    else
    {
      nonnodecells += 1;
    }
  }

  // build surface cell map

  for ( std::map<VolumeCell*, ChildCell>::iterator i=cellmap.begin(); i!=cellmap.end(); ++i )
  {
    VolumeCell * vc = i->first;
    ChildCell & cc = i->second;
    //std::set<VolumeCell*> & childset = cc.cells_;

    // find parent facets on cut surface

    std::map<Side*, std::vector<Facet*> > & facetsonsurface = cc.facetsonsurface_;

    const std::set<Facet*> & parent_facets = vc->Facets();
    for ( std::set<Facet*>::const_iterator i=parent_facets.begin();
          i!=parent_facets.end();
          ++i )
    {
      Facet * f = *i;
      if ( f->OnCutSide() )
      {
        Side * s = f->ParentSide();
        std::map<Side*, std::vector<Side*> >::iterator j = side_parent_to_child_.find( s );
        if ( j==side_parent_to_child_.end() )
          throw std::runtime_error( "unknown parent cut facet" );
        std::vector<Side*> & side_vector = j->second;
        for ( std::vector<Side*>::iterator i=side_vector.begin(); i!=side_vector.end(); ++i )
        {
          Side * cs = *i;
          facetsonsurface[cs].push_back( f );
        }
      }
    }
  }

  // emergency seed cell filling

  while ( nonnodecells > 0 )
  {
    int backup = nonnodecells;

#if 1
    for ( std::map<VolumeCell*, ChildCell>::iterator i=cellmap.begin(); i!=cellmap.end(); ++i )
    {
      //VolumeCell * vc = i->first;
      ChildCell & cc = i->second;
      //std::set<VolumeCell*> & childset = cc.cells_;
      std::map<Side*, std::vector<Facet*> > & facetsonsurface = cc.facetsonsurface_;

      // match parent and child volumes at cut surface

      for ( std::map<Side*, std::vector<Facet*> >::iterator i=facetsonsurface.begin();
            i!=facetsonsurface.end();
            ++i )
      {
        Side * child_side = i->first;
        std::vector<Facet*> & parent_facets = i->second;

        // match only possible if the side has only one facet on this volume
        // there should be no other case?

        if ( parent_facets.size()==1 )
        {
          Facet * facet = parent_facets[0];
          const std::set<VolumeCell*> & parent_cells = facet->Cells();
          std::vector<ChildCell*> parent_cell_info;
          parent_cell_info.reserve( 2 );
          for ( std::set<VolumeCell*>::const_iterator i=parent_cells.begin();
                i!=parent_cells.end();
                ++i )
          {
            VolumeCell * vc = *i;

            std::map<VolumeCell*, ChildCell>::iterator j = cellmap.find( vc );
            if ( j!=cellmap.end() )
            {
              parent_cell_info.push_back( &j->second );
            }
          }

          // If there are less than two volumes at this cut, we have a touch
          // at a boundary.

          if ( parent_cell_info.size() == 1 )
          {
            if ( not parent_cell_info[0]->done_ )
            {
              const std::vector<Facet*> & child_facets = child_side->Facets();
              for ( std::vector<Facet*>::const_iterator i=child_facets.begin();
                    i!=child_facets.end();
                    ++i )
              {
                Facet * f = *i;
                const std::set<VolumeCell*> & child_cells = f->Cells();
                if ( child_cells.size()==1 )
                {
                  VolumeCell * c = *child_cells.begin();
                  if ( done_child_cells.count( c )==0 )
                  {
                    if ( not parent_cell_info[0]->ContainsChild( c ) )
                    {
                      parent_cell_info[0]->cells_.insert( c );
                      done_child_cells.insert( c );
                    }
                  }
                }
                else if ( child_cells.size()==2 )
                {
                  // odd.
                  throw std::runtime_error( "illegal number of neighbouring volume cells: child_cells.size()==2" );
                }
                else
                {
                  std::stringstream str;
                  str << "illegal number of neighbouring volume cells: child_cells.size()==" << child_cells.size();
                  throw std::runtime_error( str.str() );
                }
              }

              ChildCell & cc = *parent_cell_info[0];
              Fill( cc.parent_, cc );
              std::copy( cc.cells_.begin(), cc.cells_.end(),
                         std::inserter( done_child_cells, done_child_cells.begin() ) );
              nonnodecells -= 1;
            }
          }
          else if ( parent_cell_info.size() == 2 )
          {
            // Only useful if one volume is already done and the other one is not.

            int doneindex = -1;
            int otherindex = -1;
            if ( parent_cell_info[0]->done_ and not parent_cell_info[1]->done_ )
            {
              doneindex = 0;
              otherindex = 1;
            }
            else if ( parent_cell_info[1]->done_ and not parent_cell_info[0]->done_ )
            {
              doneindex = 1;
              otherindex = 0;
            }

            if ( doneindex > -1 )
            {
              // There are many child facets on the child side. But those facets
              // match the one parent facet. There must be no other facet
              // outside that region.

              bool found = false;

              const std::vector<Facet*> & child_facets = child_side->Facets();
//               if ( child_facets.size()==0 )
//               {
//                 throw std::runtime_error( "Wo sind die facets?" );
//               }

              for ( std::vector<Facet*>::const_iterator i=child_facets.begin();
                    i!=child_facets.end();
                    ++i )
              {
                Facet * f = *i;
                const std::set<VolumeCell*> & child_cells = f->Cells();
                if ( child_cells.size()==2 )
                {
                  std::vector<VolumeCell*> child_cell_vector;
                  child_cell_vector.reserve( 2 );
                  child_cell_vector.assign( child_cells.begin(), child_cells.end() );
                  if ( parent_cell_info[doneindex]->ContainsChild( child_cell_vector[0] ) )
                  {
                    VolumeCell * c = child_cell_vector[1];
                    if ( done_child_cells.count( c )==0 )
                    {
                      parent_cell_info[otherindex]->cells_.insert( c );
                      done_child_cells.insert( c );
                      found = true;
                    }
                  }
                  else if ( parent_cell_info[doneindex]->ContainsChild( child_cell_vector[1] ) )
                  {
                    VolumeCell * c = child_cell_vector[0];
                    if ( done_child_cells.count( c )==0 )
                    {
                      parent_cell_info[otherindex]->cells_.insert( c );
                      done_child_cells.insert( c );
                      found = true;
                    }
                  }
                  else
                  {
                    throw std::runtime_error( "child must be part of done parent cell" );
                  }
                }
                else if ( child_cells.size()==1 )
                {
                  VolumeCell * c = *child_cells.begin();
                  if ( not parent_cell_info[doneindex]->ContainsChild( c ) )
                  {
                    if ( done_child_cells.count( c )==0 )
                    {
                      parent_cell_info[otherindex]->cells_.insert( c );
                      done_child_cells.insert( c );
                      found = true;
                    }
                  }
                }
#if 0
                else if ( child_cells.size()==0 )
                {
                  // Ignore. There is a child facet that is not connected to
                  // anything. Oh, well.
                  f->Print();
                }
#endif
                else
                {
                  std::stringstream str;
                  str << "illegal number of neighbouring volume cells: child_cells.size() == " << child_cells.size();
                  throw std::runtime_error( str.str() );
                }
              }

              if ( found )
              {
                ChildCell & cc = *parent_cell_info[otherindex];
                Fill( cc.parent_, cc );
                std::copy( cc.cells_.begin(), cc.cells_.end(),
                           std::inserter( done_child_cells, done_child_cells.begin() ) );
                nonnodecells -= 1;
              }
              else
              {
                //std::cout << "not found after " << child_facets.size() << " facets from " << child_side << "\n";
              }
            }
            else
            {
              //std::cout << "parent_cell_info[0]->done_=" << parent_cell_info[0]->done_ << "  "
              //          << "parent_cell_info[1]->done_=" << parent_cell_info[1]->done_ << "\n";
            }
          }
          else
          {
            //std::cout << "parent_cell_info.size()=" << parent_cell_info.size() << "\n";
          }
        }
        else
        {
          //std::cout << "parent_facets.size()=" << parent_facets.size() << "\n";
        }
      }
    }
#endif

    if ( nonnodecells==1 )
    {
      VolumeCell * parent_vc = NULL;
      ChildCell * child_cells = NULL;
      for ( std::map<VolumeCell*, ChildCell>::iterator i=cellmap.begin(); i!=cellmap.end(); ++i )
      {
        VolumeCell * vc = i->first;
        ChildCell & cc = i->second;
        if ( not cc.done_ )
        {
          if ( parent_vc==NULL )
          {
            parent_vc = vc;
            child_cells = &cc;
          }
          else
          {
            throw std::runtime_error( "more than one open parent cells" );
          }
        }
      }
      if ( parent_vc==NULL )
      {
        throw std::runtime_error( "no open parent cell" );
      }

      ChildCell & cc = *child_cells;
      std::set<VolumeCell*> & childset = cc.cells_;

      std::set<VolumeCell*> done_child_cells;

      for ( std::map<VolumeCell*, ChildCell>::iterator i=cellmap.begin(); i!=cellmap.end(); ++i )
      {
        //VolumeCell * vc = i->first;
        ChildCell & cc = i->second;
        std::copy( cc.cells_.begin(), cc.cells_.end(),
                   std::inserter( done_child_cells, done_child_cells.begin() ) );
      }

      const std::list<Teuchos::RCP<VolumeCell> > & all_child_cells = mesh_.VolumeCells();
      for ( std::list<Teuchos::RCP<VolumeCell> >::const_iterator i=all_child_cells.begin();
            i!=all_child_cells.end();
            ++i )
      {
        VolumeCell * child_vc = &**i;
        if ( done_child_cells.count( child_vc )==0 )
        {
          childset.insert( child_vc );
        }
      }
      if ( childset.size() > 0 )
      {
        Fill( parent_vc, cc );
        nonnodecells -= 1;
      }
      else
      {
        //throw std::runtime_error( "no child cell for open parent cell" );

        // Empty parent cell. We did not get any children. The cell is most
        // probably too small.
        cc.done_ = true;
        nonnodecells -= 1;
      }
    }

    if ( nonnodecells > 0 )
    {
      // test if there are any volume cells left

      std::set<VolumeCell*> done_child_cells;
      for ( std::map<VolumeCell*, ChildCell>::iterator i=cellmap.begin(); i!=cellmap.end(); ++i )
      {
        ChildCell & cc = i->second;
        std::copy( cc.cells_.begin(), cc.cells_.end(),
                   std::inserter( done_child_cells, done_child_cells.begin() ) );
      }

      bool found = false;
      const std::list<Teuchos::RCP<VolumeCell> > & all_child_cells = mesh_.VolumeCells();
      for ( std::list<Teuchos::RCP<VolumeCell> >::const_iterator i=all_child_cells.begin();
            i!=all_child_cells.end();
            ++i )
      {
        VolumeCell * child_vc = &**i;
        if ( done_child_cells.count( child_vc )==0 )
        {
          found = true;
          break;
        }
      }

      if ( not found )
      {
        // Done. There are a few empty parent cells. We do not mind.
        nonnodecells = 0;
        for ( std::map<VolumeCell*, ChildCell>::iterator i=cellmap.begin(); i!=cellmap.end(); ++i )
        {
          ChildCell & cc = i->second;
          cc.done_ = true;
        }
      }
    }

    if ( backup == nonnodecells )
      throw std::runtime_error( "no progress in child cell--parent cell mapping" );
  }

  for ( std::map<VolumeCell*, ChildCell>::iterator i=cellmap.begin(); i!=cellmap.end(); ++i )
  {
    VolumeCell * vc = i->first;
    ChildCell & cc = i->second;
    std::set<VolumeCell*> & childset = cc.cells_;

    if ( not cc.done_ )
    {
      // finish partly filled volume cells
      Fill( vc, cc );
    }

    RegisterNewPoints( parent_mesh, childset );
  }

  // copy volume cell position

  if ( mesh_.CreateOptions().FindPositions() )
  {
    for ( std::map<VolumeCell*, ChildCell>::iterator i=cellmap.begin(); i!=cellmap.end(); ++i )
    {
      VolumeCell * vc = i->first;
      ChildCell & cc = i->second;
      std::set<VolumeCell*> & childset = cc.cells_;

      Point::PointPosition pos = vc->Position();
//       if ( pos==Point::undecided )
//         throw std::runtime_error( "undecided volume cell" );
      if ( pos!=Point::undecided )
      for ( std::set<VolumeCell*>::iterator i=childset.begin(); i!=childset.end(); ++i )
      {
        VolumeCell * c = *i;
        c->Position( pos );
      }
    }
  }
}

void GEO::CUT::TetMeshIntersection::SeedCells( Mesh & parent_mesh,
                                               const std::set<VolumeCell*> & parent_cells,
                                               std::map<VolumeCell*, ChildCell> & cellmap,
                                               std::set<VolumeCell*> & done_child_cells )
{
  std::map<Point*, std::vector<VolumeCell*> > parent_point_cells;

  for ( std::set<VolumeCell*>::const_iterator i=parent_cells.begin();
        i!=parent_cells.end();
        ++i )
  {
    VolumeCell * vc = *i;
    ChildCell & cc = cellmap[vc];
    cc.parent_ = vc;
    std::set<VolumeCell*> & childset = cc.cells_;

    std::set<Point*> volume_points;
    vc->GetAllPoints( parent_mesh, volume_points );

    // seed cells at parent element nodes (if unique)

    for ( std::set<Point*>::iterator i=volume_points.begin(); i!=volume_points.end(); ++i )
    {
      Point * p = *i;
      if ( p->Position()!=Point::oncutsurface )
      {
        Point * np = ToChild( p );
        FindVolumeCell( np, childset );
      }
    }

    for ( std::set<Point*>::iterator i=volume_points.begin(); i!=volume_points.end(); ++i )
    {
      Point * p = *i;
      parent_point_cells[p].push_back( vc );
    }
  }

  // seed cells with unique point

  for ( std::map<Point*, std::vector<VolumeCell*> >::iterator i=parent_point_cells.begin();
        i!=parent_point_cells.end();
        ++i )
  {
    Point * p = i->first;
    std::vector<VolumeCell*> & vcs = i->second;
    if ( vcs.size()==1 )
    {
      VolumeCell * vc = vcs[0];
      Point * np = ToChild( p );
      ChildCell & cc = cellmap[vc];
      std::set<VolumeCell*> & childset = cc.cells_;
      FindVolumeCell( np, childset );
    }
  }

  // collect done cells

  for ( std::map<VolumeCell*, ChildCell>::iterator i=cellmap.begin(); i!=cellmap.end(); ++i )
  {
    ChildCell & cc = i->second;
    std::set<VolumeCell*> & childset = cc.cells_;
    std::copy( childset.begin(), childset.end(),
               std::inserter( done_child_cells, done_child_cells.begin() ) );
  }

  // look at all points of each free child volume cell and see if there is a
  // unique parent volume cell to these points

  const std::list<Teuchos::RCP<VolumeCell> > & all_child_cells = mesh_.VolumeCells();
  for ( std::list<Teuchos::RCP<VolumeCell> >::const_iterator i=all_child_cells.begin();
        i!=all_child_cells.end();
        ++i )
  {
    VolumeCell * child_vc = &**i;
    if ( done_child_cells.count( child_vc )==0 )
    {
      std::set<Point*> child_cut_points;
      child_vc->GetAllPoints( mesh_, child_cut_points );

      // Remove all points that are new in the child mesh. Those do not at
      // all help to find the parent cell.
      for ( std::set<Point*>::iterator i=child_cut_points.begin(); i!=child_cut_points.end(); )
      {
        Point * p = *i;
        if ( child_to_parent_.count( p )==0 )
        {
          child_cut_points.erase( i++ );
          //child_cut_points.clear();
          //break;
        }
        else
        {
          ++i;
        }
      }

      if ( child_cut_points.size() > 0 )
      {
        std::set<VolumeCell*> used_parent_cells;

        std::set<Point*>::iterator j = child_cut_points.begin();
        Point * p = *j;
        FindVolumeCell( ToParent( p ), used_parent_cells );

        for ( ++j; j!=child_cut_points.end(); ++j )
        {
          Point * p = *j;
          std::set<VolumeCell*> upc;
          FindVolumeCell( ToParent( p ), upc );

          std::set<VolumeCell*> intersection;
          std::set_intersection( used_parent_cells.begin(), used_parent_cells.end(),
                                 upc.begin(), upc.end(),
                                 std::inserter( intersection, intersection.begin() ) );

          std::swap( used_parent_cells, intersection );

          if ( used_parent_cells.size()==0 )
            throw std::runtime_error( "no possible parent cell" );
        }

        if ( used_parent_cells.size()==1 )
        {
          VolumeCell * parent_vc = *used_parent_cells.begin();
          ChildCell & cc = cellmap[parent_vc];
          if ( cc.done_ )
          {
            //throw std::runtime_error( "free child cell to done parent cell?" );
          }
          else
          {
            std::set<VolumeCell*> & childset = cc.cells_;
            childset.insert( child_vc );
            done_child_cells.insert( child_vc );
          }
        }
      }
      else
      {
        //throw std::runtime_error( "child cell with all new points?" );
      }
    }
  }
}

void GEO::CUT::TetMeshIntersection::Fill( Mesh & parent_mesh, Element * element, const std::set<VolumeCell*> & parent_cells, std::map<VolumeCell*, ChildCell> & cellmap )
{
  for ( std::set<VolumeCell*>::const_iterator i=parent_cells.begin();
        i!=parent_cells.end();
        ++i )
  {
    VolumeCell * parent_cell = *i;
    std::set<VolumeCell*> & childset = cellmap[parent_cell].cells_;

    std::map<Side*, std::vector<Facet*> > & facetsonsurface = cellmap[parent_cell].facetsonsurface_;

    for ( std::set<VolumeCell*>::iterator i=childset.begin(); i!=childset.end(); ++i )
    {
      VolumeCell * vc = *i;
      const std::set<IntegrationCell*> & cells = vc->IntegrationCells();
      const std::set<BoundaryCell*> & bcells = vc->BoundaryCells();

      for ( std::set<IntegrationCell*>::const_iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        IntegrationCell * ic = *i;
        const std::vector<Point*> & points = ic->Points();

        std::vector<Point*> parent_points( points );
        ToParent( parent_mesh, parent_points );

        parent_cell->NewIntegrationCell( parent_mesh, ic->Shape(), parent_points );
      }
      for ( std::set<BoundaryCell*>::const_iterator i=bcells.begin(); i!=bcells.end(); ++i )
      {
        BoundaryCell * bc = *i;
        const std::vector<Point*> & points = bc->Points();

        std::vector<Point*> parent_points( points );
        ToParent( parent_mesh, parent_points );

        Facet * parent_facet = NULL;
        Facet * child_facet = bc->GetFacet();

        if ( not child_facet->OnCutSide() )
          throw std::runtime_error( "boundary cell not on cut surface" );

        std::map<Side*, std::vector<Facet*> >::iterator j = facetsonsurface.find( child_facet->ParentSide() );
        if ( j==facetsonsurface.end() )
        {
          std::stringstream str;
          str << ( *child_facet )
              << "\non boundary cell "
              << ( *child_facet->ParentSide() )
              << "\non unknown cut surface";
          throw std::runtime_error( str.str() );
        }

        std::vector<Facet*> & facets = j->second;
        if ( facets.size()==1 )
        {
          parent_facet = facets[0];
        }
        else if ( facets.size()>1 )
        {
          throw std::runtime_error( "Multiple parts of a cut side? Hard to believe." );
        }
        else
        {
          throw std::runtime_error( "empty list bug" );
        }

        parent_cell->NewBoundaryCell( parent_mesh, bc->Shape(), parent_facet, parent_points );
      }
    }
  }
}

void GEO::CUT::TetMeshIntersection::Fill( VolumeCell * parent_cell, ChildCell & childcell )
{
  std::set<VolumeCell*> & child_cells = childcell.cells_;

  if ( child_cells.size()==0 )
  {
    throw std::runtime_error( "failed to find seed cells" );
  }

  std::set<VolumeCell*> done_child_cells;

  while ( child_cells.size() > 0 )
  {
    for ( std::set<VolumeCell*>::iterator i=child_cells.begin(); i!=child_cells.end(); ++i )
    {
      VolumeCell * vc = *i;
      if ( done_child_cells.count( vc )==0 )
      {
        child_cells.erase( vc );
        done_child_cells.insert( vc );
        const std::set<Facet*> & facets = vc->Facets();
        for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
        {
          Facet * f = *i;
          if ( not f->OnCutSide() )
          {
            VolumeCell * nc = f->Neighbor( vc );
            if ( nc!=NULL and done_child_cells.count( nc )==0 )
            {
              child_cells.insert( nc );
            }
          }
        }
        break;
      }
    }
  }

  std::swap( child_cells, done_child_cells );
  childcell.done_ = true;
}

void GEO::CUT::TetMeshIntersection::RegisterNewPoints( Mesh & parent_mesh, const std::set<VolumeCell*> & childset )
{
  for ( std::set<VolumeCell*>::iterator i=childset.begin(); i!=childset.end(); ++i )
  {
    VolumeCell * vc = *i;
    const std::set<Facet*> & facets = vc->Facets();
    for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( f->OnCutSide() )
      {
        std::set<Point*> points;
        f->AllPoints( points );
        for ( std::set<Point*>::iterator i=points.begin(); i!=points.end(); ++i )
        {
          Point * p = *i;
          if ( child_to_parent_.count( p )==0 )
          {
            Point * pp = parent_mesh.NewPoint( p->X(), NULL, NULL );
            Register( pp, p );
          }
        }
      }
    }
  }
}

void GEO::CUT::TetMeshIntersection::Status()
{
#ifdef DEBUGCUTLIBRARY
  mesh_.DumpGmsh( "tetmesh.pos" );
  cut_mesh_.DumpGmsh( "tetcutmesh.pos" );
#endif
}

void GEO::CUT::TetMeshIntersection::FindVolumeCell( Point * p, std::set<VolumeCell*> & childset )
{
  const std::set<Facet*> & facets = p->Facets();
  for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;
    const std::set<VolumeCell*> & facet_cells = f->Cells();
    std::copy( facet_cells.begin(), facet_cells.end(), std::inserter( childset, childset.begin() ) );
  }

#ifdef DEBUGCUTLIBRARY
  if ( facets.size() > 0 )
  {
    for ( std::set<VolumeCell*>::const_iterator i=childset.begin(); i!=childset.end(); ++i )
    {
      VolumeCell * vc = *i;
      if ( vc->Contains( p ) )
        return;
    }
    throw std::runtime_error( "point not contained in volume cell" );
  }
#endif
}

void GEO::CUT::TetMeshIntersection::SwapPoints( Mesh & mesh, const std::map<Point*, Point*> & pointmap, std::vector<Point*> & points )
{
  std::vector<Point*> new_points;
  new_points.reserve( points.size() );
  for ( std::vector<Point*>::iterator i=points.begin(); i!=points.end(); ++i )
  {
    Point * p = *i;
    std::map<Point*, Point*>::const_iterator j = pointmap.find( p );
    if ( j==pointmap.end() )
    {
      Point * np = mesh.NewPoint( p->X(), NULL, NULL );
      new_points.push_back( np );
    }
    else
    {
      new_points.push_back( j->second );
    }
  }
  std::swap( new_points, points );
}

void GEO::CUT::TetMeshIntersection::SwapPoints( const std::map<Point*, Point*> & pointmap, std::vector<Point*> & points )
{
  std::vector<Point*> new_points;
  new_points.reserve( points.size() );
  for ( std::vector<Point*>::iterator i=points.begin(); i!=points.end(); ++i )
  {
    Point * p = *i;
    std::map<Point*, Point*>::const_iterator j = pointmap.find( p );
    if ( j==pointmap.end() )
    {
      throw std::runtime_error( "no such point" );
    }
    new_points.push_back( j->second );
  }
  std::swap( new_points, points );
}

void GEO::CUT::TetMeshIntersection::SwapPoints( const std::map<Point*, Point*> & pointmap, std::set<Point*> & points )
{
  std::set<Point*> new_points;
  for ( std::set<Point*>::iterator i=points.begin(); i!=points.end(); ++i )
  {
    Point * p = *i;
    std::map<Point*, Point*>::const_iterator j = pointmap.find( p );
    if ( j==pointmap.end() )
    {
      throw std::runtime_error( "no such point" );
    }
    new_points.insert( j->second );
  }
  std::swap( new_points, points );
}

GEO::CUT::Point * GEO::CUT::TetMeshIntersection::SwapPoint( const std::map<Point*, Point*> & pointmap, Point* point )
{
  std::map<Point*, Point*>::const_iterator j = pointmap.find( point );
  if ( j==pointmap.end() )
  {
    //throw std::runtime_error( "no such point" );
    return NULL;
  }
  return j->second;
}

void GEO::CUT::TetMeshIntersection::Register( Point * parent_point, Point * child_point )
{
  child_to_parent_[child_point ] = parent_point;
  parent_to_child_[parent_point] = child_point;
}

void GEO::CUT::TetMeshIntersection::CopyCutSide( Side * s, Facet * f )
{
  const std::vector<Node*> & nodes = s->Nodes();
  std::vector<int> nids;
  nids.reserve( nodes.size() );
  for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    Node * n = *i;
    nids.push_back( n->Id() );
    Point * p = n->point();

    Node * new_node = cut_mesh_.GetNode( n->Id(), p->X() );
    Point * np = ToChild( p );
    if ( np != NULL )
    {
      if ( new_node->point() != np )
      {
        throw std::runtime_error( "did not catch known cut point" );
      }
    }
    else
    {
      Register( p, new_node->point() );
    }
  }

  Side * cs = cut_mesh_.GetSide( s->Id(), nids, s->Topology() );

  side_parent_to_child_[s].push_back( cs );
  //side_child_to_parent_[cs] = s;

  // Copy cut point to cut surfaces, since a second cut search could result
  // in different cut points.

  const std::vector<Edge*> & old_edges =  s->Edges();
  const std::vector<Edge*> & new_edges = cs->Edges();

  for ( std::vector<Edge*>::const_iterator ei=old_edges.begin(); ei!=old_edges.end(); ++ei )
  {
    Edge * e = *ei;
    Edge * ne = new_edges[ei-old_edges.begin()];
    const std::vector<Point*> & cutpoints = e->CutPoints();
    for ( std::vector<Point*>::const_iterator i=cutpoints.begin();
          i!=cutpoints.end();
          ++i )
    {
      Point *  p = *i;
      Point * np = Point::NewPoint( mesh_, p->X(), p->t( e ), ne, NULL );
      np->Position( Point::oncutsurface );
    }
  }

  // Copy cut points from facets. If the facets is triangulated, there is a
  // middle point that needs to be introduces as a cut point.

  std::set<Point*> points;
  f->AllPoints( points );
  for ( std::set<Point*>::iterator i=points.begin(); i!=points.end(); ++i )
  {
    Point *  p = *i;
    Point * np = ToChild( p );
    if ( np!=NULL )
    {
      np->AddSide( cs );
      np->Position( Point::oncutsurface );
    }
  }

  // Copy cut lines to cut surfaces, since a second cut search could result
  // in different cut points.

  const std::vector<Line*> & cutlines = s->CutLines();

  for ( std::vector<Line*>::const_iterator i=cutlines.begin(); i!=cutlines.end(); ++i )
  {
    Line * l = *i;
    Point * p1 = ToChild( l->BeginPoint() );
    Point * p2 = ToChild( l->EndPoint() );

    if ( p1!=NULL and p2!=NULL )
    {
      std::vector<Line*> newlines;
      mesh_.NewLine( p1, p2, cs, NULL, NULL, &newlines );

      std::set<Edge*> edges;
      p1->CommonEdge( p2, edges );
      for ( std::set<Edge*>::iterator i=edges.begin(); i!=edges.end(); ++i )
      {
        Edge * e = *i;
        const std::set<Side*> & sides = e->Sides();
        for ( std::set<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
        {
          Side * s = *i;
          for ( std::vector<Line*>::iterator i=newlines.begin(); i!=newlines.end(); ++i )
          {
            Line * nl = *i;
            nl->AddSide( s );
          }
          const std::set<Element*> & elements = s->Elements();
          for ( std::set<Element*>::const_iterator i=elements.begin(); i!=elements.end(); ++i )
          {
            Element * e = *i;
            for ( std::vector<Line*>::iterator i=newlines.begin(); i!=newlines.end(); ++i )
            {
              Line * nl = *i;
              nl->AddElement( e );
            }
          }
        }
      }
    }
  }
}


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
    cut_mesh_( options, 1, pp_, true ),
    ls_side_( 1 )
{

//   std::vector<double> lsvs;
//   if ( levelset )
//   {
//     const std::vector<Node*> & nodes = element->Nodes();
//     lsvs.reserve( nodes.size() );
//     for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
//     {
//       Node * n = *i;
//       lsvs.push_back( n->LSV() );
//     }
//   }

  std::map<Point*, Node*> point_nodes;
  if ( levelset )
  {
    const std::vector<Node*> & nodes = element->Nodes();
    for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
    {
      Node * n = *i;
      point_nodes[n->point()] = n;
    }
  }

  for ( std::vector<Point*>::const_iterator i=points.begin(); i!=points.end(); ++i )
  {
    Point * p = *i;

//     if ( levelset )
//     {
//       LINALG::Matrix<3,1> xyz( p->X() );
//       LINALG::Matrix<3,1> rst;
//       element->LocalCoordinates( xyz, rst );
//       lsv = element->Scalar( lsvs, rst );
//     }

    double lsv = 0;
    if ( levelset )
    {
      std::map<Point*, Node*>::iterator j = point_nodes.find( p );
      if ( j!=point_nodes.end() )
      {
        Node * n = j->second;
        lsv = n->LSV();
      }
    }

    Node * n = mesh_.GetNode( i - points.begin(), p->X(), lsv );

    Point * np = n->point();
    np->Position( p->Position() );
    Register( p, np );
  }

  for ( std::vector<std::vector<int> >::const_iterator i=tets.begin(); i!=tets.end(); ++i )
  {
    const std::vector<int> & tet = *i;
    unsigned id = i-tets.begin();
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

    if ( levelset )
    {
      side_parent_to_child_[s].push_back( &ls_side_ );
      continue;
    }

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

    if ( facets.size() > 1 )
      throw std::runtime_error( "more than one facet between element and cut size?" );

    if ( facets.size() == 1 )
    {
      Facet * f = *facets.begin();

      if ( f->IsTriangulated() )
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
        const std::vector<Node*> & nodes = s->Nodes();
        std::vector<int> nids;
        nids.reserve( nodes.size() );
        for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
        {
          Node * n = *i;
          nids.push_back( n->Id() );
          Point * p = n->point();

//           if ( levelset )
//           {
//             LINALG::Matrix<3,1> xyz( p->X() );
//             LINALG::Matrix<3,1> rst;
//             element->LocalCoordinates( xyz, rst );
//             lsv = element->Scalar( lsvs, rst );
//           }

          Node * new_node = cut_mesh_.GetNode( n->Id(), p->X(), 0 );
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
            Line * nl = mesh_.NewLine( p1, p2, cs, NULL, NULL );

            std::set<Edge*> edges;
            p1->CommonEdge( p2, edges );
            for ( std::set<Edge*>::iterator i=edges.begin(); i!=edges.end(); ++i )
            {
              Edge * e = *i;
              const std::set<Side*> & sides = e->Sides();
              for ( std::set<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
              {
                Side * s = *i;
                nl->AddSide( s );
                const std::set<Element*> & elements = s->Elements();
                for ( std::set<Element*>::const_iterator i=elements.begin(); i!=elements.end(); ++i )
                {
                  Element * e = *i;
                  nl->AddElement( e );
                }
              }
            }
          }
        }
      }
    }
  }

  // create nodes

  cut_mesh_.NewNodesFromPoints( nodemap );

  if ( not levelset )
  {
    // do triangulated facets (create extra cut sides)

    for ( std::vector<Facet*>::iterator i=triangulated.begin(); i!=triangulated.end(); ++i )
    {
      Facet * f = *i;
      Side * s = f->ParentSide();

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
      }
    }
  }

  Status();
}

void GEO::CUT::TetMeshIntersection::Cut( Mesh & parent_mesh, Element * element, const std::set<VolumeCell*> & parent_cells, bool levelset )
{
  mesh_.Status();

  if ( levelset )
  {
    mesh_.Cut( ls_side_ );
  }
  else
  {
    std::set<Element*> elements_done;
    cut_mesh_.Cut( mesh_, elements_done );
  }

  mesh_.Status();

  mesh_.MakeFacets();
  mesh_.MakeVolumeCells();

  std::map<VolumeCell*, ChildCell> cellmap;

  MapVolumeCells( parent_mesh, element, parent_cells, cellmap );

  if ( mesh_.CreateOptions().FindPositions() )
  {
    if ( levelset )
    {
      mesh_.FindLSNodePositions();
    }
    else
    {
      mesh_.FindNodePositions();
    }
  }

#ifdef DEBUGCUTLIBRARY
  mesh_.DumpGmsh( "mesh.pos" );
#endif

  mesh_.CreateIntegrationCells( levelset );

  Fill( parent_mesh, element, parent_cells, cellmap, levelset );
}

void GEO::CUT::TetMeshIntersection::MapVolumeCells( Mesh & parent_mesh, Element * element, const std::set<VolumeCell*> & parent_cells, std::map<VolumeCell*, ChildCell> & cellmap )
{
  int nonnodecells = 0;

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

    for ( std::set<Point*>::iterator i=volume_points.begin(); i!=volume_points.end(); ++i )
    {
      Point * p = *i;
      if ( p->Position()!=Point::oncutsurface )
      {
        Point * np = ToChild( p );
        FindVolumeCell( np, childset );
      }
    }

    if ( childset.size() > 0 )
    {
      Fill( vc, cc );
    }
    else
    {
      nonnodecells += 1;
    }
  }

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

  while ( nonnodecells > 0 )
  {
    int backup = nonnodecells;
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
          // at a boundary. This is of little use here.

          if ( parent_cell_info.size() != 2 )
          {
            continue;
          }

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

            const std::vector<Facet*> & child_facets = child_side->Facets();
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
                  parent_cell_info[otherindex]->cells_.insert( child_cell_vector[1] );
                }
                else if ( parent_cell_info[doneindex]->ContainsChild( child_cell_vector[1] ) )
                {
                  parent_cell_info[otherindex]->cells_.insert( child_cell_vector[0] );
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
                  parent_cell_info[otherindex]->cells_.insert( c );
                }
              }
              else
              {
                throw std::runtime_error( "illegal number of neighbouring volume cells" );
              }
            }

            ChildCell & cc = *parent_cell_info[otherindex];
            Fill( cc.parent_, cc );
            nonnodecells -= 1;
          }
        }
      }
    }

    if ( backup == nonnodecells )
      throw std::runtime_error( "no progress" );
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
      if ( pos==Point::undecided )
        throw std::runtime_error( "undecided volume cell" );
      for ( std::set<VolumeCell*>::iterator i=childset.begin(); i!=childset.end(); ++i )
      {
        VolumeCell * c = *i;
        c->Position( pos );
      }
    }
  }
}

void GEO::CUT::TetMeshIntersection::Fill( Mesh & parent_mesh, Element * element, const std::set<VolumeCell*> & parent_cells, std::map<VolumeCell*, ChildCell> & cellmap, bool levelset )
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
        ToParent( parent_points );

        // debug
        ic->Volume();

        parent_cell->NewIntegrationCell( parent_mesh, ic->Shape(), parent_points );
      }
      for ( std::set<BoundaryCell*>::const_iterator i=bcells.begin(); i!=bcells.end(); ++i )
      {
        BoundaryCell * bc = *i;
        const std::vector<Point*> & points = bc->Points();

        std::vector<Point*> parent_points( points );
        ToParent( parent_points );

        Facet * parent_facet = NULL;
        Facet * child_facet = bc->GetFacet();

        if ( not child_facet->OnCutSide() )
          throw std::runtime_error( "boundary cell not on cut surface" );

        std::map<Side*, std::vector<Facet*> >::iterator j = facetsonsurface.find( child_facet->ParentSide() );
        if ( j==facetsonsurface.end() )
          throw std::runtime_error( "boundary cell on unknown cut surface" );

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
//   std::set<VolumeCell*> cells;
  const std::set<Facet*> & facets = p->Facets();
  for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;
    std::set<Point*> facet_points;
    const std::set<VolumeCell*> & facet_cells = f->Cells();
    std::copy( facet_cells.begin(), facet_cells.end(), std::inserter( childset, childset.begin() ) );
  }
//   if ( cells.size()!=1 )
//     throw std::runtime_error( "expect one volume at nodal point" );
//   return *cells.begin();
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


//#include "../drt_geometry/intersection_templates.H"

#include "cut_intersection.H"
#include "cut_position.H"
#include "cut_facet.H"
#include "cut_volumecell.H"

#include <string>
#include <stack>

#include "cut_element.H"

bool GEO::CUT::Element::Cut( Mesh & mesh, Side & side )
{
  bool cut = false;

  // find nodal points inside the element
  const std::vector<Node*> side_nodes = side.Nodes();
  for ( std::vector<Node*>::const_iterator i=side_nodes.begin(); i!=side_nodes.end(); ++i )
  {
    Node * n = *i;
    Point * p = n->point();

    if ( not p->IsCut( this ) )
    {
      if ( PointInside( p ) )
      {
        p->AddElement( this );
        cut = true;
      }
    }
    else
    {
      cut = true;
    }
  }

  const std::vector<Side*> & sides = Sides();

  for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
  {
    Side * s = *i;
    FindCutPoints( mesh, *s, side );
  }

  for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
  {
    Side * s = *i;
    if ( FindCutLines( mesh, *s, side ) )
    {
      cut = true;
    }
  }

  // find lines inside the element
  const std::vector<Edge*> & side_edges = side.Edges();
  for ( std::vector<Edge*>::const_iterator i=side_edges.begin(); i!=side_edges.end(); ++i )
  {
    Edge * e = *i;
    std::vector<Point*> line;
    e->CutPointsInside( this, line );
    std::vector<Point*>::iterator i = line.begin();
    if ( i!=line.end() )
    {
      Point * bp = *i;
      for ( ++i; i!=line.end(); ++i )
      {
        Point * ep = *i;
        mesh.NewLine( bp, ep, &side, NULL, this );
        bp = ep;
      }
    }
  }

  if ( cut )
  {
    side.CreateLineSegment( mesh, this );
    cut_faces_.insert( &side );
    return true;
  }
  else
  {
    return false;
  }
}

bool GEO::CUT::Element::FindCutPoints( Mesh & mesh, Side & side, Side & other )
{
  bool cut = side.FindCutPoints( mesh, this, other );
  bool reverse_cut = other.FindCutPoints( mesh, this, side );
  return cut or reverse_cut;
}

bool GEO::CUT::Element::FindCutLines( Mesh & mesh, Side & side, Side & other )
{
  bool cut = side.FindCutLines( mesh, this, other );
  bool reverse_cut = other.FindCutLines( mesh, this, side );
  return cut or reverse_cut;
}

void GEO::CUT::Element::MakeFacets( Mesh & mesh )
{
  if ( facets_.size()==0 )
  {
    const std::vector<Side*> & sides = Sides();
    for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
    {
      Side & side = **i;
      SideElementCutFilter filter( &side, this );
      side.MakeOwnedSideFacets( mesh, filter, facets_ );
    }
    for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
    {
      Side & side = **i;
      side.MakeSideCutFacets( mesh, this, facets_ );
    }
    for ( std::set<Side*>::iterator i=cut_faces_.begin(); i!=cut_faces_.end(); ++i )
    {
      Side & cut_side = **i;
      cut_side.MakeInternalFacets( mesh, this, facets_ );
    }
  }
}

void GEO::CUT::Element::FindNodePositions()
{
  LINALG::Matrix<3,1> xyz;
  LINALG::Matrix<3,1> rst;

  const std::vector<Node*> & nodes = Nodes();
  for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    Node * n = *i;
    Point * p = n->point();
    Point::PointPosition pos = p->Position();
    if ( pos==Point::undecided )
    {
      bool done = false;
      const std::set<Facet*> & facets = p->Facets();
      for ( std::set<Facet*>::iterator i=facets.begin(); i!=facets.end(); ++i )
      {
        Facet * f = *i;
        for ( std::set<Side*>::iterator i=cut_faces_.begin(); i!=cut_faces_.end(); ++i )
        {
          Side * s = *i;

          // Only take a side that belongs to one of this points facets and
          // shares a cut edge with this point. If there are multiple cut
          // sides within the element (facets), only the close one will always
          // give the right direction.
          if ( f->IsCutSide( s ) and p->CommonCutEdge( s )!=NULL )
          {
            if ( p->IsCut( s ) )
            {
              p->Position( Point::oncutsurface );
            }
            else
            {
              p->Coordinates( xyz.A() );
              s->LocalCoordinates( xyz, rst );
              double d = rst( 2, 0 );
              if ( fabs( d ) > MINIMALTOL )
              {
                if ( d > 0 )
                {
                  p->Position( Point::outside );
                }
                else
                {
                  p->Position( Point::inside );
                }
              }
              else
              {
                // within the cut plane but not cut by the side
                break;
              }
            }
            done = true;
            break;
          }
        }
        if ( done )
          break;
      }
      if ( p->Position()==Point::undecided )
      {
        // Still undecided! No facets with cut side attached! Will be set in a
        // minute.
      }
    }
    else if ( pos==Point::outside or pos==Point::inside )
    {
      // The nodal position is already known. Set it to my facets. If the
      // facets are already set, this will not have much effect anyway. But on
      // multiple cuts we avoid unset facets this way.
      const std::set<Facet*> & facets = p->Facets();
      for ( std::set<Facet*>::iterator i=facets.begin(); i!=facets.end(); ++i )
      {
        Facet * f = *i;
        f->Position( pos );
      }
    }
  }
}

bool GEO::CUT::Element::IsCut()
{
  if ( cut_faces_.size()>0 )
  {
    return true;
  }
  for ( std::vector<Side*>::const_iterator i=Sides().begin(); i!=Sides().end(); ++i )
  {
    Side & side = **i;
    if ( side.IsCut() )
    {
      return true;
    }
  }
  return false;
}

bool GEO::CUT::Element::OnSide( Facet * f )
{
  if ( not f->HasHoles() )
  {
    return OnSide( f->Points() );
  }
  return false;
}

bool GEO::CUT::Element::OnSide( const std::vector<Point*> & facet_points )
{
  const std::vector<Node*> & nodes = Nodes();
  for ( std::vector<Point*>::const_iterator i=facet_points.begin();
        i!=facet_points.end();
        ++i )
  {
    Point * p = *i;
    if ( not p->NodalPoint( nodes ) )
    {
      return false;
    }
  }

  std::set<Point*, PointPidLess> points;
  std::copy( facet_points.begin(), facet_points.end(),
             std::inserter( points, points.begin() ) );

  for ( std::vector<Side*>::const_iterator i=Sides().begin(); i!=Sides().end(); ++i )
  {
    Side & side = **i;
    if ( side.OnSide( points ) )
    {
      return true;
    }
  }

  return false;
}


void GEO::CUT::Element::GetIntegrationCells( std::set<GEO::CUT::IntegrationCell*> & cells )
{
  for ( std::set<VolumeCell*>::iterator i=cells_.begin(); i!=cells_.end(); ++i )
  {
    VolumeCell * vc = *i;
    vc->GetIntegrationCells( cells );
  }
}

void GEO::CUT::Element::GetBoundaryCells( std::set<GEO::CUT::BoundaryCell*> & bcells )
{
  for ( std::set<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    if ( cut_faces_.count( f->ParentSide() )!= 0 )
    {
      f->GetBoundaryCells( bcells );
    }
  }
}

void GEO::CUT::Element::GetCutPoints( std::set<Point*> & cut_points )
{
  for ( std::vector<Side*>::const_iterator i=Sides().begin(); i!=Sides().end(); ++i )
  {
    Side * side = *i;

    for ( std::set<Side*>::iterator i=cut_faces_.begin(); i!=cut_faces_.end(); ++i )
    {
      Side * other = *i;
      side->GetCutPoints( this, *other, cut_points );
    }
  }
}

void GEO::CUT::Element::MakeVolumeCells( Mesh & mesh )
{
  std::map<std::pair<Point*, Point*>, std::set<Facet*> > lines;
  for ( std::set<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    f->GetLines( lines );
  }

  // Alle Facets einsammeln, die sich zu zweit eine Linie teilen. Jede Seite
  // nur ein Facet. Das sollte im Element eindeutig sein. Damit sollten die
  // Volumina erstellt werden können. Die Facets der Löcher sind mitzunehmen.

  std::set<Facet*> facets_done;

  std::set<Facet*> all_facets = facets_;
  for ( std::set<Facet*>::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    const std::set<Facet*> holes = f->Holes();
    std::copy( holes.begin(), holes.end(), std::inserter( all_facets, all_facets.begin() ) );
  }

  // We really need to do the side facets first, since there is now way to
  // decide which side facet to chose, when all we have is an cut facet. Thus
  // we need a special order here.
  //
  // This might still fail in bizarre cases. The right solution needs to
  // consider concave connections and so on.

  std::vector<Facet*> all_facets_sorted;
  all_facets_sorted.reserve( all_facets.size() );

  for ( std::vector<Node*>::iterator i=nodes_.begin(); i!=nodes_.end(); ++i )
  {
    Node * n = *i;
    Point * p = n->point();
    const std::set<Facet*> & point_facets = p->Facets();
    for ( std::set<Facet*>::const_iterator i=point_facets.begin();
          i!=point_facets.end();
          ++i )
    {
      Facet * f = *i;
      if ( all_facets.count( f ) > 0 )
      {
        all_facets.erase( f );
        all_facets_sorted.push_back( f );
      }
    }
  }

  std::copy( all_facets.begin(), all_facets.end(), std::back_inserter( all_facets_sorted ) );
  all_facets.clear();

#if 0
  int facet_count = 0;
#endif

  for ( std::vector<Facet*>::iterator i=all_facets_sorted.begin(); i!=all_facets_sorted.end(); ++i )
  {
    Facet * f = *i;
    if ( facets_done.count( f )==0 and ( OwnedSide( f->ParentSide() ) or OnSide( f ) ) )
    {
      std::set<std::vector<Facet*>::iterator> new_facets;
      std::set<Facet*> collected_facets;
      //std::set<Side*> sides_done;
      std::map<std::pair<Point*, Point*>, std::set<Facet*> > done_lines;

      new_facets.insert( std::find( all_facets_sorted.begin(), all_facets_sorted.end(), f ) );
      f->GetLines( done_lines );

      {
#if 0
        std::stringstream str;
        str << "facet" << ( facet_count++ ) << ".plot";
        std::ofstream file( str.str().c_str() );
        f->Print( file );
        file.close();
#endif
      }

      while ( new_facets.size() > 0 )
      {
        std::set<std::vector<Facet*>::iterator>::iterator first = new_facets.begin();
        Facet * f = **first;
        new_facets.erase( first );

        collected_facets.insert( f );
        Side * facet_side = f->ParentSide();

        std::map<std::pair<Point*, Point*>, std::set<Facet*> > facet_lines;
        f->GetLines( facet_lines );

        for ( std::map<std::pair<Point*, Point*>, std::set<Facet*> >::iterator i=facet_lines.begin();
              i!=facet_lines.end();
              ++i )
        {
          const std::pair<Point*, Point*> & line = i->first;
          if ( done_lines[line].size() < 2 )
          {
            std::set<Facet*> & facets = lines[line];

            Facet * found_facet = NULL;
            for ( std::set<Facet*>::iterator i=facets.begin(); i!=facets.end(); ++i )
            {
              Facet * f = *i;
              if ( collected_facets.count( f )==0 )
              {
                bool found = false;
                if ( not OwnedSide( f->ParentSide() ) ) // OnSide( f )
                {
                  found = true;
                }
                else if ( facet_side != f->ParentSide() and
                          facets_done.count( f )==0 )
                {
                  found = true;
                }

                if ( found )
                {
                  if ( found_facet==NULL )
                  {
                    found_facet = f;
                  }
                  else
                  {
                    // undecided. Ignore all matches. Hope we are able to close
                    // the volume anyway. We should be.
                    found_facet = NULL;
                    break;
                  }
                }
              }
            }
            if ( found_facet!=NULL )
            {
              new_facets.insert( std::find( all_facets_sorted.begin(), all_facets_sorted.end(), found_facet ) );
              found_facet->GetLines( done_lines );

              {
#if 0
                std::stringstream str;
                str << "facet" << ( facet_count++ ) << ".plot";
                std::ofstream file( str.str().c_str() );
                found_facet->Print( file );
                file.close();
#endif
              }

            }
          }
        }
      }

      // test for open lines in collected_facets

      std::map<std::pair<Point*, Point*>, std::set<Facet*> > volume_lines;
      for ( std::set<Facet*>::iterator i=collected_facets.begin();
            i!=collected_facets.end();
            ++i )
      {
        Facet * f = *i;
        f->GetLines( volume_lines );
      }

      for ( std::map<std::pair<Point*, Point*>, std::set<Facet*> >::iterator i=volume_lines.begin();
            i!=volume_lines.end();
            ++i )
      {
        std::set<Facet*> & facets = i->second;
        if ( facets.size()!=2 )
        {
#if 1
          std::ofstream file( "collected_facets.plot" );
          for ( std::set<Facet*>::iterator i=collected_facets.begin();
                i!=collected_facets.end();
                ++i )
          {
            Facet * f = *i;
            f->Print( file );
          }
          file.close();
#endif
          throw std::runtime_error( "not properly closed line in volume cell" );
        }
      }

      // Create new cell and remember done stuff!

      std::copy( collected_facets.begin(),
                 collected_facets.end(),
                 std::inserter( facets_done, facets_done.begin() ) );

      {
#if 0
        std::stringstream str;
        str << "volume" << cells_.size() << ".plot";
        std::ofstream file( str.str().c_str() );
        for ( std::set<Facet*>::iterator i=collected_facets.begin();
              i!=collected_facets.end();
              ++i )
        {
          Facet * f = *i;
          f->Print( file );
        }
        file.close();
#endif
      }

      cells_.insert( mesh.NewVolumeCell( collected_facets, volume_lines, this ) );

#if 0
      // create extra volumes for any holes in the wall

      for ( std::set<Facet*>::iterator i=collected_facets.begin(); i!=collected_facets.end(); ++i )
      {
        Facet * f = *i;
        if ( OwnedSide( f->ParentSide() ) )
        {
          if ( f->HasHoles() )
          {
          }
        }
      }
#endif
    }
  }

  if ( facets_done.size() != all_facets_sorted.size() )
    throw std::runtime_error( "unhandled facets" );
}

bool GEO::CUT::ConcreteElement<DRT::Element::tet4>::PointInside( Point* p )
{
  Position<DRT::Element::tet4> pos( *this, *p );
  return pos.Compute();
}

bool GEO::CUT::ConcreteElement<DRT::Element::hex8>::PointInside( Point* p )
{
  Position<DRT::Element::hex8> pos( *this, *p );
  return pos.Compute();
}

bool GEO::CUT::ConcreteElement<DRT::Element::wedge6>::PointInside( Point* p )
{
  Position<DRT::Element::wedge6> pos( *this, *p );
  return pos.Compute();
}

bool GEO::CUT::ConcreteElement<DRT::Element::pyramid5>::PointInside( Point* p )
{
  Position<DRT::Element::pyramid5> pos( *this, *p );
  return pos.Compute();
}


void GEO::CUT::ConcreteElement<DRT::Element::tet4>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position<DRT::Element::tet4> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::ConcreteElement<DRT::Element::hex8>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position<DRT::Element::hex8> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::ConcreteElement<DRT::Element::wedge6>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position<DRT::Element::wedge6> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
}

void GEO::CUT::ConcreteElement<DRT::Element::pyramid5>::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position<DRT::Element::pyramid5> pos( *this, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
//     throw std::runtime_error( "global point not within element" );
  }
  rst = pos.LocalCoordinates();
}

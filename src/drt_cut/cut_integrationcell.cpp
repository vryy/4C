
#include "cut_position2d.H"
#include "cut_integrationcell.H"
#include "cut_facet.H"
#include "cut_mesh.H"
#include "cut_boundarycell.H"
#include "cut_element.H"
#include "cut_volumecell.H"
#include "cut_projection.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#include "../drt_geometry/element_volume.H"


GEO::CUT::Point::PointPosition GEO::CUT::IntegrationCell::VolumePosition( const std::set<Facet*> & facets )
{
  Point::PointPosition pos = Point::undecided;
  for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;
    pos = std::min( pos, static_cast<Point::PointPosition>( f->PositionSideId() ) );
  }
  return pos;
}

int GEO::CUT::Hex8IntegrationCell::hex8totet4[5][4] = {
  {0, 1, 3, 4},
  {1, 2, 3, 6},
  {4, 5, 1, 6},
  {6, 7, 3, 4},
  {1, 6, 3, 4}
};

int GEO::CUT::Wedge6IntegrationCell::wedge6totet4[3][4] = {
  {0, 1, 2, 3},
  {3, 4, 1, 5},
  {1, 5, 2, 3}
};


int GEO::CUT::Pyramid5IntegrationCell::pyramid5totet4[2][4] = {
  {0, 1, 3, 4},
  {1, 2, 3, 4}
};


bool GEO::CUT::IntegrationCell::CreateCells( Mesh & mesh,
                                             VolumeCell * cell,
                                             Point::PointPosition position,
                                             const std::set<Facet*> & facets,
                                             std::set<IntegrationCell*> & integrationcells )
{
  for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;
    if ( f->HasHoles() )
    {
      return false;
    }
  }

  Element * parent = cell->ParentElement();

  switch( parent->Shape() )
  {
  case DRT::Element::hex8:
  {
    std::map<Side*, int> touched;
    for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( not f->OnCutSide() )
      {
        touched[ f->ParentSide() ] += 1;
      }
    }

    if ( touched.size() != 5 )
    {
      return false;
    }

    for ( std::map<Side*, int>::iterator i=touched.begin(); i!=touched.end(); ++i )
    {
      if ( i->second != 1 )
      {
        return false;
      }
    }

    int count = 0;
    std::set<int> uncut;
    const std::vector<Side*> & sides = parent->Sides();
    for ( std::vector<Side*>::const_iterator i=sides.begin(); i!=sides.end(); ++i )
    {
      Side * s = *i;
      if ( not s->IsCut() )
      {
        uncut.insert( count );
      }
      count += 1;
    }

    if ( uncut.size() == 2 )
    {
      double r = 0.;
      int axis = -1;

      // We use shards face numbering. Be careful.
      if ( uncut.count( 4 ) > 0 and uncut.count( 5 ) > 0 )
      {
        axis = 2;
        if ( touched.count( sides[5] ) > 0 )
        {
          r = 1;
        }
        else
        {
          r = -1;
        }
      }
      else if ( uncut.count( 0 ) > 0 and uncut.count( 2 ) > 0 )
      {
        axis = 1;
        if ( touched.count( sides[2] ) > 0 )
        {
          r = 1;
        }
        else
        {
          r = -1;
        }
      }
      else if ( uncut.count( 1 ) > 0 and uncut.count( 3 ) > 0 )
      {
        axis = 0;
        if ( touched.count( sides[1] ) > 0 )
        {
          r = 1;
        }
        else
        {
          r = -1;
        }
      }
      else
      {
        return false;
      }

      Hex8Projection projection( mesh, parent, cell, position, facets, integrationcells, axis, r );
      return true;
    }

    return false;
  }
  default:
    break;
  }
  return false;
}

bool GEO::CUT::Hex8IntegrationCell::CreateCell( Mesh & mesh,
                                                VolumeCell * cell,
                                                const std::set<Facet*> & facets,
                                                std::set<IntegrationCell*> & integrationcells )
{
  if ( facets.size()==6 )
  {
    for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( not f->Equals( DRT::Element::quad4 ) )
      {
        return false;
      }
    }

    // Must not be concave?!

    // So we have a hex8 here. Create it.
    //
    // Find two sides that do not share a common point, find positive
    // directions for both, and find the edges between those sides.

    std::set<Facet*>::const_iterator i=facets.begin();
    Facet * bot = *i;
    Facet * top = NULL;
    for ( ++i; i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( not f->Touches( bot ) )
      {
        top = f;
        break;
      }
    }
    if ( top==NULL )
    {
      throw std::runtime_error( "illegal hex8 cell" );
    }

    const std::vector<Point*> & bot_points = bot->CornerPoints();
    const std::vector<Point*> & top_points = top->CornerPoints();

    std::vector<Point*> points( 8, static_cast<Point*>( NULL ) );

    std::copy( bot_points.begin(), bot_points.end(), &points[0] );

    std::vector<Point*>::iterator end = points.begin()+4;

    for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( f!=bot and f!=top )
      {
        const std::vector<Point*> & side_points = f->CornerPoints();
        for ( std::vector<Point*>::const_iterator i=side_points.begin();
              i!=side_points.end();
              ++i )
        {
          Point * p = *i;
          std::vector<Point*>::iterator pointpos1 = std::find( points.begin(), end, p );
          if ( pointpos1!=end )
          {
            ++i;
            if ( i==side_points.end() )
            {
              throw std::runtime_error( "illegal hex8 cell" );
            }

            std::vector<Point*>::iterator pointpos2 = std::find( points.begin(), end, *i );
            if ( pointpos2==end )
            {
              i = side_points.end();
              pointpos2 = std::find( points.begin(), end, side_points[3] );
              if ( pointpos2==end )
              {
                throw std::runtime_error( "illegal hex8 cell" );
              }
              std::swap( pointpos1, pointpos2 );
            }

            unsigned pos = i - side_points.begin() - 1;
            Point * top_point2 = side_points[( pos + 2 ) % 4];
            Point * top_point1 = side_points[( pos + 3 ) % 4];

            if ( std::find( top_points.begin(), top_points.end(), top_point1 )==top_points.end() or
                 std::find( top_points.begin(), top_points.end(), top_point2 )==top_points.end() )
            {
              throw std::runtime_error( "illegal hex8 cell" );
            }

            pos = ( pointpos1 - points.begin() ) + 4;
            if ( points[pos]==NULL )
            {
              points[pos] = top_point1;
            }
            else if ( points[pos] != top_point1 )
            {
              throw std::runtime_error( "illegal hex8 cell" );
            }

            pos = ( pointpos2 - points.begin() ) + 4;
            if ( points[pos]==NULL )
            {
              points[pos] = top_point2;
            }
            else if ( points[pos] != top_point2 )
            {
              throw std::runtime_error( "illegal hex8 cell" );
            }

            break;
          }
        }
      }
    }

    // Find the geometric distance between bottom side and top side. Since we
    // can have arbitrary concave situations, all points have to be
    // checked. If we cannot decide on an orientation, we reject the cell.

    LINALG::Matrix<3,4> bot_xyze;
    LINALG::Matrix<3,4> top_xyze;

    bot->CornerCoordinates( bot_xyze.A() );
    top->CornerCoordinates( top_xyze.A() );

    int distance_counter = 0;
    for ( int i=0; i<4; ++i )
    {
      Position2d<DRT::Element::quad4> bot_distance( bot_xyze, LINALG::Matrix<3,1>( top_xyze.A()+3*i, true ) );

      bot_distance.Compute();

      if ( bot_distance.Distance() > 0 )
      {
        distance_counter += 1;
      }
      if ( bot_distance.Distance() < 0 )
      {
        distance_counter -= 1;
      }
    }

    if ( distance_counter==-4 )
    {
      std::vector<Point*> rpoints( 8 );
      rpoints[0] = points[3];
      rpoints[1] = points[2];
      rpoints[2] = points[1];
      rpoints[3] = points[0];
      rpoints[4] = points[7];
      rpoints[5] = points[6];
      rpoints[6] = points[5];
      rpoints[7] = points[4];

      for ( int i=0; i<6; ++i )
      {
        std::vector<Point*> side( 4 );
        for ( int j=0; j<4; ++j )
        {
          side[j] = rpoints[DRT::UTILS::eleNodeNumbering_hex27_surfaces[i][j]];
        }
        Facet * f = FindFacet( facets, side );
        if ( f->OnCutSide() )
          Quad4BoundaryCell::CreateCell( mesh, cell, f, side );
      }

      Point::PointPosition position = VolumePosition( facets );
      CreateCell( mesh, cell, position, rpoints, integrationcells );
      return true;
    }
    else if ( distance_counter==4 )
    {
      for ( int i=0; i<6; ++i )
      {
        std::vector<Point*> side( 4 );
        for ( int j=0; j<4; ++j )
        {
          side[j] = points[DRT::UTILS::eleNodeNumbering_hex27_surfaces[i][j]];
        }
        Facet * f = FindFacet( facets, side );
        if ( f->OnCutSide() )
          Quad4BoundaryCell::CreateCell( mesh, cell, f, side );
      }

      Point::PointPosition position = VolumePosition( facets );
      CreateCell( mesh, cell, position, points, integrationcells );
      return true;
    }
    else
    {
      return false;
    }
  }
  return false;
}

void GEO::CUT::Hex8IntegrationCell::CreateCell( Mesh & mesh,
                                                VolumeCell * cell,
                                                Point::PointPosition position,
                                                const std::vector<Point*> & points,
                                                std::set<IntegrationCell*> & integrationcells )
{
#if 0
  integrationcells.insert( mesh.NewHex8Cell( position, points, cell ) );
#else
  std::vector<Point*> tet4_points( 4 );
  for ( int i=0; i<5; ++i )
  {
    SetTetPoints( hex8totet4[i], points, tet4_points );
    integrationcells.insert( mesh.NewTet4Cell( position, tet4_points, cell ) );
  }
#endif
}

bool GEO::CUT::Tet4IntegrationCell::CreateCell( Mesh & mesh,
                                                VolumeCell * cell,
                                                const std::set<Facet*> & facets,
                                                std::set<IntegrationCell*> & integrationcells )
{
  if ( facets.size()==4 )
  {
    for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( not f->Equals( DRT::Element::tri3 ) )
      {
        return false;
      }
    }

    // So we have a tet4 here. Create it.

    std::set<Facet*>::const_iterator i=facets.begin();
    Facet * bot = *i;

    const std::vector<Point*> & bot_points = bot->CornerPoints();
    Point* top_point = NULL;

    for ( ++i; i!=facets.end(); ++i )
    {
      Facet * f = *i;

      const std::vector<Point*> & side_points = f->CornerPoints();
      for ( std::vector<Point*>::const_iterator i=side_points.begin();
            i!=side_points.end();
            ++i )
      {
        Point * p = *i;
        if ( std::find( bot_points.begin(), bot_points.end(), p )==bot_points.end() )
        {
          if ( top_point==NULL )
          {
            top_point = p;
          }
          else if ( top_point!=p )
          {
            throw std::runtime_error( "illegal tet4 cell" );
          }
        }
      }
    }

//     for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
//     {
//       Facet * f = *i;
//       f->NewTri3Cell( mesh );
//     }

    LINALG::Matrix<3,3> bot_xyze;
    LINALG::Matrix<3,1> top_xyze;

    bot->CornerCoordinates( bot_xyze.A() );
    top_point->Coordinates( top_xyze.A() );

    Position2d<DRT::Element::tri3> bot_distance( bot_xyze, top_xyze );
    bot_distance.Compute();

    if ( bot_distance.Distance() > 0 )
    {
      std::vector<Point*> points;
      points.reserve( 4 );
      std::copy( bot_points.begin(), bot_points.end(), std::back_inserter( points ) );
      points.push_back( top_point );

      for ( int i=0; i<4; ++i )
      {
        std::vector<Point*> side( 3 );
        for ( int j=0; j<3; ++j )
        {
          side[j] = points[DRT::UTILS::eleNodeNumbering_tet10_surfaces[i][j]];
        }
        Facet * f = FindFacet( facets, side );
        if ( f->OnCutSide() )
          Tri3BoundaryCell::CreateCell( mesh, cell, f, side );
      }

      Point::PointPosition position = VolumePosition( facets );
      CreateCell( mesh, cell, position, points, integrationcells );
      return true;
    }
    else if ( bot_distance.Distance() < 0 )
    {
      std::vector<Point*> points;
      points.reserve( 4 );
      std::copy( bot_points.rbegin(), bot_points.rend(), std::back_inserter( points ) );
      points.push_back( top_point );

      for ( int i=0; i<4; ++i )
      {
        std::vector<Point*> side( 3 );
        for ( int j=0; j<3; ++j )
        {
          side[j] = points[DRT::UTILS::eleNodeNumbering_tet10_surfaces[i][j]];
        }
        Facet * f = FindFacet( facets, side );
        if ( f->OnCutSide() )
          Tri3BoundaryCell::CreateCell( mesh, cell, f, side );
      }

      Point::PointPosition position = VolumePosition( facets );
      CreateCell( mesh, cell, position, points, integrationcells );
      return true;
    }
    else
    {
      throw std::runtime_error( "planar cell" );
    }
  }
  return false;
}

GEO::CUT::Tet4IntegrationCell * GEO::CUT::Tet4IntegrationCell::CreateCell( Mesh & mesh,
                                                                           VolumeCell * cell,
                                                                           Point::PointPosition position,
                                                                           const std::vector<Point*> & tet )
{
  return mesh.NewTet4Cell( position, tet, cell );
}

void GEO::CUT::Tet4IntegrationCell::CreateCell( Mesh & mesh, VolumeCell * cell, Point::PointPosition position, const std::vector<Point*> & points, std::set<IntegrationCell*> & integrationcells )
{
  integrationcells.insert( mesh.NewTet4Cell( position, points, cell ) );
}

bool GEO::CUT::Wedge6IntegrationCell::CreateCell( Mesh & mesh,
                                                  VolumeCell * cell,
                                                  const std::set<Facet*> & facets,
                                                  std::set<IntegrationCell*> & integrationcells )
{
  if ( facets.size()==5 )
  {
    std::vector<Facet*> tris;
    std::vector<Facet*> quads;

    for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( f->Equals( DRT::Element::tri3 ) )
      {
        tris.push_back( f );
      }
      else if ( f->Equals( DRT::Element::quad4 ) )
      {
        quads.push_back( f );
      }
      else
      {
        return false;
      }
    }

    if ( tris.size()!=2 or quads.size()!=3 )
    {
      return false;
    }

    Facet * bot = tris[0];
    Facet * top = tris[1];

    const std::vector<Point*> & bot_points = bot->CornerPoints();
    const std::vector<Point*> & top_points = top->CornerPoints();

    std::vector<Point*> points( 6, static_cast<Point*>( NULL ) );

    std::copy( bot_points.begin(), bot_points.end(), &points[0] );

    std::vector<Point*>::iterator end = points.begin()+3;

    for ( std::vector<Facet*>::const_iterator i=quads.begin(); i!=quads.end(); ++i )
    {
      Facet * f = *i;
      const std::vector<Point*> & side_points = f->CornerPoints();
      for ( std::vector<Point*>::const_iterator i=side_points.begin();
            i!=side_points.end();
            ++i )
      {
        Point * p = *i;
        std::vector<Point*>::iterator pointpos1 = std::find( points.begin(), end, p );
        if ( pointpos1!=end )
        {
          ++i;
          if ( i==side_points.end() )
          {
            throw std::runtime_error( "illegal wedge6 cell" );
          }

          std::vector<Point*>::iterator pointpos2 = std::find( points.begin(), end, *i );
          if ( pointpos2==end )
          {
            i = side_points.end();
            pointpos2 = std::find( points.begin(), end, side_points[3] );
            if ( pointpos2==end )
            {
              throw std::runtime_error( "illegal wedge6 cell" );
            }
            std::swap( pointpos1, pointpos2 );
          }

          unsigned pos = i - side_points.begin() - 1;
          Point * top_point2 = side_points[( pos + 2 ) % 4];
          Point * top_point1 = side_points[( pos + 3 ) % 4];

          if ( std::find( top_points.begin(), top_points.end(), top_point1 )==top_points.end() or
               std::find( top_points.begin(), top_points.end(), top_point2 )==top_points.end() )
          {
            throw std::runtime_error( "illegal wedge6 cell" );
          }

          pos = ( pointpos1 - points.begin() ) + 3;
          if ( points[pos]==NULL )
          {
            points[pos] = top_point1;
          }
          else if ( points[pos] != top_point1 )
          {
            throw std::runtime_error( "illegal wedge6 cell" );
          }

          pos = ( pointpos2 - points.begin() ) + 3;
          if ( points[pos]==NULL )
          {
            points[pos] = top_point2;
          }
          else if ( points[pos] != top_point2 )
          {
            throw std::runtime_error( "illegal wedge6 cell" );
          }

          break;
        }
      }
    }

    // Find the geometric distance between bottom side and top side. Since we
    // can have arbitrary concave situations, all points have to be
    // checked. If we cannot decide on an orientation, we reject the cell.

    LINALG::Matrix<3,3> bot_xyze;
    LINALG::Matrix<3,3> top_xyze;

    bot->CornerCoordinates( bot_xyze.A() );
    top->CornerCoordinates( top_xyze.A() );

    int distance_counter = 0;
    for ( int i=0; i<3; ++i )
    {
      Position2d<DRT::Element::tri3> bot_distance( bot_xyze, LINALG::Matrix<3,1>( top_xyze.A()+3*i, true ) );

      bot_distance.Compute();

      if ( bot_distance.Distance() > 0 )
      {
        distance_counter += 1;
      }
      if ( bot_distance.Distance() < 0 )
      {
        distance_counter -= 1;
      }
    }

//     if ( distance_counter==-3 or
//          distance_counter==3 )
//     {
//       for ( std::vector<Facet*>::const_iterator i=tris.begin(); i!=tris.end(); ++i )
//       {
//         Facet * f = *i;
//         f->NewTri3Cell( mesh );
//       }
//       for ( std::vector<Facet*>::const_iterator i=quads.begin(); i!=quads.end(); ++i )
//       {
//         Facet * f = *i;
//         f->NewQuad4Cell( mesh );
//       }
//     }

    if ( distance_counter==-3 )
    {
      std::vector<Point*> rpoints( 6 );
      rpoints[0] = points[2];
      rpoints[1] = points[1];
      rpoints[2] = points[0];
      rpoints[3] = points[5];
      rpoints[4] = points[4];
      rpoints[5] = points[3];

      for ( int i=0; i<2; ++i )
      {
        std::vector<Point*> side( 3 );
        for ( int j=0; j<3; ++j )
        {
          side[j] = rpoints[DRT::UTILS::eleNodeNumbering_wedge15_trisurfaces[i][j]];
        }
        Facet * f = FindFacet( facets, side );
        if ( f->OnCutSide() )
          Tri3BoundaryCell::CreateCell( mesh, cell, f, side );
      }
      for ( int i=0; i<3; ++i )
      {
        std::vector<Point*> side( 4 );
        for ( int j=0; j<4; ++j )
        {
          side[j] = rpoints[DRT::UTILS::eleNodeNumbering_wedge15_quadsurfaces[i][j]];
        }
        Facet * f = FindFacet( facets, side );
        if ( f->OnCutSide() )
          Quad4BoundaryCell::CreateCell( mesh, cell, f, side );
      }

      Point::PointPosition position = VolumePosition( facets );
      CreateCell( mesh, cell, position, rpoints, integrationcells );
      return true;
    }
    else if ( distance_counter==3 )
    {
      for ( int i=0; i<2; ++i )
      {
        std::vector<Point*> side( 3 );
        for ( int j=0; j<3; ++j )
        {
          side[j] = points[DRT::UTILS::eleNodeNumbering_wedge15_trisurfaces[i][j]];
        }
        Facet * f = FindFacet( facets, side );
        if ( f->OnCutSide() )
          Tri3BoundaryCell::CreateCell( mesh, cell, f, side );
      }
      for ( int i=0; i<3; ++i )
      {
        std::vector<Point*> side( 4 );
        for ( int j=0; j<4; ++j )
        {
          side[j] = points[DRT::UTILS::eleNodeNumbering_wedge15_quadsurfaces[i][j]];
        }
        Facet * f = FindFacet( facets, side );
        if ( f->OnCutSide() )
          Quad4BoundaryCell::CreateCell( mesh, cell, f, side );
      }

      Point::PointPosition position = VolumePosition( facets );
      CreateCell( mesh, cell, position, points, integrationcells );
      return true;
    }
    else
    {
      return false;
    }
  }
  return false;
}

void GEO::CUT::Wedge6IntegrationCell::CreateCell( Mesh & mesh, VolumeCell * cell, Point::PointPosition position, const std::vector<Point*> & points, std::set<IntegrationCell*> & integrationcells )
{
#if 0
  integrationcells.insert( mesh.NewWedge6Cell( position, points, cell ) );
#else
  std::vector<Point*> tet4_points( 4 );
  for ( int i=0; i<3; ++i )
  {
    SetTetPoints( wedge6totet4[i], points, tet4_points );
    integrationcells.insert( mesh.NewTet4Cell( position, tet4_points, cell ) );
  }
#endif
}

bool GEO::CUT::Pyramid5IntegrationCell::CreateCell( Mesh & mesh,
                                                    VolumeCell * cell,
                                                    const std::set<Facet*> & facets,
                                                    std::set<IntegrationCell*> & integrationcells )
{
  if ( facets.size()==5 )
  {
    std::vector<Facet*> tris;
    std::vector<Facet*> quads;

    for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
    {
      Facet * f = *i;
      if ( f->Equals( DRT::Element::tri3 ) )
      {
        tris.push_back( f );
      }
      else if ( f->Equals( DRT::Element::quad4 ) )
      {
        quads.push_back( f );
      }
      else
      {
        return false;
      }
    }

    if ( tris.size()!=4 or quads.size()!=1 )
    {
      return false;
    }

    Facet * bot = quads[0];

    const std::vector<Point*> & bot_points = bot->CornerPoints();
    Point* top_point = NULL;

    for ( std::vector<Facet*>::iterator i=tris.begin(); i!=tris.end(); ++i )
    {
      Facet * f = *i;

      const std::vector<Point*> & side_points = f->CornerPoints();
      for ( std::vector<Point*>::const_iterator i=side_points.begin();
            i!=side_points.end();
            ++i )
      {
        Point * p = *i;
        if ( std::find( bot_points.begin(), bot_points.end(), p )==bot_points.end() )
        {
          if ( top_point==NULL )
          {
            top_point = p;
          }
          else if ( top_point!=p )
          {
            // Corner point confusion. This is actually not a pyramid5.
            //
            //throw std::runtime_error( "illegal pyramid5 cell" );
            return false;
          }
        }
      }
    }

//     for ( std::vector<Facet*>::const_iterator i=tris.begin(); i!=tris.end(); ++i )
//     {
//       Facet * f = *i;
//       f->NewTri3Cell( mesh );
//     }
//     for ( std::vector<Facet*>::const_iterator i=quads.begin(); i!=quads.end(); ++i )
//     {
//       Facet * f = *i;
//       f->NewQuad4Cell( mesh );
//     }

    LINALG::Matrix<3,3> bot_xyze;
    LINALG::Matrix<3,1> top_xyze;

    bot->CornerCoordinates( bot_xyze.A() );
    top_point->Coordinates( top_xyze.A() );

    Position2d<DRT::Element::tri3> bot_distance( bot_xyze, top_xyze );
    bot_distance.Compute();

    if ( bot_distance.Distance() > 0 )
    {
      std::vector<Point*> points;
      points.reserve( 5 );
      std::copy( bot_points.begin(), bot_points.end(), std::back_inserter( points ) );
      points.push_back( top_point );

      for ( int i=0; i<4; ++i )
      {
        std::vector<Point*> side( 3 );
        for ( int j=0; j<3; ++j )
        {
          side[j] = points[DRT::UTILS::eleNodeNumbering_pyramid5_trisurfaces[i][j]];
        }
        Facet * f = FindFacet( facets, side );
        if ( f->OnCutSide() )
          Tri3BoundaryCell::CreateCell( mesh, cell, f, side );
      }
      for ( int i=0; i<1; ++i )
      {
        std::vector<Point*> side( 4 );
        for ( int j=0; j<4; ++j )
        {
          side[j] = points[DRT::UTILS::eleNodeNumbering_pyramid5_quadsurfaces[i][j]];
        }
        Facet * f = FindFacet( facets, side );
        if ( f->OnCutSide() )
          Quad4BoundaryCell::CreateCell( mesh, cell, f, side );
      }

      Point::PointPosition position = VolumePosition( facets );
      CreateCell( mesh, cell, position, points, integrationcells );
      return true;
    }
    else if ( bot_distance.Distance() < 0 )
    {
      std::vector<Point*> points;
      points.reserve( 5 );
      std::copy( bot_points.rbegin(), bot_points.rend(), std::back_inserter( points ) );
      points.push_back( top_point );

      for ( int i=0; i<4; ++i )
      {
        std::vector<Point*> side( 3 );
        for ( int j=0; j<3; ++j )
        {
          side[j] = points[DRT::UTILS::eleNodeNumbering_pyramid5_trisurfaces[i][j]];
        }
        Facet * f = FindFacet( facets, side );
        if ( f->OnCutSide() )
          Tri3BoundaryCell::CreateCell( mesh, cell, f, side );
      }
      for ( int i=0; i<1; ++i )
      {
        std::vector<Point*> side( 4 );
        for ( int j=0; j<4; ++j )
        {
          side[j] = points[DRT::UTILS::eleNodeNumbering_pyramid5_quadsurfaces[i][j]];
        }
        Facet * f = FindFacet( facets, side );
        if ( f->OnCutSide() )
          Quad4BoundaryCell::CreateCell( mesh, cell, f, side );
      }

      Point::PointPosition position = VolumePosition( facets );
      CreateCell( mesh, cell, position, points, integrationcells );
      return true;
    }
    else
    {
      throw std::runtime_error( "planar cell" );
    }
  }
  return false;
}

void GEO::CUT::Pyramid5IntegrationCell::CreateCell( Mesh & mesh, VolumeCell * cell, Point::PointPosition position, const std::vector<Point*> & points, std::set<IntegrationCell*> & integrationcells )
{
#if 0
  integrationcells.insert( mesh.NewPyramid5Cell( position, points, cell ) );
#else
  std::vector<Point*> tet4_points( 4 );
  for ( int i=0; i<2; ++i )
  {
    SetTetPoints( pyramid5totet4[i], points, tet4_points );
    integrationcells.insert( mesh.NewTet4Cell( position, tet4_points, cell ) );
  }
#endif
}

void GEO::CUT::Hex8IntegrationCell::DumpGmsh( std::ofstream & file )
{
  file << "SH(";
  for ( int i=0; i<8; ++i )
  {
    if ( i > 0 )
      file << ",";
    file << xyz_( 0, i ) << ","
         << xyz_( 1, i ) << ","
         << xyz_( 2, i );
  }
  file << "){";
  for ( int i=0; i<8; ++i )
  {
    if ( i > 0 )
      file << ",";
    file << position_;
  }
  file << "};\n";
}

void GEO::CUT::Tet4IntegrationCell::DumpGmsh( std::ofstream & file )
{
  file << "SS(";
  for ( int i=0; i<4; ++i )
  {
    if ( i > 0 )
      file << ",";
    file << xyz_( 0, i ) << ","
         << xyz_( 1, i ) << ","
         << xyz_( 2, i );
  }
  file << "){";
  for ( int i=0; i<4; ++i )
  {
    if ( i > 0 )
      file << ",";
    file << position_;
  }
  file << "};\n";
}

void GEO::CUT::Wedge6IntegrationCell::DumpGmsh( std::ofstream & file )
{
  file << "SI(";
  for ( int i=0; i<6; ++i )
  {
    if ( i > 0 )
      file << ",";
    file << xyz_( 0, i ) << ","
         << xyz_( 1, i ) << ","
         << xyz_( 2, i );
  }
  file << "){";
  for ( int i=0; i<6; ++i )
  {
    if ( i > 0 )
      file << ",";
    file << position_;
  }
  file << "};\n";
}

void GEO::CUT::Pyramid5IntegrationCell::DumpGmsh( std::ofstream & file )
{
  file << "SP(";
  for ( int i=0; i<5; ++i )
  {
    if ( i > 0 )
      file << ",";
    file << xyz_( 0, i ) << ","
         << xyz_( 1, i ) << ","
         << xyz_( 2, i );
  }
  file << "){";
  for ( int i=0; i<5; ++i )
  {
    if ( i > 0 )
      file << ",";
    file << position_;
  }
  file << "};\n";
}

double GEO::CUT::IntegrationCell::Volume() const
{
  return GEO::ElementVolume( Shape(), xyz_ );
}

int GEO::CUT::Hex8IntegrationCell::CubatureDegree( DRT::Element::DiscretizationType elementshape ) const
{
  switch ( elementshape )
  {
  case DRT::Element::hex8:
    return 3;
  case DRT::Element::hex20:
    return 4;
  case DRT::Element::hex27:
    return 4;
  default:
    throw std::runtime_error( "no rule defined for this element type" );
  }
}

int GEO::CUT::Tet4IntegrationCell::CubatureDegree( DRT::Element::DiscretizationType elementshape ) const
{
  switch ( elementshape )
  {
  case DRT::Element::hex8:
    return 11;
  case DRT::Element::hex20:
    return 12;
  case DRT::Element::hex27:
    return 14;
  case DRT::Element::tet4:
    return 5;
  case DRT::Element::tet10:
    return 7;
  default:
    throw std::runtime_error( "no rule defined for this element type" );
  }
}

int GEO::CUT::Wedge6IntegrationCell::CubatureDegree( DRT::Element::DiscretizationType elementshape ) const
{
  return 4;
}

int GEO::CUT::Pyramid5IntegrationCell::CubatureDegree( DRT::Element::DiscretizationType elementshape ) const
{
  return 4;
}

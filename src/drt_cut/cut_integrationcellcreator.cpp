
#include "cut_integrationcellcreator.H"
#include "cut_position2d.H"
#include "cut_integrationcell.H"
#include "cut_facet.H"
#include "cut_mesh.H"
#include "cut_boundarycell.H"
#include "cut_element.H"


bool GEO::CUT::IntegrationCellCreator::CreateCells( Mesh & mesh, const std::set<VolumeCell*> & cells )
{
  IntegrationCellCreator creator;

  for ( std::set<VolumeCell*>::const_iterator i=cells.begin(); i!=cells.end(); ++i )
  {
    VolumeCell * cell = *i;
    bool found = ( creator.CreateTet4Cell( mesh, cell, cell->Facets() ) or
                   creator.CreateHex8Cell( mesh, cell, cell->Facets() ) or
                   creator.CreateWedge6Cell( mesh, cell, cell->Facets() ) );

    // pyramids are not save right now.
    // Pyramid5IntegrationCell::CreateCell( mesh, this, facets_, creator ) );
    //creator.CreatePyramid5Cell( mesh, cell, cell->Facets() )

    if ( not found )
    {
      return false;
    }
  }
  creator.Execute( mesh );
  return true;
}

bool GEO::CUT::IntegrationCellCreator::CreateCell( Mesh & mesh,
                                                   DRT::Element::DiscretizationType shape,
                                                   VolumeCell * cell )
{
  IntegrationCellCreator creator;

  bool success;
  switch ( shape )
  {
  case DRT::Element::tet4:
    success = creator.CreateTet4Cell( mesh, cell, cell->Facets() );
    break;
  case DRT::Element::hex8:
    success = creator.CreateHex8Cell( mesh, cell, cell->Facets() );
    break;
  case DRT::Element::wedge6:
    success = creator.CreateWedge6Cell( mesh, cell, cell->Facets() );
    break;
  case DRT::Element::pyramid5:
    //success = creator.CreatePyramid5Cell( mesh, cell, cell->Facets() );
    success = false;
    break;
  default:
    throw std::runtime_error( "unsupported element shape" );
  }
  if ( success )
  {
    creator.Execute( mesh );
  }
  return success;
}

void GEO::CUT::IntegrationCellCreator::Execute( Mesh & mesh )
{
  for ( std::map<VolumeCell*, ic>::iterator i=cells_.begin(); i!=cells_.end(); ++i )
  {
    VolumeCell * vc = i->first;
    ic & cell = i->second;
    cell.Execute( mesh, vc );
    for ( std::vector<bc>::iterator i=cell.boundary_.begin(); i!=cell.boundary_.end(); ++i )
    {
      bc & bcell = *i;
      bcell.Execute( mesh, vc );
    }
  }
}

bool GEO::CUT::IntegrationCellCreator::CreateTet4Cell( Mesh & mesh, VolumeCell * cell, const std::set<Facet*> & facets )
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
          //Tri3BoundaryCell::CreateCell( mesh, cell, f, side );
          AddSide( cell, f, DRT::Element::tri3, side );
      }

      //cell->NewTet4Cell( mesh, points );
      Add( cell, DRT::Element::tet4, points );
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
          //Tri3BoundaryCell::CreateCell( mesh, cell, f, side );
          AddSide( cell, f, DRT::Element::tri3, side );
      }

      //cell->NewTet4Cell( mesh, points );
      Add( cell, DRT::Element::tet4, points );
      return true;
    }
    else
    {
      throw std::runtime_error( "planar cell" );
    }
  }
  return false;
}

bool GEO::CUT::IntegrationCellCreator::CreateHex8Cell( Mesh & mesh, VolumeCell * cell, const std::set<Facet*> & facets )
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
          //Quad4BoundaryCell::CreateCell( mesh, cell, f, side );
          AddSide( cell, f, DRT::Element::quad4, side );
      }

      //cell->NewHex8Cell( mesh, rpoints );
      Add( cell, DRT::Element::hex8, rpoints );
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
          //Quad4BoundaryCell::CreateCell( mesh, cell, f, side );
          AddSide( cell, f, DRT::Element::quad4, side );
      }

      //cell->NewHex8Cell( mesh, points );
      Add( cell, DRT::Element::hex8, points );
      return true;
    }
    else
    {
      return false;
    }
  }
  return false;
}

bool GEO::CUT::IntegrationCellCreator::CreateWedge6Cell( Mesh & mesh, VolumeCell * cell, const std::set<Facet*> & facets )
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
          //Tri3BoundaryCell::CreateCell( mesh, cell, f, side );
          AddSide( cell, f, DRT::Element::tri3, side );
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
          //Quad4BoundaryCell::CreateCell( mesh, cell, f, side );
          AddSide( cell, f, DRT::Element::quad4, side );
      }

      //cell->NewWedge6Cell( mesh, rpoints );
      Add( cell, DRT::Element::wedge6, rpoints );
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
          //Tri3BoundaryCell::CreateCell( mesh, cell, f, side );
          AddSide( cell, f, DRT::Element::tri3, side );
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
          //Quad4BoundaryCell::CreateCell( mesh, cell, f, side );
          AddSide( cell, f, DRT::Element::quad4, side );
      }

      //cell->NewWedge6Cell( mesh, points );
      Add( cell, DRT::Element::wedge6, points );
      return true;
    }
    else
    {
      return false;
    }
  }
  return false;
}

bool GEO::CUT::IntegrationCellCreator::CreatePyramid5Cell( Mesh & mesh, VolumeCell * cell, const std::set<Facet*> & facets )
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

    LINALG::Matrix<3,4> bot_xyze;
    LINALG::Matrix<3,1> top_xyze;

    bot->CornerCoordinates( bot_xyze.A() );
    top_point->Coordinates( top_xyze.A() );

    Position2d<DRT::Element::quad4> bot_distance( bot_xyze, top_xyze );
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
          //Tri3BoundaryCell::CreateCell( mesh, cell, f, side );
          AddSide( cell, f, DRT::Element::tri3, side );
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
          //Quad4BoundaryCell::CreateCell( mesh, cell, f, side );
          AddSide( cell, f, DRT::Element::quad4, side );
      }

      //cell->NewPyramid5Cell( mesh, points );
      Add( cell, DRT::Element::pyramid5, points );
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
          //Tri3BoundaryCell::CreateCell( mesh, cell, f, side );
          AddSide( cell, f, DRT::Element::tri3, side );
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
          //Quad4BoundaryCell::CreateCell( mesh, cell, f, side );
          AddSide( cell, f, DRT::Element::quad4, side );
      }

      //cell->NewPyramid5Cell( mesh, points );
      Add( cell, DRT::Element::pyramid5, points );
      return true;
    }
    else
    {
      throw std::runtime_error( "planar cell" );
    }
  }
  return false;
}


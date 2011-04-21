
#include "cut_projection.H"
#include "cut_volumecell.H"
#include "cut_facet.H"
#include "cut_element.H"
#include "cut_mesh.H"
#include "cut_integrationcell.H"
#include "cut_boundarycell.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"

void GEO::CUT::Hex8Projection::HorizontalCut( Mesh & mesh,
                                              Element * element,
                                              VolumeCell * cell,
                                              Point::PointPosition position,
                                              const std::set<Facet*> & facets,
                                              std::set<IntegrationCell*> & integrationcells,
                                              int axis,
                                              double r )
{
  std::set<Point*> cut_points;
  cell->GetAllPoints( mesh, cut_points );

  std::vector<Point*> points;
  points.reserve( cut_points.size() );
  points.assign( cut_points.begin(), cut_points.end() );

  //std::cout << "hex8 projection along axis " << axis << " to " << r << "\n";

  // find all inner points that need projecting

  std::set<Point*> inner;

  std::vector<Facet*> inner_facets;

  for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;
    if ( f->OnCutSide() )
    {
      inner_facets.push_back( f );
      const std::vector<Point*> & points = f->CornerPoints();
      std::copy( points.begin(), points.end(), std::inserter( inner, inner.begin() ) );
    }
  }

  std::vector<Point*> inner_points;
  std::vector<Point*> projected_points;

  inner_points.reserve( inner.size() );
  std::copy( inner.begin(), inner.end(), std::back_inserter( inner_points ) );

  // project along given axis to r

  LINALG::Matrix<3,1> rst;

  std::vector<LINALG::Matrix<3,1> > local_points;
  local_points.reserve( inner_points.size() );
  projected_points.reserve( inner_points.size() );

  for ( std::vector<Point*>::iterator i=inner_points.begin(); i!=inner_points.end(); ++i )
  {
    Point * p = *i;
    LINALG::Matrix<3,1> xyz;
    p->Coordinates( xyz.A() );
    element->LocalCoordinates( xyz, rst );

    local_points.push_back( rst );

    // projection
    rst( axis ) = r;

    // create new points
    element->GlobalCoordinates( rst, xyz );
    projected_points.push_back( mesh.NewPoint( xyz.A(), NULL, NULL ) );
    projected_points.back()->Position( position );
  }

  // create integration cells

  for ( std::vector<Facet*>::iterator i=inner_facets.begin();
        i!=inner_facets.end();
        ++i )
  {
    Facet * f = *i;
#if 1
    if ( f->Equals( DRT::Element::tri3 ) )
    {
      std::vector<Point*> points;
      points.reserve( 6 );

      const std::vector<Point*> & corner_points = f->CornerPoints();

      // find facet orientation in local coordinates
      double drs = 0;
      LINALG::Matrix<3,3> xyze;
      LINALG::Matrix<3,1> normal;
      LINALG::Matrix<2,3> deriv;
      LINALG::Matrix<2,2> metrictensor;

      double * x = xyze.A();

      int sidepos = -1;
      if ( r > 0 )
      {
        sidepos = 0;
        std::copy( corner_points.begin(), corner_points.end(), std::back_inserter( points ) );
      }

      for ( std::vector<Point*>::const_iterator i=corner_points.begin();
            i!=corner_points.end();
            ++i )
      {
        Point * p = *i;

        std::vector<Point*>::iterator pos = std::find( inner_points.begin(), inner_points.end(), p );
        if ( pos==inner_points.end() )
        {
          throw std::runtime_error( "inner point missing" );
        }

        points.push_back( projected_points[pos - inner_points.begin()] );

        const LINALG::Matrix<3,1> & rst = local_points[pos - inner_points.begin()];
        x = std::copy( rst.A(), rst.A()+3, x );
      }

      if ( r < 0 )
      {
        sidepos = 1;
        std::copy( corner_points.begin(), corner_points.end(), std::back_inserter( points ) );
      }

      DRT::UTILS::shape_function_2D_deriv1( deriv, 0., 0., DRT::Element::tri3 );
      DRT::UTILS::ComputeMetricTensorForBoundaryEle<DRT::Element::tri3>( xyze, deriv, metrictensor, drs, &normal );

      if ( normal( axis ) < 0 )
      {
        std::swap( points[1], points[2] );
        std::swap( points[4], points[5] );
      }

      std::vector<Point*> side( 3 );
      for ( int j=0; j<3; ++j )
      {
        side[j] = points[DRT::UTILS::eleNodeNumbering_wedge15_trisurfaces[sidepos][j]];
      }
      Tri3BoundaryCell::CreateCell( mesh, cell, f, side );

      cell->NewWedge6Cell( mesh, points );
    }
    else if ( f->Equals( DRT::Element::quad4 ) )
    {
      std::vector<Point*> points;
      points.reserve( 8 );

      const std::vector<Point*> & corner_points = f->CornerPoints();

      // find facet orientation in local coordinates
      double drs = 0;
      LINALG::Matrix<3,4> xyze;
      LINALG::Matrix<3,1> normal;
      LINALG::Matrix<2,4> deriv;
      LINALG::Matrix<2,2> metrictensor;

      double * x = xyze.A();

      int sidepos = -1;
      if ( r > 0 )
      {
        sidepos = 0;
        std::copy( corner_points.begin(), corner_points.end(), std::back_inserter( points ) );
      }

      for ( std::vector<Point*>::const_iterator i=corner_points.begin();
            i!=corner_points.end();
            ++i )
      {
        Point * p = *i;

        std::vector<Point*>::iterator pos = std::find( inner_points.begin(), inner_points.end(), p );
        if ( pos==inner_points.end() )
        {
          throw std::runtime_error( "inner point missing" );
        }

        points.push_back( projected_points[pos - inner_points.begin()] );

        const LINALG::Matrix<3,1> & rst = local_points[pos - inner_points.begin()];
        x = std::copy( rst.A(), rst.A()+3, x );
      }

      if ( r < 0 )
      {
        sidepos = 5;
        std::copy( corner_points.begin(), corner_points.end(), std::back_inserter( points ) );
      }

      DRT::UTILS::shape_function_2D_deriv1( deriv, 0., 0., DRT::Element::quad4 );
      DRT::UTILS::ComputeMetricTensorForBoundaryEle<DRT::Element::quad4>( xyze, deriv, metrictensor, drs, &normal );

      if ( normal( axis ) < 0 )
      {
        std::swap( points[1], points[3] );
        std::swap( points[5], points[7] );
      }

      std::vector<Point*> side( 4 );
      for ( int j=0; j<4; ++j )
      {
        side[j] = points[DRT::UTILS::eleNodeNumbering_hex27_surfaces[sidepos][j]];
      }
      Quad4BoundaryCell::CreateCell( mesh, cell, f, side );

      cell->NewHex8Cell( mesh, points );
    }
    else
#endif
    {
      std::set<Point*> points;
      const std::vector<Point*> & corner_points = f->CornerPoints();

      f->AllPoints( points );

      for ( std::vector<Point*>::const_iterator i=corner_points.begin();
            i!=corner_points.end();
            ++i )
      {
        Point * p = *i;

        std::vector<Point*>::iterator pos = std::find( inner_points.begin(), inner_points.end(), p );
        if ( pos==inner_points.end() )
        {
          throw std::runtime_error( "inner point missing" );
        }

        points.insert( projected_points[pos - inner_points.begin()] );
      }

      std::vector<Point*> cell_points;
      cell_points.reserve( points.size() );
      std::copy( points.begin(), points.end(), std::back_inserter( cell_points ) );

      // sort points that go into qhull to obtain the same result independent of
      // pointer values (compiler flags, code structure, memory usage, ...)
      std::sort( cell_points.begin(), cell_points.end(), PointPidLess() );

      std::set<Facet*> cell_facets;
      cell_facets.insert( f );
      cell->CreateTet4IntegrationCells( mesh, position, cell_points, cell_facets, true );
    }
  }
}

bool GEO::CUT::Hex8Projection::EdgeCut( Mesh & mesh,
                                        Element * element,
                                        VolumeCell * cell,
                                        Point::PointPosition position,
                                        const std::set<Facet*> & facets,
                                        std::set<IntegrationCell*> & integrationcells,
                                        int cutside1,
                                        int cutside2,
                                        int upside,
                                        int downside )
{
  Facet * cs1 = FindFacet( element, facets, cutside1 );
  Facet * cs2 = FindFacet( element, facets, cutside2 );
  Facet * us  = FindFacet( element, facets, upside );
  Facet * ds  = FindFacet( element, facets, downside );

  std::set<Facet*> myfacets1;
  std::set<Facet*> myfacets2;

  FindNeighborFacets( cs1->Points(), facets, us, ds, myfacets1 );
  FindNeighborFacets( cs2->Points(), facets, us, ds, myfacets2 );

  if ( myfacets1.size()+myfacets2.size()+2 != facets.size() )
    return false;

  std::set<Facet*> reconstruction_facets;
  std::copy( myfacets1.begin(), myfacets1.end(), std::inserter( reconstruction_facets, reconstruction_facets.begin() ) );
  std::copy( myfacets2.begin(), myfacets2.end(), std::inserter( reconstruction_facets, reconstruction_facets.begin() ) );
  reconstruction_facets.insert( us );
  reconstruction_facets.insert( ds );

  if ( reconstruction_facets != facets )
    return false;

  CreateTetMesh( mesh, cell, position, myfacets1 );
  CreateTetMesh( mesh, cell, position, myfacets2 );
  return true;
}

GEO::CUT::Facet * GEO::CUT::Hex8Projection::FindFacet( Element * element, const std::set<Facet*> & facets, int sideid )
{
  const std::vector<Side*> & sides = element->Sides();

  Facet * side_facet = NULL;
  for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;
    if ( f->ParentSide()==sides[sideid] )
    {
      if ( side_facet==NULL )
      {
        side_facet = f;
      }
      else
      {
        throw std::runtime_error( "double side facet" );
      }
    }
  }
  if ( side_facet==NULL )
    throw std::runtime_error( "facet not found" );
  return side_facet;
}

void GEO::CUT::Hex8Projection::FindNeighborFacets( const std::vector<Point*> & facet_points,
                                                   const std::set<Facet*> & facets,
                                                   Facet * us,
                                                   Facet * ds,
                                                   std::set<Facet*> & myfacets )
{
  for ( std::vector<Point*>::const_iterator i=facet_points.begin();
        i!=facet_points.end();
        ++i )
  {
    Point * p = *i;
    const std::set<Facet*> & fs = p->Facets();
    for ( std::set<Facet*>::const_iterator i=fs.begin(); i!=fs.end(); ++i )
    {
      Facet * f = *i;
      if ( f!=us and f!=ds and facets.count( f ) > 0 )
      {
        myfacets.insert( f );
      }
    }
  }
}

void GEO::CUT::Hex8Projection::CreateTetMesh( Mesh & mesh, VolumeCell * cell, Point::PointPosition position, const std::set<Facet*> & myfacets )
{
  std::set<Point*> mypoints;
  for ( std::set<Facet*>::const_iterator i=myfacets.begin(); i!=myfacets.end(); ++i )
  {
    Facet * f = *i;
    const std::vector<Point*> & ps = f->Points();
    std::copy( ps.begin(), ps.end(), std::inserter( mypoints, mypoints.begin() ) );
  }

  std::vector<Point*> points;
  points.reserve( mypoints.size() );
  points.assign( mypoints.begin(), mypoints.end() );

  // sort points that go into qhull to obtain the same result independent of
  // pointer values (compiler flags, code structure, memory usage, ...)
  std::sort( points.begin(), points.end(), PointPidLess() );

  cell->CreateTet4IntegrationCells( mesh, position, points, myfacets, true );
}


#include "cut_projection.H"
#include "cut_volumecell.H"
#include "cut_facet.H"
#include "cut_element.H"
#include "cut_mesh.H"

GEO::CUT::Hex8Projection::Hex8Projection( Mesh & mesh,
                                          Element * element,
                                          VolumeCell * cell,
                                          Point::PointPosition position,
                                          const std::set<Facet*> & facets,
                                          int axis,
                                          double r )
{
  std::set<Point*> cut_points;
  cell->GetAllPoints( mesh, cut_points );

  std::vector<Point*> points;
  points.reserve( cut_points.size() );
  points.assign( cut_points.begin(), cut_points.end() );

  std::cout << "hex8 projection along axis " << axis << " to " << r << "\n";

  // find all inner points that need projecting

  std::set<Point*> inner;

  for ( std::set<Facet*>::const_iterator i=facets.begin(); i!=facets.end(); ++i )
  {
    Facet * f = *i;
    if ( f->OnCutSide() )
    {
      inner_facets_.push_back( f );
      const std::vector<Point*> & points = f->CornerPoints();
      std::copy( points.begin(), points.end(), std::inserter( inner, inner.begin() ) );
    }
  }

  inner_points_.reserve( inner.size() );
  std::copy( inner.begin(), inner.end(), std::back_inserter( inner_points_ ) );

  // project along given axis to r

  LINALG::Matrix<3,1> rst;

  projected_points_.reserve( inner_points_.size() );

  for ( std::vector<Point*>::iterator i=inner_points_.begin(); i!=inner_points_.end(); ++i )
  {
    Point * p = *i;
    LINALG::Matrix<3,1> xyz;
    p->Coordinates( xyz.A() );
    element->LocalCoordinates( xyz, rst );

    // projection
    rst( axis ) = r;

    // create new points
    element->GlobalCoordinates( rst, xyz );
    projected_points_.push_back( mesh.NewPoint( xyz.A(), NULL, NULL ) );
    projected_points_.back()->Position( position );
  }

  // create integration cells

  for ( std::vector<Facet*>::iterator i=inner_facets_.begin();
        i!=inner_facets_.end();
        ++i )
  {
    Facet * f = *i;
#if 0
    if ( f->Equals( DRT::Element::tri3 ) )
    {
    }
    else if ( f->Equals( DRT::Element::quad4 ) )
    {
    }
    else
#endif
    {
      std::vector<Point*> points;
      const std::vector<Point*> & corner_points = f->CornerPoints();

      if ( f->IsTriangulated() )
      {
        const std::vector<std::vector<Point*> > & triangulation = f->Triangulation();
        std::set<Point*> points_set;
        for ( std::vector<std::vector<Point*> >::const_iterator i=triangulation.begin();
              i!=triangulation.end();
              ++i )
        {
          const std::vector<Point*> & tri = *i;
          std::copy( tri.begin(), tri.end(), std::inserter( points_set, points_set.begin() ) );
        }
        points.reserve( points_set.size() + corner_points.size() );
        points.assign( points_set.begin(), points_set.end() );
      }
      else
      {
        const std::vector<Point*> & facet_points = f->Points();
        points.reserve( facet_points.size() + corner_points.size() );
        points.assign( facet_points.begin(), facet_points.end() );
      }

      for ( std::vector<Point*>::const_iterator i=corner_points.begin();
            i!=corner_points.end();
            ++i )
      {
        Point * p = *i;

        std::vector<Point*>::iterator pos = std::find( inner_points_.begin(), inner_points_.end(), p );
        if ( pos==inner_points_.end() )
        {
          throw std::runtime_error( "inner point missing" );
        }

        points.push_back( projected_points_[pos - inner_points_.begin()] );
      }

      // sort points that go into qhull to obtain the same result independent of
      // pointer values (compiler flags, code structure, memory usage, ...)
      std::sort( points.begin(), points.end(), PointPidLess() );

      std::set<Facet*> cell_facets;
      cell_facets.insert( f );
      cell->CreateTet4IntegrationCells( mesh, position, points, cell_facets );
    }
  }
}

#include "cut_kernel.H"
#include "cut_point.H"
#include "cut_position2d.H"

unsigned GEO::CUT::KERNEL::FindNextCornerPoint( const std::vector<Point*> & points,
                                                LINALG::Matrix<3,1> & x1,
                                                LINALG::Matrix<3,1> & x2,
                                                LINALG::Matrix<3,1> & x3,
                                                LINALG::Matrix<3,1> & b1,
                                                LINALG::Matrix<3,1> & b2,
                                                LINALG::Matrix<3,1> & b3,
                                                unsigned i )
{
  unsigned pointsize = points.size();
  unsigned j = ( i+1 ) % pointsize;
  if ( pointsize < 3 )
  {
    return j;
  }

  points[i]->Coordinates( x1.A() );
  points[j]->Coordinates( x2.A() );

  b1.Update( 1, x2, -1, x1, 0 );

  double norm = b1.Norm2();
  if ( norm < std::numeric_limits<double>::min() )
    throw std::runtime_error( "same point in facet not supported" );

  b1.Scale( 1./norm );

  if ( b1.Norm2() < std::numeric_limits<double>::min() )
    throw std::runtime_error( "same point in facet not supported" );

  i = j;
  for ( unsigned k=2; k<pointsize; ++k )
  {
    i = ( i+1 ) % pointsize;
    Point * p = points[i];
    p->Coordinates( x3.A() );

    b2.Update( 1, x3, -1, x1, 0 );

    norm = b2.Norm2();
    if ( norm < std::numeric_limits<double>::min() )
      throw std::runtime_error( "same point in facet not supported" );

    b2.Scale( 1./norm );

    // cross product to get the normal at the point
    b3( 0 ) = b1( 1 )*b2( 2 ) - b1( 2 )*b2( 1 );
    b3( 1 ) = b1( 2 )*b2( 0 ) - b1( 0 )*b2( 2 );
    b3( 2 ) = b1( 0 )*b2( 1 ) - b1( 1 )*b2( 0 );

    if ( b3.Norm2() > PLANARTOL )
    {
      // Found. Return last node on this line.
      return ( i+pointsize-1 ) % pointsize;
    }
  }

  // All on one line. Return first and last point.
  if ( j==0 )
  {
    return 0;
  }
  else
  {
    return pointsize-1;
  }
}

void GEO::CUT::KERNEL::FindCornerPoints( const std::vector<Point*> & points,
                                         std::vector<Point*> & corner_points )
{
  LINALG::Matrix<3,1> x1;
  LINALG::Matrix<3,1> x2;
  LINALG::Matrix<3,1> x3;
  LINALG::Matrix<3,1> b1;
  LINALG::Matrix<3,1> b2;
  LINALG::Matrix<3,1> b3;

  for ( unsigned i = FindNextCornerPoint( points, x1, x2, x3, b1, b2, b3, 0 );
        true;
        i = FindNextCornerPoint( points, x1, x2, x3, b1, b2, b3, i ) )
  {
    Point * p = points[i];
    if ( corner_points.size()>0 and corner_points.front()==p )
      break;
    corner_points.push_back( p );
  }
}

bool GEO::CUT::KERNEL::IsValidQuad4( const std::vector<Point*> & points )
{
  if ( points.size()==4 )
  {
    LINALG::Matrix<3,3> xyze;
    LINALG::Matrix<3,1> xyz;
    for ( int i=0; i<4; ++i )
    {
      points[( i+0 )%4]->Coordinates( &xyze( 0, 0 ) );
      points[( i+1 )%4]->Coordinates( &xyze( 0, 1 ) );
      points[( i+2 )%4]->Coordinates( &xyze( 0, 2 ) );
      points[( i+3 )%4]->Coordinates( &xyz( 0, 0 ) );

      Position2d<DRT::Element::tri3> pos( xyze, xyz );
      if ( pos.Compute() )
      {
        return false;
      }
    }
    return true;
  }
  return false;
}

DRT::Element::DiscretizationType GEO::CUT::KERNEL::CalculateShape( const std::vector<Point*> & points,
                                                                   std::vector<Point*> & line_points )
{
  FindCornerPoints( points, line_points );

  if ( IsValidTri3( line_points ) )
  {
    return DRT::Element::tri3;
  }
  else if ( IsValidQuad4( line_points ) )
  {
    return DRT::Element::quad4;
  }

  return DRT::Element::dis_none;
}

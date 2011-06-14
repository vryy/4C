
#include <iostream>

#include "cut_kernel.H"

bool GEO::CUT::KERNEL::Intersection::operator()( const LINALG::Matrix<3,3> & tri3, const LINALG::Matrix<3,2> & line )
{
  ComputeDistance<DRT::Element::line2, DRT::Element::tri3> cd;

  LINALG::Matrix<3,1> p1( &line( 0, 0 ) );
  LINALG::Matrix<3,1> p2( &line( 0, 1 ) );

  double d1;
  double d2;

  if ( not cd( tri3, p1, d1 ) or not cd( tri3, p2, d2 ) )
  {
    throw std::runtime_error( "failed to calculate tri3--point distance" );
  }

  ComputeIntersection<DRT::Element::line2, DRT::Element::tri3> ci;

  if ( ci( tri3, line ) )
  {
    std::cout << ci.LocalCoordinates();
    return true;
  }

  return false;
}

bool GEO::CUT::KERNEL::Intersection::operator()( const LINALG::Matrix<3,4> & quad4, const LINALG::Matrix<3,2> & line )
{
  ComputeDistance<DRT::Element::line2, DRT::Element::quad4> cd;

  LINALG::Matrix<3,1> p1( &line( 0, 0 ) );
  LINALG::Matrix<3,1> p2( &line( 0, 1 ) );

  ComputeIntersection<DRT::Element::line2, DRT::Element::quad4> ci;

  if ( ci( quad4, line ) )
  {
    return true;
  }

  return false;
}

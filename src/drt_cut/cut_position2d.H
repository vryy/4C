#ifndef CUT_POSITION2D_H
#define CUT_POSITION2D_H

#include <iterator>

#include "../linalg/linalg_gauss.H"

#include "cut_side.H"
#include "cut_kernel.H"

namespace GEO
{
namespace CUT
{

/// side local coordinates of a point
template <DRT::Element::DiscretizationType elementtype,
          int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<elementtype>::numNodePerElement
>
class Position2d
{
public:

  Position2d( ConcreteSide<elementtype> & side,
              const LINALG::Matrix<3,1> & xyz )
    : px_( xyz )
  {
    side.Coordinates( xyze_ );
    Scale();
  }

  Position2d( const LINALG::Matrix<3,numNodesElement> & side_xyze,
              const LINALG::Matrix<3,1> & xyz )
    : px_( xyz ),
      xyze_( side_xyze )
  {
    Scale();
  }

  bool Compute();

  bool WithinLimits();

  const LINALG::Matrix<3,1> & LocalCoordinates() { return xsi_; }

  void LocalCoordinates( LINALG::Matrix<2,1> & rs )
  {
    rs( 0 ) = xsi_( 0 );
    rs( 1 ) = xsi_( 1 );
  }

  double Distance() const { return xsi_( 2, 0 ); }

private:

  bool LinePointDistance( const LINALG::Matrix<3,1> & p1,
                          const LINALG::Matrix<3,1> & p2,
                          const LINALG::Matrix<3,1> & p3,
                          double & distance );

  void Scale()
  {
    scale_ = 0;
    LINALG::Matrix<3,1> d;
    for ( int i=0; i<numNodesElement; ++i )
    {
      LINALG::Matrix<3,1> x1( &xyze_( 0, i ),  true );
      LINALG::Matrix<3,1> x2( &xyze_( 0, ( i+1 )%numNodesElement ),  true );
      d.Update( 1, x2, -1, x1, 0 );
      scale_ += d.Norm2();
    }
    scale_ /= numNodesElement;

    px_.Scale( 1./scale_ );
    xyze_.Scale( 1./scale_ );

    shift_ = 0;
    for ( int i=0; i<numNodesElement; ++i )
    {
      LINALG::Matrix<3,1> x1( &xyze_( 0, i ),  true );
      shift_.Update( 1, x1, 1 );
    }
    shift_.Scale( 1./numNodesElement );

    for ( int i=0; i<numNodesElement; ++i )
    {
      LINALG::Matrix<3,1> x1( &xyze_( 0, i ),  true );
      x1.Update( -1, shift_, 1 );
    }
    px_.Update( -1, shift_, 1 );
  }

  LINALG::Matrix<3,1> px_;

  LINALG::Matrix<3,1> xsi_;

  LINALG::Matrix<3, numNodesElement> xyze_;

  double scale_;

  LINALG::Matrix<3,1> shift_;
};


template <DRT::Element::DiscretizationType elementtype,
          int numNodesElement
>
bool Position2d<elementtype,numNodesElement>::Compute()
{
  KERNEL::ComputeDistance<DRT::Element::line2, elementtype> cd( xsi_ );
  //KERNEL::DebugComputeDistance<DRT::Element::line2, elementtype> cd( xsi_ );

  double dist;
  if ( not cd( xyze_, px_, dist ) )
  {
    switch ( elementtype )
    {
    case DRT::Element::quad4:
    {
      LINALG::Matrix<3,1> p1( &xyze_( 0, 0 ), true );
      LINALG::Matrix<3,1> p2( &xyze_( 0, 1 ), true );
      LINALG::Matrix<3,1> p3( &xyze_( 0, 2 ), true );
      LINALG::Matrix<3,1> p4( &xyze_( 0, 3 ), true );
      double d1, d2, d3, d4;
      if ( LinePointDistance( p1, p2, px_, d1 ) and
           LinePointDistance( p2, p3, px_, d2 ) and
           LinePointDistance( p3, p4, px_, d3 ) and
           LinePointDistance( p4, p1, px_, d4 ) )
      {
        dist = std::min( d1, std::min( d2, std::min( d3, d4 ) ) );
        // The local coordinates are wrong. But can there be right ones here?
        // Do we need them?
        xsi_ = 0;
        xsi_( 2 ) = dist;
        return true;
      }
      else
      {
        return false;
      }
    }
    case DRT::Element::tri3:
    {
      LINALG::Matrix<3,1> p1( &xyze_( 0, 0 ), true );
      LINALG::Matrix<3,1> p2( &xyze_( 0, 1 ), true );
      LINALG::Matrix<3,1> p3( &xyze_( 0, 2 ), true );
      double d1, d2, d3;
      if ( LinePointDistance( p1, p2, px_, d1 ) and
           LinePointDistance( p1, p3, px_, d2 ) and
           LinePointDistance( p2, p3, px_, d3 ) )
      {
        dist = std::min( d1, std::min( d2, d3 ) );
        // The local coordinates are wrong. But can there be right ones here?
        // Do we need them?
        xsi_ = 0;
        xsi_( 2 ) = dist;
        return true;
      }
      else
      {
        return false;
      }
    }
    default:
      throw std::runtime_error( "unsupported side type in position computation" );
    }
    throw std::runtime_error( "failed to calculate side--point distance" );
  }

  return WithinLimits();
}


template <DRT::Element::DiscretizationType elementtype,
          int numNodesElement
>
bool Position2d<elementtype,numNodesElement>::LinePointDistance( const LINALG::Matrix<3,1> & p1,
                                                                 const LINALG::Matrix<3,1> & p2,
                                                                 const LINALG::Matrix<3,1> & p3,
                                                                 double & distance )
{
  LINALG::Matrix<3,1> n;
  n.Update( 1, p2, -1, p1, 0 );

  LINALG::Matrix<3,1> m;
  m.Update( 1, p3, -1, p1, 0 );

  double a = m.Dot( n );
  double b = n.Dot( n );
  double c = a/b;

  if ( c < 0 or c > 1 )
    return false;

  LINALG::Matrix<3,1> d;
  d.Update( 1, p1, c, n, 0 );
  d.Update( -1, p3, 1 );

  distance = d.Norm2();
  return true;
}


template <DRT::Element::DiscretizationType elementtype,
          int numNodesElement
>
bool Position2d<elementtype,numNodesElement>::WithinLimits()
{
  bool intersects = false;
  switch ( elementtype )
  {
  case DRT::Element::quad4:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
  {
    intersects = ( xsi_( 0 ) >= -1-TOLERANCE and xsi_( 1 ) >= -1-TOLERANCE and
                   xsi_( 0 ) <=  1+TOLERANCE and xsi_( 1 ) <=  1+TOLERANCE );
    break;
  }
  case DRT::Element::tri3:
  case DRT::Element::tri6:
  {
    intersects = ( xsi_( 0 ) >= 0-TOLERANCE and xsi_( 1 ) >= 0-TOLERANCE and
                   xsi_( 0 ) <= 1+TOLERANCE and xsi_( 1 ) <= 1+TOLERANCE and
                   xsi_( 1 ) <= ( 1 - xsi_( 0 ) )+TOLERANCE );
    break;
  }
  default:
    throw std::runtime_error( "unsupported side type" );
  }
  return intersects and fabs( xsi_( 2 ) ) < TOLERANCE;
}

}
}

#endif

/*---------------------------------------------------------------------*/
/*!
\file cut_position2d.H

\brief compute local coords(position) on 2d element and decide if the point is inside

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

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

  /*!
  \brief Constructor
   */
  Position2d( ConcreteSide<elementtype> & side,
              const LINALG::Matrix<3,1> & xyz )
    : px_( xyz )
  {
    side.Coordinates( xyze_ );
    Scale();
    ConstructBoundingBoxSide();
  }

  /*!
  \brief Constructor
   */
  Position2d( const LINALG::Matrix<3,numNodesElement> & side_xyze,
              const LINALG::Matrix<3,1> & xyz )
    : px_( xyz ),
      xyze_( side_xyze )
  {
    Scale();
    ConstructBoundingBoxSide();
  }

  bool Compute(bool allow_dist = false);

  /*!
  \brief Returns true if the given point is within the given side. Internally local coordinates of a point is also computed
   */
  bool IsGivenPointWithinSide();

  /*!
  \brief Returns true if the local coordinate of the given is within the limiting natural coordinate system of the given Side
   */
  bool WithinLimits(bool allow_dist);

  /*!
  \brief Return the local coordinates and the distance of the given point to the side
   */
  const LINALG::Matrix<3,1> & LocalCoordinates() { return xsi_; }

  /*!
  \brief Return the local coordinates of the given point
   */
  void LocalCoordinates( LINALG::Matrix<2,1> & rs )
  {
    rs( 0 ) = xsi_( 0 );
    rs( 1 ) = xsi_( 1 );
  }

  /*!
  \brief Return the perpendicular distance between given point to the side
   */
  double Distance() const { return xsi_( 2, 0 ); }

private:

  /*!
  \brief Compute the perpendicular distance between point p3 and the line joining p1 and p2
   */
  bool LinePointDistance( const LINALG::Matrix<3,1> & p1,
                          const LINALG::Matrix<3,1> & p2,
                          const LINALG::Matrix<3,1> & p3,
                          double & distance );

  /*!
  \brief The corner points of the side and the location of the point are converted into the corresponding range of natural
  coordinate system of the element. This also gives the proper initial choice for the Newton scheme to find the local coordinate system
   */
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

  /*!
  \brief Construct bounding box over the given side (after scaling)
   */
  void ConstructBoundingBoxSide()
  {
    for ( int i=0; i<numNodesElement; ++i )
    {
      LINALG::Matrix<3,1> x1( &xyze_( 0, i ),  true );
      bbside_.AddPoint( x1 );
    }
  }

  //! Coordinate of the given point (scaled immediately after construction)
  LINALG::Matrix<3,1> px_;

  //! Local coordinate of the given point and the normal distance of the point with respect to side
  // xsi_(0) and xsi_(1)  --> Local coordinate
  // xsi_(2) --> Distance to the side
  LINALG::Matrix<3,1> xsi_;

  //! Coordinate of the given Side (scaled immediately after construction)
  LINALG::Matrix<3, numNodesElement> xyze_;

  //! Average length of side edges
  double scale_;

  //! Each node of the sides is shifted by this amount during scaling
  LINALG::Matrix<3,1> shift_;

  //! Bounding box over the given side (after scaling is performed)
  BoundingBox bbside_;
};

template <DRT::Element::DiscretizationType elementtype,
          int numNodesElement
>
bool Position2d<elementtype,numNodesElement>::Compute(bool allow_dist)
{
  // If the given point is outside the side's bounding box, no need to perform the complex calculations
  //TODO: In several places within cut, this funciton is used without properly checking whether the point is within
  // the side or not. Has to be checked thoroughly
  /*if( not allow_dist )
  {
    if( not bbside_.Within( 1.0, &px_(0,0) ))
      return false;
    xsi_ = 0.0;
  }*/

  // try to compute the local coordinates and the distance using the Newton scheme
  KERNEL::ComputeDistance<DRT::Element::line2, elementtype> cd( xsi_ );
  //KERNEL::DebugComputeDistance<DRT::Element::line2, elementtype> cd( xsi_ );

  double dist = 0.0;
  if ( not cd( xyze_, px_, dist ) ) // if newton did not converge, do a simpler check
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

  return WithinLimits(allow_dist);
 }

/*-------------------------------------------------------------------------------------------------*
 * Return true if the given point is within the given side.
 * Internally it computes the local-coordinates of the given point with respect to
 * the given sides
 *-------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType elementtype,
          int numNodesElement
>
bool Position2d<elementtype,numNodesElement>::IsGivenPointWithinSide()
{
  // If the given point is outside the side's bounding box, no need to perform the complex calculations
  if( not bbside_.Within( 1.0, &px_(0,0) ))
    return false;

  // try to compute the local coordinates and the distance using the Newton scheme
  KERNEL::ComputeDistance<DRT::Element::line2, elementtype> cd( xsi_ );
  //KERNEL::DebugComputeDistance<DRT::Element::line2, elementtype> cd( xsi_ );

  double dist = 0.0;

  if ( not cd( xyze_, px_, dist ) ) return false;

  return true;
 }


/*-------------------------------------------------------------------------------------------------------*
 * Perpendicular distance between the given point (p3) to the line joining p1 and p2
 *-------------------------------------------------------------------------------------------------------*/
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
bool Position2d<elementtype,numNodesElement>::WithinLimits(bool allow_dist)
{
  bool intersects = false;
  switch ( elementtype )
  {
  case DRT::Element::quad4:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
  {
    intersects = ( xsi_( 0 ) >= -1-REFERENCETOL and xsi_( 1 ) >= -1-REFERENCETOL and
                   xsi_( 0 ) <=  1+REFERENCETOL and xsi_( 1 ) <=  1+REFERENCETOL );
    break;
  }
  case DRT::Element::tri3:
  case DRT::Element::tri6:
  {
    intersects = ( xsi_( 0 ) >= 0-REFERENCETOL and xsi_( 1 ) >= 0-REFERENCETOL and
                   xsi_( 0 ) <= 1+REFERENCETOL and xsi_( 1 ) <= 1+REFERENCETOL and
                   xsi_( 1 ) <= ( 1 - xsi_( 0 ) )+REFERENCETOL );
    break;
  }
  default:
    throw std::runtime_error( "unsupported side type" );
  }
  if (!intersects) return false;
  else if(allow_dist) return intersects;
  else
  {
    double tol = xyze_.NormInf();
    double tol2 = px_.NormInf();
    if (tol2 > tol)
      tol = tol2;

    tol *= POSITIONTOL;
    return (intersects and fabs( xsi_( 2 ) ) < tol);
  }
}

}
}

#endif

/*---------------------------------------------------------------------*/
/*!
\file cut_position.H

\brief Class to check whether a point lies inside an element, and also to
compute local coordinates of a point with respect to the element

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#ifndef CUT_POSITION_H
#define CUT_POSITION_H

#include "cut_element.H"
#include "cut_kernel.H"

namespace GEO
{
namespace CUT
{
template <DRT::Element::DiscretizationType elementtype,
          int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<elementtype>::numNodePerElement
>
class Position
{
public:

  Position( ConcreteElement<elementtype> & element,
            Point & point )
  {
    element.Coordinates( xyze_ );
    point.Coordinates( px_.A() );
    Scale();
  }

  Position( ConcreteElement<elementtype> & element,
            const LINALG::Matrix<3,1> & xyz )
    : px_( xyz )
  {
    element.Coordinates( xyze_ );
    Scale();
  }

  Position( const LINALG::Matrix<3,numNodesElement> & xyze,
            const LINALG::Matrix<3,1> & xyz )
    : px_( xyz ),
      xyze_( xyze )
  {
    Scale();
  }

  Position( const std::vector<Node*> nodes,
            const LINALG::Matrix<3,1> & xyz )
      :px_( xyz )
  {
    if( not (nodes.size() == numNodesElement) )
      dserror( "nodes for this element are not correct\n" );

    for ( int i=0; i<numNodesElement; ++i )
    {
      Node * n = nodes[i];
      n->Coordinates( &xyze_( 0, i ) );
    }

    Scale();
  }

  /*!
  \brief Computes local coordinates of the point and returns true if it is within the element
  Uses tolerance set in cut_tolerance.H
   */
  bool Compute();

  /*!
  \brief Computes local coordinates of the point and returns true if it is within the element
   */
  bool ComputeTol(double Tol);

  /*!
  \brief Returns true if local coordinates of the point is within range defined for the element.
  Uses tolerance set in cut_tolerance.H
   */
  bool WithinLimits();

  /*!
  \brief Returns true if local coordinates of the point is within range defined for the element
   */
  bool WithinLimitsTol(double Tol);

  /*!
  \brief Returns local coordinates of the point
   */
  const LINALG::Matrix<3,1> & LocalCoordinates() { return xsi_; }

private:

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
  }

  LINALG::Matrix<3,1> px_;

  LINALG::Matrix<3,1> xsi_;

  LINALG::Matrix<3, numNodesElement> xyze_;

  double scale_;
};


template <DRT::Element::DiscretizationType elementtype,
          int numNodesElement
>
bool Position<elementtype,numNodesElement>::Compute()
{
//  KERNEL::DebugComputePosition<elementtype> cp( xsi_ );
  KERNEL::ComputePosition<elementtype> cp( xsi_ );

  if ( cp( xyze_, px_ ) )
  {
    return WithinLimits();
  }

  return false;
}

template <DRT::Element::DiscretizationType elementtype,
          int numNodesElement
>
bool Position<elementtype,numNodesElement>::ComputeTol(double Tol)
{
//  KERNEL::DebugComputePosition<elementtype> cp( xsi_ );
  KERNEL::ComputePosition<elementtype> cp( xsi_ );
  if ( cp( xyze_, px_ ) )
  {
    return WithinLimitsTol(Tol);
  }

  return false;
}


template <DRT::Element::DiscretizationType elementtype,
          int numNodesElement
>
bool Position<elementtype,numNodesElement>::WithinLimits()
{
  switch ( elementtype )
  {
  case DRT::Element::hex8:
  case DRT::Element::hex16:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
    return ( xsi_( 0 ) >= -1-REFERENCETOL and xsi_( 1 ) >= -1-REFERENCETOL and xsi_( 2 ) >= -1-REFERENCETOL and
             xsi_( 0 ) <=  1+REFERENCETOL and xsi_( 1 ) <=  1+REFERENCETOL and xsi_( 2 ) <=  1+REFERENCETOL );
  case DRT::Element::tet4:
  case DRT::Element::tet10:
    return ( xsi_( 0 ) >=   -REFERENCETOL and xsi_( 1 ) >=   -REFERENCETOL and xsi_( 2 ) >=   -REFERENCETOL and
             xsi_( 0 )+xsi_( 1 )+xsi_( 2 ) <= 1+REFERENCETOL );
  case DRT::Element::pyramid5:
    return ( xsi_( 0 ) >= -1-REFERENCETOL and xsi_( 1 ) >= -1-REFERENCETOL and xsi_( 2 ) >=   -REFERENCETOL and
             xsi_( 0 ) <=  1+REFERENCETOL and xsi_( 1 ) <=  1+REFERENCETOL and
             ( fabs( xsi_( 0 ) ) > fabs( xsi_( 1 ) ) ?
               ( xsi_( 2 )+fabs( xsi_( 0 ) ) <= 1+REFERENCETOL ) :
               ( xsi_( 2 )+fabs( xsi_( 1 ) ) <= 1+REFERENCETOL ) )
      );
  case DRT::Element::wedge6:
  case DRT::Element::wedge15:
    return ( xsi_( 0 ) >=   -REFERENCETOL and xsi_( 1 ) >=   -REFERENCETOL and xsi_( 2 ) >= -1-REFERENCETOL and
             xsi_( 2 ) <=  1+REFERENCETOL and
             xsi_( 0 )+xsi_( 1 ) <= 1+REFERENCETOL );
  default:
    throw std::runtime_error( "unsupported element type" );
  }
  return false;
}

template <DRT::Element::DiscretizationType elementtype,
          int numNodesElement
>
bool Position<elementtype,numNodesElement>::WithinLimitsTol(double Tol)
{
  switch ( elementtype )
  {
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
    return ( xsi_( 0 ) >= -1-Tol and xsi_( 1 ) >= -1-Tol and xsi_( 2 ) >= -1-Tol and
             xsi_( 0 ) <=  1+Tol and xsi_( 1 ) <=  1+Tol and xsi_( 2 ) <=  1+Tol );
  case DRT::Element::tet4:
  case DRT::Element::tet10:
    return ( xsi_( 0 ) >=   -Tol and xsi_( 1 ) >=   -Tol and xsi_( 2 ) >=   -Tol and
             xsi_( 0 )+xsi_( 1 )+xsi_( 2 ) <= 1+Tol );
  case DRT::Element::pyramid5:
    return ( xsi_( 0 ) >= -1-Tol and xsi_( 1 ) >= -1-Tol and xsi_( 2 ) >=   -Tol and
             xsi_( 0 ) <=  1+Tol and xsi_( 1 ) <=  1+Tol and
             ( fabs( xsi_( 0 ) ) > fabs( xsi_( 1 ) ) ?
               ( xsi_( 2 )+fabs( xsi_( 0 ) ) <= 1+Tol ) :
               ( xsi_( 2 )+fabs( xsi_( 1 ) ) <= 1+Tol ) )
      );
  case DRT::Element::wedge6:
  case DRT::Element::wedge15:
    return ( xsi_( 0 ) >=   -Tol and xsi_( 1 ) >=   -Tol and xsi_( 2 ) >= -1-Tol and
             xsi_( 2 ) <=  1+Tol and
             xsi_( 0 )+xsi_( 1 ) <= 1+Tol );
  default:
    throw std::runtime_error( "unsupported element type" );
  }
  return false;
}

}
}

#endif

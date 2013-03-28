#ifndef CUT_INTERSECTION_H
#define CUT_INTERSECTION_H

#include "cut_side.H"
#include "cut_mesh.H"
#include "cut_kernel.H"

namespace GEO
{
namespace CUT
{

/*!
\brief A class to calculate the intersection of an edge with a side. The core class where the all the
cut points are actually calculated
 */
template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine    = DRT::UTILS::DisTypeToNumNodePerEle<edgetype>::numNodePerElement,
          int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<sidetype>::numNodePerElement
>
class IntersectionBase
{
public:

  IntersectionBase(  LINALG::Matrix<3, numNodesSurface> xyze_surfaceElement,
                     LINALG::Matrix<3, numNodesLine> xyze_lineElement )
    : xyze_lineElement_(xyze_lineElement),
      xyze_surfaceElement_(xyze_surfaceElement),
      xsi_( true ),
      x_( true )
  {

#if 1
    scale_ = 1;
#else
    scale_ = 0;
    LINALG::Matrix<3,1> d;
    for ( int i=0; i<numNodesSurface; ++i )
    {
      LINALG::Matrix<3,1> x1( &xyze_surfaceElement_( 0, i ),  true );
      LINALG::Matrix<3,1> x2( &xyze_surfaceElement_( 0, ( i+1 )%numNodesSurface ),  true );
      d.Update( 1, x2, -1, x1, 0 );
      scale_ += d.Norm2();
    }
    scale_ /= numNodesSurface;

    xyze_lineElement_.Scale( 1./scale_ );
    xyze_surfaceElement_.Scale( 1./scale_ );
#endif
  }

  bool ComputeCurveSurfaceIntersection();

  bool AtEdge();

  bool SurfaceWithinLimits();

  bool LineWithinLimits();


private:

  LINALG::Matrix<3, numNodesLine> xyze_lineElement_;
  LINALG::Matrix<3, numNodesSurface> xyze_surfaceElement_;

  LINALG::Matrix<3,1> xsi_;
  LINALG::Matrix<3,1> x_;

  double scale_;

};



template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool IntersectionBase<edgetype, sidetype, numNodesLine, numNodesSurface>::ComputeCurveSurfaceIntersection()
{
  KERNEL::ComputeIntersection<edgetype, sidetype> ci( xsi_ );
  //KERNEL::DebugComputeIntersection<edgetype, sidetype> ci( xsi_ );


  if ( ci( xyze_surfaceElement_, xyze_lineElement_ ) )
  {
    return true;
  }

  return false;
}


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool IntersectionBase<edgetype, sidetype, numNodesLine, numNodesSurface>::AtEdge()
{
  switch ( sidetype )
  {
  case DRT::Element::quad4:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
  {
    return ( fabs( xsi_( 0 )+1 )<TOLERANCE or fabs( xsi_( 1 )+1 )<TOLERANCE or
             fabs( xsi_( 0 )-1 )<TOLERANCE or fabs( xsi_( 1 )-1 )<TOLERANCE );
  }
  case DRT::Element::tri3:
  case DRT::Element::tri6:
  {
    return ( fabs( xsi_( 0 )+0 )<TOLERANCE or fabs( xsi_( 1 )+0 )<TOLERANCE or
             fabs( xsi_( 1 )+xsi_( 0 )-1 )<TOLERANCE );
  }
  default:
    throw std::runtime_error( "unsupported side type" );
  }
}


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool IntersectionBase<edgetype, sidetype, numNodesLine, numNodesSurface>::SurfaceWithinLimits()
{
  bool intersects = false;
  switch ( sidetype )
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
  return intersects;
}


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool IntersectionBase<edgetype, sidetype, numNodesLine, numNodesSurface>::LineWithinLimits()
{
  return xsi_( 2 ) >= -1-TOLERANCE and xsi_( 2 ) <=  1+TOLERANCE;
}



/*!
\brief A class to calculate the intersection of an edge with a side. The core class where the all the
cut points are actually calculated
 */
template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine    = DRT::UTILS::DisTypeToNumNodePerEle<edgetype>::numNodePerElement,
          int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<sidetype>::numNodePerElement
>
class Intersection
{

public:

  Intersection( Mesh & mesh,
                Edge & edge,
                ConcreteSide<sidetype> & side )
    : mesh_( mesh ),
      edge_( edge ),
      side_( side ),
      xsi_( true ),
      x_( true )
  {
    edge_.Coordinates( xyze_lineElement_ );
    side_.Coordinates( xyze_surfaceElement_ );

#if 1
    scale_ = 1;
#else
    scale_ = 0;
    LINALG::Matrix<3,1> d;
    for ( int i=0; i<numNodesSurface; ++i )
    {
      LINALG::Matrix<3,1> x1( &xyze_surfaceElement_( 0, i ),  true );
      LINALG::Matrix<3,1> x2( &xyze_surfaceElement_( 0, ( i+1 )%numNodesSurface ),  true );
      d.Update( 1, x2, -1, x1, 0 );
      scale_ += d.Norm2();
    }
    scale_ /= numNodesSurface;

    xyze_lineElement_.Scale( 1./scale_ );
    xyze_surfaceElement_.Scale( 1./scale_ );
#endif
  }

  /*!
  \brief Computes the intersection points of the edge with the specified side and store the points in cuts
   */
  bool Intersect( PointSet & cuts );

  bool ComputeCurveSurfaceIntersection();

  bool AtEdge();

  bool SurfaceWithinLimits();

  bool LineWithinLimits();

  LINALG::Matrix<3,1> & FinalPoint();

private:

  double ComputeDistance( Point * p );

  bool ComputeCut( Edge * e1, Edge * e2 );

  /// add cut point that is a node to all edges and sides it touches
  void InsertCut( Node * n, PointSet & cuts )
  {
    cuts.insert( Point::InsertCut( &edge_, &side_, n ) );
  }

  /// test if we are on any edge of the cut side and add those edges
  void TestSideEdges( Point * p, const LINALG::Matrix<3,1> & xsi )
  {
    std::vector<Edge*> edges;
    TestSideEdges( p, xsi, edges );
  }

  void TestSideEdges( Point * p, const LINALG::Matrix<3,1> & xsi, std::vector<Edge*> & edges )
  {
    if ( AtEdge() )
    {
      side_.EdgeAt( xsi( 0, 0 ), xsi( 1, 0 ), edges );
      for ( std::vector<Edge*>::iterator i=edges.begin(); i!=edges.end(); ++i )
      {
        Edge * e = *i;
        p->AddEdge( e );
      }
    }
  }

  Mesh & mesh_;
  Edge & edge_;
  ConcreteSide<sidetype> & side_;

  LINALG::Matrix<3, numNodesSurface> xyze_surfaceElement_;
  LINALG::Matrix<3, numNodesLine> xyze_lineElement_;

  LINALG::Matrix<3,1> xsi_;
  LINALG::Matrix<3,1> x_;

  double scale_;
};


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::Intersect( PointSet & cuts )
{
  bool begin_within_limits = false;
  bool end_within_limits = false;

  LINALG::Matrix<3,1> begin_xsi;
  LINALG::Matrix<3,1> end_xsi;

  if ( ComputeDistance( edge_.BeginNode()->point() ) < TOLERANCE and SurfaceWithinLimits() )
  {
    begin_xsi = xsi_;
    begin_within_limits = true;
    if ( ComputeDistance( edge_.EndNode()->point() ) < TOLERANCE and SurfaceWithinLimits() )
    {
      InsertCut( edge_.BeginNode(), cuts );
      InsertCut( edge_.EndNode()  , cuts );
      TestSideEdges( edge_.BeginNode()->point(), begin_xsi );
      TestSideEdges( edge_.EndNode()  ->point(), xsi_ );
      return true;
    }
  }
  else if ( ComputeDistance( edge_.EndNode()->point() ) < TOLERANCE and SurfaceWithinLimits() )
  {
    end_xsi = xsi_;
    end_within_limits = true;
  }

  if ( ComputeCurveSurfaceIntersection() )
  {
    if ( SurfaceWithinLimits() and LineWithinLimits() )
    {
      FinalPoint();

      if ( fabs( xsi_( 2, 0 )+1 ) < TOLERANCE )
      {
        InsertCut( edge_.BeginNode(), cuts );
        TestSideEdges( edge_.BeginNode()->point(), xsi_ );
        return true;
      }
      else if ( fabs( xsi_( 2, 0 )-1 ) < TOLERANCE )
      {
        InsertCut( edge_.EndNode(), cuts );
        TestSideEdges( edge_.EndNode()->point(), xsi_ );
        return true;
      }
      else
      {
        Node * n = side_.OnNode( x_ );
        if ( n!=NULL )
        {
          InsertCut( n, cuts );
          TestSideEdges( n->point(), xsi_ );
          return true;
        }
        else
        {
          Point * p = Point::NewPoint( mesh_, x_.A(), xsi_( 2, 0 ), &edge_, &side_ );
          TestSideEdges( p, xsi_ );
          cuts.insert( p );
          return true;
        }
      }
    }
  }
  else
  {
    bool success = false;
    Point * nodalpoint = NULL;

    if ( begin_within_limits )
    {
      InsertCut( edge_.BeginNode(), cuts );

      std::vector<Edge*> cut_edges;
      nodalpoint = edge_.BeginNode()->point();
      TestSideEdges( nodalpoint, begin_xsi, cut_edges );

      success = true;
    }
    else if ( end_within_limits )
    {
      InsertCut( edge_.EndNode(), cuts );

      std::vector<Edge*> cut_edges;
      nodalpoint = edge_.EndNode()->point();
      TestSideEdges( nodalpoint, end_xsi, cut_edges );

      success = true;
    }

    // Search all side edges for a cut. We could cross the side.

    const std::vector<Edge*> & side_edges = side_.Edges();
    for ( std::vector<Edge*>::const_iterator i=side_edges.begin(); i!=side_edges.end(); ++i )
    {
      Edge * e = *i;

      if ( nodalpoint!=NULL and nodalpoint->IsCut( e ) )
      {
        // No need to do anything, since the side has a closed cycle of
        // edges. Any matching nodes will be found.
        success = true;
        continue;
      }

      PointSet cut_points;
      edge_.GetCutPoints( e, cut_points );
      if ( cut_points.size() > 0 )
      {
        // Nothing to be done. There are cut points between these edges. We
        // cannot find new ones.
      }
      else if ( ComputeCut( e, &edge_ ) )
      {
        Point * p = Point::NewPoint( mesh_, x_.A(), xsi_( 2 ), &edge_, &side_ );
        p->AddEdge( e );
        cuts.insert( p );
        success = true;
      }
      else
      {
      }
    }

    return success;
  }

  return false;
}


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::ComputeCurveSurfaceIntersection()
{
  KERNEL::ComputeIntersection<edgetype, sidetype> ci( xsi_ );
  //KERNEL::DebugComputeIntersection<edgetype, sidetype> ci( xsi_ );


  if ( ci( xyze_surfaceElement_, xyze_lineElement_ ) )
  {
    return true;
  }

  return false;
}


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
double Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::ComputeDistance( Point * p )
{
  KERNEL::ComputeDistance<edgetype, sidetype> cd( xsi_ );
  //KERNEL::DebugComputeDistance<edgetype, sidetype> cd;
  LINALG::Matrix<3,1> point( p->X() );

  double dist;
  if ( not cd( xyze_surfaceElement_, point, dist ) )
  {
    throw std::runtime_error( "failed to calculate side--point distance" );
  }

  return dist;
}


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::ComputeCut( Edge * e1, Edge * e2 )
{
  return e2->ComputeCut( e1, xsi_( 2 ), x_ );
}


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::AtEdge()
{
  switch ( sidetype )
  {
  case DRT::Element::quad4:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
  {
    return ( fabs( xsi_( 0 )+1 )<TOLERANCE or fabs( xsi_( 1 )+1 )<TOLERANCE or
             fabs( xsi_( 0 )-1 )<TOLERANCE or fabs( xsi_( 1 )-1 )<TOLERANCE );
  }
  case DRT::Element::tri3:
  case DRT::Element::tri6:
  {
    return ( fabs( xsi_( 0 )+0 )<TOLERANCE or fabs( xsi_( 1 )+0 )<TOLERANCE or
             fabs( xsi_( 1 )+xsi_( 0 )-1 )<TOLERANCE );
  }
  default:
    throw std::runtime_error( "unsupported side type" );
  }
}


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::SurfaceWithinLimits()
{
  bool intersects = false;
  switch ( sidetype )
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
  return intersects;
}


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::LineWithinLimits()
{
  return xsi_( 2 ) >= -1-TOLERANCE and xsi_( 2 ) <=  1+TOLERANCE;
}


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
LINALG::Matrix<3,1> & Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::FinalPoint()
{
  // get final point
  x_ = 0;
  LINALG::Matrix<numNodesLine,1> lineFunct;
  DRT::UTILS::shape_function_1D(lineFunct, xsi_(2), edgetype);
  for ( int inode=0; inode<numNodesLine; ++inode )
    for ( int isd=0; isd<3; ++isd )
    {
      x_(isd) += xyze_lineElement_(isd,inode) * lineFunct(inode);
    }
  x_.Scale( scale_ );
  return x_;
}

}
}

#endif

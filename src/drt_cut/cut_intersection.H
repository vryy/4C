/*---------------------------------------------------------------------*/
/*!
\file cut_intersection.H

\brief here the intersection of a (plane) surface with a line is performed

\level 2

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#ifndef CUT_INTERSECTION_H
#define CUT_INTERSECTION_H

#include "cut_side.H"
#include "cut_mesh.H"
#include "cut_kernel.H"
#include "cut_output.H"

#include <Teuchos_TimeMonitor.hpp>

//activate debug intersection output in cut_intersection
//#define DEBUG_CUTINTERSECTION_OUTPUT

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
  {
    xyze_lineElement_ = xyze_lineElement;
    xyze_surfaceElement_ = xyze_surfaceElement;
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

  IntersectionBase( Edge & edge,
                    ConcreteSide<sidetype> & side )
  :scale_(1)
{
  edge.Coordinates( xyze_lineElement_ );
  side.Coordinates( xyze_surfaceElement_ );
}

  bool ComputeCurveSurfaceIntersection(double& tolerance, int triangleid = 0 );

  bool AtEdge(const LINALG::Matrix<3,1> & xsi);

  bool SurfaceWithinLimits(double tol = REFERENCETOL, DRT::Element::DiscretizationType checksidetype = sidetype);

  bool LineWithinLimits(double tol = REFERENCETOL);


protected:

  void GetTriangle(LINALG::Matrix<3,3>& xyze_triElement, int triangleid);

  //
  static LINALG::Matrix<3, numNodesLine> xyze_lineElement_;
  static LINALG::Matrix<3, numNodesSurface> xyze_surfaceElement_;

  static LINALG::Matrix<3,1> xsi_;
  static LINALG::Matrix<3,1> x_;

  double scale_;

};

template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool IntersectionBase<edgetype, sidetype, numNodesLine, numNodesSurface>::ComputeCurveSurfaceIntersection(double& tolerance, int triangleid)
{
  TEUCHOS_FUNC_TIME_MONITOR( "ComputeCurveSurfaceIntersection" );

  bool conv;

  if (!triangleid)
  {
    KERNEL::ComputeIntersection<edgetype, sidetype> ci( xsi_ );
    //KERNEL::DebugComputeIntersection<edgetype, sidetype> ci( xsi_ );
    conv =  ci( xyze_surfaceElement_, xyze_lineElement_ );
    tolerance = ci.GetTolerance();
  }
  else
  {
    KERNEL::ComputeIntersection<edgetype,  DRT::Element::tri3> ci( xsi_ );

    LINALG::Matrix<3, 3> xyze_triElement;
    GetTriangle(xyze_triElement,triangleid);

    conv =  ci( xyze_triElement, xyze_lineElement_ );
    tolerance = ci.GetTolerance();
  }

return conv;
}

template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
//get triangle coords with triangleid = 1 or triangleid = 2
void IntersectionBase<edgetype, sidetype, numNodesLine, numNodesSurface>::GetTriangle(LINALG::Matrix<3,3>& xyze_triElement, int triangleid)
{
  if (sidetype == DRT::Element::quad4)
 {
    //here it is important that the triangle is created in the same rotation as the quad is, to get normal in the same direction and therefore the same signed distance!!!
     int n0 = 0;
     switch (triangleid)
     {
     case 1 :
       n0 = 0;
       break;
     case 2 :
       n0 = 2;
       break;
     default :
       dserror("There should be just two triangles, but triangleid = %d",triangleid);
     }

     //guess this can be done more efficient! - feel free!
     for (uint dim = 0; dim < 3; ++dim)
       for (uint nidx = 0; nidx < 3 ; nidx++)
       {
         xyze_triElement(dim,nidx) = xyze_surfaceElement_(dim,( (nidx+n0) % 4 ) );
       }
  }
  else
  {
    dserror("Cut::Intersection::ComputDistance: For Triangulation a quad4 is expected!");
  }
}

template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool IntersectionBase<edgetype, sidetype, numNodesLine, numNodesSurface>::AtEdge(const LINALG::Matrix<3,1> & xsi)
{
  switch ( sidetype )
  {
  case DRT::Element::quad4:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
  {
    return ( fabs( xsi( 0 )+1 )<REFERENCETOL or fabs( xsi( 1 )+1 )<REFERENCETOL or
             fabs( xsi( 0 )-1 )<REFERENCETOL or fabs( xsi( 1 )-1 )<REFERENCETOL );
  }
  case DRT::Element::tri3:
  case DRT::Element::tri6:
  {
    return ( fabs( xsi( 0 )+0 )<REFERENCETOL or fabs( xsi( 1 )+0 )<REFERENCETOL or
             fabs( xsi( 1 )+xsi( 0 )-1 )<REFERENCETOL );
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
bool IntersectionBase<edgetype, sidetype, numNodesLine, numNodesSurface>::SurfaceWithinLimits( double tolerance, DRT::Element::DiscretizationType checksidetype)
{
  bool intersects = false;
  switch ( checksidetype )
  {
  case DRT::Element::quad4:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
  {
    intersects = ( xsi_( 0 ) >= -1-tolerance and xsi_( 1 ) >= -1-tolerance and
                   xsi_( 0 ) <=  1+tolerance and xsi_( 1 ) <=  1+tolerance );
    break;
  }
  case DRT::Element::tri3:
  case DRT::Element::tri6:
  {
    intersects = ( xsi_( 0 ) >= 0-tolerance and xsi_( 1 ) >= 0-tolerance and
                   xsi_( 0 ) <= 1+tolerance and xsi_( 1 ) <= 1+tolerance and
                   xsi_( 1 ) <= ( 1 - xsi_( 0 ) )+tolerance );
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
bool IntersectionBase<edgetype, sidetype, numNodesLine, numNodesSurface>::LineWithinLimits(double tolerance)
{
  return xsi_( 2 ) >= -1-tolerance and xsi_( 2 ) <=  1+tolerance;
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
class Intersection : public IntersectionBase<edgetype,sidetype,numNodesLine,numNodesSurface>
{

public:

  Intersection( Mesh & mesh,
                Edge & edge,
                ConcreteSide<sidetype> & side )
    :   IntersectionBase<edgetype,sidetype,numNodesLine,numNodesSurface>(edge , side),
      mesh_( mesh ),
      edge_( edge ),
      side_( side ),
      xyze_surfaceElement_( IntersectionBase<edgetype,sidetype,numNodesLine,numNodesSurface>::xyze_surfaceElement_ ),
      xyze_lineElement_( IntersectionBase<edgetype,sidetype,numNodesLine,numNodesSurface>::xyze_lineElement_ ),
      xsi_( IntersectionBase<edgetype,sidetype,numNodesLine,numNodesSurface>::xsi_ ),
      x_( IntersectionBase<edgetype,sidetype,numNodesLine,numNodesSurface>::x_ ),
      scale_( IntersectionBase<edgetype,sidetype,numNodesLine,numNodesSurface>::scale_ )
  {  }

  /*!
  \brief Computes the intersection points of the edge with the specified side and store the points in cuts
   */
  bool Intersect( PointSet & cuts );

  LINALG::Matrix<3,1> & FinalPoint();

  bool RefinedBBOverlapCheck(int maxstep = 10);

  bool ComputeCurveSurfaceIntersection(double& tolerance, int triangleid = 0 )
  {return IntersectionBase<edgetype,sidetype,numNodesLine,numNodesSurface>::ComputeCurveSurfaceIntersection(tolerance,triangleid); }

  bool AtEdge(const LINALG::Matrix<3,1> & xsi)
  {return IntersectionBase<edgetype,sidetype,numNodesLine,numNodesSurface>::AtEdge(xsi); }

  bool SurfaceWithinLimits(double tol = REFERENCETOL, DRT::Element::DiscretizationType checksidetype = sidetype)
  {return IntersectionBase<edgetype,sidetype,numNodesLine,numNodesSurface>::SurfaceWithinLimits(tol,checksidetype); }

  bool LineWithinLimits(double tol = REFERENCETOL)
  {return IntersectionBase<edgetype,sidetype,numNodesLine,numNodesSurface>::LineWithinLimits(tol); }

protected:

  void GetTriangle(LINALG::Matrix<3,3>& xyze_triElement, int triangleid)
  {return IntersectionBase<edgetype,sidetype,numNodesLine,numNodesSurface>::GetTriangle(xyze_triElement,triangleid); }

private:

  bool ComputeDistance( LINALG::Matrix<3,1> point, double& distance, double& tolerance,  bool& zeroarea, bool signeddistance = false, int triangleid = 0 );

  bool ComputeDistance( Point * p, double& distance, double& tolerance, bool& zeroarea, bool signeddistance = false, int triangleid = 0 );

  bool ComputeCut( Edge * e1, Edge * e2, double& tolerance );

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
    if ( AtEdge(xsi) )
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

  LINALG::Matrix<3, numNodesSurface>& xyze_surfaceElement_;
  LINALG::Matrix<3, numNodesLine>& xyze_lineElement_;

  LINALG::Matrix<3,1>& xsi_;
  LINALG::Matrix<3,1>& x_;

  double& scale_;
};


template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::Intersect( PointSet & cuts )
{
  //!WARNING: Intersection just works for planes (tri3,quad4 unwarped!) with lines!!!

  //the ordering of the steps in intersection is optimized for a general case, where special cases (paralell) are not happening
 //too often. That's why we first std-Intersection and then continue with treatment of special cases. --> so minimum effort for std-cases!
//--> not optimized for the cut tests!!!

  //intersection should also return a tolerance for the evaluation of each point, which will be used to merge points in the pointpool!!!

  //1// try to find not overlapping geometries with boundigboxes to avoid a big load of work ... has to be done
    //this is here just for performance ... intersection should also be robust without that!!!

  //2// first with starting the distance to both end points of a line, go get rid of paralell cases (where intersection wouldn't converge) and also get rid
    // of cases, where the line is just on one side of the surface --> definitly no intersection!!!
    //REMARK: as for quad4 where the projected end points of the line are outside the element, we do always get reliable results (normal can flip outside the element), generally
   //the distance is computed to the two triangles (be aware of the fact, that this is just possible because we limit this function to plane(unwarped) quad4 sides!!!

  //3//take care of special intersections on the surface edges

  //4// try to calculate the intersection point directly with the Newton
    //this will basically fail if the system is conditioned badly --> means that line and plane are parallel (which should be the case anymore as it was
    //already captured in point //2//) or element is distorted or it's a quad4 and
   //the intersection point is outside the element and is not part of the interpolation space!
   //these cases should be treated seperatly later!!!

  //5// take care of other special cases --> 5.1 quad4 outside interpolation space, 5.2 distorted elements

  //6// throw dserror in case this intersection wasn't treated right --> this means there is still handling of some special cases missing
     //in the code & it does not mean that there is no intersection point

  //1 bounding box check
  bool debugbb = false;
  BoundingBox SBB = GEO::CUT::BoundingBox(side_);
  BoundingBox EBB = GEO::CUT::BoundingBox(edge_);
  if (!SBB.Within(1.0,EBB))
  {
    #ifdef DEBUG_CUTINTERSECTION_OUTPUT
      std::cout << "No Cut Point, as not even the Bounding Boxes are Overlapping! (Codeline:" << __LINE__ << ")" << std::endl;
    #endif
    debugbb = true; //at the moment we want to test as many intersections as possible, that's why bounding boxes are not used yet!
    //return false; //not even the bounding boxes are overlapping, no need to search for an intersection point!
  }

  //2 line || plane

  bool success = false;
  Point * nodalpoint = NULL;

  //all information gathered by Compute Distance is gathered here!!! - 0...begin point / 1...end point of line!!!
  bool zeroarea;
  std::vector<bool> lineendpoint_within_surfacelimits(2,false); //local coordinates are inside element
  std::vector<bool> lineendpoint_in_surface(2,false); //point is really inside element (normal distance = 0)
  std::vector<double> lineendpoint_dist(2);
  std::vector<double> lineendpoint_tol(2);
  std::vector<bool> lineendpoint_conv(2,false);
  std::vector<LINALG::Matrix<3,1> > lineendpoint_xsi(2);

  //store the intersection xsi before it is overwritten by computlineendpoint_distance for later!
  LINALG::Matrix<3,1> i_xsi;

  for (uint lp = 0; lp < 2 ; ++lp) //loop over line endpoints!
  {
    GEO::CUT::Point * actpoint = 0;
    if (lp == 0)
      actpoint = edge_.BeginNode()->point();
    else if(lp == 1)
      actpoint = edge_.EndNode()->point();
    else
      dserror("You Expect more than two Endpoints for a line?");

    if (sidetype == DRT::Element::tri3)
    {
      lineendpoint_conv[lp] = ComputeDistance( actpoint, lineendpoint_dist[lp], lineendpoint_tol[lp],zeroarea , true ) ;
      lineendpoint_within_surfacelimits[lp] = SurfaceWithinLimits();
      lineendpoint_xsi[lp] = xsi_;
    }
    else if (sidetype == DRT::Element::quad4)
    {
        bool done = false;
      //problem here is that computedistance might change the normal direction if the projected point lies outside the element, that why prefer to compute the distance onto tri3!
        std::vector<double> tri_dist(2);
        std::vector<double> tri_tol(2);
        std::vector<bool> tri_conv(2);
        uint tri;
        for (tri = 0; tri < 2; tri++)
        {
          tri_conv[tri] = ComputeDistance( actpoint, tri_dist[tri], tri_tol[tri], zeroarea , true , tri+1 );

          #ifdef DEBUGCUTLIBRARY
            if (tri == 1)//just warping check!!!
            {
              if (!(fabs(tri_dist[0] - tri_dist[1]) < tri_tol[0] or fabs(tri_dist[0] - tri_dist[1]) < tri_tol[1]))
                dserror("Expect Distance from both Triangles to be the same --> is your quad4 warped??");
            }
          #endif

          if (fabs(tri_dist[tri]) >= tri_tol[tri]) //in this case the point is not on the surface an we stop calculation of distance right now! (could also be the second triangle if first one was distorted --> very big tolerance)
          {
            //REMARK: lineendpoint_within_surfacelimits[lp] is not correctly evaluated in this case, as it is not needed!!!
            done = true;
            break;
          }
          else
          {
            lineendpoint_within_surfacelimits[lp] = SurfaceWithinLimits(tri_tol[tri], DRT::Element::tri3);
            if (lineendpoint_within_surfacelimits[lp])
            {
              if (tri_conv[tri] == false)
                dserror("You found a case where the triangulation does not converge and it is inside the element, basically removing this dserror should be ok, but please contact me so that I can have a look into it!!! Christoph Ager");
              break;
            }
          }
        }
        if (done or tri == 2) //calculate distance to triangles was enough
        {
          if (tri == 2)
          {
            lineendpoint_within_surfacelimits[lp] = false;//in this case both tri3s where outside!
            tri = 1; //set back to the last tri
          }
          lineendpoint_conv[lp] = tri_conv[tri];
          lineendpoint_dist[lp] = tri_dist[tri];
          lineendpoint_tol[lp] = tri_tol[tri];
          #ifdef DEBUG_CUTINTERSECTION_OUTPUT
            lineendpoint_xsi[lp] = xsi_; //just set for output reasons!!!
          #endif
        }
        else //calculate distance to quad4, as the projected points lies inside and we want to calc local coord (should converge in 1 newton loop (theoretically also tri3 didn't converge, but then there is anyway no hope)
        {
         // std::cout << "calculate distance from quad4 ..." << std::endl;
          lineendpoint_conv[lp] = ComputeDistance( actpoint, lineendpoint_dist[lp], lineendpoint_tol[lp], zeroarea , true ) ;
          lineendpoint_within_surfacelimits[lp] = SurfaceWithinLimits();
          lineendpoint_xsi[lp] = xsi_;
        }
        if (zeroarea)
          std::runtime_error("zeroarea handling for quad4 not implemented yet (not expected to happen here, as it should just occure in the triangulation procedure!)");
      }
      else
        dserror("Intersection Error: Other surfaces than tri3 and quad4 are not supported to avoid huge problems!!!");
  }

  if (fabs(lineendpoint_dist[0]) >= lineendpoint_tol[0] and (fabs(lineendpoint_dist[1]) >= lineendpoint_tol[1]) and lineendpoint_dist[0]*lineendpoint_dist[1]>0 and !zeroarea) //both points outside the plane on the same side --> no intersection
  {
    #ifdef DEBUG_CUTINTERSECTION_OUTPUT
      std::cout << "No Cut Point, because we detected following normal distance at the line end points: ... lineendpoint_dist[0]: " << lineendpoint_dist[0] << " lineendpoint_dist[1]: " << lineendpoint_dist[1] << "(Codeline:" << __LINE__ << ")" << std::endl;
    #endif
    return false; //there is no intersection between surface and line!!!
  }

  if (lineendpoint_conv[0] && lineendpoint_conv[1]) //for all paralell cases this should converge!!!
  {

    if ( fabs(lineendpoint_dist[0]) < lineendpoint_tol[0] and lineendpoint_within_surfacelimits[0] )
      lineendpoint_in_surface[0] = true;

    if ( fabs(lineendpoint_dist[1]) < lineendpoint_tol[1] and lineendpoint_within_surfacelimits[1] )
      lineendpoint_in_surface[1] = true;

    if ( lineendpoint_in_surface[0] and lineendpoint_in_surface[1] ) //case where both nodes of the line are inside the plane
    {
      InsertCut( edge_.BeginNode(), cuts );
      InsertCut( edge_.EndNode()  , cuts );
      TestSideEdges( edge_.BeginNode()->point(), lineendpoint_xsi[0] );
      TestSideEdges( edge_.EndNode()  ->point(), lineendpoint_xsi[1] );
      #ifdef DEBUG_CUTINTERSECTION_OUTPUT
        std::cout << "Cut points found, Begin & End Node of the Line are inside the surface: ... lineendpoint_dist[0]: " << lineendpoint_dist[0] << " lineendpoint_dist[1]: " << lineendpoint_dist[1] << " (Codeline:" << __LINE__ << ")" << std::endl;
      #endif
      if (debugbb) dserror("Bounding Boxes said that ther is no intersection point!");
      return true;
    }
    else if (lineendpoint_in_surface[0])
    {
      InsertCut( edge_.BeginNode(), cuts );
      nodalpoint = edge_.BeginNode()->point();
      TestSideEdges( nodalpoint, lineendpoint_xsi[0]);
      #ifdef DEBUG_CUTINTERSECTION_OUTPUT
        std::cout << "Cut points found, Begin Node of the Line are inside the surface: ... lineendpoint_dist[0]: " << lineendpoint_dist[0] << " lineendpoint_dist[1]: " << lineendpoint_dist[1] << " (Codeline:" << __LINE__ << ")" << std::endl;
      #endif
      success = true;
    }
    else if (lineendpoint_in_surface[1])
    {
      InsertCut( edge_.EndNode(), cuts );
      nodalpoint = edge_.EndNode()->point();
      TestSideEdges( nodalpoint, lineendpoint_xsi[1]);
      #ifdef DEBUG_CUTINTERSECTION_OUTPUT
        std::cout << "Cut points found, End Node of the Line are inside the surface: ... lineendpoint_dist[0]: " << lineendpoint_dist[0] << " lineendpoint_dist[1]: " << lineendpoint_dist[1] << " (Codeline:" << __LINE__ << ")" << std::endl;
      #endif
      success = true;
    }
  }

  //3// see if edges are cut special!!!

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
        if (cut_points.size() == 1)
          cuts.insert(*(cut_points.begin()));
        else if (cut_points.size() == 2)
        {
          PointSet::iterator cp = cut_points.begin();
          cuts.insert(*cp);
          cp++;
          cuts.insert(*cp);
        }
        else
        {
          #ifdef DEBUG_CUTINTERSECTION_OUTPUT
            for (PointSet::iterator i = cut_points.begin(); i != cut_points.end(); ++i)
              (*i)->Print(std::cout);
            edge_.BeginNode()->point()->Print(std::cout);
            edge_.EndNode()->point()->Print(std::cout);
            e->BeginNode()->point()->Print(std::cout);
            e->EndNode()->point()->Print(std::cout);
          #endif

          dserror("Two Edges have more than two cutpoint, sounds strange!");
        }
        #ifdef DEBUG_CUTINTERSECTION_OUTPUT
          std::cout << "Cut points found, by intersection of edges!" << " (Codeline:" << __LINE__ << ")" << std::endl;
        #endif
        success = true;
      }
      else
      {
        double tolerance = 0.0;
        if ( ComputeCut( e, &edge_, tolerance ))
        {
          Point * p = Point::NewPoint( mesh_, x_.A(), xsi_( 2 ), &edge_, &side_, tolerance);
          p->AddEdge( e );
          cuts.insert( p );
          #ifdef DEBUG_CUTINTERSECTION_OUTPUT
            std::cout << "Cut points found, by intersection of edges!" << " (Codeline:" << __LINE__ << ")" << std::endl;
          #endif
          success = true;
        }
      }
    }

    if (success) //parallel cases handled!!!
    {if (debugbb) dserror("Bounding Boxes said that there is no intersection point!");
      return true;
    }
    else if (lineendpoint_conv[0] and lineendpoint_conv[1] )
    {
      if (fabs(lineendpoint_dist[0]) < lineendpoint_tol[0] and fabs(lineendpoint_dist[1]) < lineendpoint_tol[1]) // parallel to the plane but No Cut Points found!!!
      {
        #ifdef DEBUG_CUTINTERSECTION_OUTPUT
          std::cout << "No Cut Point, because we detected following normal distance at the line end points: ... lineendpoint_dist[0]: " << lineendpoint_dist[0] << " lineendpoint_dist[1]: " << lineendpoint_dist[1] << " and there was no edge intersection! (Codeline:" << __LINE__ << ")" << std::endl;
        #endif
        return false;
      }
    }
    else if (zeroarea)
    {
      #ifdef DEBUG_CUTINTERSECTION_OUTPUT
        std::cout << "No Cut Point, because surface has no area and edges don't have any intersection! (Codeline:" << __LINE__ << ")" << std::endl;
      #endif
      return false;
    }

//4//
  double itol;
  bool conv = ComputeCurveSurfaceIntersection(itol);
  if ( conv ) //in newton converges we trust the result!
  {
  //  std::cout << "Intersection converged! itol: " << itol << std::endl;
    if ( SurfaceWithinLimits() and LineWithinLimits() )
    {
      FinalPoint();
      if ( fabs( xsi_( 2, 0 )+1 ) < REFERENCETOL ) //case point lies on begin node of line!!!
      {
        InsertCut( edge_.BeginNode(), cuts );
        TestSideEdges( edge_.BeginNode()->point(), xsi_ );
        #ifdef DEBUG_CUTINTERSECTION_OUTPUT
          //Basically this case should already be cought in section //2// for paralell intersections!
          std::cout << "Cut points found by intersection, Begin Node is in surface! (Codeline:" << __LINE__ << ")" << std::endl;
        #endif
          if (debugbb) dserror("Bounding Boxes said that ther is no intersection point!");
        return true;
      }
      else if ( fabs( xsi_( 2, 0 )-1 ) < REFERENCETOL ) //case point lies on end node of line!!!
      {
        InsertCut( edge_.EndNode(), cuts );
        TestSideEdges( edge_.EndNode()->point(), xsi_ );
        #ifdef DEBUG_CUTINTERSECTION_OUTPUT
          //Basically this case should already be cought in section //2// for paralell intersections!
          std::cout << "Cut points found by intersection, End Node is in surface! (Codeline:" << __LINE__ << ")" << std::endl;
        #endif
        return true;
      }
      else
      {
        Node * n = side_.OnNode( x_ ); //point lies on a node of the side???
        if ( n!=NULL ) //no is a new point
        {
          InsertCut( n, cuts );
          TestSideEdges( n->point(), xsi_ );
          #ifdef DEBUG_CUTINTERSECTION_OUTPUT
            //Basically this case should already be cought in section //2// for paralell intersections!
            std::cout << "Cut points found by intersection, with local coords: xsi: " << xsi_(0,0) << " / eta: "<< xsi_(1,0) << " / alpha: "<< xsi_(2,0) << " (Codeline:" << __LINE__ << ")" << std::endl;
          #endif
            if (debugbb) dserror("Bounding Boxes said that there is no intersection point!");
          return true;
        }
        else
        {
          Point * p = Point::NewPoint( mesh_, x_.A(), xsi_( 2, 0 ), &edge_, &side_, itol );
          TestSideEdges( p, xsi_ );
          cuts.insert( p );
          #ifdef DEBUG_CUTINTERSECTION_OUTPUT
            std::cout << "Cut points found by intersection, with local coords: xsi: " << xsi_(0,0) << " / eta: "<< xsi_(1,0) << " / alpha: "<< xsi_(2,0) << " (Codeline:" << __LINE__ << ")" << std::endl;
          #endif
            if (debugbb) dserror("Bounding Boxes said that there is no intersection point!");
          return true;
        }
      }
    }
    else
    {
      #ifdef DEBUG_CUTINTERSECTION_OUTPUT
        std::cout << "No Cut Point found by intersection, with local coords: xsi: " << xsi_(0,0) << " / eta: "<< xsi_(1,0) << " / alpha: "<< xsi_(2,0) << " and a tolerance of " << itol << " (Codeline:" << __LINE__ << ")" << std::endl;
      #endif
      return false; //there is no intersection between surface and line!!!
    }
  }
  //5// treat speacial cases ...
  else
  {
    i_xsi = xsi_;

    //5.1 point is outside the element and is not part of the interpolation space (just for quad4 try triangulation)
    if (sidetype == DRT::Element::quad4)
    {
      std::vector<double> tri_tol(2);
      std::vector<bool> tri_conv(2);
      for (uint tri = 0; tri < 2; tri++)
      {
        tri_conv[tri] = ComputeCurveSurfaceIntersection(tri_tol[tri] , tri+1 );

        if (SurfaceWithinLimits(tri_tol[tri],DRT::Element::tri3)) //we expect the quad4 to converge in case that the projected point is inside the quad4 or tri3!
        {
          std::cout << "Local Coordinates from ComputeDistance (lineendpoint_dist[0] = " << lineendpoint_dist[0] << "[tol=" << lineendpoint_tol[0] << "]) on surface: " << "xsi: " << lineendpoint_xsi[0](0,0) << " / eta: "<< lineendpoint_xsi[0](1,0) << " / alpha: "<< lineendpoint_xsi[0](2,0) << std::endl;
          std::cout << "Local Coordinates from ComputeDistance (lineendpoint_dist[1] = " << lineendpoint_dist[1] << "[tol=" << lineendpoint_tol[1] << "]) on surface: " << "xsi: " << lineendpoint_xsi[1](0,0) << " / eta: "<< lineendpoint_xsi[1](1,0) << " / alpha: "<< lineendpoint_xsi[1](2,0) << std::endl;
          std::cout << "Local Coordinates from ComputeIntersection (tol=" << itol << ") on Quad4: " << "xsi: " << i_xsi(0,0) << " / eta: "<< i_xsi(1,0) << " / alpha: "<< i_xsi(2,0) << std::endl;
          std::cout << "Local Coordinates from ComputeIntersection (tol=" << tri_tol[tri] << ") on Tri3: " << "xsi: " << xsi_(0,0) << " / eta: "<< xsi_(1,0) << " / alpha: "<< xsi_(2,0) << std::endl;
          dserror("ComputeCurveSurfaceIntersection for Quad4 didn't converge, but ComputeCurveSurfaceIntersection for triangulation %d is inside the Element!",tri+1);
        }
      }
      if (tri_conv[0] and tri_conv[1]) //both converged and are outside!!! (tested directly)
      {
        #ifdef DEBUG_CUTINTERSECTION_OUTPUT
          std::cout << "No Cut Point found by intersection with triangulated quad4! (Codeline:" << __LINE__ << ")" << std::endl;
        #endif
        return false; //point outside the interpolation space of quad4
      }
      //else
        //Looks like this is the result of an distorted element, here triangulation is not better that the distance directly from the quad4
    }

    //5.2 distorted elements!!!
    //now let's do what is possible for the rest, basically we can just find out if there is no intersection, everything else would be pure speculation!

    {
      //scale this tolerance into reference system!!!

      //calculate the tolerance in the reference coordinate system!!!
      double scaling;
      double ritol;
      {
        scaling = xyze_surfaceElement_.NormInf(); //not definitely sure if this is a choise, but at least the same as in the cut kernel!
        double linescale = xyze_lineElement_.NormInf();
        if (linescale > scaling)
          scaling = linescale;
      }
      ritol = itol/scaling;

      xsi_ = i_xsi; //set xsi_ back to the intersection xsi!!!
      if ( SurfaceWithinLimits(ritol) and LineWithinLimits(ritol) )
      {
        //we think there shouldn't be an intersection point ... to prove this we refine our bounding boxes!!!
        if (sidetype == DRT::Element::tri3)
        {
          if (!RefinedBBOverlapCheck())
          {
            #ifdef DEBUG_CUTINTERSECTION_OUTPUT
              std::cout << "No Cut Point found by Refined Bounding Boxes! (Codeline:" << __LINE__ << ")" << std::endl;
            #endif
            return false;
          }
          else
          {
            //Output of the availabe information!!!
            std::cout << "local coordinates: <plane: " << xsi_(0) << ", " << xsi_(1) << ">, <line: " << xsi_(2) << "> | relative tolerance: "  << ritol << std::endl;
            std::cout << "Local Coordinates from ComputeDistance (lineendpoint_dist[0] = " << lineendpoint_dist[0] << "[tol=" << lineendpoint_tol[0] << "], conv = " << lineendpoint_conv[0] << ") on surface: " << "xsi: " << lineendpoint_xsi[0](0,0) << " / eta: "<< lineendpoint_xsi[0](1,0) << " / alpha: "<< lineendpoint_xsi[0](2,0) << std::endl;
            std::cout << "Local Coordinates from ComputeDistance (lineendpoint_dist[1] = " << lineendpoint_dist[1] << "[tol=" << lineendpoint_tol[1] << "], conv = " << lineendpoint_conv[1] << ") on surface: " << "xsi: " << lineendpoint_xsi[1](0,0) << " / eta: "<< lineendpoint_xsi[1](1,0) << " / alpha: "<< lineendpoint_xsi[1](2,0) << std::endl;

            {
              //just to creat gmsh output of the failed intersection!!!
              KERNEL::ComputeIntersection<edgetype, sidetype> ci( xsi_ );
              ci( xyze_surfaceElement_, xyze_lineElement_ );

              std::string filename(OUTPUT::GenerateGmshOutputFilename(".intersection_CUTFAIL.pos"));
              std::ofstream file(filename.c_str());
              ci.WritetoGmsh(file);
              file.close();
            }
            //6// throw dserror! - It is not allowed to remove this dserror (don't even do that in your local version, as the cut then just generates something undefined)!!!
            dserror("CRITICAL ERROR: Your Intersection seems to be a special case, which is not treated right yet!!! --> to fix this the intersection code in baci has to be improved!!!");
          }
        }
      }
      else
      {
        #ifdef DEBUG_CUTINTERSECTION_OUTPUT
          std::cout << "No Cut Point found by intersection, with local coords: xsi: " << xsi_(0,0) << " / eta: "<< xsi_(1,0) << " / alpha: "<< xsi_(2,0) << " and a tolerance of " << itol << " (Codeline:" << __LINE__ << ")" << std::endl;
        #endif
        return false;
      }
    }
  }

  //6// throw dserror! - It is not allowed to remove this dserror (don't even do that in your local version, as the cut then just generates something undefined)!!!
  dserror("CRITICAL ERROR: Your Intersection seems to be a special case, which is not treated right yet!!! --> to fix this the intersection code in baci has to be improved!!!");
  return false;
}

template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
//compute distance to side from point
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::ComputeDistance( Point * p, double& distance, double& tolerance,  bool& zeroarea, bool signeddistance, int triangleid )
{
  LINALG::Matrix<3,1> point( p->X() );

  return ComputeDistance(point, distance, tolerance, zeroarea, signeddistance, triangleid);
}

template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
//compute distance to side from point
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::ComputeDistance( LINALG::Matrix<3,1> point, double& distance, double& tolerance,  bool& zeroarea, bool signeddistance, int triangleid)
{
  TEUCHOS_FUNC_TIME_MONITOR( "ComputeDistance" );
  bool conv;

  if (!triangleid)
  {
    KERNEL::ComputeDistance<edgetype, sidetype> cd( xsi_ );
  //KERNEL::DebugComputlineendpoint_dist[1]ance<edgetype, sidetype> cd;
    conv = cd( xyze_surfaceElement_, point, distance,signeddistance );
    tolerance = cd.GetTolerance();
    zeroarea = cd.ZeroArea();
  }
  else
  {
    KERNEL::ComputeDistance<edgetype, DRT::Element::tri3> cd( xsi_ );

    LINALG::Matrix<3, 3> xyze_triElement;
    GetTriangle(xyze_triElement,triangleid);

    conv = cd( xyze_triElement, point, distance,signeddistance );
    tolerance = cd.GetTolerance();
    zeroarea = cd.ZeroArea();
  }

  //std::cout << "dist: " << distance << "conv: " << conv << "     tolerance: " << tolerance << std::endl;

return conv;
}

template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::ComputeCut( Edge * e1, Edge * e2, double& tolerance )
{
  return e2->ComputeCut( e1, xsi_( 2 ), x_, tolerance );
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

template <DRT::Element::DiscretizationType edgetype,
          DRT::Element::DiscretizationType sidetype,
          int numNodesLine,
          int numNodesSurface
>
bool Intersection<edgetype, sidetype, numNodesLine, numNodesSurface>::RefinedBBOverlapCheck(int maxstep)
{
  if (sidetype != DRT::Element::tri3)
    dserror("RefinedBBOverlapCheck is made for distored tri3s!");

  std::vector< LINALG::Matrix<3, 1> > surfpoints;
  std::vector< LINALG::Matrix<3, 1> > linepoints;
  LINALG::Matrix<3, 1> actpoint;
  for (int nid = 0; nid < numNodesSurface; nid ++)
  {
    (side_.Nodes()[nid])->Coordinates(actpoint.A());
    surfpoints.push_back(actpoint);
  }
  //find shortest edge on the surface!
  int smallestedge;
  {
    LINALG::Matrix<3, 1> d1, d2, d3;
    d1.Update(1.0, surfpoints[1],-1.0,surfpoints[2]);
    d2.Update(1.0, surfpoints[2],-1.0,surfpoints[0]);
    d3.Update(1.0, surfpoints[0],-1.0,surfpoints[1]);
    std::vector<double> lengths;
    lengths.push_back(d1.Norm2());
    lengths.push_back(d2.Norm2());
    lengths.push_back(d3.Norm2());
    if (lengths[0] < lengths[1] && lengths[0] < lengths[2]) //0 is smallest
      smallestedge = 0;
    else if(lengths[1] < lengths[2] && lengths[1] < lengths[0]) //1 is smallest
      smallestedge = 1;
    else //2 is smallest plus all the rest
      smallestedge = 2;
  }

  LINALG::Matrix<3, 1> v1(true);
  LINALG::Matrix<3, 1> v2(true);
  v1.Update(1.0, surfpoints[(smallestedge+1)%3],-1.0,surfpoints[smallestedge]);
  v2.Update(1.0, surfpoints[(smallestedge+2)%3],-1.0,surfpoints[smallestedge]);

  for (int nid = 0; nid < numNodesLine; nid ++)
  {
    (edge_.Nodes()[nid])->Coordinates(actpoint.A());
    linepoints.push_back(actpoint);
  }
  LINALG::Matrix<3, 1> v3(true);
  v3.Update(1.0, linepoints[1],-1.0,linepoints[0]);

  bool overlap = false; //we have an overlap of the refined bounding boxes?

  int pow_fac = 4;
  int max_steps = pow(pow_fac,maxstep-1);
  double dmax_steps = 1.0/max_steps;
  int act_steps = 0;
  uint act_boxidx = 0;
  std::vector< std::vector< int > > overlappingidx;
  std::vector< std::vector< int > > newoverlappingidx;
  std::vector< int > tmpoverlappingidx;
  tmpoverlappingidx.push_back(0); //min idx with overlap for tri3
  tmpoverlappingidx.push_back(max_steps); //max idx with overlap tri3
  tmpoverlappingidx.push_back(0); //min idx with overlap for line2
  tmpoverlappingidx.push_back(max_steps); //max idx with overlap for line2
  overlappingidx.push_back(tmpoverlappingidx);

  LINALG::Matrix<3, 1> p1,p2,p3,p4,lp1,lp2;

  for (int refinestep = 1; refinestep <= maxstep; refinestep++)
  {
    newoverlappingidx.clear();
    act_boxidx = 0;
    act_steps = pow(pow_fac,(maxstep - refinestep));

    overlap = false;
    for (int surfstep = 0; surfstep < max_steps; surfstep += act_steps)
    {
      tmpoverlappingidx.clear();
      while ((act_boxidx + 1) < overlappingidx.size() and (overlappingidx[act_boxidx + 1])[0] < surfstep) act_boxidx++;
      if (!(overlappingidx[act_boxidx][0] <= surfstep and overlappingidx[act_boxidx][1] >= surfstep + act_steps )) continue;
      double alpha = surfstep*dmax_steps;
      double alphap = (surfstep+act_steps)*dmax_steps;
      p1.Update(1.0, surfpoints[smallestedge], alpha,v1);
      p2.Update(1.0, surfpoints[smallestedge], alphap,v1);
      p3.Update(1.0, surfpoints[smallestedge], alpha,v2);
      p4.Update(1.0, surfpoints[smallestedge], alphap,v2);
      BoundingBox SBB;
      SBB.AddPoint(p1);
      SBB.AddPoint(p2);
      SBB.AddPoint(p3);
      SBB.AddPoint(p4);
      for (int linestep = 0; linestep < max_steps; linestep += act_steps)
      {
        if (!(overlappingidx[act_boxidx][2] <= linestep and overlappingidx[act_boxidx][3] >= linestep + act_steps )) continue;
        double lalpha = linestep*dmax_steps;
        double lalphap = (linestep+act_steps)*dmax_steps;
        lp1.Update(1.0, linepoints[0], lalpha,v3);
        lp2.Update(1.0, linepoints[0], lalphap,v3);
        BoundingBox EBB;
        EBB.AddPoint(lp1);
        EBB.AddPoint(lp2);

        if (SBB.Within(POSITIONTOL/BOXOVERLAP,EBB))
        {
          overlap = true;
          if (!tmpoverlappingidx.size())
          {
            tmpoverlappingidx.push_back(surfstep);//min idx with overlap for tri3
            tmpoverlappingidx.push_back(surfstep+act_steps);//max idx with overlap for tri3
            tmpoverlappingidx.push_back(linestep);//min idx with overlap for line2
            tmpoverlappingidx.push_back(linestep+act_steps);//max idx with overlap for line2
          }
          else
          {
            if (tmpoverlappingidx[2] > linestep) tmpoverlappingidx[2] = linestep;
            if (tmpoverlappingidx[3] < linestep + act_steps) tmpoverlappingidx[3] = linestep + act_steps;
          }
        }
      }
      if (tmpoverlappingidx.size()) newoverlappingidx.push_back(tmpoverlappingidx);
    }
    overlappingidx = newoverlappingidx;

    std::cout << "RefinedBBOverlapCheck: Refinement Level " << refinestep << " there is " << overlap << " overlap!" << std::endl;
    if (!overlap)
      break;
  }
  return overlap;
}

}
}

//static members of Intersection base class
template <DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine, int numNodesSurface > LINALG::Matrix< 3 , numNodesLine > GEO::CUT::IntersectionBase<edgetype,sidetype,numNodesLine,numNodesSurface>::xyze_lineElement_;
template <DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine, int numNodesSurface > LINALG::Matrix< 3 , numNodesSurface > GEO::CUT::IntersectionBase<edgetype,sidetype,numNodesLine,numNodesSurface>::xyze_surfaceElement_;
template <DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine, int numNodesSurface > LINALG::Matrix< 3 , 1 > GEO::CUT::IntersectionBase<edgetype,sidetype,numNodesLine,numNodesSurface>::xsi_(true);
template <DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine, int numNodesSurface > LINALG::Matrix< 3 , 1 > GEO::CUT::IntersectionBase<edgetype,sidetype,numNodesLine,numNodesSurface>::x_(true);

#endif

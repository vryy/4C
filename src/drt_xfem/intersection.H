/*!
\file intersection.H

\brief collection of intersection tools

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/

#ifdef CCADISCRET
#ifndef INTERSECTION_H
#define INTERSECTION_H


#include "../drt_xfem/intersection_interfacepoint.H"
#include "../drt_xfem/intersection_service.H"
#include "../drt_xfem/intersection_math.H"
#include "../drt_xfem/integrationcell.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_io/io.H"
#include "tetgen.h"
#include "Teuchos_RefCountPtr.hpp"


#ifdef PARALLEL
#include <mpi.h>
#endif


#if defined(__cplusplus)
extern "C"
{
#endif
#ifdef QHULL
#include <qhull/qhull.h>
#include <qhull/qset.h>
#endif
#if defined(__cplusplus)
}
#endif


namespace XFEM
{

enum boundaryType {NONE, ONFLUID, ONSOLID};
enum SteinerType {S_SURFACE, S_EDGE, S_BOUNDARY};


class Intersection
{

  public:

  /*!
  \brief  computes the intersection of two discretizations.

  This methods localizes the interface of the
  cutter discretization within the xfem discretization. Each intersected xfem element
  is tetrahedralized with help of a Constrained Delaunay Tetrahedralization.
  The output consists of a list of qadratic tetrahedral volume integration cells and
  a list of quadratic triangular boundary integration cells
  representing the curved interface.
  */
  void computeIntersection(
      const RCP<DRT::Discretization>         xfemdis,                 ///< xfem discretization
      const RCP<DRT::Discretization>         cutterdis,               ///< cutter discretization
      const std::map<int,BlitzVec3>&         currentcutterpositions,  ///< current positions of all cutter nodes
      map< int, DomainIntCells >&			       domainintcells,          ///< domainintegrationcells for each intersected element
      map< int, BoundaryIntCells >& 		     boundaryintcells         ///< boundaryintegrationcells for each intersected element
      );



  private:

  //! offset in facet marker list
  static const int                    facetMarkerOffset_ = 11;

  //! discretization type of the xfem element
  DRT::Element::DiscretizationType     xfemDistype_;

  //! discretization type of the cutter element
  DRT::Element::DiscretizationType     cutterDistype_;

  //! number of surfaces of the xfem element
  int                                   numXFEMSurfaces_;

  //! number of corner nodes of the xfem element
  int                                   numXFEMCornerNodes_;

  //! current nodal coordinates of the XFEM element
  BlitzMat                              xyze_xfemElement_;
    
  //! counts missed points during the recovery
  int                                   countMissedPoints_;

  //! vector of nodes for each surface of a single element
  vector<vector<int> >                  eleNumberingSurfaces_;
    
  //! vector of nodes for each line of a single element
  vector<vector<int> >                  eleNumberingLines_;

  //! vector of adjacent surfaces for each node of a single element
  vector<vector<int> >                  eleNodesSurfaces_;
  
  //! vector of adjacent surfaces for each node of a single element
  vector<vector<int> >                  eleNodesLines_;

  //! vector of adjacent surfaces for each line of a single element
  vector<vector<int> >                  eleLinesSurfaces_;

  //! nodal element coordinate
  vector<vector<double> >               eleRefCoordinates_;

  //! list of interface points to be stored in the Tetgen data structure
  vector<InterfacePoint >               pointList_;

  //! list of segments to be stored in the Tetgen data structure
  vector<vector<int> >                  segmentList_;

  //! list of single points on a xfem surface to be stored in the Tetgen data structure
  vector<vector<int> >                  isolatedPointList_;

  //! list of triangles to be stored in the Tetgen data structure
  vector<vector<int> >                  triangleList_;

  //! list of triangles lying on a xfem surface to be stored in the Tetgen data structure
  map< int, vector<vector<int> > >       surfaceTriangleList_;
  
  //! list of intersecting cutter elements
  vector<DRT::Element* >                intersectingCutterElements_;

  //! list of face markers numbered in the order of the intersected surface elments
  vector<int>                           faceMarker_;
  
  /*!
  \brief initializes the private members of the xfem element
  
  \param xfemId        (in) : id of the current xfem element
  \param xfemElement   (in) : xfem element
  */
  void initializeXFEM(
    const int             xfemId,
    const DRT::Element*   xfemElement
    );

        
  /*!
  \brief collects points belonging to an interface lying inside xfem elements
         The nonlinear system of equation is solved with help of the Newton-method.

  \param xfemElement              (in)        : xfem element
  \param cutterElement            (in)        : cutter element
  \param xyze_cutterElement       (in)        : nodal coordinates of the cutter element
  \param cutterNode               (in)        : cutter node
  \param currentcutterpositions   (in)        : current nodal positions of the cutter element
  \param interfacePoints          (in/out)    : vector of interface points
  \param numTotalInterfacePoints  (in)        : total number of interface points
  \param elemId                   (in)        : xfem element id  (Debug)
  \param nodeId                   (in)        : node id
  \return true if the node lies within the xfemElement, false otherwise
  */
  bool collectInternalPoints(
    const DRT::Element*             xfemElement,             
    DRT::Element*                   cutterElement,           
    const BlitzMat&                 xyze_cutterElement,      
    const DRT::Node*                cutterNode,              
    const map<int,BlitzVec3>&       currentcutterpositions,  
    std::vector<InterfacePoint>&    interfacePoints,
    int&                            numTotalInterfacePoints,
    const int                       elemId,
    const int                       nodeId
    );


  /*!
  \brief  sets the boundary status of a boundary point
          (xfem surface, line or node)          .

  \param xfemDistype          (in)        : discretization type of xfem element
  \param xsi                  (in)        : node in element coordinates
  \param ip                   (in)        : interface point
  */
  void setBoundaryPointBoundaryStatus(
    const DRT::Element::DiscretizationType  xfemDistype,
    const BlitzVec3&                        xsi,
    InterfacePoint&                         ip
    ) const;
  
 
  /*!
  \brief  checks if an internal point
          lies on one a xfem surface, line or node or internally          .

  \param xfemDistype          (in)        : discretization type of xfem element
  \param xsi                  (in)        : node in element coordinates
  \param ip                   (in)        : interface point
  \return true if the node is lying on a boundary false otherwise
  */
  bool setInternalPointBoundaryStatus(
    const DRT::Element::DiscretizationType  xfemDistype, 
    const BlitzVec3&                        xsi,
    InterfacePoint&                         ip
    ) const;
        
        
  /*!
  \brief  checks if an intersection node
          lies on a xfem surface, line or node          .

  \param xfemElement          (in)        : xfem element
  \param surfaceElement       (in)        : xfem surface element
  \param xyze_surfaceElement  (in)        : nodal coordinates of the surface element
  \param xsi                  (in)        : node in element coordinates
  \param ip                   (in)        : interface point
  \return true if the node is lying on a boundary false otherwise
  */
  void setIntersectionPointBoundaryStatus(
    const DRT::Element*                     xfemElement,
    const DRT::Element*                     surfaceElement,
    const BlitzMat&                         xyze_surfaceElement,
    const BlitzVec3&                        xsiSurface,
    InterfacePoint&                         ip
    ) const;


  /*!
  \brief collects all intersection points between a line and a surface

  \param xfemElement              (in)        : xfem element
  \param surfaceElement           (in)        : surface element
  \param xyze_surfaceElement      (in)        : nodal coordinates of surface element
  \param lineElement              (in)        : line element
  \param interfacePointList       (out)       : nodal coordinates of line element
  \param numTotalInterfacePoints  (in)        : total number of interface points per xfem element
  \param surfaceId                (in)        : surface element id
  \param lineId                   (in)        : line element id
  \param lines                    (in)        : if lines = true intersection of cutter line with xfem element
  \param doSVD                    (in)        : compute SVD if not Cartesian surface element and linear line elments present
  */
  bool collectIntersectionPoints(
    const DRT::Element*             xfemElement,
    const DRT::Element*             surfaceElement,
    const BlitzMat&                 xyze_surfaceElement,
    const DRT::Element*             lineElement,
    const BlitzMat&                 xyze_lineElement,
    std::vector<InterfacePoint>&    interfacePoints,
    int&                            numTotalInterfacePoints,
    const int                       surfaceId,
    const int                       lineId,
    const bool                      lines,
    const bool                      doSVD
    ) const;

  
  /*!
  \brief checks if surface is Cartesian and line is linear

  \param surfaceGeoType   (in)        : geotype of surface element 
  \param lineGeoType      (in)        : geotype of line element 
  \return true if surface is Cartesian an line is linear
  */
  bool decideSVD(
      const EleGeoType surfaceGeoType,
      const EleGeoType lineGeoType
      );
       
  
  /*!
  \brief checks if a linear line lies in the plane of a Cartesian surface element

  \param surfaceElement           (in)        : surface element
  \param xyze_surfaceElement      (in)        : nodal coordinates of surface element
  \param lineElement              (in)        : line element
  */
  bool checkIfLineInSurface(
    const DRT::Element*             surfaceElement,
    const BlitzMat&                 xyze_surfaceElement,
    const BlitzMat&                 xyze_lineElement
    )const;

  
  /*!
  \brief computes a new starting points for the Newton method recursively
         in order to find all intersection points

  \param xfemElement              (in)    : xfem element
  \param surfaceElement           (in)    : surface element
  \param xyze_surfaceElement      (in)    : nodal coordinates of surface element
  \param lineElement              (in)    : line element
  \param xyze_lineElement         (in)    : nodal coordinates of line element
  \param surfaceId                (in)    : surface element id
  \param lineId                   (in)    : line element id
  \param upLimit                  (in)    : old stating vector
  \param upLimit                  (in)    : upper limit of the search interval
  \param loLimit                  (in)    : lower limit of the search interval
  \param interfacePointList       (out)   : vector of interface points
  \param lines                    (in)    : if lines = true
  \param doSVD                    (in)    : compute SVD if not Cartesian surface element and linear line elments present
  \return number of interface points
  */
  int computeNewStartingPoint(
    const DRT::Element*             xfemElement,
    const DRT::Element*             surfaceElement,
    const BlitzMat&                 xyze_surfaceElement,
    const DRT::Element*             lineElement,
    const BlitzMat&                 xyze_lineElement,
    const int                       surfaceId,
    const int                       lineId,
    const BlitzVec3&                xsiOld,
    const BlitzVec3&                upLimit,
    const BlitzVec3&                loLimit,
    std::vector<InterfacePoint>&    interfacePoints,
    const bool                      lines,
    const bool                      doSVD
    ) const;


  /*!
  \brief computes an interseticon point between a curve and a surface

      The nonlinear system of equation is solved with help of the Newton-method.

  \param surfaceElement           (in)    : surface element
  \param xyze_surfaceElement      (in)    : nodal coordinates of surface element
  \param lineElement              (in)    : line element
  \param xyze_lineElement         (in)    : nodal coordinates of line element
  \param xsi                      (in/out): starting value/vector of element coordinates
  \param upLimit                  (in)    : upper search interval boundary
  \param loLimit                  (in)    : lower search interval boundary
  \param doSVD                    (in)    : compute SVD if not Cartesian surface element and linear line elments present
  \param tol                      (in)    : tolerance
  return true if an intersection point was found, otherwise false
  */
  bool computeCurveSurfaceIntersection(
    const DRT::Element*                 surfaceElement,
    const BlitzMat&                     xyze_surfaceElement,
    const DRT::Element*                 lineElement,
    const BlitzMat&                     xyze_lineElement,
    const BlitzVec3&                    upLimit,
    const BlitzVec3&                    loLimit,
    BlitzVec3&                          xsi,
    const bool                          doSVD,
    const double                        tol
    ) const;


	/*!
	\brief adds an intersection point to the list of interface points

  \param surfaceElement           (in)    : surface element
  \param xyze_surfaceElement      (in)    : nodal coordinates of surface element
  \param lineElement              (in)    : line element
  \param xyze_lineElement         (in)    : nodal coordinates of line element
	\param xsi           			      (in)    : surface element
	\param upLimit                  (in)	  : upper limit of the search interval
	\param loLimit                  (in)	  : lower limit of the search interval
	\param interfacePointList       (out)   : vector of interface points
	\param surfaceId                (in)   	: surface element id
  \param lineId                   (in)   	: line element id
  \param lines                    (in)    : if lines = true
  \param doSVD                    (in)    : compute SVD if not Cartesian surface element and linear line elments present
	\return number of interface points
	*/
  int addIntersectionPoint(
    const DRT::Element*                 xfemElement,
    const DRT::Element*                 surfaceElement,
    const BlitzMat&                     xyze_surfaceElement,
    const DRT::Element*                 lineElement,
    const BlitzMat&                     xyze_lineElement,
    const BlitzVec3&                    xsi,
    const BlitzVec3&                    upLimit,
    const BlitzVec3&                    loLimit,
    std::vector<InterfacePoint>& 	      interfacePoints,
    const int                      	    surfaceId,
    const int                           lineId,
    const bool                          lines,
    const bool                          doSVD
    ) const;


  /*!
  \brief  creates new limits for the recursive intersection computation

  \param xsi              (in)        : intersection point
  \param upLimit          (in)        : upper limits
  \param loLimit          (in)        : lower limits
  \param upperLimits      (out)       : list of upper limits of new ranges
  \param lowerLimits      (out)       : list of lower limits of new ranges
  */
  void createNewLimits(
    const BlitzVec3&        xsi,
    const BlitzVec3&        upLimit,
    const BlitzVec3&        loLimit,
    vector<BlitzVec3>&      upperLimits,
    vector<BlitzVec3>&      lowerLimits
    ) const;

        
  /*!
  \brief  finds th id of the an xfem surface common to all interafce points.
          Returns -1 , if no common surface exists or the surface id

  \param xfemElement              (in)    : xfem element
  \param xyze_xfemElement         (in)    : nodal coordinates of xfem element
  \param cutterElement            (in)    : cutter element
  \param xyze_cutterElement       (in)    : nodal coordinates of cutter element
  \param postions                 (in)    : vector of interface point positions in pointList_
  \return                                 returns -1 or the surface id
  */
  int findCommonSurfaceID(
    const DRT::Element*     xfemElement,
    const BlitzMat&         xyze_xfemElement,
    const DRT::Element*     cutterElement,
    const BlitzMat&         xyze_cutterElement,
    const vector<int>&      positions
    );
      
  
  /*!
  \brief  finds th id of the an xfem surface common to all interafce points.
          Returns -1 , if no common surface exists or the surface id

  \param xfemElement              (in)    : xfem element
  \param xyze_xfemElement         (in)    : nodal coordinates of xfem element
  \param cutterElement            (in)    : cutter element
  \param xyze_cutterElement       (in)    : nodal coordinates of cutter element
  \param postions                 (in)    : vector of interface point positions in pointList_
  \return                                 returns true if on surface
  */
  bool checkIfCutterOnXFEMSurface(
    const DRT::Element*             xfemElement,
    const BlitzMat&                 xyze_xfemElement,
    const DRT::Element*             cutterElement,
    const BlitzMat&                 xyze_cutterElement,
    const vector<int>&              positions
    ); 
  
  
  /*!
  \brief  checks if an intersection lies on a xfem boundary

  \param xfemElement              (in)    : xfem element
  \param xyze_xfemElement         (in)    : nodal coordinates of xfem element
  \param cutterElement            (in)    : cutter element
  \param xyze_cutterElement       (in)    : nodal coordinates of cutter element
  \param interfacePoints          (in)    : vector of interface points
  \param surfId                   (in)    : Id of common surface
  */
  XFEM::boundaryType checkIfIntersectionOnXFEMBoundary(
    const DRT::Element*             xfemElement,
    const BlitzMat&                 xyze_xfemElement,
    const DRT::Element*             cutterElement,
    const BlitzMat&                 xyze_cutterElement,
    const vector<InterfacePoint>&   interfacePoints,
    const int                       surfId
    );

        
  /*!
  \brief  prepares a part of the piecewise linear complex for 
          a xfem and a cutter element

  \param xfemElement              (in)    : xfem element
  \param xyze_xfemElement         (in)    : nodal coordinates of xfem element
  \param cutterElement            (in)    : cutter element
  \param xyze_cutterElement       (in)    : nodal coordinates of cutter element
  \param interfacePoints          (in)    : list of interface points
  \param boundaryIntersectionType (in)    : boundary intersection type
  */
  void preparePLC(
    const DRT::Element*         xfemElement,
    const BlitzMat&             xyze_xfemElement,
    const DRT::Element*         cutterElement,
    const BlitzMat&             xyze_cutterElement,
    vector<InterfacePoint>&     interfacePoints,
    XFEM::boundaryType&         boundaryIntersectionType 
    );

        
  /*!
  \brief  computes the convex hull of a set of points
          points ordered according to the convex hull are stored in vertices

  \param xfemElement              (in)      : xfem element
  \param cutterElement            (in)      : cutter element
  \param xyze_cutterElement       (in)      : nodal coordinates of cutter element
  \param interfacePoints          (in)      : list of interface points
  \param vertices                 (in/out)  : ordered interface points
  \param midpoint                 (in/out)  : midpoint
  */
  void computeConvexHull(
    const DRT::Element*           xfemElement,
    const DRT::Element*           cutterElement,
    const BlitzMat&               xyze_cutterElement,
    vector<InterfacePoint>&       interfacePoints,
    vector< vector<double> >&     vertices,
    InterfacePoint&               midpoint
    );   
        
        
  /*!
  \brief  stores a part of the piecewise linear complex for 
          a xfem and a cutter element

  \param xfemElement              (in)      : xfem element
  \param cutterElement            (in)      : cutter element
  \param xyze_cutterElement       (in)      : nodal coordinates of cutter element
  \param surfId                   (in)      : id of xfem element surface
  \param vertices                 (in)      : position of points
  \param midpoint                 (in)      : midpoint
  
  */
  void storePLC(
    const DRT::Element*             xfemElement,
    const DRT::Element*             cutterElement,
    const BlitzMat&                 xyze_cutterElement,
    const int                       surfId,
    vector<int>&                    positions,
    InterfacePoint&                 midpoint
    ); 
       
  
  /*!
  \brief  completes the PLC

  */
  void completePLC(
      );
        
        
  /*!
  \brief finds the next segment of a convex hull in clockwise order

  \param vertices         (in)        : vector of facet vertices
  \param searchPoint      (in/out)    : common point of the previous and next
                                        facet/the new point of the next facet
  */
  void findNextSegment(
    vector< vector<double> >&   vertices,
    vector<double>&             searchPoint
    ) const;

  
  /*!
  \brief checks if a CDT has to be computed for the current xfem element

  \param numTotalInterfacePoints  (in)        : total number of interface points for the current xfem element
  \param boundaryIntersectionType (in)        : boundary intersection type
  */
  bool checkIfCDT(
    const int             numTotalInterfacePoints,
    const boundaryType    boundaryIntersectionType
    );
    
  
  /*!
  \brief  computes the Constrained Delaunay Tetrahedralization CDT
          of an intersected xfem element

  \param element                  (in)    : xfem element
  \param xyze_xfemElement         (in)    : nodal coordinates of xfem element
  \param currentcutterpositions   (in)    : current cutter node positions
  \param domainintcells           (out)   : list of domainintegrationcells for each intersected element
  \param boundaryintcells         (out)   : list of boundaryintegrationcells for each intersected element
  \param timestepcounter_         (out)   : time step counter
  */
  void computeCDT(
    const DRT::Element*                     xfemElement,
    const BlitzMat&                         xyze_xfemElement,
    const map<int,BlitzVec3>&               currentcutterpositions,  
    map< int, DomainIntCells >&             domainintcells,
    map< int, BoundaryIntCells >&           boundaryintcells,
    int                                     timestepcounter_
    );


  /*!
  \brief  adds the corner nodes of an xfem element to the pointlist

  */
  void startPointList();


  /*!
  \brief  stores all interface point of a curve-surface intersection 
          in a point list for the
          Constrained Delaunay Tetrahedralization (CDT) with Tetgen

  \param vertices         (in)        :   coordinates of the point to be stored
  \param interfacePoints  (in/out)    :   positions vector of interface points
  \param positions        (in/out)    :   vector of interface points
  */
  void storePointList(
    vector< vector<double> >&       vertices,
    vector<int>&                    positions,
    vector<InterfacePoint>&         interfacePoints
    );
  
  
  /*!
  \brief  stores an interface point in a point list for the
          Constrained Delaunay Tetrahedralization (CDT) with Tetgen

  \param point            (in)        :   coordinates of the point to be stored
  \param interfacePoints  (in/out)    :   vector of interface points
  \param positions        (in/out)    :   positions
  */
  void storePoint(
    const vector<double>&             point,
    vector<InterfacePoint>&           interfacePoints,
    vector<int>&                      positions
    );

  
  /*!
  \brief  stores an interface point in a point list for the
          Constrained Delaunay Tetrahedralization (CDT) with Tetgen

  \param midPoint            (in)     :   coordinates of the point to be stored
  \param positions        (in/out)    :   positions
  */ 
  void storeMidPoint(
    const InterfacePoint&             midPoint,
    vector<int>&                      positions
    );
  
  
  /*!
  \brief  rounds interface point on xfem boundary
          use TOL3 to make TetGen happy

  \param ip            (in/out)     :   interface point
  \param tol           (in)         :   tolerance
  */ 
  void roundPointOnXFEMBoundary(
        InterfacePoint&                   ip,
        const double                      tol
        );
  
  
  /*!
  \brief  stores a single point in a list of surface points if no
          segment lies in that surface

  \param interfacePoints      (in)        :   list of interface points

  */
  void storeSurfacePoints(
    const vector<InterfacePoint>&     interfacePoints
    );
        
  
  /*!
  \brief  stores a single point in a list of surface points if no
          segment lies in that surface

  \param positions            (in)        :   list of interface point positions

  */
  void storeIsolatedPoints(
    const vector<int>&                positions
    );
    
  
  /*!
  \brief  computes the midpoint of a collection of points

  \param interfacePoints            (in)        : vector of interface points
  \return returns the midpoint (type InterfacePoint)
  */
  InterfacePoint computeMidpoint(
    const vector<InterfacePoint>&     interfacePoints
    ) const;
        
  
  /*!
  \brief  computes the midpoint of a collection of points determined by 
          a pointlist position vector

  \param positions           (in)        : vector of positions
  \return returns the midpoint (type InterfacePoint)
  */
  XFEM::InterfacePoint computeMidpoint(
      const vector<int>&  positions
      ) const;
     
  
  /*!
  \brief  classifies the midpoint 

  \param surfId                     (in)        : surface Id of common surface
  \param interfacePoints            (in)        : vector of interface points
  \return returns the midpoint (type InterfacePoint)
  */
  void classifyMidpoint(
      const int         surfId,
      InterfacePoint&   midpoint
      ) const;
        
        
  /*!
  \brief  stores a segment in a segment list for the
          Constrained Delaunay Tetrahedralization (CDT) with Tetgen

  \param positions            (in)        :   positions
  */
  void storeSegments(
    const vector<int>&              positions
    );
        
  
  /*!
  \brief  store all segments for cutter elements lying with more than two
          intersection points on an xfem surface

  \param positions            (in)        :   positions
  */
  void storeSurfaceSegments(
    const vector<int>&              positions
    );

        
  /*!
  \brief  store all segments for cutter elements lying with more than two
          intersection points on an xfem surface

  \param pos1            (in)        :   position of 1. point in pointList_
  \param pos2            (in)        :   position of 2. point in pointList_
  */
  void storeSingleSegment(
    const int pos1,
    const int pos2
    );
        
  
  /*!
  \brief checks if two segment points are lying on the same xfemline

  \param position1            (in)        :   position of 1. point in pointList_
  \param position2            (in)        :   position of 2. point in pointList_
  */
  bool checkIfSegmentPointsOnSameXfemLine(const int position1,
    const int position2);
        
        
  /*!
  \brief checks if triangles have two identical points these triangles are removed from point list
  \param positions           (in)        :   point positions
  */
  bool checkIfTrianglesDegenerate(
    vector<int>&              positions
    );
      
  
  /*!
  \brief  stores a triangle facet in a triangle list for the
          Constrained Delaunay Tetrahedralization (CDT) with Tetgen

  \param positions            (in)        :   positions
  */
  void storeTriangles(
    const vector<int>               positions
    );
        
        
  /*!
  \brief  stores a triangle facet in a triangle list for the
          Constrained Delaunay Tetrahedralization (CDT) with Tetgen

  \param surfID               (in)        :   id of xfem surface
  \param positions            (in)        :   positions
  */
  void storeSurfaceTriangles(
    const int                       surfId,
    const vector<int>               positions
    );


  /*!
  \brief  stores the cutter element intersecting
          with an xfem element in a list

  \param surfaceElement   (in)    :   surface element
  */
  void storeIntersectedCutterElement(
    DRT::Element* cutterElement
    );

  
  /*!
  \brief  computes the recovery of the curved interface
          for one xfem element after a Contrained Delaunay
          Tetrahedralization

  \param xfemElement              (in)    : xfem element
  \param xyze_xfemElement         (in)    : nodal coordinates of xfem element
  \param currentcutterpositions   (in)    : current cutter node positions
  \param domainintcells           (out)   : list of domainintegrationcells for each intersected element
  \param out                      (out)   : output data structure of TetGen
  */
  void recoverCurvedInterface(
    const DRT::Element*             xfemElement,
    const BlitzMat&                 xyze_xfemElement,
    const map<int,BlitzVec3>&       currentcutterpositions,  
    map< int, BoundaryIntCells >&   boundaryintcells,
    tetgenio&                       out,
    bool 							              recovery
    );

  
 	/*!
  \brief  stores linear boundary cells

  \param xfemElement              (in)    :   xfem element
  \param xyze_xfemElement         (in)    :   nodal coordinates of xfem element
  \param currentcutterpositions   (in)    :   current cutter node positions
  \param domainintcells           (out)   :   list of domainintegrationcells for each intersected element
  \param out                      (out)   :   output data structure of TetGen
  */
  void storeIntCells(
    const DRT::Element*             xfemElement,
    const BlitzMat&                 xyze_xfemElement,
    const map<int,BlitzVec3>&       currentcutterpositions,  ///< current positions of all cutter nodes
    map< int, BoundaryIntCells >&   boundaryintcells,
    tetgenio&                       out
    );

  
  /*!
  \brief  stores linear boundary cells

  \param xfemElement              (in)    :   xfem element
  \param xyze_xfemElement         (in)    :   nodal coordinates of xfem element
  \param currentcutterpositions   (in)    :   current cutter node positions
  \param domainintcells           (out)   :   list of domainintegrationcells for each intersected element
  \param out                      (out)   :   output data structure of TetGen
  */
	void storeSurfaceIntCells(
    const bool                      higherorder,
    const DRT::Element*             xfemElement,
    const BlitzMat&                 xyze_xfemElement,
    const map<int,BlitzVec3>&       currentcutterpositions,
    BoundaryIntCells&               listBoundaryICPerElement
    );
  	
    	
  /*!
  \brief  checks if Steiner points are inserted additionally. They are lifted up such they
          are lying in a cutter element.

  \param xfemElement              (in)    :   xfem element
  \param xyze_xfemElement         (in)    :   nodal coordinates of xfem element
  \param currentcutterpositions   (in)    :   current cutter node positions
  \param out                      (out)   :   output data structure of TetGen
  */
  void liftAllSteinerPoints(
    const DRT::Element*                   xfemElement,
    const BlitzMat&                       xyze_xfemElement,
    const map<int,BlitzVec3>&             currentcutterpositions,  
    tetgenio&                             out
    );


  /*!
  \brief  collects for each Steiner point its adjacent triangular faces

  \param adjacentFacesList        (out)    :   list of adjacent faces
  \param adjacentFacemarkerList   (out)    :   list of face marker of each adjacent face
  \param out                      (in)     :   output data structure of TetGen
  */
  void locateSteinerPoints(
    vector< vector<int> >&      adjacentFacesList,
    vector< vector<int> >&      adjacentFacemarkerList,
    const tetgenio&             out
    ) const;


  /*!
  \brief  checks if the Steiner points will lie in a cutter element or on
          one of the edges of the cutter element

  \param steinerIndex             (in)     :   index of the Steiner point
  \param lineIndex                (out)    :   line index
  \param cutterIndex              (out)    :   index of the cutter element
  \param adjacentFacesList        (in)     :   list of adjacent faces
  \param adjacentFacemarkerList   (in)     :   list of face marker of each adjacent face
  \param currentcutterpositions   (in)     :   current cutter node positions
  \param edgePoint                (out)    :   point lying on the same edge as the Steiner point
  \param oppositePoint            (out)    :   point lying opposite on the edge of the Steiner point
  \param xfemElement              (in)     :   xfem element
  \param xyze_xfemElement         (in)     :   nodal coordinates of xfem element
  \param out                      (in)     :   output data structure of TetGen
  \return steiner type S_SURAFCE S_EDGE S_BOUNDARY
  */
  XFEM::SteinerType decideSteinerCase(
    const int                         steinerIndex,
    int&                              lineIndex,
    int&                              cutterIndex,
    const vector< vector<int> >&      adjacentFacesList,
    const vector< vector<int> >&      adjacentFacemarkerList,
    const map<int,BlitzVec3>&         currentcutterpositions,  
    BlitzVec&                         edgePoint,
    BlitzVec&                         oppositePoint,
    const DRT::Element*               xfemElement,
    const BlitzMat&                   xyze_xfemElement,
    const tetgenio&                   out
    ) const;


  /*!
  \brief  lifts the Steiner point into the cutter element

  \param steinerIndex             (in)     :   index of the Steiner point
  \param adjacentFacesList        (in)     :   list of adjacent faces
  \param adjacentFacemarkerList   (in)     :   list of face marker of each adjacent face
  \param currentcutterpositions   (in)     :   current cutter node positions
  \param xfemElement              (in)     :   xfem element
  \param xyze_xfemElement         (in)     :   nodal coordinates of xfem element
  \param out                      (in)     :   output data structure of TetGen
  */
  void liftSteinerPointOnSurface(
    const int                     steinerIndex,
    const vector<vector<int> >&   adjacentFacesList,
    const vector<vector<int> >&   adjacentFacemarkerList,
    const map<int,BlitzVec3>&     currentcutterpositions,  
    const DRT::Element*           xfemElement,
    const BlitzMat&               xyze_xfemElement,
    tetgenio&                     out
    );


  /*!
  \brief  lifts the Steiner point onto the edge of the cutter element
          specified by the line index

  \param steinerIndex             (in)     :   index of the Steiner point
  \param lineIndex                (in)     :   line index ??????
  \param cutterIndex              (in)     :   index of the cutter element
  \param edgePoint                (in)     :   point lying on the same edge as the Steiner point
  \param oppositePoint            (in)     :   point lying opposite on the edge of the Steiner point
  \param adjacentFacesList        (in)     :   list of adjacent faces
  \param currentcutterpositions   (in)     :   current cutter node positions
  \param xfemElement              (in)     :   xfem element
  \param xyze_xfemElement         (in)     :   nodal coordinates of xfem element
  \param out                      (in/out) :   output data structure of TetGen
  */
  void liftSteinerPointOnEdge(
    const int                         steinerIndex,
    int                               lineIndex, // TODO: ???????  why not const, if input ????????
    const int                         cutterIndex,
    BlitzVec&                         edgePoint,
    BlitzVec&                         oppositePoint,
    const vector<vector<int> >&       adjacentFacesList,
    const map<int,BlitzVec3>&         currentcutterpositions,  
    const DRT::Element*               xfemElement,
    const BlitzMat&                   xyze_xfemElement,
    tetgenio&                         out
    );


  /*!
  \brief  lifts the Steiner point on an edge that lies on the xfem boundary

  \param steinerIndex             (in)        :   index of the Steiner point
  \param adjacentFacesList        (in)        :   list of adjacent faces
  \param adjacentFacemarkerList   (in)        :   list of adjacent face marker
  \param currentcutterpositions   (in)        :   current cutter node positions
  \param xfemElement              (in)        :   xfem element
  \param xyze_xfemElement         (in)        :   nodal coordinates of xfem element
  \param out                      (in/out)    :   output data structure of TetGen
  */
  void liftSteinerPointOnBoundary(
    const int                         steinerIndex,
    const vector<vector<int> >&       adjacentFacesList,
    const vector<vector<int> >&       adjacentFacemarkerList,
    const map<int,BlitzVec3>&         currentcutterpositions,  
    const DRT::Element*               xfemElement,
    const BlitzMat&                   xyze_xfemElement,
    tetgenio&                         out
    );


  /*!
  \brief  gathers information of a tetrahedron adjacent to an
          interface

  \param tetIndex             (in)    : tetrahedron index of the tetrahedron list
  \param faceIndex            (in)    : index of the face
  \param tetraCornerIndices   (out)   : index list of tetrahedron corners
  \param order                (out)   : order the tetrahedron corner points appear
  \param out                  (in)    : output data structure of TetGen
  */
  void getTetrahedronInformation(
    const int           tetIndex,
    const int           faceIndex,
    vector<int>&        tetraCornerIndices,
    vector<int>&        order,
    const tetgenio&     out
    ) const;


  /*!
  \brief  gathers information of a tetrahedron adjacent to an
          interface

  \param tetraCornerIndices   (in)    : index list of tetrahedron corners
  \param tetraCornerIndices   (out)   : list of tetrahedron corners nodes in the current
                                        configuration of the xfem element
  \param xfemElement          (in)    : xfem element
  \param xyze_xfemElement     (in)        :   nodal coordinates of xfem element
  \param out                  (in)    : output data structure of TetGen
  */
  void getTetrahedronNodes(
    vector<BlitzVec>&       tetraCornerNodes,
    const vector<int>&      tetraCornerIndices,
    const DRT::Element*     xfemElement,
    const BlitzMat&         xyze_xfemElement,
    const tetgenio&         out
    ) const;


   /*!
   \brief  lifts a higher order point of an edge of the linearized interface
   onto the curved interface

   \param index1                   (in)    :   local index of edge point 1
   \param index2                   (in)    :   local index of edge point 2
   \param faceIndex                (in)    :   index of the face
   \param faceMarkerIndex          (in)    :   face marker
   \param globalHigherOrderIndex   (in)    :   global index of the higher-order point
   \param xfemElement              (in)    :   xfem element
   \param xyze_xfemElement         (in)    :   nodal coordinates of xfem element
   \param currentcutterpositions   (in)    :   current cutter node positions
   \param out                      (in)    :   output data structure of TetGen
   */
  void computeHigherOrderPoint(
    const int                   index1, 
    const int                   index2,
    const int                   faceIndex, 
    const int                   faceMarker,
    const int                   globalHigherOrderIndex, 
    const vector<int>&          tetraCornerIndices,
    const vector<BlitzVec>&     tetraCornerNodes,
    const DRT::Element*         xfemElement, 
    const BlitzMat&             xyze_xfemElement,
    const map<int,BlitzVec3>&   currentcutterpositions, 
    tetgenio&                   out
    );


  /*!
  \brief  returns the other two point indices belonging to a triface
          with a Steiner point

  \param out                  (in)    :   output data structure of TetGen
  \param trifaceIndex         (in)    :   index of the triface
  \param steinerPointIndex    (in)    :   index of the Steiner point
  \return vector of point indices
  */
  vector<int> getPointIndices(
    const tetgenio&   out,
    const int         trifaceIndex,
    const int         steinerPointIndex
    ) const;


  /*!
  \brief  computes the recovery of the curved interface
          for one xfem element after a Contrained Delaunay
          Tetrahedralization.

  \param xsi                      (out)   :   intersection point
  \param normal                   (in)    :   nodes of the intersection normal
  \param cutterElement            (in)    :   cutter element
  \param xyze_cutterElement       (in)    :   nodal coordinates of cutter element
  \param onBoundary               (in)    :   indicates whether the normal lies on the xfem element
                                              boundary or not
  \return true if normal intersects with cutter element or false otherwise
  */
  bool computeRecoveryNormal(
    BlitzVec3&                                  xsi,
    const vector<BlitzVec>&                     normal,
    const DRT::Element*                         cutterElement,
    const BlitzMat&                             xyze_cutterElement,
    const bool                                  onBoundary
    ) const;


  /*!
  \brief  updates the system matrix for the intersection computation
          of an average normal and the corresponding cutter element
          to lift a point of the linearized interface onto the curved interface

  \param A                    (out)       :   system matrix
  \param xsi                  (in)        :   intersection point
  \param normal               (in)        :   nodes of the intersection normal
  \param surfaceElement       (in)        :   surface element
  \param xyze_surfaceElement  (in)        :   nodal coordinates of surface element
  \param onBoundary           (in)        :   true, if intersection normal lies on the
                                              element boundary
  */
  void updateAForRCINormal(
    BlitzMat3x3&                A,
    const BlitzVec3&            xsi,
    const vector<BlitzVec>&     normal,
    const DRT::Element*         surfaceElement,
    const BlitzMat&             xyze_surfaceElement,
    const bool                  onBoundary
    ) const;


  /*!
  \brief  updates the right-hand-side for the intersection computation
          of an average normal and the corresponding cutter element
          to lift a point of the linearized interface onto the curved interface

  \param b                    (out)       :   right-hand-side
  \param xsi                  (in)        :   intersection point
  \param normal               (in)        :   nodes of the intersection normal
  \param surfaceElement       (in)        :   surface element
  \param xyze_surfaceElement  (in)        :   nodal coordinates of surface element
  \param onBoundary           (in)        :   true, if intersection normal lies on the
                                              element boundary
  */
  void updateRHSForRCINormal(
    BlitzVec3&                b,
    const BlitzVec3&          xsi,
    const vector<BlitzVec>&   normal,
    const DRT::Element*       surfaceElement,
    const BlitzMat&           xyze_surfaceElement,
    const bool                onBoundary
    ) const;


  /*!
  \brief  computes the recovery of the curved interface
          for one xfem element after a Contrained Delaunay
          Tetrahedralization.

  \param lineIndex                (in)    :   line index
  \param currentcutterpositions   (in)    :   current cutter node positions
  \param xsi                      (out)   :   intersection point
  \param plane                    (in)    :   nodes of the intersection normal
  \param cutterElement            (in)    :   cutter element
  \param xyze_cutterElement       (in)    :   nodal coordinates of cutter element

  \return true if plane intersects with the line element or false otherwise
  */
  bool computeRecoveryPlane(
    int&                                        lineIndex,
    const map<int,BlitzVec3>&                   currentcutterpositions, 
    BlitzVec3&                                  xsi,
    const vector<BlitzVec>&                     plane,
    DRT::Element*                               cutterElement,
    const BlitzMat&                             xyze_cutterElement
    ) const;


  /*!
  \brief  updates the system matrix for the intersection computation
          of an average normal and the corresponding cutter element
          to lift a point of the linearized interface onto the curved interface

  \param A                        (out)       :   system matrix
  \param xsi                      (in)        :   intersection point
  \param plane                    (in)        :   nodes of the intersection normal
  \param lineElement              (in)        :   line element
  \param xyze_lineElement         (in)        :   nodal coordinates of line element
  \param surfaceElement           (in)        :   surface element
  \param xyze_surfaceElement      (in)        :   nodal coordinates of surface element
  */
  void updateAForRCIPlane(
    BlitzMat3x3&                A,
    const BlitzVec3&            xsi,
    const vector<BlitzVec>&     plane,
    const DRT::Element*         lineElement,
    const BlitzMat&             xyze_lineElement,
    const DRT::Element*         surfaceElement,
    const BlitzMat&             xyze_surfaceElement
    ) const;


  /*!
  \brief  updates the right-hand-side for the intersection computation
          of the cutter element edge with a plane through the midpoint
          of the face edge to lift a point of the linearized interface
          onto the curved interface

  \param b                        (out)       :   right-hand-side
  \param xsi                      (in)        :   intersection point
  \param plane                    (in)        :   nodes of the intersection normal
  \param lineElement              (in)        :   line element
  \param xyze_lineElement       (in)          :   nodal coordinates of line element
  */
  void updateRHSForRCIPlane(
    BlitzVec3&                  b,
    const BlitzVec3&            xsi,
    const vector<BlitzVec>&     plane,
    const DRT::Element*         lineElement,
    const BlitzMat&             xyze_lineElement
    ) const;


  /*!
  \brief  computes the normal to the interface edge of the
          adjacent tetrahedon. In addition the normal lies within the plane
          of the tetrahedron.

  \param onBoundary               (in)    : true if edge lies on the surface of the xfem element,
                                            false otherwise
  \param index1                   (in)    : local index 1 of an interface edge
  \param index2                   (in)    : local index 2 of an interface edge
  \param oppositePointIndex       (in)    : gloabal index of the point oopoiste to the edge
  \param globalHigherOrderIndex   (in)    : global higher order index
  \param tetraCornerIndices       (in)    : indices of tetrahedron corners
  \param tetraCornerNodes         (in)    : list of tetrahedron corner nodes in
                                            current configuration of the xfem element
  \param plane                    (out)   : nodes of the normal and plane
  \param xfemElement              (in)    : xfem element
  \param out                      (in)    : output data structure of tetGen
  */
  void computeIntersectionNormalA(
    const bool                              onBoundary,
    const int                               index1,
    const int                               index2,
    const int                               oppositePointIndex,
    const int                               globalHigherOrderIndex,
    const vector<int>&                      tetraCornerIndices,
    const vector<BlitzVec>&                 tetraCornerNodes,
    vector<BlitzVec>&                       plane,
    const DRT::Element*                     xfemElement,
    const tetgenio&                         out
    ) const;


  /*!
  \brief  computes the normal to the interface edge of two
          adjacent triangular faces passing through the midpoint of that edge.

  \param index1                   (in)    : global index 1 of an interface edge
  \param index2                   (in)    : global index 2 of an interface edge
  \param faceIndex                (in)    : index of the face
  \param adjacentFaceIndex        (in)    : index of the adjacent face
  \param plane                    (out)   : nodes of the normal and plane
  \param xfemElement              (in)    : xfem element
  \param out                      (in)    : output data structure of tetGen
  */
  void computeIntersectionNormalB(
    const int                                 index1,
    const int                                 index2,
    const int                                 faceIndex,
    const int                                 adjacentFaceIndex,
    const int                                 globalHigherOrderIndex,
    vector<BlitzVec>&                         plane,
    const DRT::Element*                       xfemElement,
    const tetgenio&                           out
    ) const;


  /*!
  \brief  computes the normal to the interface edge of a triangular
          faces on the xfem boundary passing through the Steiner point of that edge.
  
  \param steinerIndex             (in)    : index of Steiner point
  \param edgeIndex                (in)    : index of the other edge point
  \param oppositeIndex            (in)    : index of the point on the boundary triangle opposite the edge
  \param plane                    (out)   : nodes of the normal and plane
  \param xfemElement              (in)    : xfem element
  \param out                      (in)    : output data structure of tetGen
  */
  void computeIntersectionNormalC(
    const int                               steinerIndex,
    const int                               edgeIndex,
    const int                               oppositeIndex,
    vector<BlitzVec>&                       plane,
    const DRT::Element*                     xfemElement,
    const tetgenio&                         out
    ) const;


  /*!
  \brief  compute the midpoint of a line described by two points

  \param p1       (in)    : point 1
  \param p2       (in)    : point 2
  \return mid point of a line
  */
  BlitzVec computeLineMidpoint(
    const BlitzVec& p1,
    const BlitzVec& p2
    ) const;


  /*!
  \brief  finds the face marker of an adjacent facet for a given edge
          of another facet

  \param edgeIndex1               (in)   : index of edge 1
  \param edgeIndex2               (in)   : index of edge 2
  \param faceMarker               (in)   : index of the given facet
  \param adjacentFaceMarker       (in)   : index of the given facet
  \param faceIndex                (in)   : index of the given facet
  \param out                      (in)   : output data structure of tetGen
  \return the index of the face marker of the adjacent facet
  */
  void findAdjacentFace(
    const int         edgeIndex1,
    const int         edgeIndex2,
    const int         faceMarker,
    int&              adjacentFaceMarker,
    const int         faceIndex,
    int&              adjacentFaceIndex,
    const tetgenio&   out
    ) const;


  /*!
  \brief  searches for the common edge of two adjacent facets

  \param faceIndex1               (in)    : index of facet 1
  \param faceIndex2               (in)    : index of facet 2
  \param adjacentFacesList        (in)    : list of adjacent facets
  \param edgepoint                (in)    : edge point
  \param oppositepoint            (in)    : opposite point
  \param out                      (in)    : output data structure of tetGen
  \return true if common edge was found, false otherwise
  */
  bool findCommonFaceEdge(
    const int                     faceIndex1,
    const int                     faceIndex2,
    const vector<int>&            adjacentFacesList,
    BlitzVec&                     edgePoint,
    BlitzVec&                     oppositePoint,
    const tetgenio&               out
    ) const;


  /*!
  \brief  searches for the common line element of two
          adjacent cutter element

  \param currentcutterpositions   (in)    : current cutter node positions
  \param faceIndex1               (in)    : index of facet 1
  \param faceIndex2               (in)    : index of facet 2
  \param lineIndex                (out)   : index of the common line element
  \param surfaceIndex             (out)   : index of the cutter element
  \return true if common line element was found, false otherwise
  */
  bool findCommonCutterLine(
    const map<int,BlitzVec3>&                       currentcutterpositions,  
    const int                                       faceIndex1,
    const int                                       faceIndex2,
    int&                                            lineIndex,
    int&                                            cutterIndex
    ) const;


  /*!
  \brief  returns the index of a point lying opposite to a face edge

  \param edgeNode1                (in)   : edge node 1 of the facet edge
  \param edgeNode2                (in)   : edge node 2 of the facet edge
  \param adjacentFaceIndex        (in)   : list of adjacent face indices
  \param out                      (in)   : output data structure of tetGen
  \return the line index of the intersection line element of the given cutter element
  */
  int findEdgeOppositeIndex(
    const int                                 edgeIndex1,
    const int                                 edgeIndex2,
    const int                                 adjacentFaceIndex,
    const tetgenio&                           out
    ) const;


  /*!
  \brief  for the recovery computation of a higher-order node
          by plane - line element intersection
          this method finds the actually intersecting line
          for a given cutter element corresponding to a given facet edge

  \param xfemElement              (in)   : xfem element
  \param xyze_xfemElement         (in)   : nodal coordinates of xfem element
  \param cutterElement            (in)   : cutter element
  \param currentcutterpositions   (in)   : current cutter node positions
  \param edgeNode1                (in)   : edge node 1 of the facet edge
  \param edgeNode2                (in)   : edge node 2 of the facet edge
  \return the line index of the intersection line element of the given cutter element
  */
  int findIntersectingSurfaceEdge(
    const DRT::Element*                 xfemElement,
    const BlitzMat&                     xyze_xfemElement,
    DRT::Element*                       cutterElement,
    const map<int,BlitzVec3>&           currentcutterpositions,  
    const BlitzVec&                     edgeNode1,
    const BlitzVec&                     edgeNode2
    ) const;


  /*!
  \brief  stores the higher order node in the pointlist of the output data
          structure of tetGen out

  \param normal                   (in)    : true if itersection of a normal with a surface was computed
                                            false otherwise
  \param globalHigherOrderIndex   (in)    : higher order index global
  \param lineIndex                (in)    : lineIndex
  \param xsi                      (in)    : coordinates of the higher order node
  \param cutterElement            (in)    : cutterelement
  \param xyze_cutterElement       (in)    : nodal coordinates of cutter element
  \param xfem element             (in)    : xfem element
  \param out                      (in)    : output data structure of tetGen
  */
  void storeHigherOrderNode(
    const bool                                  normal,
    const int                                   globalHigherOrderIndex,
    const int                                   lineIndex,
    BlitzVec3&                                  xsi,
    DRT::Element*                               cutterElement,
    const map<int,BlitzVec3>&                   currentcutterpositions,  
    const DRT::Element*                         xfemElement,
    tetgenio&                                   out
    ) const;


  /*!
  \brief  stores domain integration cells

  \param xfemElement          (in)        :   xfem element
  \param domainintcells       (out)       :   domain integration cells
  \param out                  (in)        :   tetgenio output data structure
  \param higherorder          (in)        :   boll value to decide if linear or higherorder
  */
  void addCellsToDomainIntCellsMap(
    const DRT::Element*             xfemElement,
    map< int, DomainIntCells >&     domainintcells,
    const tetgenio&                 out,
    bool 							              higherorder
    ) const;


  /*!
  \brief  stores boundary integration cells

  \param trifaceIndex             (in)        :   index of the triangular face
  \param cornerIndex              (in)        :   global index of corner node of an triangle edge
  \param globalHigherOrderIndex   (in)        :   global higher-order index the midpoint of an triangle edge
  \param faceMarker   			      (in)        :   id of cutter element
  \param currentcutterpositions   (in)        :   current cutter node positions
  \param domainCoord            	(out)       :   point coordinates of the triangular face in xfem configuration
  \param boundaryCoord            (out)       :   point coordinates of the triangular face in cutter configuration
  \param xfemElement				      (in)		    :	  xfem element
  \param xyze_xfemElement         (in)        :   nodal coordinates of xfem element
  \param out                      (in)        :   tetgenio output data structure
  */
  void addCellsToBoundaryIntCellsMap(
  	const int                        trifaceIndex,
  	const int                        cornerIndex,
  	const int                        globalHigherOrderIndex,
  	const int                        faceMarker,
  	const map<int,BlitzVec3>&        currentcutterpositions,  
  	vector<vector<double> >&         domainCoord,
  	vector<vector<double> >&         boundaryCoord,
  	const DRT::Element*              xfemElement,
  	const BlitzMat&                  xyze_xfemElement,
  	const tetgenio&                	 out
    ) const;

        
  /*!
  \brief  stores boundary integration cells

  \param trifaceIndex             (in)        :   index of the triangular face
  \param cornerIndex              (in)        :   global index of corner node of an triangle edge
  \param globalHigherOrderIndex   (in)        :   global higher-order index the midpoint of an triangle edge
  \param faceMarker               (in)        :   id of cutter element
  \param currentcutterpositions   (in)        :   current cutter node positions
  \param domainCoord              (out)       :   point coordinates of the triangular face in xfem configuration
  \param boundaryCoord            (out)       :   point coordinates of the triangular face in cutter configuration
  \param xfemElement              (in)        :   xfem element
  \param xyze_xfemElement         (in)        :   nodal coordinates of xfem element
  \param out                      (in)        :   tetgenio output data structure
  */
  void addXFEMSurfaceCellsToBoundaryIntCellsMap(
    const bool                            higherorder,
    const int                             cornerIndex,
    const int                             index1,
    const int                             index2,
    const int                             cutterPos,
    const map<int,BlitzVec3>&             currentcutterpositions,
    vector<vector<double> >&              domainCoord,
    vector<vector<double> >&              boundaryCoord,
    const DRT::Element*                   xfemElement,
    const BlitzMat&                       xyze_xfemElement
    ) const;

        
  /*!
  \brief Debugging the intersection of XAABB's (DEBUG ONLY)

  \param cutterXAABB       (in)    : XAABB of the cutting element
  \param xfemXAABB         (in)    : XAABB of the xfem element
  \param cutterElement     (in)    : cutting element
  \param xfemElement       (in)    : xfem element
  \param noC               (in)    : id of the cutting element
  \param noX               (in)    : id of the xfem element
  */
  void debugXAABBIntersection(
    const BlitzMat3x2                 cutterXAABB,
    const BlitzMat3x2                 xfemXAABB,
    const DRT::Element*               cutterElement,
    const DRT::Element*               xfemElement,
    const int                         noC,
    const int                         noX
    ) const;


  /*!
  \brief Debugging node within element (DEBUG ONLY)

  \param element          (in)    : element
  \param node             (in)    : node
  \param xsi              (in)    : element coordinates
  \param noE              (in)    : id of the element
  \param noN              (in)    : id of the node
  \param within           (in)    : true if within, false otherwise
  */
  void debugNodeWithinElement(
    const DRT::Element*               xfemElement,
    const DRT::Node*                  node,
    const BlitzVec&                   xsi,
    const int                         noE,
    const int                         noN,
    const bool                        within
    ) const;


  /*!
  \brief Debugging curve surface intersection (DEBUG ONLY)

  \param surfaceElement   (in)    : surfaceElement
  \param lineElement      (in)    : lineElement
  \param xsi              (in)    : element coordinates
  \param noSE             (in)    : id of the surface element
  \param noLE             (in)    : id of the line element
  \param within           (in)    : true if within, false otherwise
  */
  void debugCurveSurfaceIntersection(
    const DRT::Element*         surfaceElement,
    const DRT::Element*         lineElement,
    BlitzVec&                   xsi,
    int                         noSE,
    int                         noLE,
    bool                        within
    ) const;


  /*!
  \brief Debugging tetgen data structure (DEBUG ONLY)

  \param element                  (in)        : xfemElement
  */
  void debugTetgenDataStructure(
    const DRT::Element*               xfemElement
    ) const;


  /*!
  \brief 	Debugging tetgen data output (DEBUG ONLY)
  		This method prints out the piecewise linear complex (PLC)
  		and the mesh of tetrahedrons for the list of specified
  		xfem elements. Output-files can be viewed with tetview


  \param in			           (in)       : input tetgen data structure
  \param out        	     (in)       : output tetgen data structure
  \param element     	     (in)       : xfem-element
  \param elementIds  	     (in)       : vector of element-ids
  \param timestepcounter   (in)       : timestepcounter
  */
 	void debugTetgenOutput(
 	  tetgenio&             in,
 	  tetgenio&             out,
 	  const DRT::Element*   xfemElement,
    vector<int>&          elementIds,
    int                   timestepcounter
    ) const;


  /*!
  \brief  computes unique coordinates for a certain region in the
          xfem element. It is copied into the data structure of Tetgen
           (DEBUG ONLY)

  \param index            (in)        : index of the xfem element
  \param out              (in)        : tetgenio output data structure
  */
  void printTetViewOutput(
    int             index,
    tetgenio&       out
    ) const;


  /*!
  \brief  computes unique coordinates for a certain region in the
          xfem element. It is copied into the data structure of Tetgen
           (DEBUG ONLY)

  \param index            (in)        : index of the xfem element
  \param out              (in)        : tetgenio output data structure
  */
  void printTetViewOutputPLC(
    const DRT::Element*   xfemElement,
    int             index,
    tetgenio&       in
    ) const;


  /*!
  \brief  debugs face marker (DEBUG ONLY)

  \param eleId            (in)        : id of the xfem element
  \param out              (in)        : tetgenio output data structure
  */
  void debugFaceMarker(
  	const int 				 eleId,
  	tetgenio&					 out
  	) const;

   
  /*!
  \brief  debug xfem condition (DEBUG ONLY)

  \param cutterdis            (in)        : cutter discretization
  */
  void debugXFEMConditions(
    const RCP<DRT::Discretization>  cutterdis
   	) const;

   
  /*!
  \brief  writes a gmsh file with xfem element and intersection candidates (DEBUG ONLY)

  \param xfemElement                  (in)        : xfem element
  \param cutterElements               (in)        : cutter elements whose XAABB's are intersecting
  */
  void debugIntersection(
    const DRT::Element*	       xfemElement,
    set<DRT::Element*>         cutterElements
    ) const;

  
  
  /*!
  \brief  writes gmsh file with xfem element and the currently
          intersecting cutter element (DEBUG ONLY)

  \param xfemElement                  (in)        : xfem element
  \param cutterElement                (in)        : cutter elements whose XAABB's are intersecting
  \param currentcutterpositions       (in)        : current nodal position of cutter element
  */
  void debugIntersectionOfSingleElements(
    const DRT::Element*                         xfemElement,
    const DRT::Element*                         cutterElement,
    const std::map<int,BlitzVec3>&              currentcutterpositions
     ) const;
       
  
  /*!
  \brief  debugs XAABBs (DEBUG ONLY)

  \param id                 (in)        : id
  \param cutterXAABB        (in)        : cutter XAABB
  \param xfemXAABB          (in)        : xfem XABB
  */
  void debugXAABBs(
    const int							              id,
    const BlitzMat&                     cutterXAABB,
    const BlitzMat&                     xfemXAABB
     ) const;

};
}

#endif  // #ifndef INTERSECTION_H
#endif  // #ifdef CCADISCRET


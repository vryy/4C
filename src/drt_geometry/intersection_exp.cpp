/*!----------------------------------------------------------------------
\file intersection.cpp

\brief  collection of intersection tools for the computation of the
        intersection of two arbitrary discretizations

    This is a copy of the original intersection and is not meant to replace it. It's rather a test bed
    and improvement can be merged into the original intersection classes

    The class intersection handles the intersection computation of
    Cartesian, linear and quadratic discretization. The discretiazation,
    which is intersected is refered to as xfem discretization and the
    discretization acting as a cutter is called cutter discretization.
    The intersection algorithm returns a list of quadratic integration cells
    for each intersected xfem element.

    The methods are categorized as follows for clearity:
    MAIN    public method which has to be called
            from outside to perform the intersection computation
    GM      general methods
    ICS     intersection candidate search
    CLI     contruction of the linearized interface
    CDT     contrained delaunay tetrahedralization
    RCI     recovery of the curved interface
    DB      debug methods

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "intersection_exp.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 |  MAIN:   computes the interface between the xfem          u.may 06/07|
 |          discretization and the cutter discretization.               |
 |          It returns a list of intersected xfem elements              |
 |          and their integrations cells.                               |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::computeIntersection(
    const Teuchos::RCP<DRT::Discretization>        xfemdis,
    const Teuchos::RCP<DRT::Discretization>        cutterdis,
    const std::map<int,LINALG::Matrix<3,1> >&      currentcutterpositions,
    std::map< int, DomainIntCells >&               domainintcells,
    std::map< int, BoundaryIntCells >&             boundaryintcells)
{

  TEUCHOS_FUNC_TIME_MONITOR(" GEO::IntersectionExp");

  std::cout << std::endl << "GEO::IntersectionExp:";
  flush(std::cout);
  
  countMissedPoints_ = 0;
  const double t_start = ds_cputime();
  
  // initialize tree for intersection candidates search
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(*cutterdis, currentcutterpositions);
  Teuchos::RCP<GEO::SearchTree> octTree = rcp(new GEO::SearchTree(8));
  octTree->initializeTree(rootBox, *cutterdis, GEO::TreeType(GEO::OCTTREE));
  
  // stop intersection if cutterdis is empty
  if(cutterdis->NumMyColElements()==0)
  {  
    std::cout << " cutter discretization empty ";
    std::cout << endl;
    flush(cout);
    return;
  }
  
  for(int k = 0; k < xfemdis->NumMyColElements(); ++k)
  {
    // printf("eleid = %d\n", k);
    DRT::Element* xfemElement = xfemdis->lColElement(k);
    initializeXFEM(k, xfemElement);
    EleGeoType xfemGeoType = HIGHERORDER;
    //checkGeoType(xfemElement, xyze_xfemElement_, xfemGeoType);

    startPointList();
  
    // serial search
    // const std::vector<int> cutterElementIds = serialIntersectionCandidateSearch(cutterdis, currentcutterpositions, xfemElement);
    // tree search for intersection candidates
    const std::vector<int> cutterElementIds = octTree->queryIntersectionCandidates(*cutterdis, currentcutterpositions, xfemElement, xyze_xfemElement_);
  
    // debugIntersection(xfemElement, cutterElements);
    const vector<RCP<DRT::Element> > xfemElementSurfaces = xfemElement->Surfaces();
    const vector<RCP<DRT::Element> > xfemElementLines = xfemElement->Lines();

    for(vector<int>::const_iterator id = cutterElementIds.begin(); id != cutterElementIds.end(); ++id )
    {
      DRT::Element* cutterElement = cutterdis->gElement(*id);
      cutterDistype_ = cutterElement->Shape();

      if(cutterElement == NULL) dserror("cutter element is null\n");
      const LINALG::SerialDenseMatrix xyze_cutterElement(GEO::getCurrentNodalPositions(cutterElement, currentcutterpositions));
      EleGeoType cutterGeoType = HIGHERORDER;
      //checkGeoType(cutterElement, xyze_cutterElement, cutterGeoType);
      const vector<RCP<DRT::Element> > cutterElementLines = cutterElement->Lines();
      const DRT::Node*const* cutterElementNodes = cutterElement->Nodes();

      // debugIntersectionOfSingleElements(xfemElement, cutterElement, currentcutterpositions);

      vector< InterfacePoint >  interfacePoints;

      // collect internal points
      for(int m=0; m<cutterElement->NumLine() ; m++)
        collectInternalPoints(cutterElement, cutterElementNodes[m], currentcutterpositions, interfacePoints, k, m);


      // collect intersection points
      for(int m=0; m<xfemElement->NumLine() ; m++)
      {
        //printf("cutsurf = %d\t xfemline = %d\n", countCutter, m);
        const bool doSVD = decideSVD(cutterGeoType, xfemGeoType);
        const DRT::Element* xfemElementLine = xfemElementLines[m].get();
        LINALG::SerialDenseMatrix xyze_xfemElementLine(GEO::InitialPositionArray(xfemElementLine));
        
        if(collectIntersectionPoints(   cutterElement, xyze_cutterElement,
                                        xfemElementLine, xyze_xfemElementLine, 
                                        interfacePoints, 0, m, false, doSVD))
        {
          storeIntersectedCutterElement(cutterElement);
        }
      }

      for(int m=0; m<cutterElement->NumLine() ; m++)
      {
        for(int p=0; p<xfemElement->NumSurface() ; p++)
        {
          //printf("cutline = %d\t xfemsurf = %d\n", m , p);
          const bool doSVD = decideSVD(xfemGeoType, cutterGeoType);
          const DRT::Element* xfemElementSurface = xfemElementSurfaces[p].get();
          const LINALG::SerialDenseMatrix xyze_xfemElementSurface(GEO::InitialPositionArray(xfemElementSurface));
          const DRT::Element* cutterElementLine = cutterElementLines[m].get();
          const LINALG::SerialDenseMatrix xyze_cutterElementLine(GEO::getCurrentNodalPositions(cutterElementLine, currentcutterpositions));
          
          if(collectIntersectionPoints(   xfemElementSurface, xyze_xfemElementSurface,
                                          cutterElementLine, xyze_cutterElementLine, 
                                          interfacePoints, p, m, true, doSVD))
          {
            storeIntersectedCutterElement(cutterElement);
          }
        }
      }


      // order interface points
      if( interfacePoints.size() > 0)
      {
#ifdef QHULL
        preparePLC(cutterElement, xyze_cutterElement, interfacePoints);
#else
        dserror("Set QHULL flag to use XFEM intersections!!!");
#endif
      }
    }// for-loop over all cutter elements


    if(checkIfCDT())
    {
      completePLC();
      //debugTetgenDataStructure(xfemElement);
      computeCDT(xfemElement, currentcutterpositions, domainintcells, boundaryintcells);
    }

  }// for-loop over all  actdis->NumMyColElements()

  //debugDomainIntCells(domainintcells,2);
  const double t_end = ds_cputime()-t_start;
  if(countMissedPoints_ > 0)
    cout << endl << "Number of missed points during the recovery copy = " << countMissedPoints_ << endl;

  std::cout << " Success (" << t_end  <<  " secs), intersected elements: " << domainintcells.size();
  std::cout << endl;
  flush(cout);
}


/*----------------------------------------------------------------------*
 |  INIT:   initializes the private members of the           u.may 08/07|
 |          current xfem element                                        |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::initializeXFEM(
    const int             xfemLid,
    const DRT::Element*   xfemElement)
{
  xfemDistype_ = xfemElement->Shape();
  
  if(xfemLid == 0)
    xfemOldDistype_ = xfemDistype_;
  
  // copy for first element or if previous element has a different element type
  if(xfemLid == 0 || xfemOldDistype_ != xfemDistype_ )
  {
    numXFEMSurfaces_ = xfemElement->NumSurface();
    numXFEMCornerNodes_  = DRT::UTILS::getNumberOfElementCornerNodes(xfemDistype_);
        
    eleLinesSurfaces_     = DRT::UTILS::getEleNodeNumbering_lines_surfaces(xfemDistype_);
    eleNodesSurfaces_     = DRT::UTILS::getEleNodeNumbering_nodes_surfaces(xfemDistype_);
    eleNodesLines_        = DRT::UTILS::getEleNodeNumbering_nodes_lines(xfemDistype_);
    eleNumberingLines_    = DRT::UTILS::getEleNodeNumberingLines(xfemDistype_);            // can be higher order
    eleNumberingSurfaces_ = DRT::UTILS::getEleNodeNumberingSurfaces(xfemDistype_);         // can be higher order
    eleRefCoordinates_    = DRT::UTILS::getEleNodeNumbering_nodes_reference(xfemDistype_); // can be higher order
    
    xfemOldDistype_ = xfemDistype_;
  }
  
  xyze_xfemElement_ = GEO::InitialPositionArray(xfemElement);
  
  pointList_.clear();
  triangleList_.clear();

  segmentList_.clear();
  segmentList_.resize(numXFEMSurfaces_);
  surfaceTriangleList_.clear();
  
  isolatedPointList_.clear();
  isolatedPointList_.resize(numXFEMSurfaces_);

  intersectingCutterElements_.clear();
  faceMarker_.clear();
  
  // clear efficient node map
  intersectionpointmap_ = rcp(new std::map< std::vector<int>, std::vector< LINALG::Matrix<3,1> >, CompareVecInt >(CompareVecInt()) ); 
  nodeInOut_.clear();
  nodemap_.clear();
}


/*----------------------------------------------------------------------*
 |  CLI:  serial search of intersection candidates           u.may 06/07|
 |        by comparing XAABB s                                          |
 *----------------------------------------------------------------------*/
std::vector<int> GEO::IntersectionExp::serialIntersectionCandidateSearch(
    const Teuchos::RCP<DRT::Discretization>         cutterdis,
    const std::map<int,LINALG::Matrix<3,1> >&       currentcutterpositions,
    DRT::Element*                                   xfemElement)
{
  EleGeoType xfemGeoType = HIGHERORDER;
  //checkGeoType(xfemElement, xyze_xfemElement_, xfemGeoType);
  const LINALG::Matrix<3,2> xfemXAABB = computeFastXAABB(xfemElement, xyze_xfemElement_, xfemGeoType);
      
  std::vector<int> cutterElementIds;
  // search for intersection candidates
  for(int kk = 0; kk < cutterdis->NumMyColElements(); ++kk)
  {
    DRT::Element*  cutterElement = cutterdis->lColElement(kk);
    if(cutterElement == NULL) dserror("geometry does not obtain elements");

    const LINALG::SerialDenseMatrix xyze_cutterElement(GEO::getCurrentNodalPositions(cutterElement, currentcutterpositions));
    EleGeoType cutterGeoType = HIGHERORDER;
    checkGeoType(cutterElement, xyze_cutterElement, cutterGeoType);
    const LINALG::Matrix<3,2>    cutterXAABB(computeFastXAABB(cutterElement, xyze_cutterElement, cutterGeoType));

    const bool intersected = intersectionOfXAABB(cutterXAABB, xfemXAABB);

    if(intersected)
      cutterElementIds.push_back(cutterElement->Id());

  }// for-loop over all cutterdis->NumMyColElements()

  return cutterElementIds;
}



/*----------------------------------------------------------------------*
 |  CLI:    collects points that belong to the interface     u.may 06/07|
 |          and lie within an xfem element                              |
 *----------------------------------------------------------------------*/
bool GEO::IntersectionExp::collectInternalPoints(
    DRT::Element*                         cutterElement,
    const DRT::Node*                      cutterNode,
    const map<int,LINALG::Matrix<3,1> >&  currentcutterpositions,
    std::vector< InterfacePoint >&        interfacePoints,
    const int                             elemId,
    const int                             nodeId)
{
  // current nodal position
  bool nodeWithinElement = false;
  const int cunoId = cutterNode->Id();
  static LINALG::Matrix<3,1> xsi;
        
  // check in node map if already computed an if it is in or out
  if(nodeInOut_.find(cunoId) != nodeInOut_.end())
  {
    // if nodewithin
    if(nodeInOut_.find(cunoId)->second)
    {
      xsi = nodemap_.find(cunoId)->second;
      nodeWithinElement = true;
    }
  }
  else
  { 
    const LINALG::Matrix<3,1> x = currentcutterpositions.find(cunoId)->second;
    xsi = currentToVolumeElementCoordinatesExact(xfemDistype_, xyze_xfemElement_, x, TOL7);
    nodeWithinElement = nodeInOut_[cunoId] = checkPositionWithinElementParameterSpace(xsi, xfemDistype_);
    if(nodeWithinElement)
      nodemap_[cunoId] = xsi;
  }
  

  if(nodeWithinElement)
  {
    InterfacePoint ip;
    // check if node lies on the boundary of the xfem element
    setInternalPointBoundaryStatus(xsi, ip);

    // intersection coordinates in the surface
    // element element coordinate system
    ip.setCoord(DRT::UTILS::getNodeCoordinates(nodeId, cutterElement->Shape()));

    //printf("internal point\n");
    //printf("xsi = %20.16f   %20.16f   %20.16f\n", xsi(0), xsi(1), xsi(2));
    interfacePoints.push_back(ip);
    storeIntersectedCutterElement(cutterElement);
  }
  return nodeWithinElement;
}



/*----------------------------------------------------------------------*
 |  CLI:    checks if a node that lies within an element     u.may 06/07|
 |          lies on one of its surfaces or nodes                        |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::setBoundaryPointBoundaryStatus(
    const LINALG::Matrix<3,1>&              xsi,
    InterfacePoint&                         ip) const
{
  bool onSurface = false;

  vector<int> surfaces = DRT::UTILS::getSurfaces(xsi, xfemDistype_);
  const int count = surfaces.size();

  // point lies on one surface
  if(count == 1)
  {
    onSurface = true;
    ip.setPointType(SURFACE);
    ip.setSurfaceId(surfaces);
  }
  // point lies on line, which has two neighbouring surfaces
  else if(count == 2)
  {
    onSurface = true;
    ip.setPointType(LINE);
    ip.setLineId(DRT::UTILS::getLines(xsi, xfemDistype_));
    ip.setSurfaceId(surfaces);
  }
  // point lies on a node, which has three neighbouring surfaces
  else if(count == 3)
  {   
    onSurface = true;
    ip.setPointType(NODE);
    ip.setNodeId(DRT::UTILS::getNode(xsi, xfemDistype_));
    ip.setLineId(DRT::UTILS::getLines(xsi, xfemDistype_));
    ip.setSurfaceId(surfaces);
  }
  else
    dserror("not on surface !!!");
}




/*----------------------------------------------------------------------*
 |  CLI:    checks if a node that lies within an element     u.may 06/07|
 |          lies on one of its surfaces or nodes                        |
 *----------------------------------------------------------------------*/
bool GEO::IntersectionExp::setInternalPointBoundaryStatus(
    const LINALG::Matrix<3,1>&              xsi,
    InterfacePoint&                         ip) const
{
  bool onSurface = false;

  vector<int> surfaces = DRT::UTILS::getSurfaces(xsi, xfemDistype_);
  const int count = surfaces.size();

  // point lies on one surface
  if(count == 1)
  {
    onSurface = true;
    ip.setPointType(SURFACE);
    ip.setSurfaceId(surfaces);
  }
  // point lies on line, which has two neighbouring surfaces
  else if(count == 2)
  {
    onSurface = true;
    ip.setPointType(LINE);
    ip.setLineId(DRT::UTILS::getLines(xsi, xfemDistype_));
    ip.setSurfaceId(surfaces);
  }
  // point lies on a node, which has three neighbouring surfaces
  else if(count == 3)
  {   
    onSurface = true;
    ip.setPointType(NODE);
    ip.setNodeId(DRT::UTILS::getNode(xsi, xfemDistype_));
    ip.setLineId(DRT::UTILS::getLines(xsi, xfemDistype_));
    ip.setSurfaceId(surfaces);
  }
  else
  {
    onSurface = false;
    ip.setPointType(INTERNAL);
  }
  return onSurface;
}



/*----------------------------------------------------------------------*
 |  CLI:    checks if a node that lies within an element     u.may 07/08|
 |          lies on one of its surfaces or nodes                        |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::setIntersectionPointBoundaryStatus( 
    const DRT::Element*                     surfaceElement,
    const LINALG::SerialDenseMatrix&        xyze_surfaceElement,
    const LINALG::Matrix<3,1>&              xsiSurface,
    InterfacePoint&                         ip) const
{

  LINALG::Matrix<3,1> x(true);

  // surface element is an xfem surface
  elementToCurrentCoordinates(surfaceElement->Shape(), xyze_surfaceElement, xsiSurface, x);
  LINALG::Matrix<3,1> xsi = currentToVolumeElementCoordinatesExact(xfemDistype_, xyze_xfemElement_, x, TOL7); 
  vector<int> surfaces = DRT::UTILS::getSurfaces(xsi, xfemDistype_);
  const int count = surfaces.size();

  // point lies on one surface
  if(count == 1)
  {
    ip.setPointType(SURFACE);
    ip.setSurfaceId(surfaces);
  }
  // point lies on line, which has two neighbouring surfaces
  else if(count == 2)
  {
    ip.setPointType(LINE);
    ip.setLineId(DRT::UTILS::getLines(xsi, xfemDistype_));
    ip.setSurfaceId(surfaces);
  }
  // point lies on a node, which has three neighbouring surfaces
  else if(count == 3)
  {
    ip.setPointType(NODE);
    ip.setLineId(DRT::UTILS::getLines(xsi, xfemDistype_));
    ip.setNodeId(DRT::UTILS::getNode(xsi, xfemDistype_));
    ip.setSurfaceId(surfaces);
  }
  else
    dserror("not on surface !!!");

}



/*----------------------------------------------------------------------*
 |  CLI:    collects all intersection points of a line and   u.may 06/07|
 |          and a surface                                               |
 *----------------------------------------------------------------------*/
bool GEO::IntersectionExp::collectIntersectionPoints(
    const DRT::Element*               surfaceElement,
    const LINALG::SerialDenseMatrix&  xyze_surfaceElement,
    const DRT::Element*               lineElement,
    const LINALG::SerialDenseMatrix&  xyze_lineElement,
    std::vector<InterfacePoint>&      interfacePoints,
    const int                         surfaceId,
    const int                         lineId,
    const bool                        lines,
    const bool                        doSVD
) const
{

  bool intersected = true;
  static LINALG::Matrix<3,1> xsi;
  xsi.Clear();

  static LINALG::Matrix<3,1> upLimit;
  static LINALG::Matrix<3,1> loLimit;

  // for hex elements
  upLimit.PutScalar(1.0);
  loLimit.PutScalar(-1.0);

  // check XAABBs of line and surface element
  intersected = checkLineSurfaceXAABBs(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement);
    if(!intersected)
      return false;
  
  // check line lies exactly in surface
  if(!doSVD)
    if(checkIfLineInSurface(surfaceElement, xyze_surfaceElement, xyze_lineElement))
      return false;

  intersected = computeCurveSurfaceIntersection(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, 
      upLimit, loLimit, xsi, doSVD, GEO::TOL7);

  if(intersected)
  {
    //printf("intersection point\n");
    //printf("xsi = %20.16f   %20.16f   %20.16f\n", xsi(0), xsi(1), xsi(2));
    addIntersectionPoint(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit,
        interfacePoints, surfaceId, lineId, lines, doSVD);

  } 
  return intersected;
}



/*----------------------------------------------------------------------*
 |  CLI:    checks if surface element is Cartesian and       u.may 06/08|
 |          line element is linear                                      |
 *----------------------------------------------------------------------*/
bool GEO::IntersectionExp::decideSVD(
    const EleGeoType surfaceGeoType,
    const EleGeoType lineGeoType)
{
  bool doSVD = true;

  if(surfaceGeoType == CARTESIAN &&  lineGeoType != HIGHERORDER)
    doSVD = false;

  return doSVD;
}



/*----------------------------------------------------------------------*
 |  CLI:    checks if a line lies in a Cartesian surface     u.may 06/08|
 |                                                                      |
 *----------------------------------------------------------------------*/
bool GEO::IntersectionExp::checkIfLineInSurface(
    const DRT::Element*                 surfaceElement,
    const LINALG::SerialDenseMatrix&    xyze_surfaceElement,
    const LINALG::SerialDenseMatrix&    xyze_lineElement) const
{
  bool        inSurface = true;
  double      distance = GEO::LARGENUMBER;
  LINALG::Matrix<3,1>   physCoord;
  LINALG::Matrix<3,1>   normal(true);
  LINALG::Matrix<2,1>   xsi(true); 


  for(int i = 0; i < 2; i++)
  {
    for(int j = 0; j < 3; j++)
      physCoord(j) = xyze_lineElement(j,i);

    searchForNearestPointOnSurface(surfaceElement, xyze_surfaceElement, physCoord, xsi, normal, distance);
    if(fabs(distance) > GEO::TOL7)
    {
      inSurface = false;
      break;
    }
  }
  return inSurface;
}



/*----------------------------------------------------------------------*
 |  CLI:  compare line and surface XAABB's                   u.may 06/07|
 |        before CSI                                                    |
 *----------------------------------------------------------------------*/
bool GEO::IntersectionExp::checkLineSurfaceXAABBs(
    const DRT::Element*               surfaceElement,
    const LINALG::SerialDenseMatrix&  xyze_surfaceElement,
    const DRT::Element*               lineElement,
    const LINALG::SerialDenseMatrix&  xyze_lineElement) const
{
  EleGeoType lineGeoType = HIGHERORDER;
  checkGeoType(lineElement, xyze_lineElement, lineGeoType);
  const LINALG::Matrix<3,2> lineXAABB = computeFastXAABB(lineElement, xyze_lineElement,lineGeoType);
        
  EleGeoType surfaceGeoType = HIGHERORDER;
  checkGeoType(surfaceElement, xyze_surfaceElement, surfaceGeoType);
  const LINALG::Matrix<3,2>  surfaceXAABB = computeFastXAABB(surfaceElement, xyze_surfaceElement, surfaceGeoType);

  return intersectionOfXAABB(surfaceXAABB, lineXAABB);
}

    
    
/*!
\brief updates the systemmatrix at the corresponding element coordinates
       for the computation of curve surface intersections
*/
template< DRT::Element::DiscretizationType surftype,
          DRT::Element::DiscretizationType linetype >
void updateAForCSI(
    LINALG::Matrix<3,3>&              A,                   ///< system matrix
    const LINALG::Matrix<3,1>&        xsi,                 ///< vector of element coordinates (r,s,t)
    const LINALG::SerialDenseMatrix&  xyze_surfaceElement, ///< nodal positions of surface element
    const LINALG::SerialDenseMatrix&  xyze_lineElement     ///< nodal positions of line element
)
{
  const int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<surftype>::numNodePerElement;
  const int numNodesLine = DRT::UTILS::DisTypeToNumNodePerEle<linetype>::numNodePerElement;

  A.Clear();

  static LINALG::Matrix<2,numNodesSurface> surfaceDeriv1;
  DRT::UTILS::shape_function_2D_deriv1(surfaceDeriv1, xsi(0), xsi(1), surftype);
  
  for(int inode=0; inode<numNodesSurface; inode++)
    for(int isd=0; isd<3; isd++)
    {
      A(isd,0) += xyze_surfaceElement(isd,inode) * surfaceDeriv1(0,inode);
      A(isd,1) += xyze_surfaceElement(isd,inode) * surfaceDeriv1(1,inode);
    }


  static LINALG::Matrix<2, numNodesLine> lineDeriv1;
  DRT::UTILS::shape_function_1D_deriv1(lineDeriv1, xsi(2), linetype);
  
  for(int inode=0; inode<numNodesLine; inode++)
    for(int isd=0; isd<3; isd++)
      A(isd,2) -= xyze_lineElement(isd,inode) * lineDeriv1(0,inode);
  
}



/*!
\brief updates the rhs at the corresponding element coordinates
       for the computation of curve surface intersections
*/
template<DRT::Element::DiscretizationType surftype,
         DRT::Element::DiscretizationType linetype>
void updateRHSForCSI(
    LINALG::Matrix<3,1>&              b,                   ///< right-hand-side
    const LINALG::Matrix<3,1>&        xsi,                 ///< vector of element coordinates (r,s,t)
    const LINALG::SerialDenseMatrix&  xyze_surfaceElement, ///< nodal positions of surface element
    const LINALG::SerialDenseMatrix&  xyze_lineElement     ///< nodal positions of line element
    )
{
  const int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<surftype>::numNodePerElement;
  const int numNodesLine = DRT::UTILS::DisTypeToNumNodePerEle<linetype>::numNodePerElement;

  b.Clear();

  static LINALG::Matrix<numNodesSurface,1> surfaceFunct;
  DRT::UTILS::shape_function_2D(surfaceFunct, xsi(0), xsi(1), surftype);
  for(int i=0; i<numNodesSurface; i++)
  {
    for(int dim=0; dim<3; dim++)
      b(dim) -= xyze_surfaceElement(dim,i) * surfaceFunct(i);
  }

  static LINALG::Matrix<numNodesLine,1> lineFunct;
  DRT::UTILS::shape_function_1D(lineFunct, xsi(2), linetype);
  for(int i=0; i<numNodesLine; i++)
  {
    for(int dim=0; dim<3; dim++)
      b(dim) += xyze_lineElement(dim,i) * lineFunct(i);
  }
}


/*!
    \brief solves a singular system of equations

  \return true if resulting system is singular , false otherwise
*/
template<DRT::Element::DiscretizationType surftype,
         DRT::Element::DiscretizationType linetype>
bool computeSingularCSI(
    LINALG::Matrix<3,1>&                  xsi,                 ///< vector of element coordinates (r,s,t)
    const LINALG::SerialDenseMatrix&      xyze_surfaceElement, ///< nodal positions of surface element
    const LINALG::SerialDenseMatrix&      xyze_lineElement     ///< nodal positions of line element
    )
{
  bool singular = true;
  int iter = 0;
  const int maxiter = 5;
  double residual = 1.0;
  static LINALG::Matrix<3,3> A;
  static LINALG::Matrix<3,1> b;
  static LINALG::Matrix<3,1> dx;

  updateRHSForCSI<surftype,linetype>( b, xsi, xyze_surfaceElement, xyze_lineElement);

  while(residual > GEO::TOL13)
  {
    updateAForCSI<surftype,linetype>( A, xsi, xyze_surfaceElement, xyze_lineElement);

    dx.Clear();
    if(GEO::solveLinearSystemWithSVD<3>(A, b, dx, GEO::TOL14))
    {
      singular = false;
      xsi += dx;
      break;
    }

    xsi += dx;
    updateRHSForCSI<surftype,linetype>( b, xsi, xyze_surfaceElement, xyze_lineElement);
    residual = b.Norm2();
    iter++;

    if(iter >= maxiter )
    {
      singular = true;
      break;
    }
  }
  return singular;
}


/*!
\brief computes an interseticon point between a curve and a surface - templated part

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
template<DRT::Element::DiscretizationType surftype,
         DRT::Element::DiscretizationType linetype>
bool computeCurveSurfaceIntersectionT(
    const DRT::Element*               surfaceElement,
    const LINALG::SerialDenseMatrix&  xyze_surfaceElement,
    const DRT::Element*               lineElement,
    const LINALG::SerialDenseMatrix&  xyze_lineElement,
    LINALG::Matrix<3,1>&              xsi,
    const LINALG::Matrix<3,1>&        upLimit,
    const LINALG::Matrix<3,1>&        loLimit,
    const bool                        doSVD,
    const double                      tol
)
{
  if (surfaceElement->Shape() != surftype) dserror("bug in template instantiation");
  if (lineElement->Shape() != linetype) dserror("bug in template instantiation");

  bool singular = false;
  bool intersection = true;
  int iter = 0;
  const int maxiter = 20;
  double residual = 1.0;
  static LINALG::Matrix<3,3> A;
  static LINALG::Matrix<3,1> b;
  static LINALG::Matrix<3,1> dx;  dx = 0.0;


  updateRHSForCSI<surftype,linetype>( b, xsi, xyze_surfaceElement, xyze_lineElement);

  while(residual > GEO::TOL13 && !singular)
  {
    updateAForCSI<surftype,linetype>( A, xsi, xyze_surfaceElement, xyze_lineElement);
    singular = !GEO::gaussElimination<true,3>(A, b, dx, GEO::TOL14);

    if(singular && !doSVD)
    {
      intersection = false;
      break;
    }
    else if(singular && doSVD)
    {
      if(computeSingularCSI<surftype,linetype>(xsi, xyze_surfaceElement, xyze_lineElement))
      {
        intersection = false;
        iter = maxiter + 1;
      }          
      dx.Clear();
    }

    //cout << "SINGULAR << endl;
    xsi += dx;
    updateRHSForCSI<surftype,linetype>( b, xsi, xyze_surfaceElement, xyze_lineElement);
    residual = b.Norm2();
    iter++;

    //printf("xsi = %20.16f   %20.16f   %20.16f  iter = %d res = %20.16f\n", xsi(0), xsi(1), xsi(2), iter, residual  );
    // has to to be 8 , otherwise not a number is reached               // detect not a number according to IEEE NaN is not comparable to itself
    if(iter >= maxiter || GEO::SumOfFabsEntries(xsi) > GEO::TOLPLUS8 || !(xsi(0)==xsi(0))  || !(xsi(1)==xsi(1))  || !(xsi(2)==xsi(2))   )
    {
      intersection = false;
      break;
    }
  }

  //printf("xsi = %f   %f   %f  iter = %d res = %20.16f\n", xsi(0), xsi(1), xsi(2), iter, residual  );
  if(intersection)
  {
    if( (xsi(0) > (upLimit(0)+tol)) || (xsi(1) > (upLimit(1)+tol)) || (xsi(2) > (upLimit(2)+tol))  ||
        (xsi(0) < (loLimit(0)-tol)) || (xsi(1) < (loLimit(1)-tol)) || (xsi(2) < (loLimit(2)-tol)))
      intersection = false;
  }

  return intersection;
}



/*----------------------------------------------------------------------*
 |  CLI:    computes the intersection between a              u.may 06/07|
 |          curve and a surface                    (CSI)                |
 *----------------------------------------------------------------------*/
bool GEO::IntersectionExp::computeCurveSurfaceIntersection(
    const DRT::Element*               surfaceElement,
    const LINALG::SerialDenseMatrix&  xyze_surfaceElement,
    const DRT::Element*               lineElement,
    const LINALG::SerialDenseMatrix&  xyze_lineElement,
    const LINALG::Matrix<3,1>&        upLimit,
    const LINALG::Matrix<3,1>&        loLimit,
    LINALG::Matrix<3,1>&              xsi,
    const bool                        doSVD,
    const double                      tol
) const
{
  if (lineElement->Shape() == DRT::Element::line2)
  {
    switch (surfaceElement->Shape())
    {
    case DRT::Element::quad4:
      return computeCurveSurfaceIntersectionT<DRT::Element::quad4,DRT::Element::line2>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD, tol);
    case DRT::Element::quad8:
      return computeCurveSurfaceIntersectionT<DRT::Element::quad8,DRT::Element::line2>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD, tol);
    case DRT::Element::quad9:
      return computeCurveSurfaceIntersectionT<DRT::Element::quad9,DRT::Element::line2>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD, tol);
    case DRT::Element::tri3:
      return computeCurveSurfaceIntersectionT<DRT::Element::tri3 ,DRT::Element::line2>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD, tol);
    case DRT::Element::tri6:
      return computeCurveSurfaceIntersectionT<DRT::Element::tri6 ,DRT::Element::line2>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD, tol);
    default:
      dserror("template not instatiated yet");
      return false;
    };
  }
  else if (lineElement->Shape() == DRT::Element::line3)
  {
    switch (surfaceElement->Shape())
    {
    case DRT::Element::quad4:
      return computeCurveSurfaceIntersectionT<DRT::Element::quad4,DRT::Element::line3>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD, tol);
    case DRT::Element::quad8:
      return computeCurveSurfaceIntersectionT<DRT::Element::quad8,DRT::Element::line3>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD, tol);
    case DRT::Element::quad9:
      return computeCurveSurfaceIntersectionT<DRT::Element::quad9,DRT::Element::line3>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD, tol);
    case DRT::Element::tri3:
      return computeCurveSurfaceIntersectionT<DRT::Element::tri3 ,DRT::Element::line3>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD, tol);
    case DRT::Element::tri6:
      return computeCurveSurfaceIntersectionT<DRT::Element::tri6 ,DRT::Element::line3>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit, doSVD, tol);
    default:
      dserror("template not instatiated yet");
      return false;
    };
  }
  return true;
}



/*----------------------------------------------------------------------*
 |  CLI:    computes a new starting point for the            u.may 06/07|
 |          Newton-method in order to find all intersection points      |
 |          of a curve-surface intersection                             |
 *----------------------------------------------------------------------*/
int GEO::IntersectionExp::computeNewStartingPoint(
    const DRT::Element*                surfaceElement,
    const LINALG::SerialDenseMatrix&   xyze_surfaceElement,
    const DRT::Element*                lineElement,
    const LINALG::SerialDenseMatrix&   xyze_lineElement,
    const int                          surfaceId,
    const int                          lineId,
    const LINALG::Matrix<3,1>&         xsiOld,
    const LINALG::Matrix<3,1>&         upLimit,
    const LINALG::Matrix<3,1>&         loLimit,
    std::vector<InterfacePoint>&       interfacePoints,
    const bool                         lines,
    const bool                         doSVD) const
{
  bool interval = true;
  int numInterfacePoints = 0;

  if(comparePoints<3>(upLimit, loLimit))
    interval = false;

  static LINALG::Matrix<3,1> xsi;
  xsi.Update(0.5, upLimit, 0.5, loLimit);
  
  //printf("xsi = %f   %f   %f\n", fabs(xsi(0)), fabs(xsi(1)), fabs(xsi(2)) );
  //printf("lolimit = %f   %f   %f\n", fabs(loLimit(0)), fabs(loLimit(1)), fabs(loLimit(2)) );
  //printf("uplimit = %f   %f   %f\n", fabs(upLimit(0)), fabs(upLimit(1)), fabs(upLimit(2)) );

  bool intersected = computeCurveSurfaceIntersection(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, 
      upLimit, loLimit, xsi, doSVD, GEO::TOL7);

  if( comparePoints<3>(xsi, xsiOld))
    intersected = false;

  if(intersected && interval)
    numInterfacePoints = addIntersectionPoint(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit,
        interfacePoints, surfaceId, lineId, lines, doSVD);

  //printf("number of intersection points = %d\n", numInterfacePoints );
  return numInterfacePoints;
}



/*----------------------------------------------------------------------*
 |  CLI:    adds an intersection point to the 			     u.may 07/07|
 |          list of interface points                       			    |
 *----------------------------------------------------------------------*/
int GEO::IntersectionExp::addIntersectionPoint(
    const DRT::Element*               surfaceElement,
    const LINALG::SerialDenseMatrix&  xyze_surfaceElement,
    const DRT::Element*               lineElement,
    const LINALG::SerialDenseMatrix&  xyze_lineElement,
    const LINALG::Matrix<3,1>&        xsi,
    const LINALG::Matrix<3,1>&        upLimit,
    const LINALG::Matrix<3,1>&        loLimit,
    std::vector<InterfacePoint>&      interfacePoints,
    const int                         surfaceId,
    const int                         lineId,
    const bool                        lines,
    const bool                        doSVD) const
{

  InterfacePoint ip;
  
  // cutter line with xfem surface
  if(lines)
  {          
    setIntersectionPointBoundaryStatus(surfaceElement,xyze_surfaceElement,xsi, ip);
    ip.setCoord(DRT::UTILS::getLineCoordinates(lineId, xsi(2), cutterDistype_));
  }
  // xfem line with cutter surface
  else
  {
    // check if point lies on a node of the xfem line element and therefore also on
    // the xfem element
    int lineNodeId = -1;
    if(fabs(xsi(2) + 1.0) < GEO::TOL7)
      lineNodeId = 0;

    if(fabs(xsi(2) - 1.0) < GEO::TOL7)
      lineNodeId = 1; 

    if(lineNodeId > -1)
    {
      const int nodeId = eleNumberingLines_[lineId][lineNodeId];
      // point type has to be set before ids and coords are set
      ip.setPointType(NODE);
      ip.setNodeId(nodeId);
      ip.setLineId(eleNodesLines_[nodeId]);
      ip.setSurfaceId(eleNodesSurfaces_[nodeId]);
      LINALG::Matrix<3,1> coord;
      coord(0) = xsi(0);
      coord(1) = xsi(1);
      coord(2) = 0.0;
      ip.setCoord(coord);
    }
    else
    {
      ip.setPointType(LINE);
      vector<int> lineVec(1,lineId);
      ip.setLineId(lineVec);
      ip.setSurfaceId(eleLinesSurfaces_[lineId]);
      LINALG::Matrix<3,1> coord;
      coord(0) = xsi(0);
      coord(1) = xsi(1);
      coord(2) = 0.0;
      ip.setCoord(coord);
    }
  }

  vector<InterfacePoint>::const_iterator it;
  bool alreadyInList = false;
  for(it = interfacePoints.begin(); it != interfacePoints.end(); it++ )
    if(comparePoints<3>(ip.getCoord(), it->getCoord()))
    {
      //printf("alreadyinlist = true\n");
      alreadyInList = true;
      break;
    }

  int numInterfacePoints = 0;
  if(!alreadyInList)
  {
    vector< LINALG::Matrix<3,1> >  upperLimits(8, LINALG::Matrix<3,1>(true));
    vector< LINALG::Matrix<3,1> >  lowerLimits(8, LINALG::Matrix<3,1>(true));
    createNewLimits(xsi, upLimit, loLimit, upperLimits, lowerLimits);

    interfacePoints.push_back(ip);
    numInterfacePoints++;

    // recursive call
    // for linear lines and Cartesian surfaces no more than one intersection point can be expected
    if(doSVD)
      for(int i = 0; i < 8; i++)
        numInterfacePoints += computeNewStartingPoint(  surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, 
                                                        surfaceId, lineId, xsi,
                                                        upperLimits[i], lowerLimits[i], interfacePoints, lines, doSVD);

  }
  return numInterfacePoints;
}



/*----------------------------------------------------------------------*
 |  CLI:    create new ranges for the recursive              u.may 07/07|
 |          computation of all intersection points                      |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::createNewLimits(
    const LINALG::Matrix<3,1>&        xsi,
    const LINALG::Matrix<3,1>&        upLimit,
    const LINALG::Matrix<3,1>&        loLimit,
    vector< LINALG::Matrix<3,1> >&    upperLimits,
    vector< LINALG::Matrix<3,1> >&    lowerLimits) const
    {


  /*          Surface:                                Line:
   *        (-1, 1)               (1,1)
   *          0_____________________1
   *          |          s          |
   *          |         /\          |
   *          |          |          |                 4 ___________x__________ 5
   *          |          |          |              ( -1 )                    ( 1 )
   *          |          x ----> r  |
   *          |                     |
   *          |                     |
   *          |                     |
   *          2_____________________3
   *        (-1,-1)                (1,-1)
   */

  // upper left corner of surface with lower part of line
  upperLimits[0](0) = xsi(0);         lowerLimits[0](0) = loLimit(0);
  upperLimits[0](1) = upLimit(1);     lowerLimits[0](1) = xsi(1);
  upperLimits[0](2) = xsi(2);         lowerLimits[0](2) = loLimit(2);

  // upper left corner of surface with upper part of line
  upperLimits[1](0) = xsi(0);         lowerLimits[1](0) = loLimit(0);
  upperLimits[1](1) = upLimit(1);     lowerLimits[1](1) = xsi(1);
  upperLimits[1](2) = upLimit(2);     lowerLimits[1](2) = xsi(2);

  // upper right corner of surface with lower part of line
  upperLimits[2](0) = upLimit(0);     lowerLimits[2](0) = xsi(0);
  upperLimits[2](1) = upLimit(1);     lowerLimits[2](1) = xsi(1);
  upperLimits[2](2) = xsi(2);         lowerLimits[2](2) = loLimit(2);

  // upper right corner of surface with upper part of line
  upperLimits[3](0) = upLimit(0);     lowerLimits[3](0) = xsi(0);
  upperLimits[3](1) = upLimit(1);     lowerLimits[3](1) = xsi(1);
  upperLimits[3](2) = upLimit(2);     lowerLimits[3](2) = xsi(2);

  // lower right corner of surface with lower part of line
  upperLimits[4](0) = upLimit(0);     lowerLimits[4](0) = xsi(0);
  upperLimits[4](1) = xsi(1);         lowerLimits[4](1) = loLimit(1);
  upperLimits[4](2) = xsi(2);         lowerLimits[4](2) = loLimit(2);

  // lower right corner of surface with upper part of line
  upperLimits[5](0) = upLimit(0);     lowerLimits[5](0) = xsi(0);
  upperLimits[5](1) = xsi(1);         lowerLimits[5](1) = loLimit(1);
  upperLimits[5](2) = upLimit(2);     lowerLimits[5](2) = xsi(2);

  // lower left corner of surface with lower part of line
  upperLimits[6](0) = xsi(0);         lowerLimits[6](0) = loLimit(0);
  upperLimits[6](1) = xsi(1);         lowerLimits[6](1) = loLimit(1);
  upperLimits[6](2) = xsi(2);         lowerLimits[6](2) = loLimit(2);

  // lower left corner of surface with upper part of line
  upperLimits[7](0) = xsi(0);         lowerLimits[7](0) = loLimit(0);
  upperLimits[7](1) = xsi(1);         lowerLimits[7](1) = loLimit(1);
  upperLimits[7](2) = upLimit(2);     lowerLimits[7](2) = xsi(2);
}



/*----------------------------------------------------------------------*
 |  CLI:    determines the surface Id of an xfem surface, if u.may 06/08|
 |          all interface points are lying on this surface              |
 *----------------------------------------------------------------------*/
int GEO::IntersectionExp::findCommonSurfaceID(
    const DRT::Element*                 cutterElement,
    const LINALG::SerialDenseMatrix&    xyze_cutterElement,
    const vector<int>&                  positions)
{
  int surfId = -1;

  vector<vector<int> > xfemSurfPoints;
  xfemSurfPoints.resize(numXFEMSurfaces_);

  // fill data structure store all point position lying on a certain xfem surface
  for(unsigned int i = 0; i < positions.size(); i++)
  {
    const int pos = positions[i];
    for(int j = 0; j < pointList_[pos].getNumSurface(); j++)
    {
      const int surface = pointList_[pos].getSurfId()[j];
      xfemSurfPoints[surface].push_back(pos);
    }
  }
  // check if more than 2 points are lying on one xfem surface
  for(int i = 0; i < numXFEMSurfaces_; i++)
  {
    //for(unsigned int j = 0; j < xfemSurfPoints[i].size(); j++ )
    //  cout << "xfemsurf.  " << j << "        xfemSurfPoints =  "  <<  xfemSurfPoints[i][j] << endl ;

    if(xfemSurfPoints[i].size() > 2  &&  xfemSurfPoints[i].size() !=  positions.size())
    {
      for(unsigned int j = 0; j < xfemSurfPoints[i].size(); j++ )
      {
        for(int k = 0; k < 3; k++)
          cout << "point        " << pointList_[(int) xfemSurfPoints[i][j]].getCoord()(k);

        cout << endl;
      }
      cout << endl;

      printf("scenario not yet implemented\n");
    }
    else if(xfemSurfPoints[i].size() > 2 &&  xfemSurfPoints[i].size() ==  positions.size())
    {
      const bool onSurface = checkIfCutterOnXFEMSurface(cutterElement, xyze_cutterElement, positions);
      if(onSurface)
        surfId = i;
    }
  }   
  return surfId;
}



/*----------------------------------------------------------------------*
 |  CLI:    checks if the part of a cutter element           u.may 07/08|
 |          spezified by position lies on a xfem surface                |                      
 |          by checking if the cutter midpoint lies on the              |
 |          xfem surface                                                |
 *----------------------------------------------------------------------*/
bool GEO::IntersectionExp::checkIfCutterOnXFEMSurface(
    const DRT::Element*               cutterElement,
    const LINALG::SerialDenseMatrix&  xyze_cutterElement,
    const vector<int>&                positions) const
{
  bool onSurface = false;

  LINALG::Matrix<3,1> x_phys(true);

  // midpoint is computed in element coordinates of the xfem element
  InterfacePoint midpoint = computeMidpoint(positions);

  // transform to physical coordinates
  elementToCurrentCoordinates(xfemDistype_, xyze_xfemElement_, midpoint.getCoord(), x_phys);

  // check if midpoint lies on cutter element
  LINALG::Matrix<2,1> xsiCut(true);
  LINALG::Matrix<3,1> normal(true);
  double distance = 0.0;
  searchForNearestPointOnSurface(cutterElement, xyze_cutterElement, x_phys, xsiCut, normal, distance);

  if(fabs(distance) < GEO::TOL7)
    onSurface = true;

  return onSurface;
}


/*----------------------------------------------------------------------*
 |  ICS:    prepares the part of a piecewise linear complex  u.may 06/08|
 |          for a xfem and a cutter element                             |
 *----------------------------------------------------------------------*/
#ifdef QHULL
void GEO::IntersectionExp::preparePLC(
    const DRT::Element*                       cutterElement,
    const LINALG::SerialDenseMatrix&          xyze_cutterElement,
    vector<InterfacePoint>&                   interfacePoints)
{

  InterfacePoint                    midpoint; // only used when more than 2 interface points are existent
  vector< LINALG::Matrix<3,1> >     vertices;
  // for more then two interface points compute the convex hull
  // and store ordered points in vertices

  if(interfacePoints.size() > 2)
  {
    
    //compute midpoint

    // tolerance has to be twice as small than for other points because the midpoint is 
    // point by summing the other
    // points and dividing by the number of points, other wise midpoint 
    // is moved on xfem boundary even though is is still inside
    // xfem element
    midpoint = computeMidpoint(interfacePoints);
    // transform midpoint into current coordinates
    {
      static LINALG::Matrix<2,1>    eleCoordSurf;
      for(int j = 0; j < 2; j++)
        eleCoordSurf(j)  = midpoint.getCoord()(j);
      static LINALG::Matrix<3,1> curCoordVol;
      elementToCurrentCoordinates(cutterElement->Shape(), xyze_cutterElement, eleCoordSurf, curCoordVol);
      const LINALG::Matrix<3,1> eleCoordVol(currentToVolumeElementCoordinatesExact(xfemDistype_, xyze_xfemElement_, curCoordVol, TOL14));
      midpoint.setCoord(eleCoordVol);
    }    
    computeConvexHull(cutterElement, xyze_cutterElement,interfacePoints, vertices);
//    computeDelaunayTriangulation(cutterElement, xyze_cutterElement,interfacePoints, vertices, triangles);
  }
  
  // transform interface points into current coordinates
  for(vector<InterfacePoint>::iterator ipoint = interfacePoints.begin(); ipoint != interfacePoints.end(); ++ipoint)
  {
    static LINALG::Matrix<2,1>  eleCoordSurf;
    for(int j = 0; j < 2; j++)
      eleCoordSurf(j)  = ipoint->getCoord()(j);
    static LINALG::Matrix<3,1> curCoordVol;
    elementToCurrentCoordinates(cutterElement->Shape(), xyze_cutterElement, eleCoordSurf, curCoordVol);
    const LINALG::Matrix<3,1> eleCoordVol(currentToVolumeElementCoordinatesExact(xfemDistype_, xyze_xfemElement_, curCoordVol, TOL7));
    ipoint->setCoord(eleCoordVol);
  }
  
  // for 1 or 2 interface points (line segment or isolated point)
  if(interfacePoints.size() <= 2 && interfacePoints.size() > 0) 
  {
    for(vector<InterfacePoint>::const_iterator ipoint = interfacePoints.begin(); ipoint != interfacePoints.end(); ++ipoint)
    {
      // store in structure vertices
      vertices.push_back(ipoint->getCoord());
    }
  }

  if(interfacePoints.size() > 1)
  {
    // store pointList_
    vector<int> positions;
    storePointList(vertices, positions, interfacePoints);

    // find common surfID, if surfId != -1 all interface points are lying on one
    // xfem surface and have to be store in the surfaceTriangleList_ accordingly
    const int surfId = findCommonSurfaceID(cutterElement, xyze_cutterElement,  positions);

    // store part of PLC
    storePLC(cutterElement, xyze_cutterElement, surfId, positions, midpoint);
    //storePLC(cutterElement, xyze_cutterElement, surfId, surfaceInterfacePoints, vertices, midpoint);
  }

  // clear interface points
  interfacePoints.clear();
}



/*----------------------------------------------------------------------*
 |  ICS:    computes the convex hull of a set of             u.may 06/07|
 |          interface points and stores resulting points,               |
 |          segments and triangles for the use with Tetgen (CDT)        |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::computeConvexHull(
          const DRT::Element*                 cutterElement,
          const LINALG::SerialDenseMatrix&    xyze_cutterElement,
          const vector<InterfacePoint>&       interfacePoints,
          vector< LINALG::Matrix<3,1> >&      vertices) const
{
  // store coordinates in
  // points has numInterfacePoints*dim-dimensional components
  // points[0] is the first coordinate of the first point
  // points[1] is the second coordinate of the first point
  // points[dim] is the first coordinate of the second point
  coordT* coordinates = (coordT *)malloc((2*interfacePoints.size())*sizeof(coordT));
  int fill = 0;
  for(vector<InterfacePoint>::const_iterator ipoint = interfacePoints.begin(); ipoint != interfacePoints.end(); ++ipoint)
  {
    for(int j = 0; j < 2; j++)
    {
      coordinates[fill++] = ipoint->getCoord()(j);
      //printf("coord = %f\t", ipoint->getCoord()[j]);
    }
    //cout << endl;
  }

  // compute convex hull - exitcode = 0 no error
  if (qh_new_qhull(2, interfacePoints.size(), coordinates, false, "qhull ", NULL, stderr)!=0)
    dserror(" error in the computation of the convex hull (qhull error)");

//  qh_setdelaunay(2,interfacePoints.size(),coordinates);
//
//  // copy vertices out of the facet list
//  facetT* facet = qh facet_list;
//  for(int i = 0; i< qh num_facets; i++)
//  {
//    for(int j = 0; j < 3; j++) //3 points per triangle
//    {
//      double* point  = SETelemt_(facet->vertices, j, vertexT)->point;
//
//      static LINALG::Matrix<2,1>  eleCoordSurf;
//      for(int m = 0; m < 2; m++)
//        eleCoordSurf(m)  = point[m];
//      static LINALG::Matrix<3,1> curCoordVol;
//      elementToCurrentCoordinates(cutterElement->Shape(), xyze_cutterElement, eleCoordSurf, curCoordVol);
//      const LINALG::Matrix<3,1> eleCoordVol(currentToVolumeElementCoordinatesExact(xfemDistype_, xyze_xfemElement_, curCoordVol, TOL7));
//      vertices.push_back(eleCoordVol);
//    }
////
////    vector<int> triangle(3);
////    triangle[0] = vertices.size()-2;
////    triangle[1] = vertices.size()-1;
////    triangle[2] = vertices.size();
////
////    triangles.push_back(triangle);
////
//    facet = facet->next;
//  }
  
  // copy vertices out of the facet list
  facetT* facet = qh facet_list;
  for(int i = 0; i< qh num_facets; i++)
  {
    for(int j = 0; j < 2; j++)
    {
      double* point  = SETelemt_(facet->vertices, j, vertexT)->point;
      
      static LINALG::Matrix<2,1>  eleCoordSurf;
      for(int m = 0; m < 2; m++)
        eleCoordSurf(m)  = point[m];
      static LINALG::Matrix<3,1> curCoordVol;
      elementToCurrentCoordinates(cutterElement->Shape(), xyze_cutterElement, eleCoordSurf, curCoordVol);
      const LINALG::Matrix<3,1> eleCoordVol(currentToVolumeElementCoordinatesExact(xfemDistype_, xyze_xfemElement_, curCoordVol, TOL7));
      vertices.push_back(eleCoordVol);
    }
    facet = facet->next;
  }

  // for debugging if points are lying on the convex hull intersection conzinues with out any problems
  if(((int) interfacePoints.size()) != qh num_vertices)
    printf("resulting surface is concave - convex hull does not include all points\n");


  // free memory and clear vector of interface points
  qh_freeqhull(!qh_ALL);
  int curlong, totlong;           // memory remaining after qh_memfreeshort
  qh_memfreeshort (&curlong, &totlong);
  if (curlong || totlong)
    printf("qhull internal warning (main): did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);

  free(coordinates);
}
#endif //QHULL





/*----------------------------------------------------------------------*
 |  ICS:    stores the part of a piecewise linear complex    u.may 06/08|
 |          for a xfem and a cutter element                             |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::storePLC(
    const DRT::Element*               cutterElement,
    const LINALG::SerialDenseMatrix&  xyze_cutterElement,
    const int                         surfId,
    vector<int>&                      positions,
    InterfacePoint&                   midpoint)   
{ 
  
  // NOTE: position is filled in the order points appear in vertices to
  // keep the order which was determined by the convex hull computation
  const unsigned numPoints = positions.size(); 
  // store segments 
  // cutter element lies on the surface of an xfem element
  if(surfId > -1)
  {
    if(numPoints == 1)
      storeIsolatedPoints(positions);
    
    if(numPoints > 1)
    {
        // store outer triangle segments
        // possible midpoint not added to position list and point list
        storeSegments(positions);
    }
    if(numPoints > 2)
    {
      // tell midpoint on which xfem surface it lies
//      classifyMidpoint(surfId, midpoint);
//      storeMidPoint(midpoint, positions);
      // store inner segments: mipoint to outer point
      storeSurfaceSegments(positions);
      // store boundary cells immediately after
//      storeSurfaceTriangles(surfId, positions);
      storeSurfaceTrianglesNoMidPoint(surfId, positions);
    }     
  }
  else if(surfId == -1)
  {
    if(numPoints > 1)
    {
        // possible midpoint not added to position list and point list
        storeSegments( positions );
    }
    if(numPoints > 2)
    {        
//      storeMidPoint(midpoint, positions);
//      storeTriangles(positions);
      storeTrianglesNoMidPoint(positions);
    } 
    
    // this method should to be called after store segments !!!!
    // so time is saved in fillPLC
    storeIsolatedPoints(positions);
  }
  else
    dserror("surface Id is not correct");
}



/*----------------------------------------------------------------------*
 |  ICS:    stores the part of a piecewise linear complex    u.may 06/08|
 |          for a xfem and a cutter element                             |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::completePLC(
    )
{
  bool removePoint = false;
  // store isolated points
  // if any isolated point appears in the segment list remove it
  for(int i = 0; i < numXFEMSurfaces_; i++)
  {   
    vector<int> removePos;
    for(unsigned int  j = 0; j < isolatedPointList_[i].size(); j++)
    {
      removePoint = false;
      const int pointPos = isolatedPointList_[i][j]; 
      for(int ii = 0; ii < numXFEMSurfaces_; ii++)
      {
        for(unsigned int k = 0; k < segmentList_[ii].size(); k++)
          if( pointPos == segmentList_[ii][k])
          {
            removePos.push_back(j);
            removePoint = true;
            break;
          }

        if(removePoint)
          break;
      }  
    }
    // count reverse in order to be able to erase the points properly
    for(int m = (int) removePos.size()-1; m >= 0; m--)
      isolatedPointList_[i].erase(isolatedPointList_[i].begin()+removePos[m]);
  }
}



/*----------------------------------------------------------------------*
 |  ICS:    finds the next facet of a convex hull            u.may 06/07|
 |          and returns the point different form the searchpoint        |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::findNextSegment(
    vector< LINALG::Matrix<3,1> >&   vertices,
    LINALG::Matrix<3,1>&             searchPoint) const
{
  bool pointfound = false;

  if(vertices.size()==0 || !searchPoint.IsInitialized())
    dserror("one or both vectors are empty");

  for(vector< LINALG::Matrix<3,1> >::iterator it = vertices.begin(); it != vertices.end(); it=it+2 )
  {
    if(comparePoints<3>(searchPoint, *it))
    {
      pointfound = true;
      searchPoint = *(it+1);
      vertices.erase(it);
      vertices.erase(it); // remove it+ 1 because it jumped on the place of the previously removed it
      break;
    }

    if(comparePoints<3>(searchPoint, *(it+1)))
    {
      pointfound = true;
      searchPoint = *(it);
      vertices.erase(it);
      vertices.erase(it); // remove it+ 1 because it jumped on the place of the previously removed it
      break;
    }
  }
  if(!pointfound) dserror("no point found");
}


/*----------------------------------------------------------------------*
 |  CDT:    checks if a CDT has to be computed               u.may 07/08|
 |          for the current xfem element                                |
 *----------------------------------------------------------------------*/
bool GEO::IntersectionExp::checkIfCDT(
  ) const
{
  bool doCDT = true;
 
  // check if interfaces points are obtained at all
  if((int) pointList_.size() < numXFEMCornerNodes_)
    doCDT = false;
  
  // triangle lists are empty means that there are no intersecting facets within
  // the xfem element
  if(doCDT)
    if(triangleList_.empty() && surfaceTriangleList_.empty() )
      doCDT = false;
  
  return doCDT;
}


/*----------------------------------------------------------------------*
 |  CDT:    computes the Constrained Delaunay                u.may 06/07|
 |          Tetrahedralization in 3D with help of Tetgen library        |
 |  for an intersected xfem element in element configuration          |
 |  TetGen provides the method                                          |
 |  void "tetrahedralize(char *switches, tetgenio* in, tetgenio* out)"  |
 |  as an interface for its use within other codes.                     |
 |  The char switches passes all the command line switches to TetGen.   |
 |  The most important command line switches include:                   |
 |      - d     detects intersections of PLC facets                     |
 |      - p     tetrahedralizes a PLC                                   |
 |      - q     quality mesh generation                                 |
 |      - nn    writes a list of boundary faces and their adjacent      |
 |              tetrahedra to the output tetgenio data structure        |
 |      - o2    resulting tetrahedra still have linear shape but        |
 |              have a 2nd order node distribution                      |
 |      - A     assigns region attributes                               |
 |      - Q     no terminal output except errors                        |
 |      - T     sets a tolerance                                        |
 |      - V     verbose: detailed information more terminal output      |
 |      - Y     prohibits Steiner point insertion on boundaries         |
 |              very help full for later visualization                  |
 |  The data structure tetgenio* in provides Tetgen with the input      |
 |  PLC and has to filled accordingly. tetgenio* out delivers the       |
 |  output information such as the resulting tetrahedral mesh.          |
 |  These two pointers must NOT be null at any time.                    |
 |  For further information please consult the TetGen manual            |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::computeCDT(
        const DRT::Element*                     xfemElement,
        const map<int,LINALG::Matrix<3,1> >&    currentcutterpositions,
        map< int, DomainIntCells >&             domainintcells,
        map< int, BoundaryIntCells >&           boundaryintcells)
{
  const int dim = 3; 
  tetgenio in;
  tetgenio out;
  std::string switchstring = "pnnQ";    //o2 Y R
  if (xfemElement->Shape() == DRT::Element::hex8)
    switchstring = "pnnQR";
  else if (xfemElement->Shape() == DRT::Element::hex20)
    switchstring = "pnnQY";

  const int switchlen = switchstring.length();
  char* switches = new char[switchlen+1];
  switchstring.copy(switches, switchlen, 0);
  switches[switchlen] = '\0';
  
  tetgenio::facet *f;
  tetgenio::polygon *p;


  const double scalefactor =  1.0e7;
  // allocate pointlist
  in.numberofpoints = pointList_.size();
  in.pointlist = new REAL[in.numberofpoints * dim];

  // fill point list
  int fill = 0;
  for(int i = 0; i <  in.numberofpoints; i++)
    for(int j = 0; j < dim; j++)
    {
      in.pointlist[fill] = (REAL) (pointList_[i].getCoord()(j) * scalefactor);
      fill++;
    }

  
  in.pointmarkerlist = new int[in.numberofpoints];
  for(int i = 0; i < numXFEMCornerNodes_; i++)
    in.pointmarkerlist[i] = 3;    // 3 : point lying on the xfem element (corner nodes)

  for(int i = numXFEMCornerNodes_; i < in.numberofpoints; i++)
    in.pointmarkerlist[i] = 2;    // 2 : point not lying on the xfem element
                                    // changed to 3 if necessary later in the loop over
                                    // xfem surfaces for simplicity

  in.numberoffacets = numXFEMSurfaces_ + triangleList_.size();

  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];

 
 
    // loop over all xfem element surfaces
  for(int i = 0; i < numXFEMSurfaces_; i++)
  {
    f = &in.facetlist[i];
    const int nsegments = (int) (segmentList_[i].size()/2);
    const int nisoPoints = isolatedPointList_[i].size();
 
    f->numberofpolygons = 1 + nsegments + nisoPoints;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    const int numnodequad4 = 4;
    p->numberofvertices = numnodequad4;
    p->vertexlist = new int[p->numberofvertices];
    for(int ivertex = 0; ivertex < numnodequad4; ivertex ++)
        p->vertexlist[ivertex] = eleNumberingSurfaces_[i][ivertex];

    // store segments
    int count = 0;
    for(int j = 1; j < 1 + nsegments; j ++)
    {
      if(segmentList_[i].size() > 0)
      {
        p = &f->polygonlist[j];
        p->numberofvertices = 2;
        p->vertexlist = new int[p->numberofvertices];

        for(int k = 0; k < 2; k++)
        {
           p->vertexlist[k] = segmentList_[i][count];
           in.pointmarkerlist[segmentList_[i][count]] = 3;  // 3: point lying on the xfem boundary
           count++;
        }
      }
    }

    // store isolated points lying on xfem surfaces
    count = 0;
    for(int j = 1 + nsegments; j < f->numberofpolygons; j++)
    {
      if(isolatedPointList_[i].size() > 0)
      {
        p = &f->polygonlist[j];
        p->numberofvertices = 1;
        p->vertexlist = new int[p->numberofvertices];

        p->vertexlist[0] = isolatedPointList_[i][count];
        in.pointmarkerlist[isolatedPointList_[i][count]] = 3;  // 3: point lying on the xfem boundary
        count++;
      }
    }
  }

  // store triangles (tri3)
  for(int i = numXFEMSurfaces_; i < in.numberoffacets; i++)
  {
    f = &in.facetlist[i];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 3;
    p->vertexlist = new int[p->numberofvertices];
    for(int j = 0; j < p->numberofvertices; j ++)
      p->vertexlist[j] = triangleList_[i - xfemElement->NumSurface()][j];
  }

  
  // set facetmarkers
  for(int i = 0; i < in.numberoffacets; i ++)
      in.facetmarkerlist[i] = faceMarker_[i] + facetMarkerOffset_;


  {
    const std::string fbase = DRT::Problem::Instance()->OutputControlFile()->FileName()
                            + ".tetin";
    const int fbaselen = fbase.length();
    char* fbasech = new char[fbaselen+1];
    fbase.copy(fbasech, fbaselen, 0);
    fbasech[fbaselen] = '\0';
    in.save_nodes(fbasech);
    in.save_poly(fbasech);
    delete [] fbasech;
  }
  //  Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
  //  do quality mesh generation (q) with a specified quality bound
  //  (1.414), and apply a maximum volume constraint (a0.1)
  // printf("tetgen start\n");
  tetrahedralize(switches, &in, &out);
  delete [] switches;
  // printf("tetgen end\n");

  //Debug
  //vector<int> elementIds;
  //for(int i = 3; i<4; i++)
  //   elementIds.push_back(xfemElement->Id());

  //debugTetgenOutput(in, out, xfemElement, elementIds);
  //printTetViewOutputPLC( element, element->Id(), in);

  // store interface triangles (+ recovery of higher order meshes)
  fill = 0; 
  for(int i = 0; i <  out.numberofpoints; i++)
    for(int j = 0; j < dim; j++)
    {
      out.pointlist[fill] = (REAL) (out.pointlist[fill] / scalefactor);
      fill++;
    }
  
  const bool higherorder = false;
  const bool recovery = false;

  if(higherorder)
  {
    std::cout << "DO RECOVERY " << endl;
    recoverCurvedInterface(xfemElement, xyze_xfemElement_, currentcutterpositions, boundaryintcells, out, recovery);
  }
  else
  {
    storeIntCells(xfemElement, xyze_xfemElement_, currentcutterpositions, boundaryintcells, out);
  }
  // store boundaryIntCells integration cells

  //if(element->Id()==388)
  //	debugFaceMarker(element->Id(), out);

  //printTetViewOutput(element->Id(), out);
  // store domain integration cells
  addCellsToDomainIntCellsMap(xfemElement, domainintcells, out, higherorder);
}



/*----------------------------------------------------------------------*
 |  CDT:    fills the point list with the corner points      u.may 06/07|
 |          in element coordinates of the xfem element                  |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::startPointList(
    )
{
  pointList_.clear();

  for(int i = 0; i < numXFEMCornerNodes_; i++)
  {
    InterfacePoint ip;
    ip.setPointType(NODE);
    ip.setNodeId(i);
    ip.setLineId(eleNodesLines_[i]);
    ip.setSurfaceId(eleNodesSurfaces_[i]);
    ip.setCoord(eleRefCoordinates_[i]);
    pointList_.push_back(ip);
  }

  for(int i = 0; i < numXFEMSurfaces_; i++)
    faceMarker_.push_back(-1);
}


/*----------------------------------------------------------------------*
 |  CDT:    fills the point list with the points             u.may 06/08|
 |          in element coordinates of the xfem element                  |
 |          for an intersection between a single cutter and xfem element|
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::storePointList(
    vector< LINALG::Matrix<3,1> >&  vertices,
    vector<int>&                    positions,
    const vector<InterfacePoint>&   interfacePoints
    )
{  
  LINALG::Matrix<3,1>              searchPoint(true);
  
  // round points an vertices with TOL3 on xfem surface
  // has to be fixed within tetgen
  //roundPointsOnXFEMBoundary(interfacePoints, GEO::TOL3);
  //roundVerticesOnXFEMBoundary(vertices, GEO::TOL3);
  
  // store interface points in pointList_
  storePoint(vertices[0], interfacePoints, positions);
  vertices.erase(vertices.begin());

  if(interfacePoints.size() > 1)
  {
      searchPoint = vertices[0];
      storePoint(vertices[0], interfacePoints, positions );
      vertices.erase(vertices.begin());
  } 
  while(vertices.size()>2)
  {
      findNextSegment(vertices, searchPoint);
      storePoint(searchPoint, interfacePoints, positions);
  }
 
  // removeDegenerateInterfacePoints(positions);
}



/*----------------------------------------------------------------------*
 |  CDT:    stores a point within a list of points           u.may 06/07|
 |          which is to be copy to the tetgen data structure            |
 |          for the computation of the                                  |
 |          Constrained Delauney Triangulation                          |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::storePoint(
    const LINALG::Matrix<3,1>&        point,
    const vector<InterfacePoint>&     interfacePoints,
    vector<int>&                      positions)
{
  bool alreadyInList = false;

  for(vector<InterfacePoint>::const_iterator ipoint = interfacePoints.begin(); ipoint != interfacePoints.end(); ++ipoint )
  {
    if(comparePoints<3>(point, ipoint->getCoord()))
    {
      alreadyInList = false;
      int count = -1;
      for(vector<InterfacePoint>::const_iterator it = pointList_.begin(); it != pointList_.end(); ++it )
      {
        count++;
        if(comparePoints<3>(point, it->getCoord()))
        {
            alreadyInList = true;
            break;
        }
      }

      if(!alreadyInList)
      {
        pointList_.push_back(*ipoint);
        positions.push_back(pointList_.size()-1);
      }
      else
      {
        positions.push_back(count);
      }
      break;
    }
  }
}


/*----------------------------------------------------------------------*
 |  CDT:    stores a point within a list of points           u.may 06/07|
 |          which is to be copy to the tetgen data structure            |
 |          for the computation of the                                  |
 |          Constrained Delauney Triangulation                          |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::storeMidPoint(
    const InterfacePoint&             midPoint,
    vector<int>&                      positions)
{
  bool alreadyInList = false;

  int count = -1;
  for(vector<InterfacePoint>::const_iterator it = pointList_.begin(); it != pointList_.end(); ++it )
  {
      count++;
      if(comparePoints<3>(midPoint.getCoord(), it->getCoord()))
      {
          alreadyInList = true;
          break;
      }
  }

  if(!alreadyInList)
  {
      pointList_.push_back(midPoint);
      positions.push_back(pointList_.size()-1);
  }
  else
  {
      positions.push_back(count);
  }
}



/*----------------------------------------------------------------------*
 |  CDT:    rounds interface points on XFEM boundary         u.may 07/08|
 |          should be used with TOL3 to make tetgen happy               |
 *----------------------------------------------------------------------*/
bool GEO::IntersectionExp::roundPointsOnXFEMBoundary(
    vector<InterfacePoint>&           interfacePoints,
    const double                      tol)
{

  bool round = false;
  switch(xfemDistype_)
  {
    case DRT::Element::hex8:  case DRT::Element::hex20: case DRT::Element::hex27:    
    {
      for(unsigned int i = 0; i < interfacePoints.size(); i++)
      {
        LINALG::Matrix<3,1> xsi(true);
        round = false;
        for(int j = 0; j < 3; j++)
        {
          xsi = interfacePoints[i].getCoord();
          if( fabs( fabs(xsi(j)) - 1.0 ) < tol && xsi(j) < 0.0)    
          {
            xsi.PutScalar(-1.0);
            round = true;
          }
          if( fabs( fabs(xsi(j)) - 1.0) < tol &&  xsi(j) > 0.0)
          {
            xsi.PutScalar(1.0);   
            round = true;
          }   
        }
        if(round)
        {  
          setBoundaryPointBoundaryStatus(xsi, interfacePoints[i]);
          interfacePoints[i].setCoord(xsi);
        }
      }
      break;
    }
    default:
    {
      cout << DRT::DistypeToString(xfemDistype_) << endl;
      dserror("discretization type %s not yet implemented", (DRT::DistypeToString(xfemDistype_)).c_str());
    }
  }
  return round;
}




/*----------------------------------------------------------------------*
 |  CDT:    rounds interface point on XFEM boundary          u.may 07/08|
 |          should be used with TOL3 to make tetgen happy               |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::roundVerticesOnXFEMBoundary(
    vector< LINALG::Matrix<3,1> >&    vertices,
    const double                      tol)
{
  switch(xfemDistype_)
  {
    case DRT::Element::hex8:  case DRT::Element::hex20: case DRT::Element::hex27:    
    {
      for(unsigned int i = 0; i < vertices.size(); i++)
        for(int j = 0; j < 3; j++)
        {
          double coord = vertices[i](j);
          if( fabs( fabs(coord) - 1.0 ) < tol && coord < 0.0)    
          {
            vertices[i](j) = -1.0;
          }
          if( fabs( fabs(coord) - 1.0) < tol &&  coord > 0.0)
          {
            vertices[i](j) =  1.0;   
          }   
        }
      
      break;
    }
    default:
    {
      cout << DRT::DistypeToString(xfemDistype_) << endl;
      dserror("discretization type %s not yet implemented", (DRT::DistypeToString(xfemDistype_)).c_str());
    }
  }
}




/*----------------------------------------------------------------------*
 |  CDT:    computes the midpoint of a collection of         u.may 06/07|
 |          InterfacePoints                                             |
 *----------------------------------------------------------------------*/
GEO::InterfacePoint GEO::IntersectionExp::computeMidpoint(
    const vector<InterfacePoint>& interfacePoints
    ) const
{ 
   InterfacePoint ip;
   LINALG::Matrix<3,1> coord(true);

   for(unsigned int i = 0; i < interfacePoints.size(); i++)
     coord += interfacePoints[i].getCoord();

   coord.Scale(1.0/((double)interfacePoints.size()));
   ip.setPointType(INTERNAL);
   ip.setCoord(coord);

   return ip;
}



/*----------------------------------------------------------------------*
 |  CDT:    computes the midpoint of a collection of         u.may 07/08|
 |          interface points determined by a position vector            |
 |          please note: point set is default to INTERNAL               |
 *----------------------------------------------------------------------*/
GEO::InterfacePoint GEO::IntersectionExp::computeMidpoint(
    const vector<int>&  positions) const
{
  InterfacePoint ip;
  LINALG::Matrix<3,1> coord(true);

  for(unsigned int i = 0; i < positions.size() ; i++)
    coord += pointList_[positions[i]].getCoord();
 
  coord.Scale(1.0/((double) positions.size()));
  ip.setPointType(INTERNAL);
  ip.setCoord(coord);

  return ip;
}



/*----------------------------------------------------------------------*
 |  CDT:    classifies the midpoint of a collection of       u.may 07/08|
 |          InterfacePoints                                             |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::classifyMidpoint(
    const int         surfId,
    InterfacePoint&   midpoint
    ) const
{
  midpoint.setPointType(SURFACE);
  vector<int> surfaces(1, surfId);
  midpoint.setSurfaceId(surfaces);
}



/*----------------------------------------------------------------------*
 |  CDT:    stores a isolated point lying on a surface of an u.may 07/08|
 |          xfem element, if it is not a segment point                  |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::storeIsolatedPoints(
    const vector<int>&                positions)
{

  bool noIsolatedPoint = false;

  for(unsigned int i = 0; i < positions.size(); i++)
  {  
    // NODE type interface points don't need to be stored, these points are the corner
    // points
    noIsolatedPoint = false;
    const InterfacePoint ip = pointList_[positions[i]];

    if(ip.getPointType() == LINE || ip.getPointType() == SURFACE)
    {
      int countEnd = 0;
      for(int j = 0; j < ip.getNumSurface(); j++)
      {
        // check if point poistion is already stored in the segment list of this surface
        const int surf_j = ip.getSurfId()[j];
        if(segmentList_[surf_j].size() > 0)
        {
          vector<int>::iterator itSegment = find ( segmentList_[surf_j].begin(), 
              segmentList_[surf_j].end(), positions[i]);

          if( (itSegment) != segmentList_[surf_j].end())
            noIsolatedPoint = true;
        }

        if(noIsolatedPoint)
          break;

        if(isolatedPointList_[surf_j].size() > 0)
        {
          // check if point poistion is already stored in isolated point list
          vector<int>::iterator itPoint = find (  isolatedPointList_[surf_j].begin(), 
              isolatedPointList_[surf_j].end(), 
              positions[i]);

          if(itPoint == isolatedPointList_[surf_j].end())
            countEnd++;
        }
        else
          countEnd++;

      }
      // store only on one surface even if the point lies on a line
      if(countEnd == ip.getNumSurface())
        isolatedPointList_[ip.getSurfId()[0]].push_back(positions[i]);
    } 
  }
}



/*----------------------------------------------------------------------*
 |  CDT:    stores a single segment within a list of         u.may 07/08|
 |          segments, which is to be copied to the tetgen data structure|
 |          for the computation of the                                  |
 |          Constrained Delauney Triangulation                          |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::storeSingleSegment(
    const int pos1,
    const int pos2)
{
  if(!checkIfSegmentPointsOnSameXfemLine(pos1, pos2) && pos1 != pos2)
  {   
    // loops over all NODE LINE and SURFACE type points
    for(int j = 0; j < pointList_[pos1].getNumSurface(); j++ )
    {
      for(int k = 0; k < pointList_[pos2].getNumSurface(); k++ )
      {
        const int surf1 = pointList_[pos1].getSurfId()[j];
        const int surf2 = pointList_[pos2].getSurfId()[k];

        if( (surf1 == surf2) )
        {
          bool alreadyInList = false;

          for(unsigned int is = 0 ; is < segmentList_[surf1].size() ; is = is + 2)
          {
            if( (segmentList_[surf1][is] == pos1  &&  segmentList_[surf1][is+1] == pos2)  ||
                (segmentList_[surf1][is] == pos2  &&  segmentList_[surf1][is+1] == pos1) )
            {
              alreadyInList = true;
              break;
            }
          }

          if(!alreadyInList)
          {
            segmentList_[surf1].push_back(pos1);
            segmentList_[surf1].push_back(pos2);
          }
        }
      }
    }
  }
}



/*----------------------------------------------------------------------*
 |  CDT:    stores segments within a list of segments       u.may 06/07|
 |          which is to be copied to the tetgen data structure          |
 |          for the computation of the                                  |
 |          Constrained Delauney Triangulation                          |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::storeSegments(
    const vector<int>&              positions)
{
  // midpoint is not yet addded to position list !!!
  for(unsigned int i = 0; i < positions.size(); i++ )
  {
    const int pos1 = positions[i];
    int pos2 = 0;
    if(pos1 ==  positions[positions.size()-1])
      pos2 = positions[0];
    else
      pos2 = positions[i+1];

    // if both point are lying on the same surface but not on the same line store segment
    // if not already stored
    storeSingleSegment(pos1, pos2);
  }
}



/*----------------------------------------------------------------------*
 |  CDT:    stores segments within a list of segments       u.may 06/07|
 |          which is to be copied to the tetgen data structure          |
 |          for the computation of the                                  |
 |          Constrained Delauney Triangulation                          |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::storeSurfaceSegments(
    const vector<int>&              positions)
{

  // If more than two points resulting from an intersection of a single
  // cutter element and a single xfem element require special treatment
  // to keep the tetgen data structure as simple as possible
  // In this case the cutter element lies partially on one xfem surface.
  // The mid point is computed and the resulting triangles are stored as
  // segments.

  // last entry in position vector corresponds to the midpoint
  // store outer segments 
  const int pos1 = positions[positions.size()-1];  // midpoint
  for(unsigned int i = 0; i < positions.size()-1; i++ )
  {
    const int pos2 = positions[i];                   // all other points 

    // if both point are lying on the same surface but not on the same line store segment
    // if not already stored
    storeSingleSegment(pos1, pos2);
  }
}



/*----------------------------------------------------------------------*
 |  CDT:    check if two segment points are on the same line u.may 07/08|
 *----------------------------------------------------------------------*/
bool GEO::IntersectionExp::checkIfSegmentPointsOnSameXfemLine(
    const int position1,
    const int position2)
{
  bool onSameLine = false;

  if((pointList_[position1].getPointType() == NODE || pointList_[position1].getPointType() == LINE) && 
     (pointList_[position2].getPointType() == NODE || pointList_[position2].getPointType() == LINE))
  { 
    for(int i = 0; i < pointList_[position1].getNumLine(); i++)
    {  
      const int lineId = pointList_[position1].getLineId()[i];
      for(int j = 0; j < pointList_[position2].getNumLine(); j++)
      {
        if(lineId == pointList_[position2].getLineId()[j])
        {
          onSameLine = true;
          break;
        }
      }
      if(onSameLine)
        break;
    }
  }

  return onSameLine;
}


/*----------------------------------------------------------------------*
 |  CDT:    removes equal interface points from              u.may 07/08|
 |          positions                                                   |  
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::removeDegenerateInterfacePoints(
    vector<int>&              positions)
{
  vector<int> removePoints;

  for(int i = 0; i < (int) positions.size()-1; i++)
    for(unsigned int j = i+1; j < positions.size(); j++)
      if(positions[i] == positions[j])
      {
        removePoints.push_back(i);
        break;
      }

  for(int i = removePoints.size()-1; i >= 0; i--)
    positions.erase(positions.begin()+removePoints[i]);

}


/*----------------------------------------------------------------------*
 |  CDT:    stores a triangle within a list of trianles      u.may 06/07|
 |          which is to be copy to the tetgen data structure            |
 |          for the computation of the                                  |
 |          Constrained Delauney Triangulation                          |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::storeTriangles(
    const vector<int>&               positions)
{
  vector<int> triangle(3,0);

  // store midpoint before on last entry in positions and point list
  for(unsigned int i = 0; i < positions.size()-2; i++ )
  {
    triangle[0] = positions[i];
    triangle[1] = positions[i+1];
    triangle[2] = positions[positions.size()-1];      // add midpoint

    triangleList_.push_back(triangle);
    faceMarker_.push_back(intersectingCutterElements_.size()-1);
  }

  // last point and first point
  triangle[0] = positions[positions.size()-2];
  triangle[1] = positions[0];
  triangle[2] = positions[positions.size()-1];          // add midpoint

  triangleList_.push_back(triangle);
  faceMarker_.push_back(intersectingCutterElements_.size()-1);
}


/*----------------------------------------------------------------------*
 |  CDT:    stores a triangle within a list of trianles      u.may 06/07|
 |          which is to be copy to the tetgen data structure            |
 |          for the computation of the                                  |
 |          Constrained Delauney Triangulation                          |
 |          no midpoint was added to positions                          |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::storeTrianglesNoMidPoint(
    const vector<int>&               positions)
{
  vector<int> triangle(3,0);

  // store midpoint before on last entry in positions and point list
  for(unsigned int i = 0; i < positions.size()-2; i++ )
  {
    triangle[0] = positions[0];
    triangle[1] = positions[i+1];
    triangle[2] = positions[i+2];

    triangleList_.push_back(triangle);
    faceMarker_.push_back(intersectingCutterElements_.size()-1);
  }
}



/*----------------------------------------------------------------------*
 |  CDT:    stores a triangle within a list of triangles     u.may 06/08|
 |          only lying on the xfem surfaces                             |
 |          which is to be copy to the tetgen data structure            |
 |          for the computation of the                                  |
 |          Constrained Delauney Triangulation                          |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::storeSurfaceTriangles(
    const int                       surfId,
    const vector<int>               positions)
{
    vector<int> triangle(3,0);
    vector<vector<int> > triangleList;
    
    // store midpoint before in last entry in positions and point list
    for(unsigned int i = 0; i < positions.size()-2; i++ )
    {
        triangle[0] = positions[i];
        triangle[1] = positions[i+1];
        triangle[2] = positions[positions.size()-1];      // add midpoint

        triangleList.push_back(triangle); 
    }

    // last point to first point
    triangle[0] = positions[positions.size()-2];
    triangle[1] = positions[0];
    triangle[2] = positions[positions.size()-1];      // add midpoint

    triangleList.push_back(triangle);
    
    // current element is always the last element in vector
    const int cutterPosition = intersectingCutterElements_.size()-1;
    surfaceTriangleList_.insert(make_pair(cutterPosition,triangleList));
}


/*----------------------------------------------------------------------*
 |  CDT:    stores a triangle within a list of triangles     u.may 06/08|
 |          only lying on the xfem surfaces                             |
 |          which is to be copy to the tetgen data structure            |
 |          for the computation of the                                  |
 |          Constrained Delauney Triangulation                          |
 |          no midpoint was added to positions                          |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::storeSurfaceTrianglesNoMidPoint(
    const int                       surfId,
    const vector<int>               positions)
{
    vector<int> triangle(3,0);
    vector<vector<int> > triangleList;
    
    for(unsigned int i = 0; i < positions.size()-2; i++ )
    {
        triangle[0] = positions[0];
        triangle[1] = positions[i+1];
        triangle[2] = positions[i+2];

        triangleList.push_back(triangle); 
    }
    
    // current element is always the last element in vector
    const int cutterPosition = intersectingCutterElements_.size()-1;
    surfaceTriangleList_.insert(make_pair(cutterPosition,triangleList));
}


/*----------------------------------------------------------------------*
 |  RCI:    stores a pointer to each intersecting            u.may 08/07|
 |          cutter element  used for the recovery of curved interface   |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::storeIntersectedCutterElement(
        DRT::Element*  cutterElement)
{
    bool alreadyInList = false;

    vector<DRT::Element*>::const_iterator eleiter;
    for(eleiter = intersectingCutterElements_.begin(); eleiter != intersectingCutterElements_.end(); ++eleiter)
        if(*eleiter == cutterElement)
        {
            alreadyInList = true;
            break;
        }

    if(!alreadyInList)
        intersectingCutterElements_.push_back(cutterElement);
}



/*----------------------------------------------------------------------*
 |  RCI:    recovers the curved interface after the          u.may 08/07|
 |          Constrained Delaunay Tetrahedralization                     |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::recoverCurvedInterface(
        const DRT::Element*                     xfemElement,
        const LINALG::SerialDenseMatrix&        xyze_xfemElement,
        const map<int,LINALG::Matrix<3,1> >&    currentcutterpositions,
        map< int, BoundaryIntCells >&           boundaryintcells,
        tetgenio&                               out,
        bool							                      recovery
        )
{
    BoundaryIntCells                     			listBoundaryICPerElement;

    // list of point markers , if already visited = 1 , if not = 0
    vector<int> visitedPointIndexList(out.numberofpoints,0);

    // store triface lying completey on xfem surfaces
    if(!surfaceTriangleList_.empty())
        storeSurfaceIntCells(true,xfemElement,xyze_xfemElement, currentcutterpositions,listBoundaryICPerElement);
    
    // lifts all corner points into the curved interface
    if(recovery)
    	liftAllSteinerPoints(xfemElement, xyze_xfemElement, currentcutterpositions, out);

    for(int i=0; i<out.numberoftrifaces; i++)
    {
        // run over all faces not lying in on of the xfem element planes
        const int faceMarker = out.trifacemarkerlist[i] - facetMarkerOffset_;
        vector<vector<double> > domainCoord(6, std::vector<double>(3,0.0));
        vector<vector<double> > boundaryCoord(6, std::vector<double>(3,0.0));

        if(faceMarker > -1)
        {
            const int tetIndex = out.adjtetlist[i*2];
            //printf("tetIndex = %d\n", tetIndex);
            vector<int>             order(3,0);
            vector<int>             tetraCornerIndices(4,0);
            vector<LINALG::Matrix<3,1> >        tetraCornerNodes(4, LINALG::Matrix<3,1>(true));
            getTetrahedronInformation(tetIndex, i, tetraCornerIndices, order, out );
            getTetrahedronNodes(tetraCornerNodes, tetraCornerIndices, xfemElement, xyze_xfemElement, out);

            // run over each triface
            for(int index1 = 0; index1 < 3 ;index1++)
            {
                int index2 = index1+1;
                if(index2 > 2) index2 = 0;

                //printf("edgeIndex1 = %d\n", out.tetrahedronlist[tetIndex*out.numberofcorners+order[index1]]);
                //printf("edgeIndex2 = %d\n", out.tetrahedronlist[tetIndex*out.numberofcorners+order[index2]]);

                const int localHigherOrderIndex = DRT::UTILS::getHigherOrderIndex(order[index1], order[index2], DRT::Element::tet10);
                //printf("localHOindex = %d\n", localHigherOrderIndex);
                const int globalHigherOrderIndex = out.tetrahedronlist[tetIndex*out.numberofcorners+localHigherOrderIndex];
                //printf("globalHOindex = %d\n", globalHigherOrderIndex);
                if(visitedPointIndexList[globalHigherOrderIndex]== 0 && recovery)
                {
                    visitedPointIndexList[globalHigherOrderIndex] = 1;

                    computeHigherOrderPoint(    index1, index2, i, faceMarker, globalHigherOrderIndex,
                                                tetraCornerIndices, tetraCornerNodes, xfemElement, xyze_xfemElement, currentcutterpositions, out);
                }

                // store boundary integration cells
                addCellsToBoundaryIntCellsMap(	i, index1, globalHigherOrderIndex, faceMarker, currentcutterpositions,
                						                    domainCoord, boundaryCoord, xfemElement, xyze_xfemElement, out);
            }

            /*for(int ii=0; ii < 6; ii++)
            	for(int jj=0; jj < 3; jj++)
            		printf("boundary = %f\n",boundaryCoord[ii][jj]);
            	*/
            const int ele_gid = intersectingCutterElements_[faceMarker]->Id();
            listBoundaryICPerElement.push_back(BoundaryIntCell(DRT::Element::tri6, ele_gid, domainCoord, boundaryCoord));

        }

    }

    boundaryintcells.insert(make_pair(xfemElement->Id(),listBoundaryICPerElement));

    intersectingCutterElements_.clear();
}



/*----------------------------------------------------------------------*
 |  RCI:    store linear boundary and integration cells      u.may 03/08|
 |                               										                    |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::storeIntCells(
    const DRT::Element*                   xfemElement,
    const LINALG::SerialDenseMatrix&      xyze_xfemElement,
    const map<int,LINALG::Matrix<3,1> >&  currentcutterpositions,
    map< int, BoundaryIntCells >&         boundaryintcells,
    tetgenio&                             out)
{

  BoundaryIntCells  listBoundaryICPerElement;

  // store cells completey lying on xfem boundaries
  // no lifting necessayr if -Y switch is applied and/or volume element is Cartesian
  if(!surfaceTriangleList_.empty())
    storeSurfaceIntCells(false, xfemElement, xyze_xfemElement, currentcutterpositions,listBoundaryICPerElement);

  // lifts all corner points into the curved interface
  //liftAllSteinerPoints(xfemElement, xyze_xfemElement, currentcutterpositions, out);

  for(int i=0; i<out.numberoftrifaces; i++)
  {
    // run over all faces not lying in on of the xfem element planes
    const int faceMarker = out.trifacemarkerlist[i] - facetMarkerOffset_;
    vector<vector<double> > domainCoord(3, std::vector<double>(3,0.0));
    vector<vector<double> > boundaryCoord(3, std::vector<double>(3,0.0));

    if(faceMarker > -1)
    {
      // run over each triface
      for(int index1 = 0; index1 < 3 ;index1++)
      {
          // store boundary integration cells
          const int globalHigherOrderIndex = -1; // means no quadratic points (tri3 instead of tri6)
          addCellsToBoundaryIntCellsMap(	i, index1, globalHigherOrderIndex, faceMarker, currentcutterpositions,
          								domainCoord, boundaryCoord, xfemElement, xyze_xfemElement, out);
      }
      const int ele_gid = intersectingCutterElements_[faceMarker]->Id();
      listBoundaryICPerElement.push_back(BoundaryIntCell(DRT::Element::tri3, ele_gid, domainCoord, boundaryCoord));
    }
  }
  boundaryintcells.insert(make_pair(xfemElement->Id(),listBoundaryICPerElement));
  intersectingCutterElements_.clear();
}



/*----------------------------------------------------------------------*
 |  RCI:    store linear boundary and integration cells      u.may 07/08|
 |          lying on xfem surfaces                                      |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::storeSurfaceIntCells(
    const bool                            higherorder,
    const DRT::Element*                   xfemElement,
    const LINALG::SerialDenseMatrix&      xyze_xfemElement,
    const map<int,LINALG::Matrix<3,1> >&  currentcutterpositions,
    BoundaryIntCells&                     listBoundaryICPerElement)
{
  vector<vector<double> > domainCoord(3, std::vector<double>(3,0.0));
  vector<vector<double> > boundaryCoord(3, std::vector<double>(3,0.0));

  for(map< int, vector<vector<int> > >::iterator mapit = surfaceTriangleList_.begin(); mapit != surfaceTriangleList_.end(); ++mapit )
  {
    const int cutterPosition = mapit->first;
    const int ele_gid = intersectingCutterElements_[cutterPosition]->Id();
    for(vector<vector<int> >::iterator vecit = (mapit->second).begin(); vecit != (mapit->second).end(); ++vecit )
    {
      const vector<int> triface = (*vecit);
      for(int index1 = 0; index1 < 3; index1++)
      {
        int index2 = index1+1;
        if(index2 > 2) index2 = 0;

        addXFEMSurfaceCellsToBoundaryIntCellsMap( higherorder, index1, triface[index1], triface[index2], cutterPosition, 
                                                  currentcutterpositions, domainCoord, boundaryCoord, 
                                                  xfemElement, xyze_xfemElement);
      }
      listBoundaryICPerElement.push_back(BoundaryIntCell(DRT::Element::tri3, ele_gid, domainCoord, boundaryCoord));
    }
  }
}


/*----------------------------------------------------------------------*
 |  RCI:    checks if all tetrahedra corner points are lying u.may 09/07|
 |          in a surface element                                        |
 |          if not corner points is recovered on the surface element    |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::liftAllSteinerPoints(
        const DRT::Element*                     xfemElement,
        const LINALG::SerialDenseMatrix&        xyze_xfemElement,
        const map<int,LINALG::Matrix<3,1> >&    currentcutterpositions,
        tetgenio&                               out)
{
    LINALG::Matrix<3,1> edgePoint(true);     
    LINALG::Matrix<3,1> oppositePoint(true);  
    vector< vector<int> > adjacentFacesList;
    vector< vector<int> > adjacentFacemarkerList;
    
    locateSteinerPoints(adjacentFacesList, adjacentFacemarkerList, out);
   
    if(!adjacentFacesList.empty())
    {
        // run over all Steiner points
        for(unsigned int i=0; i<adjacentFacesList.size(); i++)
        {
            int lineIndex = -1;
            int cutterIndex = -1;
            const SteinerType caseSteiner = decideSteinerCase(  i, lineIndex, cutterIndex, adjacentFacesList, adjacentFacemarkerList, currentcutterpositions,
                                                         edgePoint, oppositePoint, xfemElement, xyze_xfemElement, out);
           
            switch(caseSteiner)
            {
                case S_SURFACE:
                {
                    liftSteinerPointOnSurface(i, adjacentFacesList, adjacentFacemarkerList, currentcutterpositions, xfemElement, xyze_xfemElement, out);
                    break;
                }
                case S_EDGE:
                {   
                    liftSteinerPointOnEdge( i, lineIndex, cutterIndex, edgePoint, oppositePoint,
                                            adjacentFacesList, currentcutterpositions, xfemElement, xyze_xfemElement, out);
                    break;
                }
                case S_BOUNDARY:
                {
                    liftSteinerPointOnBoundary( i, adjacentFacesList, adjacentFacemarkerList, currentcutterpositions, xfemElement, xyze_xfemElement, out);
                    break;
                }
                default:
                    dserror("case of lifting Steiner point does not exist");
            }
        }
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    stores for each Steiner point its adjacent       u.may 11/07|
 |          faces and face markers                                      |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::locateSteinerPoints(
    vector< vector<int> >&     adjacentFacesList,
    vector< vector<int> >&     adjacentFacemarkerList,
    const tetgenio&             out
    ) const
{

    for(int i = 0; i < out.numberoftrifaces; i++)
    {
        if( (out.trifacemarkerlist[i]-facetMarkerOffset_) > -1)
        {
            for (int j = 0; j < 3; j++)
            {
                const int pointIndex = out.trifacelist[i*3+j];

                // check if point is a Steiner point
                if(out.pointmarkerlist[pointIndex] != 2 && out.pointmarkerlist[pointIndex] != 3 )
                {
                    bool alreadyInList = false;
                    const vector<int> pointIndices = getPointIndices(out, i, j);

                    for(unsigned int k = 0; k < adjacentFacesList.size(); k++)
                        if(adjacentFacesList[k][0] == pointIndex)
                        {
                            alreadyInList = true;

                            adjacentFacesList[k].push_back(pointIndices[0]);
                            adjacentFacesList[k].push_back(pointIndices[1]);
                            adjacentFacemarkerList[k].push_back( out.trifacemarkerlist[i]-facetMarkerOffset_ );
                            break;
                        }


                    if(!alreadyInList)
                    {
                        vector<int> adjacentFaces(3);
                        adjacentFaces[0] = pointIndex;
                        adjacentFaces[1] = pointIndices[0];
                        adjacentFaces[2] = pointIndices[1];

                        adjacentFacesList.push_back(adjacentFaces);

                        vector<int> adjacentFacemarkers(1);
                        adjacentFacemarkers[0] = out.trifacemarkerlist[i]-facetMarkerOffset_ ;
                        adjacentFacemarkerList.push_back(adjacentFacemarkers);
                    }
                }
            }
        }
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    checks if the Steiner points lies within         u.may 11/07|
 |          the cutter element or on one of its edges                   |
 *----------------------------------------------------------------------*/
GEO::SteinerType GEO::IntersectionExp::decideSteinerCase(
        const int                             steinerIndex,
        int&                                  lineIndex,
        int&                                  cutterIndex,
        const vector< vector<int> >&          adjacentFacesList,
        const vector< vector<int> >&          adjacentFacemarkerList,
        const map<int,LINALG::Matrix<3,1> >&  currentcutterpositions,
        LINALG::Matrix<3,1>&                  edgePoint,
        LINALG::Matrix<3,1>&                  oppositePoint,
        const DRT::Element*                   xfemElement,
        const LINALG::SerialDenseMatrix&      xyze_xfemElement,
        const tetgenio&                       out
        ) const
{
    const int pointIndex = adjacentFacesList[steinerIndex][0];

    static LINALG::Matrix<3,1> x;
    for(int k=0; k<3; k++)
        x(k)   = out.pointlist[pointIndex*3 + k];

    static LINALG::Matrix<3,1>    xsi;
    // check exact TODO
    xsi = currentToVolumeElementCoordinatesExact(xfemElement->Shape(), xyze_xfemElement_, x, TOL7);
    //currentToVolumeElementCoordinates(xfemElement, x, xsi);

    InterfacePoint emptyIp;
    if(setInternalPointBoundaryStatus(xsi, emptyIp))
        out.pointmarkerlist[pointIndex] = 3;    // on xfem boundary
    else
        out.pointmarkerlist[pointIndex] = 2;    // not on xfem boundary

    bool normalSteiner = true;
    for(unsigned int j = 0; j < adjacentFacemarkerList[steinerIndex].size(); j++ )
    {
        for(unsigned int k = 0; k < adjacentFacemarkerList[steinerIndex].size(); k++ )
        {
            if(adjacentFacemarkerList[steinerIndex][j] != adjacentFacemarkerList[steinerIndex][k])
            {
                //printf("a = %d and b = %d\n", adjacentFacemarkerList[steinerIndex][j], adjacentFacemarkerList[steinerIndex][k]);
                if(findCommonFaceEdge(  j, k, adjacentFacesList[steinerIndex], edgePoint, oppositePoint, out))
                {
                    if(!findCommonCutterLine(  currentcutterpositions, adjacentFacemarkerList[steinerIndex][j], adjacentFacemarkerList[steinerIndex][k],
                                               lineIndex, cutterIndex))
                                               dserror("no common line element found\n");

                    normalSteiner = false;
                }
            }
            if(!normalSteiner)      break;
        }
        if(!normalSteiner)      break;
    }

    SteinerType caseSteiner = S_SURFACE;
    if(!normalSteiner)
        caseSteiner = S_EDGE;
    if(out.pointmarkerlist[pointIndex] == 3)
        caseSteiner = S_BOUNDARY;

    return caseSteiner;
}



/*----------------------------------------------------------------------*
 |  RCI:    lifts Steiner points lying                       u.may 11/07|
 |          within a cutter element                                     |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::liftSteinerPointOnSurface(
        const int                             steinerIndex,
        const vector< vector<int> >&          adjacentFacesList,
        const vector< vector<int> >&          adjacentFacemarkerList,
        const map<int,LINALG::Matrix<3,1> >&  currentcutterpositions,
        const DRT::Element*                   xfemElement,
        const LINALG::SerialDenseMatrix&      xyze_xfemElement,
        tetgenio&                             out
        )
{
    // get Steiner point coordinates
    LINALG::Matrix<3,1>   Steinerpoint;
    for(int j=0; j<3; ++j)
    {
        Steinerpoint(j) = out.pointlist[adjacentFacesList[steinerIndex][0]*3 + j];
    }
    elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement, Steinerpoint);

    LINALG::Matrix<3,1>  averageNormal(true);

    const int length = (int) ( ( (double) (adjacentFacesList[steinerIndex].size()-1))*0.5 );
    vector< LINALG::Matrix<3,1> > normals;
    normals.reserve(length);

    for(int j = 0; j < length; ++j)
    {
        const int pointIndex1 = adjacentFacesList[steinerIndex][1 + 2*j];
        const int pointIndex2 = adjacentFacesList[steinerIndex][1 + 2*j + 1];

        LINALG::Matrix<3,1>    p1;
        LINALG::Matrix<3,1>    p2;
        for(int k = 0; k < 3; ++k)
        {
            p1(k) =  out.pointlist[pointIndex1*3 + k];
            p2(k) =  out.pointlist[pointIndex2*3 + k];
        }
        elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement, p1);
        elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement, p2);

        LINALG::Matrix<3,1> n1(true); 
        n1.Update(1.0, p1, -1.0, Steinerpoint);
        LINALG::Matrix<3,1> n2(true); 
        n2.Update(1.0, p2, -1.0, Steinerpoint);

        LINALG::Matrix<3,1> normal = computeCrossProduct(n1, n2);
        //normalizeVectorInPLace(normal);
        normal.Scale(1.0/normal.Norm2());
        averageNormal += normal;

        normals.push_back(normal);
    }

    averageNormal.Scale(1.0/((double)length));

    const int faceMarker = adjacentFacemarkerList[steinerIndex][0];
    DRT::Element* cutterElement =  intersectingCutterElements_[faceMarker];
    const LINALG::SerialDenseMatrix xyze_cutterElement(GEO::getCurrentNodalPositions(cutterElement, currentcutterpositions));

    LINALG::Matrix<3,1> xsi(true);
    vector< LINALG::Matrix<3,1> > plane;
    LINALG::Matrix<3,1> normalNode(true);
    normalNode.Update(1.0, Steinerpoint, 1.0, averageNormal);
    plane.push_back(normalNode);
    normalNode.Update(1.0, Steinerpoint, -1.0, averageNormal);
    plane.push_back(normalNode);
    
    bool intersected = computeRecoveryNormal( xsi, plane, cutterElement, xyze_cutterElement, false);
    if(intersected)
    {
        storeHigherOrderNode(   true, adjacentFacesList[steinerIndex][0], -1, xsi,
            cutterElement, currentcutterpositions, xfemElement, out);
    }
    else
    {
      intersected = false;
      // loop over all individual normals
      vector< LINALG::Matrix<3,1> >::const_iterator normalptr;
      for(normalptr = normals.begin(); normalptr != normals.end(); ++normalptr )
      {
        vector< LINALG::Matrix<3,1> > plane;
        normalNode.Update(1.0, Steinerpoint, 1.0, (*normalptr));
        plane.push_back(normalNode);
        normalNode.Update(1.0, Steinerpoint, -1.0, (*normalptr));
        plane.push_back(normalNode);
        intersected = computeRecoveryNormal( xsi, plane, cutterElement, xyze_cutterElement, false);
        if(intersected)
        {
          storeHigherOrderNode( true, adjacentFacesList[steinerIndex][0], -1, xsi,
              cutterElement, currentcutterpositions, xfemElement, out);
          break;
        }
      }
      if(!intersected)
      {
          countMissedPoints_++;
          //printf("STEINER POINT NOT LIFTED in liftSteinerPointOnSurface()!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      }
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    lifts Steiner points lying                       u.may 11/07|
 |          on the edge of a  cutter element                            |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::liftSteinerPointOnEdge(
    const int                             steinerIndex,
    int                                   lineIndex,     // no const because line index is changed in computeRecoveryPlane !!
    const int                             cutterIndex,
    LINALG::Matrix<3,1>&                  edgePoint,
    LINALG::Matrix<3,1>&                  oppositePoint,
    const vector< vector<int> >&          adjacentFacesList,
    const map<int,LINALG::Matrix<3,1> >&  currentcutterpositions,
    const DRT::Element*                   xfemElement,
    const LINALG::SerialDenseMatrix&      xyze_xfemElement,
    tetgenio&                             out)
{
    // get Steiner point coordinates
    LINALG::Matrix<3,1>    Steinerpoint(true);
    for(int j=0; j<3; j++)
        Steinerpoint(j) = out.pointlist[adjacentFacesList[steinerIndex][0]*3 + j];

    elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement, Steinerpoint);
    elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement, edgePoint);
    elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement, oppositePoint);

    LINALG::Matrix<3,1> r1; r1.Update(1.0, edgePoint, -1.0, Steinerpoint);
    LINALG::Matrix<3,1> r2; r2.Update(1.0, oppositePoint, -1.0, Steinerpoint);
    
    LINALG::Matrix<3,1> n1 = computeCrossProduct( r1, r2);
    LINALG::Matrix<3,1> n2 = computeCrossProduct( r1, n1);

    //normalizeVectorInPLace(n1);
    //normalizeVectorInPLace(n2);  
    n1.Scale(1.0/n1.Norm2());
    n2.Scale(1.0/n2.Norm2());

    vector< LINALG::Matrix<3,1> > plane;
    LINALG::Matrix<3,1> planeNode(true);
    planeNode.Update(1.0, Steinerpoint, 1.0, n1);
    plane.push_back(planeNode);
    
    planeNode.Update(1.0,Steinerpoint, -1.0, n1);
    plane.push_back(planeNode);
    
    planeNode.Update(1.0,plane[1], 1.0, n2);
    plane.push_back(planeNode);
    
    planeNode.Update(1.0,plane[0], 1.0, n2);
    plane.push_back(planeNode);
    
    static LINALG::Matrix<3,1> xsi;
    xsi = 0.0;
    DRT::Element* cutterElement = intersectingCutterElements_[cutterIndex];
    const LINALG::SerialDenseMatrix xyze_cutterElement(GEO::getCurrentNodalPositions(cutterElement, currentcutterpositions));  
    const bool intersected = computeRecoveryPlane( lineIndex, currentcutterpositions, xsi, plane, cutterElement, xyze_cutterElement);
   
    if(intersected)
    {
        storeHigherOrderNode( false, adjacentFacesList[steinerIndex][0], lineIndex, xsi,
                              cutterElement, currentcutterpositions, xfemElement, out);
    }
    else
    {
        countMissedPoints_++;
        //printf("STEINER POINT NOT LIFTED in liftSteinerPointOnEdge()!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    }
}


/*----------------------------------------------------------------------*
 |  RCI:    lifts Steiner points lying                       u.may 11/07|
 |          on the boundary of the xfem element                         |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::liftSteinerPointOnBoundary(
        const int                             steinerIndex,
        const vector< vector<int> >&          adjacentFacesList,
        const vector< vector<int> >&          adjacentFacemarkerList,
        const map<int,LINALG::Matrix<3,1> >&  currentcutterpositions,
        const DRT::Element*                   xfemElement,
        const LINALG::SerialDenseMatrix&      xyze_xfemElement,
        tetgenio&                             out
        )
{
    int edgeIndex = 0;
    int oppositeIndex = 0;
    int facemarkerIndex = 0;

    // find egde on boundary
    for(unsigned int i = 1; i < adjacentFacesList[steinerIndex].size(); i++ )
    {
        edgeIndex = adjacentFacesList[steinerIndex][i];
        if(out.pointmarkerlist[edgeIndex] == 3)
        {
            facemarkerIndex = (int) (i+1)/2;
            break;
        }
    }

    int faceIndex = adjacentFacemarkerList[steinerIndex][facemarkerIndex];
    // locate triangle on boundary
    bool oppositeFound = false;
    for(int i = 0; i < out.numberoftrifaces; i++)
    {
        if( (out.trifacemarkerlist[i]-facetMarkerOffset_) == -1)
        {
            int countIndex = 0;
            for(int j = 0; j < 3; j++)
            {
                const int index = out.trifacelist[i*3+j];
                if(index == steinerIndex || index == edgeIndex)
                    countIndex++;
            }

            if(countIndex == 2)
            {
                for(int j = 0; j < 3; j++)
                {
                    const int index = out.trifacelist[i*3+j];
                    if(index != steinerIndex && index != edgeIndex)
                    {
                        oppositeIndex = index;
                        oppositeFound = true;
                        break;
                    }
                }
            }
        }
        if(!oppositeFound)
            break;
    }
    
    
    // compute normal through Steiner point lying in the boundary triangle
    vector<LINALG::Matrix<3,1> > plane;
    computeIntersectionNormalC( adjacentFacesList[steinerIndex][0], edgeIndex, oppositeIndex, plane,
                               xfemElement, out);

   
    // compute intersection normal on boundary
    static LINALG::Matrix<3,1> xsi;
    DRT::Element* cutterElement = intersectingCutterElements_[faceIndex];
    const LINALG::SerialDenseMatrix xyze_cutterElement(GEO::getCurrentNodalPositions(cutterElement, currentcutterpositions));   
    const bool intersected = computeRecoveryNormal(xsi, plane, cutterElement, xyze_cutterElement, true);
   
    if(intersected)
    {
        storeHigherOrderNode(   true, adjacentFacesList[steinerIndex][0], -1, xsi,
                                intersectingCutterElements_[faceIndex], currentcutterpositions, xfemElement, out);
    }
    else
    {
        countMissedPoints_++;
        //printf("STEINER POINT NOT LIFTED in liftSteinerPointOnBoundary()!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    }

}



/*----------------------------------------------------------------------*
 |  RCI:    returns information of the tetrahedra            u.may 08/07|
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::getTetrahedronInformation(
    const int           tetIndex,
    const int           faceIndex,
    vector<int>&        tetraCornerIndices,
    vector<int>&        order,
    const tetgenio&     out) const
{

    // store boundary face node indices
    for(int j=0; j<3; j++)
        tetraCornerIndices[j] = out.trifacelist[faceIndex*3+j];

    // store node index opposite to the boundary face of the tetrahedron
    for(int j=0; j<4; j++)
    {
        const int nodeIndex = out.tetrahedronlist[tetIndex*out.numberofcorners+j];
        if(nodeIndex != tetraCornerIndices[0] && nodeIndex != tetraCornerIndices[1] && nodeIndex != tetraCornerIndices[2])
        {
            tetraCornerIndices[3] = out.tetrahedronlist[tetIndex*out.numberofcorners+j];
            break;
        }
    }

    for(int j=0; j<4; j++)
    {
        const int nodeIndex = out.tetrahedronlist[tetIndex*out.numberofcorners+j];
        for(int k=0; k<3; k++)
            if(nodeIndex == tetraCornerIndices[k])
            {
                order[k] = j;
                break;
            }
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    collects the tetrahedron corner nodes            u.may 09/07|
 |          transforms them into current coordinates                    |
 |          of the xfem-element                                         |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::getTetrahedronNodes(
        vector<LINALG::Matrix<3,1> >&           tetraCornerNodes,
        const vector<int>&                      tetraCornerIndices,
        const DRT::Element*                     xfemElement,
        const LINALG::SerialDenseMatrix&        xyze_xfemElement,
        const tetgenio&                         out
        ) const
{

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 3; j++)
            tetraCornerNodes[i](j) = out.pointlist[tetraCornerIndices[i]*3+j];

        elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement, tetraCornerNodes[i]);
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    lifts the higher-order point of an edge of the   u.may 09/07|
 |          linearized interface onto the curved interface              |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::computeHigherOrderPoint(
        const int                                 index1,
        const int                                 index2,
        const int                                 faceIndex,
        const int                                 faceMarker,
        const int                                 globalHigherOrderIndex,
        const vector<int>&                        tetraCornerIndices,
        const vector<LINALG::Matrix<3,1> >&       tetraCornerNodes,
        const DRT::Element*                       xfemElement,
        const LINALG::SerialDenseMatrix&          xyze_xfemElement,
        const map<int,LINALG::Matrix<3,1> >&      currentcutterpositions,
        tetgenio&                                 out
        )
{

    bool                                    intersected             = false;
    bool                                    intersectionNormal      = true;
    int                                     lineIndex               = -1;
    int                                     adjacentFaceMarker      = -1;
    int                                     adjacentFaceIndex       = -1;
    static LINALG::Matrix<3,1> xsi;
    xsi.Clear();

    findAdjacentFace(  tetraCornerIndices[index1], tetraCornerIndices[index2],
                       faceMarker, adjacentFaceMarker, faceIndex, adjacentFaceIndex, out);


    // edge lies within the xfem element
    if(adjacentFaceMarker  > -1)
    {
        vector<LINALG::Matrix<3,1> >  plane;
        computeIntersectionNormalB(   tetraCornerIndices[index1], tetraCornerIndices[index2], faceIndex,
                                      adjacentFaceIndex, globalHigherOrderIndex, plane, xfemElement, out);

        // higher order node lies within the cutter element
        if(adjacentFaceMarker == faceMarker)
        {
            DRT::Element* cutterElement = intersectingCutterElements_[faceMarker];
            const LINALG::SerialDenseMatrix xyze_cutterElement(GEO::getCurrentNodalPositions(cutterElement, currentcutterpositions));
            intersected = computeRecoveryNormal(xsi, plane, cutterElement, xyze_cutterElement, false);
            intersectionNormal = true;
        }
        // higher-order point lies on one of the boundary lines of the cutter element
        else if(adjacentFaceMarker != faceMarker)
        {
            int cutterIndex        = -1;
            findCommonCutterLine(currentcutterpositions, faceMarker, adjacentFaceMarker, lineIndex, cutterIndex);

            if(lineIndex != -1)
            {
                DRT::Element* cutterElement = intersectingCutterElements_[cutterIndex];
                const LINALG::SerialDenseMatrix xyze_cutterElement(GEO::getCurrentNodalPositions(cutterElement, currentcutterpositions));
                intersected = computeRecoveryPlane( lineIndex, currentcutterpositions, xsi, plane, cutterElement, xyze_cutterElement);
                intersectionNormal = false;
            }
        }
    }
    // edge lies on the surface of the xfem element
    else if(adjacentFaceMarker  == -1)
    {
        const int oppositeIndex = findEdgeOppositeIndex(  tetraCornerIndices[index1], tetraCornerIndices[index2],
                                                    adjacentFaceIndex, out);

        //printf("oppo = %d\n", oppositeIndex);

        DRT::Element* cutterElement = intersectingCutterElements_[faceMarker];
        const LINALG::SerialDenseMatrix xyze_cutterElement(GEO::getCurrentNodalPositions(cutterElement, currentcutterpositions));
        
        vector<LINALG::Matrix<3,1> >  plane;

        computeIntersectionNormalA(  true,  index1, index2, oppositeIndex, globalHigherOrderIndex,
                                    tetraCornerIndices, tetraCornerNodes, plane, xfemElement, out);

        intersected = computeRecoveryNormal(xsi, plane, cutterElement, xyze_cutterElement, true);
        intersectionNormal = true;

        if(!intersected)
        {
            printf("REFERNCE DOMAIN");
            lineIndex = findIntersectingSurfaceEdge(    xfemElement, xyze_xfemElement, cutterElement, currentcutterpositions,
                                                        tetraCornerNodes[index1], tetraCornerNodes[index2]);
            if(lineIndex != -1)
            {
                intersected = computeRecoveryPlane( lineIndex, currentcutterpositions, xsi, plane, cutterElement, xyze_cutterElement);
                intersectionNormal = false;
            }
        }
    }


    if(intersected)
    {
        storeHigherOrderNode(   intersectionNormal, globalHigherOrderIndex, lineIndex,
                                xsi, intersectingCutterElements_[faceMarker], currentcutterpositions, xfemElement, out);
    }
    else
    {
        countMissedPoints_++;
        std::cout << "faceMarker = " << faceMarker << endl;
        printf("NO INTERSECTION POINT FOUND!!!!! adjacentFaceMarker = %d\n", adjacentFaceMarker);
    }

}



/*----------------------------------------------------------------------*
 |  RCI:    returns the other two point indices belonging    u.may 09/07|
 |          to a triface that obtains a Steiner point                   |
 *----------------------------------------------------------------------*/
vector<int> GEO::IntersectionExp::getPointIndices(
        const tetgenio&   out,
        const int         trifaceIndex,
        const int         steinerPointIndex
        ) const
{

    int count = 0;
    vector<int> pointIndices(2);

    for(int i = 0; i < 3; i++)
        if(i != steinerPointIndex)
            pointIndices[count++] = out.trifacelist[trifaceIndex*3+i];

    return pointIndices;
}



/*----------------------------------------------------------------------*
 |  RCI:    computes the intersection between a              u.may 08/07|
 |          line and a surface                                          |
 *----------------------------------------------------------------------*/
bool GEO::IntersectionExp::computeRecoveryNormal(
    LINALG::Matrix<3,1>&                      xsi,
    const vector<LINALG::Matrix<3,1> >&       normal,
    const DRT::Element*                       cutterElement,
    const LINALG::SerialDenseMatrix&          xyze_cutterElement,
    const bool                                onBoundary
    ) const
{
    bool                        intersection = true;
    int                         iter = 0;
    int                         countSingular = 0;
    const int                   maxiter = 20;
    double                      residual = 1.0;
    static LINALG::Matrix<3,3> A;
    static LINALG::Matrix<3,1> b;
    static LINALG::Matrix<3,1> dx;
    b.Clear();
    dx.Clear();

    xsi.Clear();
    updateRHSForRCINormal( b, xsi, normal, cutterElement, xyze_cutterElement, onBoundary);

    while(residual > GEO::TOL13)
    {
        updateAForRCINormal( A, xsi, normal, cutterElement, xyze_cutterElement, onBoundary);

        if(!solveLinearSystemWithSVD<3>(A, b, dx, GEO::TOL13))
            countSingular++;

        if(countSingular > 10)
        {
            intersection = false;
            break;
        }

        xsi += dx;
        //printf("dx0 = %20.16f\t, dx1 = %20.16f\t, dx2 = %20.16f\n", dx(0), dx(1), dx(2));
        if(iter >= maxiter || GEO::SumOfFabsEntries(xsi) > GEO::TOLPLUS8)
        {
            intersection = false;
            break;
        }

        updateRHSForRCINormal( b, xsi, normal, cutterElement, xyze_cutterElement, onBoundary);
        residual = b.Norm2();
        iter++;
   
    }

    if( (xsi(0) > ( 1.0 + TOL7))  || (xsi(1) > ( 1.0 + TOL7)) ||    // line coordinate may be bigger than 1
        (xsi(0) < (-1.0 - TOL7))  || (xsi(1) < (-1.0 - TOL7)) ) 
    {
        //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi(0), xsi(1), xsi(2), residual, TOL14);
        // cout << "false" << endl; 
        intersection = false;
    }

    return intersection;
}



/*----------------------------------------------------------------------*
 |  RCI:    updates the systemmatrix for the                 u.may 08/07|
 |          computation of a curve-surface intersection                 |
 |          for the recovery of the curved interface                    |
 |          (line-surface intersection)                                 |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::updateAForRCINormal(
    LINALG::Matrix<3,3>&                        A,
    const LINALG::Matrix<3,1>&                  xsi,
    const vector<LINALG::Matrix<3,1> >&         normal,
    const DRT::Element*                         surfaceElement,
    const LINALG::SerialDenseMatrix&            xyze_surfaceElement,
    const bool                                  onBoundary
    ) const
{
    const int numNodesSurface = surfaceElement->NumNode();

    A = 0.0;
    LINALG::SerialDenseMatrix surfaceDeriv1(2,numNodesSurface);
    DRT::UTILS::shape_function_2D_deriv1(surfaceDeriv1, xsi(0), xsi(1), surfaceElement->Shape());

    if(!onBoundary)
    {
        dsassert(normal.size() >= 2, "mismatch in length");
        for(int i=0; i<numNodesSurface; i++)
        {
            for(int dim=0; dim<3; dim++)
            {
                A(dim,0) += xyze_surfaceElement(dim,i) * surfaceDeriv1(0, i);
                A(dim,1) += xyze_surfaceElement(dim,i) * surfaceDeriv1(1, i);
            }
        }

        for(int dim=0; dim<3; dim++)
            A(dim,2) -= 0.5*( -normal[0](dim) + normal[1](dim));
    }
    else
    {
        const int numNodesLine = 3;
        LINALG::SerialDenseMatrix lineDeriv1(1,numNodesLine);
        DRT::UTILS::shape_function_1D_deriv1(lineDeriv1, xsi(2), DRT::Element::line3);

        for(int i=0; i<numNodesSurface; i++)
        {
            for(int dim=0; dim<3; dim++)
            {
                A(dim,0) += xyze_surfaceElement(dim,i) * surfaceDeriv1(0,i);
                A(dim,1) += xyze_surfaceElement(dim,i) * surfaceDeriv1(1,i);
            }
        }

        for(int i = 0; i < numNodesLine; i++)
            for(int dim=0; dim<3; dim++)
            {
                int index = i;
                if(i > 1)   index = 4;
                A(dim,2) -= normal[index](dim) * lineDeriv1(0,i);
            }
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    updates the right-hand-side for the              u.may 08/07|
 |          computation of a curve-surface intersection                 |
 |          for the recovery of the curved surface (RCI)                |
 |          (line-surface intersection)                                 |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::updateRHSForRCINormal(
    LINALG::Matrix<3,1>&                  b,
    const LINALG::Matrix<3,1>&            xsi,
    const vector<LINALG::Matrix<3,1> >&   normal,
    const DRT::Element*                   surfaceElement,
    const LINALG::SerialDenseMatrix&      xyze_surfaceElement,
    const bool                            onBoundary
    ) const
{
    const int numNodesSurface = surfaceElement->NumNode();
    LINALG::SerialDenseVector surfaceFunct(numNodesSurface);
    DRT::UTILS::shape_function_2D(surfaceFunct, xsi(0), xsi(1), surfaceElement->Shape());

    b = 0.0;
    if(!onBoundary)
    {
        dsassert(normal.size() >= 2, "mismatch in length");
        for(int i=0; i<numNodesSurface; i++)
        {
            for(int dim=0; dim<3; dim++)
                b(dim) -= xyze_surfaceElement(dim,i) * surfaceFunct(i);
        }

        for(int dim=0; dim<3; dim++)
            b(dim) += 0.5*(normal[0](dim)*(1.0 - xsi(2)) + normal[1](dim)*(1.0 + xsi(2)));
    }
    else
    {
        const int numNodesLine = 3;
        LINALG::SerialDenseVector lineFunct(numNodesLine);
        DRT::UTILS::shape_function_1D(lineFunct, xsi(2), DRT::Element::line3);

        for(int i=0; i<numNodesSurface; i++)
        {
            for(int dim=0; dim<3; dim++)
                b(dim) -= xyze_surfaceElement(dim,i) * surfaceFunct(i);
        }


        for(int i=0; i<numNodesLine; i++)
            for(int dim=0; dim<3; dim++)
            {
                int index = i;
                if(i > 1)   index = 4;

                b(dim) += normal[index](dim) * lineFunct(i);
            }
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    computes the intersection between a              u.may 08/07|
 |          curve and a plane                      RCI                  |
 *----------------------------------------------------------------------*/
bool GEO::IntersectionExp::computeRecoveryPlane(
        int&                                        lineIndex,
        const map<int,LINALG::Matrix<3,1> >&        currentcutterpositions,
        LINALG::Matrix<3,1> &                       xsi,
        const vector<LINALG::Matrix<3,1> >&         plane,
        DRT::Element*                               cutterElement,
        const LINALG::SerialDenseMatrix&            xyze_cutterElement
        ) const
{
  const int     numLines = cutterElement->NumLine();

  // run over all lines (curves)
  int     begin, end;
  if(lineIndex == -1)
  {
      begin = 0;
      end = numLines;
  }
  else
  {
      begin = lineIndex;
      end   = lineIndex + 1;
  }

  bool    intersection = true;
  for(int i = begin; i < end; i++)
  {
    int                         iter = 0;
    const int                   maxiter = 20;
    double                      residual = 1.0;
    const vector<RCP<DRT::Element> > cutterElementLines = cutterElement->Lines();
    DRT::Element*               lineElement = (cutterElementLines[i]).get();
    const LINALG::SerialDenseMatrix xyze_lineElement(GEO::getCurrentNodalPositions(lineElement, currentcutterpositions));
    
    static LINALG::Matrix<3,3> A;
    static LINALG::Matrix<3,1> b;
    static LINALG::Matrix<3,1> dx;  dx.Clear();

    intersection = true;
    xsi.Clear();

    updateRHSForRCIPlane( b, xsi, plane, lineElement, xyze_lineElement);

    while(residual > GEO::TOL13)
    {
      updateAForRCIPlane( A, xsi, plane, lineElement, xyze_lineElement, cutterElement, xyze_cutterElement);

      if(!gaussElimination<true,3>(A, b, dx, GEO::TOL13))
      {
          intersection = false;
          break;
      }

      if(iter >= maxiter)
      {
          intersection = false;
          break;
      }

      xsi += dx;
      
      updateRHSForRCIPlane( b, xsi, plane, lineElement, xyze_lineElement);
      residual = b.Norm2();
      iter++;
      //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi(0), xsi(1), xsi(2), residual, TOL14);
    }

    if( (xsi(2) > ( 1.0 + TOL7))  || (xsi(2) < ( -1.0 - TOL7)) )     // planes coordinate may be bigger than 1
    {   //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi(0), xsi(1), xsi(2), residual, TOL14);
        intersection = false;
    }

    if(intersection)
    {
        lineIndex = begin;
        break;
    }
  }

  return intersection;
}



/*----------------------------------------------------------------------*
 |  RCI:    updates the systemmatrix for the                 u.may 09/07|
 |          computation of a curve-surface intersection                 |
 |          for the recovery of the curved interface                    |
 |          (curve-plane intersection)                                  |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::updateAForRCIPlane(
        LINALG::Matrix<3,3>&                  A,
        const LINALG::Matrix<3,1>&            xsi,
        const vector<LINALG::Matrix<3,1> >&   plane,
        const DRT::Element*                   lineElement,
        const LINALG::SerialDenseMatrix&      xyze_lineElement,
        const DRT::Element*                   surfaceElement,
        const LINALG::SerialDenseMatrix&      xyze_surfaceElement
        ) const
{
    const int numNodesLine = lineElement->NumNode();
    const int numNodesSurface = 4;

    static LINALG::Matrix<2,numNodesSurface> surfaceDeriv;
    DRT::UTILS::shape_function_2D_deriv1(surfaceDeriv, xsi(0),  xsi(1), DRT::Element::quad4);
    static LINALG::SerialDenseMatrix lineDeriv(1,numNodesLine);
    DRT::UTILS::shape_function_1D_deriv1(lineDeriv, xsi(2), lineElement->Shape());

    dsassert((int)plane.size() >= numNodesSurface, "plane array has to have size numNodesSurface ( = 4)!");

    A.Clear();
    for(int dim=0; dim<3; dim++)
        for(int i=0; i<numNodesSurface; i++)
        {
            A(dim,0) += plane[i](dim) * surfaceDeriv(0,i);
            A(dim,1) += plane[i](dim) * surfaceDeriv(1,i);
        }

    for(int i=0; i<numNodesLine; i++)
    {
        for(int dim=0; dim<3; dim++)
            A(dim,2) -= xyze_lineElement(dim,i) * lineDeriv(0,i);
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    updates the right-hand-side for the              u.may 09/07|
 |          computation of a curve-surface intersection                 |
 |          for the recovery of the curved surface (RCI)                |
 |          (curve-plane intersection)                                  |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::updateRHSForRCIPlane(
    LINALG::Matrix<3,1>&                  b,
    const LINALG::Matrix<3,1>&            xsi,
    const vector<LINALG::Matrix<3,1> >&   plane,
    const DRT::Element*                         lineElement,
    const LINALG::SerialDenseMatrix&            xyze_lineElement
    ) const
{
    const int numNodesLine    = lineElement->NumNode();
    const int numNodesSurface = 4;

    LINALG::Matrix<numNodesSurface,1> surfaceFunct;
    DRT::UTILS::shape_function_2D(surfaceFunct, xsi(0), xsi(1), DRT::Element::quad4 );
    LINALG::SerialDenseVector lineFunct(numNodesLine);
    DRT::UTILS::shape_function_1D(lineFunct, xsi(2), lineElement->Shape());

    dsassert((int)plane.size() >= numNodesSurface, "plane array has to have size numNodesSurface ( = 4)!");

    b.Clear();
    for(int dim=0; dim<3; dim++)
        for(int i=0; i<numNodesSurface; i++)
            b(dim) -= plane[i](dim) * surfaceFunct(i);

    for(int i=0; i<numNodesLine; i++)
        for(int dim=0; dim<3; dim++)
            b(dim) +=  xyze_lineElement(dim,i)  * lineFunct(i);
}



/*----------------------------------------------------------------------*
 |  RCI:    computes the normal to the interface edge of     u.may 08/07|
 |          the tetrahedron facet lying within this facet               |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::computeIntersectionNormalA(
    const bool                              onBoundary,
    const int                               index1,
    const int                               index2,
    const int                               oppositePointIndex,
    const int                               globalHigherOrderIndex,
    const vector<int>&                      tetraCornerIndices,
    const vector<LINALG::Matrix<3,1> >&     tetraCornerNodes,
    vector<LINALG::Matrix<3,1> >&           plane,
    const DRT::Element*                     xfemElement,
    const tetgenio&                         out) const
{

    LINALG::Matrix<3,1>  p1;
    LINALG::Matrix<3,1>  p2;
    LINALG::Matrix<3,1>  p3;


    if(!onBoundary)
    {
        for(int i=0; i<3; i++)
        {
            p1(i) = tetraCornerNodes[3](i);
            p2(i) = tetraCornerNodes[index1](i);
            p3(i) = tetraCornerNodes[index2](i);
        }
    }
    else
    {
        for(int i=0; i<3; i++)
        {
            p1(i) = out.pointlist[oppositePointIndex*3          + i];
            p2(i) = out.pointlist[tetraCornerIndices[index1]*3  + i];
            p3(i) = out.pointlist[tetraCornerIndices[index2]*3  + i];
        }
    }

    // compute direction vectors of the plane
    LINALG::Matrix<3,1>  r1; r1.Update(1.0, p1, -1.0, p2);
    LINALG::Matrix<3,1>  r2; r2.Update(1.0, p3, -1.0, p2);

    // normal of the plane
    LINALG::Matrix<3,1>  n = computeCrossProduct(r1, r2);
    //normalizeVectorInPLace(n);
    n.Scale(1.0/n.Norm2());

    // direction vector of the intersection line
    LINALG::Matrix<3,1>  r = computeCrossProduct(n, r2);
    //normalizeVectorInPLace(r);
    r.Scale(1.0/r.Norm2());

    // computes the start point of the line
    LINALG::Matrix<3,1>  midpoint(true);

    if(!onBoundary)
      midpoint = computeLineMidpoint(p2, p3);
    else
    {
        for(int i = 0; i < 3; i++)
          midpoint(i) = out.pointlist[globalHigherOrderIndex*3+i];
    }

    // compute nodes of the normal to the interface edge of the tetrahedron
    plane.clear();
    plane.reserve(5);
    LINALG::Matrix<3,1> planeNode;
    planeNode.Update(1.0, midpoint,  1.0, r);
    plane.push_back(planeNode);
    planeNode.Update(1.0, midpoint, -1.0, r);
    plane.push_back(planeNode);
    planeNode.Update(1.0, plane[1],  1.0, n);
    plane.push_back(planeNode);
    planeNode.Update(1.0, plane[0],  1.0, n); // check sign
    plane.push_back(planeNode);


    if(onBoundary)
    {
        for(int i = 0; i < 4; i++)
            elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement_, plane[i]);

        elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement_, midpoint);
        plane.push_back(midpoint);
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    computes the normal to the interface edge of     u.may 11/07|
 |          two adjacent triangular faces                               |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::computeIntersectionNormalB(
        const int                                 index1,
        const int                                 index2,
        const int                                 faceIndex,
        const int                                 adjacentFaceIndex,
        const int                                 globalHigherOrderIndex,
        vector<LINALG::Matrix<3,1> >&             plane,
        const DRT::Element*                       xfemElement,
        const tetgenio&                           out
        ) const
{

    int oppositePointIndex = -1;
    int adjacentOppositePointIndex = -1;

    for(int i = 0; i < 3; i++)
        if((out.trifacelist[faceIndex*3+i] != index1) && (out.trifacelist[faceIndex*3+i] != index2 ))
        {
            oppositePointIndex = out.trifacelist[faceIndex*3+i];
            break;
        }

    for(int i = 0; i < 3; i++)
        if((out.trifacelist[adjacentFaceIndex*3+i] != index1) && (out.trifacelist[adjacentFaceIndex*3+i] != index2 ))
        {
            adjacentOppositePointIndex = out.trifacelist[adjacentFaceIndex*3+i];
            break;
        }

    // compute averageNormal of two faces
    LINALG::Matrix<3,1> p1;
    LINALG::Matrix<3,1> p2;
    LINALG::Matrix<3,1> p3;
    LINALG::Matrix<3,1> p4;

    for(int i=0; i<3; i++)
    {
        p1(i) = out.pointlist[index1*3 + i];
        p2(i) = out.pointlist[index2*3 + i];
        p3(i) = out.pointlist[oppositePointIndex*3 + i];
        p4(i) = out.pointlist[adjacentOppositePointIndex*3 + i];
    }

    elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement_, p1);
    elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement_, p2);
    elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement_, p3);
    elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement_, p4);

    LINALG::Matrix<3,1> r1; r1.Update(1.0, p1, -1.0, p2);
    LINALG::Matrix<3,1> r2; r2.Update(1.0, p3, -1.0, p2);
    LINALG::Matrix<3,1> r3; r3.Update(1.0, p4, -1.0, p2);

    LINALG::Matrix<3,1> n1 = computeCrossProduct(r2, r1);
    LINALG::Matrix<3,1> n2 = computeCrossProduct(r1, r3);

    LINALG::Matrix<3,1> averageNormal; averageNormal.Update(1.0, n1, 1.0, n2);
    LINALG::Matrix<3,1> rPlane = computeCrossProduct(n1, r1);

//    for(int i = 0; i < 3; i++)
//        averageNormal(i) = 0.5*averageNormal(i);
    averageNormal.Scale(0.5);

    //normalizeVectorInPLace(averageNormal);
    //normalizeVectorInPLace(rPlane);
    
    averageNormal.Scale(1.0/averageNormal.Norm2());
    rPlane.Scale(1.0/rPlane.Norm2());

    LINALG::Matrix<3,1> midpoint;
    for(int i = 0; i < 3; i++)
      midpoint(i) = out.pointlist[globalHigherOrderIndex*3+i];

    elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement_, midpoint);

    // compute nodes of the normal to the interface edge of the tetrahedron
    plane.clear();
    plane.reserve(4);
    LINALG::Matrix<3,1> planeNode;
    planeNode.Update(1.0, midpoint,  1.0, averageNormal);
    plane.push_back(planeNode);
    planeNode.Update(1.0, midpoint, -1.0, averageNormal);
    plane.push_back(planeNode);
    planeNode.Update(1.0, plane[1],  1.0, rPlane);
    plane.push_back(planeNode);
    planeNode.Update(1.0, plane[0],  1.0, rPlane); //TODO check signs
    plane.push_back(planeNode);
    

   /* for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 3 ; j++)
            printf("plane = %f\t", plane[i](j));

        printf("\n");
    }
    */
}



/*----------------------------------------------------------------------*
 |  RCI:    computes the normal to the interface edge of     u.may 08/07|
 |          the tetrahedron facet lying within this facet               |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::computeIntersectionNormalC(
    const int                               steinerIndex,
    const int                               edgeIndex,
    const int                               oppositeIndex,
    vector<LINALG::Matrix<3,1> >&           plane,
    const DRT::Element*                     xfemElement,
    const tetgenio&                         out) const
{

    LINALG::Matrix<3,1>  p1;
    LINALG::Matrix<3,1>  p2;
    LINALG::Matrix<3,1>  p3;


    for(int i=0; i<3; i++)
    {
        p1(i) = out.pointlist[oppositeIndex*3  + i];
        p2(i) = out.pointlist[steinerIndex*3   + i];
        p3(i) = out.pointlist[edgeIndex*3      + i];
    }

    // compute direction vectors of the plane
    LINALG::Matrix<3,1>  r1;  r1.Update(1.0, p1, -1.0,  p2);
    LINALG::Matrix<3,1>  r2;  r2.Update(1.0, p3, -1.0,  p2);

    // normal of the plane
    LINALG::Matrix<3,1>  n = computeCrossProduct(r1, r2);
    //normalizeVectorInPLace(n);
    n.Scale(1.0/n.Norm2());
    // direction vector of the intersection line
    LINALG::Matrix<3,1>  r = computeCrossProduct(n, r2);
    //normalizeVectorInPLace(r);
    r.Scale(1.0/r.Norm2());

    // compute nodes of the normal to the interface edge of the tetrahedron
    plane.clear();
    plane.reserve(5);
    LINALG::Matrix<3,1> planeNode;
    planeNode.Update(1.0, p2,  1.0, r);
    plane.push_back(planeNode);
    planeNode.Update(1.0, p2, -1.0, r);
    plane.push_back(planeNode);
    planeNode.Update(1.0, plane[1],  1.0, n);
    plane.push_back(planeNode);
    planeNode.Update(1.0, plane[0],  1.0, n); //TODO check signs
    plane.push_back(planeNode);


    for(int i = 0; i < 4; i++)
      elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement_, plane[i]);
    
    elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement_, p2);
   
    plane.push_back(p2);
}



/*----------------------------------------------------------------------*
 |  RCI:    computes the midpoint of a line                  u.may 08/07|
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,1> GEO::IntersectionExp::computeLineMidpoint(
    const LINALG::Matrix<3,1>& p1,
    const LINALG::Matrix<3,1>& p2) const
{
    LINALG::Matrix<3,1> midpoint;
    midpoint.Update(0.5, p1, 0.5, p2);
    
    //for(int i=0; i<3; i++)
    //    midpoint(i) = (p1(i) + p2(i))*0.5;

    return midpoint;
}



/*----------------------------------------------------------------------*
 |  RCI:    searches for the face marker                     u.may 10/07|
 |          of a facet adjacent to a given edge of                      |
 |          of a given facet                                            |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::findAdjacentFace(
    const int         edgeIndex1,
    const int         edgeIndex2,
    const int         faceMarker,
    int&              adjacentFaceMarker,
    const int         faceIndex,
    int&              adjacentFaceIndex,
    const tetgenio&   out
    ) const
{

    bool    faceMarkerFound     = false;

    for(int i=0; i<out.numberoftrifaces; i++)
    {
        adjacentFaceMarker = out.trifacemarkerlist[i] - facetMarkerOffset_;
        adjacentFaceIndex = i;

        if((adjacentFaceMarker > -2) && (faceIndex != adjacentFaceIndex))
        {
            int countPoints = 0;
            for(int j = 0; j < 3 ; j++)
            {
                const int pointIndex = out.trifacelist[ i*3 + j ];
                if(pointIndex == edgeIndex1 || pointIndex == edgeIndex2)
                    countPoints++;
            }

            if(countPoints == 2)
                faceMarkerFound = true;
        }
        if(faceMarkerFound)
            break;

    }

    if(!faceMarkerFound)
        adjacentFaceMarker = -2;

}



/*----------------------------------------------------------------------*
 |  RCI:    finds the global index of the point opposite     u.may 08/07|
 |          to an edge in the adjacent triangular face                  |
 *----------------------------------------------------------------------*/
int GEO::IntersectionExp::findEdgeOppositeIndex(
        const int                                 edgeIndex1,
        const int                                 edgeIndex2,
        const int                                 adjacentFaceIndex,
        const tetgenio&                           out
        ) const
{
    int oppositePointIndex = -1;

    for(int i = 0; i < 3; i++)
        if((out.trifacelist[adjacentFaceIndex*3+i] != edgeIndex1) && (out.trifacelist[adjacentFaceIndex*3+i] != edgeIndex2 ))
        {
            oppositePointIndex = out.trifacelist[adjacentFaceIndex*3+i];
            break;
        }

    return oppositePointIndex;
}



/*----------------------------------------------------------------------*
 |  RCI:    searchs for the common edge of two               u.may 10/07|
 |          adjacent facets                                             |
 *----------------------------------------------------------------------*/
bool GEO::IntersectionExp::findCommonFaceEdge(
        const int                           faceIndex1,
        const int                           faceIndex2,
        const vector<int>&                  adjacentFacesList,
        LINALG::Matrix<3,1>&                edgePoint,
        LINALG::Matrix<3,1>&                oppositePoint,
        const tetgenio&                     out
        ) const
{
    bool edgeFound = false;

    for(int i = 0; i < 2; i++)
    {
        for(int j = 0; j < 2; j++)
        {
            if(adjacentFacesList[faceIndex1*2+i+1] == adjacentFacesList[faceIndex2*2+j+1])
            {
                int index = 0;
                if(i == 0)  index = 1;

                for(int k = 0; k < 3; k++)
                {
                    edgePoint(k)        = out.pointlist[adjacentFacesList[faceIndex1*2+i+1]*3 + k];
                    oppositePoint(k)    = out.pointlist[adjacentFacesList[faceIndex1*2+index+1]*3 + k];
                }
                edgeFound = true;
                break;
            }
        }
        if(edgeFound)
            break;
    }

    return edgeFound;
}



/*----------------------------------------------------------------------*
 |  RCI:    searches for the common edge of two              u.may 10/07|
 |          adjacent cutter elements                                    |
 |          corresponding to the common face edge                       |
 |          of face 1 and facet 2                                       |
 *----------------------------------------------------------------------*/
bool GEO::IntersectionExp::findCommonCutterLine(
    const map<int,LINALG::Matrix<3,1> >&            currentcutterpositions,
    const int                                       faceIndex1,
    const int                                       faceIndex2,
    int&                                            lineIndex,
    int&                                            cutterIndex
    ) const
{
    bool comparison = false;
    // get line arrays, which are computed onthe fly within the cutter elements
    const vector<RCP<DRT::Element> > lines1 = intersectingCutterElements_[faceIndex1]->Lines();
    const vector<RCP<DRT::Element> > lines2 = intersectingCutterElements_[faceIndex2]->Lines();

    const int numLines1 = intersectingCutterElements_[faceIndex1]->NumLine();
    const int numLines2 = intersectingCutterElements_[faceIndex2]->NumLine();
    const int numNodes  = lines2[0]->NumNode();

    for(int i = 0; i < numLines1; i++)
    {
        for(int j = 0; j < numLines2; j++)
        {
            comparison = true;
            for(int k  = 0; k < numNodes; k++)
            {
                const DRT::Node* node1 = lines1[i]->Nodes()[k];
                const DRT::Node* node2 = lines2[j]->Nodes()[k];
                LINALG::Matrix<3,1> pos1 = currentcutterpositions.find(node1->Id())->second;
                LINALG::Matrix<3,1> pos2 = currentcutterpositions.find(node2->Id())->second;
                if(!comparePoints<3>(pos1, pos2))
                {
                    comparison = false;
                    break;
                }
            }

            if(!comparison)
            {
                comparison = true;
                for(int k  = 0; k < numNodes; k++)
                {
                    if(k==2)
                    {
                        const DRT::Node* node1 = lines1[i]->Nodes()[k];
                        const DRT::Node* node2 = lines2[j]->Nodes()[k];
                        LINALG::Matrix<3,1> pos1 = currentcutterpositions.find(node1->Id())->second;
                        LINALG::Matrix<3,1> pos2 = currentcutterpositions.find(node2->Id())->second;
                        if(!comparePoints<3>(pos1, pos2))
                            comparison = false;
                    }
                    else
                    {
                        const DRT::Node* node1 = lines1[i]->Nodes()[k];
                        const DRT::Node* node2 = lines2[j]->Nodes()[1-k];
                        LINALG::Matrix<3,1> pos1 = currentcutterpositions.find(node1->Id())->second;
                        LINALG::Matrix<3,1> pos2 = currentcutterpositions.find(node2->Id())->second;
                        if(!comparePoints<3>(pos1, pos2))
                            comparison = false;
                    }
                }
            }

            if(comparison)
            {
                lineIndex    =  i;
                cutterIndex =  faceIndex1;
                break;
            }
        }
        if(comparison)
            break;
    }
    return comparison;
}



/*----------------------------------------------------------------------*
 |  RCI:    for the recovery of a higher-order node          u.may 10/07|
 |          by a plane - line element intersection                      |
 |          this method finds the line element of the given cutter      |
 |          element intersecting the plane                              |
 |          checking if the edge nodes of the correspondning            |
 |          facet edge are lying on this line element                   |
 *----------------------------------------------------------------------*/
int GEO::IntersectionExp::findIntersectingSurfaceEdge(
        const DRT::Element*                       xfemElement,
        const LINALG::SerialDenseMatrix&          xyze_xfemElement,
        DRT::Element*                             cutterElement,
        const map<int,LINALG::Matrix<3,1> >&      currentcutterpositions,
        const LINALG::Matrix<3,1>&                edgeNode1,
        const LINALG::Matrix<3,1>&                edgeNode2
        ) const
{
    dserror("to be improved");
    int lineIndex = -1;
    static LINALG::Matrix<3,1> x1;
    static LINALG::Matrix<3,1> x2;

    LINALG::Matrix<3,1> node1 = edgeNode1;
    LINALG::Matrix<3,1> node2 = edgeNode2;

    elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement, node1);
    elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement, node2);

    x1(0) = node1(0);
    x2(0) = node2(0);

    vector<RCP<DRT::Element> > lines = cutterElement->Lines();
//    for(int i = 0; i < cutterElement->NumLine(); i++)
//    {
//        const DRT::Element*  lineElement = lines[i];
//
//        LINALG::Matrix<1,1> xsi1(1);
//        LINALG::Matrix<1,1> xsi2(1);
//        currentToVolumeElementCoordinates( lineElement, xyze_lineElement, x1, xsi1);  // TODO: will crash
//        currentToVolumeElementCoordinates( lineElement, xyze_lineElement, x2, xsi2);
//
//        const bool check1 = checkPositionWithinElementParameterSpace( xsi1, lineElement->Shape());
//        const bool check2 = checkPositionWithinElementParameterSpace( xsi2, lineElement->Shape());
//        if( check1  &&  check2 )
//        {
//            lineIndex = i;  //countIndex;
//            break;
//        }
//    }

    /*
    for(int i = 0; i < lines[lineIndex]->NumNode(); i++ )
    {
        lines[lineIndex]->Nodes()[i]->Print(cout);
        printf("\n");
    }
    */
    return lineIndex;
}



/*----------------------------------------------------------------------*
 |  RCI:    stores the higher-order node in the pointlist    u.may 08/07|
 |          at the place of the linear node                             |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::storeHigherOrderNode(
    const bool                                  normal,
    const int                                   globalHigherOrderIndex,
    const int                                   lineIndex,
    LINALG::Matrix<3,1>&                        xsi,
    DRT::Element*                               cutterElement,
    const map<int,LINALG::Matrix<3,1> >&        currentcutterpositions,
    const DRT::Element*                         xfemElement,
    tetgenio&                                   out
    ) const
{
    static LINALG::Matrix<3,1> curr;

    if(normal)
    {
        LINALG::Matrix<2,1> xsiSurf;
        xsiSurf(0) = xsi(0);
        xsiSurf(1) = xsi(1);
        
        const LINALG::SerialDenseMatrix xyze_cutterElement(GEO::getCurrentNodalPositions(cutterElement, currentcutterpositions));
        elementToCurrentCoordinates(cutterElement->Shape(), xyze_cutterElement, xsiSurf, curr);
    }
    else
    {
        LINALG::Matrix<1,1> xsiLine;
        xsiLine(0) = xsi(2);
        const vector<RCP<DRT::Element> > cutterElementLines = cutterElement->Lines();
        const DRT::Element* lineele = cutterElementLines[lineIndex].get();
        const LINALG::SerialDenseMatrix xyze_lineElement(GEO::getCurrentNodalPositions(lineele, currentcutterpositions));
        elementToCurrentCoordinates(lineele->Shape(), xyze_lineElement, xsiLine, curr);
    }
    xsi = currentToVolumeElementCoordinatesExact(xfemElement->Shape(), xyze_xfemElement_, curr, TOL7);

    //printf("xsiold0 = %20.16f\t, xsiold1 = %20.16f\t, xsiold2 = %20.16f\n", out.pointlist[globalHigherOrderIndex*3], out.pointlist[globalHigherOrderIndex*3+1], out.pointlist[globalHigherOrderIndex*3+2]);

    for(int i = 0; i < 3; i++)
        out.pointlist[globalHigherOrderIndex*3+i]   = xsi(i);

    //printf("xsi0    = %20.16f\t, xsi1    = %20.16f\t, xsi2    = %20.16f\n", xsi(0), xsi(1), xsi(2));
    //printf("\n");
}



/*----------------------------------------------------------------------*
 |  RCI:    stores domain integration cells                  u.may 11/07|
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::addCellsToDomainIntCellsMap(
        const DRT::Element*             xfemElement,
        map< int, DomainIntCells >&     domainintcells,
        const tetgenio&                 out,
        bool 							              higherorder
        ) const
{

    // store domain integration cells
    DomainIntCells   						listDomainICPerElement;
    DRT::Element::DiscretizationType		distype = DRT::Element::tet4;
    if(higherorder)							distype = DRT::Element::tet10;
    
    const int numTetNodes = DRT::UTILS::getNumberOfElementNodes(distype);
    if (out.numberofcorners < numTetNodes)
      dserror("you fool, you need to turn on quadratic tets with tetgen -o2 switch!");
    
    for(int i=0; i<out.numberoftetrahedra; i++ )
    {
        vector< vector<double> > tetrahedronCoord;
        for(int j = 0; j < numTetNodes; j++)
        {
            vector<double> tetnodes(3);
            for(int k = 0; k < 3; k++)
                tetnodes[k] = out.pointlist[out.tetrahedronlist[i*out.numberofcorners+j]*3+k];

            tetrahedronCoord.push_back(tetnodes);
        }
        listDomainICPerElement.push_back(DomainIntCell(distype, tetrahedronCoord));
    }
    domainintcells.insert(make_pair(xfemElement->Id(),listDomainICPerElement));
}


/*----------------------------------------------------------------------*
 |  RCI:    stores boundary integration cells                u.may 11/07|
 |																		|
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::addCellsToBoundaryIntCellsMap(
    const int                         		trifaceIndex,
    const int                         		cornerIndex,
    const int                         		globalHigherOrderIndex,
    const int                         		faceMarker,
    const map<int,LINALG::Matrix<3,1> >&  currentcutterpositions,
    vector<vector<double> >&              domainCoord,
    vector<vector<double> >&              boundaryCoord,
    const DRT::Element*						        xfemElement,
    const LINALG::SerialDenseMatrix&      xyze_xfemElement,
    const tetgenio&               			  out
    ) const
{

    //vector<double> trinodes(3);

//    printf("i length = %d\n", domainCoord.size());
//    printf("j length = %d\n", domainCoord[0].size());
//    for(int i = 0; i < 6; i++)
//    	for(int j = 0; j < 3; j++)
//    		printf("vector = %f\n ", domainCoord[i][j]);


    // store corner node
    {
    static LINALG::Matrix<3,1> eleCoordDomainCorner;
    for(int k = 0; k < 3; k++)
        eleCoordDomainCorner(k) = out.pointlist[(out.trifacelist[trifaceIndex*3+cornerIndex])*3+k];

    domainCoord[cornerIndex][0] = eleCoordDomainCorner(0);
    domainCoord[cornerIndex][1] = eleCoordDomainCorner(1);
    domainCoord[cornerIndex][2] = eleCoordDomainCorner(2);

    static LINALG::Matrix<3,1> physCoordCorner;
    elementToCurrentCoordinates(xfemElement->Shape(), xyze_xfemElement, eleCoordDomainCorner, physCoordCorner);

    //domainCoord.push_back(trinodes);

    const DRT::Element* cutterElement = intersectingCutterElements_[faceMarker];
    const LINALG::SerialDenseMatrix xyze_cutterElement(GEO::getCurrentNodalPositions(cutterElement, currentcutterpositions));
    
    LINALG::Matrix<2,1> eleCoordBoundaryCorner;
//    cout << *cutterElement << endl;
//    cout << "xyze_cutterElement: " << xyze_cutterElement << endl;
    CurrentToSurfaceElementCoordinates(cutterElement->Shape(), xyze_cutterElement, physCoordCorner, eleCoordBoundaryCorner);

//    std::cout << "physcood = " << physCoordCorner << "    ";
//    std::cout << "elecood = " << eleCoordBoundaryCorner << endl;
//    std::cout << "cornerIndex = " << cornerIndex << endl;
    boundaryCoord[cornerIndex][0] = eleCoordBoundaryCorner(0);
    boundaryCoord[cornerIndex][1] = eleCoordBoundaryCorner(1);
    boundaryCoord[cornerIndex][2] = 0.0;

    //boundaryCoord.push_back(trinodes);
    }

    if(globalHigherOrderIndex > -1)
    {
        // store higher order node
        static LINALG::Matrix<3,1> eleCoordDomaninHO;
        for(int k = 0; k < 3; k++)
            eleCoordDomaninHO(k) = out.pointlist[globalHigherOrderIndex*3+k];



        domainCoord[cornerIndex+3][0] = eleCoordDomaninHO(0);
        domainCoord[cornerIndex+3][1] = eleCoordDomaninHO(1);
        domainCoord[cornerIndex+3][2] = eleCoordDomaninHO(2);

        static LINALG::Matrix<3,1> physCoordHO;
        elementToCurrentCoordinates(xfemElement->Shape(), xyze_xfemElement, eleCoordDomaninHO, physCoordHO);

        //domainCoord.push_back(trinodes);

        const DRT::Element* cutterElement = intersectingCutterElements_[faceMarker];
        const LINALG::SerialDenseMatrix xyze_cutterElement(GEO::getCurrentNodalPositions(cutterElement, currentcutterpositions));
        LINALG::Matrix<2,1> eleCoordBoundaryHO;
        CurrentToSurfaceElementCoordinates(cutterElement->Shape(), xyze_cutterElement, physCoordHO, eleCoordBoundaryHO);

//        std::cout << "physcood = " << physCoordHO << "    ";
//        std::cout << "elecood = " << eleCoordBoundaryHO << endl;

        boundaryCoord[cornerIndex+3][0] = eleCoordBoundaryHO(0);
        boundaryCoord[cornerIndex+3][1] = eleCoordBoundaryHO(1);
        boundaryCoord[cornerIndex+3][2] = 0.0;

        //boundaryCoord.push_back(trinodes);
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    stores boundary integration cells                u.may 07/08|
 |          of cells lying on xfem surfaces                             |
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::addXFEMSurfaceCellsToBoundaryIntCellsMap(
    const bool                            higherorder,
    const int                             cornerIndex,
    const int                             index1,
    const int                             index2,
    const int                             cutterPos,
    const map<int,LINALG::Matrix<3,1> >&  currentcutterpositions,
    vector<vector<double> >&              domainCoord,
    vector<vector<double> >&              boundaryCoord,
    const DRT::Element*                   xfemElement,
    const LINALG::SerialDenseMatrix&      xyze_xfemElement
    ) const
{
    // store corner node
    {
    static LINALG::Matrix<3,1> eleCoordDomainCorner;
    for(int k = 0; k < 3; k++)
        eleCoordDomainCorner(k) = pointList_[index1].getCoord()(k);

    domainCoord[cornerIndex][0] = eleCoordDomainCorner(0);
    domainCoord[cornerIndex][1] = eleCoordDomainCorner(1);
    domainCoord[cornerIndex][2] = eleCoordDomainCorner(2);

    static LINALG::Matrix<3,1> physCoordCorner;
    elementToCurrentCoordinates(xfemElement->Shape(), xyze_xfemElement, eleCoordDomainCorner, physCoordCorner);

    const DRT::Element* cutterElement = intersectingCutterElements_[cutterPos];
    const LINALG::SerialDenseMatrix xyze_cutterElement(GEO::getCurrentNodalPositions(cutterElement, currentcutterpositions));
    
    LINALG::Matrix<2,1> eleCoordBoundaryCorner;
    CurrentToSurfaceElementCoordinates(cutterElement->Shape(), xyze_cutterElement, physCoordCorner, eleCoordBoundaryCorner);

    boundaryCoord[cornerIndex][0] = eleCoordBoundaryCorner(0);
    boundaryCoord[cornerIndex][1] = eleCoordBoundaryCorner(1);
    boundaryCoord[cornerIndex][2] = 0.0;

    }

    if(higherorder)
    {
      // store higher order node
      static LINALG::Matrix<3,1> eleCoordDomaninCorner1;
      static LINALG::Matrix<3,1> eleCoordDomaninCorner2;
      static LINALG::Matrix<3,1> eleCoordDomaninHO;
      for(int k = 0; k < 3; k++)
      {
        eleCoordDomaninCorner1(k) = pointList_[index1].getCoord()(k);
        eleCoordDomaninCorner2(k) = pointList_[index2].getCoord()(k);
        eleCoordDomaninHO(k) = 0.5*(eleCoordDomaninCorner1(0) + eleCoordDomaninCorner2(0));
      }   

      domainCoord[cornerIndex+3][0] = eleCoordDomaninHO(0);
      domainCoord[cornerIndex+3][1] = eleCoordDomaninHO(1);
      domainCoord[cornerIndex+3][2] = eleCoordDomaninHO(2);

      static LINALG::Matrix<3,1> physCoordHO;
      elementToCurrentCoordinates(xfemElement->Shape(), xyze_xfemElement, eleCoordDomaninHO, physCoordHO);

      const DRT::Element* cutterElement = intersectingCutterElements_[cutterPos];
      const LINALG::SerialDenseMatrix xyze_cutterElement(GEO::getCurrentNodalPositions(cutterElement, currentcutterpositions));
      LINALG::Matrix<2,1> eleCoordBoundaryHO;
      CurrentToSurfaceElementCoordinates(cutterElement->Shape(), xyze_cutterElement, physCoordHO, eleCoordBoundaryHO);

      boundaryCoord[cornerIndex+3][0] = eleCoordBoundaryHO(0);
      boundaryCoord[cornerIndex+3][1] = eleCoordBoundaryHO(1);
      boundaryCoord[cornerIndex+3][2] = 0.0;
    }
}


/*----------------------------------------------------------------------*
 |  DB:     Debug only                                       u.may 06/07|
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::debugXAABBIntersection(
        const LINALG::Matrix<3,2>& cutterXAABB,
        const LINALG::Matrix<3,2>& xfemXAABB,
        const DRT::Element* cutterElement,
        const DRT::Element* xfemElement,
        const int noC,
        const int noX
        ) const
{
	cout << endl;
    std::cout << "===============================================================" << endl;
	cout << "Debug Intersection of XAABB's" << endl;
	cout << "===============================================================" << endl;
	cout << endl;
	cout << "CUTTER ELEMENT " << noC << " :" << endl;
	cout << endl;
	for(int jE = 0; jE < cutterElement->NumNode(); jE++)
	{
		cutterElement->Nodes()[jE]->Print(cout);
		cout << endl;
	}
	cout << endl;
    std::cout << endl;
    std::cout << "CUTTER XAABB: "<< "                     " << "XFEM XAABB: " << endl;
    std::cout << endl;
    std::cout <<  "minX = " << cutterXAABB(0,0) << "      " << "maxX = " << cutterXAABB(0,1) << "      " ;
    std::cout <<  "minX = " << xfemXAABB(0,0)   << "      " << "maxX = " << xfemXAABB(0,1)   << endl;
    std::cout <<  "minY = " << cutterXAABB(1,0) << "      " << "maxY = " << cutterXAABB(1,1) << "      " ;
    std::cout <<  "minY = " << xfemXAABB(1,0)   << "      " << "maxY = " << xfemXAABB(1,1)   << endl;
    std::cout <<  "minZ = " << cutterXAABB(2,0) << "      " << "maxZ = " << cutterXAABB(2,1) << "      " ;
    std::cout <<  "minZ = " << xfemXAABB(2,0)   << "      " << "maxZ = " << xfemXAABB(2,1) << endl;
    std::cout << endl;
    std::cout << endl;
	cout << "XFEM ELEMENT " << noX << " :" << endl;
	cout << endl;
	for(int jE = 0; jE < xfemElement->NumNode(); jE++)
	{
		xfemElement->Nodes()[jE]->Print(cout);
		cout << endl;
	}
    std::cout << endl;
    std::cout << endl;
    std::cout << "CUTTER XAABB: "<< "                     " << "XFEM XAABB: " << endl;
    std::cout << endl;
    std::cout <<  "minX = " << cutterXAABB(0,0) << "      " << "maxX = " << cutterXAABB(0,1) << "      " ;
    std::cout <<  "minX = " << xfemXAABB(0,0)   << "      " << "maxX = " << xfemXAABB(0,1)   << endl;
    std::cout <<  "minY = " << cutterXAABB(1,0) << "      " << "maxY = " << cutterXAABB(1,1) << "      " ;
    std::cout <<  "minY = " << xfemXAABB(1,0)   << "      " << "maxY = " << xfemXAABB(1,1)   << endl;
    std::cout <<  "minZ = " << cutterXAABB(2,0) << "      " << "maxZ = " << cutterXAABB(2,1) << "      " ;
    std::cout <<  "minZ = " << xfemXAABB(2,0)   << "      " << "maxZ = " << xfemXAABB(2,1) << endl;
	cout << endl;
	cout << endl;
    std::cout << "===============================================================" << endl;
    std::cout << "End Debug Intersection of XAABB's" << endl;
	cout << "===============================================================" << endl;
	cout << endl; std::cout << endl; std::cout << endl;


}



/*----------------------------------------------------------------------*
 |  DB:     Debug only                                       u.may 06/07|
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::debugNodeWithinElement(
        const DRT::Element*         xfemElement,
        const DRT::Node*            node,
        const LINALG::Matrix<3,1>&  xsi,
        const int                   noE,
        const int                   noN,
        const bool                  within
        ) const
{
    dserror("fix");
    const int numnodes = xfemElement->NumNode();
    vector<int> actParams(1,0);
    LINALG::SerialDenseVector funct(numnodes);
    LINALG::Matrix<3,1> x(true);
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;

    params.set("action","calc_Shapefunction");
    actParams[0] = numnodes;

    dserror("we don't use Evaluate anymore, so thius function does not make sence!");
    //element->Evaluate(params, dummyDis, actParams, emptyM, emptyM , funct, xsi, emptyV);

    for(int dim=0; dim<3; dim++)
        for(int i=0; i<numnodes; i++)
        {
            x(dim) += xfemElement->Nodes()[i]->X()[dim] * funct(i);
        }

    std::cout << endl;
    std::cout << "===============================================================" << endl;
    std::cout << "Debug Node within element" << endl;
    std::cout << "===============================================================" << endl;
    std::cout << endl;
    std::cout << "ELEMENT " << noE << " :" << endl;
    std::cout << endl;
/*    for(int jE = 0; jE < element->NumNode(); jE++)
    {
        element->Nodes()[jE]->Print(cout);
        std::cout << endl;
    }
*/
    std::cout << endl;
    std::cout << endl;
    std::cout << "NODE " << noN << " :" << endl;
    std::cout << endl;
        node->Print(cout);
    std::cout << endl;
    std::cout << endl;
    std::cout << "XSI :";
    std::cout << "   r = " << xsi(0) << "     s = " <<  xsi(1) << "     t = "  << xsi(2) << endl;
    std::cout << endl;
    std::cout << endl;
    std::cout << "CURRENT COORDINATES :";
    std::cout << "   x = " << x(0) << "     y = " <<  x(1) << "     z = "  << x(2) << endl;
    std::cout << endl;
    std::cout << endl;
    if(within) std::cout << "NODE LIES WITHIN ELEMENT" << endl;
    else            std::cout << "NODE DOES NOT LIE WITHIN ELEMENT" << endl;
    std::cout << endl;
    std::cout << endl;
    std::cout << "===============================================================" << endl;
    std::cout << "End Debug Node within element" << endl;
    std::cout << "===============================================================" << endl;
    std::cout << endl; std::cout << endl; std::cout << endl;
}



/*----------------------------------------------------------------------*
 |  DB:     Debug only                                       u.may 06/07|
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::debugTetgenDataStructure(
        const DRT::Element*               xfemElement) const
{
    std::cout << endl;
    std::cout << "===============================================================" << endl;
    std::cout << "Debug Tetgen Data Structure " << endl;
    std::cout << "===============================================================" << endl;
    std::cout << endl;
    std::cout << "POINT LIST " << " :" << endl;
    std::cout << endl;
    LINALG::Matrix<3,1> xsi(true);
    for(unsigned int i = 0; i < pointList_.size(); i++)
    {
        for(int j = 0; j< 3; j++)
        {
            xsi(j) = pointList_[i].getCoord()(j);
        }
        elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement_, xsi);

        std::cout << i << ".th point:   ";
        for(int j = 0; j< 3; j++)
        {
            //cout << setprecision(10) << pointList_[i].coord[j] << "\t";
             printf("%20.16f\t", pointList_[i].getCoord()(j) );
        }
        std::cout << endl;
        std::cout << endl;

      /*  for(int j = 0; j< 3; j++)
        {
            std::cout << xsi[j] << "\t";
        }
        std::cout << endl;
        std::cout << endl;*/
    }
    std::cout << endl;
    std::cout << endl;

    std::cout << endl;
    std::cout << "SEGMENT LIST " << " :" << endl;
    std::cout << endl;
    for(unsigned int i = 0; i < segmentList_.size(); i++)
    {
        std::cout << i << ".th segment:   ";
        int count = 0;
        for(unsigned int j = 0; j < segmentList_[i].size(); j++)
                std::cout << segmentList_[i][count++] << "\t";

        count = 0;
        for(unsigned int j = 0; j < isolatedPointList_[i].size(); j++)
                std::cout << isolatedPointList_[i][count++] << "\t";

        std::cout << endl;
        std::cout << endl;
    }
    std::cout << endl;
    std::cout << endl;

    std::cout << endl;
    std::cout << "TRIANGLE LIST " << " :" << endl;
    std::cout << endl;
    for(unsigned int i = 0; i < triangleList_.size(); i++)
    {
        std::cout << i << ".th triangle:   ";
        for(int j = 0; j< 3; j++)
        {
            std::cout << triangleList_[i][j] << "\t";
        }
        std::cout << endl;
        std::cout << endl;
    }
    std::cout << endl;
    std::cout << endl;

    std::cout << "===============================================================" << endl;
    std::cout << "Debug Tetgen Data Structure" << endl;
    std::cout << "===============================================================" << endl;
    std::cout << endl; std::cout << endl; std::cout << endl;
}



/*----------------------------------------------------------------------*
 |  Debug only                                               u.may 06/07|
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::debugTetgenOutput( 	
    tetgenio& in,
    tetgenio& out,
    const DRT::Element*   xfemElement,
    vector<int>& elementIds,
    int timestepcounter) const
{
	std::string tetgenIn = "tetgenPLC";
	std::string tetgenOut = "tetgenMesh";
	char tetgenInId[30];
	char tetgenOutId[30];

	for(unsigned int i = 0; i < elementIds.size(); i++)
	{
		if(xfemElement->Id()== elementIds[i])
		{
			// change filename
			sprintf(tetgenInId,"%s%d%d", tetgenIn.c_str(), elementIds[i], timestepcounter);
			sprintf(tetgenOutId,"%s%d%d", tetgenOut.c_str(), elementIds[i],timestepcounter);

			// write piecewise linear complex
			in.save_nodes(tetgenInId);
	    	in.save_poly(tetgenInId);

	    	// write tetrahedron mesh
	    	out.save_elements(tetgenOutId);
	    	out.save_nodes(tetgenOutId);
	    	out.save_faces(tetgenOutId);

	    	cout << "Saving tetgen output for the " <<  elementIds[i] << ".xfem element" << endl;
            flush(cout);
	    }
    }
}



/*----------------------------------------------------------------------*
 |  DB:     Debug only                                       u.may 09/07|
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::printTetViewOutput(
    int             index,
    tetgenio&       out) const
{

    FILE *outFile;
    char filename[100];


    sprintf(filename, "tetgenMesh%d.node", index);

    outFile = fopen(filename, "w");
    fprintf(outFile, "%d  %d  %d  %d\n", out.numberofpoints, out.mesh_dim,
          out.numberofpointattributes, out.pointmarkerlist != NULL ? 1 : 0);
    for (int i = 0; i < out.numberofpoints; i++)
    {
        fprintf(outFile, "%d  %.16g  %.16g  %.16g", i, out.pointlist[i * 3], out.pointlist[i * 3 + 1], out.pointlist[i * 3 + 2]);

        for (int j = 0; j < out.numberofpointattributes; j++)
        {
            fprintf(outFile, "  %.16g", out.pointattributelist[i * out.numberofpointattributes + j]);
        }
        if (out.pointmarkerlist != NULL)
        {
            fprintf(outFile, "  %d", out.pointmarkerlist[i]);
        }
        fprintf(outFile, "\n");
    }
    fclose(outFile);
}



/*----------------------------------------------------------------------*
 |  DB:     Debug only                                       u.may 09/07|
 *----------------------------------------------------------------------*/
void GEO::IntersectionExp::printTetViewOutputPLC(
        const DRT::Element*   xfemElement,
        int                   index,
        tetgenio&             in) const
{

    FILE *outFile;
    char filename[100];
    LINALG::Matrix<3,1> xsi(true);

    sprintf(filename, "tetgenPLC%d.node", index);

    outFile = fopen(filename, "w");
    fprintf(outFile, "%d  %d  %d  %d\n", in.numberofpoints, in.mesh_dim,
          in.numberofpointattributes, in.pointmarkerlist != NULL ? 1 : 0);
    for (int i = 0; i < in.numberofpoints; i++)
    {

        for(int j = 0; j < 3; j++)
            xsi(j) = in.pointlist[i*3 + j];

        elementToCurrentCoordinatesInPlace(xfemElement->Shape(), xyze_xfemElement_, xsi);

        fprintf(outFile, "%d  %.16g  %.16g  %.16g", i, xsi(0), xsi(1), xsi(2));

        for (int j = 0; j < in.numberofpointattributes; j++)
        {
            fprintf(outFile, "  %.16g", in.pointattributelist[i * in.numberofpointattributes + j]);
        }
        if (in.pointmarkerlist != NULL)
        {
            fprintf(outFile, "  %d", in.pointmarkerlist[i]);
        }
        fprintf(outFile, "\n");
    }
    fclose(outFile);
}


void GEO::IntersectionExp::debugFaceMarker(
		const int       eleId,
		tetgenio&			  out) const
{

  ofstream f_system("element_faceMarker.pos");
  f_system << "View \" Face Markers" << " \" {" << endl;

  for(int iface=0; iface<out.numberoftrifaces; iface++)
  {
    int trifaceMarker = out.trifacemarkerlist[iface] - facetMarkerOffset_;

    if(trifaceMarker > -2)
    {
      LINALG::Matrix<3,3> triface;
      for(int inode = 0; inode < 3; inode++)
        for(int isd = 0; isd < 3; isd++)
          triface(isd,inode) = out.pointlist[out.trifacelist[iface*3+inode]*3 + isd];

      f_system << IO::GMSH::cellWithScalarToString(DRT::Element::tri3, trifaceMarker, triface) << endl;
    }
  }
  f_system << "};" << endl;
  f_system.close();
}


void GEO::IntersectionExp::debugXFEMConditions(
	const RCP<DRT::Discretization>  cutterdis) const
{

	vector< DRT::Condition * >	xfemConditions;
	cutterdis->GetCondition ("XFEMCoupling", xfemConditions);

	ofstream f_system("element_xfemconditions.pos");
	f_system << "View \" XFEM conditions" << " \" {" << endl;

	for(unsigned int i=0; i<xfemConditions.size(); i++)
    {
        map<int, RCP<DRT::Element > >  geometryMap = xfemConditions[i]->Geometry();
        map<int, RCP<DRT::Element > >::iterator iterGeo;
        //if(geometryMap.size()==0)   printf("geometry does not obtain elements\n");
      	//printf("size of %d.geometry map = %d\n",i, geometryMap.size());

        for(iterGeo = geometryMap.begin(); iterGeo != geometryMap.end(); iterGeo++ )
        {
            const DRT::Element*  cutterElement = iterGeo->second.get();
            f_system << IO::GMSH::elementAtInitialPositionToString(i, cutterElement) << endl;
        }
    }

	f_system << "};" << endl;
	f_system.close();
}


void GEO::IntersectionExp::debugIntersection(
    const DRT::Element* xfemElement,
    set<DRT::Element*> cutterElements) const
{
  int count = 0;
  ofstream f_system("intersection.pos");
  f_system << "View \" Intersection" << " \" {" << endl;
  
  f_system << IO::GMSH::elementAtInitialPositionToString(0, xfemElement) << endl;
  
  for(set< DRT::Element* >::iterator i = cutterElements.begin(); i != cutterElements.end(); ++i )
    f_system << IO::GMSH::elementAtInitialPositionToString(count++, (*i)) << endl;
  
  f_system << "};" << endl;
  f_system.close();
}



void GEO::IntersectionExp::debugIntersectionOfSingleElements(
    const DRT::Element*                         xfemElement,
    const DRT::Element*                         cutterElement,
    const std::map<int,LINALG::Matrix<3,1> >&   currentcutterpositions) const
{
  ofstream f_system("intersectionOfSingleElements.pos");
  f_system << "View \" IntersectionOfSingleElements" << " \" {" << endl;
  
  f_system << IO::GMSH::elementAtInitialPositionToString(0, xfemElement) << endl;
  f_system << IO::GMSH::elementAtCurrentPositionToString(1, cutterElement, currentcutterpositions) << endl;
  
  f_system << "};" << endl;
  f_system.close();
}

//! return Gmsh representation of a bounding box
static std::string XAABBToString(
    const double                    scalar,
    const vector<vector<double> >&  XAABB)
{
  const DRT::Element::DiscretizationType distype = DRT::Element::hex8;
  const int numnode = IO::GMSH::distypeToGmshNumNode(distype);
  
  stringstream pos_array_string;
  pos_array_string << "S" << IO::GMSH::distypeToGmshElementHeader(distype)
      << "(";
  for (int i = 0; i<numnode; ++i)
  {
    pos_array_string << scientific << XAABB[i][0] << ",";
    pos_array_string << scientific << XAABB[i][1] << ",";
    pos_array_string << scientific << XAABB[i][2];
    if (i < numnode-1)
    {
      pos_array_string << ",";
    }
  };
  pos_array_string << ")";
  // values
  pos_array_string << "{";
  for (int i = 0; i<numnode; ++i)
  {
    pos_array_string << scientific << scalar;
    if (i < numnode-1)
    {
      pos_array_string << ",";
    }
  };
  pos_array_string << "};";
  return pos_array_string.str();
}


void GEO::IntersectionExp::debugXAABBs(
    const int id,
    const LINALG::SerialDenseMatrix& cutterXAABB,
    const LINALG::SerialDenseMatrix& xfemXAABB) const
{
  char filename[100];
  sprintf(filename, "element_XAABB%d.pos", id);
  
  ofstream f_system(filename);
  f_system << "View \" XAABB " << " \" {" << endl;
  std::vector<std::vector<double> > nodes(8, vector<double>(3, 0.0));
  
  //cutterXAABB
  nodes[0][0] = cutterXAABB(0, 0);
  nodes[0][1] = cutterXAABB(1, 0);
  nodes[0][2] = cutterXAABB(2, 0); // node 0
  nodes[1][0] = cutterXAABB(0, 1);
  nodes[1][1] = cutterXAABB(1, 0);
  nodes[1][2] = cutterXAABB(2, 0); // node 1
  nodes[2][0] = cutterXAABB(0, 1);
  nodes[2][1] = cutterXAABB(1, 1);
  nodes[2][2] = cutterXAABB(2, 0); // node 2
  nodes[3][0] = cutterXAABB(0, 0);
  nodes[3][1] = cutterXAABB(1, 1);
  nodes[3][2] = cutterXAABB(2, 0); // node 3
  nodes[4][0] = cutterXAABB(0, 0);
  nodes[4][1] = cutterXAABB(1, 0);
  nodes[4][2] = cutterXAABB(2, 1); // node 4
  nodes[5][0] = cutterXAABB(0, 1);
  nodes[5][1] = cutterXAABB(1, 0);
  nodes[5][2] = cutterXAABB(2, 1); // node 5
  nodes[6][0] = cutterXAABB(0, 1);
  nodes[6][1] = cutterXAABB(1, 1);
  nodes[6][2] = cutterXAABB(2, 1); // node 6
  nodes[7][0] = cutterXAABB(0, 0);
  nodes[7][1] = cutterXAABB(1, 1);
  nodes[7][2] = cutterXAABB(2, 1); // node 7

  f_system << XAABBToString(id+1, nodes) << endl;
  
  //xfemXAABB
  nodes[0][0] = xfemXAABB(0, 0);
  nodes[0][1] = xfemXAABB(1, 0);
  nodes[0][2] = xfemXAABB(2, 0); // node 0
  nodes[1][0] = xfemXAABB(0, 1);
  nodes[1][1] = xfemXAABB(1, 0);
  nodes[1][2] = xfemXAABB(2, 0); // node 1
  nodes[2][0] = xfemXAABB(0, 1);
  nodes[2][1] = xfemXAABB(1, 1);
  nodes[2][2] = xfemXAABB(2, 0); // node 2
  nodes[3][0] = xfemXAABB(0, 0);
  nodes[3][1] = xfemXAABB(1, 1);
  nodes[3][2] = xfemXAABB(2, 0); // node 3
  nodes[4][0] = xfemXAABB(0, 0);
  nodes[4][1] = xfemXAABB(1, 0);
  nodes[4][2] = xfemXAABB(2, 1); // node 4
  nodes[5][0] = xfemXAABB(0, 1);
  nodes[5][1] = xfemXAABB(1, 0);
  nodes[5][2] = xfemXAABB(2, 1); // node 5
  nodes[6][0] = xfemXAABB(0, 1);
  nodes[6][1] = xfemXAABB(1, 1);
  nodes[6][2] = xfemXAABB(2, 1); // node 6
  nodes[7][0] = xfemXAABB(0, 0);
  nodes[7][1] = xfemXAABB(1, 1);
  nodes[7][2] = xfemXAABB(2, 1); // node 7

  f_system << XAABBToString(0, nodes) << endl;
  
  f_system << "};" << endl;
  f_system.close();
}



#endif  // #ifdef CCADISCRET

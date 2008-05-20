/*!----------------------------------------------------------------------
\file intersection.cpp

\brief  collection of intersection tools for the computation of the
        intersection of two arbitrary discretizations

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
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "intersection.H"

#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_utils_local_connectivity_matrices.H"
#include "../drt_io/io_gmsh.H"

#ifdef PARALLEL
#include <mpi.h>
#endif


/*----------------------------------------------------------------------*
 |  MAIN:   computes the interface between the xfem          u.may 06/07|
 |          discretization and the cutter discretization.               |
 |          It returns a list of intersected xfem elements              |
 |          and their integrations cells.                               |
 *----------------------------------------------------------------------*/
void XFEM::Intersection::computeIntersection( const Teuchos::RCP<DRT::Discretization>  xfemdis,
                                        const Teuchos::RCP<DRT::Discretization>        cutterdis,
                                        const std::map<int,BlitzVec3>&                 currentcutterpositions,
                                        std::map< int, DomainIntCells >&               domainintcells,
                                        std::map< int, BoundaryIntCells >&             boundaryintcells
                                        )
{

    static int timestepcounter_ = -1;
    timestepcounter_++;

    bool xfemIntersection;
    

    countMissedPoints_ = 0;

    const double t_start = ds_cputime();

    //debugXFEMConditions(cutterdis);
    // obtain vector of pointers to all xfem conditions of the cutter discretization
    vector< DRT::Condition * >      xfemConditions;
    cutterdis->GetCondition ("XFEMCoupling", xfemConditions);

    if(xfemConditions.size()==0)
        cout << "number of XFEMCoupling conditions = 0" << endl;

    //  k < xfemdis->NumMyColElements()
    for(int k = 0; k < xfemdis->NumMyColElements(); ++k)
    {
        xfemIntersection = false;
        DRT::Element* xfemElement = xfemdis->lColElement(k);
        initializeXFEM(k, xfemElement);

        std::set< DRT::Element* >        cutterElements;

        // initial positions, since the xfem element does not move
        const BlitzMat xyze_xfemElement(DRT::UTILS::PositionArrayBlitz(xfemElement));
        
        const BlitzMat3x2 xfemXAABB = computeFastXAABB(xfemElement, xyze_xfemElement);

        startPointList();

        for(vector<DRT::Condition*>::const_iterator conditer = xfemConditions.begin(); conditer!=xfemConditions.end(); ++conditer)
        {
            DRT::Condition* xfemCondition = *conditer;
            const map<int, RCP<DRT::Element > > geometryMap = xfemCondition->Geometry();
            map<int, RCP<DRT::Element > >::const_iterator iterGeo;
            //if(geometryMap.size()==0)   printf("geometry does not obtain elements\n");
          	//printf("size of %d.geometry map = %d\n",i, geometryMap.size());

            for(iterGeo = geometryMap.begin(); iterGeo != geometryMap.end(); ++iterGeo )
            {
                DRT::Element*  cutterElement = iterGeo->second.get();
                if(cutterElement == NULL) dserror("geometry does not obtain elements");
                
                // fill current positions into an array
                const BlitzMat xyze_cutterElement(getCurrentNodalPositions(cutterElement, currentcutterpositions));
                
                const BlitzMat3x2    cutterXAABB(computeFastXAABB(cutterElement, xyze_cutterElement));

                const bool intersected = intersectionOfXAABB(cutterXAABB, xfemXAABB);

                if(intersected)
                {
                  cutterElements.insert(cutterElement);
                }
            }// for-loop over all geometryMap
		    }// for-loop over all xfemConditions

        const DRT::Element*const* xfemElementSurfaces = xfemElement->Surfaces();
        const DRT::Element*const* xfemElementLines = xfemElement->Lines();

        for(set< DRT::Element* >::iterator i = cutterElements.begin(); i != cutterElements.end(); ++i )
        {
            DRT::Element* cutterElement = (*i);
            if(cutterElement == NULL) dserror("cutter element is null\n");
            const BlitzMat xyze_cutterElement(getCurrentNodalPositions(cutterElement, currentcutterpositions));
            const DRT::Element*const* cutterElementLines = cutterElement->Lines();
            const DRT::Node*const* cutterElementNodes = cutterElement->Nodes();

            // number of internal points
            int numInternalPoints = 0;
            // number of surface points
            int numBoundaryPoints = 0;
            vector< InterfacePoint >  interfacePoints;

            // collect internal points
            for(int m=0; m<cutterElement->NumLine() ; m++)
            {
                collectInternalPoints( xfemElement, cutterElement, xyze_cutterElement, cutterElementNodes[m], currentcutterpositions,
                                        interfacePoints, numInternalPoints, numBoundaryPoints, k, m);
            }

            // collect intersection points
            for(int m=0; m<xfemElement->NumLine() ; m++)
            {
                const DRT::Element* xfemElementLine = xfemElementLines[m];
                const BlitzMat xyze_xfemElementLine(DRT::UTILS::PositionArrayBlitz(xfemElementLine));
                if(collectIntersectionPoints(   cutterElement, xyze_cutterElement,
                                                xfemElementLine, xyze_xfemElementLine,
                                                interfacePoints, numBoundaryPoints, 0, m, false, xfemIntersection))
                {
                    storeIntersectedCutterElement(cutterElement);
                }
            }

            for(int m=0; m<cutterElement->NumLine() ; m++)
            {
                for(int p=0; p<xfemElement->NumSurface() ; p++)
                {
                    const DRT::Element* xfemElementSurface = xfemElementSurfaces[p];
                    const BlitzMat xyze_xfemElementSurface(DRT::UTILS::PositionArrayBlitz(xfemElementSurface));
                    const DRT::Element* cutterElementLine = cutterElementLines[m];
                    const BlitzMat xyze_cutterElementLine(getCurrentNodalPositions(cutterElementLine, currentcutterpositions));
                    if(collectIntersectionPoints(   xfemElementSurface, xyze_xfemElementSurface,
                                                    cutterElementLine, xyze_cutterElementLine,
                                                    interfacePoints, numBoundaryPoints,
                                                    p, m, true, xfemIntersection))
                    {
                        storeIntersectedCutterElement(cutterElement);
                    }
                }
            }

            // order interface points
            if( interfacePoints.size() > 0)
            {
#ifdef QHULL
                computeConvexHull( xfemElement, xyze_xfemElement, cutterElement, xyze_cutterElement, interfacePoints, numInternalPoints, numBoundaryPoints);
#else
                dserror("Set QHULL flag to use XFEM intersections!!!");
#endif
            }
        }// for-loop over all cutter elements

        if(xfemIntersection)
        {
            //debugTetgenDataStructure(xfemElement);
            computeCDT(xfemElement, xyze_xfemElement, currentcutterpositions, domainintcells, boundaryintcells, timestepcounter_);
        }

    }// for-loop over all  actdis->NumMyColElements()

    //debugDomainIntCells(domainintcells,2);
    const double t_end = ds_cputime()-t_start;
    cout << endl;
    if(countMissedPoints_ > 0)
    	cout << "Number of missed points during the recovery copy = " << countMissedPoints_ << endl;

    cout << "Intersection computed sucessfully in " << t_end  <<  " secs";
    cout << endl << endl;
}


/*----------------------------------------------------------------------*
 |  INIT:   initializes the private members of the           u.may 08/07|
 |          current xfem element                                        |
 *----------------------------------------------------------------------*/
void XFEM::Intersection::initializeXFEM(
    const int             xfemId,
    const DRT::Element*   xfemElement)
{
    const DRT::Element::DiscretizationType xfemDistype = xfemElement->Shape();

    numXFEMSurfaces_ = xfemElement->NumSurface();

    numXFEMCornerNodes_  = DRT::UTILS::getNumberOfElementCornerNodes(xfemDistype);
    
    xyze_xfemElement_    = DRT::UTILS::PositionArrayBlitz(xfemElement);

    eleLinesSurfaces_     = DRT::UTILS::getEleNodeNumbering_lines_surfaces(xfemDistype);
    eleNodesSurfaces_     = DRT::UTILS::getEleNodeNumbering_nodes_surfaces(xfemDistype);
    eleNumberingSurfaces_ = DRT::UTILS::getEleNodeNumberingSurfaces(xfemDistype);
    eleRefCoordinates_    = DRT::UTILS::getEleNodeNumbering_nodes_reference(xfemDistype);

    pointList_.clear();
    triangleList_.clear();

    segmentList_.clear();
    segmentList_.resize(numXFEMSurfaces_);
    surfacePointList_.clear();
    surfacePointList_.resize(numXFEMSurfaces_);

    intersectingCutterElements_.clear();
    faceMarker_.clear();
}


/*----------------------------------------------------------------------*
 |  CLI:    collects points that belong to the interface     u.may 06/07|
 |          and lie within an xfem element                              |
 *----------------------------------------------------------------------*/
bool XFEM::Intersection::collectInternalPoints(
    const DRT::Element*             xfemElement,
    DRT::Element*                   cutterElement,
    const BlitzMat&                 xyze_cutterElement,
    const DRT::Node*                cutterNode,
    const map<int,BlitzVec3>&       currentcutterpositions,
    std::vector< InterfacePoint >&  interfacePoints,
    int&                            numInternalPoints,
    int&                            numBoundaryPoints,
    const int                       elemId,
    const int                       nodeId)
{
    // current nodal position
    const BlitzVec3 x = currentcutterpositions.find(cutterNode->Id())->second;

    static BlitzVec3 xsi;
    currentToVolumeElementCoordinates(xfemElement, x, xsi);
    const bool nodeWithinElement = checkPositionWithinElementParameterSpace(xsi, xfemElement->Shape());
    // debugNodeWithinElement(xfemElement,cutterNode, xsi, elemId ,nodeId, nodeWithinElement);

    if(nodeWithinElement)
    {
        InterfacePoint ip;
        //debugNodeWithinElement(xfemElement,cutterNode,xsi,elemId ,nodeId, nodeWithinElement);

        numInternalPoints++;

        // check if node lies on the boundary of the xfem element
        if(setInterfacePointBoundaryStatus(xfemElement->Shape(), xsi, ip))
            numBoundaryPoints++;

        // intersection coordinates in the surface
        // element element coordinate system
        DRT::UTILS::getNodeCoordinates(nodeId, ip.coord, cutterElement->Shape());

        interfacePoints.push_back(ip);

        storeIntersectedCutterElement(cutterElement);
    }

    return nodeWithinElement;
}



/*----------------------------------------------------------------------*
 |  CLI:    checks if a node that lies within an element     u.may 06/07|
 |          lies on one of its surfaces or nodes                        |
 *----------------------------------------------------------------------*/
bool XFEM::Intersection::setInterfacePointBoundaryStatus(
        const DRT::Element::DiscretizationType  xfemDistype,
        const BlitzVec3&                        xsi,
        InterfacePoint&                         ip
        ) const
{
    bool onSurface = false;
    const int count = DRT::UTILS::getSurfaces(xsi, ip.surfaces, xfemDistype);

    // point lies on one surface
    if(count == 1)
    {
        onSurface = true;
        ip.nsurf = count;
        ip.pType = SURFACE;
    }
    // point lies on line, which has two neighbouring surfaces
    else if(count == 2)
    {
        onSurface = true;
        ip.nsurf = count;
        ip.pType = LINE;
    }
    // point lies on a node, which has three neighbouring surfaces
    else if(count == 3)
    {
        onSurface = true;
        ip.nsurf = count;
        ip.pType = NODE;
    }
    else
    {
        onSurface = false;
        ip.nsurf = 0;
        ip.pType = INTERNAL;
    }

    return onSurface;
}



/*----------------------------------------------------------------------*
 |  CLI:    collects all intersection points of a line and   u.may 06/07|
 |          and a surface                                               |
 *----------------------------------------------------------------------*/
bool XFEM::Intersection::collectIntersectionPoints(
    const DRT::Element*             surfaceElement,
    const BlitzMat&                 xyze_surfaceElement,
    const DRT::Element*             lineElement,
    const BlitzMat&                 xyze_lineElement,
    std::vector<InterfacePoint>&    interfacePoints,
    const int                       numBoundaryPoints,
    const int                       surfaceId,
    const int                       lineId,
    const bool                      lines,
    bool&                           xfemIntersection
    ) const
{
    static BlitzVec3 xsi;
    xsi = 0.0;

    static BlitzVec3 upLimit;
    static BlitzVec3 loLimit;

    // for hex elements
    upLimit  =  1.0;
    loLimit  = -1.0;

    const bool intersected = computeCurveSurfaceIntersection(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit);

    if(intersected)
        addIntersectionPoint( surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit,
                              interfacePoints, surfaceId, lineId, lines);


    // in this case a node of this line lies on the facet of the xfem element
    // but there is no intersection within the element
    if(!((int) interfacePoints.size() == numBoundaryPoints))
        xfemIntersection = true;

    return intersected;
}


/*!
\brief updates the systemmatrix at the corresponding element coordinates
       for the computation of curve surface intersections

\param A                 (out)      : system matrix
\param xsi               (in)       : vector of element coordinates
\param surfaceElement    (in)       : surface element
\param lineElement       (in)       : line element
*/
template<DRT::Element::DiscretizationType surftype,
         DRT::Element::DiscretizationType linetype>
void updateAForCSI(
        BlitzMat3x3&                      A,
        const BlitzVec3&                  xsi,
        const DRT::Element*               surfaceElement,
        const BlitzMat&                   xyze_surfaceElement,
        const DRT::Element*               lineElement,
        const BlitzMat&                   xyze_lineElement
        )
{
    const int numNodesSurface = DRT::UTILS::getNumberOfElementNodes<surftype>();
    const int numNodesLine = DRT::UTILS::getNumberOfElementNodes<linetype>();

    A = 0.0;

    static BlitzMat surfaceDeriv1(2, numNodesSurface);
    DRT::UTILS::shape_function_2D_deriv1(surfaceDeriv1, xsi(0), xsi(1), surftype);
    for(int inode=0; inode<numNodesSurface; inode++)
    {
        for(int isd=0; isd<3; isd++)
        {
            A(isd,0) += xyze_surfaceElement(isd,inode) * surfaceDeriv1(0,inode);
            A(isd,1) += xyze_surfaceElement(isd,inode) * surfaceDeriv1(1,inode);
        }
    }

    static BlitzMat lineDeriv1(2, numNodesLine);
    DRT::UTILS::shape_function_1D_deriv1(lineDeriv1, xsi(2), linetype);
    for(int inode=0; inode<numNodesLine; inode++)
    {
        for(int isd=0; isd<3; isd++)
        {
            A(isd,2) -= xyze_lineElement(isd,inode) * lineDeriv1(0,inode);
        }
    }
}



/*!
\brief updates the rhs at the corresponding element coordinates
       for the computation of curve surface intersections

\param b                 (out)      : right-hand-side
\param xsi               (in)       : vector of element coordinates
\param surfaceElement    (in)       : surface element
\param lineElement       (in)       : line element
*/
template<DRT::Element::DiscretizationType surftype,
         DRT::Element::DiscretizationType linetype>
void updateRHSForCSI(
        BlitzVec3&         b,
        const BlitzVec3&   xsi,
        const DRT::Element*               surfaceElement,
        const BlitzMat&                   xyze_surfaceElement,
        const DRT::Element*               lineElement,
        const BlitzMat&                   xyze_lineElement
        )
{
    const int numNodesSurface = DRT::UTILS::getNumberOfElementNodes<surftype>();
    const int numNodesLine = DRT::UTILS::getNumberOfElementNodes<linetype>();

    b = 0.0;

    static BlitzVec surfaceFunct(numNodesSurface);
    DRT::UTILS::shape_function_2D(surfaceFunct, xsi(0), xsi(1), surftype);
    for(int i=0; i<numNodesSurface; i++)
    {
        for(int dim=0; dim<3; dim++)
            b(dim) -= xyze_surfaceElement(dim,i) * surfaceFunct(i);
    }

    static BlitzVec lineFunct(numNodesLine);
    DRT::UTILS::shape_function_1D(lineFunct, xsi(2), linetype);
    for(int i=0; i<numNodesLine; i++)
    {
       for(int dim=0; dim<3; dim++)
            b(dim) += xyze_lineElement(dim,i) * lineFunct(i);
    }
}


/*!
    \brief solves a singular system of equations

\param xsi              (in/out)    : vector of element domain coordinates
\param lineElement      (in)        : line element
\param surfaceElement   (in)        : surface element
return true if resulting system is singular , false otherwise
*/
template<DRT::Element::DiscretizationType surftype,
         DRT::Element::DiscretizationType linetype>
bool computeSingularCSI(
        BlitzVec3&                  xsi,
        const DRT::Element*         surfaceElement,
        const BlitzMat&             xyze_surfaceElement,
        const DRT::Element*         lineElement,
        const BlitzMat&             xyze_lineElement
        )
{
    bool singular = false;
    int iter = 0;
    const int maxiter = 5;
    double residual = 1.0;
    static BlitzMat3x3 A;
    static BlitzVec3   b;
    static BlitzVec3   dx;

    updateRHSForCSI<surftype,linetype>( b, xsi, surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement);

    while(residual > XFEM::TOL14)
    {
        updateAForCSI<surftype,linetype>( A, xsi, surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement);

        if(XFEM::solveLinearSystemWithSVD<3>(A, b, dx))
        {
            singular = false;
            xsi += dx;
            break;
        }

        xsi += dx;
        updateRHSForCSI<surftype,linetype>( b, xsi, surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement);
        residual = XFEM::Norm2(b);
        iter++;

        if(iter >= maxiter )
        {
            singular = true;
            break;
        }
    }
    return singular;
}


template<DRT::Element::DiscretizationType surftype,
         DRT::Element::DiscretizationType linetype>
bool computeCurveSurfaceIntersectionT(
    const DRT::Element*               surfaceElement,
    const BlitzMat&                   xyze_surfaceElement,
    const DRT::Element*               lineElement,
    const BlitzMat&                   xyze_lineElement,
    BlitzVec3&                        xsi,
    const BlitzVec3&                  upLimit,
    const BlitzVec3&                  loLimit
    )
{
    bool intersection = true;
    int iter = 0;
    const int maxiter = 30;
    double residual = 1.0;
    static BlitzMat3x3 A;
    static BlitzVec3   b;
    static BlitzVec3   dx;

    updateRHSForCSI<surftype,linetype>( b, xsi, surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement);

    while(residual > XFEM::TOL14)
    {
        updateAForCSI<surftype,linetype>( A, xsi, surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement);

        if(!XFEM::gaussElimination<true, 3, 1>(A, b, dx))
        {
            if(computeSingularCSI<surftype,linetype>(xsi, surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement))
            {
                intersection = false;
                break;
            }
            dx = 0.0;
            iter++;
            printf("SINGULAR\n");
        }

        xsi += dx;
        updateRHSForCSI<surftype,linetype>( b, xsi, surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement);
        residual = XFEM::Norm2(b);
        iter++;

        if(iter >= maxiter )
        {
            intersection = false;
            break;
        }
    }

    if(intersection)
    {
        if( (xsi(0) > (upLimit(0)+XFEM::TOL7)) || (xsi(1) > (upLimit(1)+XFEM::TOL7)) || (xsi(2) > (upLimit(2)+XFEM::TOL7))  ||
            (xsi(0) < (loLimit(0)-XFEM::TOL7)) || (xsi(1) < (loLimit(1)-XFEM::TOL7)) || (xsi(2) < (loLimit(2)-XFEM::TOL7)))
                intersection = false;
    }

    return intersection;
}


/*----------------------------------------------------------------------*
 |  CLI:    computes the intersection between a              u.may 06/07|
 |          curve and a surface                    (CSI)                |
 *----------------------------------------------------------------------*/
bool XFEM::Intersection::computeCurveSurfaceIntersection(
    const DRT::Element*               surfaceElement,
    const BlitzMat&                   xyze_surfaceElement,
    const DRT::Element*               lineElement,
    const BlitzMat&                   xyze_lineElement,
    BlitzVec3&                        xsi,
    const BlitzVec3&                  upLimit,
    const BlitzVec3&                  loLimit
    ) const
{
    if (lineElement->Shape() == DRT::Element::line2)
    {
        switch (surfaceElement->Shape())
        {
        case DRT::Element::quad4:
            return computeCurveSurfaceIntersectionT<DRT::Element::quad4,DRT::Element::line2>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit);
        case DRT::Element::quad8:
            return computeCurveSurfaceIntersectionT<DRT::Element::quad8,DRT::Element::line2>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit);
        case DRT::Element::quad9:
            return computeCurveSurfaceIntersectionT<DRT::Element::quad9,DRT::Element::line2>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit);
        case DRT::Element::tri3:
            return computeCurveSurfaceIntersectionT<DRT::Element::tri3 ,DRT::Element::line2>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit);
        case DRT::Element::tri6:
            return computeCurveSurfaceIntersectionT<DRT::Element::tri6 ,DRT::Element::line2>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit);
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
            return computeCurveSurfaceIntersectionT<DRT::Element::quad4,DRT::Element::line3>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit);
        case DRT::Element::quad8:
            return computeCurveSurfaceIntersectionT<DRT::Element::quad8,DRT::Element::line3>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit);
        case DRT::Element::quad9:
            return computeCurveSurfaceIntersectionT<DRT::Element::quad9,DRT::Element::line3>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit);
        case DRT::Element::tri3:
            return computeCurveSurfaceIntersectionT<DRT::Element::tri3 ,DRT::Element::line3>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit);
        case DRT::Element::tri6:
            return computeCurveSurfaceIntersectionT<DRT::Element::tri6 ,DRT::Element::line3>(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit);
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
int XFEM::Intersection::computeNewStartingPoint(
    const DRT::Element*                surfaceElement,
    const BlitzMat&                    xyze_surfaceElement,
    const DRT::Element*                lineElement,
    const BlitzMat&                    xyze_lineElement,
    const int                          surfaceId,
    const int                          lineId,
    const BlitzVec3&                   xsiOld,
    const BlitzVec3&                   upLimit,
    const BlitzVec3&                   loLimit,
    std::vector<InterfacePoint>&       interfacePoints,
    const bool                         lines
    ) const
{
    bool interval = true;
	int numInterfacePoints = 0;
	static BlitzVec3 xsi;

    //printf("xsi = %f   %f   %f\n", fabs(xsi(0)), fabs(xsi(1)), fabs(xsi(2)) );
    //printf("lolimit = %f   %f   %f\n", fabs(loLimit(0)), fabs(loLimit(1)), fabs(loLimit(2)) );
    //printf("uplimit = %f   %f   %f\n", fabs(upLimit(0)), fabs(upLimit(1)), fabs(upLimit(2)) );

    if(comparePoints<3>(upLimit, loLimit))
        interval = false;

    xsi = upLimit;
    xsi += loLimit;
    xsi *= 0.5;

	bool intersected = computeCurveSurfaceIntersection(surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit);

    if( comparePoints<3>(xsi, xsiOld))
        intersected = false;

	if(intersected && interval)
   		numInterfacePoints = addIntersectionPoint(	surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, xsi, upLimit, loLimit,
   													interfacePoints, surfaceId, lineId, lines);

   	//printf("number of intersection points = %d\n", numInterfacePoints );
   	return numInterfacePoints;
}



/*----------------------------------------------------------------------*
 |  CLI:    adds an intersection point to the 			     u.may 07/07|
 |          list of interface points                       			    |
 *----------------------------------------------------------------------*/
int XFEM::Intersection::addIntersectionPoint(
    const DRT::Element*             surfaceElement,
    const BlitzMat&                 xyze_surfaceElement,
    const DRT::Element*             lineElement,
    const BlitzMat&                 xyze_lineElement,
    const BlitzVec3&                xsi,
    const BlitzVec3&                upLimit,
    const BlitzVec3&                loLimit,
    std::vector<InterfacePoint>& 	interfacePoints,
    const int                       surfaceId,
    const int                       lineId,
    const bool 					    lines
    ) const
{

	int numInterfacePoints = 0;

 	InterfacePoint ip;
    if(lines)
    {
        ip.nsurf = 1;
        ip.surfaces[0] = surfaceId;

        DRT::UTILS::getLineCoordinates(lineId, xsi(2), ip.coord, surfaceElement->Shape());
    }
    else
    {
        ip.nsurf = 2;
        ip.surfaces[0] = eleLinesSurfaces_[lineId][0];
        ip.surfaces[1] = eleLinesSurfaces_[lineId][1];
        ip.coord[0] = xsi(0);
        ip.coord[1] = xsi(1);
    }

    ip.coord[2] = 0.0;
    ip.pType = INTERSECTION;

    vector<InterfacePoint>::iterator it;
    bool alreadyInList = false;
    for(it = interfacePoints.begin(); it != interfacePoints.end(); it++ )
        if(comparePoints<3>(ip.coord, it->coord))
        {
            //printf("alreadyinlist = true\n");
            alreadyInList = true;
            break;

        }

    if(!alreadyInList)
    {
        vector< BlitzVec3 >  upperLimits(8, BlitzVec3(0.0));
        vector< BlitzVec3 >  lowerLimits(8, BlitzVec3(0.0));
        createNewLimits(xsi, upLimit, loLimit, upperLimits, lowerLimits);

        interfacePoints.push_back(ip);
        numInterfacePoints++;

        // recursive call
        for(int i = 0; i < 8; i++)
            numInterfacePoints += computeNewStartingPoint(
                                        surfaceElement, xyze_surfaceElement, lineElement, xyze_lineElement, surfaceId, lineId, xsi,
                                        upperLimits[i], lowerLimits[i], interfacePoints, lines);

    }
	return numInterfacePoints;
}


/*----------------------------------------------------------------------*
 |  CLI:    create new ranges for the recursive              u.may 07/07|
 |          computation of all intersection points                      |
 *----------------------------------------------------------------------*/
void XFEM::Intersection::createNewLimits(
    const BlitzVec3&        xsi,
    const BlitzVec3&        upLimit,
    const BlitzVec3&        loLimit,
    vector< BlitzVec3 >&    upperLimits,
    vector< BlitzVec3 >&    lowerLimits) const
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
 |  ICS:    computes the convex hull of a set of             u.may 06/07|
 |          interface points and stores resulting points,               |
 |          segments and triangles for the use with Tetgen (CDT)        |
 *----------------------------------------------------------------------*/
#ifdef QHULL
void XFEM::Intersection::computeConvexHull(
        const DRT::Element*     xfemElement,
        const BlitzMat&         xyze_xfemElement,
        const DRT::Element*     cutterElement,
        const BlitzMat&         xyze_cutterElement,
        vector<InterfacePoint>& interfacePoints,
        const int               numInternalPoints,
        const int               numBoundaryPoints)
{

    vector<int>                 positions;
    vector<double>              searchPoint(3,0);
    vector<double>              vertex(3,0);
    vector< vector<double> >    vertices;
    InterfacePoint              midpoint;


    if(interfacePoints.size() > 2)
    {
        midpoint = computeMidpoint(interfacePoints);
        // transform it into current coordinates
        {
            static BlitzVec2    eleCoordSurf;
            for(int j = 0; j < 2; j++)
                eleCoordSurf(j)  = midpoint.coord[j];
            static BlitzVec3 curCoordVol;
            elementToCurrentCoordinates(cutterElement, xyze_cutterElement, eleCoordSurf, curCoordVol);
            const BlitzVec3 eleCoordVol(currentToVolumeElementCoordinatesExact(xfemElement, curCoordVol));
            for(int j = 0; j < 3; j++)
                midpoint.coord[j] = eleCoordVol(j);
        }

        // store coordinates in
        // points has numInterfacePoints*dim-dimensional components
        // points[0] is the first coordinate of the first point
        // points[1] is the second coordinate of the first point
        // points[dim] is the first coordinate of the second point
        coordT* coordinates = (coordT *)malloc((2*interfacePoints.size())*sizeof(coordT));
        int fill = 0;
        for(vector<InterfacePoint>::iterator ipoint = interfacePoints.begin(); ipoint != interfacePoints.end(); ++ipoint)
        {
            for(int j = 0; j < 2; j++)
            {
                coordinates[fill++] = ipoint->coord[j];
               // printf("coord = %f\t", ipoint->coord[j]);
            }
            // printf("\n");
            // transform interface points into current coordinates
            {
                static BlitzVec2  eleCoordSurf;
                for(int j = 0; j < 2; j++)
                    eleCoordSurf(j)  = ipoint->coord[j];
                static BlitzVec3 curCoordVol;
                elementToCurrentCoordinates(cutterElement, xyze_cutterElement, eleCoordSurf, curCoordVol);
                const BlitzVec3 eleCoordVol(currentToVolumeElementCoordinatesExact(xfemElement, curCoordVol));
                for(int j = 0; j < 3; j++)
                    ipoint->coord[j] = eleCoordVol(j);
            }

        }

        // compute convex hull - exitcode = 0 no error
        if (qh_new_qhull(2, interfacePoints.size(), coordinates, false, "qhull ", NULL, stderr)!=0)
            dserror(" error in the computation of the convex hull (qhull error)");

        if(((int) interfacePoints.size()) != qh num_vertices)
            dserror("resulting surface is concave - convex hull does not include all points");

        // copy vertices out of the facet list
        facetT* facet = qh facet_list;
        for(int i = 0; i< qh num_facets; i++)
        {
            for(int j = 0; j < 2; j++)
            {
                double* point  = SETelemt_(facet->vertices, j, vertexT)->point;
                for(int k = 0; k < 2; k++)
                    vertex[k] = point[k];
                {
                    static BlitzVec2  eleCoordSurf;
                    for(int m = 0; m < 2; m++)
                        eleCoordSurf(m)  = vertex[m];
                    static BlitzVec3 curCoordVol;
                    elementToCurrentCoordinates(cutterElement, xyze_cutterElement, eleCoordSurf, curCoordVol);
                    const BlitzVec3 eleCoordVol(currentToVolumeElementCoordinatesExact(xfemElement, curCoordVol));
                    for(int m = 0; m < 3; m++)
                        vertex[m] = eleCoordVol(m);
                }

                vertices.push_back(vertex);
            }
            facet = facet->next;
        }

        // free memory and clear vector of interface points
        qh_freeqhull(!qh_ALL);
        int curlong, totlong;           // memory remaining after qh_memfreeshort
        qh_memfreeshort (&curlong, &totlong);
        if (curlong || totlong)
            printf("qhull internal warning (main): did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);

        free(coordinates);

    }
    else if(interfacePoints.size() <= 2 && interfacePoints.size() > 0) // ??? == 1 ???
    {
        for(vector<InterfacePoint>::iterator ipoint = interfacePoints.begin(); ipoint != interfacePoints.end(); ++ipoint)
        {
            // transform interface points into current coordinates
            {
                static BlitzVec2  eleCoordSurf;
                for(int j = 0; j < 2; j++)
                    eleCoordSurf(j)  = ipoint->coord[j];
                static BlitzVec3 curCoordVol;
                elementToCurrentCoordinates(cutterElement, xyze_cutterElement, eleCoordSurf, curCoordVol);
                const BlitzVec3 eleCoordVol(currentToVolumeElementCoordinatesExact(xfemElement, curCoordVol));
                for(int j = 0; j < 3; j++)
                {
                    ipoint->coord[j] = eleCoordVol(j);
                    vertex[j] = eleCoordVol(j);
                }
            }
            vertices.push_back(vertex);
        }
    }
    else
        dserror("collection of interface points is empty");

    storePoint(vertices[0], interfacePoints, positions);
    vertices.erase(vertices.begin());

    if(interfacePoints.size() > 1)
    {
        // store points, segments and triangles for the computation of the
        // Constrained Delaunay Tetrahedralization with Tetgen
        searchPoint = vertices[0];
        storePoint(vertices[0], interfacePoints, positions );
        vertices.erase(vertices.begin());
    }


    while(vertices.size()>2)
    {
        findNextSegment(vertices, searchPoint);
        storePoint(searchPoint, interfacePoints, positions);
    }


    storeSurfacePoints(interfacePoints);

    // cutter element lies on the surface of an xfem element
    if(numInternalPoints == numBoundaryPoints && numInternalPoints != 0)
    {
        if(numBoundaryPoints > 1)
            storeSegments(positions);

    }
    else
    {
        if(interfacePoints.size() > 1)
            storeSegments( positions );

        if(interfacePoints.size() > 2)
        {
            pointList_.push_back(midpoint);
            storeTriangles(positions);
        }
    }
    interfacePoints.clear();

}
#endif //QHULL



/*----------------------------------------------------------------------*
 |  ICS:    finds the next facet of a convex hull            u.may 06/07|
 |          and returns the point different form the searchpoint        |
 *----------------------------------------------------------------------*/
void XFEM::Intersection::findNextSegment(
    vector< vector<double> >&   vertices,
    vector<double>&             searchPoint) const
{
    vector< vector<double> >::iterator it;
    bool pointfound = false;

    if(vertices.size()==0 || searchPoint.size()==0)
        dserror("one or both vectors are empty");

    for(it = vertices.begin(); it != vertices.end(); it=it+2 )
    {
        if(comparePoints<3>(searchPoint, *it))
        {
            pointfound = true;
            searchPoint = *(it+1);
            vertices.erase(it);
            vertices.erase(it); // remove it+ 1
            break;
        }

        if(comparePoints<3>(searchPoint, *(it+1)))
        {
            pointfound = true;
            searchPoint = *(it);
            vertices.erase(it);
            vertices.erase(it); // remove it+ 1
            break;
        }
    }
    if(!pointfound) dserror("no point found");
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
void XFEM::Intersection::computeCDT(
        const DRT::Element*                     xfemElement,
        const BlitzMat&                         xyze_xfemElement,
        const map<int,BlitzVec3>&               currentcutterpositions,
        map< int, DomainIntCells >&	            domainintcells,
        map< int, BoundaryIntCells >&           boundaryintcells,
        int                                     timestepcounter_)
{
    int dim = 3;
    int nsegments = 0;
    int nsurfPoints = 0;
    tetgenio in;
    tetgenio out;
    char switches[] = "pnnQ";    //o2 Y
    tetgenio::facet *f;
    tetgenio::polygon *p;


    // allocate pointlist
    in.numberofpoints = pointList_.size();
    in.pointlist = new REAL[in.numberofpoints * dim];

    // fill point list
    int fill = 0;
    for(int i = 0; i <  in.numberofpoints; i++)
        for(int j = 0; j < dim; j++)
        {
            in.pointlist[fill] = (REAL) pointList_[i].coord[j];
            fill++;
        }


    in.pointmarkerlist = new int[in.numberofpoints];
    for(int i = 0; i < numXFEMCornerNodes_; i++)
        in.pointmarkerlist[i] = 3;    // 3 : point not lying on the xfem element

    for(int i = numXFEMCornerNodes_; i < in.numberofpoints; i++)
        in.pointmarkerlist[i] = 2;    // 2 : point not lying on the xfem element


    in.numberoffacets = numXFEMSurfaces_ + triangleList_.size();

    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];


    // loop over all xfem element surfaces
    for(int i = 0; i < numXFEMSurfaces_; i++)
    {
        f = &in.facetlist[i];
        nsegments = (int) (segmentList_[i].size()/2);
        nsurfPoints = surfacePointList_[i].size();

        f->numberofpolygons = 1 + nsegments + nsurfPoints;
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];
        const int numnodequad4 = 4;
        p->numberofvertices = numnodequad4;
        p->vertexlist = new int[p->numberofvertices];
        for(int ivertex = 0; ivertex < numnodequad4; ivertex ++)
            p->vertexlist[ivertex] = eleNumberingSurfaces_[i][ivertex];


        int count = 0;
        for(int j = 1; j < 1 + nsegments; j ++)
        {
            if(segmentList_[i].size() > 0)
            {
                p = &f->polygonlist[j];
                p->numberofvertices = 2;
                p->vertexlist = new int[p->numberofvertices];

                for(int k = 0; k < 2; k ++)
                {
                   p->vertexlist[k] = segmentList_[i][count];
                   in.pointmarkerlist[segmentList_[i][count]] = 3;  // 3: point lying on the xfem boundary
                   count++;
                }
            }
        }

        count = 0;
        for(int j = 1 + nsegments; j < f->numberofpolygons; j++)
        {
            if(surfacePointList_[i].size() > 0)
            {
                p = &f->polygonlist[j];
                p->numberofvertices = 1;
                p->vertexlist = new int[p->numberofvertices];

                p->vertexlist[0] = surfacePointList_[i][count];
                in.pointmarkerlist[surfacePointList_[i][count]] = 3;  // 3: point lying on the xfem boundary
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


    //in.save_nodes("tetin");
    //in.save_poly("tetin");
    //  Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
    //  do quality mesh generation (q) with a specified quality bound
    //  (1.414), and apply a maximum volume constraint (a0.1)
    tetrahedralize(switches, &in, &out);


    //Debug
    //vector<int> elementIds;
    //for(int i = 388; i<389; i++)
    //    elementIds.push_back(element->Id());

    //debugTetgenOutput(in, out, xfemElement, elementIds, timestepcounter_);
    //printTetViewOutputPLC( element, element->Id(), in);

    // store interface triangles (+ recovery of higher order meshes)
    const bool higherorder = false;
    const bool recovery = false;

    if(higherorder)
    {
      cout << "lifting of Stinerpoints" << endl;
      dserror("hu");
    	recoverCurvedInterface(xfemElement, xyze_xfemElement, currentcutterpositions, boundaryintcells, out, recovery);
    }
   	else
   	{
   		storeIntCells(xfemElement, xyze_xfemElement, currentcutterpositions, boundaryintcells, out);
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
void XFEM::Intersection::startPointList(
    )
{
    InterfacePoint ip;

    pointList_.clear();

    for(int i = 0; i < numXFEMCornerNodes_; i++)
    {
        ip.nsurf = 3;

        // change for other element types
        for(int j = 0; j < 3; j++)
        {
           ip.coord[j]      =   eleRefCoordinates_[i][j];
           ip.surfaces[j]   =   eleNodesSurfaces_[i][j];
           ip.pType         =   NODE;
        }
        pointList_.push_back(ip);
    }

    for(int i = 0; i < numXFEMSurfaces_; i++)
        faceMarker_.push_back(-1);

}



/*----------------------------------------------------------------------*
 |  CDT:    stores a point within a list of points           u.may 06/07|
 |          which is to be copy to the tetgen data structure            |
 |          for the computation of the                                  |
 |          Constrained Delauney Triangulation                          |
 *----------------------------------------------------------------------*/
void XFEM::Intersection::storePoint(
    const vector<double>&             point,
    const vector<InterfacePoint>&     interfacePoints,
    vector<int>&                      positions)
{
    bool alreadyInList = false;

    for(vector<InterfacePoint>::const_iterator ipoint = interfacePoints.begin(); ipoint != interfacePoints.end(); ++ipoint )
    {
        if(comparePoints<3>(point, ipoint->coord))
        {
            alreadyInList = false;
            int count = 0;
            for(vector<InterfacePoint>::const_iterator it = pointList_.begin(); it != pointList_.end(); ++it )
            {
                if(comparePoints<3>(point, it->coord))
                {
                    alreadyInList = true;
                    break;
                }
                count++;
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
 |  CDT:    computes the midpoint of a collection of         u.may 06/07|
 |          InterfacePoints                                             |
 *----------------------------------------------------------------------*/
XFEM::InterfacePoint XFEM::Intersection::computeMidpoint(
    const vector<InterfacePoint>& interfacePoints
    ) const
{

     int n = interfacePoints.size();
     InterfacePoint ip;

     ip.nsurf = 0;

     for(int i = 0; i < 3 ; i++)
        ip.coord[i] = 0.0;

     for(int i = 0; i < n ; i++)
        for(int j = 0; j < 3 ; j++)
         ip.coord[j] += interfacePoints[i].coord[j];

     for(int i = 0; i < 3 ; i++)
       ip.coord[i] = ip.coord[i]/(double)n;

     return ip;
}



/*----------------------------------------------------------------------*
 |  CDT:    stores a single point lying on a surface of an   u.may 06/07|
 |          xfem element, if no segments are lying in                   |
 |          that surface                                                |
 *----------------------------------------------------------------------*/
void XFEM::Intersection::storeSurfacePoints(
    const vector<InterfacePoint>&     interfacePoints)
{

    for(unsigned int i = 0; i < interfacePoints.size(); i++)
    {
        bool singlePoint = true;
        if(interfacePoints[i].pType == SURFACE || interfacePoints[i].pType == LINE)
        {
            for(unsigned int j = 0; j < interfacePoints.size(); j++)
            {
                if((interfacePoints[j].pType != INTERNAL) && (i != j))
                {
                    for(int k = 0; k < interfacePoints[i].nsurf; k++)
                    {
                        for(int l = 0; l < interfacePoints[j].nsurf; l++)
                        {
                            const int surf1 = interfacePoints[i].surfaces[k];
                            const int surf2 = interfacePoints[j].surfaces[l];

                            if(surf1 == surf2)
                            {
                                singlePoint = false;
                                break;
                            }
                        }
                        if(!singlePoint)
                            break;
                    }
                }
                if(!singlePoint)
                        break;
            }
        }
        else
            singlePoint = false;


        if(singlePoint)
        {
            for(unsigned int jj = numXFEMCornerNodes_; jj < pointList_.size(); jj++ )
            {
                if(comparePoints<3>(interfacePoints[i].coord, pointList_[jj].coord))
                {
                    bool alreadyInList = false;
                    for(int kk = 0; kk < numXFEMSurfaces_; kk++)
                        for(unsigned int ll = 0; ll < surfacePointList_[kk].size(); ll++)
                            if(surfacePointList_[kk][ll] == (int) jj)
                            {
                                alreadyInList = true;
                                break;
                            }

                    if(!alreadyInList)
                        surfacePointList_[interfacePoints[i].surfaces[0]].push_back(jj);

                    break;
                }
            }
        }
    }
}



/*----------------------------------------------------------------------*
 |  CDT:    stores a segment within a list of segments       u.may 06/07|
 |          which is to be copied to the tetgen data structure          |
 |          for the computation of the                                  |
 |          Constrained Delauney Triangulation                          |
 *----------------------------------------------------------------------*/
void XFEM::Intersection::storeSegments(
    const vector<int>&              positions)
{

    for(unsigned int i = 0; i < positions.size(); i++ )
    {
        const int pos1 = positions[i];
        int pos2 = 0;
        if(pos1 ==  positions[positions.size()-1])
            pos2 = positions[0];
        else
            pos2 = positions[i+1];


        for(int j = 0; j < pointList_[pos1].nsurf; j++ )
            for(int k = 0; k < pointList_[pos2].nsurf; k++ )
            {
                const int surf1 = pointList_[pos1].surfaces[j];
                const int surf2 = pointList_[pos2].surfaces[k];

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


/*----------------------------------------------------------------------*
 |  CDT:    stores a triangle within a list of trianles      u.may 06/07|
 |          which is to be copy to the tetgen data structure            |
 |          for the computation of the                                  |
 |          Constrained Delauney Triangulation                          |
 *----------------------------------------------------------------------*/
void XFEM::Intersection::storeTriangles(
    const vector<int>               positions)
{
    vector<int> triangle(3,0);

    for(unsigned int i = 0; i < positions.size()-1; i++ )
    {
        triangle[0] = positions[i];
        triangle[1] = positions[i+1];
        triangle[2] = pointList_.size()-1;

        triangleList_.push_back(triangle);
        faceMarker_.push_back(intersectingCutterElements_.size()-1);
    }

    triangle[0] = positions[positions.size()-1];
    triangle[1] = positions[0];
    triangle[2] = pointList_.size()-1;

    triangleList_.push_back(triangle);
    faceMarker_.push_back(intersectingCutterElements_.size()-1);
}


/*----------------------------------------------------------------------*
 |  RCI:    stores a pointer to each intersecting            u.may 08/07|
 |          cutter element  used for the recovery of curved interface   |
 *----------------------------------------------------------------------*/
void XFEM::Intersection::storeIntersectedCutterElement(
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
void XFEM::Intersection::recoverCurvedInterface(
        const DRT::Element*             xfemElement,
        const BlitzMat&                 xyze_xfemElement,
        const map<int,BlitzVec3>&       currentcutterpositions,
        map< int, BoundaryIntCells >&   boundaryintcells,
        tetgenio&                       out,
        bool							recovery
        )
{
    BoundaryIntCells                     			listBoundaryICPerElement;

    // list of point markers , if already visited = 1 , if not = 0
    vector<int> visitedPointIndexList(out.numberofpoints,0);
//    for(int i = 0; i<out.numberofpoints; ++i)
//        visitedPointIndexList[i] = 0;

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
            vector<BlitzVec>        tetraCornerNodes(4, BlitzVec(3));
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
 |                               										|
 *----------------------------------------------------------------------*/
void XFEM::Intersection::storeIntCells(
		const DRT::Element*             xfemElement,
		const BlitzMat&                 xyze_xfemElement,
		const map<int,BlitzVec3>&       currentcutterpositions,
		map< int, BoundaryIntCells >&   boundaryintcells,
		tetgenio&                       out)
{

	BoundaryIntCells                     			listBoundaryICPerElement;

    // lifts all corner points into the curved interface
    //liftAllSteinerPoints(xfemElement, out);

    for(int i=0; i<out.numberoftrifaces; i++)
    {
        // run over all faces not lying in on of the xfem element planes
        const int faceMarker = out.trifacemarkerlist[i] - facetMarkerOffset_;
        vector<vector<double> > domainCoord(3, std::vector<double>(3,0.0));
        vector<vector<double> > boundaryCoord(3, std::vector<double>(3,0.0));

        if(faceMarker > -1)
        {
            //const int tetIndex = out.adjtetlist[i*2];
            //printf("tetIndex = %d\n", tetIndex);
            //vector<int>             order(3,0);
            //vector<int>             tetraCornerIndices(4,0);
            //vector<BlitzVec>        tetraCornerNodes(4, BlitzVec(3));
            //getTetrahedronInformation(tetIndex, i, tetraCornerIndices, order, out );
            //getTetrahedronNodes(tetraCornerNodes, tetraCornerIndices, xfemElement, out);

            // run over each triface
            for(int index1 = 0; index1 < 3 ;index1++)
            {
                //int index2 = index1+1;
                //if(index2 > 2) index2 = 0;

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
 |  RCI:    checks if all tetrahedra corner points are lying u.may 09/07|
 |          in a surface element                                        |
 |          if not corner points is recovered on the surface element    |
 *----------------------------------------------------------------------*/
void XFEM::Intersection::liftAllSteinerPoints(
        const DRT::Element*                             xfemElement,
        const BlitzMat&                                 xyze_xfemElement,
        const map<int,BlitzVec3>&                       currentcutterpositions,
        tetgenio&                                       out)
{
    BlitzVec edgePoint(3);      edgePoint = 0.0;
    BlitzVec oppositePoint(3);  oppositePoint = 0.0;
    vector< vector<int> > adjacentFacesList;
    vector< vector<int> > adjacentFacemarkerList;

    locateSteinerPoints(adjacentFacesList, adjacentFacemarkerList, out);

    if(!adjacentFacesList.empty())
    {
        // run over all Steiner points
        for(unsigned int i=0; i<adjacentFacesList.size(); i++)
        {
            int lineIndex = -1, cutterIndex = -1;
            const int  caseSteiner = decideSteinerCase(  i, lineIndex, cutterIndex, adjacentFacesList, adjacentFacemarkerList, currentcutterpositions,
                                                   edgePoint, oppositePoint, xfemElement, xyze_xfemElement, out);
            switch(caseSteiner)
            {
                case 1:
                {
                    liftSteinerPointOnSurface(i, adjacentFacesList, adjacentFacemarkerList, currentcutterpositions, xfemElement, xyze_xfemElement, out);
                    break;
                }
                case 2:
                {
                    liftSteinerPointOnEdge( i, lineIndex, cutterIndex, edgePoint, oppositePoint,
                                        adjacentFacesList, currentcutterpositions, xfemElement, xyze_xfemElement, out);
                    break;
                }
                case 3:
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
void XFEM::Intersection::locateSteinerPoints(
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
int XFEM::Intersection::decideSteinerCase(
        const int                       steinerIndex,
        int&                            lineIndex,
        int&                            cutterIndex,
        const vector< vector<int> >&   adjacentFacesList,
        const vector< vector<int> >&   adjacentFacemarkerList,
        const map<int,BlitzVec3>&       currentcutterpositions,
        BlitzVec&                       edgePoint,
        BlitzVec&                       oppositePoint,
        const DRT::Element*             xfemElement,
        const BlitzMat&                 xyze_xfemElement,
        const tetgenio&                 out
        ) const
{
    const int pointIndex = adjacentFacesList[steinerIndex][0];

    static BlitzVec3    x;
    for(int k=0; k<3; k++)
        x(k)   = out.pointlist[pointIndex*3 + k];

    static BlitzVec3    xsi;
    currentToVolumeElementCoordinates(xfemElement, x, xsi);

    InterfacePoint emptyIp;
    if(setInterfacePointBoundaryStatus(xfemElement-> Shape(), xsi, emptyIp))
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

    int caseSteiner = 0;
    if(normalSteiner)
        caseSteiner = 1;
    else
        caseSteiner = 2;
    if(out.pointmarkerlist[pointIndex] == 3)
        caseSteiner = 3;

    return caseSteiner;
}



/*----------------------------------------------------------------------*
 |  RCI:    lifts Steiner points lying                       u.may 11/07|
 |          within a cutter element                                     |
 *----------------------------------------------------------------------*/
void XFEM::Intersection::liftSteinerPointOnSurface(
        const int                       steinerIndex,
        const vector< vector<int> >&    adjacentFacesList,
        const vector< vector<int> >&    adjacentFacemarkerList,
        const map<int,BlitzVec3>&       currentcutterpositions,
        const DRT::Element*             xfemElement,
        const BlitzMat&                 xyze_xfemElement,
        tetgenio&                       out
        )
{
    // get Steiner point coordinates
    BlitzVec    Steinerpoint(3);
    for(int j=0; j<3; ++j)
    {
        Steinerpoint(j) = out.pointlist[adjacentFacesList[steinerIndex][0]*3 + j];
    }
    elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement, Steinerpoint);

    BlitzVec    averageNormal(3);
    averageNormal = 0.0;

    const int length = (int) ( ( (double) (adjacentFacesList[steinerIndex].size()-1))*0.5 );
    vector<BlitzVec> normals;
    normals.reserve(length);

    for(int j = 0; j < length; ++j)
    {
        const int pointIndex1 = adjacentFacesList[steinerIndex][1 + 2*j];
        const int pointIndex2 = adjacentFacesList[steinerIndex][1 + 2*j + 1];

        BlitzVec    p1(3);
        BlitzVec    p2(3);
        for(int k = 0; k < 3; ++k)
        {
            p1(k) =  out.pointlist[pointIndex1*3 + k];
            p2(k) =  out.pointlist[pointIndex2*3 + k];
        }
        elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement, p1);
        elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement, p2);

        const BlitzVec n1(p1 - Steinerpoint);
        const BlitzVec n2(p2 - Steinerpoint);

        BlitzVec normal(computeCrossProduct( n1, n2));
        normalizeVectorInPLace(normal);

        averageNormal += normal;

        normals.push_back(normal);
    }

    // compute average normal
    averageNormal /= ((double)length);

    const int faceMarker = adjacentFacemarkerList[steinerIndex][0];
    DRT::Element* cutterElement =  intersectingCutterElements_[faceMarker];
    const BlitzMat xyze_cutterElement(getCurrentNodalPositions(cutterElement, currentcutterpositions));

    static BlitzVec3 xsi;
    vector<BlitzVec> plane;
    plane.push_back(BlitzVec(Steinerpoint + averageNormal));
    plane.push_back(BlitzVec(Steinerpoint - averageNormal));
    const bool intersected = computeRecoveryNormal( xsi, plane, cutterElement, xyze_cutterElement, false);
    if(intersected)
    {
        storeHigherOrderNode(   true, adjacentFacesList[steinerIndex][0], -1, xsi,
            cutterElement, currentcutterpositions, xfemElement, out);
    }
    else
    {
        // loop over all individual normals
        bool intersected = false;
        vector<BlitzVec>::const_iterator normalptr;
        for(normalptr = normals.begin(); normalptr != normals.end(); ++normalptr )
        {
            vector<BlitzVec> plane;
            plane.push_back(BlitzVec(Steinerpoint + (*normalptr)));
            plane.push_back(BlitzVec(Steinerpoint - (*normalptr)));
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
            printf("STEINER POINT NOT LIFTED in liftSteinerPointOnSurface()\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        }
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    lifts Steiner points lying                       u.may 11/07|
 |          on the edge of a  cutter element                            |
 *----------------------------------------------------------------------*/
void XFEM::Intersection::liftSteinerPointOnEdge(
    const int                         steinerIndex,
    int                               lineIndex,     // TODO:??????????????? why is this an argument to this function?
    const int                         cutterIndex,
    BlitzVec&                         edgePoint,
    BlitzVec&                         oppositePoint,
    const vector< vector<int> >&     adjacentFacesList,
    const map<int,BlitzVec3>&         currentcutterpositions,
    const DRT::Element*               xfemElement,
    const BlitzMat&                   xyze_xfemElement,
    tetgenio&                         out)
{
    // get Steiner point coordinates
    BlitzVec    Steinerpoint(3);
    for(int j=0; j<3; j++)
        Steinerpoint(j) = out.pointlist[adjacentFacesList[steinerIndex][0]*3 + j];

    elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement, Steinerpoint);
    elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement, edgePoint);
    elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement, oppositePoint);

    const BlitzVec r1(edgePoint - Steinerpoint);
    const BlitzVec r2(oppositePoint - Steinerpoint);

    BlitzVec n1(computeCrossProduct( r1, r2));
    BlitzVec n2(computeCrossProduct( r1, n1));

    normalizeVectorInPLace(n1);
    normalizeVectorInPLace(n2);

    vector<BlitzVec> plane;
    plane.push_back(BlitzVec(Steinerpoint + n1));
    plane.push_back(BlitzVec(Steinerpoint - n1));
    plane.push_back(BlitzVec(plane[1] + n2));
    plane.push_back(BlitzVec(plane[0] + n2));

    static BlitzVec3 xsi;
    xsi = 0.0;
    DRT::Element* cutterElement = intersectingCutterElements_[cutterIndex];
    const BlitzMat xyze_cutterElement(getCurrentNodalPositions(cutterElement, currentcutterpositions));
    const bool intersected = computeRecoveryPlane( lineIndex, currentcutterpositions, xsi, plane, cutterElement, xyze_cutterElement);

    if(intersected)
    {
        storeHigherOrderNode(   false, adjacentFacesList[steinerIndex][0], lineIndex, xsi,
            cutterElement, currentcutterpositions, xfemElement, out);
    }
    else
    {
        countMissedPoints_++;
        printf("STEINER POINT NOT LIFTED in liftSteinerPointOnEdge()\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    }
}


/*----------------------------------------------------------------------*
 |  RCI:    lifts Steiner points lying                       u.may 11/07|
 |          on the boundary of the xfem element                         |
 *----------------------------------------------------------------------*/
void XFEM::Intersection::liftSteinerPointOnBoundary(
        const int                         steinerIndex,
        const vector< vector<int> >&      adjacentFacesList,
        const vector< vector<int> >&      adjacentFacemarkerList,
        const map<int,BlitzVec3>&         currentcutterpositions,
        const DRT::Element*               xfemElement,
        const BlitzMat&                   xyze_xfemElement,
        tetgenio&                         out
        )
{
    int edgeIndex = 0;
    int oppositeIndex = 0;
    int facemarkerIndex = 0;

    // find egde on boundary
    bool edgeFound = false;
    for(unsigned int i = 1; i < adjacentFacesList[steinerIndex].size(); i++ )
    {
        edgeIndex = adjacentFacesList[steinerIndex][i];
        if(out.pointmarkerlist[edgeIndex] == 3)
        {
            edgeFound = true;  // TODO:  why is this calculated?
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
    vector<BlitzVec> plane;
    computeIntersectionNormalC( adjacentFacesList[steinerIndex][0], edgeIndex, oppositeIndex, plane,
                               xfemElement, out);

    // compute intersection normal on boundary
    static BlitzVec3 xsi;
    DRT::Element* cutterElement = intersectingCutterElements_[faceIndex];
    const BlitzMat xyze_cutterElement(getCurrentNodalPositions(cutterElement, currentcutterpositions));
    const bool intersected = computeRecoveryNormal(xsi, plane, cutterElement, xyze_cutterElement, true);

    if(intersected)
    {
        storeHigherOrderNode(   true, adjacentFacesList[steinerIndex][0], -1, xsi,
                                intersectingCutterElements_[faceIndex], currentcutterpositions, xfemElement, out);
    }
    else
    {
        countMissedPoints_++;
        printf("STEINER POINT NOT LIFTED in liftSteinerPointOnBoundary()\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    }

}



/*----------------------------------------------------------------------*
 |  RCI:    returns information of the tetrahedra            u.may 08/07|
 *----------------------------------------------------------------------*/
void XFEM::Intersection::getTetrahedronInformation(
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
void XFEM::Intersection::getTetrahedronNodes(
        vector<BlitzVec>&                       tetraCornerNodes,
        const vector<int>&                      tetraCornerIndices,
        const DRT::Element*                     xfemElement,
        const BlitzMat&                         xyze_xfemElement,
        const tetgenio&                         out
        ) const
{

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 3; j++)
            tetraCornerNodes[i](j) = out.pointlist[tetraCornerIndices[i]*3+j];

        elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement, tetraCornerNodes[i]);
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    lifts the higher-order point of an edge of the   u.may 09/07|
 |          linearized interface onto the curved interface              |
 *----------------------------------------------------------------------*/
void XFEM::Intersection::computeHigherOrderPoint(
        const int                                 index1,
        const int                                 index2,
        const int                                 faceIndex,
        const int                                 faceMarker,
        const int                                 globalHigherOrderIndex,
        const vector<int>&                        tetraCornerIndices,
        const vector<BlitzVec>&                   tetraCornerNodes,
        const DRT::Element*                       xfemElement,
        const BlitzMat&                           xyze_xfemElement,
        const map<int,BlitzVec3>&                 currentcutterpositions,
        tetgenio&                           out
        )
{

    bool                                    intersected             = false;
    bool                                    intersectionNormal      = true;
    int                                     lineIndex               = -1;
    int                                     adjacentFaceMarker      = -1;
    int                                     adjacentFaceIndex       = -1;
    static BlitzVec3 xsi;
    xsi = 0.0;

    findAdjacentFace(  tetraCornerIndices[index1], tetraCornerIndices[index2],
                       faceMarker, adjacentFaceMarker, faceIndex, adjacentFaceIndex, out);

    // TODO: something seems weird about the following if else construct, sometimes it goes through, somethimes not...
    // suggestion: for every if, there should be an else, otherwise, values could become undefined...

    // edge lies within the xfem element
    if(adjacentFaceMarker  > -1)
    {
        vector<BlitzVec>        plane;
        computeIntersectionNormalB(  tetraCornerIndices[index1], tetraCornerIndices[index2], faceIndex,
                                    adjacentFaceIndex, globalHigherOrderIndex, plane, xfemElement, out);

        // higher order node lies within the cutter element
        if(adjacentFaceMarker == faceMarker)
        {
            DRT::Element* cutterElement = intersectingCutterElements_[faceMarker];
            const BlitzMat xyze_cutterElement(getCurrentNodalPositions(cutterElement, currentcutterpositions));
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
                const BlitzMat xyze_cutterElement(getCurrentNodalPositions(cutterElement, currentcutterpositions));
                intersected = computeRecoveryPlane( lineIndex, currentcutterpositions, xsi, plane, cutterElement, xyze_cutterElement);
                intersectionNormal = false;
            }
            else
            {
                dserror("what do we do here? Implement something?");
            }
        }
        else
        {
            dserror("should we ever arive here?");
        }
    }
    // edge lies on the surface of the xfem element
    else if(adjacentFaceMarker  == -1)
    {
        const int oppositeIndex = findEdgeOppositeIndex(  tetraCornerIndices[index1], tetraCornerIndices[index2],
                                                    adjacentFaceIndex, out);

        //printf("oppo = %d\n", oppositeIndex);

        DRT::Element* cutterElement = intersectingCutterElements_[faceMarker];
        const BlitzMat xyze_cutterElement(getCurrentNodalPositions(cutterElement, currentcutterpositions));
        
        vector<BlitzVec>     plane;

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
            else
            {
                dserror("What do we do here?");
            }
        }
    }
    else // (adjacentFaceMarker < -1) should be impossible
    {
        dserror("bug in adjacentFaceMarker numbering?");
    }


    if(intersected)
    {
        storeHigherOrderNode(   intersectionNormal, globalHigherOrderIndex, lineIndex,
                                xsi, intersectingCutterElements_[faceMarker], currentcutterpositions, xfemElement, out);
    }
    else
    {
        countMissedPoints_++;
        cout << "faceMarker = " << faceMarker << endl;
        dserror("NO INTERSECTION POINT FOUND!!!!! adjacentFaceMarker = %d\n", adjacentFaceMarker);
    }

}



/*----------------------------------------------------------------------*
 |  RCI:    returns the other two point indices belonging    u.may 09/07|
 |          to a triface that obtains a Steiner point                   |
 *----------------------------------------------------------------------*/
vector<int> XFEM::Intersection::getPointIndices(
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
bool XFEM::Intersection::computeRecoveryNormal(
    BlitzVec3&                                  xsi,
    const vector<BlitzVec>&                     normal,
    const DRT::Element*                         cutterElement,
    const BlitzMat&                             xyze_cutterElement,
    const bool                                  onBoundary
    ) const
{
    bool                        intersection = true;
    int                         iter = 0;
    int                         countSingular = 0;
    const int                   maxiter = 50;
    double                      residual = 1.0;
    static BlitzMat3x3 A;
    static BlitzVec3   b;
    static BlitzVec3   dx;
    b = 0.0;
    dx = 0.0;

    xsi = 0.0;
    updateRHSForRCINormal( b, xsi, normal, cutterElement, xyze_cutterElement, onBoundary);

    while(residual > TOL14)
    {
        updateAForRCINormal( A, xsi, normal, cutterElement, xyze_cutterElement, onBoundary);

        if(!solveLinearSystemWithSVD<3>(A, b, dx))
            countSingular++;

        if(countSingular > 5)
        {
            intersection = false;
            break;
        }

        xsi += dx;
        //printf("dx0 = %20.16f\t, dx1 = %20.16f\t, dx2 = %20.16f\n", dx(0), dx(1), dx(2));
        if(iter >= maxiter)
        {
            intersection = false;
            break;
        }

        updateRHSForRCINormal( b, xsi, normal, cutterElement, xyze_cutterElement, onBoundary);
        residual = Norm2(b);
        iter++;

        //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi(0), xsi(1), xsi(2), residual, TOL14);
    }

    if( (fabs(xsi(0))-1.0) > TOL7  || (fabs(xsi(1))-1.0) > TOL7 )    // line coordinate may be bigger than 1
    {
        //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi(0), xsi(1), xsi(2), residual, TOL14);
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
void XFEM::Intersection::updateAForRCINormal(
    BlitzMat3x3&                                A,
    const BlitzVec3&                            xsi,
    const vector<BlitzVec>&                     normal,
    const DRT::Element*                         surfaceElement,
    const BlitzMat&                             xyze_surfaceElement,
    const bool                                  onBoundary
    ) const
{
    const int numNodesSurface = surfaceElement->NumNode();

    A = 0.0;
    BlitzMat surfaceDeriv1(2,numNodesSurface);
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
        BlitzMat lineDeriv1(1,numNodesLine);
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
void XFEM::Intersection::updateRHSForRCINormal(
    BlitzVec3&                  b,
    const BlitzVec3&            xsi,
    const vector<BlitzVec>&     normal,
    const DRT::Element*                         surfaceElement,
    const BlitzMat&                             xyze_surfaceElement,
    const bool                                  onBoundary
    ) const
{
    const int numNodesSurface = surfaceElement->NumNode();
    BlitzVec surfaceFunct(numNodesSurface);
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
        BlitzVec lineFunct(numNodesLine);
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
bool XFEM::Intersection::computeRecoveryPlane(
        int&                                        lineIndex,
        const map<int,BlitzVec3>&                   currentcutterpositions,
        BlitzVec3&                                  xsi,
        const vector<BlitzVec>&                     plane,
        DRT::Element*                               cutterElement,
        const BlitzMat&                             xyze_cutterElement
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
        const int                   maxiter = 50;
        double                      residual = 1.0;
        DRT::Element*               lineElement = cutterElement->Lines()[i];
        const BlitzMat xyze_lineElement(getCurrentNodalPositions(lineElement, currentcutterpositions));
        
        static BlitzMat3x3 A;
        static BlitzVec3   b;
        static BlitzVec3   dx;  dx = 0.0;

        intersection = true;
        xsi = 0.0;

        updateRHSForRCIPlane( b, xsi, plane, lineElement, xyze_lineElement);

        while( residual > TOL14 )
        {
            updateAForRCIPlane( A, xsi, plane, lineElement, xyze_lineElement, cutterElement, xyze_cutterElement);

            if(!gaussElimination<true, 3, 1>(A, b, dx))
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
            residual = Norm2(b);
            iter++;
        }

        if( (fabs(xsi(2))-1.0) > TOL7 )     // planes coordinate may be bigger than 1
        {   printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi(0), xsi(1), xsi(2), residual, TOL14);
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
void XFEM::Intersection::updateAForRCIPlane(
        BlitzMat3x3&                A,
        const BlitzVec3&            xsi,
        const vector<BlitzVec>&     plane,
        const DRT::Element*         lineElement,
        const BlitzMat&             xyze_lineElement,
        const DRT::Element*         surfaceElement,
        const BlitzMat&             xyze_surfaceElement
        ) const
{
    const int numNodesLine = lineElement->NumNode();
    const int numNodesSurface = 4;

    BlitzMat surfaceDeriv(2,numNodesSurface);
    DRT::UTILS::shape_function_2D_deriv1(surfaceDeriv, xsi(0),  xsi(1), DRT::Element::quad4);
    BlitzMat lineDeriv(1,numNodesLine);
    DRT::UTILS::shape_function_1D_deriv1(lineDeriv, xsi(2), lineElement->Shape());

    dsassert((int)plane.size() >= numNodesSurface, "plane array has to have size numNodesSurface ( = 4)!");

    A = 0.0;
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
void XFEM::Intersection::updateRHSForRCIPlane(
    BlitzVec3&                  b,
    const BlitzVec3&            xsi,
    const vector<BlitzVec>&     plane,
    const DRT::Element*         lineElement,
    const BlitzMat&             xyze_lineElement
    ) const
{
    const int numNodesLine    = lineElement->NumNode();
    const int numNodesSurface = 4;

    BlitzVec surfaceFunct(numNodesSurface);
    DRT::UTILS::shape_function_2D(surfaceFunct, xsi(0), xsi(1), DRT::Element::quad4 );
    BlitzVec lineFunct(numNodesLine);
    DRT::UTILS::shape_function_1D(lineFunct, xsi(2), lineElement->Shape());

    dsassert((int)plane.size() >= numNodesSurface, "plane array has to have size numNodesSurface ( = 4)!");

    b = 0.0;
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
void XFEM::Intersection::computeIntersectionNormalA(
    const bool                              onBoundary,
    const int                               index1,
    const int                               index2,
    const int                               oppositePointIndex,
    const int                               globalHigherOrderIndex,
    const vector<int>&                      tetraCornerIndices,
    const vector<BlitzVec>&                 tetraCornerNodes,
    vector<BlitzVec>&                       plane,
    const DRT::Element*                     xfemElement,
    const tetgenio&                         out) const
{

    BlitzVec  p1(3);
    BlitzVec  p2(3);
    BlitzVec  p3(3);


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
    const BlitzVec  r1(p1 - p2);
    const BlitzVec  r2(p3 - p2);

    // normal of the plane
    BlitzVec  n(computeCrossProduct(r1, r2));
    normalizeVectorInPLace(n);

    // direction vector of the intersection line
    BlitzVec  r(computeCrossProduct(n, r2));
    normalizeVectorInPLace(r);

    // computes the start point of the line
    BlitzVec  m(3);

    if(!onBoundary)
        m = computeLineMidpoint(p2, p3);
    else
    {
        for(int i = 0; i < 3; i++)
            m(i) = out.pointlist[globalHigherOrderIndex*3+i];
    }

    // compute nodes of the normal to the interface edge of the tetrahedron
    plane.clear();
    plane.reserve(5);
    plane.push_back(BlitzVec(m + r));
    plane.push_back(BlitzVec(m - r));
    plane.push_back(BlitzVec(plane[1] + n));
    plane.push_back(BlitzVec(plane[0] + n));


    if(onBoundary)
    {
        for(int i = 0; i < 4; i++)
            elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement_, plane[i]);

        elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement_, m);
        plane.push_back(BlitzVec(m));
    }
}



/*----------------------------------------------------------------------*
 |  RCI:    computes the normal to the interface edge of     u.may 11/07|
 |          two adjacent triangular faces                               |
 *----------------------------------------------------------------------*/
void XFEM::Intersection::computeIntersectionNormalB(
        const int                                 index1,
        const int                                 index2,
        const int                                 faceIndex,
        const int                                 adjacentFaceIndex,
        const int                                 globalHigherOrderIndex,
        vector<BlitzVec>&                         plane,
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
    BlitzVec p1(3);
    BlitzVec p2(3);
    BlitzVec p3(3);
    BlitzVec p4(3);

    for(int i=0; i<3; i++)
    {
        p1(i) = out.pointlist[index1*3 + i];
        p2(i) = out.pointlist[index2*3 + i];
        p3(i) = out.pointlist[oppositePointIndex*3 + i];
        p4(i) = out.pointlist[adjacentOppositePointIndex*3 + i];
    }

    elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement_, p1);
    elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement_, p2);
    elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement_, p3);
    elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement_, p4);

    const BlitzVec r1(p1 - p2);
    const BlitzVec r2(p3 - p2);
    const BlitzVec r3(p4 - p2);

    BlitzVec n1(computeCrossProduct(r2, r1));
    BlitzVec n2(computeCrossProduct(r1, r3));

    BlitzVec averageNormal(n1 + n2);
    BlitzVec rPlane(computeCrossProduct(n1, r1));

//    for(int i = 0; i < 3; i++)
//        averageNormal(i) = 0.5*averageNormal(i);
    averageNormal *= 0.5;

    normalizeVectorInPLace(averageNormal);
    normalizeVectorInPLace(rPlane);

    BlitzVec m(3);
    for(int i = 0; i < 3; i++)
        m(i) = out.pointlist[globalHigherOrderIndex*3+i];

    elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement_, m);

    // compute nodes of the normal to the interface edge of the tetrahedron
    plane.clear();
    plane.reserve(4);
    plane.push_back(BlitzVec(m + averageNormal));
    plane.push_back(BlitzVec(m - averageNormal));
    plane.push_back(BlitzVec(plane[1] + rPlane));
    plane.push_back(BlitzVec(plane[0] + rPlane));

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
void XFEM::Intersection::computeIntersectionNormalC(
    const int                               steinerIndex,
    const int                               edgeIndex,
    const int                               oppositeIndex,
    vector<BlitzVec>&                       plane,
    const DRT::Element*                     xfemElement,
    const tetgenio&                         out) const
{

    BlitzVec  p1(3);
    BlitzVec  p2(3);
    BlitzVec  p3(3);


    for(int i=0; i<3; i++)
    {
        p1(i) = out.pointlist[oppositeIndex*3  + i];
        p2(i) = out.pointlist[steinerIndex*3   + i];
        p3(i) = out.pointlist[edgeIndex*3      + i];
    }

    // compute direction vectors of the plane
    BlitzVec  r1(p1 - p2);
    BlitzVec  r2(p3 - p2);

    // normal of the plane
    BlitzVec  n(computeCrossProduct(r1, r2));
    normalizeVectorInPLace(n);

    // direction vector of the intersection line
    BlitzVec  r(computeCrossProduct(n, r2));
    normalizeVectorInPLace(r);

    // compute nodes of the normal to the interface edge of the tetrahedron
    plane.clear();
    plane.reserve(5);
    plane.push_back(BlitzVec(p2 + r));
    plane.push_back(BlitzVec(p2 - r));
    plane.push_back(BlitzVec(plane[1] + n));
    plane.push_back(BlitzVec(plane[0] + n));

    for(int i = 0; i < 4; i++)
        elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement_, plane[i]);

    elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement_, p2);
    plane.push_back(BlitzVec(p2));
}



/*----------------------------------------------------------------------*
 |  RCI:    computes the midpoint of a line                  u.may 08/07|
 *----------------------------------------------------------------------*/
BlitzVec XFEM::Intersection::computeLineMidpoint(
    const BlitzVec& p1,
    const BlitzVec& p2) const
{
    BlitzVec midpoint(3);

    for(int i=0; i<3; i++)
        midpoint(i) = (p1(i) + p2(i))*0.5;
//    midpoint = 0.5*(p1 + p2);

    return midpoint;
}



/*----------------------------------------------------------------------*
 |  RCI:    searches for the face marker                     u.may 10/07|
 |          of a facet adjacent to a given edge of                      |
 |          of a given facet                                            |
 *----------------------------------------------------------------------*/
void XFEM::Intersection::findAdjacentFace(
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
int XFEM::Intersection::findEdgeOppositeIndex(
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
bool XFEM::Intersection::findCommonFaceEdge(
        const int                           faceIndex1,
        const int                           faceIndex2,
        const vector<int>&                  adjacentFacesList,
        BlitzVec&                           edgePoint,
        BlitzVec&                           oppositePoint,
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
bool XFEM::Intersection::findCommonCutterLine(
    const map<int,BlitzVec3>&                       currentcutterpositions,
    const int                                       faceIndex1,
    const int                                       faceIndex2,
    int&                                            lineIndex,
    int&                                            cutterIndex
    ) const
{
    bool comparison = false;
    // get line arrays, which are computed onthe fly within the cutter elements
    const DRT::Element*const* lines1 = intersectingCutterElements_[faceIndex1]->Lines();
    const DRT::Element*const* lines2 = intersectingCutterElements_[faceIndex2]->Lines();

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
                BlitzVec3 pos1 = currentcutterpositions.find(node1->Id())->second;
                BlitzVec3 pos2 = currentcutterpositions.find(node2->Id())->second;
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
                        BlitzVec3 pos1 = currentcutterpositions.find(node1->Id())->second;
                        BlitzVec3 pos2 = currentcutterpositions.find(node2->Id())->second;
                        if(!comparePoints<3>(pos1, pos2))
                            comparison = false;
                    }
                    else
                    {
                        const DRT::Node* node1 = lines1[i]->Nodes()[k];
                        const DRT::Node* node2 = lines2[j]->Nodes()[1-k];
                        BlitzVec3 pos1 = currentcutterpositions.find(node1->Id())->second;
                        BlitzVec3 pos2 = currentcutterpositions.find(node2->Id())->second;
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
int XFEM::Intersection::findIntersectingSurfaceEdge(
        const DRT::Element*                       xfemElement,
        const BlitzMat&                           xyze_xfemElement,
        DRT::Element*                             cutterElement,
        const map<int,BlitzVec3>&                 currentcutterpositions,
        const BlitzVec&                           edgeNode1,
        const BlitzVec&                           edgeNode2
        ) const
{
    dserror("to be improved by Ursula");
    int lineIndex = -1;
    static BlitzVec3 x1;
    static BlitzVec3 x2;

    BlitzVec node1 = edgeNode1;
    BlitzVec node2 = edgeNode2;

    elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement, node1);
    elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement, node2);

    x1(0) = node1(0);
    x2(0) = node2(0);

    DRT::Element** lines = cutterElement->Lines();
//    for(int i = 0; i < cutterElement->NumLine(); i++)
//    {
//        const DRT::Element*  lineElement = lines[i];
//
//        BlitzVec xsi1(1);
//        BlitzVec xsi2(1);
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
void XFEM::Intersection::storeHigherOrderNode(
    const bool                                  normal,
    const int                                   globalHigherOrderIndex,
    const int                                   lineIndex,
    BlitzVec3&                                  xsi,
    DRT::Element*                               cutterElement,
    const map<int,BlitzVec3>&                   currentcutterpositions,
    const DRT::Element*                         xfemElement,
    tetgenio&                                   out
    ) const
{
    static BlitzVec3 curr;

    if(normal)
    {
        BlitzVec2 xsiSurf;
        xsiSurf(0) = xsi(0);
        xsiSurf(1) = xsi(1);
        const BlitzMat xyze_cutterElement(getCurrentNodalPositions(cutterElement, currentcutterpositions));
        elementToCurrentCoordinates(cutterElement, xyze_cutterElement, xsiSurf, curr);
    }
    else
    {
        BlitzVec1 xsiLine;
        xsiLine(0) = xsi(2);
        const DRT::Element* lineele = cutterElement->Lines()[lineIndex];
        const BlitzMat xyze_lineElement(getCurrentNodalPositions(lineele, currentcutterpositions));
        elementToCurrentCoordinates(lineele, xyze_lineElement, xsiLine, curr);
    }
    xsi = currentToVolumeElementCoordinatesExact(xfemElement, curr);

    //printf("xsiold0 = %20.16f\t, xsiold1 = %20.16f\t, xsiold2 = %20.16f\n", out.pointlist[index*3], out.pointlist[index*3+1], out.pointlist[index*3+2]);

    for(int i = 0; i < 3; i++)
        out.pointlist[globalHigherOrderIndex*3+i]   = xsi(i);

    //printf("xsi0    = %20.16f\t, xsi1    = %20.16f\t, xsi2    = %20.16f\n", xsi(0), xsi(1), xsi(2));
    //printf("\n");
}



/*----------------------------------------------------------------------*
 |  RCI:    stores domain integration cells                  u.may 11/07|
 *----------------------------------------------------------------------*/
void XFEM::Intersection::addCellsToDomainIntCellsMap(
        const DRT::Element*             xfemElement,
        map< int, DomainIntCells >&     domainintcells,
        const tetgenio&                 out,
        bool 							higherorder
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
void XFEM::Intersection::addCellsToBoundaryIntCellsMap(
    const int                         		trifaceIndex,
    const int                         		cornerIndex,
    const int                         		globalHigherOrderIndex,
    const int                         		faceMarker,
    const map<int,BlitzVec3>&             currentcutterpositions,
    vector<vector<double> >&              domainCoord,
    vector<vector<double> >&              boundaryCoord,
    const DRT::Element*						        xfemElement,
    const BlitzMat&                       xyze_xfemElement,
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
    static BlitzVec3 eleCoordDomainCorner;
    for(int k = 0; k < 3; k++)
        eleCoordDomainCorner(k) = out.pointlist[(out.trifacelist[trifaceIndex*3+cornerIndex])*3+k];

    domainCoord[cornerIndex][0] = eleCoordDomainCorner(0);
    domainCoord[cornerIndex][1] = eleCoordDomainCorner(1);
    domainCoord[cornerIndex][2] = eleCoordDomainCorner(2);

    static BlitzVec3 physCoordCorner;
    elementToCurrentCoordinates(xfemElement, xyze_xfemElement, eleCoordDomainCorner, physCoordCorner);

    //domainCoord.push_back(trinodes);

    const DRT::Element* cutterElement = intersectingCutterElements_[faceMarker];
    const BlitzMat xyze_cutterElement(getCurrentNodalPositions(cutterElement, currentcutterpositions));
    
    BlitzVec2 eleCoordBoundaryCorner;
    CurrentToSurfaceElementCoordinates(cutterElement, xyze_cutterElement, physCoordCorner, eleCoordBoundaryCorner);

//    cout << "physcood = " << physCoord << "    ";
//    cout << "elecood = " << eleCoord << endl;
//    cout << "cornerIndex = " << cornerIndex << endl;
    boundaryCoord[cornerIndex][0] = eleCoordBoundaryCorner(0);
    boundaryCoord[cornerIndex][1] = eleCoordBoundaryCorner(1);
    boundaryCoord[cornerIndex][2] = 0.0;

    //boundaryCoord.push_back(trinodes);
    }

    if(globalHigherOrderIndex > -1)
    {
	    // store higher order node
        static BlitzVec3 eleCoordDomaninHO;
	    for(int k = 0; k < 3; k++)
	        eleCoordDomaninHO(k) = out.pointlist[globalHigherOrderIndex*3+k];



	    domainCoord[cornerIndex+3][0] = eleCoordDomaninHO(0);
	    domainCoord[cornerIndex+3][1] = eleCoordDomaninHO(1);
	    domainCoord[cornerIndex+3][2] = eleCoordDomaninHO(2);

	    static BlitzVec3 physCoordHO;
	    elementToCurrentCoordinates(xfemElement, xyze_xfemElement, eleCoordDomaninHO, physCoordHO);

	    //domainCoord.push_back(trinodes);

	    const DRT::Element* cutterElement = intersectingCutterElements_[faceMarker];
	    const BlitzMat xyze_cutterElement(getCurrentNodalPositions(cutterElement, currentcutterpositions));
	    BlitzVec2 eleCoordBoundaryHO;
	    CurrentToSurfaceElementCoordinates(cutterElement, xyze_cutterElement, physCoordHO, eleCoordBoundaryHO);

	//    cout << "physcood = " << physCoord << "    ";
	//    cout << "elecood = " << eleCoord << endl;

	    boundaryCoord[cornerIndex+3][0] = eleCoordBoundaryHO(0);
	    boundaryCoord[cornerIndex+3][1] = eleCoordBoundaryHO(1);
	    boundaryCoord[cornerIndex+3][2] = 0.0;

	    //boundaryCoord.push_back(trinodes);
    }
}



/*----------------------------------------------------------------------*
 |  DB:     Debug only                                       u.may 06/07|
 *----------------------------------------------------------------------*/
void XFEM::Intersection::debugXAABBIntersection(
        const BlitzMat3x2 cutterXAABB,
        const BlitzMat3x2 xfemXAABB,
        const DRT::Element* cutterElement,
        const DRT::Element* xfemElement,
        const int noC,
        const int noX
        ) const
{
	cout << endl;
    cout << "===============================================================" << endl;
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
    cout << endl;
    cout << "CUTTER XAABB: "<< "                     " << "XFEM XAABB: " << endl;
    cout << endl;
    cout <<  "minX = " << cutterXAABB(0,0) << "      " << "maxX = " << cutterXAABB(0,1) << "      " ;
    cout <<  "minX = " << xfemXAABB(0,0)   << "      " << "maxX = " << xfemXAABB(0,1)   << endl;
    cout <<  "minY = " << cutterXAABB(1,0) << "      " << "maxY = " << cutterXAABB(1,1) << "      " ;
    cout <<  "minY = " << xfemXAABB(1,0)   << "      " << "maxY = " << xfemXAABB(1,1)   << endl;
    cout <<  "minZ = " << cutterXAABB(2,0) << "      " << "maxZ = " << cutterXAABB(2,1) << "      " ;
    cout <<  "minZ = " << xfemXAABB(2,0)   << "      " << "maxZ = " << xfemXAABB(2,1) << endl;
    cout << endl;
    cout << endl;
	cout << "XFEM ELEMENT " << noX << " :" << endl;
	cout << endl;
	for(int jE = 0; jE < xfemElement->NumNode(); jE++)
	{
		xfemElement->Nodes()[jE]->Print(cout);
		cout << endl;
	}
    cout << endl;
    cout << endl;
    cout << "CUTTER XAABB: "<< "                     " << "XFEM XAABB: " << endl;
    cout << endl;
    cout <<  "minX = " << cutterXAABB(0,0) << "      " << "maxX = " << cutterXAABB(0,1) << "      " ;
    cout <<  "minX = " << xfemXAABB(0,0)   << "      " << "maxX = " << xfemXAABB(0,1)   << endl;
    cout <<  "minY = " << cutterXAABB(1,0) << "      " << "maxY = " << cutterXAABB(1,1) << "      " ;
    cout <<  "minY = " << xfemXAABB(1,0)   << "      " << "maxY = " << xfemXAABB(1,1)   << endl;
    cout <<  "minZ = " << cutterXAABB(2,0) << "      " << "maxZ = " << cutterXAABB(2,1) << "      " ;
    cout <<  "minZ = " << xfemXAABB(2,0)   << "      " << "maxZ = " << xfemXAABB(2,1) << endl;
	cout << endl;
	cout << endl;
    cout << "===============================================================" << endl;
    cout << "End Debug Intersection of XAABB's" << endl;
	cout << "===============================================================" << endl;
	cout << endl; cout << endl; cout << endl;


}



/*----------------------------------------------------------------------*
 |  DB:     Debug only                                       u.may 06/07|
 *----------------------------------------------------------------------*/
void XFEM::Intersection::debugNodeWithinElement(
        const DRT::Element* xfemElement,
        const DRT::Node* node,
        const BlitzVec& xsi,
        const int noE,
        const int noN,
        const bool within
        ) const
{
    const int numnodes = xfemElement->NumNode();
    vector<int> actParams(1,0);
    BlitzVec funct(numnodes);
    BlitzMat emptyM;
    BlitzVec emptyV;
    BlitzVec x(3);
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;

    params.set("action","calc_Shapefunction");
    actParams[0] = numnodes;

    dserror("we don't use Evaluate anymore, so thius function does not make sence!");
    //element->Evaluate(params, dummyDis, actParams, emptyM, emptyM , funct, xsi, emptyV);

    for(int dim=0; dim<3; dim++)
        x(dim) = 0.0;

    for(int dim=0; dim<3; dim++)
        for(int i=0; i<numnodes; i++)
        {
            x(dim) += xfemElement->Nodes()[i]->X()[dim] * funct(i);
        }

    cout << endl;
    cout << "===============================================================" << endl;
    cout << "Debug Node within element" << endl;
    cout << "===============================================================" << endl;
    cout << endl;
    cout << "ELEMENT " << noE << " :" << endl;
    cout << endl;
/*    for(int jE = 0; jE < element->NumNode(); jE++)
    {
        element->Nodes()[jE]->Print(cout);
        cout << endl;
    }
*/
    cout << endl;
    cout << endl;
    cout << "NODE " << noN << " :" << endl;
    cout << endl;
        node->Print(cout);
    cout << endl;
    cout << endl;
    cout << "XSI :";
    cout << "   r = " << xsi(0) << "     s = " <<  xsi(1) << "     t = "  << xsi(2) << endl;
    cout << endl;
    cout << endl;
    cout << "CURRENT COORDINATES :";
    cout << "   x = " << x[0] << "     y = " <<  x[1] << "     z = "  << x[2] << endl;
    cout << endl;
    cout << endl;
    if(within) cout << "NODE LIES WITHIN ELEMENT" << endl;
    else            cout << "NODE DOES NOT LIE WITHIN ELEMENT" << endl;
    cout << endl;
    cout << endl;
    cout << "===============================================================" << endl;
    cout << "End Debug Node within element" << endl;
    cout << "===============================================================" << endl;
    cout << endl; cout << endl; cout << endl;
}



/*----------------------------------------------------------------------*
 |  DB:     Debug only                                       u.may 06/07|
 *----------------------------------------------------------------------*/
void XFEM::Intersection::debugTetgenDataStructure(
        const DRT::Element*               xfemElement) const
{
    cout << endl;
    cout << "===============================================================" << endl;
    cout << "Debug Tetgen Data Structure " << endl;
    cout << "===============================================================" << endl;
    cout << endl;
    cout << "POINT LIST " << " :" << endl;
    cout << endl;
    BlitzVec xsi(3);
    for(unsigned int i = 0; i < pointList_.size(); i++)
    {
        for(int j = 0; j< 3; j++)
        {
            xsi[j] = pointList_[i].coord[j];
        }
        elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement_, xsi);

        cout << i << ".th point:   ";
        for(int j = 0; j< 3; j++)
        {
            //cout << setprecision(10) << pointList_[i].coord[j] << "\t";
             printf("%20.16f\t", pointList_[i].coord[j] );
        }
        cout << endl;
        cout << endl;

      /*  for(int j = 0; j< 3; j++)
        {
            cout << xsi[j] << "\t";
        }
        cout << endl;
        cout << endl;*/
    }
    cout << endl;
    cout << endl;

    cout << endl;
    cout << "SEGMENT LIST " << " :" << endl;
    cout << endl;
    for(unsigned int i = 0; i < segmentList_.size(); i++)
    {
        cout << i << ".th segment:   ";
        int count = 0;
        for(unsigned int j = 0; j < segmentList_[i].size(); j++)
                cout << segmentList_[i][count++] << "\t";

        count = 0;
        for(unsigned int j = 0; j < surfacePointList_[i].size(); j++)
                cout << surfacePointList_[i][count++] << "\t";

        cout << endl;
        cout << endl;
    }
    cout << endl;
    cout << endl;

    cout << endl;
    cout << "TRIANGLE LIST " << " :" << endl;
    cout << endl;
    for(unsigned int i = 0; i < triangleList_.size(); i++)
    {
        cout << i << ".th triangle:   ";
        for(int j = 0; j< 3; j++)
        {
            cout << triangleList_[i][j] << "\t";
        }
        cout << endl;
        cout << endl;
    }
    cout << endl;
    cout << endl;

    cout << "===============================================================" << endl;
    cout << "Debug Tetgen Data Structure" << endl;
    cout << "===============================================================" << endl;
    cout << endl; cout << endl; cout << endl;
}



/*----------------------------------------------------------------------*
 |  Debug only                                               u.may 06/07|
 *----------------------------------------------------------------------*/
void XFEM::Intersection::debugTetgenOutput( 	tetgenio& in,
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
void XFEM::Intersection::printTetViewOutput(
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
void XFEM::Intersection::printTetViewOutputPLC(
        const DRT::Element*   xfemElement,
        int                   index,
        tetgenio&             in) const
{

    FILE *outFile;
    char filename[100];
    BlitzVec xsi(3);

    sprintf(filename, "tetgenPLC%d.node", index);

    outFile = fopen(filename, "w");
    fprintf(outFile, "%d  %d  %d  %d\n", in.numberofpoints, in.mesh_dim,
          in.numberofpointattributes, in.pointmarkerlist != NULL ? 1 : 0);
    for (int i = 0; i < in.numberofpoints; i++)
    {

        for(int j = 0; j < 3; j++)
            xsi[j] = in.pointlist[i*3 + j];

        elementToCurrentCoordinatesInPlace(xfemElement, xyze_xfemElement_, xsi);

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


void XFEM::Intersection::debugFaceMarker(
		const int 						eleId,
		tetgenio&						out) const
{

  ofstream f_system("element_faceMarker.pos");
  f_system << "View \" Face Markers" << " \" {" << endl;

  for(int iface=0; iface<out.numberoftrifaces; iface++)
  {
    int trifaceMarker = out.trifacemarkerlist[iface] - facetMarkerOffset_;

    if(trifaceMarker > -2)
    {
      BlitzMat triface(3,3);
      for(int inode = 0; inode < 3; inode++)
        for(int isd = 0; isd < 3; isd++)
          triface(isd,inode) = out.pointlist[out.trifacelist[iface*3+inode]*3 + isd];

      f_system << IO::GMSH::cellWithScalarToString(DRT::Element::tri3, trifaceMarker, triface) << endl;
    }
  }
  f_system << "};" << endl;
  //f_system << IO::GMSH::getConfigString(2);
  f_system.close();
}


void XFEM::Intersection::debugXFEMConditions(
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
            f_system << IO::GMSH::elementToString(i, cutterElement) << endl;
        }
    }

	f_system << "};" << endl;
	f_system.close();
}


void XFEM::Intersection::debugIntersection(
		const DRT::Element*			xfemElement,
		vector<DRT::Element*>       cutterElements
		) const
{

	ofstream f_system("intersection.pos");
	f_system << "View \" Intersection" << " \" {" << endl;

	f_system << IO::GMSH::elementToString(0, xfemElement) << endl;

	for(unsigned int i=0; i<cutterElements.size(); i++)
		f_system << IO::GMSH::elementToString(i+1, cutterElements[i]) << endl;


	f_system << "};" << endl;
	f_system.close();
}

static string XAABBToString(const double scalar, const vector<vector<double> >& XAABB)
{
    const DRT::Element::DiscretizationType distype = DRT::Element::hex8;
    const int numnode = IO::GMSH::distypeToGmshNumNode(distype);

    stringstream pos_array_string;
    pos_array_string << "S" << IO::GMSH::distypeToGmshElementHeader(distype) << "(";
    for (int i = 0; i<numnode;++i)
    {
        pos_array_string << scientific << XAABB[i][0] << ",";
        pos_array_string << scientific << XAABB[i][1] << ",";
        pos_array_string << scientific << XAABB[i][2];
        if (i < numnode-1){
            pos_array_string << ",";
        }
    };
    pos_array_string << ")";
    // values
    pos_array_string << "{";
    for (int i = 0; i<numnode;++i)
    {
        pos_array_string << scientific << scalar;
        if (i < numnode-1){
            pos_array_string << ",";
        }
    };
    pos_array_string << "};";
    return pos_array_string.str();
}


void XFEM::Intersection::debugXAABBs(
	const int							id,
	const BlitzMat&     cutterXAABB,
	const BlitzMat&     xfemXAABB) const
{
	char filename[100];
	sprintf(filename, "element_XAABB%d.pos", id);

	ofstream f_system(filename);
	f_system << "View \" XAABB " << " \" {" << endl;
	std::vector<std::vector<double> > nodes(8, vector<double> (3, 0.0));

	//cutterXAABB
	nodes[0][0] = cutterXAABB(0,0);	nodes[0][1] = cutterXAABB(1,0);	nodes[0][2] = cutterXAABB(2,0);	// node 0
	nodes[1][0] = cutterXAABB(0,1);	nodes[1][1] = cutterXAABB(1,0);	nodes[1][2] = cutterXAABB(2,0);	// node 1
	nodes[2][0] = cutterXAABB(0,1);	nodes[2][1] = cutterXAABB(1,1);	nodes[2][2] = cutterXAABB(2,0);	// node 2
	nodes[3][0] = cutterXAABB(0,0);	nodes[3][1] = cutterXAABB(1,1);	nodes[3][2] = cutterXAABB(2,0);	// node 3
	nodes[4][0] = cutterXAABB(0,0);	nodes[4][1] = cutterXAABB(1,0);	nodes[4][2] = cutterXAABB(2,1);	// node 4
	nodes[5][0] = cutterXAABB(0,1);	nodes[5][1] = cutterXAABB(1,0);	nodes[5][2] = cutterXAABB(2,1);	// node 5
	nodes[6][0] = cutterXAABB(0,1);	nodes[6][1] = cutterXAABB(1,1);	nodes[6][2] = cutterXAABB(2,1);	// node 6
	nodes[7][0] = cutterXAABB(0,0);	nodes[7][1] = cutterXAABB(1,1);	nodes[7][2] = cutterXAABB(2,1);	// node 7

	f_system << XAABBToString(id+1, nodes) << endl;

	//xfemXAABB
	nodes[0][0] = xfemXAABB(0,0);	nodes[0][1] = xfemXAABB(1,0);	nodes[0][2] = xfemXAABB(2,0);	// node 0
	nodes[1][0] = xfemXAABB(0,1);	nodes[1][1] = xfemXAABB(1,0);	nodes[1][2] = xfemXAABB(2,0);	// node 1
	nodes[2][0] = xfemXAABB(0,1);	nodes[2][1] = xfemXAABB(1,1);	nodes[2][2] = xfemXAABB(2,0);	// node 2
	nodes[3][0] = xfemXAABB(0,0);	nodes[3][1] = xfemXAABB(1,1);	nodes[3][2] = xfemXAABB(2,0);	// node 3
	nodes[4][0] = xfemXAABB(0,0);	nodes[4][1] = xfemXAABB(1,0);	nodes[4][2] = xfemXAABB(2,1);	// node 4
	nodes[5][0] = xfemXAABB(0,1);	nodes[5][1] = xfemXAABB(1,0);	nodes[5][2] = xfemXAABB(2,1);	// node 5
	nodes[6][0] = xfemXAABB(0,1);	nodes[6][1] = xfemXAABB(1,1);	nodes[6][2] = xfemXAABB(2,1);	// node 6
	nodes[7][0] = xfemXAABB(0,0);	nodes[7][1] = xfemXAABB(1,1);	nodes[7][2] = xfemXAABB(2,1);	// node 7

	f_system << XAABBToString(0, nodes) << endl;

	f_system << "};" << endl;
	f_system.close();
}


#endif  // #ifdef CCADISCRET

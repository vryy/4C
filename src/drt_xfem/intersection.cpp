/*!----------------------------------------------------------------------
\file intersection.cpp

\brief collection of intersection tools for the interface determination of two meshes

    ML      math library
    GM      general methods
    ICS     intersection candidate search
    CLI     contruction of the linearized interface
    CDT     contrained delaunay triangulation
    RQI     recovery of the curved interface
    DB      debug methods
    
<pre>
Maintainer: Ursula Mayer
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "intersection.H"
#include "../drt_xfem/intersection_math.H"
#include "../drt_xfem/integrationcell.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_element.H"
//#include "../drt_f3/fluid3.H"

#ifdef PARALLEL
#include <mpi.h>
#endif 

using namespace XFEM;
using namespace DRT::Utils;

/*----------------------------------------------------------------------*
 |  MAIN:   computes the interface between the xfem          u.may 06/07|
 |          discretization and the cutter discretization.               |
 |          It returns a list of intersected xfem elements              |
 |          and their integrations cells.                               |
 *----------------------------------------------------------------------*/  
void Intersection::computeIntersection( const RefCountPtr<DRT::Discretization>  xfemdis,
                                        const RefCountPtr<DRT::Discretization>  cutterdis,  
                                        map< int, vector <Integrationcell> >&   integrationcellList)
{
    bool intersected =                      false;
    bool xfemIntersection =                 false;
    bool store = false;
    
    vector< DRT::Condition * >              xfemConditions;
    vector< InterfacePoint >                interfacePoints;
    
    map< int, RefCountPtr<DRT::Element > >  geometryMap;
    DRT::Element*                           cutterElement; 
    DRT::Element*                           xfemElement; 
    Epetra_SerialDenseMatrix                cutterXAABB;
    Epetra_SerialDenseMatrix                xfemXAABB;
   
    higherorder_ = false;
   
    
    // obtain vector of pointers to all xfem conditions of all cutter discretizations
    cutterdis->GetCondition ("XFEMCoupling", xfemConditions);
    //cout << endl << "number of xfem conditions =" << xfemConditions.size() << endl; cout << endl;
        
    if(xfemConditions.size()==0)
        dserror("number of fsi xfem conditions = 0");
        
    /* for(unsigned int i=0; i<xfemConditions.size(); i++)  cout << *xfemConditions[i]; */   

    /*cout << endl; 
    xfemElement = xfemdis->gElement(xfemdis->NumMyRowElements()-1);
    for(int jE = 0; jE < xfemElement->NumNode(); jE++)
    {
        xfemElement->Nodes()[jE]->Print(cout);
        cout << endl;
    }
    cout << endl;
    */
       
    //  k < xfemdis->NumMyRowElements()-1
    for(int k = 0; k < xfemdis->NumMyRowElements(); k++)
    {
        xfemIntersection = false;
        
        xfemElement = xfemdis->gElement(k);
        initializeXFEM(xfemElement);  
        xfemXAABB = computeFastXAABB(xfemElement);
            
        startPointList();
        
        //int countE = 0; 
        //xfemConditions.size()
        for(unsigned int i=0; i<xfemConditions.size(); i++)
        {
            geometryMap = xfemConditions[i]->Geometry();
            if(geometryMap.size()==0)   dserror("geometry does not obtain elements");
            printf("size of %d.geometry map = %d\n",i, geometryMap.size());
            
           
            for(unsigned int j=0; j<geometryMap.size(); j++)
            { 
                cutterElement = geometryMap.find(j)->second.get();
                if(cutterElement == NULL) dserror("geometry does not obtain elements");
            
                initializeCutter(cutterElement);  
                //printf("cutterElementcount = %d\n", countE++);
            
                cutterXAABB = computeFastXAABB(cutterElement);                                 
                intersected = intersectionOfXAABB(cutterXAABB, xfemXAABB);    
                //debugXAABBIntersection( cutterXAABB, xfemXAABB, cutterElement, xfemElement, i, k);
                                     
                if(intersected)
                {
                    
                    // collect internal points
                    numInternalPoints_= 0; 
                    numBoundaryPoints_ = 0;

                    for(int m=0; m<cutterElement->NumLine() ; m++)                    
                        collectInternalPoints( xfemElement, cutterElement, cutterElement->Nodes()[m],
                                                interfacePoints, k, m);
                      
                    // collect intersection points                                   
                    for(int m=0; m<xfemElement->NumLine() ; m++) 
                        if(collectIntersectionPoints(   cutterElement, xfemElement->Lines()[m],
                                                        interfacePoints, 0, m, false, xfemIntersection))                                         
                            storeIntersectedCutterElement(cutterElement); 
                    
                    
                    for(int m=0; m<cutterElement->NumLine() ; m++)                                              
                        for(int p=0; p<xfemElement->NumSurface() ; p++) 
                            if(collectIntersectionPoints(   xfemElement->Surfaces()[p], 
                                                            cutterElement->Lines()[m], interfacePoints,
                                                            p, m, true, xfemIntersection)) 
                                storeIntersectedCutterElement(cutterElement); 
                        
                        
                    // order interface points          
                    if(interfacePoints.size()!=0)
                    {
#ifdef QHULL
                        computeConvexHull( xfemElement, cutterElement, interfacePoints);
#else
                        dserror("Set QHULL flag to use XFEM intersections!!!");
#endif
                    }
                    
                    interfacePoints.clear();     
                }// if intersected
            }// for-loop over all geometryMap.size()
        }// for-loop over all xfemConditions.size() 
        
        if(xfemIntersection)
        {                                                                          
            //debugTetgenDataStructure(xfemElement);
            computeCDT(xfemElement, cutterElement, integrationcellList);
        }
 
  
    }// for-loop over all  actdis->NumMyRowElements()
    
    
    //debugIntegrationcells(integrationcellList,2);
    printf("\n");
    printf("Intersection computed sucessfully\n");
    printf("\n");
}



/*----------------------------------------------------------------------*
 |  INIT:   initializes the private members of the           u.may 08/07|
 |          current xfem element                                        |
 *----------------------------------------------------------------------*/
void Intersection::initializeXFEM(  
    DRT::Element* xfemElement)
{ 
    xfemDistype_ = xfemElement->Shape();
    
    numXFEMSurfaces_ = xfemElement->NumSurface(); 
    numXFEMLines_ = xfemElement->NumLine(); 
   
    numXFEMCornerNodes_  = getNumberOfElementCornerNodes(xfemDistype_);
    
    pointList_.clear();
    segmentList_.clear();                 
    surfacePointList_.clear();                
    triangleList_.clear();
    
    intersectedSurfaces_.clear();
    faceMarker_.clear();
    
    eleLinesSurfaces_.clear();
    eleNodesSurfaces_.clear();
    eleNumberingSurfaces_.clear();
    eleRefCoordinates_.clear(); 
     
    segmentList_.resize(numXFEMSurfaces_);
    surfacePointList_.resize(numXFEMSurfaces_);    
    
    eleLinesSurfaces_ = getEleNodeNumbering_lines_surfaces(xfemDistype_);
    eleNodesSurfaces_ = getEleNodeNumbering_nodes_surfaces(xfemDistype_);
    eleNumberingSurfaces_ = getEleNodeNumberingSurfaces(xfemDistype_);
    eleRefCoordinates_ = getEleNodeNumbering_nodes_reference(xfemDistype_);
}



/*----------------------------------------------------------------------*
 |  INIT:   initializes the private members of the           u.may 08/07|
 |          current cutter element                                        |
 *----------------------------------------------------------------------*/
void Intersection::initializeCutter(  
    DRT::Element* cutterElement)
{  
    cutterDistype_ = cutterElement->Shape();
}



/*----------------------------------------------------------------------*
 | GM:      transforms a node in reference coordinates       u.may 07/07|
 |          into current coordinates                                    |
 *----------------------------------------------------------------------*/  
void Intersection::referenceToCurrentCoordinates(   
    DRT::Element* element, 
    Epetra_SerialDenseVector& xsi)
{
    const int numNodes = element->NumNode();
    vector<int> actParams(1,0);
    Epetra_SerialDenseVector funct(numNodes);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;   
    params.set("action","calc_Shapefunction");
           
    actParams[0] = numNodes;            
    element->Evaluate(params, dummyDis, actParams, emptyM, emptyM , funct, xsi, emptyV);           
        
    for(int i=0; i<3; i++)
    {
        xsi[i] = 0.0;
        for(int j=0; j<numNodes; j++)
            xsi[i] += element->Nodes()[j]->X()[i] * funct(j);
    }
    
    /*
    for(int dim=0; dim<3; dim++)
    {
        xsi[dim] = 0.0;
        for(int j=0; j<numNodes; j++)
            xsi[dim] += element->Nodes()[j]->X()[dim] * funct(j);
    }
    */
}



/*----------------------------------------------------------------------*
 | GM:      transforms a node in reference coordinates       u.may 07/07|
 |          into current coordinates                                    |
 *----------------------------------------------------------------------*/  
void Intersection::referenceToCurrentCoordinates(   
    DRT::Element*                               surfaceElement, 
    Epetra_SerialDenseVector&                   xsi,
    const vector<Epetra_SerialDenseVector>&     surfaceNodes)
{
    const int numNodes = surfaceElement->NumNode();
    vector<int> actParams(1,0);
    Epetra_SerialDenseVector funct(numNodes);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;   
    params.set("action","calc_Shapefunction");
           
    actParams[0] = numNodes;            
    surfaceElement->Evaluate(params, dummyDis, actParams, emptyM, emptyM , funct, xsi, emptyV);           
        
    for(int dim=0; dim<3; dim++)
    {
        xsi[dim] = 0.0;
        for(int i=0; i<numNodes; i++)
            xsi[dim] += surfaceNodes[i][dim] * funct(i);
    }
}



/*----------------------------------------------------------------------*
 | GM:      transforms a node in reference coordinates       u.may 07/07|
 |          into current coordinates                                    |
 *----------------------------------------------------------------------*/  
void Intersection::referenceToCurrentCoordinates(  
    Epetra_SerialDenseVector&                   xsi,
    const vector<Epetra_SerialDenseVector>&     plane)
{
    const int numNodes = 4;
    Epetra_SerialDenseVector funct(numNodes);
    
    shape_function_2D(funct, xsi[0], xsi[1], DRT::Element::quad4);
        
    for(int dim=0; dim<3; dim++)
    {
        xsi[dim] = 0.0;
        for(int i=0; i<numNodes; i++)
            xsi[dim] += plane[i][dim] * funct(i);
    }
}



/*----------------------------------------------------------------------*
 | GM:  transforms a node in current coordinates            u.may 07/07 |
 |      into reference coordinates                                      |
 *----------------------------------------------------------------------*/  
void Intersection::currentToReferenceCoordinates(   
    DRT::Element*               element, 
    Epetra_SerialDenseVector&   xsi)
{
    bool nodeWithinElement;
    Epetra_SerialDenseVector x(3);
    
    for(int i = 0; i < 3; i++)
    {
        x[i] = xsi[i];
        xsi[i] = 0.0;
    }
    
    nodeWithinElement = checkNodeWithinElement(element, x, xsi);
    
/*    if(!nodeWithinElement)
    {
        printf("xsi = %f \t %f \t %f\n", xsi[0], xsi[1], xsi[2]);
    }
    //    dserror("node not within element");
*/        
        
    // rounding 1 and -1 to be exact for the CDT
    for(int j = 0; j < 3; j++)
    {
        if( fabs((fabs(xsi[j])-1.0)) < TOL7_ &&  xsi[j] < 0)    xsi[j] = -1.0;
        if( fabs((fabs(xsi[j])-1.0)) < TOL7_ &&  xsi[j] > 0)    xsi[j] =  1.0;      
    }  
}     



/*----------------------------------------------------------------------*
 |  GM:     compares two nodes                               u.may 06/07|
 |          overloaded method:                                          |
 |          double*  and  double*                                       |
 *----------------------------------------------------------------------*/  
bool Intersection::comparePoints(    
    const double*     point1,
    const double*     point2,
    const int         length)
{   
    bool equal = true;
             
    for(int i = 0; i < length; i++)
        if(fabs(point1[i] - point2[i]) > TOL7_)
        {
            equal = false;
            break;
        }
  
    return equal;
}



/*----------------------------------------------------------------------*
 |  GM:     compares two nodes                               u.may 06/07|
 |          overloaded method:  vector<double>  and double*             |
 *----------------------------------------------------------------------*/  
bool Intersection::comparePoints(   
    const vector<double>&     point1,
    const double*             point2)
{   
    bool equal = true;
    
    if(point2 == NULL)
        dserror("array is NULL");
        
    for(unsigned int i = 0; i < point1.size() ; i++)
        if(fabs(point1[i] - point2[i]) > TOL7_)
        {
            equal = false;
            break;
        }
    
    return equal;
}



/*----------------------------------------------------------------------*
 |  GM:     compares two nodes                               u.may 06/07|
 |          overloaded method:  vector<double>  and  vector<double>     |
 *----------------------------------------------------------------------*/  
bool Intersection::comparePoints(   
    const vector<double>& point1,
    const vector<double>& point2)
{   
    bool equal = true;
    
    if(point1.size() != point2.size())
        dserror("arrays of nodes need to have the same length");
             
    for(unsigned int i = 0; i < point1.size() ; i++)
        if(fabs(point1[i] - point2[i]) > TOL7_)
        {
            equal = false;
            break;
        }
  
    return equal;
}



/*----------------------------------------------------------------------*
 |  GM:     compares two nodes                               u.may 06/07|
 |          overloaded method:  DRT::Node*  and  DRT::Node*             |
 *----------------------------------------------------------------------*/  
bool Intersection::comparePoints( 
    const DRT::Node*     point1,
    const DRT::Node*     point2)
{
    bool equal = true;
    
             
    for(unsigned int i = 0; i < 3 ; i++)
        if(fabs(point1->X()[i] - point2->X()[i]) > TOL7_)
        {
            equal = false;
            break;
        }
  
    return equal;

}



/*----------------------------------------------------------------------*
 |  GM:     compares two nodes                               u.may 06/07|
 |          overloaded method:                                          |
 |          Epetra_SerialDenseVector  and  Epetra_SerialDenseVector     |
 *----------------------------------------------------------------------*/  
bool Intersection::comparePoints(    
    const Epetra_SerialDenseVector&     point1,
    const Epetra_SerialDenseVector&     point2)
{   
    bool equal = true;
    
    if(point1.Length() != point2.Length())
        dserror("arrays of nodes need to have the same length");
             
    for(int i = 0; i < point1.Length() ; i++)
        if(fabs(point1[i] - point2[i]) > TOL7_)
        {
            equal = false;
            break;
        }
  
    return equal;
}



/*----------------------------------------------------------------------*
 |  GM:     checks if a certain element is a                 u.may 06/07|
 |          volume  element with help of the discretization type        |
 *----------------------------------------------------------------------*/  
bool Intersection::checkIfVolumeElement(
    DRT::Element* element)
{
    bool isVolume = false;
    DRT::Element::DiscretizationType distype = element->Shape();
    
    if( distype == DRT::Element::hex8  ||
        distype == DRT::Element::hex20 ||
        distype == DRT::Element::hex27 ||
        distype == DRT::Element::tet4  ||
        distype == DRT::Element::tet10  )
    {
        isVolume = true;        
    }
    return isVolume;
}



/*----------------------------------------------------------------------*
 |  GM:     checks if a certain element is a                 u.may 06/07|
 |          with help of the discretization type                        |
 *----------------------------------------------------------------------*/  
bool Intersection::checkIfSurfaceElement(
    DRT::Element* element)
{
    bool isSurface = false;
    DRT::Element::DiscretizationType distype = element->Shape();
    
    if( distype == DRT::Element::quad4 ||
        distype == DRT::Element::quad8 ||
        distype == DRT::Element::quad9 ||
        distype == DRT::Element::tri3  ||
        distype == DRT::Element::tri6  )
    {
        isSurface = true;       
    }
    return isSurface;
}



/*----------------------------------------------------------------------*
 |  GM:     checks if a certain element is a                 u.may 06/07|
 |          line element with help of the discretization type           |
 *----------------------------------------------------------------------*/  
bool Intersection::checkIfLineElement(
    DRT::Element* element)
{
    bool isLine = false;
    DRT::Element::DiscretizationType distype = element->Shape();
    
    if( distype == DRT::Element::line2 ||
        distype == DRT::Element::line3 )
    {
        isLine = true;      
    }
    return isLine;
}



/*----------------------------------------------------------------------*
 |  GM:     returns the order of the element                 u.may 06/07|
 *----------------------------------------------------------------------*/  
int Intersection::getDimension(
    DRT::Element* element)
{
    int order = 0;
    DRT::Element::DiscretizationType distype = element->Shape();
    
    if( distype == DRT::Element::line2 ||
        distype == DRT::Element::line3 )
    {
        order = 1;      
    }
    if( distype == DRT::Element::quad4 ||
        distype == DRT::Element::quad8 ||
        distype == DRT::Element::quad9 ||
        distype == DRT::Element::tri3  ||
        distype == DRT::Element::tri6  )
    {
        order = 2;      
    }
    else if(    distype == DRT::Element::hex8  ||
        distype == DRT::Element::hex20 ||
        distype == DRT::Element::hex27 ||
        distype == DRT::Element::tet4  ||
        distype == DRT::Element::tet10  )
    {
        order = 3;      
    }
    
    return order;
}



/*----------------------------------------------------------------------*
 |  ICS:    computes an extended axis-aligned bounding box   u.may 06/07|
 |          XAABB for a given element                                   |
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix Intersection::computeFastXAABB( 
    DRT::Element* element)
{
    double  maxDistance; 
    Epetra_SerialDenseMatrix XAABB(3, 2);
    
	for(int dim=0; dim<3; dim++)
	{
		XAABB(dim, 0) = element->Nodes()[0]->X()[dim] - TOL7_;
   		XAABB(dim, 1) = element->Nodes()[0]->X()[dim] + TOL7_;
	}
    
    for(int i=1; i<element->NumNode(); i++)
        for(int dim=0; dim<3; dim++)
		{
            XAABB(dim, 0) = std::min( XAABB(dim, 0), element->Nodes()[i]->X()[dim] - TOL7_);
			XAABB(dim, 1) = std::max( XAABB(dim, 1), element->Nodes()[i]->X()[dim] + TOL7_);
		}
 
    maxDistance = fabs(XAABB(0,1) - XAABB(0,0));
 	for(int dim=1; dim<3; dim++)
	   maxDistance = std::max(maxDistance, fabs(XAABB(dim,1)-XAABB(dim,0)) );
	
    // subtracts half of the maximal distance to minX, minY, minZ
    // adds half of the maximal distance to maxX, maxY, maxZ 
	for(int dim=0; dim<3; dim++)
	{
		XAABB(dim, 0) = XAABB(dim, 0) - 0.5*maxDistance;
		XAABB(dim, 1) = XAABB(dim, 1) + 0.5*maxDistance;
	}	
	
	/*
    printf("\n");
	printf("axis-aligned bounding box:\n minX = %f\n minY = %f\n minZ = %f\n maxX = %f\n maxY = %f\n maxZ = %f\n", 
			  XAABB(0,0), XAABB(1,0), XAABB(2,0), XAABB(0,1), XAABB(1,1), XAABB(2,1));
	printf("\n");
	*/
    
	return XAABB;
}



/*----------------------------------------------------------------------*
 |  ICS:    checks if a node is within an XAABB               u.may 06/07|
 *----------------------------------------------------------------------*/
bool Intersection::isNodeWithinXAABB(    
    const std::vector<double>&         node,
    const Epetra_SerialDenseMatrix&    XAABB)
{
	bool isWithin = true;
    double diffMin = 0;
    double diffMax = 0;
    const double tol = TOL7_;
	
	for (int dim=0; dim<3; dim++)
	{
        diffMin = node[dim] - XAABB(dim,0);
        diffMax = XAABB(dim,1) - node[dim];
        
   	    if((diffMin < -tol)||(diffMax < -tol)) //check again !!!!!      
            isWithin = false;
    }
		 	
	return isWithin;
}



/*----------------------------------------------------------------------*
 |  ICS:    checks if two XAABB's intersect                  u.may 06/07|
 *----------------------------------------------------------------------*/
bool Intersection::intersectionOfXAABB(  
    const Epetra_SerialDenseMatrix&     cutterXAABB, 
    const Epetra_SerialDenseMatrix&     xfemXAABB)
{
	
  /*====================================================================*/
  /* bounding box topology*/
  /*--------------------------------------------------------------------*/
  /* parameter coordinates (x,y,z) of nodes
   * node 0: (minX, minY, minZ)
   * node 1: (maxX, minY, minZ)
   * node 2: (maxX, maxY, minZ)
   * node 3: (minX, maxY, minZ)
   * node 4: (minX, minY, maxZ)
   * node 5: (maxX, minY, maxZ)
   * node 6: (maxX, maxY, maxZ)
   * node 7: (minX, maxY, maxZ)
   * 
   *                      z
   *                      |           
   *             4========|================7
   *           //|        |               /||
   *          // |        |              //||
   *         //  |        |             // ||
   *        //   |        |            //  ||
   *       //    |        |           //   ||
   *      //     |        |          //    ||
   *     //      |        |         //     ||
   *     5=========================6       ||
   *    ||       |        |        ||      ||
   *    ||       |        o--------||---------y
   *    ||       |       /         ||      ||
   *    ||       0------/----------||------3
   *    ||      /      /           ||     //
   *    ||     /      /            ||    //
   *    ||    /      /             ||   //
   *    ||   /      /              ||  //
   *    ||  /      /               || //
   *    || /      x                ||//
   *    ||/                        ||/
   *     1=========================2
   *
   */
  /*====================================================================*/
	
    bool intersection =  false;
	std::vector < std::vector <double> > nodes(8, vector<double> (3, 0.0));
	
	nodes[0][0] = cutterXAABB(0,0);	nodes[0][1] = cutterXAABB(1,0);	nodes[0][2] = cutterXAABB(2,0);	// node 0	
	nodes[1][0] = cutterXAABB(0,1);	nodes[1][1] = cutterXAABB(1,0);	nodes[1][2] = cutterXAABB(2,0);	// node 1
	nodes[2][0] = cutterXAABB(0,1);	nodes[2][1] = cutterXAABB(1,1);	nodes[2][2] = cutterXAABB(2,0);	// node 2
	nodes[3][0] = cutterXAABB(0,0);	nodes[3][1] = cutterXAABB(1,1);	nodes[3][2] = cutterXAABB(2,0);	// node 3
	nodes[4][0] = cutterXAABB(0,0);	nodes[4][1] = cutterXAABB(1,0);	nodes[4][2] = cutterXAABB(2,1);	// node 4
	nodes[5][0] = cutterXAABB(0,1);	nodes[5][1] = cutterXAABB(1,0);	nodes[5][2] = cutterXAABB(2,1);	// node 5
	nodes[6][0] = cutterXAABB(0,1);	nodes[6][1] = cutterXAABB(1,1);	nodes[6][2] = cutterXAABB(2,1);	// node 6
	nodes[7][0] = cutterXAABB(0,0);	nodes[7][1] = cutterXAABB(1,1);	nodes[7][2] = cutterXAABB(2,1);	// node 7
	
	for (int i = 0; i < 8; i++)
		if(isNodeWithinXAABB(nodes[i], xfemXAABB))
		{
			intersection = true;
			break;
		}
	
		
	if(!intersection)
	{
		nodes[0][0] = xfemXAABB(0,0);	nodes[0][1] = xfemXAABB(1,0);	nodes[0][2] = xfemXAABB(2,0);	// node 0	
		nodes[1][0] = xfemXAABB(0,1);	nodes[1][1] = xfemXAABB(1,0);	nodes[1][2] = xfemXAABB(2,0);	// node 1
		nodes[2][0] = xfemXAABB(0,1);	nodes[2][1] = xfemXAABB(1,1);	nodes[2][2] = xfemXAABB(2,0);	// node 2
		nodes[3][0] = xfemXAABB(0,0);	nodes[3][1] = xfemXAABB(1,1);	nodes[3][2] = xfemXAABB(2,0);	// node 3
		nodes[4][0] = xfemXAABB(0,0);	nodes[4][1] = xfemXAABB(1,0);	nodes[4][2] = xfemXAABB(2,1);	// node 4
		nodes[5][0] = xfemXAABB(0,1);	nodes[5][1] = xfemXAABB(1,0);	nodes[5][2] = xfemXAABB(2,1);	// node 5
		nodes[6][0] = xfemXAABB(0,1);	nodes[6][1] = xfemXAABB(1,1);	nodes[6][2] = xfemXAABB(2,1);	// node 6
		nodes[7][0] = xfemXAABB(0,0);	nodes[7][1] = xfemXAABB(1,1);	nodes[7][2] = xfemXAABB(2,1);	// node 7
	
		for (int i = 0; i < 8; i++)
			if(isNodeWithinXAABB(nodes[i], cutterXAABB))
			{
				intersection = true;
				break;
			}
	}	
	return intersection;
}



/*----------------------------------------------------------------------*
 |  CLI:    collects points that belong to the interface     u.may 06/07|
 |          and lie within an xfem element                              |
 *----------------------------------------------------------------------*/
bool Intersection::collectInternalPoints(   
    DRT::Element*                   element,
    DRT::Element*                   surfaceElement,
    DRT::Node*                      node,
    std::vector< InterfacePoint >&  interfacePoints,
    const int                       elemId,
    const int                       nodeId)
{
    Epetra_SerialDenseVector xsi(3);
    Epetra_SerialDenseVector x(3);
       
    x[0] = node->X()[0];
    x[1] = node->X()[1];
    x[2] = node->X()[2];
    
    bool nodeWithinElement = checkNodeWithinElement(element, x, xsi);
   // debugNodeWithinElement(element,node,xsi,elemId ,nodeId, nodeWithinElement);  
    
    if(nodeWithinElement)
    {   
        InterfacePoint ip;
        //debugNodeWithinElement(element,node,xsi,elemId ,nodeId, nodeWithinElement);  
          
        numInternalPoints_++;
        // check if node lies on the boundary of the xfem element
        if(checkIfOnBoundary(xsi, ip))  numBoundaryPoints_++;
                                   
        
        // intersection coordinates in the surface 
        // element reference coordinate system
        getNodeCoordinates(nodeId, ip.coord, surfaceElement->Shape());
                          
        interfacePoints.push_back(ip);
        
        // for RQI do if decision
        storeIntersectedCutterElement(surfaceElement); 
    }
      
    return nodeWithinElement;
}



/*----------------------------------------------------------------------*
 |  CLI:    updates the Jacobi matrix for the computation    u.may 06/07|
 |          if a node is in a given element                             |
 *----------------------------------------------------------------------*/
void Intersection::updateAForNWE(   
    const int                   dim,
    Epetra_SerialDenseMatrix&   A,
    Epetra_SerialDenseVector&   xsi,
    DRT::Element*               element)                                                  
{	
    const int numNodes = element->NumNode();
    vector<int> actParams(1,0);
    Epetra_SerialDenseMatrix deriv1(dim, numNodes);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
   
    
    params.set("action","calc_ShapeDeriv1");
    actParams[0] = numNodes;  
       
    element->Evaluate(params, dummyDis, actParams, deriv1, emptyM , emptyV, xsi, emptyV);       
    
    
    for(int i=0; i<dim; i++)
        for(int j=0; j<dim; j++)
            A[i][j] = 0.0;
        
        
    for(int i=0; i<dim; i++)
        for(int k=0; k<dim; k++)
            for(int j=0; j<numNodes; j++)
                A[i][k] += element->Nodes()[j]->X()[i] * deriv1[j][k];
        
}



/*----------------------------------------------------------------------*
 |  CLI:    updates the rhs for the computation if a         u.may 06/07|
 |          node is in a given element                                  |
 *----------------------------------------------------------------------*/
void Intersection::updateRHSForNWE( 
    const int                           dim,
    Epetra_SerialDenseVector&           b,
    Epetra_SerialDenseVector&           xsi,
    const Epetra_SerialDenseVector&     x,
    DRT::Element*                       element)                                                  
{
    int numNodes = element->NumNode();
    vector<int> actParams(1,0);
    Epetra_SerialDenseVector funct(numNodes);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
      
    params.set("action","calc_Shapefunction");
    actParams[0] = numNodes;   
    
      
    element->Evaluate(params, dummyDis, actParams, emptyM, emptyM , funct, xsi, emptyV);        
    
    for(int i=0; i<dim; i++)
        b[i] = 0.0;
    
    for(int i=0; i<dim; i++)
        for(int j=0; j<numNodes; j++)
            b[i] += (-1.0) * element->Nodes()[j]->X()[i] * funct[j];
        
     for(int i=0; i<dim; i++)
        b[i] += x[i];
}



/*----------------------------------------------------------------------*
 |  CLI:    checks if a node is within a given element       u.may 06/07|	
 *----------------------------------------------------------------------*/
bool Intersection::checkNodeWithinElement(	
    DRT::Element*                       element,
    const Epetra_SerialDenseVector&     x,
    Epetra_SerialDenseVector&           xsi)
{

	bool nodeWithinElement = true;
    int iter = 0;
    int dim = getDimension(element);
    const int maxiter = 500;
    double residual = 1.0;
    
    Epetra_SerialDenseMatrix A(dim,dim);
    Epetra_SerialDenseVector b(dim);
    
    Epetra_SerialDenseMatrix A_old(dim,dim);
    Epetra_SerialDenseVector b_old(dim);
    Epetra_SerialDenseVector dx(dim);
    
    for(int i = 0; i<dim; i++)
        xsi[i] = 0.0;

    dx = xsi;
            
    updateRHSForNWE( dim, b, xsi, x, element);
   
    while(residual > TOL14_)
    { 
        updateAForNWE( dim, A, xsi, element);
        A_old = A;
        b_old = b;
        //!gaussElimination(A, b, dx, true, dim, 1)
        //!solveLinearSystemWithSVD(A, b, dx, dim)
        if(!gaussElimination(A, b, dx, true, dim, 1))
        {
            printf("MATRIX SINGULAR\n");
            nodeWithinElement = false;
            break;
        }   
        
        Epetra_SerialDenseVector b_new(dim);
        
        /*for(int i = 0 ; i < dim ; i++)
        {
            b_new[i] = 0.0;
            for(int j = 0 ; j < dim ; j++)
                b_new[i] += A_old[i][j]*dx[j];
        }
        printf("x = %20.16f\t, x = %20.16f\t, x = %20.16f\n", x[0],x[1],x[2]);
        printf("dx = %20.16f\t, dx = %20.16f\t, dx = %20.16f\n", dx[0],dx[1],dx[2]);
        printf("b = %20.16f\t, b = %20.16f\t, b = %20.16f\n", b_old[0],b_old[1],b_old[2]);
        printf("b1 = %20.16f\t, b1 = %20.16f\t, b1 = %20.16f\n", b_new[0],b_new[1],b_new[2]);
        printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\n", xsi[0],xsi[1],xsi[2]);
        
        printf("\n");
        */
        xsi = addTwoVectors(xsi,dx);
        
        if(iter >= maxiter)
        {   
            printf("ITERATION > maxiter\n");
            nodeWithinElement = false;
            break;
        }   
        
            
        updateRHSForNWE(dim, b, xsi, x, element);
        residual = b.Norm2();
        iter++; 
    }
    
    //printf("iter = %d\n", iter);
    //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi[0],xsi[1],xsi[2], residual, TOL14_);
    
    for(int i=0; i<dim; i++)
        if( (fabs(xsi[i])-1.0) > TOL7_)    // || (fabs(xsi[1])-1.0) > TOL7_  || (fabs(xsi[2])-1.0) > TOL7_ )   
        {    
            //printf(" i.%d  xsi = %f\n",i,  xsi[i]);
            nodeWithinElement = false;
            break;
        }
        
    return nodeWithinElement;
}



/*----------------------------------------------------------------------*
 |  CLI:    checks if a node that lies within an element     u.may 06/07|
 |          lies on one of its surfaces or nodes                        |                                           
 *----------------------------------------------------------------------*/
bool Intersection::checkIfOnBoundary( 
    Epetra_SerialDenseVector&       xsi,    
    InterfacePoint&                 ip)
{
    bool onSurface = false;
    int count = 0;
      
    count = getSurfaces(xsi, ip.surfaces, xfemDistype_);
        
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
bool Intersection::collectIntersectionPoints(  
    DRT::Element*                   surfaceElement,
    DRT::Element*                   lineElement,
    std::vector<InterfacePoint>&    interfacePoints, 
    const int                       surfaceId,
    const int                       lineId,
    const bool                      lines,
    bool&                           xfemIntersection)
{
    bool intersected = false;
    Epetra_SerialDenseVector xsi(3);
    Epetra_SerialDenseVector upLimit(3);
    Epetra_SerialDenseVector loLimit(3);
    int numInterfacePoints = 0;
  
    if(!checkIfSurfaceElement(surfaceElement))
        dserror("surface element has to be a surface element\n");
    if(!checkIfLineElement(lineElement))
        dserror("line element has to be a line element\n");
  
  
    for(int i = 0; i < 3; i++)
    {
       xsi[i] =  0.0; 
       upLimit[i]  =  1.0;  
       loLimit[i]  = -1.0;
       // extend for triangle
    }
    
    intersected = computeCurveSurfaceIntersection(surfaceElement, lineElement, xsi, upLimit, loLimit);
                                        
    if(intersected) 
        numInterfacePoints = addIntersectionPoint( surfaceElement, lineElement,xsi, upLimit, loLimit, 
                                                   interfacePoints, surfaceId, lineId, lines);
      
    
    // in this case a node of this line lies on the facet of the xfem element
    // but there is no intersection within the element                                          
    if(!((int) interfacePoints.size() == numBoundaryPoints_)) 
        xfemIntersection = true;
      
    return intersected;
}



/*----------------------------------------------------------------------*
 |  CLI:    computes the intersection between a              u.may 06/07|
 |          curve and a surface                    (CSI)                |
 *----------------------------------------------------------------------*/
bool Intersection::computeCurveSurfaceIntersection( 
    DRT::Element*               surfaceElement,
    DRT::Element*               lineElement,
    Epetra_SerialDenseVector&   xsi,
    Epetra_SerialDenseVector&   upLimit,
    Epetra_SerialDenseVector&   loLimit)
{
    int iter = 0;
    int maxiter = 500;
    bool intersection = true;
    double residual = 1.0;
    Epetra_SerialDenseMatrix A(3,3);
    Epetra_SerialDenseVector b(3);
    Epetra_SerialDenseVector dx(3);
   
    Epetra_SerialDenseMatrix A_old(3,3);
    Epetra_SerialDenseVector b_old(3);
    Epetra_SerialDenseVector b_gauss(3);
    dx = xsi;
 
    updateRHSForCSI( b, xsi, surfaceElement, lineElement);
  
                 
    while(residual > TOL14_)
    {   
        updateAForCSI( A, xsi, surfaceElement, lineElement);
   
        A_old = A;
        b_old = b;
        
      /*if(!gaussElimination(A, b, dx, true, 3, 1))
        {
            intersection = false;
            break;
        } 
        */
        if(!solveLinearSystemWithSVD(A, b, dx, 3))
        {
            intersection = false;
            break;
        } 
   
      /*printf("\n");  
        printf("Intersection\n");  
        printf("=============================================================================\n");  
        printf("\n");  
        printf("dx = %20.16f\t, dx = %20.16f\t, dx = %20.16f\n", dx[0], dx[1], dx[2]);
        
       for(int i = 0; i < 3; i++ )
       {
            b_gauss[i] = 0.0;
            for(int  k = 0; k < 3; k++ )
            {  
                b_gauss[i] += A_old[i][k]*dx[k];
            }
            printf("b = %f\t  b_gauss = %f\n", b_old[i], b_gauss[i]);
        }
        printf("\n");
        printf("\n");  
        printf("=============================================================================\n");  
        printf("\n");  
        */
        
        xsi = addTwoVectors(xsi,dx);
       
        if(iter >= maxiter)
        {   
            intersection = false;
            break;
        }       
       
        updateRHSForCSI( b, xsi, surfaceElement, lineElement);
        residual = b.Norm2(); 
        iter++;
    } 
    
    
    if( (xsi[0] > upLimit[0]+TOL7_) || (xsi[1] > upLimit[1]+TOL7_) || (xsi[2] > upLimit[2]+TOL7_)  || 
        (xsi[0] < loLimit[0]-TOL7_) || (xsi[1] < loLimit[1]-TOL7_) || (xsi[2] < loLimit[2]-TOL7_)) 
            intersection = false;
        
        
    return intersection;
}



/*----------------------------------------------------------------------*
 |  CLI:    updates the systemmatrix for the computation     u.may 06/07| 
 |          of a curve-surface intersection (CSI)                       |
 *----------------------------------------------------------------------*/
void Intersection::updateAForCSI(  	Epetra_SerialDenseMatrix&   A,
                                    Epetra_SerialDenseVector&   xsi,
                                    DRT::Element*               surfaceElement,
                                    DRT::Element*               lineElement)        											
{	
	const int numNodesSurface = surfaceElement->NumNode();
   	const int numNodesLine = lineElement->NumNode();
	vector<int> actParams(1,0);
	Epetra_SerialDenseMatrix surfaceDeriv1(2,numNodesSurface);
	Epetra_SerialDenseMatrix lineDeriv1(1,numNodesLine);
	Epetra_SerialDenseMatrix emptyM;
	Epetra_SerialDenseVector emptyV;
    Epetra_SerialDenseVector xsiSurface(2);
    Epetra_SerialDenseVector xsiLine(1);
	DRT::Discretization dummyDis("dummy discretization", null);
	ParameterList params;
	 
    xsiSurface(0) = xsi[0]; // r-coordinate surface
    xsiSurface(1) = xsi[1]; // s-coordinate surface
    xsiLine(0) = xsi[2];    // r-coordinate line
	actParams[0] = numNodesSurface;     
	
    for(unsigned int i=0; i<3; i++)
        for(unsigned int j=0; j<3; j++)
            A[i][j] = 0.0;
	
	params.set("action","calc_ShapeDeriv1");
	surfaceElement->Evaluate(params, dummyDis, actParams, surfaceDeriv1, emptyM , emptyV, xsiSurface, emptyV);			
    //printf("num of nodes surfaces = %d\n", numOfNodesSurface);
       	
    for(int dim=0; dim<3; dim++)
   	    for(int i=0; i<numNodesSurface; i++)
		{
			A[dim][0] += surfaceElement->Nodes()[i]->X()[dim] * surfaceDeriv1(0,i);
			A[dim][1] += surfaceElement->Nodes()[i]->X()[dim] * surfaceDeriv1(1,i);
		}
		
   
    actParams[0] = numNodesLine;
    lineElement->Evaluate(params, dummyDis, actParams, lineDeriv1, emptyM , emptyV, xsiLine, emptyV);

    for(int dim=0; dim<3; dim++)
        for(int i=0; i<numNodesLine; i++)
        {
            A[dim][2] +=  (-1) * lineElement->Nodes()[i]->X()[dim] * lineDeriv1(0,i);
        }
        
}



/*----------------------------------------------------------------------*
 |  CLI:    updates the right-hand-side for the              u.may 06/07|
 |          computation of                                              |
 |          a curve-surface intersection (CSI)                          |
 *----------------------------------------------------------------------*/
void Intersection::updateRHSForCSI( 
    Epetra_SerialDenseVector&   b,
    Epetra_SerialDenseVector&   xsi,
    DRT::Element*               surfaceElement,
    DRT::Element*               lineElement)        											
{
    int numNodesSurface = surfaceElement->NumNode();
    int numNodesLine = lineElement->NumNode();
    vector<int> actParams(1,0);
    Epetra_SerialDenseVector surfaceFunct(numNodesSurface);
    Epetra_SerialDenseVector lineFunct(numNodesLine);
    Epetra_SerialDenseMatrix emptyML(1,numNodesLine);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    Epetra_SerialDenseVector xsiSurface(2);
    Epetra_SerialDenseVector xsiLine(1);
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
    
    
    xsiSurface(0) = xsi[0]; // r-coordinate surface
    xsiSurface(1) = xsi[1]; // s-coordinate surface
    xsiLine(0) = xsi[2];    // r-coordinate line
    params.set("action","calc_Shapefunction");
    		
   	for(unsigned int i=0; i<3; i++)   
		b[i] = 0.0;
    
   	actParams[0] = numNodesSurface;            
	surfaceElement->Evaluate(params, dummyDis, actParams, emptyM, emptyM , surfaceFunct, xsiSurface, emptyV);			
  
	actParams[0] = numNodesLine;     
	lineElement->Evaluate(params, dummyDis, actParams, emptyM, emptyML, lineFunct, xsiLine, emptyV);
		
   	for(int dim=0; dim<3; dim++)
   		for(int i=0; i<numNodesSurface; i++)
		{
			b[dim] += (-1) * surfaceElement->Nodes()[i]->X()[dim] * surfaceFunct(i);
		}
	
	for(int dim=0; dim<3; dim++)
   		for(int i=0; i<numNodesLine; i++)
		{
			b[dim] += lineElement->Nodes()[i]->X()[dim] * lineFunct(i);
		}
}



/*----------------------------------------------------------------------*
 |  CLI:    computes a new starting point for the            u.may 06/07| 
 |          Newton-method in order to find all intersection points      | 
 |          of a curve-surface intersection                             |
 *----------------------------------------------------------------------*/  
int Intersection::computeNewStartingPoint(
    DRT::Element*                surfaceElement,
    DRT::Element*                lineElement,
    const int                    surfaceId,
    const int                    lineId,
    Epetra_SerialDenseVector&    upLimit,
    Epetra_SerialDenseVector&    loLimit,
    std::vector<InterfacePoint>& interfacePoints,
    const bool                   lines)
{	
	bool intersected = false;
    bool interval = true;
	int numInterfacePoints = 0;
	Epetra_SerialDenseVector xsi(3);
	
    //printf("xsi2 = %f   %f   %f\n", fabs(xsi[0]), fabs(xsi[1]), fabs(xsi[2]) ); 
    //printf("lolimit = %f   %f   %f\n", fabs(loLimit[0]), fabs(loLimit[1]), fabs(loLimit[2]) ); 
    //printf("uplimit = %f   %f   %f\n", fabs(upLimit[0]), fabs(upLimit[1]), fabs(upLimit[2]) ); 
    
    if(comparePoints(upLimit, loLimit))
        interval = false;
        
	for(int i = 0; i < 3; i++)
		xsi[i] = (double) (( upLimit[i] + loLimit[i] )/2.0);
       
 
	intersected = computeCurveSurfaceIntersection(surfaceElement, lineElement, xsi, upLimit, loLimit);
   
    
    if( comparePoints(xsi,upLimit))     intersected = false;
    if( comparePoints(xsi,loLimit))     intersected = false;
       							 	
	if(intersected && interval)		
   		numInterfacePoints = addIntersectionPoint(	surfaceElement, lineElement,xsi, upLimit, loLimit, 
   													interfacePoints, surfaceId, lineId, lines);
   	
   	//printf("number of intersection points = %d\n", numInterfacePoints );
   	return numInterfacePoints;    						
}



/*----------------------------------------------------------------------*
 |  CLI:    adds an intersection point to the 			     u.may 07/07|
 |          list of interface points                       			    |
 *----------------------------------------------------------------------*/  
int Intersection::addIntersectionPoint(
    DRT::Element*                	surfaceElement,
    DRT::Element*                	lineElement,
    Epetra_SerialDenseVector&		xsi,
    Epetra_SerialDenseVector&     	upLimit,
    Epetra_SerialDenseVector&     	loLimit,
    std::vector<InterfacePoint>& 	interfacePoints,
    const int                       surfaceId,
    const int                       lineId,							  
    const bool 					    lines)
{

	int numInterfacePoints = 0;
	
 	InterfacePoint ip;
    if(lines)
    {   
        ip.nsurf = 1;
        ip.surfaces[0] = surfaceId;
         
        getLineCoordinates(lineId, xsi[2], ip.coord, surfaceElement->Shape());
    }
    else
    {
        ip.nsurf = 2;
        ip.surfaces[0] = eleLinesSurfaces_[lineId][0];
        ip.surfaces[1] = eleLinesSurfaces_[lineId][1];
        ip.coord[0] = xsi[0]; 
        ip.coord[1] = xsi[1]; 
    }
    
    ip.coord[2] = 0.0; 
    ip.pType = INTERSECTION;
    
    vector<InterfacePoint>::iterator it;
    bool alreadyInList = false;
    for(it = interfacePoints.begin(); it != interfacePoints.end(); it++ )  
        if(comparePoints(ip.coord, it->coord,3))   
        {   
            //printf("alreadyinlist = true\n");
            alreadyInList = true;
            break;
            
        }
      
    if(!alreadyInList)      
    {
        //printf("alreadyinlist = false\n");
        interfacePoints.push_back(ip);  
        numInterfacePoints++;
    }   
    
    // recursive call 8 times !!!!
    numInterfacePoints =  numInterfacePoints +
      computeNewStartingPoint(	surfaceElement, lineElement, surfaceId, lineId, 
      							upLimit, xsi, interfacePoints, lines) +
      computeNewStartingPoint(	surfaceElement, lineElement, surfaceId, lineId, 
      							xsi, loLimit, interfacePoints, lines);

	return numInterfacePoints;
}



/*----------------------------------------------------------------------*
 |  ICS:    computes the convex hull of a set of             u.may 06/07|
 |          interface points and stores resulting points,               |
 |          segments and triangles for the use with Tetgen (CDT)        |
 *----------------------------------------------------------------------*/
#ifdef QHULL
void Intersection::computeConvexHull(   
    DRT::Element*           xfemElement,
    DRT::Element*           surfaceElement,
    vector<InterfacePoint>& interfacePoints)
{
    
    double*                     point;
    vector<int>                 positions;
    vector<double>              searchPoint(3,0);
    vector<double>              vertex(3,0);
    vector< vector<double> >    vertices;  
    Epetra_SerialDenseVector    curCoord(3);  
    InterfacePoint              midpoint;
    
           
    if(!checkIfSurfaceElement(surfaceElement))
   		dserror("surface element has to be a surface element\n");
       
           
    if(interfacePoints.size() != 0)
    {    
        if(interfacePoints.size() > 2)  
        {
           
            midpoint = computeMidpoint(interfacePoints);
            // transform it into current coordinates
            for(int j = 0; j < 2; j++)      curCoord[j]  = midpoint.coord[j];            
            referenceToCurrentCoordinates(surfaceElement, curCoord);    
            currentToReferenceCoordinates(xfemElement, curCoord);    
            for(int j = 0; j < 3; j++)      midpoint.coord[j] = curCoord[j]; 
         
            // store coordinates in 
            // points has numInterfacePoints*dim-dimensional components
            // points[0] is the first coordinate of the first point
            // points[1] is the second coordinate of the first point
            // points[dim] is the first coordinate of the second point                     
            coordT* coordinates = (coordT *)malloc((2*interfacePoints.size())*sizeof(coordT));
            int fill = 0;
            for(unsigned int i = 0; i < interfacePoints.size(); i++)
            {
                for(int j = 0; j < 2; j++)
                {
                    coordinates[fill++] = interfacePoints[i].coord[j]; 
                   // printf("coord = %f\t", interfacePoints[i].coord[j]);
                }
                // printf("\n");
                // transform interface points into current coordinates
                for(int j = 0; j < 2; j++)      
                    curCoord[j]  = interfacePoints[i].coord[j];   
                              
                referenceToCurrentCoordinates(surfaceElement, curCoord);  
                currentToReferenceCoordinates(xfemElement, curCoord);
                         
                for(int j = 0; j < 3; j++)         
                    interfacePoints[i].coord[j] = curCoord[j];
                    
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
                    point  = SETelemt_(facet->vertices, j, vertexT)->point;
                    for(int k = 0; k < 2; k++)  
                        vertex[k] = point[k];
                                    
                    for(int m = 0; m < 2; m++)      curCoord[m]  = vertex[m];           
                    // surface reference coordinates to current coordinates       
                    referenceToCurrentCoordinates(surfaceElement, curCoord); 
                    // current coordinates to xfem element reference coordinates
                    currentToReferenceCoordinates(xfemElement, curCoord);    
                    for(int m = 0; m < 3; m++)      vertex[m] = curCoord[m];
                                                    
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
        else
        {       
            for(unsigned int i = 0; i < interfacePoints.size(); i++)
            {
                // transform interface points into current coordinates
                for(int j = 0; j < 2; j++)         
                    curCoord[j]  = interfacePoints[i].coord[j];   
                
                // surface reference coordinates to current coordinates       
                referenceToCurrentCoordinates(surfaceElement, curCoord); 
                // current coordinates to xfem element reference coordinates
                currentToReferenceCoordinates(xfemElement, curCoord);   
                  
                for(int j = 0; j < 3; j++)      
                {   
                    interfacePoints[i].coord[j] = curCoord[j];
                    vertex[j] = curCoord[j];
                }
                vertices.push_back(vertex);
            }              
        }  
    
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
        
        int countWhile = 0;
        while(vertices.size()>2)
        {                    
            findNextSegment(vertices, searchPoint);
            storePoint(searchPoint, interfacePoints, positions);
            countWhile++;           
        } 
        vertices.clear();
       
       
        storeSurfacePoints(interfacePoints);
            
        // cutter element lies on the surface of an xfem element
        if(numInternalPoints_ == numBoundaryPoints_ && numInternalPoints_ != 0)
        {
            if(numBoundaryPoints_ > 1)              
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
}
#endif //QHULL
 


/*----------------------------------------------------------------------*
 |  ICS:    finds the next facet of a convex hull            u.may 06/07|
 |          and returns the point different form the searchpoint        |
 *----------------------------------------------------------------------*/  
void Intersection::findNextSegment(   
    vector< vector<double> >&   vertices, 
    vector<double>&             searchPoint)
{     
    vector< vector<double> >::iterator it;
    bool pointfound = false;
    
    if(vertices.size()==0 || searchPoint.size()==0)
        dserror("one or both vectors are empty");   
    
    for(it = vertices.begin(); it != vertices.end(); it=it+2 )
    {      
        if(comparePoints(searchPoint, *it))
        {
            pointfound = true;
            searchPoint = *(it+1);              
            vertices.erase(it);
            vertices.erase(it); // remove it+ 1
            break; 
        }
               
        if(comparePoints(searchPoint, *(it+1)))
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
 |  for an intersected xfem element in reference configuration          |
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
void Intersection::computeCDT(  
    DRT::Element*               			element,
    DRT::Element*               			cutterElement,
    map< int, vector <Integrationcell> >&	integrationcellList)
{
    int dim = 3;
    int nsegments = 0; 
    int nsurfPoints = 0;
    tetgenio in;
    tetgenio out;
    char switches[] = "pnno2QY";    //o2
    tetgenio::facet *f;
    tetgenio::polygon *p;
    double regionCoordinates[6];
    

    // allocate pointlist
    in.numberofpoints = pointList_.size();
    in.pointlist = new REAL[in.numberofpoints * dim];
       
    // fill point list
    int fill = 0;
    for(int i = 0; i <  in.numberofpoints; i++)
        for(int j = 0; j < dim; j++)  
            in.pointlist[fill++] = (REAL) pointList_[i].coord[j]; 
 
 
    in.pointmarkerlist = new int[in.numberofpoints];   
    for(int i = 0; i < numXFEMCornerNodes_; i++)
        in.pointmarkerlist[i] = 3;    // 3 : point not lying on the xfem element
        
    for(int i = numXFEMCornerNodes_; i < in.numberofpoints; i++)
        in.pointmarkerlist[i] = 2;    // 2 : point not lying on the xfem element

   
    if(triangleList_.size()>0)      in.numberoffacets = numXFEMSurfaces_ + triangleList_.size(); 
    else                            in.numberoffacets = numXFEMSurfaces_;   
      
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    
    
    // loop over all xfem element surfaces
    for(int i = 0; i < numXFEMSurfaces_; i++)
    {
        f = &in.facetlist[i];
        if(segmentList_[i].size() > 0)          nsegments = (int) (segmentList_[i].size()/2);
        else                                    nsegments = 0;
        if(surfacePointList_[i].size() > 0)     nsurfPoints = surfacePointList_[i].size();
        else                                    nsurfPoints = 0;
        f->numberofpolygons = 1 + nsegments + nsurfPoints; 
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];
        p->numberofvertices = 4;
        p->vertexlist = new int[p->numberofvertices];
        for(int j = 0; j < 4; j ++)
            p->vertexlist[j] = eleNumberingSurfaces_[i][j];
           
      
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
    
    // store triangles
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
        for(int j = 0; j < 3; j ++)
            p->vertexlist[j] = triangleList_[i - element->NumSurface()][j];    
    }
      
        
    // set facetmarkers
    for(int i = 0; i < in.numberoffacets; i ++)
        in.facetmarkerlist[i] = faceMarker_[i] + facetMarkerOffset_;   
        
        
	// specify regions
	bool regions = false;	
	if(regions)
	{
		//computeRegionCoordinates(element,cutterElement,regionCoordinates);
	    fill = 0;
	    int read = 0;
	    in.numberofregions = 2;
	    in.regionlist = new REAL[in.numberofregions*5];
	  	//for(int i = 0; i < in.numberofregions; i++)
	    //{
	    	// store coordinates 
	    	//for(int j = 0; j < 3; j++)
	    	//	in.regionlist[fill++] = regionCoordinates[read++];
            in.regionlist[fill++] = 0.9;
            in.regionlist[fill++] = 0.9;
            in.regionlist[fill++] = 0.9;
	    		
	    	// store regional attribute (switch A)   i=0 cutter i=1 fluid
	    	in.regionlist[fill++] = 0;
	    	
	    	// store volume constraint (switch a)   
	    	in.regionlist[fill++] = 0.0;
            
            in.regionlist[fill++] = -1.0;
            in.regionlist[fill++] = -1.0;
            in.regionlist[fill++] = -1.0;
                
            // store regional attribute (switch A)   i=0 cutter i=1 fluid
            in.regionlist[fill++] = 1;
            
            // store volume constraint (switch a)   
            in.regionlist[fill++] = 0.0;
		//}
	}
	   
    in.save_nodes("tetin");   
    in.save_poly("tetin");   
    //  Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
    //  do quality mesh generation (q) with a specified quality bound
    //  (1.414), and apply a maximum volume constraint (a0.1)
    tetrahedralize(switches, &in, &out); 
 
    
    //Debug
    vector<int> elementIds;
    for(int i = 0; i<numXFEMCornerNodes_; i++)
        elementIds.push_back(i);
    
    debugTetgenOutput(in, out, element, elementIds);
    //printTetViewOutputPLC( element, element->Id(), in);
    // recover curved interface for higher order meshes
    bool curvedInterface = true;
    if(curvedInterface)
        recoverCurvedInterface(element, out);
 
    //printTetViewOutput(element->Id(), out);
   
    
   
    // store integrationcells
    vector<double> tetnodes(3);
    vector< vector<double> > tetrahedronCoord;
    vector< Integrationcell > listperElement;
    
    for(int i=0; i<out.numberoftetrahedra; i++ )
    {   
        for(int j = 0; j < out.numberofcorners; j++)
        {
            for(int dim = 0; dim < 3; dim++)
                tetnodes[dim] = out.pointlist[out.tetrahedronlist[i*out.numberofcorners+j]*3+dim];
         
            tetrahedronCoord.push_back(tetnodes);    
        }
        Integrationcell cell(i, tetrahedronCoord);
        tetrahedronCoord.clear();
        listperElement.push_back(cell);                 
    }
    integrationcellList.insert(make_pair(element->Id(),listperElement));
    
    // clear vectors
    listperElement.clear();
}



/*----------------------------------------------------------------------*
 |  CDT:    fills the point list with the corner points      u.may 06/07|
 |          in reference coordinates of the xfem element                |
 *----------------------------------------------------------------------*/  
void Intersection::startPointList(
    )
{
    InterfacePoint ip;
        
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
void Intersection::storePoint(  
    const vector<double>&       point, 
    vector<InterfacePoint>&     interfacePoints, 
    vector<int>&                positions)
{
    int count = 0;
    bool alreadyInList = false; 
    vector<InterfacePoint>::iterator it1;
        

    for(unsigned int i = 0; i < interfacePoints.size(); i++ )
    {   
        if(comparePoints(point, interfacePoints[i].coord))
        {         
            alreadyInList = false;
            count = 0;
            for(it1 = pointList_.begin(); it1 != pointList_.end(); it1++ )  
            {
                if(comparePoints(point, it1->coord))   
                {   
                    alreadyInList = true;
                    break;
                }
                count++;
            }  
            
            if(!alreadyInList) 
            {               
                pointList_.push_back(interfacePoints[i]);
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
InterfacePoint Intersection::computeMidpoint(
    const vector<InterfacePoint>& interfacePoints)
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
void Intersection::storeSurfacePoints(  
    vector<InterfacePoint>&     interfacePoints)
{   
    
    int surf1;
    int surf2;
    bool singlePoint = true;
    bool alreadyInList = false;
    vector<InterfacePoint>::iterator it1;
   
   
    for(unsigned int i = 0; i < interfacePoints.size(); i++)
    {
        singlePoint = true;
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
                            surf1 = interfacePoints[i].surfaces[k];
                            surf2 = interfacePoints[j].surfaces[l];
                        
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
                if(comparePoints(interfacePoints[i].coord, pointList_[jj].coord, 3)) 
                {
                    alreadyInList = false;
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
void Intersection::storeSegments(   
    const vector<int>&              positions)
{
    int pos1 = 0;
    int pos2 = 0;
    int surf1 = 0;
    int surf2 = 0;
    bool alreadyInList = false;
    
    
    for(unsigned int i = 0; i < positions.size()-1; i++ )
    {
        pos1 = positions[i];
        pos2 = positions[i+1];
        
        for(int j = 0; j < pointList_[pos1].nsurf; j++ )  
            for(int k = 0; k < pointList_[pos2].nsurf; k++ ) 
            {
                surf1 = pointList_[pos1].surfaces[j];
                surf2 = pointList_[pos2].surfaces[k];
             
                if( (surf1 == surf2) ) 
                { 
                    alreadyInList = false;
                    
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
      
    pos1 = positions[positions.size()-1];
    pos2 = positions[0]; 
    
    for(int j = 0; j < pointList_[pos1].nsurf; j++ )  
        for(int k = 0; k < pointList_[pos2].nsurf; k++ ) 
        {
            surf1 = pointList_[pos1].surfaces[j];
            surf2 = pointList_[pos2].surfaces[k];

            if((surf1 == surf2) ) 
            { 
                alreadyInList = false;
                for(unsigned int is = 0 ; is < segmentList_[surf1].size(); is = is + 2)
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
    


/*----------------------------------------------------------------------*
 |  CDT:    stores a triangle within a list of trianles      u.may 06/07|
 |          which is to be copy to the tetgen data structure            |
 |          for the computation of the                                  |
 |          Constrained Delauney Triangulation                          |
 *----------------------------------------------------------------------*/  
void Intersection::storeTriangles(  
    const vector<int>               positions)
{
    vector<int> triangle(3,0);
    
    for(unsigned int i = 0; i < positions.size()-1; i++ )
    {
        triangle[0] = positions[i];
        triangle[1] = positions[i+1];
        triangle[2] = pointList_.size()-1;
        
        triangleList_.push_back(triangle);
        faceMarker_.push_back(intersectedSurfaces_.size()-1);
    }
    
    triangle[0] = positions[positions.size()-1];
    triangle[1] = positions[0];
    triangle[2] = pointList_.size()-1;
        
    triangleList_.push_back(triangle);
    faceMarker_.push_back(intersectedSurfaces_.size()-1);
}



/*----------------------------------------------------------------------*
 |  RQI:    recovers the curved interface after the          u.may 08/07|
 |          Constrained Delaunay Tetrahedralization                     |
 *----------------------------------------------------------------------*/  
void Intersection::recoverCurvedInterface(
    DRT::Element*   element, 
    tetgenio&       out
    )
{
    bool currentDomain;
    bool intersected = false;
    int index1, index2;
    int tetIndex, faceMarkerIndex, localHigherOrderIndex, globalHigherOrderIndex;
    vector<int> order(3,0);
    vector<int> tetraCornerIndices(4,0); 
    Epetra_SerialDenseVector node2(3);
    Epetra_SerialDenseVector node1(3);
    Epetra_SerialDenseVector xsi(3);  
    vector < Epetra_SerialDenseVector > plane(5, Epetra_SerialDenseVector(3)); 
    vector < vector <Epetra_SerialDenseVector> > surfaceNodes;
    vector < Epetra_SerialDenseVector >  tetraCornerNodes(4, Epetra_SerialDenseVector(3) );
    
    // list of point markers , if already visited = 1 , if not = 0
    int* visitedPointIndexList = new int[out.numberofpoints];      
    for(int i = 0; i<out.numberofpoints; i++)
        visitedPointIndexList[i] = 0;
        
        
    // transform surface elements into the ref system of the xfem element
    transformSurfaceNodes(element, surfaceNodes);    
        
    // lifts all corner points into the curved interface
    liftAllCornerPoints(out, element, surfaceNodes);
     
    int countTri = 0;
   
    for(int i=0; i<out.numberoftrifaces; i++)
    {
        printf("tri face marker = %d\n", out.trifacemarkerlist[i]);   
        printf("num intersected surfaces = %d\n", intersectedSurfaces_.size()); 
        // run over all faces not lying in on of the xfem element planes
        faceMarkerIndex = out.trifacemarkerlist[i] - facetMarkerOffset_;
        
        if(faceMarkerIndex > -1)
        {        
            //printf("faceMarkerIndex = %d\n", faceMarkerIndex);
            /*for(int k = 0; k < 3; k++)
                printf("triface list = %d\t", out.trifacelist[ i*3 + k ]);
            
            printf("\n");
            */    
            
            tetIndex = out.adjtetlist[i*2];
            //printf("tetIndex = %d\n", tetIndex); 

            countTri++;   
            getTetrahedronInformation(out, tetraCornerIndices, order, tetIndex, i);
            getTetrahedronNodes(out, element, tetraCornerNodes, tetraCornerIndices);
            
         
            // run over each triface
            for(index1 = 0; index1 < 3 ;index1++)
            {                   
                index2 = index1+1;
                if(index2 > 2) index2 = 0;
                
                //printf("edgeIndex1 = %d\n", out.tetrahedronlist[tetIndex*out.numberofcorners+order[index1]]);
                //printf("edgeIndex2 = %d\n", out.tetrahedronlist[tetIndex*out.numberofcorners+order[index2]]);
                
                localHigherOrderIndex = getHigherOrderIndex(order[index1], order[index2], DRT::Element::tet10); 
                //printf("localHOindex = %d\n", localHigherOrderIndex);
                globalHigherOrderIndex = out.tetrahedronlist[tetIndex*out.numberofcorners+localHigherOrderIndex];
                //printf("globalHOindex = %d\n", globalHigherOrderIndex);
                if(visitedPointIndexList[globalHigherOrderIndex]== 0)
                {                            
                    visitedPointIndexList[globalHigherOrderIndex] = 1; 
                         
                    if(facetOnExternalBoundary(index1, index2, 3, tetIndex, out))       currentDomain = false;
                    else                                                                currentDomain = true;
                       
                    computeIntersectionNormal(  element, index1, index2, globalHigherOrderIndex, out, 
                                                tetraCornerIndices, order, tetraCornerNodes, plane, currentDomain);
                
                    
                    intersected = computeRecoveryNormal(surfaceNodes[faceMarkerIndex], intersectedSurfaces_[faceMarkerIndex], 
                                                        plane, xsi, currentDomain);
              
                    
                    if(intersected)
                    {   
                        printf("\n");
                        printf("INTERSECTED NORMAL\n");
                        printf("\n");
                        storeHigherOrderNode(   globalHigherOrderIndex, xsi, surfaceNodes[faceMarkerIndex], plane, -1,
                                                intersectedSurfaces_[faceMarkerIndex], element, out, true, true);
                        
                    }
                    else
                    {
                        int lineIndex = -1;
                        int surfaceIndex = -1;
                        int faceMarkerIndex2;
                        // find other facemarker
                        // determine common lineElement - the intersection line

                        if(currentDomain)
                        {
                            printf("CURRENT DOMAIN\n");
                            faceMarkerIndex2 = findAdjacentFaceMarker(  tetraCornerIndices[index1], 
                                                                        tetraCornerIndices[index2], 
                                                                        faceMarkerIndex, out); 
                            findCommonSurfaceEdge(faceMarkerIndex, faceMarkerIndex2, surfaceIndex, lineIndex); 
                        }
                        else
                        { 
                            printf("REFERNCE DOMAIN\n");
                            surfaceIndex = faceMarkerIndex;
                            lineIndex = findIntersectingSurfaceEdge(element, intersectedSurfaces_[surfaceIndex],
                                                                    tetraCornerNodes[index1], 
                                                                    tetraCornerNodes[index2]);
                            printf("lineIndex = %d\n", lineIndex);
                            
                        }
                        
                        intersected = computeRecoveryPlane( surfaceNodes[surfaceIndex], intersectedSurfaces_[surfaceIndex], 
                                                            lineIndex, plane, xsi, true);
                       
                        if(intersected)
                        {   
                            printf("\n");
                            printf("INTERSECTED PLANE\n");
                            printf("\n");
                            storeHigherOrderNode(   globalHigherOrderIndex, xsi, surfaceNodes[faceMarkerIndex], plane,
                                                    lineIndex, intersectedSurfaces_[faceMarkerIndex], 
                                                    element, out, true, false);
                        
                        }
                        else dserror("NO INTERSECTION POINT FOUND\n");
                                
                    }       
                }
            }
        }
    }
    
 
    //printf("tri = %d\n", countTri);
    delete [] visitedPointIndexList;
    visitedPointIndexList = (int *) NULL;
    intersectedSurfaces_.clear();
}



/*----------------------------------------------------------------------*
 |  RQI:    checks if all tetrahedra corner points are lying u.may 09/07|
 |          in a surface element                                        |   
 |          if not corner points is recovered on the surface element    | 
 *----------------------------------------------------------------------*/
void Intersection::liftAllCornerPoints(
    tetgenio&                                       out,
    DRT::Element*                                   xfemElement,
    vector< vector<Epetra_SerialDenseVector> >&     surfaceNodes)
{
    bool alreadyInList = false;
    bool currentDomain = true;
    int pointIndex, lineIndex, surfaceIndex;
    unsigned int size;
    vector< vector <int> > adjacentFacesList;
    vector< vector <int> > adjacentFacemarkerList;
    Epetra_SerialDenseVector x(3);
    Epetra_SerialDenseVector xsi(3);
    
    for(int i = 0; i < out.numberoftrifaces; i++)
    {   
        if( (out.trifacemarkerlist[i]-facetMarkerOffset_) > -1)
        {  
            for (int j = 0; j < 3; j++)
            {
                pointIndex = out.trifacelist[i*3+j];
        
                alreadyInList = false;
                
                // check if point is a Steiner point
                if(out.pointmarkerlist[pointIndex] != 2 && out.pointmarkerlist[pointIndex] != 3 )
                {
                    //printf("pointIndex = %d\n", pointIndex);
                    size = adjacentFacesList.size();
                    vector<int> pointIndices(2);
                    pointIndices = getPointIndices(out, i, j);
                    //printf("pointIndex2 = %d\n", pointIndices[0]);
                    //printf("pointIndex3 = %d\n", pointIndices[1]);
                    
                    for(unsigned int k = 0; k < size; k++)
                    {
                        if(adjacentFacesList[k][0] == pointIndex)
                        {
                            alreadyInList = true;
                            
                            adjacentFacesList[k].push_back(pointIndices[0]);
                            adjacentFacesList[k].push_back(pointIndices[1]);
                            adjacentFacemarkerList[k].push_back( out.trifacemarkerlist[i]-facetMarkerOffset_ );
                            break;
                        }
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
    
    bool intersected = false;
    bool normalSteiner = true;
    int length = 0;
    int pointIndex1;
    int pointIndex2;
    Epetra_SerialDenseVector Steinerpoint(3);
    Epetra_SerialDenseVector p1(3);
    Epetra_SerialDenseVector p2(3);
    Epetra_SerialDenseVector node(3);
    Epetra_SerialDenseVector normal(3);
    Epetra_SerialDenseVector averageNormal(3); 
    vector<Epetra_SerialDenseVector> normals; 
    vector<Epetra_SerialDenseVector> plane(4); 
    Epetra_SerialDenseVector n1(3); 
    Epetra_SerialDenseVector n2(3);
    Epetra_SerialDenseVector r1(3); 
    Epetra_SerialDenseVector r2(3);
    Epetra_SerialDenseVector edgePoint(3);
    Epetra_SerialDenseVector oppositePoint(3);
    DRT::Element* lineElement;
     
    if(!adjacentFacesList.empty())
    {
        printf("number of steiner points = %d\n",adjacentFacesList.size());
        for(unsigned int i=0; i<adjacentFacesList.size(); i++)
        {
            // check if steiner point is already on the surface element
            pointIndex = adjacentFacesList[i][0];
            for(int k=0; k<3; k++)
            {
                x[k]   = out.pointlist[pointIndex*3 + k];
                xsi[k] = 0.0;
            }
       
            InterfacePoint emptyIp;
            checkNodeWithinElement(xfemElement, x, xsi);
            if(checkIfOnBoundary(xsi, emptyIp))     out.pointmarkerlist[pointIndex] = 3;    // on xfem boundary
            else                                    out.pointmarkerlist[pointIndex] = 2;    // not on xfem boundary
            
            normalSteiner = true;
            
            for(unsigned int j = 0; j < adjacentFacemarkerList[i].size(); j++ )
            {
                for(unsigned int k = 0; k < adjacentFacemarkerList[i].size(); k++ )
                {
                    if(adjacentFacemarkerList[i][j] != adjacentFacemarkerList[i][k])
                    {
                        //printf("a = %d and b = %d\n", adjacentFacemarkerList[i][j], adjacentFacemarkerList[i][k]);
                        if(findCommonFaceEdge(  j, k, 
                                                adjacentFacesList[i], edgePoint, oppositePoint, out))
                        {
                            if(!findCommonSurfaceEdge(  adjacentFacemarkerList[i][j], adjacentFacemarkerList[i][k],
                                                    surfaceIndex,  lineIndex))
                              
                                                    dserror("no common line element found\n");
                            normalSteiner = false;
                        }
                    }
                    if(!normalSteiner)      break;
                }
                if(!normalSteiner)      break;
            }
            
            // get Steiner point coordinates 
            for(int j=0; j<3; j++)
            {
                averageNormal[j] = 0.0;
                Steinerpoint[j] = out.pointlist[adjacentFacesList[i][0]*3 + j];
                //printf("face marker = %d\n", adjacentFacesList[i][0]);
                printf("Steinerpoint = %10.8f\t", Steinerpoint[j]);
            }
            printf("\n");
            
            if(normalSteiner)
            {
                if(!normals.empty())    normals.clear();
                // TODO check if on boundary if yes set marker list = 3
                referenceToCurrentCoordinates(xfemElement, Steinerpoint);   
           
                // printf("\n");
                length = (int) ( ( (double) (adjacentFacesList[i].size()-1))*0.5 );
                //printf("length = %d\n", length);
            
                for(int j = 0; j < length; j++)
                {
                    pointIndex1 = adjacentFacesList[i][1 + 2*j];
                    pointIndex2 = adjacentFacesList[i][1 + 2*j + 1];
                    for(int k = 0; k < 3; k++)
                    {
                        p1[k] =  out.pointlist[pointIndex1*3 + k];
                        p2[k] =  out.pointlist[pointIndex2*3 + k];
                    }
                
                    referenceToCurrentCoordinates(xfemElement, p1);   
                    referenceToCurrentCoordinates(xfemElement, p2);   
                    n1 = subtractsTwoVectors(p1, Steinerpoint);
                    n2 = subtractsTwoVectors(p2, Steinerpoint);
                
                    normal = computeCrossProduct( n1, n2);
                    normal = normalizeVector(normal);
                
                    for(int k=0; k<3; k++)
                        averageNormal[k] += normal[k];
                        
                     normals.push_back(normal);
                }
            
                // compute average normal
                for(int j=0; j<3; j++)
                    averageNormal[j] = averageNormal[j]/( (double)length);
                   
                
                                  
                plane[0] = addTwoVectors(Steinerpoint, averageNormal);               
                plane[1] = subtractsTwoVectors(Steinerpoint, averageNormal);
         
                int faceMarker = adjacentFacemarkerList[i][0];
                intersected = computeRecoveryNormal(surfaceNodes[faceMarker], intersectedSurfaces_[faceMarker], 
                                                    plane, xsi, true);
                if(intersected)
                {
                    printf("LIFT STEINER POINT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                    storeHigherOrderNode(   adjacentFacesList[i][0], xsi, surfaceNodes[faceMarker], plane,
                                            -1, intersectedSurfaces_[faceMarker], xfemElement, out, true, true);
                }
                else
                {
                    for(unsigned int j = 0; j < normals.size(); j++ )
                    {
                        plane[0] = addTwoVectors(Steinerpoint, normals[j]);               
                        plane[1] = subtractsTwoVectors(Steinerpoint, normals[j]);
                        intersected = computeRecoveryNormal(    surfaceNodes[faceMarker], intersectedSurfaces_[faceMarker], 
                                                                plane, xsi, true);
                        if(intersected)
                        {
                            printf("LIFT STEINER POINT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                            storeHigherOrderNode(   adjacentFacesList[i][0], xsi, surfaceNodes[faceMarker], plane,
                                                    -1, intersectedSurfaces_[faceMarker], xfemElement, out, true, true);
                            break;
                        }
                    }
                    if(!intersected) dserror("STEINER POINT NOT LIFTED");
                }
            }
            else
            {
                referenceToCurrentCoordinates(xfemElement, Steinerpoint);   
                referenceToCurrentCoordinates(xfemElement, edgePoint);   
                referenceToCurrentCoordinates(xfemElement, oppositePoint);  
                 
                r1 = subtractsTwoVectors( edgePoint, Steinerpoint);
                r2 = subtractsTwoVectors( oppositePoint, Steinerpoint);
            
                n1 = computeCrossProduct( r1, r2);
                n2 = computeCrossProduct( r1, n1);
                
                n1 = normalizeVector(n1);
                n2 = normalizeVector(n2);
            
                plane[0] = addTwoVectors(Steinerpoint, n1);               
                plane[1] = subtractsTwoVectors(Steinerpoint, n1);
                plane[2] = addTwoVectors(plane[1], n2);
                plane[3] = addTwoVectors(plane[0], n2);
                
                intersected = computeRecoveryPlane( surfaceNodes[surfaceIndex], intersectedSurfaces_[surfaceIndex], 
                                                    lineIndex, plane, xsi, true);
                if(intersected)
                {   
                    storeHigherOrderNode(   adjacentFacesList[i][0], xsi, surfaceNodes[surfaceIndex], plane,
                                            lineIndex, intersectedSurfaces_[surfaceIndex], xfemElement, out, true, false);
                }
                else
                    dserror("STEINER POINT NOT LIFTED");   
       
            }
        }
    }
}
            


/*----------------------------------------------------------------------*
 |  RQI:    returns the other two point indices belonging    u.may 09/07|
 |          to a triface that obtains a Steiner point                   |           
 *----------------------------------------------------------------------*/
vector<int> Intersection::getPointIndices(
    tetgenio&   out, 
    int         trifaceIndex, 
    int         steinerPointIndex)
{
    
    int count = 0;
    vector<int> pointIndices(2);
    
    for(int i = 0; i < 3; i++)
        if(i != steinerPointIndex)
            pointIndices[count++] = out.trifacelist[trifaceIndex*3+i];
           
    return pointIndices;
}



/*----------------------------------------------------------------------*
 |  RQI:    computes the intersection between a              u.may 08/07|
 |          line  and a surface                    RQI                  |           
 *----------------------------------------------------------------------*/
bool Intersection::computeRecoveryNormal( 
    vector <Epetra_SerialDenseVector>&          surfaceNodes,
    DRT::Element*                               surfaceElement,
    const vector<Epetra_SerialDenseVector>&     plane,
    Epetra_SerialDenseVector&                   xsi,
    const bool                                  current)
{
    int iter = 0;
    int countSingular = 0;
    int maxiter = 100;
    bool intersection = true;
    double residual = 1.0;
    Epetra_SerialDenseMatrix A(3,3);
    Epetra_SerialDenseVector b(3);
    Epetra_SerialDenseVector dx(3);
    
    for(int i = 0; i < 3; i++)
        xsi[i]= 0.0; 
        
    
    dx = xsi;
    
    
    updateRHSForRQINormal( b, xsi, surfaceNodes, surfaceElement, plane, current);
                                
    while(residual > TOL14_)
    {   
        updateAForRQINormal( A, xsi, plane, surfaceNodes, surfaceElement, current);
         
        if(!solveLinearSystemWithSVD(A, b, dx, 3))
        {  
            countSingular++;
            printf("MATRIX SINGULAR step %d\n", countSingular );
        } 
        
        if(countSingular > 5)
        {
            //printf("MATRIX SINGULAR\n");
            intersection = false;
            break;  
        }
        xsi = addTwoVectors(xsi,dx);
        //printf("dx0 = %20.16f\t, dx1 = %20.16f\t, dx2 = %20.16f\n", dx[0], dx[1], dx[2]);
        if(iter >= maxiter)
        {   
            printf("ITERATION > 20\n");
            intersection = false;
            break;
        }       
        
        updateRHSForRQINormal( b, xsi, surfaceNodes, surfaceElement, plane, current); 
        residual = b.Norm2(); 
        iter++;
        
        //printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi[0], xsi[1], xsi[2], residual, TOL14_);
    } 
    
    if( (fabs(xsi[0])-1.0) > TOL7_  || (fabs(xsi[1])-1.0) > TOL7_ )    // line coordinate may be bigger than 1
        intersection = false;
        
    return intersection;
}




/*----------------------------------------------------------------------*
 |  RQI:    updates the systemmatrix for the                 u.may 08/07|
 |          computation of a curve-surface intersection                 |
 |          for the recovery of the curved interface                    |
 |          (line-surface intersection)                                 |   
 *----------------------------------------------------------------------*/
void Intersection::updateAForRQINormal(   
    Epetra_SerialDenseMatrix&                   A,
    const Epetra_SerialDenseVector&             xsi,
    const vector<Epetra_SerialDenseVector>&     plane,
    const vector <Epetra_SerialDenseVector>&    surfaceNodes,
    DRT::Element*                               surfaceElement,
    const bool                                  current)                                                 
{   
    const int numNodesSurface = surfaceNodes.size();
    const int numNodesLine = 3;
    vector<int> actParams(1,0);
    Epetra_SerialDenseMatrix surfaceDeriv1(2,numNodesSurface);
    Epetra_SerialDenseMatrix lineDeriv1(1,numNodesLine);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    Epetra_SerialDenseVector xsiSurface(2);
    Epetra_SerialDenseVector xsiLine(1);
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
     
    xsiSurface(0) = xsi[0]; // r-coordinate surface
    xsiSurface(1) = xsi[1]; // s-coordinate surface
    xsiLine(0)    = xsi[2]; // r-coordinate line
    actParams[0] = numNodesSurface;     
    
    for(unsigned int i=0; i<3; i++)
        for(unsigned int j=0; j<3; j++)
            A[i][j] = 0.0;
    
    params.set("action","calc_ShapeDeriv1");
    surfaceElement->Evaluate(params, dummyDis, actParams, surfaceDeriv1, emptyM , emptyV, xsiSurface, emptyV);          
    //printf("num of nodes surfaces = %d\n", numOfNodesSurface);
        
    if(current)
    {
        for(int dim=0; dim<3; dim++)
        {
            for(int i=0; i<numNodesSurface; i++)
            {
                A[dim][0] += surfaceElement->Nodes()[i]->X()[dim] * surfaceDeriv1[i][0];
                A[dim][1] += surfaceElement->Nodes()[i]->X()[dim] * surfaceDeriv1[i][1];
            }
            A[dim][2] = (-1.0) * ( plane[0][dim] * (-0.5) + plane[1][dim] * 0.5 );  
        }
    }
    else
    {
        shape_function_1D_deriv1(lineDeriv1, xsiLine[0], DRT::Element::line3);
     
        for(int dim=0; dim<3; dim++)
        {
            for(int i=0; i<numNodesSurface; i++)
            {
                A[dim][0] += surfaceElement->Nodes()[i]->X()[dim] * surfaceDeriv1[i][0];
                A[dim][1] += surfaceElement->Nodes()[i]->X()[dim] * surfaceDeriv1[i][1];
            }
            
            for(int i = 0; i < 3; i++)
            {
                int index = i;
                if(i > 1)   index = 4;
                A[dim][2] += (-1.0) * plane[index][dim] * lineDeriv1[i][0]; 
            } 
        }
    }
}


/*
 
 
 for(int dim=0; dim<3; dim++)
        {
            for(int i=0; i<numNodesSurface; i++)
            {
                A[dim][0] += surfaceNodes[i][dim] * surfaceDeriv1[i][0];
                A[dim][1] += surfaceNodes[i][dim] * surfaceDeriv1[i][1];
            }
            A[dim][2] = (-1.0) * ( plane[0][dim] * (-0.5) + plane[1][dim] * 0.5 );  
        }
 
 
 */



/*----------------------------------------------------------------------*
 |  RQI:    updates the right-hand-side for the              u.may 08/07|
 |          computation of a curve-surface intersection                 |     
 |          for the recovery of the curved surface (RQI)                |
 |          (line-surface intersection)                                 |   
 *----------------------------------------------------------------------*/
void Intersection::updateRHSForRQINormal( 
    Epetra_SerialDenseVector&                   b,
    Epetra_SerialDenseVector&                   xsi,    
    const vector <Epetra_SerialDenseVector>&    surfaceNodes,
    DRT::Element*                               surfaceElement,
    const vector<Epetra_SerialDenseVector>&     plane,
    const bool                                  current)                                                    
{
    const int numNodesSurface = surfaceNodes.size();
    const int numNodesLine = 3;
    vector<int> actParams(1,0);
    Epetra_SerialDenseVector surfaceFunct(numNodesSurface);
    Epetra_SerialDenseVector lineFunct(numNodesLine);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    Epetra_SerialDenseVector xsiSurface(2);
    Epetra_SerialDenseVector xsiLine(1);
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
    
    
    xsiSurface(0) = xsi[0]; // r-coordinate surface
    xsiSurface(1) = xsi[1]; // s-coordinate surface
    xsiLine(0)    = xsi[2]; // r-coordinate line
    params.set("action","calc_Shapefunction");
            
    for(unsigned int i=0; i<3; i++)   
        b[i] = 0.0;
    
    actParams[0] = numNodesSurface;            
    surfaceElement->Evaluate(params, dummyDis, actParams, emptyM, emptyM , surfaceFunct, xsiSurface, emptyV);           
  
    if(current)
    {
        for(int dim=0; dim<3; dim++)
            for(int i=0; i<numNodesSurface; i++)
            {
                b[dim] += (-1.0) * surfaceElement->Nodes()[i]->X()[dim] * surfaceFunct[i];
            }
            
        for(int dim=0; dim<3; dim++)
            b[dim] += plane[0][dim] * 0.5*(1.0 - xsi[2]) + plane[1][dim] * 0.5*(1.0 + xsi[2]);
    }
    else
    {
        shape_function_1D( lineFunct, xsiLine[0], DRT::Element::line3);
       
        for(int dim=0; dim<3; dim++)
            for(int i=0; i<numNodesSurface; i++)
            {
                b[dim] += (-1.0) * surfaceElement->Nodes()[i]->X()[dim] * surfaceFunct[i];
            }
        
        
        for(int dim=0; dim<3; dim++)
            for(int i=0; i<numNodesLine; i++)
            {
                int index = i;
                if(i > 1)   index = 4;
                
                b[dim] += plane[index][dim] * lineFunct[i]; 
            } 
    }
          
    
/*
for(int dim=0; dim<3; dim++)
            for(int i=0; i<numNodesSurface; i++)
            {
                b[dim] += (-1.0) * surfaceNodes[i][dim] * surfaceFunct[i];
            }
*/

  
}



/*----------------------------------------------------------------------*
 |  RQI:    stores a pointer to each intersecting            u.may 08/07|
 |          cutter element  used for the recovery of curved interface   |     
 *----------------------------------------------------------------------*/  
void Intersection::storeIntersectedCutterElement(
    DRT::Element* surfaceElement)
{
    bool alreadyInList = false;
    vector< DRT::Element* >::iterator iter;
  
    for(unsigned int i = 0; i < intersectedSurfaces_.size(); i++)
        if(intersectedSurfaces_[i] == surfaceElement)
        {
            alreadyInList = true;
            break; 
        }
      
    if(!alreadyInList)
        intersectedSurfaces_.push_back(surfaceElement);
     
}   



/*----------------------------------------------------------------------*
 |  RQI:    returns information of the tetrahedra            u.may 08/07|
 *----------------------------------------------------------------------*/  
void Intersection::getTetrahedronInformation(   
    const tetgenio&     out, 
    vector<int>&        tetraCornerIndices,
    vector<int>&        order,
    const int           tetIndex,
    const int           index)
{            
      
    int nodeIndex = 0;         
    
    // store boundary face node indices
    for(int j=0; j<3; j++)
        tetraCornerIndices[j] = out.trifacelist[index*3+j];
   
    // store node index opposite to the boundary face of the tetrahedron
    for(int j=0; j<4; j++)
    {
        nodeIndex = out.tetrahedronlist[tetIndex*out.numberofcorners+j];
        if(nodeIndex != tetraCornerIndices[0] && nodeIndex != tetraCornerIndices[1] && nodeIndex != tetraCornerIndices[2])
        {
            tetraCornerIndices[3] = out.tetrahedronlist[tetIndex*out.numberofcorners+j];
            break;
        }
    }
    
    for(int j=0; j<4; j++)
    {
        nodeIndex = out.tetrahedronlist[tetIndex*out.numberofcorners+j];
        for(int k=0; k<3; k++)
            if(nodeIndex == tetraCornerIndices[k])
            {
                order[k] = j;
                break;
            }
    }
}


/*----------------------------------------------------------------------*
 |  RQI:    collects the tetrahedron corner nodes            u.may 09/07|
 |          transforms them into current coordinates                    |
 |          of the xfem-element                                         |
 *----------------------------------------------------------------------*/  
void Intersection::getTetrahedronNodes(
    const tetgenio&                         out, 
    DRT::Element*                           xfemElement,
    vector<Epetra_SerialDenseVector>&       tetraCornerNodes,
    vector<int>&                            tetraCornerIndices)
{

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 3; j++)
            tetraCornerNodes[i][j] = out.pointlist[tetraCornerIndices[i]*3+j];
        
        referenceToCurrentCoordinates(xfemElement, tetraCornerNodes[i]);   
    }
}



/*----------------------------------------------------------------------*
 |  RQI:    computes the normal to the interface edge of     u.may 08/07|
 |          the tetrahedron facet lying within this facet               |
 *----------------------------------------------------------------------*/  
void Intersection::computeIntersectionNormal(    
    const int                           index1,
    const int                           index2,
    const int                           globalHigherOrderIndex, 
    const tetgenio&                     out,
    const vector<int>&                  tetraCornerIndices,
    vector<Epetra_SerialDenseVector>&   plane)
{            
    
    Epetra_SerialDenseVector  p1(3);
    Epetra_SerialDenseVector  p2(3);
    Epetra_SerialDenseVector  p3(3);        
    Epetra_SerialDenseVector  m(3);
    Epetra_SerialDenseVector  n(3);
    Epetra_SerialDenseVector  r(3);
    Epetra_SerialDenseVector  r1(3);
    Epetra_SerialDenseVector  r2(3);
    
    
    for(int i=0; i<3; i++)
    {
        p1[i] = out.pointlist[tetraCornerIndices[3]*3 + i];
        p2[i] = out.pointlist[tetraCornerIndices[index1]*3 + i];
        p3[i] = out.pointlist[tetraCornerIndices[index2]*3 + i];                
    }
    
    
/*    for(int i=0; i<3; i++)
        printf("p1 = %f\t", p1[i]);
        
    printf("\n"); printf("\n");
    
    for(int i=0; i<3; i++)
        printf("p2 = %f\t", p2[i]);
        
    printf("\n"); printf("\n");
    
    for(int i=0; i<3; i++)
        printf("p3 = %f\t", p3[i]);
        
    printf("\n"); printf("\n");    
    
 */                
    // compute direction vectors of the plane 
    r1 = subtractsTwoVectors(p1, p2);
/*
    for(int i=0; i<3; i++)
        printf("r1 = %f\t", r1[i]);
        
    printf("\n"); printf("\n");       
 */   
    r2 = subtractsTwoVectors(p3, p2);
   // r2 = normalizeVector(r2);
/*    for(int i=0; i<3; i++)
        printf("r2 = %f\t", r2[i]);
        
    printf("\n"); printf("\n");  
*/    
    // normal of the plane
    n = computeCrossProduct(r1, r2);
    n = normalizeVector(n);
 /*   for(int i=0; i<3; i++)
        printf("n = %f\t", n[i]);
        
    printf("\n"); printf("\n");
*/     
    // direction vector of the intersection line
    r = computeCrossProduct(n, r2);  
    r = normalizeVector(r);
 
 /*   for(int i=0; i<3; i++)
        printf("r = %f\t", r[i]);
        
    printf("\n"); printf("\n");
    */                 
    // computes the start point of the line
    m = computeLineMidpoint(p2, p3);
   // for(int i=0; i<3; i++)
   //     printf("m = %f\t", m[i]);
        
   // printf("\n"); printf("\n");      
    
    for(int i = 0; i < 3; i++)
        m[i] = out.pointlist[globalHigherOrderIndex*3+i];

   /* for(int i=0; i<3; i++)
        printf("m1 = %f\t", m[i]);
        
    printf("\n"); printf("\n");      
   */
    
    //dis = subtractsTwoVectors(p1, m);
    
    //dis = normalizeVector(dis);
    //for(int i=0; i<3; i++)
    //    printf("dis  = %f\t", dis[i]);
        
   // printf("\n"); printf("\n");      
    
   // absVal = computeAbsoluteValue(dis);
    
    // printf("absVal = %f\n",absVal);

    //absVal = absVal*100;
    //absDirection = computeScalarVectorMultiplication(absVal, r);
    //for(int i=0; i<3; i++)
    //    printf("absDirection  = %f\t", absDirection[i]);
        
    //printf("\n"); printf("\n");   
    
    // compute nodes of the normal to the interface edge of the tetrahedron
    plane[0] = addTwoVectors(m, r);               
    plane[1] = subtractsTwoVectors(m, r);
    plane[2] = addTwoVectors(plane[1], n);               
    plane[3] = addTwoVectors(plane[0], n);
}



/*----------------------------------------------------------------------*
 |  RQI:    computes the normal to the interface edge of     u.may 08/07|
 |          the tetrahedron facet lying within this facet               |
 *----------------------------------------------------------------------*/  
void Intersection::computeIntersectionNormal( 
    DRT::Element*                           xfemElement,   
    const int                               index1,
    const int                               index2,
    const int                               globalHigherOrderIndex, 
    const tetgenio&                         out,
    const vector<int>&                      tetraCornerIndices,
    const vector<int>&                      order,
    const vector<Epetra_SerialDenseVector>& tetraCornerNodes,
    vector<Epetra_SerialDenseVector>&       plane,
    const bool                              current)
{            
    
    Epetra_SerialDenseVector  p1(3);
    Epetra_SerialDenseVector  p2(3);
    Epetra_SerialDenseVector  p3(3);        
    Epetra_SerialDenseVector  m(3);
    Epetra_SerialDenseVector  n(3);
    Epetra_SerialDenseVector  r(3);
    Epetra_SerialDenseVector  r1(3);
    Epetra_SerialDenseVector  r2(3);
    
    
    if(current)
    {
        for(int i=0; i<3; i++)
        {    
            p1[i] = tetraCornerNodes[3][i];
            p2[i] = tetraCornerNodes[index1][i];
            p3[i] = tetraCornerNodes[index2][i];             
        }
    }
    else
    {
        for(int i=0; i<3; i++)
        {
            p1[i] = out.pointlist[tetraCornerIndices[3]*3       + i];
            p2[i] = out.pointlist[tetraCornerIndices[index1]*3  + i];
            p3[i] = out.pointlist[tetraCornerIndices[index2]*3  + i];                       
        }
    }
 
/*    
    for(int i=0; i<3; i++)
        printf("p1 = %f\t", p1[i]);
        
    printf("\n"); printf("\n");
    
    for(int i=0; i<3; i++)
        printf("p2 = %f\t", p2[i]);
        
    printf("\n"); printf("\n");
    
    for(int i=0; i<3; i++)
        printf("p3 = %f\t", p3[i]);
        
    printf("\n"); printf("\n");    

*/
                 
    // compute direction vectors of the plane 
    r1 = subtractsTwoVectors(p1, p2);
/*
    for(int i=0; i<3; i++)
        printf("r1 = %f\t", r1[i]);
        
    printf("\n"); printf("\n");       
 */   
    r2 = subtractsTwoVectors(p3, p2);
   // r2 = normalizeVector(r2);
/*    for(int i=0; i<3; i++)
        printf("r2 = %f\t", r2[i]);
        
    printf("\n"); printf("\n");  
*/    
    // normal of the plane
    n = computeCrossProduct(r1, r2);
    n = normalizeVector(n);
 /*   for(int i=0; i<3; i++)
        printf("n = %f\t", n[i]);
        
    printf("\n"); printf("\n");
*/     
    // direction vector of the intersection line
    r = computeCrossProduct(n, r2);  
    r = normalizeVector(r);
 
 /*   for(int i=0; i<3; i++)
        printf("r = %f\t", r[i]);
        
    printf("\n"); printf("\n");
    */                 
    // computes the start point of the line
    m = computeLineMidpoint(p2, p3);
    
    if(current)
        m = computeLineMidpoint(p2, p3);
    else
    {
        for(int i = 0; i < 3; i++)
            m[i] = out.pointlist[globalHigherOrderIndex*3+i];
    }
    
/*    for(int i=0; i<3; i++)
        printf("m = %f\t", m[i]);
        
    printf("\n"); printf("\n");
*/   
    //dis = subtractsTwoVectors(p1, m);
    
    //dis = normalizeVector(dis);
    //for(int i=0; i<3; i++)
    //    printf("dis  = %f\t", dis[i]);
        
   // printf("\n"); printf("\n");      
    
   // absVal = computeAbsoluteValue(dis);
    
    // printf("absVal = %f\n",absVal);

    //absVal = absVal*100;
    //absDirection = computeScalarVectorMultiplication(absVal, r);
    //for(int i=0; i<3; i++)
    //    printf("absDirection  = %f\t", absDirection[i]);
        
    //printf("\n"); printf("\n");   
    
    // compute nodes of the normal to the interface edge of the tetrahedron
    plane[0] = addTwoVectors(m, r);               
    plane[1] = subtractsTwoVectors(m, r);
    plane[2] = addTwoVectors(plane[1], n);               
    plane[3] = addTwoVectors(plane[0], n);
    
    
    if(!current)
    {
        for(int i = 0; i < 4; i++)
            referenceToCurrentCoordinates(xfemElement, plane[i]);
        
        referenceToCurrentCoordinates(xfemElement, m);
          
        plane[4] = m;
      
    }
    
    
/*    for(int i = 0; i < 2; i++)
    {
        for(int j = 0; j < 3; j++)
            printf("plane = %f\n", plane[i][j]);
     
        printf("\n");       
    }
*/    

}



/*----------------------------------------------------------------------*
 |  RQI:    computes the midpoint of a line                  u.may 08/07|
 *----------------------------------------------------------------------*/  
Epetra_SerialDenseVector Intersection::computeLineMidpoint( 
    const Epetra_SerialDenseVector& p1, 
    const Epetra_SerialDenseVector& p2)
{
    Epetra_SerialDenseVector midpoint(3);
    
    for(int i=0; i<3; i++)
        midpoint[i] = (p1[i] + p2[i])*0.5;
        
    return midpoint;
}



/*----------------------------------------------------------------------*
 |  RQI:    transforms the nodes of the intersecting         u.may 08/07|
 |          surface elements into the reference configuration           |
 |          of the xfem element                                         |
 *----------------------------------------------------------------------*/  
void Intersection::transformSurfaceNodes(   
    DRT::Element*                                   element,
    vector < vector <Epetra_SerialDenseVector> >&   surfaceNodes)
{

    Epetra_SerialDenseVector coord(3);
    vector <Epetra_SerialDenseVector> nodesPerElement;
    
    
    surfaceNodes.clear();
    
    for(unsigned int i = 0; i < intersectedSurfaces_.size(); i++ )
    {  
        for(int j = 0; j < intersectedSurfaces_[i]->NumNode(); j++ )
        {
            for(int k = 0; k < 3; k++)
                coord[k] = intersectedSurfaces_[i]->Nodes()[j]->X()[k];
            
            currentToReferenceCoordinates(element, coord);
            nodesPerElement.push_back(coord);
        }
        surfaceNodes.push_back(nodesPerElement);
        nodesPerElement.clear();
    }
    
    
    
/*    for(unsigned int i = 0; i < intersectedSurfaces_.size(); i++ )
    {  
        for(int j = 0; j < intersectedSurfaces_[i]->NumNode(); j++ )
        {
            for(int k = 0; k < 3; k++)
                coord[k] = surfaceNodes[i][j][k];
    
            referenceToCurrentCoordinates(element, coord);
            
            intersectedSurfaces_[i]->Nodes()[j]->Print(cout);
            printf("\t"); 
            for(int k = 0; k < 3; k++)
                printf("coord = %5.2f\t", coord[k]);
                
            
            printf("\n");
            
        }
        printf("\n");
    }
    
    
    
    for(int i = 0; i < element->NumNode(); i++)
    {
        printf("%2.2d. node\t",i);
        for(int j = 0; j < 3; j++)
        {
            coord[j] = eleRefCoordinates_[i][j];
            printf("refCoord = %5.2f\t", coord[j]);
        }
        printf("\n");
        referenceToCurrentCoordinates(element, coord);
        element->Nodes()[i]->Print(cout);
        printf("\n");
        
        printf("\t");
        for(int j = 0; j < 3; j++) 
            printf("refcoord = %3.2f\t",coord[j]);
            
        printf("\n");
       
        currentToReferenceCoordinates(element, coord);
        
        printf("%2.2d. node\t",i);
        for(int j = 0; j < 3; j++)
            printf("curCoord = %19.16f\t", coord[j]);
            
        printf("\n");
        printf("\n");
    }
    printf("\n");
*/     
    
}




/*----------------------------------------------------------------------*
 |  RQI:    stores the higher-order node in the pointlist    u.may 08/07|
 |          at the place of the linear node                             |
 *----------------------------------------------------------------------*/  
void Intersection::storeHigherOrderNode( 
    const int                                   index,
    Epetra_SerialDenseVector&                   xsi, 
    const vector <Epetra_SerialDenseVector>&    surfaceNodes,
    const vector <Epetra_SerialDenseVector>&    plane,
    int                                         lineIndex, 
    DRT::Element*                               surfaceElement, 
    DRT::Element*                               xfemElement, 
    tetgenio&                                   out,
    const bool                                  current,
    const bool                                  normal )
{
    
    Epetra_SerialDenseVector xsiLine(3);
    for(int i = 0; i < 3; i++)
        xsiLine[i] = 0.0;
        
    if(current)
    {
        if(normal)     referenceToCurrentCoordinates(surfaceElement, xsi);
        else           
        {   
            xsiLine[0] = xsi[2];
            referenceToCurrentCoordinates(surfaceElement->Lines()[lineIndex], xsiLine);
            for(int i=0; i<3; i++)
                xsi[i] = xsiLine[i];
        }       
        currentToReferenceCoordinates(xfemElement, xsi);
    }
    else
    {
        if(normal)  referenceToCurrentCoordinates(surfaceElement, xsi, surfaceNodes);
        else        referenceToCurrentCoordinates(xsi, plane); 
    }
    
    printf("xsiold0 = %20.16f\t, xsiold1 = %20.16f\t, xsiold2 = %20.16f\n", out.pointlist[index*3], out.pointlist[index*3+1], out.pointlist[index*3+2]);
    
    for(int i = 0; i < 3; i++)
        out.pointlist[index*3+i]   = xsi[i];  
   
    printf("xsi0    = %20.16f\t, xsi1    = %20.16f\t, xsi2    = %20.16f\n", xsi[0], xsi[1], xsi[2]);
    printf("\n");
    
}



/*----------------------------------------------------------------------*
 |  RQI:    computes the intersection between a              u.may 08/07|
 |          curve and a plane                      RQI                  |           
 *----------------------------------------------------------------------*/
bool Intersection::computeRecoveryPlane( 
    vector <Epetra_SerialDenseVector>&          surfaceNodes,
    DRT::Element*                               surfaceElement,
    int&                                        lineIndex,
    const vector<Epetra_SerialDenseVector>&     plane,
    Epetra_SerialDenseVector&                   xsi,
    const bool                                  current)
{
    int                                 iter = 0;
    int                                 begin, end;
    int                                 index1, index2;
    int                                 maxiter = 100;
    int                                 numLines = surfaceElement->NumLine();
    int                                 numNodesLine = surfaceElement->Lines()[0]->NumNode();
    bool                                intersection = true;
    double                              residual = 1.0;
    Epetra_SerialDenseMatrix            A(3,3);
    Epetra_SerialDenseVector            b(3);
    Epetra_SerialDenseVector            dx(3);
    vector <Epetra_SerialDenseVector>   curve(numNodesLine, Epetra_SerialDenseVector(3));
    vector <vector<int> >               nodeIndices(numLines, vector<int>(numNodesLine));
    DRT::Element*                       lineElement = surfaceElement->Lines()[0];
   
   
    if(!current)
        // get node indices for each line of the surface element discretization type 
        for(int i = 0; i < numLines; i++)
            for(int j = 0; j < numNodesLine; j++)
                nodeIndices = getEleNodeNumberingLines( surfaceElement->Shape() );
      
            
    printf("numlines = %d\n", numLines );
    // run over all lines (curves)
    
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
    
    for(int i = begin; i < end; i++)
    {
        intersection = true;
        residual = 1.0;
        iter = 0;
        printf("num = %d\n", i );
        // starting value equals (0,0,0)
        for(int j = 0; j < 3; j++)
            xsi[j]= 0.0; 
        
        dx = xsi;
        
        if(current)
        {
            lineElement = surfaceElement->Lines()[i];
            lineElement->Print(cout);
            printf("\n");
        }   
        else
        {
            index1 = nodeIndices[i][0];
            index2 = nodeIndices[i][1];
                
            for(int dim = 0; dim < 3; dim++)
            {
                curve[0][dim] =  surfaceNodes[index1][dim];
                curve[1][dim] =  surfaceNodes[index2][dim];
                if(numNodesLine > 2)
                    curve[2][dim] =  surfaceNodes[ nodeIndices[i][2] ][dim];  
                    //surfaceNodes[getHigherOrderIndex(index1, index2, surfaceElement->Shape())][dim];  
            }    
        } 
        
     /*   for(int k = 0; k < 3; k++)
        {
            for(int j = 0; j < 3; j++)
                printf("curve = %f\n", curve[k][j]);
                
            printf("\n");       
        }
        
        for(int k = 0; k < 4; k++)
        {
            for(int j = 0; j < 3; j++)
                printf("plane = %f\n", plane[k][j]);
                
            printf("\n");       
        }
        */
        updateRHSForRQIPlane( b, xsi, curve, lineElement, plane, current);                
        while( residual > TOL14_ )
        {   
            updateAForRQIPlane( A, xsi, plane, curve, surfaceElement, lineElement, current);
            
            if(!gaussElimination(A, b, dx, true, 3, 1))
            {
                printf("MATRIX SINGULAR\n");
                intersection = false;
                break;  
            } 
             
            if(iter >= maxiter)
            {   
                printf("ITERATION > 20\n");
                intersection = false;
                break;
            }       
        
            xsi = addTwoVectors(xsi, dx);
            //printf("dx0 = %20.16f\t, dx1 = %20.16f\t, dx2 = %20.16f\n", dx[0], dx[1], dx[2]);
            updateRHSForRQIPlane( b, xsi, curve, lineElement, plane, current);
            residual = b.Norm2(); 
            iter++;
        
            printf("xsi0 = %20.16f\t, xsi1 = %20.16f\t, xsi2 = %20.16f\t, res = %20.16f\t, tol = %20.16f\n", xsi[0], xsi[1], xsi[2], residual, TOL14_);
        } 
    
        if( (fabs(xsi[2])-1.0) > TOL7_ )     // planes coordinate may be bigger than 1
            intersection = false;
            
        if(intersection)
        {   
            lineIndex = begin;
            break;
        }
    }
           
    return intersection;
}



/*----------------------------------------------------------------------*
 |  RQI:    updates the systemmatrix for the                 u.may 09/07|
 |          computation of a curve-surface intersection                 |
 |          for the recovery of the curved interface                    |
 |          (curve-plane intersection)                                  |   
 *----------------------------------------------------------------------*/
void Intersection::updateAForRQIPlane(   
    Epetra_SerialDenseMatrix&                   A,
    const Epetra_SerialDenseVector&             xsi,
    const vector<Epetra_SerialDenseVector>&     plane,
    const vector<Epetra_SerialDenseVector>&     curve,
    DRT::Element*                               surfaceElement,
    DRT::Element*                               lineElement,
    const bool                                  current)                                                 
{   
    const int numNodesLine = curve.size();
    const int numNodesSurface = 4;
    vector<int> actParams(1,0);
    Epetra_SerialDenseMatrix surfaceDeriv(2, numNodesSurface);
    Epetra_SerialDenseMatrix lineDeriv(1,numNodesLine);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    Epetra_SerialDenseVector xsiSurface(2);
    Epetra_SerialDenseVector xsiLine(1);
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
     
    xsiSurface[0] = xsi[0]; // r-coordinate surface
    xsiSurface[1] = xsi[1]; // s-coordinate surface
    xsiLine[0]    = xsi[2]; // r-coordinate line
    
    for(unsigned int i=0; i<3; i++)
        for(unsigned int j=0; j<3; j++)
            A[i][j] = 0.0;
    
    
    params.set("action","calc_ShapeDeriv1");
    //actParams[0] = surfaceElement->NumNode();  
    //surfaceElement->Evaluate(params, dummyDis, actParams, surfaceDeriv, emptyM , emptyV, xsiSurface, emptyV);   
    shape_function_2D_deriv1(surfaceDeriv,  xsiSurface[0],  xsiSurface[1], DRT::Element::quad4);
       
    actParams[0] = numNodesLine;    
    lineElement->Evaluate(params, dummyDis, actParams, lineDeriv, emptyM , emptyV, xsiLine, emptyV);               
    
      
    for(int dim=0; dim<3; dim++)
    {
        for(int i=0; i<numNodesSurface; i++)
        {
            A[dim][0] += plane[i][dim] * surfaceDeriv[i][0];
            A[dim][1] += plane[i][dim] * surfaceDeriv[i][1];
        }
    }
        
    if(current)
    {
        for(int dim=0; dim<3; dim++)
        {
            for(int i=0; i<numNodesLine; i++)
            {
                A[dim][2] += (-1.0) * lineElement->Nodes()[i]->X()[dim] * lineDeriv[i][0];   
            }
        }        
    }
    else
    {
        for(int dim=0; dim<3; dim++)
        {
            for(int i=0; i<numNodesLine; i++)
            {
                A[dim][2] += (-1.0) * curve[i][dim] * lineDeriv[i][0];   
            }
        }             
    }
}



/*----------------------------------------------------------------------*
 |  RQI:    updates the right-hand-side for the              u.may 09/07|
 |          computation of a curve-surface intersection                 |     
 |          for the recovery of the curved surface (RQI)                |
 |          (curve-plane intersection)                                  |   
 *----------------------------------------------------------------------*/
void Intersection::updateRHSForRQIPlane( 
    Epetra_SerialDenseVector&                   b,
    Epetra_SerialDenseVector&                   xsi,    
    const vector <Epetra_SerialDenseVector>&    curve,
    DRT::Element*                               lineElement,
    const vector<Epetra_SerialDenseVector>&     plane,
    const bool                                  current)                                                    
{
    const int numNodesLine    = curve.size();
    const int numNodesSurface = 4;
    vector<int> actParams(1,0);
    Epetra_SerialDenseVector surfaceFunct(numNodesSurface);
    Epetra_SerialDenseVector lineFunct(numNodesLine);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    Epetra_SerialDenseVector xsiSurface(2);
    Epetra_SerialDenseVector xsiLine(1);
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
    
    xsiSurface[0] = xsi[0]; // r-coordinate surface
    xsiSurface[1] = xsi[1]; // s-coordinate surface
    xsiLine[0]    = xsi[2]; // r-coordinate line
    params.set("action","calc_Shapefunction");
            
    for(unsigned int i=0; i<3; i++)   
        b[i] = 0.0;
    
    // shape function for the plane 
    shape_function_2D( surfaceFunct, xsiSurface[0], xsiSurface[1], DRT::Element::quad4 ); 
    // shape function for the curve
    actParams[0] = numNodesLine;   
    lineElement->Evaluate(params, dummyDis, actParams, emptyM, emptyV , lineFunct, xsiLine, emptyV);    
  
    for(int dim=0; dim<3; dim++)
        for(int i=0; i<numNodesSurface; i++)
        {
            b[dim] += (-1.0) * plane[i][dim] * surfaceFunct[i];
        }
        
    if(current)
    {
        for(int dim=0; dim<3; dim++)
            for(int i=0; i<numNodesLine; i++)
            {
                b[dim] +=  lineElement->Nodes()[i]->X()[dim]  * lineFunct[i];        
            }
    }
    else
    {
        for(int dim=0; dim<3; dim++)
            for(int i=0; i<numNodesLine; i++)
            {
                b[dim] += curve[i][dim] * lineFunct[i];        
            }
    }
    
}



/*----------------------------------------------------------------------*
 |  RQI:    checks if the considered tetrahedron facet       u.may 09/07|
 |          lies on a xfem boundary                                     |     
 *----------------------------------------------------------------------*/
bool Intersection::facetOnExternalBoundary(
    int         index1,
    int         index2,
    int         cornerIndex,
    int         tetIndex,
    tetgenio&   out)
{
    
    bool onBoundary = true;
    int  globalIndex;
    vector<int> localIndex(3,0);
    
    localIndex[0] = index1; 
    localIndex[1] = index2; 
    localIndex[2] = cornerIndex; 
       
    for(int i = 0; i < 3; i++)
    {
        globalIndex = out.tetrahedronlist[tetIndex*out.numberofcorners + localIndex[i]];
        //printf("marker = %d\t", out.pointmarkerlist[globalIndex]);
        if(out.pointmarkerlist[globalIndex] != 3)
        {
            onBoundary = false;
            break;
        }
    } 
    printf("\n");
     
    return onBoundary;
}




bool Intersection::findCommonFaceEdge( 
    int                                 faceIndex1, 
    int                                 faceIndex2, 
    vector<int>&                        adjacentFacesList,
    Epetra_SerialDenseVector&           edgePoint,
    Epetra_SerialDenseVector&           oppositePoint,
    tetgenio&                           out)
{
    bool edgeFound = false;
    int index;
    
    for(int i = 0; i < 2; i++)
    {
        for(int j = 0; j < 2; j++)
        {
            printf("pointIndex = %d, pointIndex = %d\n ", adjacentFacesList[faceIndex1*2+i+1], adjacentFacesList[faceIndex2*2+j+1]);
            if(adjacentFacesList[faceIndex1*2+i+1] == adjacentFacesList[faceIndex2*2+j+1])
            {
                if(i == 0)  index = 1;
                if(i == 1)  index = 0;
                
                for(int k = 0; k < 3; k++)
                {
                    edgePoint[k]        = out.pointlist[adjacentFacesList[faceIndex1*2+i+1]*3 + k];
                    oppositePoint[k]    = out.pointlist[adjacentFacesList[faceIndex1*2+index+1]*3 + k];
                }
                printf("edgeFound\n");
                edgeFound = true;
                break;
            }   
        }
        if(edgeFound)
            break; 
        else
            printf("edge not Found\n");      
    }
    
    return edgeFound;
}
                            
                    
                    
        
bool Intersection::findCommonSurfaceEdge(  
    int                                             faceIndex1, 
    int                                             faceIndex2,
    int&                                            surfaceIndex,
    int&                                            lineIndex)
{
    bool comparison = false;
    int numLines1 = intersectedSurfaces_[faceIndex1]->NumLine(); 
    int numLines2 = intersectedSurfaces_[faceIndex2]->NumLine(); 
    int numNodes  = intersectedSurfaces_[faceIndex2]->Lines()[0]->NumNode();
    DRT::Node* node1;
    DRT::Node* node2;
    DRT::Element*                                   lineElement2; 
    
    
/*    for(unsigned int i = 0; i < intersectedSurfaces_.size(); i++)
    {
        for(int j = 0; j < numLines1; j++)
        {
            intersectedSurfaces_[i]->Lines()[j]->Print(cout); 
            printf("\n");
        } 
        printf("\n");  
    }
*/    

    printf("helloooo\n");
    for(int i = 0; i < numLines1; i++)
    {
        for(int j = 0; j < numLines2; j++)
        {
            comparison = true;
            for(int k  = 0; k < numNodes; k++)
            {
                node1 = intersectedSurfaces_[faceIndex1]->Lines()[i]->Nodes()[k];
                node2 = intersectedSurfaces_[faceIndex2]->Lines()[j]->Nodes()[k];
                if(!comparePoints(node1, node2))
                    comparison = false;
            }  
             
            if(!comparison)
            {
                comparison = true;
                for(int k  = 0; k < numNodes; k++)
                {
                    if(k==2)
                    {
                        node1 = intersectedSurfaces_[faceIndex1]->Lines()[i]->Nodes()[k];
                        node2 = intersectedSurfaces_[faceIndex2]->Lines()[j]->Nodes()[k];
                    }
                    else
                    {
                        node1 = intersectedSurfaces_[faceIndex1]->Lines()[i]->Nodes()[k];
                        node2 = intersectedSurfaces_[faceIndex2]->Lines()[j]->Nodes()[1-k];
                    }
                    
                    if(!comparePoints(node1, node2))
                        comparison = false;
                }   
            }
            
            if(comparison)
            {
                lineIndex    =  i;
                surfaceIndex =  faceIndex1;
                printf("faceIndex1 = %d, lineindex = %d\n", faceIndex1, i);
                printf("faceIndex2 = %d, lineindex = %d\n", faceIndex2, j);
                intersectedSurfaces_[faceIndex2]->Lines()[j]->Print(cout);
                printf("\n");
                break;
            }
        }
        if(comparison) 
            break;
    }
    return comparison;
}



int Intersection::findAdjacentFaceMarker(
    int         edgeIndex1, 
    int         edgeIndex2, 
    int         faceMarkerIndex,
    tetgenio&   out)
{

    bool faceMarkerFound = false;
    int faceMarkerIndex2 = 0;
    int pointIndex = 0;
    int countPoints = 0;
    
    
    printf("edgeIndex1 = %d\n", edgeIndex1 );
    printf("edgeIndex2 = %d\n", edgeIndex2 );
    
    for(int i=0; i<out.numberoftrifaces; i++)
    {
        faceMarkerIndex2 = out.trifacemarkerlist[i] - facetMarkerOffset_;
        printf("faceMarkerIndex = %d\n", faceMarkerIndex2);
        if(faceMarkerIndex2 > -1)
        {        
            countPoints = 0;
            for(int j = 0; j < 3 ; j++)
            {
                pointIndex = out.trifacelist[ i*3 + j ];
                printf("pointIndex = %d\t", pointIndex);
                if(pointIndex == edgeIndex1 || pointIndex == edgeIndex2 )
                    countPoints++;
                    
                printf("\n");
            }
            
            if(countPoints == 2 && faceMarkerIndex != faceMarkerIndex2)
                faceMarkerFound = true;
            
        }
        if(faceMarkerFound)
            break;
    }
    
    printf("faceMarkerIndex found = %d\n", faceMarkerIndex2);
    return faceMarkerIndex2;
}


int Intersection::findIntersectingSurfaceEdge(
    DRT::Element*                       xfemElement,
    DRT::Element*                       surfaceElement,
    Epetra_SerialDenseVector            edgeNode1,
    Epetra_SerialDenseVector            edgeNode2)
{

    int lineIndex = -1;
    int countIndex = -1;
    int numLines = surfaceElement->NumLine();
    DRT::Element* lineElement;
    Epetra_SerialDenseVector xsi1(1);
    Epetra_SerialDenseVector xsi2(1);
    Epetra_SerialDenseVector x1(1);
    Epetra_SerialDenseVector x2(1);
      
    referenceToCurrentCoordinates(xfemElement, edgeNode1);
    referenceToCurrentCoordinates(xfemElement, edgeNode2);
    
    
    x1[0] = edgeNode1[0];
    x2[0] = edgeNode2[0];
    
    bool check1;
    bool check2;
    
    
    for(int i = 0; i < numLines; i++)
    {
        lineElement = surfaceElement->Lines()[i];
        
        check1 = checkNodeWithinElement( lineElement, x1, xsi1);
        printf("xsi1 = %f\n", xsi1[0]);
        check2 = checkNodeWithinElement( lineElement, x2, xsi2);
        printf("xsi2 = %f\n", xsi2[0]);
        countIndex++;
        if( check1  &&  check2 )
        {   
            lineIndex = countIndex;
            printf("LINEINDEX = %d\n", lineIndex);
            break;
        }
    }
    
    for(int i = 0; i < surfaceElement->Lines()[lineIndex]->NumNode(); i++ )
    {
        surfaceElement->Lines()[lineIndex]->Nodes()[i]->Print(cout); 
        printf("\n");  
    }
    return lineIndex;
}


/*----------------------------------------------------------------------*
 |  DB:     Debug only                                       u.may 06/07|
 *----------------------------------------------------------------------*/  
void Intersection::debugXAABBIntersection( 	Epetra_SerialDenseMatrix cutterXAABB, 
											Epetra_SerialDenseMatrix xfemXAABB,
											DRT::Element* cutterElement,
				 							DRT::Element* xfemElement,
				 							int noC,
				 							int noX)
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
void Intersection::debugNodeWithinElement(  DRT::Element* element,
                                            DRT::Node* node,
                                            Epetra_SerialDenseVector& xsi,
                                            int noE,
                                            int noN,
                                            bool within)
{
    int numnodes = element->NumNode();
    vector<int> actParams(1,0);
    Epetra_SerialDenseVector funct(numnodes);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    Epetra_SerialDenseVector x(3);
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
      
    params.set("action","calc_Shapefunction");
    actParams[0] = numnodes;   
      
    element->Evaluate(params, dummyDis, actParams, emptyM, emptyM , funct, xsi, emptyV);        
    
    for(int dim=0; dim<3; dim++)
        x(dim) = 0.0;
    
    for(int dim=0; dim<3; dim++)
        for(int i=0; i<numnodes; i++)
        {
            x(dim) += element->Nodes()[i]->X()[dim] * funct(i);
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
    cout << "   r = " << xsi[0] << "     s = " <<  xsi[1] << "     t = "  << xsi[2] << endl;
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
void Intersection::debugTetgenDataStructure(    DRT::Element*               element)
{    
    cout << endl;
    cout << "===============================================================" << endl;
    cout << "Debug Tetgen Data Structure " << endl;
    cout << "===============================================================" << endl;
    cout << endl;
    cout << "POINT LIST " << " :" << endl;
    cout << endl;
    Epetra_SerialDenseVector xsi(3);
    for(unsigned int i = 0; i < pointList_.size(); i++)
    {
        for(int j = 0; j< 3; j++)
        {
            xsi[j] = pointList_[i].coord[j];
        }
        referenceToCurrentCoordinates(element, xsi);
        
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
void Intersection::debugTetgenOutput( 	tetgenio& in,
										tetgenio& out, 
    									DRT::Element * element,
    									vector<int>& elementIds)
{
	char* tetgenIn = "tetgenPLC";
	char* tetgenOut = "tetgenMesh";
	char tetgenInId[30];
	char tetgenOutId[30];
		
	for(unsigned int i = 0; i < elementIds.size(); i++)
	{
		if(element->Id()== elementIds[i])
		{
			// change filename
			sprintf(tetgenInId,"%s%d", tetgenIn, elementIds[i]);
			sprintf(tetgenOutId,"%s%d", tetgenOut, elementIds[i]);
			
			// write piecewise linear complex
			in.save_nodes(tetgenInId);
	    	in.save_poly(tetgenInId);
	      
	    	// write tetrahedron mesh
	    	out.save_elements(tetgenOutId);
	    	out.save_nodes(tetgenOutId);
	    	out.save_faces(tetgenOutId);
	    	
	    	printf("Saving tetgen output for the %d.xfem element\n", elementIds[i]);
	    }
    }
}





/*----------------------------------------------------------------------*
 |  DB:     Debug only                                       u.may 06/07|
 *----------------------------------------------------------------------*/  
void Intersection::debugIntegrationcells(	map< int, vector <Integrationcell> >&	integrationcellList,
											int id)
{    
    cout << endl;
    cout << "===============================================================" << endl;
    cout << "Debug RESULTING INTEGRATION CELLS " << endl;
    cout << "===============================================================" << endl;
    cout << endl;
    
    map<int,vector <Integrationcell> >::iterator iter;
  	for( iter = integrationcellList.begin(); iter != integrationcellList.end(); ++iter ) 
    {
		cout << "XFEM ELEMENT " << iter->first << " :" << endl;	
		cout << endl;
		
		if(iter->first == id)
		{
			cout << "NUMBER OF INTEGRATIONCELLS : " << iter->second.size() << endl;	
			cout << endl;
			for(unsigned int j = 0; j < iter->second.size(); j++ )
			{
				cout << "IC " << j << ":  " << endl;
				cout << endl;	
				for(unsigned int k = 0; k < iter->second[j].GetCoord().size(); k++)
				{
					//cout << k << "\t";
					for(unsigned int m = 0; m < iter->second[j].GetCoord()[k].size(); m++)
					{
						cout << iter->second[j].GetCoord()[k][m] << "\t";	
					}
					cout << endl;  
				}
				cout << endl;	
			}
			cout << endl;	cout << endl;
		}	
	}		
       
    cout << endl;
    cout << endl;
    cout << endl;
    cout << "===============================================================" << endl;
    cout << "Debug RESULTING INTEGRATION CELLS" << endl;
    cout << "===============================================================" << endl;
    cout << endl; cout << endl; cout << endl;
}



/*----------------------------------------------------------------------*
 |  DB:     Debug only                                       u.may 09/07|
 *----------------------------------------------------------------------*/  
void Intersection::printTetViewOutput(
    int             index,
    tetgenio&       out)
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
void Intersection::printTetViewOutputPLC(
    DRT::Element*   xfemElement,
    int             index,
    tetgenio&       in)
{    
    
    FILE *outFile;
    char filename[100];
    Epetra_SerialDenseVector xsi(3);

    sprintf(filename, "tetgenPLC%d.node", index);
   
    outFile = fopen(filename, "w");
    fprintf(outFile, "%d  %d  %d  %d\n", in.numberofpoints, in.mesh_dim,
          in.numberofpointattributes, in.pointmarkerlist != NULL ? 1 : 0);
    for (int i = 0; i < in.numberofpoints; i++) 
    {
        
        for(int j = 0; j < 3; j++)
            xsi[j] = in.pointlist[i*3 + j];
        
        referenceToCurrentCoordinates(xfemElement, xsi);
        
        fprintf(outFile, "%d  %.16g  %.16g  %.16g", i, xsi[0], xsi[1], xsi[2]);
    
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




/*----------------------------------------------------------------------*
 |  computes the coordinates of a point within a region for  u.may 07/07|
 |  the Tetgen data structure       (only for visualization)            |
 *----------------------------------------------------------------------*/  
// dependency on Fluid3Surface. This cannot be used for the standard fluid tests
// so switch on, if needed, but don't commit
#if 0
void Intersection::computeRegionCoordinates(    
    DRT::Element*  xfemElement,
    DRT::Element*  cutterElement,
    double* regionCoordinates)
{

    bool inStored = false;
    bool outStored = false;
    bool nodeWithin = false;
    int fill = 0;
    int numCornerPoints = 8;
    DRT::Element*  cutterVolume;
    Epetra_SerialDenseVector xsi(3);
    
    
    const DRT::Element::ElementType etype = cutterElement->Type();
  
    switch(etype)
    {
        case DRT::Element::element_fluid3surface:
        {       
            DRT::Elements::Fluid3Surface* fluid3surface = (DRT::Elements::Fluid3Surface *) (cutterElement);
            
            if(fluid3surface == NULL) 
                dserror("fluid3surface is NULL");
                
            cutterVolume = fluid3surface->GetParent();
            
            if(cutterVolume == NULL) 
                dserror("cuttervolume is NULL");
            break;
        }
        default:
            dserror("elementtype not yet handled");
    }

    for(int i = 0; i < numCornerPoints; i++)
    {
        Epetra_SerialDenseVector x;       
        for(int ii = 0; ii < 3; ii++)
            x[ii] = xfemElement->Nodes()[i]->X()[ii];
            
        nodeWithin = checkNodeWithinElement(cutterVolume, x, xsi);
                                
        if(nodeWithin && !inStored)
        {
            for(int j = 0; j < 3; j++)
                regionCoordinates[fill++] = xfemElement->Nodes()[i]->X()[j];
            
            inStored = true;
        }
        
        if(!nodeWithin && !outStored)
        {
            for(int j = 0; j < 3; j++)
                regionCoordinates[fill++] = xfemElement->Nodes()[i]->X()[j];
    
            outStored = true;
        }
        
        if(inStored && outStored)
        {
            break;
        }
    }
}
#endif




#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET

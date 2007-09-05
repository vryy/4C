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

#ifdef D_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "intersection.H"
#include "../drt_xfem/integrationcell.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_element.H"


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
    int numInternalPoints =                 0;
    int numSurfacePoints =                  0;
    vector< DRT::Condition * >              xfemConditions;
    vector< InterfacePoint >                interfacePoints;
    vector <vector< InterfacePoint > >      interfacePointCollection;
    map< int, RefCountPtr<DRT::Element > >  geometryMap;
    DRT::Element*                           cutterElement; 
    DRT::Element*                           xfemElement; 
    Epetra_SerialDenseMatrix                cutterXAABB;
    Epetra_SerialDenseMatrix                xfemXAABB;
   
    higherorder_ = false;
    // put in initialize routine
    initialize(xfemdis, cutterdis);
    vector< InterfacePoint >                pointList;
    vector< vector<int> >                   segmentList(numXFEMSurfaces_);                 
    vector< vector<int> >                   surfacePointList(numXFEMSurfaces_);                
    vector< vector<int> >                   triangleList;
    
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
    for(int k=0; k<xfemdis->NumMyRowElements()-1; k++)
    {
        xfemElement = xfemdis->gElement(k);
        xfemXAABB = computeFastXAABB(xfemElement);
        xfemIntersection = false;
        
        pointList.clear();
        startPointList(pointList);
        
        int countE = 0; 
        //xfemConditions.size()
        for(unsigned int i=0; i<xfemConditions.size(); i++)
        {
            geometryMap = xfemConditions[i]->Geometry();
            if(geometryMap.size()==0)   dserror("geometry does not obtain elements");
            // printf("size of %d.geometry map = %d\n",i, geometryMap.size());
            
            //geometryMap.size()
           
            for(unsigned int j=0; j<geometryMap.size(); j++)
            { 
                cutterElement = geometryMap.find(j)->second.get();
                if(cutterElement == NULL) dserror("geometry does not obtain elements");
            
                //printf("cutterElementcount = %d\n", countE++);
            
                cutterXAABB = computeFastXAABB(cutterElement);                                 
                intersected = intersectionOfXAABB(cutterXAABB, xfemXAABB);    
                //debugXAABBIntersection( cutterXAABB, xfemXAABB, cutterElement, xfemElement, i, k);
                                     
                if(intersected)
                {
                    // collect internal points
                    numInternalPoints= 0; 
                    numSurfacePoints = 0;
                    for(int m=0; m<cutterElement->NumLine() ; m++)                    
                        collectInternalPoints( xfemElement, cutterElement, cutterElement->Nodes()[m],
                                                interfacePoints, k, m, numInternalPoints, numSurfacePoints);
                    
                    /*printf("MAIN 1  PRINT surfaceslist\n");
                    printf("\n");   
                    for(unsigned int i = 0; i < intersectedSurfaces_.size(); i++)
                    {          
                        printf("count = %d\n", i);
                        intersectedSurfaces_[i]->Print(cout);
                        printf("\n");       
                        printf("count = %d\n", i);
                    }
                    printf("\n");   
                    */
                      
                    // collect intersection points                                   
                    for(int m=0; m<xfemElement->NumLine() ; m++) 
                        collectIntersectionPoints(  cutterElement, xfemElement->Lines()[m],
                                                    interfacePoints, numInternalPoints, numSurfacePoints, 
                                                    0, m, false, xfemIntersection);                                         
                     
                    
                    for(int m=0; m<cutterElement->NumLine() ; m++)                                              
                        for(int p=0; p<xfemElement->NumSurface() ; p++) 
                            collectIntersectionPoints(  xfemElement->Surfaces()[p], cutterElement->Lines()[m],
                                                        interfacePoints, numInternalPoints, numSurfacePoints, 
                                                        p, m, true, xfemIntersection);  
                       
                                
                    if(interfacePoints.size()!=0)
                        computeConvexHull(  xfemElement, cutterElement, interfacePoints,  
                                            pointList, surfacePointList, segmentList, triangleList,
                                            numInternalPoints, numSurfacePoints);    
                        
                    interfacePoints.clear();     
                }// if intersected
            }// for-loop over all geometryMap.size()
        }// for-loop over all xfemConditions.size() 
        
        if(xfemIntersection)
        {                                                                          
            //debugTetgenDataStructure(xfemElement, pointList, surfacePointList, segmentList, triangleList);
            computeCDT(xfemElement, cutterElement, pointList, surfacePointList, segmentList, triangleList, integrationcellList);
        }
        
        interfacePointCollection.clear();  
        pointList.clear();
        for(int i=0; i<numXFEMSurfaces_; i++)
        {
            segmentList[i].clear();
            surfacePointList[i].clear();
        }
        triangleList.clear(); 
    }// for-loop over all  actdis->NumMyRowElements()
    
    
    //debugIntegrationcells(integrationcellList,3);
    printf("\n");
    printf("Intersection computed sucessfully\n");
    printf("\n");
}



/*----------------------------------------------------------------------*
 |  INIT:   initializes the private members of the           u.may 08/07|
 |          interface class                                             |
 *----------------------------------------------------------------------*/
void Intersection::initialize(  
    const RefCountPtr<DRT::Discretization> xfemdis,
    const RefCountPtr<DRT::Discretization> cutterdis)
{

    DRT::Element* xfemElement;
    DRT::Element::DiscretizationType xfemDistype;
    
    xfemElement = xfemdis->gElement(0);
    xfemDistype = xfemElement->Shape();
    
    numXFEMSurfaces_ = xfemElement->NumSurface(); 
    numXFEMLines_ = xfemElement->NumLine(); 
   
    
    getNumberOfElementCornerNodes(numXFEMCornerNodes_, xfemDistype);
    getEleNodeNumberingSurfaces(eleNodeNumberingSurfaces_, xfemDistype);
    getEleNodeNumbering_lines_surfaces(eleNodeNumberingLinesSurfaces_, xfemDistype);
    getEleNodeNumbering_nodes_surfaces(eleNodeNumberingNodesSurfaces_, xfemDistype);
    getEleNodeNumbering_nodes_reference(eleNodeNumberingNodesReference_, xfemDistype);
   
}



/*----------------------------------------------------------------------*
 |  ML:     adds two Epetra_SerialDenseVector                u.may 06/07|
 *----------------------------------------------------------------------*/
Epetra_SerialDenseVector Intersection::addTwoVectors(   
    Epetra_SerialDenseVector&   v1,
    Epetra_SerialDenseVector&   v2)
{   
    if(v1.Length() != v2.Length())
        dserror("both vectors need to have the same size\n"); 

    for(int i = 0; i < v1.Length(); i++)
        v1[i] = v1[i] + v2[i];
 
    return v1;
}
	
    
   
/*----------------------------------------------------------------------*
 |  ML:     adds two vector<double>                          u.may 06/07|
 *----------------------------------------------------------------------*/
vector<double> Intersection::addTwoVectors(
    vector<double>&   v1,
    vector<double>&   v2)
{   
    if(v1.size() != v2.size())
        dserror("both vectors need to have the same size\n"); 

    for(unsigned int i = 0; i < v1.size(); i++)
        v1[i] = v1[i] + v2[i];
 
    return v1;
}
    
    

/*----------------------------------------------------------------------*
 |  ML:     subtracts one Epetra_SerialDenseVector from   u.may 06/07   |
 |          another Epetra_SerialDenseVector.                           |
 |          The result is stored in v1                                  |
 *----------------------------------------------------------------------*/
Epetra_SerialDenseVector Intersection::subtractsTwoVectors( 
    Epetra_SerialDenseVector& v1,
    Epetra_SerialDenseVector& v2)
{   
    if(v1.Length() != v2.Length())
        dserror("both vectors need to have the same size\n"); 

    for(int i = 0; i < v1.Length(); i++)
        v1[i] = v1[i] - v2[i];
 
    return v1;
}

    

/*----------------------------------------------------------------------*
 |  ML :    subtracts one vector<double> from another        u.may 06/07|
 |          vector<double> . The result is stored in v1.                |
 *----------------------------------------------------------------------*/
vector<double> Intersection::subtractsTwoVectors(   
    vector <double>& v1,
    vector <double>& v2)
{   
    if(v1.size() != v2.size())
        dserror("both vectors need to have the same size\n"); 

    for(unsigned int i = 0; i < v1.size(); i++)
        v1[i] = v1[i] - v2[i];
 
    return v1;
}



/*----------------------------------------------------------------------*
 |  ML:     computes the cross product                       u.may 08/07|
 |          of 2 Epetra_SerialDenseVector c = a x b                     |
 *----------------------------------------------------------------------*/  
Epetra_SerialDenseVector Intersection::computeCrossProduct(
    const Epetra_SerialDenseVector& a,
    const Epetra_SerialDenseVector& b)
{
    Epetra_SerialDenseVector c(3);
   
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    
    return c;
}



/*----------------------------------------------------------------------*
 |  ML:     computes the absolute value (L2-norm)            u.may 08/07|
 |          of an Epetra_SerialDenseVector                              |
 *----------------------------------------------------------------------*/  
double Intersection::computeAbsoluteValue(   
    Epetra_SerialDenseVector& a)
{
    double absoluteValue = 0.0;
   
    for(int i = 0; i < a.Length(); i++)
         absoluteValue += a[i]*a[i];
        
    absoluteValue = sqrt(absoluteValue);
    
    return absoluteValue;
}



/*----------------------------------------------------------------------*
 |  ML:     multiplies a scalar with a vector                u.may 08/07|
 *----------------------------------------------------------------------*/  
Epetra_SerialDenseVector Intersection::computeScalarVectorMultiplication(   
    const double                s, 
    Epetra_SerialDenseVector&   v)
{
   
    for(int i = 0; i < v.Length(); i++)
        v[i] = s*v[i];
    
    return v;
}



/*----------------------------------------------------------------------*
 |  ML:     computes a Gaussian Elimination for a linear     u.may 06/07|
 |          system of equations                                         |
 *----------------------------------------------------------------------*/
bool Intersection::gaussElimination(    
    Epetra_SerialDenseMatrix&   A,
    Epetra_SerialDenseVector&   b,
    Epetra_SerialDenseVector&   x,
    const bool                  do_piv,
    const int                   dim,    
    const int                   order)
{
    bool solution = true;
    int pivot;
    double tmp[4];    // check tmp and size of A with dim
    double det;
   

    //printf("A[0][] = %8e %8e %8e | %8e\n", A(0,0),A(0,1),A(0,2),b[0]);
    //printf("A[1][] = %8e %8e %8e | %8e\n", A(1,0),A(1,1),A(1,2),b[1]);
    //printf("A[2][] = %8e %8e %8e | %8e\n", A(2,0),A(2,1),A(2,2),b[2]);

    if (!do_piv) 
    {
        for (int k=0;k<dim;k++)
        {            
            A(k,k)=1./A(k,k);

            for (int i=k+1;i<dim;i++)
            {
                A(i,k) = A(i,k) * A(k,k);
                x[i] = A(i,k);

                for (int j=k+1;j<dim;j++)
                    A(i,j) = A(i,j) - A(i,k) * A(k,j);               
            }

            for (int i=k+1;i<dim;i++)
                b[i]=b[i]-x[i]*b[k];           
        }
    }
    else 
    {
        for (int k=0;k<dim;k++)
        {
            pivot = k;
            // search for pivot element 
            for (int i=k+1;i<dim;i++)
                pivot = (fabs(A(pivot,pivot)) < fabs(A(i,k))) ? i : pivot;
            
            // copy pivot row to current row 
            if (pivot != k) 
            {
                for (int j=0;j<dim;j++)
                    tmp[j] = A(pivot,j);

                tmp[dim] = b[pivot];

                for (int j=0;j<dim;j++)
                    A(pivot,j) = A(k,j);

                b[pivot] = b[k];

                for (int j=0;j<dim;j++)
                    A(k,j) = tmp[j];

                b[k] = tmp[dim];
            }

            A(k,k) = 1./A(k,k);
            //printf("inf_diag = %8e\n", A(k,k));
            //fflush(NULL);

            for (int i=k+1;i<dim;i++)
            {
                A(i,k) = A(i,k) * A(k,k);
                x[i] = A(i,k);

                for (int j=k+1;j<dim;j++)
                    A(i,j) = A(i,j) - A(i,k) * A(k,j);               
            }

            for (int i=k+1;i<dim;i++)
                b[i]=b[i]-x[i]*b[k];
         
            //printf("A[0][] = %8e %8e %8e | %8e\n", A(0,0),A(0,1),A(0,2),b[0]);
            //printf("A[1][] = %8e %8e %8e | %8e\n", A(1,0),A(1,1),A(1,2),b[1]);
            //printf("A[2][] = %8e %8e %8e | %8e\n", A(2,0),A(2,1),A(2,2),b[2]);
        }
    }

    det = (1.0/A(0,0))*(1.0/A(1,1))*(1.0/A(2,2));
    if(fabs(det) < TOL7_ && order == 1)
    {
        solution = false;
        //printf("matrix is singular A1 = %f, A2 = %f, A3 = %f, det = %f\n ", 1/A(0,0), 1/A(1,1), 1/A(2,2),det );
    }
    // backward substitution 
    x[dim-1]=b[dim-1]*A(dim-1,dim-1);

    for (int i=dim-2;i>=0;i--)
    {
        for (int j=dim-1;j>i;j--)
            b[i]=b[i]-A(i,j)*x[j];
        
        x[i]=b[i]*A(i,i);
    }

    //for (i=0;i<dim;i++)
    //  printf("%8e ",x[i]);
    //printf("\n");
    return solution;
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
        
    for(int dim=0; dim<3; dim++)
    {
        xsi(dim) = 0.0;
        for(int i=0; i<numNodes; i++)
            xsi(dim) += element->Nodes()[i]->X()[dim] * funct(i);
    }
}



/*----------------------------------------------------------------------*
 | GM:  transforms a node in current coordinates            u.may 07/07 |
 |      into reference coordinates                                      |
 *----------------------------------------------------------------------*/  
void Intersection::currentToReferenceCoordinates(   
    DRT::Element* element, 
    Epetra_SerialDenseVector& xsi)
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
    std::vector<double>&         node,
    Epetra_SerialDenseMatrix&    XAABB)
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
    Epetra_SerialDenseMatrix& cutterXAABB, 
    Epetra_SerialDenseMatrix& xfemXAABB)
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
    const int                       nodeId,
    int&                            numInternalPoints,
    int&                            numSurfacePoints)
{
    Epetra_SerialDenseVector xsi(3);
    Epetra_SerialDenseVector x(3);
       
    x[0] = node->X()[0];
    x[1] = node->X()[1];
    x[2] = node->X()[2];
    
    bool nodeWithinElement = checkNodeWithinElement(element, x, xsi);
    //debugNodeWithinElement(element,node,xsi,elemId ,nodeId, nodeWithinElement);  
    
    if(nodeWithinElement)
    {   
        numInternalPoints++; 
        InterfacePoint ip;
        //debugNodeWithinElement(element,node,xsi,elemId ,nodeId, nodeWithinElement);  
          
        // check if set neighbouring surfaces
        checkIfOnSurfaceAndNode(element, xsi, ip, numSurfacePoints);
        
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
    Epetra_SerialDenseMatrix&   A,
    Epetra_SerialDenseVector&   xsi,
    DRT::Element*               element)                                                  
{	
    const int numNodes = element->NumNode();
    vector<int> actParams(1,0);
    Epetra_SerialDenseMatrix deriv1(3, numNodes);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
    
    if(!checkIfVolumeElement(element))
   		dserror("element has to be a volume element\n");
    
    params.set("action","calc_ShapeDeriv1");
    actParams[0] = numNodes;  
       
    element->Evaluate(params, dummyDis, actParams, deriv1, emptyM , emptyV, xsi, emptyV);       
    
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            A(i,j) = 0.0;
            
    for(int dim=0; dim<3; dim++)
        for(int j=0; j<numNodes; j++)
        {
            A(dim, 0) += element->Nodes()[j]->X()[dim] * deriv1(0,j);
            A(dim, 1) += element->Nodes()[j]->X()[dim] * deriv1(1,j);
            A(dim, 2) += element->Nodes()[j]->X()[dim] * deriv1(2,j);
        }
}



/*----------------------------------------------------------------------*
 |  CLI:    updates the rhs for the computation if a         u.may 06/07|
 |          node is in a given element                                  |
 *----------------------------------------------------------------------*/
void Intersection::updateRHSForNWE( 
    Epetra_SerialDenseVector&   b,
    Epetra_SerialDenseVector&   xsi,
    Epetra_SerialDenseVector&   x,
    DRT::Element*               element)                                                  
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
    
    if(!checkIfVolumeElement(element))
   		dserror("element has to be a volume element\n");
      
    element->Evaluate(params, dummyDis, actParams, emptyM, emptyM , funct, xsi, emptyV);        
    
    for(int dim=0; dim<3; dim++)
        b(dim) = 0.0;
    
    for(int dim=0; dim<3; dim++)
        for(int i=0; i<numNodes; i++)
        {
            b(dim) += (-1) * element->Nodes()[i]->X()[dim] * funct(i);
        }
        
    b = addTwoVectors(b,x);
}




/*----------------------------------------------------------------------*
 |  CLI:    checks if a node is within a given element       u.may 06/07|	
 *----------------------------------------------------------------------*/
bool Intersection::checkNodeWithinElement(	
    DRT::Element* element,
    Epetra_SerialDenseVector& x,
    Epetra_SerialDenseVector& xsi)
{

	bool nodeWithinElement = true;
    int iter = 0;
    const int maxiter = 500;
    double residual = 1.0;
    Epetra_SerialDenseMatrix A(3,3);
    Epetra_SerialDenseVector b(3);
    Epetra_SerialDenseVector dx(3);
  
    xsi[0] = 0.0; xsi[1] = 0.0; xsi[2] = 0.0;
    dx = xsi;
        
    updateRHSForNWE( b, xsi, x, element);
    
    while(residual > TOL14_)
    { 
        updateAForNWE( A, xsi, element);
        
        if(!gaussElimination(A, b, dx, true, 3, 1))
        {
            nodeWithinElement = false;
            break;
        }   
        xsi = addTwoVectors(xsi,dx);
        
        if(iter >= maxiter)
        {   
            nodeWithinElement = false;
            break;
        }       
        updateRHSForNWE(b, xsi, x, element);
        residual = b.Norm2();
        iter++;
    }
    
    if( (fabs(xsi[0])-1.0) > TOL7_  || (fabs(xsi[1])-1.0) > TOL7_  || (fabs(xsi[2])-1.0) > TOL7_ )   
        nodeWithinElement = false;
    
    return nodeWithinElement;
}



/*----------------------------------------------------------------------*
 |  CLI:    checks if a node that lies within an element     u.may 06/07|
 |          lies on one of its surfaces or nodes                        |                                           
 *----------------------------------------------------------------------*/
bool Intersection::checkIfOnSurfaceAndNode( 
    DRT::Element*                   element,
    Epetra_SerialDenseVector&       xsi,    
    InterfacePoint&                 ip,
    int&                            numSurfacePoints)
{
    bool onSurface = false;
    int count = 0;
      
    count = getSurfaces(xsi, ip.surfaces, element->Shape());
        
    if(count > 0)
    {
        onSurface = true;
        ip.nsurf = count;
        ip.pType = surfaceP;
        numSurfacePoints++;
    }
    else
    {
        onSurface = false;
        ip.nsurf = 0;
        ip.pType = internalP;
    }
    
    return onSurface;
}



/*----------------------------------------------------------------------*
 |  CLI:    collects all intersection points of a line and   u.may 06/07|
 |          and a surface                                               |
 *----------------------------------------------------------------------*/  
void Intersection::collectIntersectionPoints(   
    DRT::Element*                   surfaceElement,
    DRT::Element*                   lineElement,
    std::vector<InterfacePoint>&    interfacePoints, 
    const int                       numInternalPoints,
    const int                       numSurfacePoints,
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
    }
    
    intersected = computeCurveSurfaceIntersection(surfaceElement, lineElement, xsi, upLimit, loLimit);
                                        
    if(intersected) 
        numInterfacePoints = addIntersectionPoint(  surfaceElement, lineElement,xsi, upLimit, loLimit, 
                                                    interfacePoints, surfaceId, lineId, lines);
      
    
    // in this case a node of this line lies on the facet of the xfem element
    // but there is no intersection within the element                                          
    if(!((int) interfacePoints.size() == numSurfacePoints)) 
        xfemIntersection = true;
      
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
    dx = xsi;
    
    updateRHSForCSI( b, xsi, surfaceElement, lineElement);
                                
    while(residual > TOL14_)
    {   
        updateAForCSI( A, xsi, surfaceElement, lineElement);
        if(!gaussElimination(A, b, dx, true, 3, 1))
        {
            intersection = false;
            break;  
        } 
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
            A(i,j) = 0.0;
	
	params.set("action","calc_ShapeDeriv1");
	surfaceElement->Evaluate(params, dummyDis, actParams, surfaceDeriv1, emptyM , emptyV, xsiSurface, emptyV);			
    //printf("num of nodes surfaces = %d\n", numOfNodesSurface);
       	
    for(int dim=0; dim<3; dim++)
   	    for(int i=0; i<numNodesSurface; i++)
		{
			A(dim,0) += surfaceElement->Nodes()[i]->X()[dim] * surfaceDeriv1(0,i);
			A(dim,1) += surfaceElement->Nodes()[i]->X()[dim] * surfaceDeriv1(1,i);
		}
		
   
    actParams[0] = numNodesLine;
    lineElement->Evaluate(params, dummyDis, actParams, lineDeriv1, emptyM , emptyV, xsiLine, emptyV);

    for(int dim=0; dim<3; dim++)
        for(int i=0; i<numNodesLine; i++)
        {
            A(dim,2) +=  (-1) * lineElement->Nodes()[i]->X()[dim] * lineDeriv1(0,i);
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
			b(dim) += (-1) * surfaceElement->Nodes()[i]->X()[dim] * surfaceFunct(i);
		}
	
	for(int dim=0; dim<3; dim++)
   		for(int i=0; i<numNodesLine; i++)
		{
			b(dim) += lineElement->Nodes()[i]->X()[dim] * lineFunct(i);
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
        ip.surfaces[0] = eleNodeNumberingLinesSurfaces_[lineId][0];
        ip.surfaces[1] = eleNodeNumberingLinesSurfaces_[lineId][1];
        ip.coord[0] = xsi[0]; 
        ip.coord[1] = xsi[1]; 
    }
    
    ip.coord[2] = 0.0; 
    ip.pType = intersectionP;
    
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
void Intersection::computeConvexHull(   
    DRT::Element*           xfemElement,
    DRT::Element*           surfaceElement,
    vector<InterfacePoint>& interfacePoints,
    vector<InterfacePoint>& pointList,
    vector< vector<int> >&  surfacePointList,
    vector< vector<int> >&  segmentList,
    vector< vector<int> >&  triangleList,
    const int               numInternalPoints,
    const int               numSurfacePoints )
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
                    coordinates[fill++] = interfacePoints[i].coord[j]; 

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
    
        storePoint(vertices[0], interfacePoints, positions, pointList );
        vertices.erase(vertices.begin());
        
        if(interfacePoints.size() > 1)
        {
            // store points, segments and triangles for the computation of the
            // Constrained Delaunay Tetrahedralization with Tetgen  
            searchPoint = vertices[0];   
            storePoint(vertices[0], interfacePoints, positions, pointList );
            vertices.erase(vertices.begin());
        }
        
        int countWhile = 0;
        while(vertices.size()>2)
        {                    
            findNextSegment(vertices, searchPoint);
            storePoint(searchPoint, interfacePoints, positions, pointList);
            countWhile++;           
        } 
        vertices.clear();
       

        // only one point on a plane
        if(numSurfacePoints  == 1)
            storeSurfacePoints(pointList, interfacePoints, surfacePointList);
            
            
        // cutter element lies on the surface of an xfem element
        if(numInternalPoints == numSurfacePoints && numInternalPoints != 0)
        {
            if(numSurfacePoints > 1)              
                storeSegments(pointList, positions, segmentList);           
           
        }
        else
        {
            if(interfacePoints.size() > 1)
                storeSegments(pointList, positions, segmentList);
           
            if(interfacePoints.size() > 2)
            {
                pointList.push_back(midpoint);
                storeTriangles(pointList, positions, triangleList);
            }
        }
        interfacePoints.clear();
    }
}



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
    vector< InterfacePoint >&   			pointList,
    vector< vector<int> >&                  surfacePointList,
    vector< vector<int> >&      			segmentList,
    vector< vector<int> >&      			triangleList,
    map< int, vector <Integrationcell> >&	integrationcellList)
{
    int dim = 3;
    int nsegments = 0; 
    int nsurfPoints = 0;
    tetgenio in;
    tetgenio out;
    char switches[] = "pnno2QY";    
    tetgenio::facet *f;
    tetgenio::polygon *p;
    double regionCoordinates[6];
    

    // allocate pointlist
    in.numberofpoints = pointList.size();
    in.pointlist = new REAL[in.numberofpoints * dim];
       
    // fill point list
    int fill = 0;
    for(int i = 0; i <  in.numberofpoints; i++)
        for(int j = 0; j < dim; j++)  
            in.pointlist[fill++] = (REAL) pointList[i].coord[j]; 
 
    // allocate facetlist
    if(triangleList.size()>0)       in.numberoffacets = numXFEMSurfaces_ + triangleList.size(); 
    else                            in.numberoffacets = numXFEMSurfaces_;   
      
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    
    
    // loop over all xfem element surfaces
    for(int i = 0; i < numXFEMSurfaces_; i++)
    {
        f = &in.facetlist[i];
        if(segmentList[i].size() > 0)           nsegments = (int) (segmentList[i].size()/2);
        else                                    nsegments = 0;
        if(surfacePointList[i].size() > 0)      nsurfPoints = surfacePointList[i].size();
        else                                    nsurfPoints = 0;
        f->numberofpolygons = 1 + nsegments + nsurfPoints; 
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];
        p->numberofvertices = 4;
        p->vertexlist = new int[p->numberofvertices];
        for(int j = 0; j < 4; j ++)
            p->vertexlist[j] = eleNodeNumberingSurfaces_[i][j];
           
      
        int count = 0;
        for(int j = 1; j < 1 + nsegments; j ++)
        {
            if(segmentList[i].size() > 0)
            {             
                p = &f->polygonlist[j];
                p->numberofvertices = 2;
                p->vertexlist = new int[p->numberofvertices];
            
                for(int k = 0; k < 2; k ++)
                   p->vertexlist[k] = segmentList[i][count++];   
            }
        } 
        
        count = 0;
        for(int j = 1 + nsegments; j < f->numberofpolygons; j++)
        {
            if(surfacePointList[i].size() > 0)
            {             
                p = &f->polygonlist[j];
                p->numberofvertices = 1;
                p->vertexlist = new int[p->numberofvertices];
            
                p->vertexlist[0] = surfacePointList[i][count++];   
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
            p->vertexlist[j] = triangleList[i - element->NumSurface()][j];    
    }
  
        
    // set facetmarkers external boundary facets are marked with -1 
    // and internal boundary facets with -2
    for(int i = 0; i < numXFEMSurfaces_; i ++)
        in.facetmarkerlist[i] = -1;   
    
    for(int i = numXFEMSurfaces_; i < in.numberoffacets; i ++)
        in.facetmarkerlist[i] = -2;   
        
        
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
 
  
    // recover curved interface for higher order meshes
    bool curvedInterface = true;
    if(curvedInterface)
        recoverCurvedInterface(element, out);
 
 
    // store integrationcells
    vector<double> tetnodes(3);
    vector< vector<double> > tetrahedronCoord;
    vector< Integrationcell > listperElement;
    
    for(int i=0; i<out.numberoftetrahedra; i++ )
    {   
        for(int j = 0; j < out.numberofcorners; j++)
        {
            for(int dim = 0; dim < 3; dim++)
                tetnodes[dim] = out.pointlist[out.tetrahedronlist[i*out.numberofcorners+j]*dim+dim];
         
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
    vector< InterfacePoint >&  pointList)
{
    InterfacePoint ip;
        
    for(int i = 0; i < numXFEMCornerNodes_; i++)
    {
        ip.nsurf = 3;
        
        // change for other element types
        for(int j = 0; j < 3; j++) 
        { 
           ip.coord[j] = eleNodeNumberingNodesReference_[i][j]; 
           ip.surfaces[j] = eleNodeNumberingNodesSurfaces_[i][j]; 
        }
        pointList.push_back(ip);               
    }
}



/*----------------------------------------------------------------------*
 |  CDT:    stores a point within a list of points           u.may 06/07|
 |          which is to be copy to the tetgen data structure            |
 |          for the computation of the                                  |
 |          Constrained Delauney Triangulation                          |
 *----------------------------------------------------------------------*/  
void Intersection::storePoint(  
    vector<double>&             point, 
    vector<InterfacePoint>&     interfacePoints, 
    vector<int>&                positions, 
    vector<InterfacePoint>&     pointList )
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
            for(it1 = pointList.begin(); it1 != pointList.end(); it1++ )  
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
                pointList.push_back(interfacePoints[i]);
                positions.push_back(pointList.size()-1);
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
    vector<InterfacePoint>& interfacePoints)
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
    vector<InterfacePoint>&     pointList, 
    vector<InterfacePoint>&     interfacePoints,
    vector< vector<int> >&      surfacePointList)
{   
    int count = -1;
    vector<InterfacePoint>::iterator it1;
    
    for(unsigned int i = 0; i < interfacePoints.size(); i++)
    {
        if(interfacePoints[i].pType == surfaceP)
        {  
            count = -1;
            for(it1 = pointList.begin(); it1 != pointList.end(); it1++ )  
            {
                count++;
                if(comparePoints(interfacePoints[i].coord, it1->coord, 3)) 
                {
                    surfacePointList[interfacePoints[i].surfaces[0]].push_back(count);          
                    break;
                }
            }
            break;
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
    vector<InterfacePoint>&     pointList, 
    vector<int>&                positions, 
    vector< vector<int> >&      segmentList)
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
        
        for(int j = 0; j < pointList[pos1].nsurf; j++ )  
            for(int k = 0; k < pointList[pos2].nsurf; k++ ) 
            {
                surf1 = pointList[pos1].surfaces[j];
                surf2 = pointList[pos2].surfaces[k];
             
                if( (surf1 == surf2) &&  (pointList[pos1].nsurf < 2 || pointList[pos2].nsurf < 2) ) 
                { 
                    alreadyInList = false;
                    
                    for(unsigned int is = 0 ; is < segmentList[surf1].size() ; is = is + 2)
                    {
                        if( (segmentList[surf1][is] == pos1  &&  segmentList[surf1][is+1] == pos2)  ||
                            (segmentList[surf1][is] == pos2  &&  segmentList[surf1][is+1] == pos1) )
                        {
                            alreadyInList = true;
                            break;
                        }
                    }
                    
                    if(!alreadyInList)
                    {
                        segmentList[surf1].push_back(pos1);
                        segmentList[surf1].push_back(pos2);
                    }
                }
            }
    }
      
    pos1 = positions[positions.size()-1];
    pos2 = positions[0]; 
    
    for(int j = 0; j < pointList[pos1].nsurf; j++ )  
        for(int k = 0; k < pointList[pos2].nsurf; k++ ) 
        {
            surf1 = pointList[pos1].surfaces[j];
            surf2 = pointList[pos2].surfaces[k];

            if((surf1 == surf2) && (pointList[pos1].nsurf < 2 || pointList[pos2].nsurf < 2) ) 
            { 
                alreadyInList = false;
                for(unsigned int is = 0 ; is < segmentList[surf1].size(); is = is + 2)
                {
                    if( (segmentList[surf1][is] == pos1  &&  segmentList[surf1][is+1] == pos2)  ||
                        (segmentList[surf1][is] == pos2  &&  segmentList[surf1][is+1] == pos1) )
                    {
                        alreadyInList = true;
                        break;
                    }
                }
                
                if(!alreadyInList)
                {
                    segmentList[surf1].push_back(pos1);
                    segmentList[surf1].push_back(pos2);
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
    vector<InterfacePoint>&     pointList, 
    vector<int>                 positions, 
    vector< vector<int> >&      triangleList)
{
    vector<int> triangle(3,0);
    
    for(unsigned int i = 0; i < positions.size()-1; i++ )
    {
        triangle[0] = positions[i];
        triangle[1] = positions[i+1];
        triangle[2] = pointList.size()-1;
        
        triangleList.push_back(triangle);
    }
    
    triangle[0] = positions[positions.size()-1];
    triangle[1] = positions[0];
    triangle[2] = pointList.size()-1;
        
    triangleList.push_back(triangle);
}



/*----------------------------------------------------------------------*
 |  RQI:    recovers the curved interface after the          u.may 08/07|
 |          Constrained Delaunay Tetrahedralization                     |
 *----------------------------------------------------------------------*/  
void Intersection::recoverCurvedInterface(
    DRT::Element*  element, 
    tetgenio& out
    )
{
    bool intersected  = false;
    int index1, index2;
    int higherOrderIndex;
    vector<int> order(3,0);
    vector<int> tetraCornerIndices(4,0); 
    Epetra_SerialDenseVector node2(3);
    Epetra_SerialDenseVector node1(3);
    Epetra_SerialDenseVector xsi(3);
    vector < vector <Epetra_SerialDenseVector> > surfaceNodes;
    
    // list of point markers , if already visited = 1 , if not = 0
    int* visitedPointIndexList = new int[out.numberofpoints];      
    for(int i = 0; i<out.numberofpoints; i++)
        visitedPointIndexList[i] = 0;
        
    for(int i=0; i<out.numberoftrifaces; i++)
    {
        //printf("tri face marker = %d\n", out.trifacemarkerlist[i]);     
        // run over all faces not lying in on of the xfem element planes
        if(out.trifacemarkerlist[i] == -2)
        {           
            getTetrahedraInformation(out, tetraCornerIndices, order, i);
             
            // run over each triface
            for(int j=0; j<3 ;j++)
            {                   
                index1 = j;
                index2 = j+1;
                if(index2 > 2) index2 = 0;
                
                higherOrderIndex = getHigherOrderTetIndex(order[index1], order[index2]); 
                if(visitedPointIndexList[higherOrderIndex] == 0)
                {                                              
                    computeIntersectionNormal(index1, index2, out, tetraCornerIndices, node1, node2);
                    
                    // transform surface elements into the ref system of the xfem element
                    transformSurfaceNodes(element, surfaceNodes);
                   
                    // compute Newton for intersection
                    for(unsigned int k = 0; k < intersectedSurfaces_.size(); k++)
                    {    
                        intersected = computeRecovery(surfaceNodes[k], intersectedSurfaces_[k], node1, node2, xsi);
                        if(intersected)
                        {
                            storeHigherOrderNode(higherOrderIndex, xsi, intersectedSurfaces_[k], out);
                            visitedPointIndexList[higherOrderIndex] = 1;
                            break;
                        }
                    }
                }
            }
        }
    }
    delete [] visitedPointIndexList;
    intersectedSurfaces_.clear();
}



/*----------------------------------------------------------------------*
 |  RQI:    computes the intersection between a              u.may 08/07|
 |          curve and a surface                    RQI                  |           
 *----------------------------------------------------------------------*/
bool Intersection::computeRecovery( 
    vector <Epetra_SerialDenseVector>&  surfaceNodes,
    DRT::Element*                       surfaceElement,
    Epetra_SerialDenseVector&           node1,
    Epetra_SerialDenseVector&           node2,
    Epetra_SerialDenseVector&           xsi)
{
    int iter = 0;
    int maxiter = 500;
    bool intersection = true;
    double residual = 1.0;
    Epetra_SerialDenseMatrix A(3,3);
    Epetra_SerialDenseVector b(3);
    Epetra_SerialDenseVector dx(3);
    dx = xsi;
    
    updateRHSForRQI( b, xsi, surfaceNodes, surfaceElement, node1, node2);
                                
    while(residual > TOL14_)
    {   
        updateAForRQI( A, xsi, surfaceNodes, surfaceElement);
        if(!gaussElimination(A, b, dx, true, 3, 1))
        {
            intersection = false;
            break;  
        } 
        xsi = addTwoVectors(xsi,dx);
        
        if(iter >= maxiter)
        {   
            intersection = false;
            break;
        }       
       
        updateRHSForRQI( b, xsi, surfaceNodes, surfaceElement, node1, node2);
        residual = b.Norm2(); 
        iter++;
    } 
    
    if( (fabs(xsi[0])-1.0) > TOL7_  || (fabs(xsi[1])-1.0) > TOL7_  || (fabs(xsi[2])-1.0) > TOL7_ )   
        intersection = false;
        
    return intersection;
}



/*----------------------------------------------------------------------*
 |  RQI:    updates the systemmatrix for the                 u.may 06/07|
 |          computation of a curve-surface intersection                 |
 |          for the recovery of the curved interface                    |
 *----------------------------------------------------------------------*/
void Intersection::updateAForRQI(   
    Epetra_SerialDenseMatrix&           A,
    Epetra_SerialDenseVector&           xsi,
    vector <Epetra_SerialDenseVector>&  surfaceNodes,
    DRT::Element*                       surfaceElement)                                                 
{   
    const int numNodesSurface = surfaceNodes.size();
    vector<int> actParams(1,0);
    Epetra_SerialDenseMatrix surfaceDeriv1(2,numNodesSurface);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    Epetra_SerialDenseVector xsiSurface(2);
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
     
    xsiSurface(0) = xsi[0]; // r-coordinate surface
    xsiSurface(1) = xsi[1]; // s-coordinate surface
    actParams[0] = numNodesSurface;     
    
    for(unsigned int i=0; i<3; i++)
        for(unsigned int j=0; j<3; j++)
            A(i,j) = 0.0;
    
    params.set("action","calc_ShapeDeriv1");
    surfaceElement->Evaluate(params, dummyDis, actParams, surfaceDeriv1, emptyM , emptyV, xsiSurface, emptyV);          
    //printf("num of nodes surfaces = %d\n", numOfNodesSurface);
        
    for(int dim=0; dim<3; dim++)
        for(int i=0; i<numNodesSurface; i++)
        {
            A(dim,0) += surfaceNodes[i][dim] * surfaceDeriv1(0,i);
            A(dim,1) += surfaceNodes[i][dim] * surfaceDeriv1(1,i);
        }
        
    //A(dim,2) = 0 for lines
}



/*----------------------------------------------------------------------*
 |  RQI:    updates the right-hand-side for the              u.may 08/07|
 |          computation of a curve-surface intersection                 |     
 |          for the recovery of the curved surface (RQI)                |
 *----------------------------------------------------------------------*/
void Intersection::updateRHSForRQI( 
    Epetra_SerialDenseVector&           b,
    Epetra_SerialDenseVector&           xsi,    
    vector <Epetra_SerialDenseVector>&  surfaceNodes,
    DRT::Element*                       surfaceElement,
    Epetra_SerialDenseVector&           node1,  
    Epetra_SerialDenseVector&           node2)                                                    
{
    int numNodesSurface = surfaceNodes.size();
    vector<int> actParams(1,0);
    Epetra_SerialDenseVector surfaceFunct(numNodesSurface);
    Epetra_SerialDenseMatrix emptyM;
    Epetra_SerialDenseVector emptyV;
    Epetra_SerialDenseVector xsiSurface(2);
    DRT::Discretization dummyDis("dummy discretization", null);
    ParameterList params;
    
    
    xsiSurface(0) = xsi[0]; // r-coordinate surface
    xsiSurface(1) = xsi[1]; // s-coordinate surface
    params.set("action","calc_Shapefunction");
            
    for(unsigned int i=0; i<3; i++)   
        b[i] = 0.0;
    
    actParams[0] = numNodesSurface;            
    surfaceElement->Evaluate(params, dummyDis, actParams, emptyM, emptyM , surfaceFunct, xsiSurface, emptyV);           
  
    for(int dim=0; dim<3; dim++)
        for(int i=0; i<numNodesSurface; i++)
        {
            b(dim) += (-1) * surfaceNodes[i][dim] * surfaceFunct(i);
        }
        
    for(int dim=0; dim<3; dim++)
        b(dim) += node1[dim] * 0.5*(1.0 - xsi[2]) + node2[dim] * 0.5*(1.0 + xsi[2]);
        
}



/*----------------------------------------------------------------------*
 |  RQI:    stores a pointer to each intersecting            u.may 07/07|
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
void Intersection::getTetrahedraInformation(   
    tetgenio& out, 
    vector<int>& tetraCornerIndices,
    vector<int>& order,
    const int index)
{            
      
    int nodeIndex = 0;         
    int tetIndex = out.adjtetlist[index*2];
    
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
    
    // find order    
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
 |  RQI:    returns the index of a higher order              u.may 08/07|
 |          tetrahedron node index lying between                        |
 |          two specified corner node indices                           |
 *----------------------------------------------------------------------*/  
int Intersection::getHigherOrderTetIndex(
    const int index1, 
    const int index2)
{

    int cornerIndex = 0;
    
    if(index1 == 0 && index2 == 1 || index1 == 1 && index2 == 0 )   cornerIndex = 4;
    else if(index1 == 1 && index2 == 2 || index1 == 2 && index2 == 1 )   cornerIndex = 5;
    else if(index1 == 2 && index2 == 0 || index1 == 0 && index2 == 2 )   cornerIndex = 6;
    else if(index1 == 0 && index2 == 3 || index1 == 3 && index2 == 0 )   cornerIndex = 7;
    else if(index1 == 1 && index2 == 3 || index1 == 3 && index2 == 1 )   cornerIndex = 8;
    else if(index1 == 2 && index2 == 3 || index1 == 3 && index2 == 2 )   cornerIndex = 9;
    else dserror("no valid tetrahedron edge found");

    return cornerIndex;
}



/*----------------------------------------------------------------------*
 |  RQI:    computes the normal to the interface edge of     u.may 08/07|
 |          the tetrahedron facet lying within this facet               |
 *----------------------------------------------------------------------*/  
void Intersection::computeIntersectionNormal(    
    const int                   index1,
    const int                   index2,
    tetgenio&                   out,
    vector<int>&                tetraCornerIndices,
    Epetra_SerialDenseVector&   node1,
    Epetra_SerialDenseVector&   node2)
{            
      
    double                    absVal = 0.0;   
    Epetra_SerialDenseVector  p1(3);
    Epetra_SerialDenseVector  p2(3);
    Epetra_SerialDenseVector  p3(3);        
    Epetra_SerialDenseVector  m(3);
    Epetra_SerialDenseVector  n(3);
    Epetra_SerialDenseVector  r(3);
    Epetra_SerialDenseVector  r1(3);
    Epetra_SerialDenseVector  r2(3);
    Epetra_SerialDenseVector  dis(3);
    Epetra_SerialDenseVector  absDirection(3); 
    
    
    for(int i=0; i<3; i++)
    {
        p1[i] = out.pointlist[tetraCornerIndices[3]*3 + i];
        p2[i] = out.pointlist[tetraCornerIndices[index1]*3 + i];
        p3[i] = out.pointlist[tetraCornerIndices[index2]*3 + i];                
    }
    
                         
    // compute direction vectors of the plane 
    r1 = subtractsTwoVectors(p1, p2);
    r2 = subtractsTwoVectors(p3, p2);
    // normal of the plane
    n = computeCrossProduct(r1, r2);
    // direction vector of the intersection line
    r = computeCrossProduct(n, r2);              
    // computes the start point of the line
    m = computeLineMidpoint(p2, p3);


    dis = subtractsTwoVectors(p1, m);
    absVal = computeAbsoluteValue(dis);

    absDirection = computeScalarVectorMultiplication(absVal, r);

    // compute nodes of the normal to the interface edge of the tetrahedron
    node1 = addTwoVectors(m, absDirection);               
    node2 = subtractsTwoVectors(m, absDirection);
}



/*----------------------------------------------------------------------*
 |  RQI:    computes the midpoint of a line                  u.may 08/07|
 *----------------------------------------------------------------------*/  
Epetra_SerialDenseVector Intersection::computeLineMidpoint( 
    Epetra_SerialDenseVector& p1, 
    Epetra_SerialDenseVector& p2)
{
    Epetra_SerialDenseVector midpoint(3);
    
    for(int i=0; i<3; i++)
        midpoint[i] = (p1[i] + p2[i])/2.0;
        
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
}




/*----------------------------------------------------------------------*
 |  RQI:    stores the higher-order node in the pointlist    u.may 08/07|
 |          at the place of the linear node                             |
 *----------------------------------------------------------------------*/  
void Intersection::storeHigherOrderNode( 
    int                         index,
    Epetra_SerialDenseVector&   xsi, 
    DRT::Element*               surfaceElement, 
    tetgenio&                   out )
{
    referenceToCurrentCoordinates(surfaceElement, xsi);
    for(int i = 0; i < 3; i++)
        out.pointlist[index*3+i]   = xsi[i];  
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
    for(int jE = 0; jE < element->NumNode(); jE++)
    {
        element->Nodes()[jE]->Print(cout);
        cout << endl;
    }
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
void Intersection::debugTetgenDataStructure(    DRT::Element*               element,
                                                vector< InterfacePoint >&   pointList,
                                                vector< vector<int> >&      surfacePointList,
                                                vector< vector<int> >&      segmentList,
                                                vector< vector<int> >&      triangleList)
{    
    cout << endl;
    cout << "===============================================================" << endl;
    cout << "Debug Tetgen Data Structure " << endl;
    cout << "===============================================================" << endl;
    cout << endl;
    cout << "POINT LIST " << " :" << endl;
    cout << endl;
    Epetra_SerialDenseVector xsi(3);
    for(unsigned int i = 0; i < pointList.size(); i++)
    {
        for(int j = 0; j< 3; j++)
        {
            xsi[j] = pointList[i].coord[j];
        }
        referenceToCurrentCoordinates(element, xsi);
        
        cout << i << ".th point:   ";
        for(int j = 0; j< 3; j++)
        {
            //cout << setprecision(10) << pointList[i].coord[j] << "\t"; 
             printf("%20.16f\t", pointList[i].coord[j] );
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
    for(unsigned int i = 0; i < segmentList.size(); i++)
    {
        cout << i << ".th segment:   ";
        int count = 0;
        for(unsigned int j = 0; j < segmentList[i].size(); j++)
                cout << segmentList[i][count++] << "\t";
        
        count = 0;
        for(unsigned int j = 0; j < surfacePointList[i].size(); j++)
                cout << surfacePointList[i][count++] << "\t";
                
        cout << endl;
        cout << endl;
    }
    cout << endl;
    cout << endl;
    
    cout << endl;
    cout << "TRIANGLE LIST " << " :" << endl;
    cout << endl;
    for(unsigned int i = 0; i < triangleList.size(); i++)
    {
        cout << i << ".th triangle:   ";
        for(int j = 0; j< 3; j++)
        {
            cout << triangleList[i][j] << "\t";
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
					cout << k << "\t";
					for(unsigned int m = 0; m < iter->second[j].GetCoord()[k].size(); m++)
					{
						cout << iter->second[j].GetCoord()[k][m] << "\t";	
					}
					cout << "\t";	
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
 |  computes the coordinates of a point within a region for  u.may 07/07|
 |  the Tetgen data structure       (only for visualization)            |
 *----------------------------------------------------------------------*/  
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




#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_XFEM



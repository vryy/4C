/*!----------------------------------------------------------------------
\file intersection.cpp

\brief collection of intersection tool 

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
#include "../drt_lib/drt_element.H"
#include "../drt_f3/fluid3.H"

#include "../headers/definitions.h"     //remove

using namespace XFEM;


/*----------------------------------------------------------------------*
 |  adds two Epetra_SerialDenseVector                        u.may 06/07|
 *----------------------------------------------------------------------*/
Epetra_SerialDenseVector Intersection::addTwoVectors(   Epetra_SerialDenseVector&   v1,
                                                        Epetra_SerialDenseVector&   v2)
{   
    if(v1.Length() != v2.Length())
        dserror("both vectors need to have the same size\n"); 

    for(int i = 0; i < v1.Length(); i++)
        v1[i] = v1[i] + v2[i];
 
    return v1;
}
	
    

/*----------------------------------------------------------------------*
 |  subtracts one vector<double> from another vector<double> u.may 06/07|
 |  The result is stored in v1                                          |
 *----------------------------------------------------------------------*/
std::vector<double> Intersection::subtractsTwoVectors(	std::vector <double>& v1, 
														std::vector <double>& v2)
{   
    if(v1.size() != v2.size())
        dserror("both vectors need to have the same size\n"); 

    for(unsigned int i = 0; i < v1.size(); i++)
        v1[i] = v1[i] - v2[i];
 
    return v1;
}



/*----------------------------------------------------------------------*
 |  computes an approximate axis-aligned bounding box        u.may 06/07|
 |  for all types of elements                                           |
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix Intersection::computeFastAABB(DRT::Element* element)
{
    double  maxDistance; 
    const double TOL = EPS8;
    Epetra_SerialDenseMatrix AABB(3, 2);
    
	for(int dim=0; dim<3; dim++)
	{
		AABB(dim, 0) = element->Nodes()[0]->X()[dim] - TOL;
   		AABB(dim, 1) = element->Nodes()[0]->X()[dim] + TOL;
	}
    
    for(int i=1; i<element->NumNode(); i++)
        for(int dim=0; dim<3; dim++)
		{
            AABB(dim, 0) = std::min( AABB(dim, 0), element->Nodes()[i]->X()[dim] - TOL);
			AABB(dim, 1) = std::max( AABB(dim, 1), element->Nodes()[i]->X()[dim] + TOL);
		}
 
    maxDistance = fabs(AABB(0,1) - AABB(0,0));
 	for(int dim=1; dim<3; dim++)
	   maxDistance = std::max(maxDistance, fabs(AABB(dim,1)-AABB(dim,0)) );
	
    // subtracts half of the maximal distance to minX, minY, minZ
    // adds half of the maximal distance to maxX, maxY, maxZ 
	for(int dim=0; dim<3; dim++)
	{
		AABB(dim, 0) = AABB(dim, 0) - 0.5*maxDistance;
		AABB(dim, 1) = AABB(dim, 1) + 0.5*maxDistance;
	}	
	
	/*
    printf("\n");
	printf("axis-aligned bounding box:\n minX = %f\n minY = %f\n minZ = %f\n maxX = %f\n maxY = %f\n maxZ = %f\n", 
			  AABB(0,0), AABB(1,0), AABB(2,0), AABB(0,1), AABB(1,1), AABB(2,1));
	printf("\n");
	*/
    
	return AABB;
}



/*----------------------------------------------------------------------*
 |  checks if a node is within an AABB                       u.may 06/07|
 *----------------------------------------------------------------------*/
bool Intersection::isNodeWithinAABB(    std::vector<double>&         node,
                                        Epetra_SerialDenseMatrix&    AABB)
{
	bool isWithin = true;
    const double TOL = -EPS8;
    double diffMin = 0;
    double diffMax = 0;
	
	for (int dim=0; dim<3; dim++)
	{
        diffMin = node[dim] - AABB(dim,0);
        diffMax = AABB(dim,1) - node[dim];
        
   	    if((diffMin < TOL)||(diffMax < TOL)) //check again !!!!!
            isWithin = false;
    }
		 	
	return isWithin;
}



/*----------------------------------------------------------------------*
 |  checks if two AABB's intersect                           u.may 06/07|
 *----------------------------------------------------------------------*/
bool Intersection::intersectionOfAABB(  Epetra_SerialDenseMatrix& cutterAABB, 
										Epetra_SerialDenseMatrix& xfemAABB)
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
	
	nodes[0][0] = cutterAABB(0,0);	nodes[0][1] = cutterAABB(1,0);	nodes[0][2] = cutterAABB(2,0);	// node 0	
	nodes[1][0] = cutterAABB(0,1);	nodes[1][1] = cutterAABB(1,0);	nodes[1][2] = cutterAABB(2,0);	// node 1
	nodes[2][0] = cutterAABB(0,1);	nodes[2][1] = cutterAABB(1,1);	nodes[2][2] = cutterAABB(2,0);	// node 2
	nodes[3][0] = cutterAABB(0,0);	nodes[3][1] = cutterAABB(1,1);	nodes[3][2] = cutterAABB(2,0);	// node 3
	nodes[4][0] = cutterAABB(0,0);	nodes[4][1] = cutterAABB(1,0);	nodes[4][2] = cutterAABB(2,1);	// node 4
	nodes[5][0] = cutterAABB(0,1);	nodes[5][1] = cutterAABB(1,0);	nodes[5][2] = cutterAABB(2,1);	// node 5
	nodes[6][0] = cutterAABB(0,1);	nodes[6][1] = cutterAABB(1,1);	nodes[6][2] = cutterAABB(2,1);	// node 6
	nodes[7][0] = cutterAABB(0,0);	nodes[7][1] = cutterAABB(1,1);	nodes[7][2] = cutterAABB(2,1);	// node 7
	
	for (int i = 0; i < 8; i++)
		if(isNodeWithinAABB(nodes[i], xfemAABB))
		{
			intersection = true;
			break;
		}
	
		
	if(!intersection)
	{
		nodes[0][0] = xfemAABB(0,0);	nodes[0][1] = xfemAABB(1,0);	nodes[0][2] = xfemAABB(2,0);	// node 0	
		nodes[1][0] = xfemAABB(0,1);	nodes[1][1] = xfemAABB(1,0);	nodes[1][2] = xfemAABB(2,0);	// node 1
		nodes[2][0] = xfemAABB(0,1);	nodes[2][1] = xfemAABB(1,1);	nodes[2][2] = xfemAABB(2,0);	// node 2
		nodes[3][0] = xfemAABB(0,0);	nodes[3][1] = xfemAABB(1,1);	nodes[3][2] = xfemAABB(2,0);	// node 3
		nodes[4][0] = xfemAABB(0,0);	nodes[4][1] = xfemAABB(1,0);	nodes[4][2] = xfemAABB(2,1);	// node 4
		nodes[5][0] = xfemAABB(0,1);	nodes[5][1] = xfemAABB(1,0);	nodes[5][2] = xfemAABB(2,1);	// node 5
		nodes[6][0] = xfemAABB(0,1);	nodes[6][1] = xfemAABB(1,1);	nodes[6][2] = xfemAABB(2,1);	// node 6
		nodes[7][0] = xfemAABB(0,0);	nodes[7][1] = xfemAABB(1,1);	nodes[7][2] = xfemAABB(2,1);	// node 7
	
		for (int i = 0; i < 8; i++)
			if(isNodeWithinAABB(nodes[i], cutterAABB))
			{
				intersection = true;
				break;
			}
	}	
	return intersection;
}



/*----------------------------------------------------------------------*
 |  computes a Gaussian Elimination for a linear system of equations    |
 |                                                           u.may 06/07|
 *----------------------------------------------------------------------*/
bool Intersection::gaussElimination(    Epetra_SerialDenseMatrix&   A,
                                        Epetra_SerialDenseVector&   b,
                                        Epetra_SerialDenseVector&   x,
                                        bool                        do_piv,
                                        const int                   dim,
                                        const int 					order)
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
    if(fabs(det) < EPS8 && order == 1)
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
 |  updates the Jacobi matrix for the computation if a node is in       |
 |  a given element                                          u.may 06/07|
 *----------------------------------------------------------------------*/
void Intersection::updateAForNWE(   Epetra_SerialDenseMatrix&   A,
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
 |  updates the rhs for the computation if a node is in                 |
 |  a given element                                          u.may 06/07|
 *----------------------------------------------------------------------*/
void Intersection::updateRHSForNWE( Epetra_SerialDenseVector&   b,
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
 |  checks if a node is within a given element               u.may 06/07|								|
 *----------------------------------------------------------------------*/
bool Intersection::checkNodeWithinElement(	DRT::Element* element,
 											Epetra_SerialDenseVector& x,
                                            Epetra_SerialDenseVector& xsi)
{

	bool nodeWithinElement = true;
    int iter = 0;
    const int maxiter = 500;
    double residual = 1.0;
    const double TOL = EPS8;
    Epetra_SerialDenseMatrix A(3,3);
    Epetra_SerialDenseVector b(3);
    Epetra_SerialDenseVector dx(3);
  
    xsi[0] = 0.0; xsi[1] = 0.0; xsi[2] = 0.0;
    dx = xsi;
        
    updateRHSForNWE( b, xsi, x, element);
    
    while(residual > TOL)
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
    
    if( (fabs(xsi[0])-1.0) > TOL  || (fabs(xsi[1])-1.0) > TOL  || (fabs(xsi[2])-1.0) > TOL )   
        nodeWithinElement = false;
    
    return nodeWithinElement;
}



/*----------------------------------------------------------------------*
 |  checks if a node that lies within an element lies on     u.may 06/07|
 |  one of its surfaces or nodes                                        |
 *----------------------------------------------------------------------*/
bool Intersection::checkIfOnSurfaceAndNode( DRT::Element*                   element,
                                            DRT::Node*                      node,
                                            Epetra_SerialDenseVector&       xsi, 
                                            InterfacePoint&                 ip,
                                            int&                            numSurfacePoints)
{
    bool onSurface = false;
    const double TOL = EPS8; 
    int count = 0;
    
    for(int i = 0; i < 3; i++)
        if(fabs(fabs(xsi[i])-1.0) < TOL)
            count++;
     
    
    if(count == 1)  
    {
        onSurface = true;
        numSurfacePoints++;
        ip.nsurf = 1;
        ip.pType = surfaceP;
        
       
        // only for hex
        if(     fabs(xsi[0]-1.0) < TOL)                     ip.surfaces[0] = 2;        
        else if(fabs(xsi[0]+1.0) < TOL)                     ip.surfaces[0] = 4;        
        else if(fabs(xsi[1]-1.0) < TOL)                     ip.surfaces[0] = 3;        
        else if(fabs(xsi[1]+1.0) < TOL)                     ip.surfaces[0] = 1;        
        else if(fabs(xsi[2]-1.0) < TOL)                     ip.surfaces[0] = 5;        
        else if(fabs(xsi[2]+1.0) < TOL)                     ip.surfaces[0] = 0;        
        else dserror("node does not lie a surface");
            
    }
    else if(count == 2)
    {   
        onSurface = true;
        numSurfacePoints++;
        ip.nsurf = 2;
        ip.pType = surfaceP;
        int countSurf = 0;
        
        // only hex
        if(fabs(xsi[0]-1.0) < TOL)      ip.surfaces[countSurf++] = 2;        
        if(fabs(xsi[0]+1.0) < TOL)      ip.surfaces[countSurf++] = 4;        
        if(fabs(xsi[1]-1.0) < TOL)      ip.surfaces[countSurf++] = 3;        
        if(fabs(xsi[1]+1.0) < TOL)      ip.surfaces[countSurf++] = 1;        
        if(fabs(xsi[2]-1.0) < TOL)      ip.surfaces[countSurf++] = 5;        
        if(fabs(xsi[2]+1.0) < TOL)      ip.surfaces[countSurf++] = 0;        
        
    }
    // coincident with an element node
    else if(count == 3)
    {
        onSurface = true;
        numSurfacePoints++;
        ip.pType = surfaceP;
        
        for(int i = 0; i < element->NumNode(); i++)
            if(comparePoints(node->X(), element->Nodes()[i]->X(), 3))
            {
                ip.nsurf = 3;
                for(int j = 0; j < 3; j++)
                    ip.surfaces[j] = eleNodeNumbering_hex27_nodes_surfaces[i][j];
                   
                break;
            }
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
 |  collects points that belong to the interface and lie 	 u.may 06/07|
 |	within an xfem element             									|
 *----------------------------------------------------------------------*/
bool Intersection::collectInternalPoints(  	DRT::Element*                   element,
                                            DRT::Node*                      node,
                                            std::vector< InterfacePoint >&  interfacePoints,
                                            int                             elemId,
                                            int                             nodeId,
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
        checkIfOnSurfaceAndNode(element, node, xsi, ip, numSurfacePoints);
        
        
        switch(nodeId)
        {
            case 0:
            {
                ip.coord[0] = -1.0; 
                ip.coord[1] = -1.0;
                break;
            }
            case 1:
            {
                ip.coord[0] =  1.0; 
                ip.coord[1] = -1.0;
                break;
            }
            case 2:
            {
                ip.coord[0] =  1.0; 
                ip.coord[1] =  1.0;
                break;
            }
            case 3:
            {
                ip.coord[0] = -1.0; 
                ip.coord[1] =  1.0;
                break;
            }
            default:
                dserror("node number not correct");    
        }
        ip.coord[2] = 0.0; 
        
        interfacePoints.push_back(ip);
    }
      
    return nodeWithinElement;
}



/*----------------------------------------------------------------------*
 |  updates the systemmatrix for the computation of          u.may 06/07|
 |  a curve-surface intersection (CSI)                                  |
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
 |  updates the right-hand-side for the computation of       u.may 06/07|
 |  a curve-surface intersection (CSI)                                  |
 *----------------------------------------------------------------------*/
void Intersection::updateRHSForCSI( Epetra_SerialDenseVector&   b,
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
 |  computes the intersection between a curve and            u.may 06/07|
 |  a surface                    (CSI)                                  |
 *----------------------------------------------------------------------*/
bool Intersection::computeCurveSurfaceIntersection( DRT::Element*               surfaceElement,
                                                    DRT::Element*               lineElement,
                                                    Epetra_SerialDenseVector&   xsi,
                                                    Epetra_SerialDenseVector&   upLimit,
                                                    Epetra_SerialDenseVector&   loLimit)
{
    int iter = 0;
	int maxiter = 500;
	bool intersection = true;
    double residual = 1.0;
    const double TOL = EPS8;
    Epetra_SerialDenseMatrix A(3,3);
    Epetra_SerialDenseVector b(3);
    Epetra_SerialDenseVector dx(3);
    dx = xsi;
    
	updateRHSForCSI( b, xsi, surfaceElement, lineElement);
             					
	while(residual > 1e-14)
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
    
    if( (xsi[0] > upLimit[0]+TOL) || (xsi[1] > upLimit[1]+TOL) || (xsi[2] > upLimit[2]+TOL)  || 
        (xsi[0] < loLimit[0]-TOL) || (xsi[1] < loLimit[1]-TOL) || (xsi[2] < loLimit[2]-TOL)) 
    	intersection = false;
      	
	return intersection;
}



/*----------------------------------------------------------------------*
 |  collects all intersection points of a line and           u.may 06/07|
 |  a surface                                                           |
 *----------------------------------------------------------------------*/  
void Intersection::collectIntersectionPoints(   DRT::Element*                   surfaceElement,
					 							DRT::Element*                   lineElement,
					 							std::vector<InterfacePoint>&    interfacePoints,
	                                            int                            	numInternalPoints,
                                                int                             numSurfacePoints,
	                                           	int                             surfaceId,
	        									int                             lineId,
	                                            bool                            lines,
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
        numInterfacePoints = addIntersectionPoint(	surfaceElement, lineElement,xsi, upLimit, loLimit, 
        											interfacePoints, surfaceId, lineId, lines);
      
    
    // in this case a node of this line lies on the facet of the xfem element
    // but there is no intersection within the element											
    if(!((int) interfacePoints.size() == numSurfacePoints)) 
        xfemIntersection = true;
      
}



/*----------------------------------------------------------------------*
 |  computes a new starting point for the Newton-method      u.may 06/07|
 |  in order to find all intersection points                            |
 |  of a curve-surface intersection                                     |
 *----------------------------------------------------------------------*/  
int Intersection::computeNewStartingPoint(	DRT::Element*                surfaceElement,
											DRT::Element*                lineElement,
                                          	int                          surfaceId,
        									int                          lineId,
											Epetra_SerialDenseVector&    upLimit,
     										Epetra_SerialDenseVector&    loLimit,
     										std::vector<InterfacePoint>& interfacePoints,
                                            bool                         lines)
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
 |  adds an intersection point to the 					     u.may 07/07|
 |  list of interface points                       						|
 *----------------------------------------------------------------------*/  
int Intersection::addIntersectionPoint(	DRT::Element*                	surfaceElement,
										DRT::Element*                	lineElement,
										Epetra_SerialDenseVector&		xsi,
										Epetra_SerialDenseVector&     	upLimit,
										Epetra_SerialDenseVector&     	loLimit,
										std::vector<InterfacePoint>& 	interfacePoints,
										int                      		surfaceId,
			        					int                       		lineId,							
										bool 							lines)
{

	int numInterfacePoints = 0;
	
	
 	InterfacePoint ip;
    if(lines)
    {   
        ip.nsurf = 1;
        ip.surfaces[0] = surfaceId;
         
        // change minus sign if you change the line numbering   
        switch(lineId)
        {
            case 0:
            {
                ip.coord[0] = xsi[2]; 
                ip.coord[1] = -1.0;
                break;
            }
            case 1:
            {
                ip.coord[0] =  1.0; 
                ip.coord[1] = xsi[2];
                break;
            }
            case 2:
            {
                ip.coord[0] =  -xsi[2];    // changed
                ip.coord[1] =  1.0;
                break;
            }
            case 3:
            {
                ip.coord[0] = -1.0; 
                ip.coord[1] =  xsi[2];
                break;
            }
            default:
                dserror("node number not correct");    
        }           
    }
    else
    {
        ip.nsurf = 2;
        // change for all distypes
        ip.surfaces[0] = eleNodeNumbering_hex27_lines_surfaces[lineId][0];
        ip.surfaces[1] = eleNodeNumbering_hex27_lines_surfaces[lineId][1];
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
    
    numInterfacePoints =  numInterfacePoints +
      computeNewStartingPoint(	surfaceElement, lineElement, surfaceId, lineId, 
      							upLimit, xsi, interfacePoints, lines) +
      computeNewStartingPoint(	surfaceElement, lineElement, surfaceId, lineId, 
      							xsi, loLimit, interfacePoints, lines);

	return numInterfacePoints;
}



/*----------------------------------------------------------------------*
 |  transforms a node in reference coordinates            u.may 06/07|
 |  into current coordinates                                            |
 *----------------------------------------------------------------------*/  
void Intersection::referenceToCurrentCoordinates(   DRT::Element* element, 
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
 |  transforms a node in current coordinates            u.may 06/07|
 |  into reference coordinates                                            |
 *----------------------------------------------------------------------*/  
void Intersection::currentToReferenceCoordinates(   DRT::Element* element, 
                                                    Epetra_SerialDenseVector& xsi)
{
    bool nodeWithinElement;
    double TOL = EPS8;
    Epetra_SerialDenseVector x(3);
    
    for(int i = 0; i < 3; i++)
    {
        x[i] = xsi[i];
        xsi[i] = 0.0;
    }
    
    nodeWithinElement = checkNodeWithinElement(element, x, xsi);
    
    if(!nodeWithinElement)
        dserror("node not within element");
        
        
    // rounding 1 and -1 to be exact for the CDT
    for(int j = 0; j < 3; j++)
    {
        if( fabs((fabs(xsi[j])-1.0)) < TOL &&  xsi[j] < 0)    xsi[j] = -1.0;
        if( fabs((fabs(xsi[j])-1.0)) < TOL &&  xsi[j] > 0)    xsi[j] =  1.0;      
    }  
}     



/*----------------------------------------------------------------------*
 |  compares two nodes                                       u.may 06/07|
 |  overloaded method:                                                  |
 |  double*  and  double*                                               |
 *----------------------------------------------------------------------*/  
bool Intersection::comparePoints(    const double*     point1,
                                     const double*     point2,
                                     const int         length)
{   
    bool equal = true;
             
    for(int i = 0; i < length; i++)
        if(fabs(point1[i] - point2[i]) > EPS8)
        {
            equal = false;
            break;
        }
  
    return equal;
}



/*----------------------------------------------------------------------*
 |  compares two nodes                                       u.may 06/07|
 |  overloaded method:  vector<double>  and double*                     |
 *----------------------------------------------------------------------*/  
bool Intersection::comparePoints(   const vector<double>&     point1,
                                    const double*             point2)
{   
    bool equal = true;
    
    if(point2 == NULL)
        dserror("array is NULL");
        
    for(unsigned int i = 0; i < point1.size() ; i++)
        if(fabs(point1[i] - point2[i]) > EPS8)
        {
            equal = false;
            break;
        }
    
    return equal;
}



/*----------------------------------------------------------------------*
 |  compares two nodes                                       u.may 06/07|
 |  overloaded method:  vector<double>  and  vector<double>             |
 *----------------------------------------------------------------------*/  
bool Intersection::comparePoints(	const vector<double>& point1,
                                    const vector<double>& point2)
{   
    bool equal = true;
    
    if(point1.size() != point2.size())
        dserror("arrays of nodes need to have the same length");
             
    for(unsigned int i = 0; i < point1.size() ; i++)
        if(fabs(point1[i] - point2[i]) > EPS8)
        {
            equal = false;
            break;
        }
  
    return equal;
}



/*----------------------------------------------------------------------*
 |  compares two nodes                                       u.may 06/07|
 |  overloaded method:                                                  |
 |  Epetra_SerialDenseVector  and  Epetra_SerialDenseVector             |
 *----------------------------------------------------------------------*/  
bool Intersection::comparePoints(    const Epetra_SerialDenseVector&     point1,
                                     const Epetra_SerialDenseVector&     point2)
{   
    bool equal = true;
    
    if(point1.Length() != point2.Length())
        dserror("arrays of nodes need to have the same length");
             
    for(int i = 0; i < point1.Length() ; i++)
        if(fabs(point1[i] - point2[i]) > EPS8)
        {
            equal = false;
            break;
        }
  
    return equal;
}



/*----------------------------------------------------------------------*
 |  checks if a certain element is a volume  element         u.may 06/07|
 |  with help of the discretization type                                |
 *----------------------------------------------------------------------*/  
bool Intersection::checkIfVolumeElement(DRT::Element* element)
{
	bool isVolume = false;
	DRT::Element::DiscretizationType distype = element->Shape();
	
	if(	distype == DRT::Element::hex8  ||
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
 |  checks if a certain element is a surface element         u.may 06/07|
 |  with help of the discretization type                                |
 *----------------------------------------------------------------------*/  
bool Intersection::checkIfSurfaceElement(DRT::Element* element)
{
	bool isSurface = false;
	DRT::Element::DiscretizationType distype = element->Shape();
	
	if(	distype == DRT::Element::quad4 ||
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
 |  checks if a certain element is a line element            u.may 06/07|
 |  with help of the discretization type                                |
 *----------------------------------------------------------------------*/  
bool Intersection::checkIfLineElement(DRT::Element* element)
{
	bool isLine = false;
	DRT::Element::DiscretizationType distype = element->Shape();
	
	if(	distype == DRT::Element::line2 ||
		distype == DRT::Element::line3 )
	{
		isLine = true;		
	}
	return isLine;
}



/*----------------------------------------------------------------------*
 |  returns the order of the element			             u.may 06/07|
 *----------------------------------------------------------------------*/  
int Intersection::computeOrder(DRT::Element* element)
{
	int order = 0;
	DRT::Element::DiscretizationType distype = element->Shape();
	
	if(	distype == DRT::Element::line2 ||
		distype == DRT::Element::line3 )
	{
		order = 1;		
	}
	if(	distype == DRT::Element::quad4 ||
		distype == DRT::Element::quad8 ||
		distype == DRT::Element::quad9 ||
		distype == DRT::Element::tri3  ||
		distype == DRT::Element::tri6  )
	{
		order = 2;		
	}
	else if(	distype == DRT::Element::hex8  ||
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
 |  finds the next facet of a convex hull                    u.may 06/07|
 |  and returns the point different form the searchpoint                |
 *----------------------------------------------------------------------*/  
void Intersection::findNextFacet(   vector< vector<double> >&   vertices, 
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
 |  stores a point within a list of points                   u.may 06/07|
 |  which is to be copy to the tetgen data structure                    |
 |  for the computation of the Constrained Delauney Triangulation       |
 *----------------------------------------------------------------------*/  
void Intersection::storePoint(  vector<double>&             point, 
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
 |  stores a point lying on a surface of an                  u.may 06/07|
 |  xfem element                                                        |
 *----------------------------------------------------------------------*/  
void Intersection::storeSurfacePoints(  vector<InterfacePoint>&     pointList, 
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
 |  stores a segment within a list of segments               u.may 06/07|
 |  which is to be copied to the tetgen data structure                  |
 |  for the computation of the Constrained Delauney Triangulation       |
 *----------------------------------------------------------------------*/  
void Intersection::storeSegments(   vector<InterfacePoint>&     pointList, 
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
 |  stores a triangle within a list of trianles              u.may 06/07|
 |  which is to be copy to the tetgen data structure                    |
 |  for the computation of the Constrained Delauney Triangulation       |
 *----------------------------------------------------------------------*/  
void Intersection::storeTriangles(  vector<InterfacePoint>&     pointList, 
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
 |  stores a triangles lying within a xfem facet             u.may 06/07|
 |  in a polygon list                                                   |
 |  for the computation of the Constrained Delauney Triangulation       |
 *----------------------------------------------------------------------*/  
void Intersection::storePolygons(   vector<InterfacePoint>&     pointList, 
                                    vector<int>                 positions, 
                                    int                         surface,
                                    vector< vector<int> >&      polygonSurfaceList)
{
    
    for(unsigned int i = 0; i < positions.size()-1; i++ )
    {
        polygonSurfaceList[surface].push_back( positions[i] );
        polygonSurfaceList[surface].push_back( positions[i+1] ) ;
        polygonSurfaceList[surface].push_back( pointList.size()-1 );
    }
    
    polygonSurfaceList[surface].push_back( positions[positions.size()-1] );
    polygonSurfaceList[surface].push_back( positions[0] );
    polygonSurfaceList[surface].push_back( pointList.size()-1 );
}



/*----------------------------------------------------------------------*
 |  computes the midpoint of a collection of InterfacePoints u.may 06/07|
 *----------------------------------------------------------------------*/  
InterfacePoint Intersection::computeMidpoint(vector<InterfacePoint>& interfacePoints)
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
 |  computes surface id for a cutter element lying in a      u.may 06/07|
 |  facet of an xfemelement                                             |
 *----------------------------------------------------------------------*/  
int  Intersection::determineSurfaceId(  InterfacePoint& midpoint)
{

    int surfaceId = 0;
    double TOL = EPS8;
       
    if(     fabs(midpoint.coord[0]-1.0) < TOL)                     surfaceId = 2;        
    else if(fabs(midpoint.coord[0]+1.0) < TOL)                     surfaceId = 4;        
    else if(fabs(midpoint.coord[1]-1.0) < TOL)                     surfaceId = 3;        
    else if(fabs(midpoint.coord[1]+1.0) < TOL)                     surfaceId = 1;        
    else if(fabs(midpoint.coord[2]-1.0) < TOL)                     surfaceId = 5;        
    else if(fabs(midpoint.coord[2]+1.0) < TOL)                     surfaceId = 0;        
    else dserror("node does not lie a surface");
    
    return surfaceId;
}



/*----------------------------------------------------------------------*
 |  computes the coordinates of a point within a region for  u.may 07/07|
 |  the Tetgen data structure		(only for debugging)				|
 *----------------------------------------------------------------------*/  
void Intersection::computeRegionCoordinates(	DRT::Element*  xfemElement,
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



/*----------------------------------------------------------------------*
 |  computes the convex hull of a set of interface points    u.may 06/07|
 |  and stores resulting points, segments and triangles                 |
 |  for the use with Tetgen (CDT)                                       |
 *----------------------------------------------------------------------*/  
void Intersection::computeConvexHull(   DRT::Element*           xfemElement,
                                        DRT::Element*           surfaceElement,
                                        vector<InterfacePoint>& interfacePoints,
                                        vector<InterfacePoint>& pointList,
                                        vector< vector<int> >&  surfacePointList,
                                        vector< vector<int> >&  segmentList,
                                        vector< vector<int> >&  triangleList,
                                        int                     numInternalPoints,
                                        int                     numSurfacePoints )
{
    double* point;
    vector<int> positions;
    vector<double> searchPoint(3,0);
    vector<double> vertex(3,0);
    vector< vector<double> > vertices;  
    Epetra_SerialDenseVector curCoord(3);  
    InterfacePoint midpoint;
    
           
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
            findNextFacet(vertices, searchPoint);
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
 |  fills the point list with the corner points              u.may 06/07|
 |  in reference coordinates of the xfem element                        |
 *----------------------------------------------------------------------*/  
void Intersection::startPointList(vector< InterfacePoint >&  pointList)
{
    int cornerpoints = 8;
    InterfacePoint ip;
        
    for(int i = 0; i < cornerpoints; i++)
    {
        ip.nsurf = 3;
        
        // change for other element types
        for(int j = 0; j < 3; j++) 
        { 
           ip.coord[j] = eleNodeNumbering_hex27_nodes_reference[i][j]; 
           ip.surfaces[j] = eleNodeNumbering_hex27_nodes_surfaces[i][j]; 
        }
        pointList.push_back(ip);               
    }
}



/*----------------------------------------------------------------------*
 |  computes the Constrained Delaunay Tetrahedralization     u.may 06/07|
 |  in 3D with help of Tetgen library for an intersected                |
 |  xfem element                                                        |
 *----------------------------------------------------------------------*/  
void Intersection::computeCDT(  DRT::Element*               			element,
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
    char switches[] = "po2QY";    //YA for regions Q quiet
    tetgenio::facet *f;
    tetgenio::polygon *p;
    double regionCoordinates[6];
    

    // set points
    in.numberofpoints = pointList.size();
    in.pointlist = new REAL[in.numberofpoints * dim];
       
    // fill point list
    int fill = 0;
    for(int i = 0; i <  in.numberofpoints; i++)
        for(int j = 0; j < dim; j++)  
            in.pointlist[fill++] = (REAL) pointList[i].coord[j]; 
 
    // set facets
    if(triangleList.size()>0)       in.numberoffacets = 6 + triangleList.size(); 
    else                            in.numberoffacets = 6;   
      
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    
    
    // loop over all xfem element surfaces
    for(int i = 0; i < element->NumSurface(); i++)
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
            p->vertexlist[j] = eleNodeNumbering_hex27_surfaces[i][j];
           
      
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
    for(int i = element->NumSurface(); i < in.numberoffacets; i++)
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
  
        
    // set facetmarkers
    for(int i = 0; i < in.numberoffacets; i ++)
        in.facetmarkerlist[i] = 0;   


	// specify regions
	bool regions = false;	
	if(regions)
	{
		computeRegionCoordinates(element,cutterElement,regionCoordinates);
	    fill = 0;
	    int read = 0;
	    in.numberofregions = 2;
	    in.regionlist = new REAL[in.numberofregions*5];
	  	for(int i = 0; i < in.numberofregions; i++)
	    {
	    	// store coordinates 
	    	for(int j = 0; j < 3; j++)
	    		in.regionlist[fill++] = regionCoordinates[read++];
	    		
	    	// store regional attribute (switch A)   i=0 cutter i=1 fluid
	    	in.regionlist[fill++] = i;
	    	
	    	// store volume constraint (switch a)   
	    	in.regionlist[fill++] = 0.0;
		}
	}
	      
    //  Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
    //  do quality mesh generation (q) with a specified quality bound
    //  (1.414), and apply a maximum volume constraint (a0.1)
    tetrahedralize(switches, &in, &out); 
  
    //Debug
    vector<int> elementIds;
    for(int i = 0; i<8; i++)
        elementIds.push_back(i);
    
    debugTetgenOutput(in, out, element, elementIds);
 
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
 |  computes the intersection between two different          u.may 06/07|
 |  discretizations and returns a list of intersected xfem elements     |
 |  and their integrations cell                                         |
 *----------------------------------------------------------------------*/  
void Intersection::computeIntersection( RefCountPtr<DRT::Discretization> actdis, 
                                        map< int, vector <Integrationcell> >&	integrationcellList)
{
    bool intersected =                      false;
    bool xfemIntersection =                 false;
    int numInternalPoints = 				0;
    int numSurfacePoints =                  0;
    vector< DRT::Condition * >              xfemConditions;
    vector< InterfacePoint >                interfacePoints;
    vector <vector< InterfacePoint > >      interfacePointCollection;
    map< int, RefCountPtr<DRT::Element > >  geometryMap;
    DRT::Element*                           cutterElement; 
    DRT::Element*                           xfemElement; 
    Epetra_SerialDenseMatrix                cutterAABB;
    Epetra_SerialDenseMatrix                xfemAABB;
    vector< InterfacePoint >                pointList;
    vector< vector<int> >                   segmentList(6);                 // adjust for all element types
    vector< vector<int> >                   surfacePointList(6);                 // adjust for all element types
    vector< vector<int> >                   triangleList;
    
    
        
    // obtain vector of pointers to all xfem conditions
    actdis->GetCondition ("XFEMCoupling", xfemConditions);
    //cout << endl << "number of xfem conditions =" << xfemConditions.size() << endl; cout << endl;
        
    if(xfemConditions.size()==0)
        dserror("number of fsi xfem conditions = 0");
        
    /* for(unsigned int i=0; i<xfemConditions.size(); i++)  cout << *xfemConditions[i]; */   

       
    //  k<actdis->NumMyRowElements()-1
    for(int k=0; k<actdis->NumMyRowElements()-1; k++)
    {
        xfemElement = actdis->gElement(k);
        xfemAABB = computeFastAABB(xfemElement);
        xfemIntersection = false;
        
        pointList.clear();
        startPointList(pointList);
        
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
            
                cutterAABB = computeFastAABB(cutterElement);                                 
                intersected = intersectionOfAABB(cutterAABB, xfemAABB);    
                //debugAABBIntersection( cutterAABB, xfemAABB, cutterElement, xfemElement, i, k);
                                     
                if(intersected)
                {
                	// collect internal points
                	numInternalPoints= 0; 
                    numSurfacePoints = 0;
                    for(int m=0; m<cutterElement->NumLine() ; m++)                    
                        collectInternalPoints( xfemElement, cutterElement->Nodes()[m],
                                                interfacePoints, k, m, numInternalPoints, numSurfacePoints);
                    
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
        for(unsigned int i  = 0; i < 6; i++)
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
 |  Debug only                                               u.may 06/07|
 *----------------------------------------------------------------------*/  
void Intersection::debugAABBIntersection( 	Epetra_SerialDenseMatrix cutterAABB, 
											Epetra_SerialDenseMatrix xfemAABB,
											DRT::Element* cutterElement,
				 							DRT::Element* xfemElement,
				 							int noC,
				 							int noX)
{
	cout << endl;
    cout << "===============================================================" << endl;
	cout << "Debug Intersection of AABB's" << endl;
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
    cout << "CUTTER AABB: "<< "                     " << "XFEM AABB: " << endl;
    cout << endl;
    cout <<  "minX = " << cutterAABB(0,0) << "      " << "maxX = " << cutterAABB(0,1) << "      " ;
    cout <<  "minX = " << xfemAABB(0,0)   << "      " << "maxX = " << xfemAABB(0,1)   << endl;  
    cout <<  "minY = " << cutterAABB(1,0) << "      " << "maxY = " << cutterAABB(1,1) << "      " ;
    cout <<  "minY = " << xfemAABB(1,0)   << "      " << "maxY = " << xfemAABB(1,1)   << endl;   
    cout <<  "minZ = " << cutterAABB(2,0) << "      " << "maxZ = " << cutterAABB(2,1) << "      " ;
    cout <<  "minZ = " << xfemAABB(2,0)   << "      " << "maxZ = " << xfemAABB(2,1) << endl;
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
    cout << "CUTTER AABB: "<< "                     " << "XFEM AABB: " << endl;
    cout << endl;
    cout <<  "minX = " << cutterAABB(0,0) << "      " << "maxX = " << cutterAABB(0,1) << "      " ;
    cout <<  "minX = " << xfemAABB(0,0)   << "      " << "maxX = " << xfemAABB(0,1)   << endl;  
    cout <<  "minY = " << cutterAABB(1,0) << "      " << "maxY = " << cutterAABB(1,1) << "      " ;
    cout <<  "minY = " << xfemAABB(1,0)   << "      " << "maxY = " << xfemAABB(1,1)   << endl;   
    cout <<  "minZ = " << cutterAABB(2,0) << "      " << "maxZ = " << cutterAABB(2,1) << "      " ;
    cout <<  "minZ = " << xfemAABB(2,0)   << "      " << "maxZ = " << xfemAABB(2,1) << endl;
	cout << endl;
	cout << endl;    
    cout << "===============================================================" << endl;
    cout << "End Debug Intersection of AABB's" << endl;
	cout << "===============================================================" << endl;
	cout << endl; cout << endl; cout << endl;


}



/*----------------------------------------------------------------------*
 |  Debug only                                               u.may 06/07|
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
 |  Debug only                                               u.may 06/07|
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
 |  Debug only                                               u.may 06/07|
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



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_XFEM



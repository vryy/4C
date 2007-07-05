/*!
\file intersection.cpp

\brief collection of intersection tool 

<pre>
Maintainer: Ursula Mayer
</pre>
*/

#ifdef XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "intersection.H"
#include "../drt_xfem/integrationcell.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"

#include "../headers/definitions.h"     //remove



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
                                        const int                   dim)
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
    if(fabs(det) < EPS8)
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
 |  checks if a node is in a given element                   u.may 06/07|
 *----------------------------------------------------------------------*/
bool Intersection::checkNodeWithinElement(  DRT::Element*                   element,
                                            DRT::Node*                      node,
                                            std::vector< InterfacePoint >&  interfacePoints,
                                            int                             elemId,
                                            int                             nodeId,
                                            int&                            numInterfacePoints)
{
    bool nodeWithinElement = true;
    int iter = 0;
    const int maxiter = 500;
    double residual = 1.0;
    const double TOL = EPS8;
    Epetra_SerialDenseMatrix A(3,3);
    Epetra_SerialDenseVector b(3);
    Epetra_SerialDenseVector dx(3);
    Epetra_SerialDenseVector x(3);
    Epetra_SerialDenseVector xsi(3);
       
    x[0] = node->X()[0];
    x[1] = node->X()[1];
    x[2] = node->X()[2];
    xsi[0] = 0.0; xsi[1] = 0.0; xsi[2] = 0.0;
    dx = xsi;
        
    updateRHSForNWE( b, xsi, x, element);
    
    while(residual > TOL)
    { 
        updateAForNWE( A, xsi, element);
        
        if(!gaussElimination(A, b, dx, true, 3))
        {
            nodeWithinElement = false;
            break;
        }   
        xsi = addTwoVectors(xsi,dx);
        //cout << "xsi: x = " << xsi[0] << "  y = " <<  xsi[1] << "  z = "  << xsi[2] << endl;
      
        if( (fabs(xsi[0]) > (1.0+TOL) ) || (fabs(xsi[1]) > (1.0+TOL) ) || (fabs(xsi[2]) > (1.0+TOL) ) )   
        {
            nodeWithinElement = false;
            break;
        }
        
        if(iter >= maxiter)
        {   
            nodeWithinElement = false;
            break;
        }       
        updateRHSForNWE(b, xsi, x, element);
        residual = b.Norm2();
        iter++;
    }
    
    if(nodeWithinElement)
    {
        numInterfacePoints++; 
        InterfacePoint ip;
        //debugNodeWithinElement(element,node,xsi,elemId ,nodeId, nodeWithinElement);  
        
        ip.coord[2] = 0.0;
        ip.nsurf = 0;
        // include what happens if it lies on a node or surfaces
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
                                                    Epetra_SerialDenseVector&   E1,
                                                    Epetra_SerialDenseVector&   E2)
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
	
    if(surfaceElement == NULL )     dserror("pointer to surface element is NULL");  
    if(lineElement == NULL)         dserror("pointer to line element is NULL");
    
	updateRHSForCSI( b, xsi, surfaceElement, lineElement);
             					
	while(residual > TOL)
   	{   
   		updateAForCSI( A, xsi, surfaceElement, lineElement);
  		if(!gaussElimination(A, b, dx, true, 3))
        {
            intersection = false;
            break;  
        } 
     	xsi = addTwoVectors(xsi,dx);
      
    	if( (fabs(xsi[0]) > E1[0]) || (fabs(xsi[1]) > E1[1]) || (fabs(xsi[2]) > E1[2])) 
    	{
      		intersection = false;
      		break;
		}
        
        if(iter >= maxiter)
        {   
            intersection = false;
            break;
        }       
       
		updateRHSForCSI( b, xsi, surfaceElement, lineElement);
     	residual = b.Norm2();      
     	iter++;
    } 
    return intersection;
}



/*----------------------------------------------------------------------*
 |  collects all intersection points of a line and           u.may 06/07|
 |  a surface                                                           |
 *----------------------------------------------------------------------*/  
void Intersection::collectInterfacePoints(  DRT::Element*                   surfaceElement,
				 							DRT::Element*                   lineElement,
				 							std::vector<InterfacePoint>&    interfacePoints,
                                            int&                            numInterfacePoints,
                                            int                             numSurface,
                                            int                             numLine,
                                            bool                            lines)
{

	bool intersected = false;
    const double TOL = EPS8;
	Epetra_SerialDenseVector xsi(3);
	Epetra_SerialDenseVector E1(3);
   	Epetra_SerialDenseVector E2(3);
    InterfacePoint ip;
	
	xsi[0] = 0.0; xsi[1] = 0.0; xsi[2] = 0.0;
	E1[0] = 1.0 + TOL;  E1[1] = 1.0 + TOL ;  E1[2] = 1.0 + TOL;
	E2[0] = -1.0 - TOL; E2[1] = -1.0 - TOL;  E2[2] = -1.0 - TOL;

	intersected = computeCurveSurfaceIntersection( surfaceElement, lineElement, xsi, E1, E2);
             							 	
	if(intersected)		
   	{	    
   		// include that intersectionpoint lies on a node or surface
        if(lines)
        {   
            ip.nsurf = 1;
            ip.surfaces[0] = numSurface;
             
            // change minus sign if you change the line numbering   
            switch(numLine)
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
                    dserror("node number not correct");    // changed
            }           
        }
        else
        {
            ip.nsurf = 2;
            // change for all distypes
            ip.surfaces[0] = eleNodeNumbering_hex27_lines_surfaces[numLine][0];
            ip.surfaces[1] = eleNodeNumbering_hex27_lines_surfaces[numLine][1];
            ip.coord[0] = xsi[0]; 
            ip.coord[1] = xsi[1]; 
        }
        
        ip.coord[2] = 0.0; 
        interfacePoints.push_back(ip);  
   		numInterfacePoints++; //+
  // 		computeNewStartingPoint(surfaceElement, lineElement, E1, xsi, interfacePoints) +
  // 		computeNewStartingPoint(surfaceElement, lineElement, xsi, E2, interfacePoints);
   	}      							 	
}



/*----------------------------------------------------------------------*
 |  computes a new statring point for the Newton-method      u.may 06/07|
 |  in order to find all intersection points                            |
 |  of a curve-surface intersection                                     |
 *----------------------------------------------------------------------*/  
int Intersection::computeNewStartingPoint(	DRT::Element*                surfaceElement,
											DRT::Element*                lineElement,
											Epetra_SerialDenseVector     E1,
     										Epetra_SerialDenseVector     E2,
     										std::vector<InterfacePoint>& interfacePoints)
{	
    // to be tested
	bool intersected = false;
	int numInterfacePoints = 0;
	Epetra_SerialDenseVector xsi(3);
	
	for(int i = 0; i<3; i++)
		xsi[i] = (double) (( E1[i] + E2[i] )/2.0);
		
	// intersected = computeCurveSurfaceIntersection(surfaceElement, lineElement, xsi, E1, E2);
             							 	
	if(intersected)		
   	{
   		printf("hello2.0\n");
        //Intersection
   		//interfacePoints.push_back(xsi);        
   		//numInterfacePoints = 1 +
   		//computeNewStartingPoint(surfaceElement, lineElement, E1, xsi, interfacePoints) +
   		//computeNewStartingPoint(surfaceElement, lineElement, xsi, E2, interfacePoints);
   		printf("hello2.1\n");
   	}      
   
   	printf("number of intersection points = %d\n", numInterfacePoints );
   	return numInterfacePoints;    						
}



/*----------------------------------------------------------------------*
 |  transforms a elememt in reference coordinates            u.may 06/07|
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
        xsi[dim] = 0.0;
        for(int i=0; i<numNodes; i++)
            xsi(dim) += element->Nodes()[i]->X()[dim] * funct(i);
    }
}


/*----------------------------------------------------------------------*
 |  compares two nodes                                       u.may 06/07|
 |  overloaded method:  vector<double>  and double*                     |
 *----------------------------------------------------------------------*/  
bool Intersection::compareNodes(    vector<double>&     node1,
                                    double*             node2)
{   
    bool equal = true;
    
    if(node2 == NULL)
        dserror("array is NULL");
        
    for(unsigned int i = 0; i < node1.size() ; i++)
        if(fabs(node1[i] - node2[i]) > EPS8)
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
bool Intersection::compareNodes(    vector<double>& node1,
                                    vector<double>& node2)
{   
    bool equal = true;
    
    if(node1.size() != node2.size())
        dserror("arrays of nodes need to have the same length");
             
    for(unsigned int i = 0; i < node1.size() ; i++)
        if(fabs(node1[i] - node2[i]) > EPS8)
        {
            equal = false;
            break;
        }
  
    return equal;
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
        if(compareNodes(searchPoint, *it))
        {
            pointfound = true;
            searchPoint = *(it+1);              
            vertices.erase(it);
            vertices.erase(it); // remove it+ 1
            break; 
        }
               
        if(compareNodes(searchPoint, *(it+1)))
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
    vector<InterfacePoint>::iterator it2;
        
    for(unsigned int i = 0; i < interfacePoints.size(); i++ )
        if(compareNodes(point, interfacePoints[i].coord))
        {         
            alreadyInList = false;
            count = 0;
            for(it2 = pointList.begin(); it2 != pointList.end(); it2++ )  
            {
                if(compareNodes(point, it2->coord))   
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



/*----------------------------------------------------------------------*
 |  stores a segment within a list of segments               u.may 06/07|
 |  which is to be copy to the tetgen data structure                    |
 |  for the computation of the Constrained Delauney Triangulation       |
 *----------------------------------------------------------------------*/  
void Intersection::storeSegments(   vector<InterfacePoint>&     pointList, 
                                    vector<int>                 positions, 
                                    vector< vector<int> >&      segmentList)
{
    int pos1 = 0;
    int pos2 = 0;
    
    // change + 8 number of corner points
    for(unsigned int i = 0; i < positions.size()-2; i++ )
    {
        pos1 = positions[i];
        pos2 = positions[i+1];
        
        for(int j = 0; j < pointList[pos1].nsurf; j++ )  
            for(int k = 0; k < pointList[pos2].nsurf; k++ ) 
                if(pointList[pos1].surfaces[j] == pointList[pos2].surfaces[k]) 
                { 
                    segmentList[pointList[pos1].surfaces[j]].push_back(pos1+8);
                    segmentList[pointList[pos1].surfaces[j]].push_back(pos2+8);
                }
    }
      
    pos1 = positions[positions.size()-2];
    pos2 = positions[0]; 
    
    for(int j = 0; j < pointList[pos1].nsurf; j++ )  
        for(int k = 0; k < pointList[pos2].nsurf; k++ ) 
            if(pointList[pos1].surfaces[j] == pointList[pos2].surfaces[k]) 
            { 
                segmentList[pointList[pos1].surfaces[j]].push_back(pos1+8);
                segmentList[pointList[pos1].surfaces[j]].push_back(pos2+8);
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
    
    for(unsigned int i = 0; i < positions.size()-2; i++ )
    {
        triangle[0] = positions[i]+8;
        triangle[1] = positions[i+1]+8;
        triangle[2] = pointList.size()-1+8;
        
        triangleList.push_back(triangle);
    }
    
    triangle[0] = positions[positions.size()-2]+8;
    triangle[1] = positions[0]+8;
    triangle[2] = pointList.size()-1+8;
        
    triangleList.push_back(triangle);
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
 |  computes the convex hull of a set of interface points    u.may 06/07|
 |  and stores resulting points, segments and triangles                 |
 |  for the use with Tetgen (CDT)                                       |
 *----------------------------------------------------------------------*/  
void Intersection::computeConvexHull(   DRT::Element*           surfaceElement,
                                        vector<InterfacePoint>& interfacePoints,
                                        vector<InterfacePoint>& pointList,
                                        vector< vector<int> >&  segmentList,
                                        vector< vector<int> >&  triangleList )
{
    double* point;
    vector<int> positions;
    vector<double> searchPoint(3,0);
    vector<double> vertex(3,0);
    vector< vector<double> > vertices;  
    Epetra_SerialDenseVector curCoord(3);   
    
       
    // compute midpoint in reference coordinates
    InterfacePoint midpoint = computeMidpoint(interfacePoints);
    // transform it into current coordinates
    for(int j = 0; j < 2; j++)      curCoord[j]  = midpoint.coord[j];            
    referenceToCurrentCoordinates(surfaceElement, curCoord);        
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
        for(int j = 0; j < 2; j++)      curCoord[j]  = interfacePoints[i].coord[j];           
        referenceToCurrentCoordinates(surfaceElement, curCoord);       
        for(int j = 0; j < 3; j++)      interfacePoints[i].coord[j] = curCoord[j];
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
            {
                vertex[k] = point[k];
            }      
     
            for(int m = 0; m < 2; m++)      curCoord[m]  = vertex[m];           
            referenceToCurrentCoordinates(surfaceElement, curCoord);       
            for(int m = 0; m < 3; m++)      
            {
                vertex[m] = curCoord[m];
                //printf("vertex2 = %f   ", vertex[m]); 
            }    
               
            
           
            vertices.push_back(vertex);
        }
        facet = facet->next;
    }
    
    // store points, segments and triangles for the computation of the
    // Constrained Delaunay Tetrahedralization with Tetgen  
    searchPoint = vertices[1];  
    storePoint(vertices[0], interfacePoints, positions, pointList );
    storePoint(vertices[1], interfacePoints, positions, pointList );
    vertices.erase(vertices.begin());
    vertices.erase(vertices.begin());
    
    while(!vertices.empty())
    {                    
        findNextFacet(vertices, searchPoint);
        storePoint(searchPoint, interfacePoints, positions, pointList);
    } 
   
    storeSegments(pointList, positions, segmentList);
    pointList.push_back(midpoint);
    storeTriangles(pointList, positions, triangleList);
   
    // free memory and clear vector of interface points
    qh_freeqhull(!qh_ALL);
    int curlong, totlong;           // memory remaining after qh_memfreeshort 
    qh_memfreeshort (&curlong, &totlong);
    if (curlong || totlong) 
        printf("qhull internal warning (main): did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);
    
    free(coordinates);
    interfacePoints.clear();
}



/*----------------------------------------------------------------------*
 |  computes the Constrained Delaunay Tetrahedralization     u.may 06/07|
 |  in 3D with help of Tetgen library for an intersected                |
 |  xfem element                                                        |
 *----------------------------------------------------------------------*/  
void Intersection::computeCDT(  DRT::Element*               element,
                                vector< InterfacePoint >&   pointList,
                                vector< vector<int> >&      segmentList,
                                vector< vector<int> >&      triangleList,
                                vector < vector <Integrationcell> >& integrationcellList)
{
    int dim = 3;
    int cornerpoints = 8;
    int nsegments = 0; 
    tetgenio in;
    tetgenio out;
    char switches[] = "po2";
    tetgenio::facet *f;
    tetgenio::polygon *p;
        
    
    // set points
    in.numberofpoints = cornerpoints + pointList.size();
    in.pointlist = new REAL[in.numberofpoints * dim];
    
    
    // fill point list
    int fill = 0;
    for(int i = 0; i < cornerpoints; i++)
        for(int j = 0; j < dim; j++)  
            in.pointlist[fill++] = element->Nodes()[i]->X()[j]; 
            
    for(int i = cornerpoints; i < in.numberofpoints; i++)
        for(int j = 0; j < dim; j++)  
            in.pointlist[fill++] = pointList[i-cornerpoints].coord[j];
       
       
    // set facets
    if(triangleList.size()>0)       in.numberoffacets = 6 + triangleList.size(); 
    else                            in.numberoffacets = 6;   
      
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    
    
    // loop over all xfem element surfaces
    for(int i = 0; i < element->NumSurface(); i++)
    {
        f = &in.facetlist[i];
        if(segmentList[i].size() > 0)       nsegments = (int) (segmentList[i].size()/2);
        else                                nsegments = 0;
        f->numberofpolygons = 1 + nsegments; 
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];
        p->numberofvertices = 4;
        p->vertexlist = new int[p->numberofvertices];
        for(int j = 0; j < 4; j ++)
        {
             p->vertexlist[j] = eleNodeNumbering_hex27_surfaces[i][j];
             //printf("surfaces = %d\t", eleNodeNumbering_hex27_surfaces[i][j]);
        }  
        //printf("\n");  
        int count = 0;
        for(int j = 1; j < f->numberofpolygons; j ++)
        {
            if(segmentList[i].size() > 0)
            {             
                p = &f->polygonlist[j];
                p->numberofvertices = 2;
                p->vertexlist = new int[p->numberofvertices];
            
                for(int k = 0; k < 2; k ++)
                {
                   p->vertexlist[k] = segmentList[i][count];
                   //printf("segment = %d\t",segmentList[i][count] );
                   count++;
                }
                //printf("\n");
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
        {
             p->vertexlist[j] = triangleList[i - element->NumSurface()][j];
             //printf("triangle = %d\t",triangleList[i - element->NumSurface()][j] );
        }
        //printf("\n");
    }
        
        
    // set facetmarkers
    for(int i = 0; i < in.numberoffacets; i ++)
        in.facetmarkerlist[i] = 0;   

    in.save_nodes("tetgen");
    in.save_poly("tetgen");
    
    //  Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
    //  do quality mesh generation (q) with a specified quality bound
    //  (1.414), and apply a maximum volume constraint (a0.1)
    tetrahedralize(switches, &in, &out, NULL, NULL); 
    
    out.save_elements("tetgenout");
    
    vector<double> tetnodes(3);
    vector< vector<double> > tetrahedronCoord;
    vector< Integrationcell > listperElement;
    
    for(int i=0; i<out.numberoftetrahedra; i++ )
    {   
        for(int j = 0; j < out.numberofcorners; j++)
        {
            for(int dim = 0; dim < 3; dim++)
                tetnodes[dim] = out.pointlist[out.tetrahedronlist[i*10+j]*3+dim];
         
            tetrahedronCoord.push_back(tetnodes);    
        }
        Integrationcell cell(i, tetrahedronCoord);
        listperElement.push_back(cell);                 
    }
    integrationcellList.push_back(listperElement);
    
    
    // clear vectors
    pointList.clear();
    for(unsigned int i  = 0; i < segmentList.size(); i++)
        segmentList[i].clear();
    triangleList.clear();
    listperElement.clear();
}



/*----------------------------------------------------------------------*
 |  computes the intersection between two different          u.may 06/07|
 |  discretizations and returns a list of intersected xfem elements     |
 |  and their integrations cell                                         |
 *----------------------------------------------------------------------*/  
void Intersection::computeIntersection( RefCountPtr<DRT::Discretization> actdis, 
                                        vector < vector <Integrationcell> >& integrationcellList)
{
    bool intersected =                      false;
    int numInterfacePoints =                0;
    vector< DRT::Condition * >              xfemConditions;
    vector< InterfacePoint >                interfacePoints;
    map< int, RefCountPtr<DRT::Element > >  geometryMap;
    DRT::Element*                           cutterElement; 
    DRT::Element*                           xfemElement; 
    Epetra_SerialDenseMatrix                cutterAABB;
    Epetra_SerialDenseMatrix                xfemAABB;
    vector< InterfacePoint >                pointList;
    vector< vector<int> >                   segmentList(6);
    vector <vector<int> >                   triangleList;
    
    
        
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
        
        //xfemConditions.size()
        for(unsigned int i=0; i<xfemConditions.size(); i++)
        {
            geometryMap = xfemConditions[i]->Geometry();
            if(geometryMap.size()==0)   dserror("geometry does not obtain elements");
            // printf("size of %d.geometry map = %d\n",i, geometryMap.size());
            
            //geometryMap.size()
            for(unsigned int j=0; j<geometryMap.size(); j++)
            {
                numInterfacePoints = 0;  
                cutterElement = geometryMap.find(j)->second.get();
                        
                if(cutterElement == NULL) dserror("geometry does not obtain elements");
            
                cutterAABB = computeFastAABB(cutterElement);                                 
                intersected = intersectionOfAABB(cutterAABB, xfemAABB);    
                //debugAABBIntersection( cutterAABB, xfemAABB, cutterElement, xfemElement, i, k);
                                     
                if(intersected)
                {
                    for(int m=0; m<cutterElement->NumLine() ; m++)
                    {                    
                        checkNodeWithinElement( xfemElement, cutterElement->Nodes()[m],
                                                interfacePoints, k, m, numInterfacePoints);
                                                                
                        for(int p=0; p<xfemElement->NumSurface() ; p++)       
                            collectInterfacePoints( xfemElement->Surfaces()[p], cutterElement->Lines()[m],
                                                    interfacePoints, numInterfacePoints, p, m, true);  
                    }
                                      
                    for(int m=0; m<xfemElement->NumLine() ; m++) 
                        collectInterfacePoints( cutterElement, xfemElement->Lines()[m],
                                                interfacePoints, numInterfacePoints, 0, m, false);                                         
                                       
                    if(numInterfacePoints!= 0)
                        computeConvexHull(  cutterElement, interfacePoints,  pointList,
                                            segmentList, triangleList);                      
                }// if intersected
            }// for-loop over all geometryMap.size()
        }// for-loop over all xfemConditions.size() 
        //debugTetgenDataStructure(pointList, segmentList, triangleList);
        computeCDT(xfemElement, pointList, segmentList, triangleList, integrationcellList);
    }// for-loop over all  actdis->NumMyRowElements()
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
void Intersection::debugTetgenDataStructure(    vector< InterfacePoint >&   pointList,
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
    for(unsigned int i = 0; i< pointList.size(); i++)
    {
        cout << i+8 << ".th point:   ";
        for(int j = 0; j< 3; j++)
        {
            cout << pointList[i].coord[j] << "\t";
        }
        cout << endl;
        cout << endl;
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
void Intersection::debugIntegrationcells(vector < vector <Integrationcell> >& integrationcellList)
{    
    cout << endl;
    cout << "===============================================================" << endl;
    cout << "Debug RESULTING INTEGRATION CELLS " << endl;
    cout << "===============================================================" << endl;
    cout << endl;
    
    for(unsigned int i = 0; i < integrationcellList.size(); i++ )
	{
		cout << "XFEM ELEMENT " << i << " :" << endl;	
		cout << endl;
		cout << "NUMBER OF INTEGRATIONCELLS : " << integrationcellList[i].size() << endl;	
		cout << endl;
		for(unsigned int j = 0; j < integrationcellList[i].size(); j++ )
		{
			cout << "IC " << j << ":  " << endl;
			cout << endl;	
			for(unsigned int k = 0; k < integrationcellList[i][j].GetCoord().size(); k++)
			{
				for(unsigned int m = 0; m < integrationcellList[i][j].GetCoord()[k].size(); m++)
				{
					cout << integrationcellList[i][j].GetCoord()[k][m] << "   ";	
				}
				cout << endl;	
			}
			cout << endl;	
		}
		cout << endl;	cout << endl;	
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
#endif  // #ifdef XFEM



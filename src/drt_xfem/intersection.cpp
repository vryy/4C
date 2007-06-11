/*!----------------------------------------------------------------------
\file drt_intersection.cpp

\brief collection of intersection tool 

<pre>
Maintainer: Ursula Mayer
</pre>

*----------------------------------------------------------------------*/


#ifdef D_FLUID3_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "intersection.h"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_element.H"
#include "../drt_f3_xfem/fluid3_xfem.H"
#include "Epetra_SerialDenseVector.h"




/*
 * adds two vectors v1 + v2
 */
 
Epetra_SerialDenseVector Intersection::addTwoVectors(	Epetra_SerialDenseVector& v1, 
																		Epetra_SerialDenseVector& v2)
{   
   if(v1.Length() != v2.Length())
       dserror("both vectors need to have the same size\n"); 

   for(int i = 0; i < v1.Length(); i++)
     v1[i] = v1[i] + v2[i];
 
   return v1;
}
	

/*
 * subtracts two vectors v1 - v2
 */
std::vector<double> Intersection::subtractsTwoVectors(	std::vector <double>& v1, 
																			std::vector <double>& v2)
{   
   if(v1.size() != v2.size())
      dserror("both vectors need to have the same size\n"); 

   for(unsigned int i = 0; i < v1.size(); i++)
     v1[i] = v1[i] - v2[i];
 
   return v1;
}


/*
 * computes a rough axis-aligned bounding box
 */
Epetra_SerialDenseMatrix Intersection::computeFastAABB(DRT::Element* element)
{
   double  maxdistance; 
   double tol = 1e-10;
   // minX, minY, minZ, 
   // maxX, maxY, maxZ;
   Epetra_SerialDenseMatrix AABB(3, 2);
   
      
   if(element->Nodes()==NULL)
   	dserror("element does not contain nodes \n");
    
	for(int dim=0; dim<3; dim++)
	{
		AABB(dim, 0) = element->Nodes()[0]->X()[dim]-tol;
   	AABB(dim, 1) = element->Nodes()[0]->X()[dim]+tol;
	}
    
   for(int i=1; i<element->NumNode(); i++)
   	for(int dim=0; dim<3; dim++)
		{
			AABB(dim, 0) = std::min( AABB(dim, 0), element->Nodes()[i]->X()[dim]-tol);
			AABB(dim, 1) = std::max( AABB(dim, 1), element->Nodes()[i]->X()[dim]+tol);
		}
 
 	maxdistance = fabs(AABB(0,1) - AABB(0,0));
 	for(int dim=1; dim<3; dim++)
		maxdistance = std::max(maxdistance, fabs(AABB(dim,1)-AABB(dim,0)) );
	
	for(int dim=0; dim<3; dim++)
	{
		AABB(dim, 0) = AABB(dim, 0) - 0.5*maxdistance;
		AABB(dim, 1) = AABB(dim, 1) + 0.5*maxdistance;
	}	
	
	printf("\n");
	printf("axis-aligned bounding box:\n minX = %f\n minY = %f\n minZ = %f\n maxX = %f\n maxY = %f\n maxZ = %f\n", 
			  AABB(0,0), AABB(1,0), AABB(2,0), AABB(0,1), AABB(1,1), AABB(2,1));
	printf("\n");
	
	return AABB;
}



/*
 *  checks is a node is within a axis aligned bounding box
 */
bool Intersection::isNodeWithinAABB( 	std::vector<double> node, 
													Epetra_SerialDenseMatrix AABB)
{
	bool isWithin = true;
	
	for (int dim=0; dim<3; dim++)
	{
   	if((node[dim]<AABB(dim,0))||(node[dim]>AABB(dim,1)))
     		isWithin = false;
   }
		 	
	return isWithin;
}


/*
 *  checks if two axis-aligned bounding boxes intersect
 */
bool Intersection::intersectionOfAABB(	Epetra_SerialDenseMatrix cutterAABB, 
													Epetra_SerialDenseMatrix xfemAABB)
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



/*
 * computes a Gaussian Elimination for a linear system of equations
 */
void Intersection::gauss_elimination( 	Epetra_SerialDenseMatrix& A,
                         					Epetra_SerialDenseVector& b,
                						  		Epetra_SerialDenseVector& x,
                								bool do_piv,
                								const int dim)
{
   int i,j,k;
   double tmp[4];    // check tmp and size of A with dim
   int pivot;

    //printf("A[0][] = %8e %8e %8e | %8e\n", A(0,0),A(0,1),A(0,2),b[0]);
    //printf("A[1][] = %8e %8e %8e | %8e\n", A(1,0),A(1,1),A(1,2),b[1]);
    //printf("A[2][] = %8e %8e %8e | %8e\n", A(2,0),A(2,1),A(2,2),b[2]);

   if (!do_piv) 
   {
      for (k=0;k<dim;k++)
      {            
         A(k,k)=1./A(k,k);

         for (i=k+1;i<dim;i++)
         {
            A(i,k) = A(i,k) * A(k,k);
            x[i] = A(i,k);

            for (int j=k+1;j<dim;j++)
            {
               A(i,j) = A(i,j) - A(i,k) * A(k,j);
            }
         }

         for (i=k+1;i<dim;i++)
         {
            b[i]=b[i]-x[i]*b[k];
         }
      }
   }
   else 
   {
      for (k=0;k<dim;k++)
      {
         pivot = k;
         // search for pivot element 
         for (i=k+1;i<dim;i++)
         {
            pivot = (fabs(A(pivot,pivot)) < fabs(A(i,k))) ? i : pivot;
         }
         // copy pivot row to current row 
         if (pivot != k) 
         {
            for (j=0;j<dim;j++)
               tmp[j] = A(pivot,j);

            tmp[dim] = b[pivot];

            for (j=0;j<dim;j++)
               A(pivot,j) = A(k,j);

            b[pivot] = b[k];

            for (j=0;j<dim;j++)
               A(k,j) = tmp[j];

            b[k] = tmp[dim];
         }

         A(k,k) = 1./A(k,k);
         //printf("inf_diag = %8e\n", A(k,k));
         //fflush(NULL);

         for (i=k+1;i<dim;i++)
         {
            A(i,k) = A(i,k) * A(k,k);
            x[i] = A(i,k);

            for (j=k+1;j<dim;j++)
            {
               A(i,j) = A(i,j) - A(i,k) * A(k,j);
            }
         }

         for (i=k+1;i<dim;i++)
            b[i]=b[i]-x[i]*b[k];
         
         //printf("A[0][] = %8e %8e %8e | %8e\n", A(0,0),A(0,1),A(0,2),b[0]);
         //printf("A[1][] = %8e %8e %8e | %8e\n", A(1,0),A(1,1),A(1,2),b[1]);
         //printf("A[2][] = %8e %8e %8e | %8e\n", A(2,0),A(2,1),A(2,2),b[2]);
      }
   }

   // backward substitution 
   x[dim-1]=b[dim-1]*A(dim-1,dim-1);

   for (i=dim-2;i>=0;i--)
   {
      for (j=dim-1;j>i;j--)
      {
         b[i]=b[i]-A(i,j)*x[j];
      }
      x[i]=b[i]*A(i,i);
   }

    //for (i=0;i<dim;i++)
    //    printf("%8e ",x[i]);
    //printf("\n");
}



void Intersection::updateForCurveSurfaceIntersection( Epetra_SerialDenseMatrix& A,
				 														Epetra_SerialDenseVector& b,
             														Epetra_SerialDenseVector& xsi,
             														DRT::Element* surfaceElement,
				 														DRT::Element* lineElement)        											
{
	int numOfNodesSurface = surfaceElement->NumNode();
   int numOfNodesLine = lineElement->NumNode();
	vector<int> actParams(1,0);
	Epetra_SerialDenseMatrix surfaceDeriv1(2,numOfNodesSurface);
	Epetra_SerialDenseMatrix lineDeriv1(1,numOfNodesLine);
	Epetra_SerialDenseVector surfaceFunct(numOfNodesSurface);
	Epetra_SerialDenseVector lineFunct(numOfNodesLine);
	Epetra_SerialDenseMatrix emptyM;
	Epetra_SerialDenseVector emptyV;
	DRT::Discretization dummyDis("dummy discretization", null);
	ParameterList params;
	
	params.set("action","calc_ShapefunctDeriv1Deriv2");
	actParams[0] = numOfNodesSurface;     
	
	surfaceElement->Evaluate(params, dummyDis, actParams, surfaceDeriv1, emptyM , surfaceFunct, xsi, emptyV);		
	
   for(int dim=0; dim<3; dim++)
   	for(int i=0; i<numOfNodesSurface; i++)
		{
			A(dim, 0) += surfaceElement->Nodes()[i]->X()[dim] * surfaceDeriv1(0,i);
			A(dim, 1) += surfaceElement->Nodes()[i]->X()[dim] * surfaceDeriv1(1,i);
		}
		
	
	actParams[0] = numOfNodesLine;     
	lineElement->Evaluate(params, dummyDis, actParams, lineDeriv1, emptyM, lineFunct, xsi, emptyV);
	
	for(int dim=0; dim<3; dim++)
   	for(int i=0; i<numOfNodesLine; i++)
		{
			A(dim, 2) +=  (-1) * lineElement->Nodes()[i]->X()[dim] * lineDeriv1(0,i);
		}
	
	
   for(int dim=0; dim<3; dim++)
   	for(int i=0; i<numOfNodesSurface; i++)
		{
			b(dim) += (-1) * surfaceElement->Nodes()[i]->X()[dim] * surfaceFunct(i);
		}
	
	for(int dim=0; dim<3; dim++)
   	for(int i=0; i<numOfNodesLine; i++)
		{
			b(dim) += lineElement->Nodes()[i]->X()[dim] * lineFunct(i);
		}
}



void Intersection::updateForPointWithinElement( Epetra_SerialDenseMatrix& A,
													 			Epetra_SerialDenseVector& b,
													 			Epetra_SerialDenseVector& xsi,
													 			Epetra_SerialDenseVector& x,
													 			DRT::Element* element)         											
{
	int numnodes = element->NumNode();
	vector<int> actParams(1,0);
	Epetra_SerialDenseMatrix deriv1(3, numnodes);
	Epetra_SerialDenseVector funct(numNodes);
	Epetra_SerialDenseMatrix emptyM;
	Epetra_SerialDenseVector emptyV;
	DRT::Discretization dummyDis("dummy discretization", null);
	ParameterList params;
	
	
	params.set("action","calc_ShapefunctDeriv1Deriv2");
	actParams[0] = numnodes;     
	
	element->Evaluate(params, dummyDis, actParams, deriv1, emptyM , funct, xsi, emptyV);		
	
   for(int dim=0; dim<3; dim++)
   	for(int i=0; i<numnodes; i++)
		{
			A(dim, 0) += element->Nodes()[i]->X()[dim] * deriv1(0,i);
			A(dim, 1) += element->Nodes()[i]->X()[dim] * deriv1(1,i);
			A(dim, 2) += element->Nodes()[i]->X()[dim] * deriv1(2,i);
		}
	
	
   for(int dim=0; dim<3; dim++)
   	for(int i=0; i<numNodes; i++)
		{
			b(dim) += (-1) * element->Nodes()[i]->X()[dim] * funct(i);
		}
		
	b = addTwoVectors(b,x);
}



bool Intersection::computeCurveSurfaceIntersection( 	DRT::Element* surfaceElement,
				 														DRT::Element* lineElement,
             														Epetra_SerialDenseVector& xsi,
             														Epetra_SerialDenseVector& E1,
             														Epetra_SerialDenseVector& E2)
{
	int iter = 0;
	int maxiter = 500;
	bool intersection = true;
   double residual = 1.0;
   Epetra_SerialDenseMatrix A(3,3);
   Epetra_SerialDenseVector b(3);
   Epetra_SerialDenseVector dx(3);
   dx = xsi;
	
	updateForCurveSurfaceIntersection( A, b, xsi, surfaceElement, lineElement);
	
	while(residual > 1e-8 && iter > maxiter)
   {
      
      gauss_elimination(A, b, dx, false, 3);       
      xsi = addTwoVectors(xsi,dx);
          
      if( (fabs(xsi[0]) > E1[0]) || (fabs(xsi[1]) > E1[1]) || (fabs(xsi[2]) > E1[2]) ||
      	 (fabs(xsi[0]) > E2[0]) || (fabs(xsi[1]) > E2[1]) || (fabs(xsi[2]) > E2[2]) )   
      {
      	intersection = false;
      	break;
     	}

      updateForCurveSurfaceIntersection( A, b, xsi, surfaceElement, lineElement);
      residual = b.Norm2();
      iter++;
   }
   return intersection;
}


bool Intersection::checkPointWithinElement( DRT::Element* element,
				 										  Epetra_SerialDenseVector& xsi,
             										  Epetra_SerialDenseVector& x)
{
	bool pointWithinElement = true;
   double residual = 1.0;
   Epetra_SerialDenseMatrix A(3,3);
   Epetra_SerialDenseVector b(3);
   Epetra_SerialDenseVector dx(3);
   dx = xsi;
	
	updateForPointWithinElement( A, b, xsi, x, element);

	while(residual > 1e-8)
   { 
      gauss_elimination(A, b, dx, false, 3);       
      xsi = addTwoVectors(xsi,dx);
      
      if( (fabs(xsi[0]) > 1.0) || (fabs(xsi[1]) > 1.0) || (fabs(xsi[2]) > 1.0) )     // eve + TOL
      {
      	pointWithinElement = false;
      	break;
     	}

      updateForPointWithinElement( A, b, xsi, x, element);
      residual = b.Norm2();
   }
   return pointWithinElement;
}


   
int Intersection::collectIntersectionPoints(	DRT::Element* surfaceElement,
				 											DRT::Element* lineElement,
				 											std::set< Epetra_SerialDenseVector& >& intersectionPointSet)
{

	bool intersected;
	int numIntersectionPoints = 0;
	Epetra_SerialDenseVector xsi(3);
	Epetra_SerialDenseVector E1(3);
   Epetra_SerialDenseVector E2(3);
	
	xsi[0] = 0.0; xsi[1] = 0.0; xsi[2] = 0.0;
	E1[0] = 1.0;  E1[1] = 1.0;  E1[2] = 1.0;
	E2[0] = -1.0; E2[1] = -1.0; E2[2] = -1.0;


	intersected = computeCurveSurfaceIntersection( surfaceElement, lineElement, xsi, E1, E2);
             							 	
	if(intersected == true)		
   {
   	// IntersectionPointSet.add(xsi);
   	computeNewStartingPoint(surfaceElement, lineElement, E1, xsi, intersectionPointSet);
   	computeNewStartingPoint(surfaceElement, lineElement, xsi, E2, intersectionPointSet);
   }          							 	
             							 	
	return numIntersectionPoints;         							 	
}



void computeNewStartingPoint(	DRT::Element* surfaceElement,
										DRT::Element* lineElement,
										Epetra_SerialDenseVector& E1,
     									Epetra_SerialDenseVector& E2,
     									std::set< Epetra_SerialDenseVector& >& intersectionPointSet)
{	
	
	bool intersected;
	Epetra_SerialDenseVector xsi(3);
	
	for(int i = 0; i<3; i++)
		xsi[i] = ( E1[i] + E2[i] )/2;
		
	intersected = computeCurveSurfaceIntersection( emptyDis, surfaceElement, lineElement, xsi,
             							 					  E1, E2, true, 3);
             							 	
	if(intersected == true)		
   {
   	// add xsi to intersection list IntersectionPointSet.add(xsi);
   	computeNewStartingPoint(surfaceElement, lineElement, E1, xsi, intersectionPointSet);
   	computeNewStartingPoint(surfaceElement, lineElement, xsi, E2, intersectionPointSet);
   }          						
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3_XFEM



/*

get geometry (DRT::Condition) of solid elements

get Discretization of fluid elements

loop over all surface solid elements

  compute the boundingbox of the solid element 
  // read wihich type of bounding boxes (Jean Hamman, Huber)

  loop over all fluid elemenets // replace by a tree search later
    compute the bounding of the fluid element
    if the bounding boxes intersect

      fluid element stores pointer to solid   // map of vectors possible 
      a pointer to the fluid element is stored in a list of possible xfem fluid elements
    end
  end
end             
        


loop over all possible xfem fluid elements

  loop over all intersecting surface elements:
    loop over all lines of the surface element
         check if endpoints are in or outside the fluid element
            store the points inside the fluid element in parameter coordinates

      loop over all fluid surfaces
                compute all curve-surface intersections in the surface parent domain
                use modified Newton for ill-conditioned nonlinear systems
                (if surfaces have different parametrization do an normalization of the tangent vectors)

                check if intersection points are within the element domain
                    store them in parameter coordinates in a list
                    compute a convex hull
                    check if lines are in one plane
                        if not 
                            add new points
                    store the lines in the linearized discretization or set first

                create an xfem element
                do CDT
                recover quadratic lines,k compute the midpoint of each line and relax the midpoint back onto the element
            end
end


*/

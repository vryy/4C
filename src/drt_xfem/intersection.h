/*!----------------------------------------------------------------------
\file drt_intersection.H

\brief collection of intersection tools

<pre>
Maintainer: Ursula Mayer
</pre>

*----------------------------------------------------------------------*/


#ifdef D_FLUID3_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef INTERSECTION_H
#define INTERSECTION_H

#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Vector.h"


#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_f3_xfem/fluid3_xfem.H"

#include <vector>
#include <set>
#include <cmath>
#include <algorithm>


using namespace std;
using namespace Teuchos;


class Intersection
{
	
   public:	
		
		/*!
	  	\brief Returns the sum of two double vectors
	
	  	\param v1 (in) : arbitrary double vector
	  	\param v2 (in) : arbitrary double vector
	  	\return sum of two double vectors
	  	*/
		Epetra_SerialDenseVector addTwoVectors(	Epetra_SerialDenseVector& v1,
	                           						Epetra_SerialDenseVector& v2);
		
		
		/*!
	  	\brief Returns the difference of two double vectors v1 - v2
	  	\param v1 (in) : arbitrary double vector		
	  	\param v2 (in) : arbitrary double vector
	  	\return difference of two double vectors v1 - v2
	  	*/
		std::vector<double> subtractsTwoVectors(	std::vector <double>& v1, 
																std::vector <double>& v2);
																	
	  	/*!
	  	\brief computes a rough overestimating bounding box for an element
	
	   \param element 		 (in) 	   : element (DRT::Element*)
		\return bounding box for an element
	  	*/
		Epetra_SerialDenseMatrix computeFastAABB(	DRT::Element* element);
		                      
	                           
		/*!
	  	\brief checks is a node is within a axis aligned bounding box
	
	  	\param node (in) 		: node to checked
	  	\param AABB (in) 		: axis aligned bouning box
	  	\return isWithin = true if node is within the AABB or false otherwise
	  	*/
		bool isNodeWithinAABB(	std::vector<double> node, 
										Epetra_SerialDenseMatrix AABB);
										
													
		/*!
	  	\brief checks if two axis aligned bounding boxes intersect
	
	  	\param cutterAABB (in) 		: node to checked
	  	\param xfemAABB   (in) 		: axis aligned bouning box
	  	\return intersection = true if the AABB's intersect or false otherwise
	  	*/
		bool intersectionOfAABB(	Epetra_SerialDenseMatrix cutterAABB, 
											Epetra_SerialDenseMatrix xfemAABB);

																		
		/*!
	  	\brief computes a Gaussian elimination of linear system of equations
	
	  	\param A (in) 		: system matrix
	  	\param b (in) 		: right-hand-side
	  	\param x (in) 		: solution vector
	  	\param do_piv (in)	: do_piv = true does pivoting, do_piv = false does not do pivoting
	  	\param dim 			: dimension of the matrix
	  	*/
	   void gauss_elimination( 	Epetra_SerialDenseMatrix& A, 
	                           	Epetra_SerialDenseVector& b,
	                           	Epetra_SerialDenseVector& x,
	                           	bool do_piv,
	                           	const int dim);
										
			
		/*!
	  	\brief updates the systemmatrix at the corresponding reference coordinates 
	
		\param A  				 (out)		: system matrix
		\param b  				 (out)		: right-hand-side
		\param xsi				 (in)			: vector of reference coordinates
	  	\param surfaceElement (in) 		: surface element
	  	\param lineElement    (in) 		: line element
	  	*/							
		void updateForCurveSurfaceIntersection(	Epetra_SerialDenseMatrix& A,
				 												Epetra_SerialDenseVector& b,
             												Epetra_SerialDenseVector& xsi,
             												DRT::Element* surfaceElement,			
	  	 														DRT::Element* lineElement);
	  	 														
	  	 														
      /*!
	  	\brief updates the system matrix at the corresponding reference coordinates for the 
	  			 computation if a certain point in current coordinates lies within an element 
	
		\param A  				 (out)		: system matrix
		\param b  				 (out)		: right-hand-side	  	
	  	\param xsi				 (in)			: vector of reference coordinates
	  	\param x  				 (in)			: point in current coordinates
	  	\param element 		 (in) 		: surface element
	  	\return systemmatrix
	  	*/			      																						
      void updateForPointWithinElement( Epetra_SerialDenseMatrix& A,
													 Epetra_SerialDenseVector& b,
													 Epetra_SerialDenseVector& xsi,
													 Epetra_SerialDenseVector& x,
													 DRT::Element* element) ;        
      
		
		/*!
	  	\brief computes an interseticon point between a curve and a surface
	  		
	  		The nonlinear system of equation is solved with help of the Newton-method.
	
	  	\param surfaceElement (in) 		: surface element
	  	\param lineElement    (in) 		: line element
	  	\param xsi				 (in/out)	: starting value/vector of reference coordinates
	  	\update 					 (in)			: if true system matrix is updated at each iteration step
	  	return intersection = true if an intersection point was found, otherwise false	
	  	*/
   	bool computeCurveSurfaceIntersection(  DRT::Element* surfaceElement,
				 											DRT::Element* lineElement,
             											Epetra_SerialDenseVector& xsi,
             											Epetra_SerialDenseVector& E1,
             											Epetra_SerialDenseVector& E2);
      
      /*!
	  	\brief checks if a point in current coordinates lies within a certain element
	  		
	  		The nonlinear system of equation is solved with help of the Newton-method.
	
	  	\param element 		 (in) 	   : element
	  	\param xsi				 (in/out)	: starting value/vector of reference coordinates
	  	\param x 				 (in)	      : point in current coordinates
	  	\update 					 (in)			: if true system matrix is updated at each iteration step
	  	\dim 						 (in)			: dimension of systemmatrix
	  	return pointWithinElement = true if the point lies within the element, false otherwise
	  	*/     									
		bool checkPointWithinElement(		DRT::Element* element,
				 									Epetra_SerialDenseVector& xsi,
             									Epetra_SerialDenseVector& x);
  
             									
          
      /*!
	  	\brief collects all intersection points between a line and a surface
	  		
	  	\param surfaceElement       (in)    : surface element
	  	\param lineElement          (in)    : line element
	  	\param intersectionPointSet (out)	: set of intersection points
	  	return number of intersection points 
	  	*/     				 									
		int collectIntersectionPoints(	DRT::Element* surfaceElement,
				 									DRT::Element* lineElement,
				 									std::set< Epetra_SerialDenseVector& >& intersectionPointSet);
				 									
       
      /*!
	  	\brief computes a new starting points for the Newton method recursively
	  			 in order to find all intersection points
	
	  	\param surfaceElement       (in) : surface element
	  	\param lineElement          (in) : line element
	  	\param E1                   (in)	: first endpoint
	  	\param E2                   (in)	: second endpoint
	  	\param intersectionPointSet (out): set of intersection points
	  	*/     				 	            												
		void computeNewStartingPoint(		DRT::Element* surfaceElement,
													DRT::Element* lineElement,
													Epetra_SerialDenseVector& E1,
     												Epetra_SerialDenseVector& E2,
     												std::set< Epetra_SerialDenseVector& >& intersectionPointSet);
   
};

#endif  // #ifndef INTERSECTION_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3_XFEM


/*!
\file drt_intersection.H

\brief collection of intersection tools

<pre>
Maintainer: Ursula Mayer
</pre>
*/

#ifdef XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef INTERSECTION_H
#define INTERSECTION_H


#include "../drt_xfem/integrationcell.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_element.H"
#include "../drt_f3/fluid3.H"
#include "/home/axelchen/include/tetgen.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RefCountPtr.hpp"

#if defined(__cplusplus)
extern "C"
{
#endif
#include <qhull/qhull.h>
#include <qhull/qset.h>
#if defined(__cplusplus)
}
#endif

using namespace std;
using namespace DRT;
using namespace DRT::Utils;
using namespace Teuchos;


typedef struct
{
  int           nsurf;
  int           surfaces[3];    /*!< surfaces*/
  double        coord[3];       /*!< coordinates */
} InterfacePoint;



class Intersection
{

       public:

                                                       
        /*!
        \brief computes the intersection between to different fields
    
        \param actdis               (in)    : discretization
        \param integrationcellList  (out)   : list of integrationcells for each intersected element
        */            
        void computeIntersection(   RefCountPtr<DRT::Discretization>        actdis, 
                                    map< int, vector <Integrationcell> >&	integrationcellList);
		
        
        
        /*!
        \brief Debugging integration cells (DEBUG ONLY)
    
        \param integrationcellList                (out)       : list of integration cells
        */                                  
       	void debugIntegrationcells(	map< int, vector <Integrationcell> >&	integrationcellList,
									int id);  
        
        private:
        
        
		/*!
	  	\brief Returns the sum of two Epetra_SerialDenseVectors
	
	  	\param v1 (in) : arbitrary Epetra_SerialDenseVector
	  	\param v2 (in) : arbitrary Epetra_SerialDenseVector
	  	\return sum of two Epetra_SerialDenseVector
	  	*/
        Epetra_SerialDenseVector addTwoVectors(     Epetra_SerialDenseVector&   v1,
                                                    Epetra_SerialDenseVector&   v2);
		
		
		/*!
	  	\brief Returns the difference of two double vectors v1 - v2
	  	\param v1 (in) : arbitrary double vector		
	  	\param v2 (in) : arbitrary double vector
	  	\return difference of two double vectors v1 - v2
	  	*/
        std::vector<double> subtractsTwoVectors(    std::vector <double>&   v1,
                                                    std::vector <double>&   v2);
			
            														
	  	/*!
	  	\brief Computes a rough overestimating axis-aligned bounding 
               box for an element (AABB)
	
	   	\param element 		 (in) 	   : element (DRT::Element*)
		\return axis-aligned bounding box for an element (AABB)
	  	*/
		Epetra_SerialDenseMatrix computeFastAABB( DRT::Element*   element);
		                      
	                           
		/*!
	  	\brief Checks if a node is within an axis-aligned bounding box (AABB)
	
	  	\param node (in) 		: node to checked
	  	\param AABB (in) 		: axis-aligned bouning box
	  	\return true if node is within the AABB or false otherwise
	  	*/
        bool isNodeWithinAABB(  std::vector<double>&         node,
                                Epetra_SerialDenseMatrix&    AABB);
										
													
		/*!
	  	\brief checks if two axis aligned bounding boxes intersect
	
	  	\param cutterAABB (in) 		: AABB of the cutting element
	  	\param xfemAABB   (in) 		: AABB of the xfem element
	  	\return true if the AABB's intersect or false otherwise
	  	*/
        bool intersectionOfAABB(    Epetra_SerialDenseMatrix&    cutterAABB,
                                    Epetra_SerialDenseMatrix&    xfemAABB);

																		
		/*!
	  	\brief computes a Gaussian elimination for a linear system of equations
	
	  	\param A        (in)    : system matrix
	  	\param b        (in)    : right-hand-side
	  	\param x        (out)   : solution vector
	  	\param do_piv   (in)	: do_piv = true does pivoting, do_piv = false does not do pivoting
	  	\param dim 	    (in)	: dimension of the matrix
        \return true if matrix is not singular , false if matrix is singular
	  	*/
        bool gaussElimination(  Epetra_SerialDenseMatrix&   A,
                                Epetra_SerialDenseVector&   b,
                                Epetra_SerialDenseVector&   x,
                                bool                        do_piv,
                                const int                   dim);
			
            							
        /*!
        \brief updates the system matrix at the corresponding reference coordinates for the 
               computation if a node in current coordinates lies within an element 
    
        \param A                 (out)      : system matrix
        \param xsi               (in)       : vector of reference coordinates
        \param element           (in)       : element (volume)
        */                                                                                                      
        void updateAForNWE( Epetra_SerialDenseMatrix&   A,
                            Epetra_SerialDenseVector&   xsi,
                            DRT::Element*               element);
          
          
        /*!
        \brief updates the rhs at the corresponding reference coordinates for the 
               computation whether a node in current coordinates lies within an element 
    
        \param b                 (out)      : right-hand-side       
        \param xsi               (in)       : vector of reference coordinates
        \param x                 (in)       : node in current coordinates
        \param element           (in)       : element (volume)
        */                               
        void updateRHSForNWE(   Epetra_SerialDenseVector&   b,
                                Epetra_SerialDenseVector&   xsi,
                                Epetra_SerialDenseVector&   x,
                                DRT::Element*               element);         
        
        
        /*!
        \brief collects points belonging to an interface lying inside xfem elements         
               The nonlinear system of equation is solved with help of the Newton-method.
    
        \param element              (in)        : element (surface)
        \param node                 (in)        : node in current coordinates
        \param interfacePoints      (in/out)    : vector of interface points
        \param nodeId               (in)        : node id
        \param numInterfacePoints   (in/out)    : number of interface points
        return true if the node lies within the element, false otherwise
        */                                      
		bool collectInternalPoints(	DRT::Element*                   element,
                                  	DRT::Node*                      node,
                                	std::vector<InterfacePoint>&    interfacePoints,
                                	int                             elemId,
                                 	int                             nodeId,
                                   	int&                            numInterfacePoints);
			
		/*!
        \brief checks if a node in current coordinates lies within a certain element           
               The nonlinear system of equation is solved with help of the Newton-method.
    
        \param element              (in)        : element (surface)
        \param node                 (in)        : node in current coordinates
        */	
		bool checkNodeWithinElement(	DRT::Element* element,
 										DRT::Node*  node);	
        
		/*!
	  	\brief updates the systemmatrix at the corresponding reference coordinates 
               for the computation of curve surface intersections
	
		\param A  				 (out)		: system matrix
		\param xsi				 (in)	    : vector of reference coordinates
	  	\param surfaceElement    (in) 		: surface element
	  	\param lineElement       (in) 		: line element
	  	*/							
        void updateAForCSI( Epetra_SerialDenseMatrix&   A,
                            Epetra_SerialDenseVector&   xsi,
                            DRT::Element*               surfaceElement,
                            DRT::Element*               lineElement); 
		
        	            
        /*!
        \brief updates the rhs at the corresponding reference coordinates 
               for the computation of curve surface intersections
   
        \param b                 (out)      : right-hand-side
        \param xsi               (in)       : vector of reference coordinates
        \param surfaceElement    (in)       : surface element
        \param lineElement       (in)       : line element
        */                          	 					
        void updateRHSForCSI(   Epetra_SerialDenseVector&   b,
                                Epetra_SerialDenseVector&   xsi,
                                DRT::Element*               surfaceElement,
                                DRT::Element*               lineElement);       		


		/*!
	  	\brief computes an interseticon point between a curve and a surface
	  		
	  		The nonlinear system of equation is solved with help of the Newton-method.
	
	  	\param surfaceElement   (in) 		: surface element
	  	\param lineElement      (in) 		: line element
	  	\param xsi              (in/out)	: starting value/vector of reference coordinates
        \param upLimit          (in)        : upper search interval boundary
        \param loLimit          (int)       : lower search interval boundary	  	
	  	return true if an intersection point was found, otherwise false	
	  	*/
        bool computeCurveSurfaceIntersection(   DRT::Element*               surfaceElement,
                                                DRT::Element*               lineElement,
                                                Epetra_SerialDenseVector&   xsi,
                                                Epetra_SerialDenseVector&   upLimit,
                                                Epetra_SerialDenseVector&   loLimit);
        
                
        /*!
	  	\brief collects all intersection points between a line and a surface
	  		
	  	\param surfaceElement           (in)        : surface element
	  	\param lineElement              (in)        : line element
	  	\param interfacePointList       (out)       : vector of interface points
        \param numInterfacePoints       (in/out)    : number of interface points
        \param surfaceId                (in)        : surface element id
        \param lineId                   (in)        : line element id
        \param lines                    (in)        : if lines = true
	  	*/     				 									
        void collectIntersectionPoints(	DRT::Element*                   surfaceElement,
        								DRT::Element*                   lineElement,
        								std::vector<InterfacePoint>&    interfacePoints,
        								int&                            numInterfacePoints,
        								int                             surfaceId,
        								int                             lineId,
        								bool                            lines);
		
        	      
        /*!
	  	\brief computes a new starting points for the Newton method recursively
	  		   in order to find all intersection points
	
	  	\param surfaceElement           (in)    : surface element
	  	\param lineElement              (in)    : line element
	  	\param surfaceId                (in)   	: surface element id
        \param lineId                   (in)   	: line element id
	  	\param upLimit                  (in)	: upper limit of the search interval 
	  	\param loLimit                  (in)	: lower limit of the search interval
	  	\param interfacePointList       (out)   : vector of interface points
	  	\param lines                    (in)	: if lines = true
	  	\return number of interface points
	  	*/     				 	            												
        int computeNewStartingPoint(	DRT::Element*                surfaceElement,
										DRT::Element*                lineElement,
                                      	int                          surfaceId,
    									int                          lineId,
										Epetra_SerialDenseVector&    upLimit,
 										Epetra_SerialDenseVector&    loLimit,
 										std::vector<InterfacePoint>& interfacePoints,
                                        bool                         lines);
  
  
  		/*!
	  	\brief adds an intersection point to the list of interface points
	
		\param surfaceElement           (in)    : surface element
	  	\param lineElement              (in)    : line element
	  	\param xsi           			(in)    : surface element
	  	\param upLimit                  (in)	: upper limit of the search interval 
	  	\param loLimit                  (in)	: lower limit of the search interval 
	  	\param interfacePointList       (out)   : vector of interface points
	  	\param surfaceId                (in)   	: surface element id
        \param lineId                   (in)   	: line element id
		\param lines                    (in)  	: if lines = true        
	  	\return number of interface points
	  	*/     	
  		int addIntersectionPoint(	DRT::Element*                	surfaceElement,
									DRT::Element*                	lineElement,
  									Epetra_SerialDenseVector&		xsi,
									Epetra_SerialDenseVector&     	upLimit,
									Epetra_SerialDenseVector&     	loLimit,
									std::vector<InterfacePoint>& 	interfacePoints,
									int                      		surfaceId,
		        					int                       		lineId,							
									bool 							lines);
									
									
        /*!
        \brief computes a new starting points for the Newton method recursively
               in order to find all intersection points
    
        \param element          (in)    : surface element
        \param xsi              (in)    : vector of reference coordinates
        */
        void referenceToCurrentCoordinates( DRT::Element*               element, 
                                            Epetra_SerialDenseVector&   xsi);
  
  
        /*!
        \brief compares two points (overloaded method)
    
        \param point1       (in)    : first point   (vector<double>)
        \param point2       (in)    : second point  (double*)
        \return true if both points equal each other, false otherwise
        */
        bool comparePoints( 	vector<double>&     point1,
        						double*             point2);
          
                            
        /*!
        \brief compares two points (overloaded method)
    
        \param point1       (in)    : first point   (vector<double>)
        \param point2       (in)    : second point  (vector<double>)
        \return true if both points equal each other, false otherwise
        */                    
        bool comparePoints(	vector<double>&     point1,
                            vector<double>&     point2);
                            
           
      	/*!
        \brief compares two points (overloaded method)
    
        \param point1       (in)    : first point   (Epetra_SerialDenseVector)
        \param point2       (in)    : second point  (Epetra_SerialDenseVector)
        \return true if both points equal each other, false otherwise
        */                      
        bool comparePoints(	Epetra_SerialDenseVector&     point1,
        					Epetra_SerialDenseVector&     point2);
        
        
  		/*!
        \brief finds the next facet of a convex hull in clockwise order
    
        \param vertices         (in)        : vector of facet vertices
        \param searchPoint      (in/out)    : common point of the previous and next facet/the new point of the next facet
        */          
        void findNextFacet( vector< vector<double> >&   vertices,
                            vector<double>&             searchPoint);
  
          
        /*!
        \brief  stores an interface point in a point list for the 
                Constrained Delaunay Tetrahedralization (CDT) with Tetgen
    
        \param point            (in)        :   coordinates of the point to be stored 
        \param interfacePoints  (in/out)    :   vector of interface points
        \param positions        (in/out)    :   positions
        \param pointList        (in/out)    :   point list, data structure for CDT with Tetgen
        */                          
        void storePoint(    vector<double>&             point, 
                            vector<InterfacePoint>&     interfacePoints, 
                            vector<int>&                positions, 
                            vector<InterfacePoint>&     pointList );
  
   
        /*!
        \brief  stores a segment in a segment list for the 
                Constrained Delaunay Tetrahedralization (CDT) with Tetgen
    
        \param pointList            (in)        :   coordinates of the point to be stored 
        \param positions            (in)        :   positions
        \param segmentList          (in/out)    :   segment list, data structure for CDT with Tetgen
        */            
        void storeSegments( vector<InterfacePoint>&     pointList, 
                            vector<int>                 positions, 
                            vector< vector<int> >&      segmentList);
          
                            
        /*!
        \brief  stores a triangle facet in a triangle list for the 
                Constrained Delaunay Tetrahedralization (CDT) with Tetgen
    
        \param pointList            (in)        :   coordinates of the point to be stored 
        \param positions            (in)        :   positions
        \param triangleList         (in/out)    :   triangle list, data structure for CDT with Tetgen
        */                         
        void storeTriangles(    vector<InterfacePoint>&     pointList, 
                                vector<int>                 positions, 
                                vector< vector<int> >&      triangleList);
         
                   
        /*!
        \brief  computes the midpoint of a collection of points
    
        \param interfacePoints            (in)        : vector of interface points
        \return returns the midpoint (type InterfacePoint)
        */                                 
        InterfacePoint computeMidpoint( vector<InterfacePoint>&     interfacePoints);          					      
        
       
       	/*!
        \brief  computes unique coordinates for a certain region in the 
        		xfem element. It is copied into the data structure of Tetgen
        		   
        \param xfemElement           (in)        : xfem element
        \param cutterElement         (in)        : cutter element
        \param regioncoordinates	 (out)		 : coordinates of the region points
        */    
		void computeRegionCoordinates(	DRT::Element*  xfemElement,
 										DRT::Element*  cutterElement,
 										double* regionCoordinates);
 										
 										       
        /*!
        \brief  computes the convex hull of a set of points
    
        \param surfaceElement           (in)        : surface element
        \param interfacePoints          (in)        : list of interface points
        \param pointList                (out)       : list of points
        \param segmentList              (out)       : list of segments
        \param triangleList             (out)       : list of triangle facets
        */           
        void computeConvexHull( DRT::Element*               surfaceElement,
                                vector<InterfacePoint>&     interfacePoints, 
                                vector< InterfacePoint >&   pointList,
                                vector< vector<int> >&      segmentList,
                                vector< vector<int> >&      triangleList );
    
    
        /*!
        \brief  computes the convex hull of a set of points
    
        \param element                  (in)        : element
        \param cutterElement            (in)        : cutter element
        \param pointList                (out)       : list of points
        \param segmentList              (out)       : list of segments
        \param triangleList             (out)       : list of triangle facets
        */       
        void computeCDT(    DRT::Element*               element,
        					DRT::Element*  				cutterElement,
                            vector< InterfacePoint >&   pointList,
                            vector< vector<int> >&      segmentList,
                            vector< vector<int> >&      triangleList,
                            map< int, vector <Integrationcell> >&	integrationcellList);
    
          
        /*!
        \brief Debugging the intersection of AABB's (DEBUG ONLY)
    
        \param cutterAABB       (in)    : AABB of the cutting element
        \param xfemAABB         (in)    : AABB of the xfem element
        \param cutterElement    (in)    : cutting element
        \param xfemElement      (in)    : xfem element
        \param noC              (in)    : id of the cutting element
        \param noX              (in)    : id of the xfem element
        */                        
        void debugAABBIntersection( Epetra_SerialDenseMatrix    cutterAABB,
                                    Epetra_SerialDenseMatrix    xfemAABB,
                                    DRT::Element*               cutterElement,
                                    DRT::Element*               xfemElement,
                                    int                         noC,
                                    int                         noX);
                                    
          
        /*!
        \brief Debugging node within element (DEBUG ONLY)
    
        \param element          (in)    : element
        \param node             (in)    : node
        \param xsi              (in)    : reference coordinates
        \param noE              (in)    : id of the element
        \param noN              (in)    : id of the node
        \param within           (in)    : true if within, false otherwise
        */                                     
        void debugNodeWithinElement(    DRT::Element*               element,
                                        DRT::Node*                  node,
                                        Epetra_SerialDenseVector&   xsi,
                                        int                         noE,
                                        int                         noN,
                                        bool                        within);
                                        
                                        
        /*!
        \brief Debugging curve surface intersection (DEBUG ONLY)
    
        \param surfaceElement   (in)    : surfaceElement
        \param lineElement      (in)    : lineElement 
        \param xsi              (in)    : reference coordinates
        \param noSE             (in)    : id of the surface element
        \param noLE             (in)    : id of the line element
        \param within           (in)    : true if within, false otherwise
        */                                     
        void debugCurveSurfaceIntersection(     DRT::Element*               surfaceElement,
                                                DRT::Element*               lineElement,
                                                Epetra_SerialDenseVector&   xsi,
                                                int                         noSE,
                                                int                         noLE,
                                                bool                        within);
    
      
        /*!
        \brief Debugging tetgen data structure (DEBUG ONLY)
    
        \param pointList                (out)       : list of points
        \param segmentList              (out)       : list of segments
        \param triangleList             (out)       : list of triangle facets
        */                             
        void debugTetgenDataStructure(  vector< InterfacePoint >&   pointList,
                                        vector< vector<int> >&      segmentList,
                                        vector< vector<int> >&      triangleList);  
                                                                                                                            
};




#endif  // #ifndef INTERSECTION_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef XFEM


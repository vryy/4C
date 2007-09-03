/*!
\file drt_intersection.H

\brief collection of intersection tools

<pre>
Maintainer: Ursula Mayer
</pre>
*/

#ifdef D_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef INTERSECTION_H
#define INTERSECTION_H


#include "../drt_xfem/integrationcell.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_element.H"
#include "../drt_f3/fluid3.H"
#include "tetgen.h"
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

namespace XFEM
{

enum pointType {internalP, surfaceP, intersectionP};

enum lineType {inPlane , inElement, inPlaneCurved, inElementCurved};

typedef struct
{
  int           nsurf;
  int           surfaces[3];    /*!< surfaces*/
  double        coord[3];       /*!< coordinates */
  pointType     pType;
} InterfacePoint;



class Intersection
{

       public:

                                                       
        /*!
        \brief computes the intersection between to different fields
    
        \param xfemdis              (in)    : xfem discretization
        \param cutterdis            (in)    : cutter discretization
        \param integrationcellList  (out)   : list of integrationcells for each intersected element
        */            
        void computeIntersection(   
            const RefCountPtr<DRT::Discretization>        xfemdis, 
            const RefCountPtr<DRT::Discretization>        cutterdis, 
            map< int, vector <Integrationcell> >&	integrationcellList
            );
		
        
        
        private:
        
        
        bool higherorder_;
        
        
        //! number of surfaces of the xfem element
        int numXFEMSurfaces_;
        
        
        //! number of lines of the xfem element
        int numXFEMLines_;
        
        
        //! number of corner nodes of the xfem element
        int numXFEMCornerNodes_;
        
        
        //! tolerance 1e-14 used for residuals in all Newton-methods
        static const double TOL14_ = 1e-14;
        
        
        //! tolerance 1e-7 
        static const double TOL7_ = 1e-7;
        
        
        //! lists the node indices for each surface of the xfem element
        vector< vector<int> > eleNodeNumberingSurfaces_;
        
        
        //! list of all surfaces adjacent to an xfem element line
        vector< vector<int> > eleNodeNumberingLinesSurfaces_;
       
       
        //! //! list of all surfaces adjacent to an xfem element node
        vector< vector<int> > eleNodeNumberingNodesSurfaces_; 
        
        
        //! list of all nodes of the xfem element in reference coordinates
        vector< vector<double> > eleNodeNumberingNodesReference_; 
        
        
        //! list of intersected surface elements
        vector< DRT::Element* >          intersectedSurfaces_; 
        
        
        /*!
        \brief initializes the private members
    
        \param xfemdis   (in) : discretization of the xfem mesh
        \param cutterdis (in) : discretization of the cutter mesh
        */
        void initialize(  
            const RefCountPtr<DRT::Discretization> xfemdis,
            const RefCountPtr<DRT::Discretization> cutterdis
            );
        
                
		/*!
	  	\brief Returns the sum of two Epetra_SerialDenseVectors
	
	  	\param v1 (in) : arbitrary Epetra_SerialDenseVector
	  	\param v2 (in) : arbitrary Epetra_SerialDenseVector
	  	\return sum of two Epetra_SerialDenseVectors
	  	*/
        Epetra_SerialDenseVector addTwoVectors(     
            Epetra_SerialDenseVector&   v1,
            Epetra_SerialDenseVector&   v2
            );		
		
        
        /*!
        \brief Returns the sum of two vector<double>
    
        \param v1 (in) :    arbitrary vector<double>
        \param v2 (in) :    arbitrary vector<double>
        \return sum of two vector<double>
        */
        vector<double> addTwoVectors(   
            vector<double>&   v1,
            vector<double>&   v2
            );
        
                                              
        /*!
        \brief Returns the difference of Epetra_SerialDenseVectors v1 - v2
        \param v1 (in) : Epetra_SerialDenseVector      
        \param v2 (in) : Epetra_SerialDenseVector
        \return difference of Epetra_SerialDenseVectors v1 - v2
        */
        Epetra_SerialDenseVector subtractsTwoVectors(   
            Epetra_SerialDenseVector&   v1,
            Epetra_SerialDenseVector&   v2
            );       
        
                                            
		/*!
	  	\brief Returns the difference of two vector<double> v1 - v2
	  	\param v1 (in) : arbitrary vector<double>		
	  	\param v2 (in) : arbitrary vector<double>
	  	\return difference of two vector<double> v1 - v2
	  	*/
        vector<double> subtractsTwoVectors( 
            vector <double>&   v1,
            vector <double>&   v2
            );
		
        	     
        /*!
        \brief computes the cross product of two Epetra_SerialDenseVectors a x b
        \param a (in) : arbitrary Epetra_SerialDenseVector     
        \param b (in) : arbitrary Epetra_SerialDenseVector
        \return cross product of two Epetra_SerialDenseVectors a x b
        */
        Epetra_SerialDenseVector computeCrossProduct(   
            const Epetra_SerialDenseVector& a,
            const Epetra_SerialDenseVector& b
            );   
        
                                                                                               
        /*!
        \brief compute absolute value (L2-norm) of a Epetra_SerialDenseVector
        \param a    (in) : arbitrary Epetra_SerialDenseVector  
        \return absolute value
        */                                            
        double computeAbsoluteValue( 
            Epetra_SerialDenseVector&   a
            );
        
                       
        /*!
        \brief compute a scalar vector multiplication
        \param s (in) : scalar  
        \param v (in) : Epetra_SerialDenseVector
        \return Epetra_SerialDenseVector resulting from scalar vector multiplication
        */       
        Epetra_SerialDenseVector computeScalarVectorMultiplication( 
            const double                s, 
            Epetra_SerialDenseVector&   v
            ); 
        
                																
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
                                const bool                  do_piv,
                                const int                   dim,
                                const int 					order);
		
        
        /*!
        \brief transforms a point in reference coordinates to a point
               in current coordinates with respect to a given element
    
        \param element          (in)        : element
        \param xsi              (in/out)    : vector of reference coordinates/
                                              vector of current coordinates
        */
        void referenceToCurrentCoordinates( 
            DRT::Element*               element, 
            Epetra_SerialDenseVector&   xsi
            );
 
 
        /*!
        \brief transforms a point in current coordinates to a point
               in reference coordinates with respect to a given element
    
        \param element          (in)        : element
        \param xsi              (in/out)    : vector of current coordinates/
                                              vector of reference coordinates
        */
        void currentToReferenceCoordinates( 
            DRT::Element*               element,                     
            Epetra_SerialDenseVector&   xsi
            );
 
            
        /*!
        \brief compares two points (overloaded method)
    
        \param point1       (in)    : first point   (vector<double>)
        \param point2       (in)    : second point  (double*)
        \return true if both points equal each other, false otherwise
        */
        bool comparePoints( 
            const vector<double>&     point1,
            const double*             point2
            );
 
             
        /*!
        \brief compares two points (overloaded method)
    
        \param point1       (in)    : first point   (double*)
        \param point2       (in)    : second point  (double*)
        \param lenght       (in)    : length of arrays
        \return true if both points equal each other, false otherwise
        */
        bool comparePoints( 
            const double*     point1,
            const double*     point2,
            const int         length
            ); 
 
                                    
        /*!
        \brief compares two points (overloaded method)
    
        \param point1       (in)    : first point   (vector<double>)
        \param point2       (in)    : second point  (vector<double>)
        \return true if both points equal each other, false otherwise
        */                    
        bool comparePoints( 
            const vector<double>&     point1,
            const vector<double>&     point2
            );
 
                            
        /*!
        \brief compares two points (overloaded method)
    
        \param point1       (in)    : first point   (Epetra_SerialDenseVector)
        \param point2       (in)    : second point  (Epetra_SerialDenseVector)
        \return true if both points equal each other, false otherwise
        */                      
        bool comparePoints( 
            const Epetra_SerialDenseVector&     point1,
            const Epetra_SerialDenseVector&     point2
            );
 
        
        /*!
        \brief checks if element is a volume element
    
        \param element       (in)    : element
        \return true if it is a volume element, false otherwise
        */      
        bool checkIfVolumeElement(
            DRT::Element* element
            );
        
        
        /*!
        \brief checks if element is a surface element
    
        \param element       (in)    : element
        \return true if it is a surface element, false otherwise
        */     
        bool checkIfSurfaceElement(
            DRT::Element* element
            );
        
        
        /*!
        \brief checks if element is a line element
    
        \param element       (in)    : element
        \return true if it is a line element, false otherwise
        */     
        bool checkIfLineElement(
            DRT::Element* element
            );
           
           
        /*!
        \brief returns the order of the element
    
        \param element       (in)    : element
        \return the order of the element
        */    
        int getDimension(
            DRT::Element* element
            );     
        
        
        /*!
        \brief Computes a rough overestimating extended  
               axis-aligned bounding box for an element (XAABB)
    
        \param  element          (in)      : element (DRT::Element*)
        \return extended axis-aligned bounding box  (XAABB) for an element
        */
        Epetra_SerialDenseMatrix computeFastXAABB( 
            DRT::Element*   element
            );
        
                                                
        /*!
        \brief Checks if a node is within an axis-aligned bounding box (XAABB)
    
        \param node (in)        : node to be checked
        \param XAABB (in)       : extended axis-aligned bouning box
        \return true if node is within the XAABB or false otherwise
        */
        bool isNodeWithinXAABB(  
            std::vector<double>&         node,
            Epetra_SerialDenseMatrix&    XAABB);
        
                                                                                        
        /*!
        \brief checks if two extended axis aligned bounding boxes intersect
    
        \param cutterXAABB (in)         : XAABB of the cutting element
        \param xfemXAABB   (in)         : XAABB of the xfem element
        \return true if the XAABB's intersect or false otherwise
        */
        bool intersectionOfXAABB(    
            Epetra_SerialDenseMatrix&    cutterXAABB,
            Epetra_SerialDenseMatrix&    xfemXAABB
            );
        
        	       
        /*!
        \brief collects points belonging to an interface lying inside xfem elements         
               The nonlinear system of equation is solved with help of the Newton-method.
    
        \param element              (in)        : element (surface)
        \param node                 (in)        : node in current coordinates
        \param interfacePoints      (in/out)    : vector of interface points
        \param nodeId               (in)        : node id
        \param numInternalPoints    (in/out)    : number of interface points
        \param numSurfacePoints     (in/out)    : number points lying on a surface or node of an xfem-element
        \return true if the node lies within the element, false otherwise
        */                                      
        bool collectInternalPoints( 
            DRT::Element*                   element,
            DRT::Element*                   surfaceElement,
            DRT::Node*                      node,
            std::vector<InterfacePoint>&    interfacePoints,
            const int                       elemId,
            const int                       nodeId,
            int&                            numInternalPoints,
            int&                            numSurfacePoints);
        
                       							
        /*!
        \brief updates the system matrix at the corresponding reference coordinates for the 
               computation if a node in current coordinates lies within an element 
    
        \param A                 (out)      : system matrix
        \param xsi               (in)       : vector of reference coordinates
        \param element           (in)       : element (volume)
        */                                                                                                      
        void updateAForNWE( 
            Epetra_SerialDenseMatrix&   A,
            Epetra_SerialDenseVector&   xsi,
            DRT::Element*               element
            );
        
          
        /*!
        \brief updates the rhs at the corresponding reference coordinates for the 
               computation whether a node in current coordinates lies within an element 
    
        \param b                 (out)      : right-hand-side       
        \param xsi               (in)       : vector of reference coordinates
        \param x                 (in)       : node in current coordinates
        \param element           (in)       : element (volume)
        */                               
        void updateRHSForNWE(   
            Epetra_SerialDenseVector&   b,
            Epetra_SerialDenseVector&   xsi,
            Epetra_SerialDenseVector&   x,
            DRT::Element*               element
            );         
        
            
		/*!
        \brief checks if a node in current coordinates lies within a certain element           
               The nonlinear system of equation is solved with help of the Newton-method.
    
        \param element              (in)        : element (surface)
        \param node                 (in)        : node in current coordinates
        \param xsi                  (in)        : node in reference coordinates
        */	
		bool checkNodeWithinElement(	
            DRT::Element* element,
            Epetra_SerialDenseVector& x,
            Epetra_SerialDenseVector& xsi
            );	       
        
        
        /*!
        \brief checks if a node that is within a node lies on one of its surfaces or node          .
    
        \param element              (in)        : element (surface)
        \param xsi                  (in)        : node in reference coordinates
        \param ip                   (in)        : interface point
        \param numSurfacePoints     (in)        : number of points lying on a surface and a node
        \return true if the node is lying on a surface or a node, false otherwise
        */  
        bool checkIfOnSurfaceAndNode(   
            DRT::Element*                   element,
            Epetra_SerialDenseVector&       xsi, 
            InterfacePoint&                 ip,
            int&                            numSurfacePoints);
        
              
        /*!
        \brief collects all intersection points between a line and a surface
            
        \param surfaceElement           (in)        : surface element
        \param lineElement              (in)        : line element
        \param interfacePointList       (out)       : vector of interface points
        \param numInternalPoints        (in)        : number of internal points
        \param numSurfacePoints         (in)        : number of surface points
        \param surfaceId                (in)        : surface element id
        \param lineId                   (in)        : line element id
        \param lines                    (in)        : if lines = true
        */                                                      
        void collectIntersectionPoints( 
            DRT::Element*                   surfaceElement,
            DRT::Element*                   lineElement,
            std::vector<InterfacePoint>&    interfacePoints,
            const int                       numInternalPoints,
            const int                       numSurfacePoints,
            const int                       surfaceId,
            const int                       lineId,
            const bool                      lines,
            bool&                           xfemIntersection
            );
        
         
        /*!
        \brief computes a new starting points for the Newton method recursively
               in order to find all intersection points
    
        \param surfaceElement           (in)    : surface element
        \param lineElement              (in)    : line element
        \param surfaceId                (in)    : surface element id
        \param lineId                   (in)    : line element id
        \param upLimit                  (in)    : upper limit of the search interval 
        \param loLimit                  (in)    : lower limit of the search interval
        \param interfacePointList       (out)   : vector of interface points
        \param lines                    (in)    : if lines = true
        \return number of interface points
        */                                                                                  
        int computeNewStartingPoint(    
            DRT::Element*                   surfaceElement,
            DRT::Element*                   lineElement,
            const int                       surfaceId,
            const int                       lineId,
            Epetra_SerialDenseVector&       upLimit,
            Epetra_SerialDenseVector&       loLimit,
            std::vector<InterfacePoint>&    interfacePoints,
            const bool                      lines
            );
        
              
        /*!
        \brief computes an interseticon point between a curve and a surface
            
            The nonlinear system of equation is solved with help of the Newton-method.
    
        \param surfaceElement   (in)        : surface element
        \param lineElement      (in)        : line element
        \param xsi              (in/out)    : starting value/vector of reference coordinates
        \param upLimit          (in)        : upper search interval boundary
        \param loLimit          (int)       : lower search interval boundary        
        return true if an intersection point was found, otherwise false 
        */
        bool computeCurveSurfaceIntersection(   
            DRT::Element*               surfaceElement,
            DRT::Element*               lineElement,
            Epetra_SerialDenseVector&   xsi,
            Epetra_SerialDenseVector&   upLimit,
            Epetra_SerialDenseVector&   loLimit
            );
        
                                                           
		/*!
	  	\brief updates the systemmatrix at the corresponding reference coordinates 
               for the computation of curve surface intersections
	
		\param A  				 (out)		: system matrix
		\param xsi				 (in)	    : vector of reference coordinates
	  	\param surfaceElement    (in) 		: surface element
	  	\param lineElement       (in) 		: line element
	  	*/							
        void updateAForCSI( 
            Epetra_SerialDenseMatrix&   A,
            Epetra_SerialDenseVector&   xsi,
            DRT::Element*               surfaceElement,
            DRT::Element*               lineElement
            ); 
		
          	            
        /*!
        \brief updates the rhs at the corresponding reference coordinates 
               for the computation of curve surface intersections
   
        \param b                 (out)      : right-hand-side
        \param xsi               (in)       : vector of reference coordinates
        \param surfaceElement    (in)       : surface element
        \param lineElement       (in)       : line element
        */                          	 					
        void updateRHSForCSI(   
            Epetra_SerialDenseVector&   b,
            Epetra_SerialDenseVector&   xsi,
            DRT::Element*               surfaceElement,
            DRT::Element*               lineElement);       		
        
            
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
  		int addIntersectionPoint(	
            DRT::Element*                	surfaceElement,
            DRT::Element*                	lineElement,
  		    Epetra_SerialDenseVector&		xsi,
            Epetra_SerialDenseVector&     	upLimit,
            Epetra_SerialDenseVector&     	loLimit,
            std::vector<InterfacePoint>& 	interfacePoints,
            const int                      	surfaceId,
            const int                       lineId,
            const bool 					    lines
            );
		
        				   							
        
        /*!
        \brief  computes the convex hull of a set of points
    
        \param surfaceElement           (in)        : surface element
        \param interfacePoints          (in)        : list of interface points
        \param pointList                (out)       : list of points
        \param segmentList              (out)       : list of segments
        \param triangleList             (out)       : list of triangle facets
        */           
        void computeConvexHull( 
            DRT::Element*               xfemElement,
            DRT::Element*               surfaceElement,
            vector<InterfacePoint>&     interfacePoints, 
            vector< InterfacePoint >&   pointList,
            vector< vector<int> >&      surfacePointList,
            vector< vector<int> >&      segmentList,
            vector< vector<int> >&      triangleList,
            const int                   numInternalPoints,
            const int                   numSurfacePoints 
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
            );
  
        
        /*!
        \brief  computes the Constrained Delaunay Tetrahedralization CDT
                of an intersected xfem element
    
        \param element                  (in)        : element
        \param cutterElement            (in)        : cutter element
        \param pointList                (out)       : list of points
        \param segmentList              (out)       : list of segments
        \param triangleList             (out)       : list of triangle facets
        */       
        void computeCDT(    
            DRT::Element*               element,
            DRT::Element*               cutterElement,
            vector< InterfacePoint >&   pointList,
            vector< vector<int> >&      surfacePointList,
            vector< vector<int> >&      segmentList,
            vector< vector<int> >&      triangleList,
            map< int, vector <Integrationcell> >&   integrationcellList
            );
    
    
        /*!
        \brief  adds the corner nodes of an xfem element to the pointlist
    
        \param pointList            (in/out)        : list of interface points
        */          
        void startPointList(
            vector< InterfacePoint >&  pointList
            );
            
            
        /*!
        \brief  stores an interface point in a point list for the 
                Constrained Delaunay Tetrahedralization (CDT) with Tetgen
    
        \param point            (in)        :   coordinates of the point to be stored 
        \param interfacePoints  (in/out)    :   vector of interface points
        \param positions        (in/out)    :   positions
        \param pointList        (in/out)    :   point list, data structure for CDT with Tetgen
        */                          
        void storePoint(    
            vector<double>&             point, 
            vector<InterfacePoint>&     interfacePoints, 
            vector<int>&                positions, 
            vector<InterfacePoint>&     pointList 
            );
  
         
        /*!
        \brief  stores a single point in a list of surface points if no 
                segment lies in that surface
    
        \param pointList            (in)        :   coordinates of the point to be stored 
        \param interfacePoints      (in)        :   list of interface points
        \param surfacePointList     (in/out)    :   list of single points lying of an surface
                                                    of an xfem element
        */      
        void storeSurfacePoints(    
            vector<InterfacePoint>&     pointList, 
            vector<InterfacePoint>&     interfacePoints,
            vector< vector<int> >&      surfacePointList
            );
            
                           
        /*!
        \brief  stores a segment in a segment list for the 
                Constrained Delaunay Tetrahedralization (CDT) with Tetgen
    
        \param pointList            (in)        :   coordinates of the point to be stored 
        \param positions            (in)        :   positions
        \param segmentList          (in/out)    :   segment list, data structure for CDT with Tetgen
        */            
        void storeSegments( 
            vector<InterfacePoint>&     pointList, 
            vector<int>&                positions, 
            vector< vector<int> >&      segmentList
            );
          
          
        /*!
        \brief  computes the midpoint of a collection of points
    
        \param interfacePoints            (in)        : vector of interface points
        \return returns the midpoint (type InterfacePoint)
        */                                 
        InterfacePoint computeMidpoint( 
            vector<InterfacePoint>&     interfacePoints
            );       
            
                               
        /*!
        \brief  stores a triangle facet in a triangle list for the 
                Constrained Delaunay Tetrahedralization (CDT) with Tetgen
    
        \param pointList            (in)        :   coordinates of the point to be stored 
        \param positions            (in)        :   positions
        \param triangleList         (in/out)    :   triangle list, data structure for CDT with Tetgen
        */                         
        void storeTriangles(    
            vector<InterfacePoint>&     pointList, 
            vector<int>                 positions, 
            vector< vector<int> >&      triangleList
            );
         
  
        /*!
        \brief  computes the recovery of the curved interface
                for one xfem element after a Contrained Delaunay 
                Tetrahedralization
    
        \param element              (in)        :   xfem element
        \param out                  (in/out)    :   output data structure of TetGen
        */      
        void recoverCurvedInterface(
            DRT::Element* element, 
            tetgenio& out
            );
 
 
        /*!
        \brief  computes the recovery of the curved interface
                for one xfem element after a Contrained Delaunay 
                Tetrahedralization. The surface element is
                transformed into the reference system of the
                xfem element, since the intersection is carried
                out in the reference configuration of the xfem element
    
        \param surfaceNodes             (in)    :   nodes of the surface element transformed 
                                                    into the reference system of the xfem element
        \param surfaceElement           (in)    :   surface element
        \param node1                    (in)    :   1.node of normal
        \param node2                    (in)    :   2.node of normal
        \param xsi                      (out)   :   intersection point
        \return true if normal intersects with curved surface element or false otherwise
        */      
        bool computeRecovery(  
            vector <Epetra_SerialDenseVector>&  surfaceNodes,
            DRT::Element*                       surfaceElement,
            Epetra_SerialDenseVector&           node1,
            Epetra_SerialDenseVector&           node2,
            Epetra_SerialDenseVector&           xsi
            );        
        
        
        /*!
        \brief  computes the recovery of the curved interface
                for one xfem element after a Contrained Delaunay 
                Tetrahedralization
    
        \param A                (in)    :   xfem element
        \param xsi              (in)    :   intersection point
        \param surfaceNodes     (in)    :   surface nodes transformed into the
                                            reference system of the xfem element
        \param surfaceElement   (in)    :   surface element
        */      
        void updateAForRQI( 
            Epetra_SerialDenseMatrix&           A,
            Epetra_SerialDenseVector&           xsi,
            vector <Epetra_SerialDenseVector>&  surfaceNodes,
            DRT::Element*                       surfaceElement
            );
         
        
        /*!
        \brief  computes the recovery of the curved interface
                for one xfem element after a Contrained Delaunay 
                Tetrahedralization
    
        \param b                (in)    :   right-hand-side
        \param xsi              (out)   :   intersection point
        \param surfaceNodes     (in)    :   surface nodes transformed into the
                                            reference system of the xfem element
        \param surfaceElement   (in)    :   surface element
        \param node1            (in)    :   1.node of a normal
        \param node2            (in)    :   2.node of a normal
        */                         
        void updateRHSForRQI(   
            Epetra_SerialDenseVector&           b,
            Epetra_SerialDenseVector&           xsi,
            vector <Epetra_SerialDenseVector>&  surfaceNodes,
            DRT::Element*                       surfaceElement,
            Epetra_SerialDenseVector&           node1,
            Epetra_SerialDenseVector&           node2
            );
       
    
        /*!
        \brief  stores the surface element actually intersecting 
                with an xfem element in a list
    
        \param surfaceElement   (in)    :   surface element
        */                  
        void storeIntersectedCutterElement(
            DRT::Element* surfaceElement
            );  
         
         
        /*!
        \brief  transforms the intersecting surface elements into the reference
                system of the xfem elements
    
        \param surfaceElement   (in)    :   element
        \param surfaceNodes     (out)   :   surface nodes transformed into the
                                            reference system of the xfem element
        */        
        void transformSurfaceNodes( 
            DRT::Element*                                   element,
            vector < vector <Epetra_SerialDenseVector> >&   surfaceNodes
            );
    
       
        /*!
        \brief  returns the index of the higher order node of an
                tetrahedra edge spezified by their end point indices
                
        \param index1            (in)        : index of an end point node
        \param index2            (in)        : index of an end point node
        \return index of the higher order node
        */
        int getHigherOrderTetIndex(
            const int index1, 
            const int index2
            ); 
            
            
        /*!
        \brief  gathers information of a tetrahedron adjacent to an 
                interface 
                
        \param out                  (in)    : output data structure of TetGen
        \param tetraCornerIndices   (out)   : index list of tetrahedron corners
        \param order                (out)   : order the tetrahedron corner points appear
        \param index                (in)    : tetrahedron index of the tetrahedron list
        */
        void getTetrahedraInformation(   
            tetgenio& out, 
            vector<int>& tetraCornerIndices,
            vector<int>& order,
            const int index
            );          
        
        
        /*!
        \brief  compute the midpoint of a line described by two points
                
        \param p1       (in)    : point 1
        \param p2       (in)    : point 2
        \return mid point of a line
        */    
        Epetra_SerialDenseVector computeLineMidpoint(  
            Epetra_SerialDenseVector& p1,
            Epetra_SerialDenseVector& p2
            );      


        /*!
        \brief  computes the normal to the interface edge of the
                adjacent tetrahedon. In addition the normal lies within the plane
                of the tetrahedron. 
                
        \param index1               (in)    : index 1 of an interface edge
        \param index2               (in)    : index 2 of an interface edge
        \param out                  (in)    : out
        \param tetraCornerIndices   (in)    : indices of tetrahedron corners
        \param node1                (in)    : node 1 of the normal
        \param node2                (in)    : node 2 of the normal
        */    
        void computeIntersectionNormal(    
            const int                   index1,
            const int                   index2,
            tetgenio&                   out,
            vector<int>&                tetraCornerIndices,
            Epetra_SerialDenseVector&   node1,
            Epetra_SerialDenseVector&   node2
            );
               
        
        /*!
        \brief  computes the normal to the interface edge of the
                adjacent tetrahedon. In addition the normal lies within the plane
                of the tetrahedron. 
                
        \param index                (in)    : index 1 of an interface edge
        \param xsi                  (in)    : coordinates of the higher order node 
        \param surface element      (in)    : surface element
        \param out                  (in)    : output data structure of tetGen
        */                              
        void storeHigherOrderNode(   
            const int                   index,
            Epetra_SerialDenseVector&   xsi, 
            DRT::Element*               surfaceElement, 
            tetgenio&                   out 
            );
          
          
        /*!
        \brief Debugging the intersection of XAABB's (DEBUG ONLY)
    
        \param cutterXAABB       (in)    : XAABB of the cutting element
        \param xfemXAABB         (in)    : XAABB of the xfem element
        \param cutterElement    (in)    : cutting element
        \param xfemElement      (in)    : xfem element
        \param noC              (in)    : id of the cutting element
        \param noX              (in)    : id of the xfem element
        */                        
        void debugXAABBIntersection( 
            Epetra_SerialDenseMatrix    cutterXAABB,
            Epetra_SerialDenseMatrix    xfemXAABB,
            DRT::Element*               cutterElement,
            DRT::Element*               xfemElement,
            int                         noC,
            int                         noX
            );
                                    
          
        /*!
        \brief Debugging node within element (DEBUG ONLY)
    
        \param element          (in)    : element
        \param node             (in)    : node
        \param xsi              (in)    : reference coordinates
        \param noE              (in)    : id of the element
        \param noN              (in)    : id of the node
        \param within           (in)    : true if within, false otherwise
        */                                     
        void debugNodeWithinElement(    
            DRT::Element*               element,
            DRT::Node*                  node,
            Epetra_SerialDenseVector&   xsi,
            int                         noE,
            int                         noN,
            bool                        within
            );
                                        
                                        
        /*!
        \brief Debugging curve surface intersection (DEBUG ONLY)
    
        \param surfaceElement   (in)    : surfaceElement
        \param lineElement      (in)    : lineElement 
        \param xsi              (in)    : reference coordinates
        \param noSE             (in)    : id of the surface element
        \param noLE             (in)    : id of the line element
        \param within           (in)    : true if within, false otherwise
        */                                     
        void debugCurveSurfaceIntersection(     
            DRT::Element*               surfaceElement,    
            DRT::Element*               lineElement,
            Epetra_SerialDenseVector&   xsi,
            int                         noSE,
            int                         noLE,
            bool                        within
            );
    
      
        /*!
        \brief Debugging tetgen data structure (DEBUG ONLY)
        
        \param element                  (in)        : element
        \param pointList                (in)        : list of points
        \param segmentList              (in)        : list of segments
        \param triangleList             (in)        : list of triangle facets
        */                             
        void debugTetgenDataStructure(  
            DRT::Element*               element,
            vector< InterfacePoint >&   pointList,
            vector< vector<int> >&      surfacePointList,
            vector< vector<int> >&      segmentList,
            vector< vector<int> >&      triangleList
            );  
                                        
          
        /*!
        \brief 	Debugging tetgen data output (DEBUG ONLY)
        		This method prints out the piecewise linear complex (PLC)
        		and the mesh of tetrahedrons for the list of specified
        		xfem elements. Output-files can be viewed with tetview
        		
    
        \param in			(in)       : input tetgen data structure
        \param out        	(in)       : output tetgen data structure
        \param element     	(in)       : xfem-element
        \param elementIds  	(in)       : vector of element-ids
        */                                   
       	void debugTetgenOutput(	
            tetgenio& in,
			tetgenio& out, 
		    DRT::Element* element,
            vector<int>& elementIds
            );
														
						
		/*!
        \brief Debugging integration cells (DEBUG ONLY)
    
        \param integrationcellList                (out)       : list of integration cells
        */                                  
       	void debugIntegrationcells(
            map< int, vector <Integrationcell> >&	integrationcellList,
            int id
            );  
            
       
        /*!
        \brief  computes unique coordinates for a certain region in the 
                xfem element. It is copied into the data structure of Tetgen
                   
        \param xfemElement           (in)        : xfem element
        \param cutterElement         (in)        : cutter element
        \param regioncoordinates     (out)       : coordinates of the region points
        */    
        void computeRegionCoordinates(  
            DRT::Element*  xfemElement,
            DRT::Element*  cutterElement,
            double* regionCoordinates
            );     
                                               
};


};


#endif  // #ifndef INTERSECTION_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_XFEM


/*!
\file intersection.H

\brief collection of intersection tools

<pre>
Maintainer: Ursula Mayer
</pre>
*/

#ifdef CCADISCRET
#ifndef INTERSECTION_H
#define INTERSECTION_H


#include "../drt_xfem/intersection_math.H"
#include "../drt_xfem/integrationcell.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_element.H"
#include "../io/io_drt.H"
#include "../io/gmsh.H"
#include "tetgen.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RefCountPtr.hpp"

 
#ifdef PARALLEL
#include <mpi.h>
#endif


#if defined(__cplusplus)
extern "C"
{
#endif
#ifdef QHULL
#include <qhull/qhull.h>
#include <qhull/qset.h>
#endif
#if defined(__cplusplus)
}
#endif

using namespace std;
using namespace DRT;
using namespace DRT::Utils;
using namespace Teuchos;

namespace XFEM
{


enum pointType {INTERNAL, SURFACE, LINE, NODE, INTERSECTION };



typedef struct
{
  int           nsurf;
  int           surfaces[3];    /*!< adjacent surfaces */
  double        coord[3];       /*!< point coordinates */
  pointType     pType;          /*!< point type */
} InterfacePoint;



class Intersection
{

       public:

                                                       
        /*!
        \brief  computes the intersection of two discretizations.
        
        This methods localizes the interface of the
        cutter discretization within the xfem discretization. Each intersected xfem element
        is tetrahedralized with help of a Constrained Delaunay Tetrahedralization. 
        The output consists of a list of qadratic tetrahedral volume integration cells and 
        a list of quadratic triangular boundary integration cells 
        representing the curved interface.
               
    
        \param xfemdis              (in)    : xfem discretization
        \param cutterdis            (in)    : cutter discretization
        \param domainintcells       (out)   : list of domainintegrationcells for each intersected element
        \param boundaryintcells     (out)   : list of boundaryintegrationcells for each intersected element
        */            
        void computeIntersection(   
            const RefCountPtr<DRT::Discretization>        	xfemdis, 
            const RefCountPtr<DRT::Discretization>        	cutterdis, 
            map< int, DomainIntCells >&						domainintcells,
            map< int, BoundaryIntCells >& 					boundaryintcells
            );
		
        
        
        private:
        
        //! offset in facet marker list
        static const int                    facetMarkerOffset_ = 11;
              
        //! discretization type of the xfem element
        DRT::Element::DiscretizationType    xfemDistype_;      
        
        //! discretization type of the previously handled xfem element
        DRT::Element::DiscretizationType    xfemOldDistype_;      
        
        //! discretization type of the cutter element
        DRT::Element::DiscretizationType    cutterDistype_;
            
        //! number of surfaces of the xfem element
        int                                 numXFEMSurfaces_;       
        
        //! number of lines of the xfem element
        int                                 numXFEMLines_;
            
        //! number of corner nodes of the xfem element
        int                                 numXFEMCornerNodes_;
        
        //! number of internal points
        int                                 numInternalPoints_;
    
        //! number of surface points
        int                                 numBoundaryPoints_;
        
        //! counts missed points during the recovery
        int                                 countMissedPoints_;
        
         //! vector of nodes for each surface of a single element
        vector< vector <int> >              eleNumberingSurfaces_;
        
        //! vector of adjacent surfaces for each node of a single element
        vector< vector <int> >              eleNodesSurfaces_;
        
        //! vector of adjacent surfaces for each line of a single element
        vector< vector <int> >              eleLinesSurfaces_;
        
        //! nodal element coordinate
        vector< vector <double> >           eleRefCoordinates_;
        
        //! list of interface points to be stored in the Tetgen data structure
        vector< InterfacePoint >            pointList_;
    
        //! list of segments to be stored in the Tetgen data structure
        vector< vector<int> >               segmentList_;                 
    
        //! list of single points on a xfem surface to be stored in the Tetgen data structure
        vector< vector<int> >               surfacePointList_;                
    
        //! list of triangles to be stored in the Tetgen data structure
        vector< vector<int> >               triangleList_;
       
        //! list of intersecting cutter elements
        vector< DRT::Element* >             intersectingCutterElements_; 
        
        //! list of face markers numbered in the order of the intersected surface elments
        vector<int>                         faceMarker_; 
        
        
        /*!
        \brief initializes the the intersection computation
    
        \param xfemdis   (in) : xfem discretization
        */
        void initializeIntersection(
            const RefCountPtr<DRT::Discretization>  xfemdis
            );
        
        /*!
        \brief initializes the private members of the xfem element
    
        \param xfemId        (in) : id of the current xfem element
        \param xfemElement   (in) : xfem element
        */
        void initializeXFEM(  
            int             xfemId,
            DRT::Element*   xfemElement
            );
            
            
        /*!
        \brief initializes the private members of the cutter element
    
        \param cutterElement   (in) : cutter element
        */
        void initializeCutter(  
            DRT::Element* cutterElement
            );
           		
        
        /*!
        \brief transforms a point in element coordinates to a point
               in current coordinates with respect to a given element
    
        \param element          (in)        : element
        \param xsi              (in/out)    : vector of element coordinates/
                                              vector of current coordinates
        */
        void elementToCurrentCoordinates( 
            DRT::Element*               element, 
            Epetra_SerialDenseVector&   xsi
            );
 
 
        /*!
        \brief transforms a point in element coordinates to a point
               in current coordinates with respect to a given surfaceElement
               based on given nodal coordinates
    
        \param surfacElement    (in)        : surfaceElement
        \param xsi              (in/out)    : vector of element coordinates/
                                              vector of current coordinates
        \param surfaceNodes     (in)        : nodal coordinates of the surface element
        */
        void elementToCurrentCoordinates(   
            DRT::Element*                               surfaceElement, 
            Epetra_SerialDenseVector&                   xsi,
            const vector <Epetra_SerialDenseVector>&    surfaceNodes
            );
    
        
        /*!
        \brief transforms a point in element coordinates to a point
               in current coordinates with respect to a given plane
               based on given nodal coordinates
    
        \param xsi          (in/out)    : vector of element coordinates/
                                          vector of current coordinates
        \param plane        (in)        : quad4 surface element
        */
        void elementToCurrentCoordinates(  
            Epetra_SerialDenseVector&                   xsi,
            const vector<Epetra_SerialDenseVector>&     plane
            );


        /*!
        \brief transforms a point in current coordinates to a point
               in element coordinates with respect to a given element
    
        \param element          (in)        : element
        \param xsi              (in/out)    : vector of current coordinates/
                                              vector of element coordinates
        */
        void currentToElementCoordinates( 
            DRT::Element*               element,                     
            Epetra_SerialDenseVector&   xsi
            );
 
        
        /*!
       	\brief transforms a point in current coordinates to a point
               in element coordinates with respect to a given element
            
        \param element          (in)        : element
        \param xsi              (in/out)    : vector of current coordinates/
                                              vector of element coordinates
        */
        void currentToElementCoordinates(   
        	DRT::Element*               element, 
        	vector<double>&   			xsiVector
          	);
    
        
        /*!
       	\brief transforms a point in current coordinates to a point
               in element coordinates with respect to a given element
            
        \param element          (in)        : element
        \param xsi              (in/out)    : vector of current coordinates/
                                              vector of element coordinates
        */
        void elementToCurrentCoordinates(   
            DRT::Element* 	element, 
            vector<double>& xsi
            );
        
        
        /*!
        \brief computes a global numbering for cutter elements, which are boundary
               elements with xfem conditions
    
        \param cutterdis            (in)    : cutter discretization
        \param xfemConditions       (in)    : list of boundary elements with xfem conditions
        \param conditionEleCount    (out)   : global element counter
        */
    	void adjustCutterElementNumbering( 
			const RefCountPtr<DRT::Discretization>& 	cutterdis,  
			vector< DRT::Condition * >&   				xfemConditions,
			vector<int>& 								conditionEleCount
			);
    
    
        /*!
        \brief packs data for parallelization
    
        \param cutterdis            (in)        : cutter discretization
        \param conditionSend        (in)        : list of boundary elements with xfem conditions
        \param lengthSend           (in)        : size of packed data
        \param nodeSetSizeSend      (in)        : size of node set
        \param nodeVectorSend       (in)        : node vector
        \param cutterDataSend       (out)       : cutter data
        */
        void packData(
            const RefCountPtr<DRT::Discretization>      cutterdis,
            vector<int>&                                conditionSend, 
            vector<int>&                                lengthSend, 
            int&                                        nodeSetSizeSend,  
            vector<int>&                                nodeVectorSend,
            vector<char>&                               cutterDataSend 
            );
        
        
        /*!
        \brief unpacks nodes after a MPI send
    
        \param index            (in)        : start index in cutter data to read the nodes
        \param cutterDataRecv   (in)        : received cutter data
        \param nodeVectorRecv   (in)        : received node vector
        \param nodeMap          (out)       : node map
        */
        void unpackNodes(
            int                                     index, 
            vector<char>&                           cutterDataRecv, 
            vector<int>&                            nodeVectorRecv, 
            map< int, RefCountPtr<DRT::Node> >&     nodeMap 
            );
    
    
        /*!
        \brief sends and recieves cutter element data
    
        \param xfemConditions       (in)        : list of cutter elements with xfem conditions
        \param conditionEleCount    (in)        : global cutter element count
        \param cutterElementMap     (in)        : map of cutter elements
        \param cutterNodeMap        (out)       : map of nodes for cutter elements
        \param xfemCutterIdMap      (out)       : map of cutter ids
        \param xfemdis              (in)        : xfem discretization
        \param cutterdis            (in)        : cutter discretization
        */
    	void getCutterElementsInParallel(	
			vector< DRT::Condition * >&      		xfemConditions,
			vector<int>&							conditionEleCount, 
			map< int, vector< DRT::Element* > >&   	cutterElementMap,
			map< int, RefCountPtr<DRT::Node> >&  	cutterNodeMap,
			map< int, set<int> >& 					xfemCutterIdMap,
			const RefCountPtr<DRT::Discretization>& xfemdis,
		   	const RefCountPtr<DRT::Discretization>& cutterdis  
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
        \brief compares two points (overloaded method)
    
        \param point1       (in)    : first point   (DRT::Node*)
        \param point2       (in)    : second point  (DRT::Node*)
        \return true if both points equal each other, false otherwise
        */                      
        bool comparePoints( 
            const DRT::Node*     point1,
            const DRT::Node*     point2
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
        \brief Computes a rough overestimating extended  
               axis-aligned bounding box for an element (XAABB)
    
        \param  element          (in)      : element (DRT::Element*)
        \return extended axis-aligned bounding box  (XAABB) for an element
        */
        Epetra_SerialDenseMatrix computeFastXAABB( 
            DRT::Element*   element
            );
        
                                                
        /*!
        \brief Checks if a node is within an 
               extended axis-aligned bounding box (XAABB)
    
        \param node (in)        : node to be checked
        \param XAABB (in)       : extended axis-aligned bouning box
        \return true if node is within the XAABB or false otherwise
        */
        bool isNodeWithinXAABB(  
            const std::vector<double>&         node,
            const Epetra_SerialDenseMatrix&    XAABB
            );
        
         
        bool isLineWithinXAABB(    
            const std::vector<double>&         node1,
            const std::vector<double>&         node2,
            const Epetra_SerialDenseMatrix&    XAABB
            );
        
        
        /*!
        \brief checks if two extended axis aligned bounding boxes intersect
    
        \param cutterXAABB (in)         : XAABB of the cutting element
        \param xfemXAABB   (in)         : XAABB of the xfem element
        \return true if the XAABB's intersect or false otherwise
        */
        bool intersectionOfXAABB(    
            const Epetra_SerialDenseMatrix&    cutterXAABB,
            const Epetra_SerialDenseMatrix&    xfemXAABB
            );
        
        	       
        /*!
        \brief collects points belonging to an interface lying inside xfem elements         
               The nonlinear system of equation is solved with help of the Newton-method.
    
        \param xfemElement          (in)        : xfem element 
        \param cutterElement        (in)        : cutter element 
        \param node                 (in)        : node in current coordinates (Debug)
        \param interfacePoints      (in/out)    : vector of interface points
        \param elemId               (in)        : xfem element id  (Debug)
        \param nodeId               (in)        : node id
        \return true if the node lies within the element, false otherwise
        */                                      
        bool collectInternalPoints( 
            DRT::Element*                   xfemElement,
            DRT::Element*                   cutterElement,
            DRT::Node*                      node,
            std::vector<InterfacePoint>&    interfacePoints,
            const int                       elemId,
            const int                       nodeId
            );
        
                       							
        /*!
        \brief updates the system matrix at the corresponding element coordinates for the 
               computation if a node in current coordinates lies within an element 
    
        \param dim               (in)       : dimension of the problem
        \param A                 (out)      : system matrix
        \param xsi               (in)       : vector of element coordinates
        \param element           (in)       : element 
        */                                                                                                      
        void updateAForNWE( 
            const int                   dim,
            Epetra_SerialDenseMatrix&   A,
            Epetra_SerialDenseVector&   xsi,
            DRT::Element*               element
            );
        
          
        /*!
        \brief updates the rhs at the corresponding element coordinates for the 
               computation whether a node in current coordinates lies within an element 
    
        \param dim               (in)       : dimension of the problem
        \param b                 (out)      : right-hand-side       
        \param xsi               (in)       : vector of element coordinates
        \param x                 (in)       : node in current coordinates
        \param element           (in)       : element
        */                               
        void updateRHSForNWE(  
            const int                           dim, 
            Epetra_SerialDenseVector&           b,
            Epetra_SerialDenseVector&           xsi,
            const Epetra_SerialDenseVector&     x,
            DRT::Element*                       element
            );         
        
            
		/*!
        \brief checks if a node in current coordinates lies within a certain element           
               The nonlinear system of equation is solved with help of the Newton-method.
    
        \param element              (in)        : element 
        \param x                    (in)        : node in current coordinates
        \param xsi                  (in)        : node in element coordinates
        */	
		bool checkNodeWithinElement(	
            DRT::Element*                       element,
            const Epetra_SerialDenseVector&     x,
            Epetra_SerialDenseVector&           xsi
            );	       
        
        
        /*!
        \brief  checks if a node that lies within an element and if it 
                lies on one of its surfaces, lines or nodes          .
    
        \param xsi                  (in)        : node in element coordinates
        \param ip                   (in)        : interface point
        \return true if the node is lying on a boundary false otherwise
        */  
        bool checkIfOnBoundary(   
            Epetra_SerialDenseVector&       xsi, 
            InterfacePoint&                 ip);
        
              
        /*!
        \brief collects all intersection points between a line and a surface
            
        \param surfaceElement           (in)        : surface element
        \param lineElement              (in)        : line element
        \param interfacePointList       (out)       : vector of interface points
        \param surfaceId                (in)        : surface element id
        \param lineId                   (in)        : line element id
        \param lines                    (in)        : if lines = true
        */                                                      
        bool collectIntersectionPoints( 
            DRT::Element*                   surfaceElement,
            DRT::Element*                   lineElement,
            std::vector<InterfacePoint>&    interfacePoints,
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
            Epetra_SerialDenseVector&       xsiOld,
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
        \param xsi              (in/out)    : starting value/vector of element coordinates
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
            \brief solves a singular system of equations
    
        \param xsi              (in/out)    : vector of element domain coordinates
        \param lineElement      (in)        : line element
        \param surfaceElement   (in)        : surface element
        return true if resulting system is singular , false otherwise
        */                                               
        bool computeSingularCSI(
            Epetra_SerialDenseVector&   xsi,
            DRT::Element*               lineElement,
            DRT::Element*               surfaceElement
            );
            
                 
		/*!
	  	\brief updates the systemmatrix at the corresponding element coordinates 
               for the computation of curve surface intersections
	
		\param A  				 (out)		: system matrix
		\param xsi				 (in)	    : vector of element coordinates
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
        \brief updates the rhs at the corresponding element coordinates 
               for the computation of curve surface intersections
   
        \param b                 (out)      : right-hand-side
        \param xsi               (in)       : vector of element coordinates
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
        \brief  creates new limits for the recursive intersection computation
    
        \param xsi              (in)        : intersection point
        \param upLimit          (in)        : upper limits
        \param loLimit          (in)        : lower limits
        \param upperLimits      (out)       : list of upper limits of new ranges
        \param lowerLimits      (out)       : list of lower limits of new ranges
        */              
        void createNewLimits(
            Epetra_SerialDenseVector&               xsi, 
            Epetra_SerialDenseVector&               upLimit,
            Epetra_SerialDenseVector&               loLimit,
            vector< Epetra_SerialDenseVector >&     upperLimits, 
            vector< Epetra_SerialDenseVector >&     lowerLimits
            );		   							
        
        
        /*!
        \brief  computes the convex hull of a set of points
    
        \param surfaceElement           (in)        : surface element
        \param interfacePoints          (in)        : list of interface points
        */           
        void computeConvexHull( 
            DRT::Element*               xfemElement,
            DRT::Element*               surfaceElement,
            vector<InterfacePoint>&     interfacePoints
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
    
        \param element              (in)    : xfem element
        \param cutterElement        (in)    : cutter element
        \param domainintcells       (out)   : list of domainintegrationcells for each intersected element
        \param boundaryintcells     (out)   : list of boundaryintegrationcells for each intersected element
        */       
        void computeCDT(    
            DRT::Element*                           element,
            DRT::Element*                           cutterElement,
            map< int, DomainIntCells >&             domainintcells,
            map< int, BoundaryIntCells >&           boundaryintcells
            );
    
    
        /*!
        \brief  adds the corner nodes of an xfem element to the pointlist
    
        */          
        void startPointList();
            
            
        /*!
        \brief  stores an interface point in a point list for the 
                Constrained Delaunay Tetrahedralization (CDT) with Tetgen
    
        \param point            (in)        :   coordinates of the point to be stored 
        \param interfacePoints  (in/out)    :   vector of interface points
        \param positions        (in/out)    :   positions
 
        */                          
        void storePoint(    
            const vector<double>&       point, 
            vector<InterfacePoint>&     interfacePoints, 
            vector<int>&                positions
            );
  
         
        /*!
        \brief  stores a single point in a list of surface points if no 
                segment lies in that surface
  
        \param interfacePoints      (in)        :   list of interface points
        
        */      
        void storeSurfacePoints(    
            vector<InterfacePoint>&     interfacePoints
            );
          
            
        /*!
        \brief  computes the midpoint of a collection of points
    
        \param interfacePoints            (in)        : vector of interface points
        \return returns the midpoint (type InterfacePoint)
        */                                 
        InterfacePoint computeMidpoint( 
            const vector<InterfacePoint>&     interfacePoints
            );       
            
                             
        /*!
        \brief  stores a segment in a segment list for the 
                Constrained Delaunay Tetrahedralization (CDT) with Tetgen
    
        \param positions            (in)        :   positions
        */            
        void storeSegments( 
            const vector<int>&              positions
            );
          
                               
        /*!
        \brief  stores a triangle facet in a triangle list for the 
                Constrained Delaunay Tetrahedralization (CDT) with Tetgen
    
        \param positions            (in)        :   positions
        */                         
        void storeTriangles(    
            const vector<int>               positions
            );
         

        /*!
        \brief  stores the cutter element intersecting 
                with an xfem element in a list
    
        \param surfaceElement   (in)    :   surface element
        */                  
        void storeIntersectedCutterElement(
            DRT::Element* surfaceElement
            );  
         
  
        /*!
        \brief  computes the recovery of the curved interface
                for one xfem element after a Contrained Delaunay 
                Tetrahedralization
    
        \param xfemElement              (in)        :   xfem element
        \param out                  (in/out)    :   output data structure of TetGen
        */      
        void recoverCurvedInterface(
            DRT::Element*                   xfemElement, 
            map< int, BoundaryIntCells >&   boundaryintcells,
            tetgenio&                       out
            );
           
            
        /*!
        \brief  checks if Steiner points are inserted additionally. They are lifted up such they
                are lying in a cutter element.
    
        \param xfemElement          (out)    :   xfem element
        \param out                  (out)    :   output data structure of TetGen
        */     
        void liftAllSteinerPoints(
            DRT::Element*                                   xfemElement,
            tetgenio&                                       out
            
            );      
           
           
        /*!
        \brief  collects for each Steiner point its adjacent triangular faces
    
        \param adjacentFacesList        (out)    :   list of adjacent faces
        \param adjacentFacemarkerList   (out)    :   list of face marker of each adjacent face
        \param out                      (in)     :   output data structure of TetGen
        */     
        void locateSteinerPoints(
            vector< vector <int> >&     adjacentFacesList,
            vector< vector <int> >&     adjacentFacemarkerList,
            tetgenio&                   out
            );
         
         
        /*!
        \brief  checks if the Steiner points will lie in a cutter element or on
                one of the edges of the cutter element
    
        \param steinerIndex             (in)     :   index of the Steiner point
        \param lineIndex                (out)    :   line index
        \param cutterIndex              (out)    :   index of the cutter element
        \param adjacentFacesList        (in)     :   list of adjacent faces
        \param adjacentFacemarkerList   (in)     :   list of face marker of each adjacent face
        \param edgePoint                (out)    :   point lying on the same edge as the Steiner point  
        \param oppositePoint            (out)    :   point lying opposite on the edge of the Steiner point  
        \param xfemElement              (in)     :   xfem element  
        \param out                      (in)     :   output data structure of TetGen
        \return true the point lies in the cutter element and false if on an edge
        */     
        int decideSteinerCase(  
            int                         steinerIndex,
            int&                        lineIndex, 
            int&                        cutterIndex,  
            vector< vector <int> >&     adjacentFacesList,
            vector< vector <int> >&     adjacentFacemarkerList,
            Epetra_SerialDenseVector&   edgePoint,
            Epetra_SerialDenseVector&   oppositePoint,
            DRT::Element*               xfemElement, 
            tetgenio&                   out
            );   
    
    
        /*!
        \brief  lifts the Steiner point into the cutter element
    
        \param steinerIndex             (in)     :   index of the Steiner point
        \param adjacentFacesList        (in)     :   list of adjacent faces
        \param adjacentFacemarkerList   (in)     :   list of face marker of each adjacent face
        \param xfemElement              (in)     :   xfem element  
        \param out                      (in)     :   output data structure of TetGen
        */     
        void liftSteinerPointOnSurface(
            int                         steinerIndex,
            vector< vector <int> >&     adjacentFacesList,
            vector< vector <int> >&     adjacentFacemarkerList,
            DRT::Element*               xfemElement, 
            tetgenio&                   out
            ); 
        
        
        /*!
        \brief  lifts the Steiner point onto the edge of the cutter element
                specified by the line index
    
        \param steinerIndex             (in)     :   index of the Steiner point
        \param lineIndex                (in)     :   line index
        \param cutterIndex              (in)     :   index of the cutter element
        \param edgePoint                (in)     :   point lying on the same edge as the Steiner point  
        \param oppositePoint            (in)     :   point lying opposite on the edge of the Steiner point 
        \param adjacentFacesList        (in)     :   list of adjacent faces 
        \param xfemElement              (in)     :   xfem element  
        \param out                      (in)     :   output data structure of TetGen
        */     
        void liftSteinerPointOnEdge(
            int                         steinerIndex,
            int                         lineIndex,
            int                         cutterIndex,   
            Epetra_SerialDenseVector&   edgePoint,
            Epetra_SerialDenseVector&   oppositePoint,
            vector< vector <int> >&     adjacentFacesList,
            DRT::Element*               xfemElement, 
            tetgenio&                   out
            );
          
        
        /*!
        \brief  lifts the Steiner point on an edge that lies on the xfem boundary
    
        \param steinerIndex             (in)     :   index of the Steiner point 
        \param adjacentFacesList        (in)     :   list of adjacent faces 
        \param adjacentFacemarkerList   (in)     :   list of adjacent face marker
        \param xfemElement              (in)     :   xfem element  
        \param out                      (in)     :   output data structure of TetGen
        */       
        void liftSteinerPointOnBoundary( 
            int                         steinerIndex, 
            vector< vector <int> >&     adjacentFacesList, 
            vector< vector <int> >&     adjacentFacemarkerList, 
            DRT::Element*               xfemElement, 
            tetgenio&                   out
            );
            
             
        /*!
        \brief  gathers information of a tetrahedron adjacent to an 
                interface 
                
        \param tetIndex             (in)    : tetrahedron index of the tetrahedron list
        \param faceIndex            (in)    : index of the face        
        \param tetraCornerIndices   (out)   : index list of tetrahedron corners
        \param order                (out)   : order the tetrahedron corner points appear
        
        \param out                  (in)    : output data structure of TetGen
        */
        void getTetrahedronInformation(   
            const int           tetIndex,
            const int           faceIndex,  
            vector<int>&        tetraCornerIndices,
            vector<int>&        order,
            const tetgenio&     out
            );          
        
        
        /*!
        \brief  gathers information of a tetrahedron adjacent to an 
                interface 
                
        \param tetraCornerIndices   (in)    : index list of tetrahedron corners
        \param tetraCornerIndices   (out)   : list of tetrahedron corners nodes in the current 
                                              configuration of the xfem element
        \param xfemElement          (in)    : xfem element
        \param out                  (in)    : output data structure of TetGen
        */
        void getTetrahedronNodes(
            vector<Epetra_SerialDenseVector>&       tetraCornerNodes,
            vector<int>&                            tetraCornerIndices,
            DRT::Element*                           xfemElement,
            const tetgenio&                         out
            ); 
        
                
        /*!
        \brief  lifts a higher order point of an edge of the linearized interface
                onto the curved interface
    
        \param index1                   (in)    :   local index of edge point 1
        \param index2                   (in)    :   local index of edge point 2
        \param faceIndex                (in)    :   index of the face 
        \param faceMarkerIndex          (in)    :   face marker
        \param globalHigherOrderIndex   (in)    :   global index of the higher-order point
        output xfemElement              (in)    :   xfem element
        \param out                      (in)    :   output data structure of TetGen
        */     
        void computeHigherOrderPoint(    
            int                                 index1, 
            int                                 index2, 
            int                                 faceIndex, 
            int                                 faceMarker, 
            int                                 globalHigherOrderIndex, 
            vector<int>&                        tetraCornerIndices, 
            vector<Epetra_SerialDenseVector>&   tetraCornerNodes, 
            DRT::Element*                       xfemElement, 
            tetgenio&                           out
            );
    
    
        /*!
        \brief  returns the other two point indices belonging to a triface 
                with a Steiner point
    
        \param out                  (in)    :   output data structure of TetGen
        \param trifaceIndex         (in)    :   index of the triface
        \param steinerPointIndex    (in)    :   index of the Steiner point
        \return vector of point indices
        */         
        vector<int> getPointIndices(
            tetgenio&   out, 
            int         trifaceIndex, 
            int         steinerPointIndex
            );
            
              
        /*!
        \brief  computes the recovery of the curved interface
                for one xfem element after a Contrained Delaunay 
                Tetrahedralization. 
    
        \param xsi                      (out)   :   intersection point
        \param normal                   (in)    :   nodes of the intersection normal
        \param cutterElement            (in)    :   cutter element
        \param onBoundary               (in)    :   indicates whether the normal lies on the xfem element 
                                                    boundary or not
        \return true if normal intersects with cutter element or false otherwise
        */      
        bool computeRecoveryNormal( 
            Epetra_SerialDenseVector&                   xsi,
            const vector<Epetra_SerialDenseVector>&     normal,
            DRT::Element*                               cutterElement,
            const bool                                  onBoundary
            );
        
        
        /*!
        \brief  updates the system matrix for the intersection computation
                of an average normal and the corresponding cutter element
                to lift a point of the linearized interface onto the curved interface 
    
        \param A                (out)       :   system matrix
        \param xsi              (in)        :   intersection point
        \param normal           (in)        :   nodes of the intersection normal
        \param surfaceElement   (in)        :   surface element
        \param onBoundary       (in)        :   true, if intersection normal lies on the
                                                element boundary
        */      
        void updateAForRCINormal(   
            Epetra_SerialDenseMatrix&                   A,
            const Epetra_SerialDenseVector&             xsi,
            const vector<Epetra_SerialDenseVector>&     normal,
            DRT::Element*                               surfaceElement,
            const bool                                  onBoundary
            );         
         
        
        /*!
        \brief  updates the right-hand-side for the intersection computation
                of an average normal and the corresponding cutter element
                to lift a point of the linearized interface onto the curved interface 
                
        \param b                (out)       :   right-hand-side
        \param xsi              (in)        :   intersection point
        \param normal           (in)        :   nodes of the intersection normal
        \param surfaceElement   (in)        :   surface element
        \param onBoundary       (in)        :   true, if intersection normal lies on the
                                                element boundary
        */                         
        void updateRHSForRCINormal( 
            Epetra_SerialDenseVector&                   b,
            Epetra_SerialDenseVector&                   xsi,    
            const vector<Epetra_SerialDenseVector>&     normal,
            DRT::Element*                               surfaceElement,
            const bool                                  onBoundary
            );           
              
        
        /*!
        \brief  computes the recovery of the curved interface
                for one xfem element after a Contrained Delaunay 
                Tetrahedralization. 
                
        \param lineIndex                (in)    :   line index
        \param xsi                      (out)   :   intersection point
        \param plane                    (in)    :   nodes of the intersection normal
        \param surfaceElement           (in)    :   surface element
        
        \return true if plane intersects with the line element or false otherwise
        */      
        bool computeRecoveryPlane( 
            int&                                        lineIndex,
            Epetra_SerialDenseVector&                   xsi,
            const vector<Epetra_SerialDenseVector>&     plane,
            DRT::Element*                               surfaceElement
            );
            
         
        /*!
        \brief  updates the system matrix for the intersection computation
                of an average normal and the corresponding cutter element
                to lift a point of the linearized interface onto the curved interface 
    
        \param A                (out)       :   system matrix
        \param xsi              (in)        :   intersection point
        \param plane            (in)        :   nodes of the intersection normal
        \param lineElement      (in)        :   line element
        \param surfaceElement   (in)        :   surface element
        */         
        void updateAForRCIPlane(   
            Epetra_SerialDenseMatrix&                   A,
            const Epetra_SerialDenseVector&             xsi,
            const vector<Epetra_SerialDenseVector>&     plane,
            DRT::Element*                               lineElement,
            DRT::Element*                               surfaceElement
            );
            
            
        /*!
        \brief  updates the right-hand-side for the intersection computation
                of the cutter element edge with a plane through the midpoint 
                of the face edge to lift a point of the linearized interface 
                onto the curved interface 
                
        \param b                (out)       :   right-hand-side
        \param xsi              (in)        :   intersection point
        \param plane            (in)        :   nodes of the intersection normal
        \param lineElement      (in)        :   line element
        */                       
        void updateRHSForRCIPlane( 
            Epetra_SerialDenseVector&                   b,
            Epetra_SerialDenseVector&                   xsi,    
            const vector <Epetra_SerialDenseVector>&    plane,
            DRT::Element*                               lineElement
            );
            
                           
        /*!
        \brief  computes the normal to the interface edge of the
                adjacent tetrahedon. In addition the normal lies within the plane
                of the tetrahedron. 
              
        \param onBoundary               (in)    : true if edge lies on the surface of the xfem element, 
                                                  false otherwise
        \param index1                   (in)    : local index 1 of an interface edge
        \param index2                   (in)    : local index 2 of an interface edge
        \param oppositePointIndex       (in)    : gloabal index of the point oopoiste to the edge
        \param globalHigherOrderIndex   (in)    : global higher order index
        \param tetraCornerIndices       (in)    : indices of tetrahedron corners
        \param tetraCornerNodes         (in)    : list of tetrahedron corner nodes in 
                                                  current configuration of the xfem element
        \param plane                    (out)   : nodes of the normal and plane 
        \param xfemElement              (in)    : xfem element
        \param out                      (in)    : output data structure of tetGen
        */    
        void computeIntersectionNormal( 
            const bool                              onBoundary,
            const int                               index1,
            const int                               index2,
            const int                               oppositePointIndex,
            const int                               globalHigherOrderIndex, 
            const vector<int>&                      tetraCornerIndices,
            const vector<Epetra_SerialDenseVector>& tetraCornerNodes,
            vector<Epetra_SerialDenseVector>&       plane,
            DRT::Element*                           xfemElement,
            const tetgenio&                         out
            );
            

        /*!
        \brief  computes the normal to the interface edge of two
                adjacent triangular faces passing through the midpoint of that edge.             
              
        \param index1                   (in)    : global index 1 of an interface edge
        \param index2                   (in)    : global index 2 of an interface edge
        \param faceIndex                (in)    : index of the face
        \param adjacentFaceIndex        (in)    : index of the adjacent face
        \param plane                    (out)   : nodes of the normal and plane 
        \param xfemElement              (in)    : xfem element
        \param out                      (in)    : output data structure of tetGen
        */    
        void computeIntersectionNormal(  
            int                                 index1,
            int                                 index2, 
            int                                 faceIndex,
            int                                 adjacentFaceIndex,
            int                                 globalHigherOrderIndex,
            vector<Epetra_SerialDenseVector>&   plane,
            DRT::Element*                       xfemElement,
            tetgenio&                           out
            );    
   
   
        /*!
        \brief  computes the normal to the interface edge of a triangular
                faces on the xfem boundary passing through the Steiner point of that edge.             
              
        \param steinerIndex             (in)    : index of Steiner point
        \param edgeIndex                (in)    : index of the other edge point
        \param oppositeIndex            (in)    : index of the point on the boundary triangle opposite the edge
        \param plane                    (out)   : nodes of the normal and plane 
        \param xfemElement              (in)    : xfem element
        \param out                      (in)    : output data structure of tetGen
        */    
        void computeIntersectionNormal( 
            const int                               steinerIndex,
            const int                               edgeIndex,
            const int                               oppositeIndex,
            vector<Epetra_SerialDenseVector>&       plane,
            DRT::Element*                           xfemElement,
            const tetgenio&                         out
            );
            
            
        /*!
        \brief  compute the midpoint of a line described by two points
                
        \param p1       (in)    : point 1
        \param p2       (in)    : point 2
        \return mid point of a line
        */    
        Epetra_SerialDenseVector computeLineMidpoint(  
            const Epetra_SerialDenseVector& p1,
            const Epetra_SerialDenseVector& p2
            );      
           
           
        /*!
        \brief  finds the face marker of an adjacent facet for a given edge 
                of another facet
                
        \param edgeIndex1               (in)   : index of edge 1
        \param edgeIndex2               (in)   : index of edge 2
        \param faceMarkerIndex          (in)   : index of the given facet
        \param out                      (in)   : output data structure of tetGen
        \return the index of the face marker of the adjacent facet
        */         
        void findAdjacentFace(
            const int       edgeIndex1, 
            const int       edgeIndex2, 
            const int       faceMarker,
            int&            adjacentFaceMarker,
            int             faceIndex,
            int&            adjacentFaceIndex,
            tetgenio&       out
            );   
            
            
        /*!
        \brief  searches for the common edge of two adjacent facets
                
        \param faceIndex1               (in)    : index of facet 1
        \param faceIndex2               (in)    : index of facet 2
        \param adjacentFacesList        (in)    : list of adjacent facets
        \param edgepoint                (in)    : edge point
        \param oppositepoint            (in)    : opposite point
        \param out                      (in)    : output data structure of tetGen
        \return true if common edge was found, false otherwise
        */          
        bool findCommonFaceEdge( 
            const int                           faceIndex1, 
            const int                           faceIndex2, 
            vector<int>&                        adjacentFacesList,
            Epetra_SerialDenseVector&           edgePoint,
            Epetra_SerialDenseVector&           oppositePoint,
            tetgenio&                           out
            );
            
      
        /*!
        \brief  searches for the common line element of two 
                adjacent cutter element
                
        \param faceIndex1               (in)    : index of facet 1
        \param faceIndex2               (in)    : index of facet 2
        \param lineIndex                (out)   : index of the common line element
        \param surfaceIndex             (out)   : index of the cutter element
        \return true if common line element was found, false otherwise
        */          
        bool findCommonCutterLine(  
            const int                                       faceIndex1, 
            const int                                       faceIndex2,
            int&                                            lineIndex,
            int&                                            cutterIndex
            );
        
        
       /*!
        \brief  returns the index of a point lying opposite to a face edge
                
        \param edgeNode1                (in)   : edge node 1 of the facet edge
        \param edgeNode2                (in)   : edge node 2 of the facet edge
        \param adjacentFaceIndex        (in)   : list of adjacent face indices
        \param out                      (in)   : output data structure of tetGen
        \return the line index of the intersection line element of the given cutter element
        */              
       int findEdgeOppositeIndex( 
            int                                 edgeIndex1,
            int                                 edgeIndex2, 
            int                                 adjacentFaceIndex,
            tetgenio&                           out
            );
    
             
        /*!
        \brief  for the recovery computation of a higher-order node
                by plane - line element intersection
                this method finds the actually intersecting line
                for a given cutter element corresponding to a given facet edge
                
        \param xfemElement              (in)   : xfem element
        \param cutterElement            (in)   : cutter element
        \param edgeNode1                (in)   : edge node 1 of the facet edge
        \param edgeNode2                (in)   : edge node 2 of the facet edge
        \return the line index of the intersection line element of the given cutter element
        */             
        int findIntersectingSurfaceEdge(
            DRT::Element*                       xfemElement,
            DRT::Element*                       cutterElement,
            Epetra_SerialDenseVector&           edgeNode1,
            Epetra_SerialDenseVector&           edgeNode2
            );
    
    
        /*!
        \brief  stores the higher order node in the pointlist of the output data
                structure of tetGen out
                
        \param normal                   (in)    : true if itersection of a normal with a surface was computed
                                                  false otherwise   
        \param globalHigherOrderIndex   (in)    : lineIndex
        \param lineIndex                (in)    : lineIndex
        \param xsi                      (in)    : coordinates of the higher order node 
        \param surface element          (in)    : surface element
        \param xfem element             (in)    : xfem element
        \param out                      (in)    : output data structure of tetGen
        */                              
        void storeHigherOrderNode( 
            const bool                                  normal,
            const int                                   globalHigherOrderIndex,
            const int                                   lineIndex, 
            Epetra_SerialDenseVector&                   xsi, 
            DRT::Element*                               surfaceElement, 
            DRT::Element*                               xfemElement, 
            tetgenio&                                   out
            );
        
        
        /*!
        \brief  stores domain integration cells
    
        \param xfemElement          (in)        :   xfem element
        \param domainintcells       (out)       :   domain integration cells
        \param out                  (in)        :   tetgenio output data structure 
        */    
        void storeDomainIntCells(
            DRT::Element*                   xfemElement,
            map< int, DomainIntCells >&     domainintcells,
            tetgenio&                       out
            ); 
             
              
        /*!
        \brief  stores boundary integration cells
    
        \param trifaceIndex             (in)        :   index of the triangular face
        \param cornerIndex              (in)        :   global index of corner node of an triangle edge
        \param globalHigherOrderIndex   (in)        :   global higher-order index the midpoint of an triangle edge
        \param faceMarker   			(in)        :   id of cutter element
        \param domainCoord            	(out)       :   point coordinates of the triangular face in xfem configuration
        \param boundaryCoord            (out)       :   point coordinates of the triangular face in cutter configuration
        \param xfemElement				(in)		:	xfem element
        \param out                      (in)        :   tetgenio output data structure 
        */      
        void storeBoundaryIntCells(
        	const int                         		trifaceIndex,
        	const int                         		cornerIndex, 
        	const int                         		globalHigherOrderIndex, 
        	const int                         		faceMarker, 
		    std::vector< std::vector<double> >&   	domainCoord, 
		    std::vector< std::vector<double> >&   	boundaryCoord, 
		    DRT::Element*							xfemElement,
		    tetgenio&                   			out
            );
            
                    
        /*!
        \brief Debugging the intersection of XAABB's (DEBUG ONLY)
    
        \param cutterXAABB       (in)    : XAABB of the cutting element
        \param xfemXAABB         (in)    : XAABB of the xfem element
        \param cutterElement     (in)    : cutting element
        \param xfemElement       (in)    : xfem element
        \param noC               (in)    : id of the cutting element
        \param noX               (in)    : id of the xfem element
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
        \param xsi              (in)    : element coordinates
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
        \param xsi              (in)    : element coordinates
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
        
        */                             
        void debugTetgenDataStructure(  
            DRT::Element*               element
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
    
        \param domainintcells                (out)       : list of integration cells
        */                                  
       	void debugDomainIntCells(
            map< int, DomainIntCells >&	domainintcells,
            int id
            );  
            

        /*!
        \brief  computes unique coordinates for a certain region in the 
                xfem element. It is copied into the data structure of Tetgen
                   
        \param index            (in)        : index of the xfem element
        \param out              (in)        : tetgenio output data structure
        */        
        void printTetViewOutput(
            int             index,
            tetgenio&       out
            ); 
            
            
        /*!
        \brief  computes unique coordinates for a certain region in the 
                xfem element. It is copied into the data structure of Tetgen
                   
        \param index            (in)        : index of the xfem element
        \param out              (in)        : tetgenio output data structure
        */           
        void printTetViewOutputPLC(
            DRT::Element*   xfemElement,
            int             index,
            tetgenio&       in
            );   
        
        
       void debugFaceMarker(
        	const int 						eleId,
        	tetgenio&						out
        	);
       
       void debugXFEMConditions(
       	const RefCountPtr<DRT::Discretization>  cutterdis
       	);
       
       void debugIntersection(
       		DRT::Element*			xfemElement,
       		vector <DRT::Element*> 	cutterElements
       		);
       
       void debugXAABBs(
        const int							id,
       	const Epetra_SerialDenseMatrix&     cutterXAABB, 
       	const Epetra_SerialDenseMatrix&     xfemXAABB
       	);
                                               
};


};


#endif  // #ifndef INTERSECTION_H
#endif  // #ifdef CCADISCRET

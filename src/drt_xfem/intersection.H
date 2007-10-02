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
#include "../drt_f3/fluid3.H"
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

 //{INTERNAL, LINE, SURFACE, NODE };
enum pointType {INTERNAL, SURFACE, LINE, NODE, INTERSECTION };



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
        
        //! tolerance 1e-14 used for residuals in all Newton-methods
        static const double                 TOL14_ = 1e-14;
        
        //! tolerance 1e-7 
        static const double                 TOL7_ = 1e-7;
        
        //! offset in facet marker list
        static const int                    facetMarkerOffset_ = 11;
              
        //! discretization type of the xfem element
        DRT::Element::DiscretizationType    xfemDistype_;      
        
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
        
         //! vector of nodes for each surface
        vector< vector <int> >              eleNumberingSurfaces_;
        
        //! vector of adjacent surfaces for each node
        vector< vector <int> >              eleNodesSurfaces_;
        
        //! vector of adjacent surfaces for each line
        vector< vector <int> >              eleLinesSurfaces_;
        
        //! nodal reference coordinate
        vector< vector <double> >           eleRefCoordinates_;
        
        //! list of interface points to be stored in the Tetgen data structure
        vector< InterfacePoint >            pointList_;
    
        //! list of segments to be stored in the Tetgen data structure
        vector< vector<int> >               segmentList_;                 
    
        //! list of single points on a xfem surface to be stored in the Tetgen data structure
        vector< vector<int> >               surfacePointList_;                
    
        //! list of triangles to be stored in the Tetgen data structure
        vector< vector<int> >               triangleList_;
       
        //! list of intersected surface elements
        vector< DRT::Element* >             intersectedSurfaces_; 
        
        //! list of facet marker numbered in the order of the intersected surface elments
        vector<int>                         faceMarker_; 
        
        
        /*!
        \brief initializes the private members of the xfem element
    
        \param xfemElement   (in) : xfem element
        */
        void initializeXFEM(  
            DRT::Element* xfemElement
            );
            
            
        /*!
        \brief initializes the private members of the cutter element
    
        \param cutterElement   (in) : cutter element
        */
        void initializeCutter(  
            DRT::Element* cutterElement
            );
           		
        
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
        \brief transforms a point in reference coordinates to a point
               in current coordinates with respect to a given surfaceElement
               based on given nodal coordinates
    
        \param element          (in)        : surfaceElement
        \param xsi              (in/out)    : vector of reference coordinates/
                                              vector of current coordinates
        \param surfaceNodes     (in)        : nodal coordinates of the surface element
        */
        void referenceToCurrentCoordinates(   
            DRT::Element*                               surfaceElement, 
            Epetra_SerialDenseVector&                   xsi,
            const vector <Epetra_SerialDenseVector>&    surfaceNodes);
    
    
        void referenceToCurrentCoordinates(  
            Epetra_SerialDenseVector&                   xsi,
            const vector<Epetra_SerialDenseVector>&     plane
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
            const std::vector<double>&         node,
            const Epetra_SerialDenseMatrix&    XAABB);
        
                                                                                        
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
    
        \param element              (in)        : element (surface)
        \param node                 (in)        : node in current coordinates
        \param interfacePoints      (in/out)    : vector of interface points
        \param nodeId               (in)        : node id
        \return true if the node lies within the element, false otherwise
        */                                      
        bool collectInternalPoints( 
            DRT::Element*                   element,
            DRT::Element*                   surfaceElement,
            DRT::Node*                      node,
            std::vector<InterfacePoint>&    interfacePoints,
            const int                       elemId,
            const int                       nodeId);
        
                       							
        /*!
        \brief updates the system matrix at the corresponding reference coordinates for the 
               computation if a node in current coordinates lies within an element 
    
        \param dim               (in)       : dimension of the problem
        \param A                 (out)      : system matrix
        \param xsi               (in)       : vector of reference coordinates
        \param element           (in)       : element (volume)
        */                                                                                                      
        void updateAForNWE( 
            const int                   dim,
            Epetra_SerialDenseMatrix&   A,
            Epetra_SerialDenseVector&   xsi,
            DRT::Element*               element
            );
        
          
        /*!
        \brief updates the rhs at the corresponding reference coordinates for the 
               computation whether a node in current coordinates lies within an element 
    
        \param dim               (in)       : dimension of the problem
        \param b                 (out)      : right-hand-side       
        \param xsi               (in)       : vector of reference coordinates
        \param x                 (in)       : node in current coordinates
        \param element           (in)       : element (volume)
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
    
        \param element              (in)        : element (surface)
        \param node                 (in)        : node in current coordinates
        \param xsi                  (in)        : node in reference coordinates
        */	
		bool checkNodeWithinElement(	
            DRT::Element*                       element,
            const Epetra_SerialDenseVector&     x,
            Epetra_SerialDenseVector&           xsi
            );	       
        
        
        /*!
        \brief  checks if a node that lies within an xfem-element 
                lies on one of its surfaces, lines or nodes          .
    
        \param xsi                  (in)        : node in reference coordinates
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
    
        \param element                  (in)        : element
        \param cutterElement            (in)        : cutter element
        */       
        void computeCDT(    
            DRT::Element*                           element,
            DRT::Element*                           cutterElement,
            map< int, vector <Integrationcell> >&   integrationcellList
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
        \brief  computes the recovery of the curved interface
                for one xfem element after a Contrained Delaunay 
                Tetrahedralization
    
        \param element              (in)        :   xfem element
        \param out                  (in/out)    :   output data structure of TetGen
        */      
        void recoverCurvedInterface(
            DRT::Element*   element, 
            tetgenio&       out
            );
           
            
        /*!
        \brief  checks if Steiner points are inserted additionally. If the Steiner points
                are not lying in a surface element already they are lifted up such they
                are lying in a surface element.
    
        \param out                  (out)    :   output data structure of TetGen
        */     
        void liftAllCornerPoints(
            tetgenio&                                       out,
            DRT::Element*                                   xfemElement,
            vector< vector<Epetra_SerialDenseVector> >&     surfaceNodes
            );      
           
    
        /*!
        \brief  returns the other two point indices belonging to a triface 
                with a Steiner point
    
        \param out                  (out)    :   output data structure of TetGen
        \param trifaceIndex         (out)    :   index of the triface
        \param steinerPointIndex    (out)    :   index of the Steiner point
        */         
        vector<int> getPointIndices(
            tetgenio&   out, 
            int         trifaceIndex, 
            int         steinerPointIndex
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
        \param plane                    (in)    :   plane nodes
        \param xsi                      (out)   :   intersection point
        \return true if normal intersects with curved surface element or false otherwise
        */      
        bool computeRecoveryNormal(  
            vector <Epetra_SerialDenseVector>&          surfaceNodes,
            DRT::Element*                               surfaceElement,
            const vector<Epetra_SerialDenseVector>&     plane,
            Epetra_SerialDenseVector&                   xsi,
            const bool                                  current
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
        void updateAForRQINormal( 
            Epetra_SerialDenseMatrix&                   A,
            const Epetra_SerialDenseVector&             xsi,
            const vector<Epetra_SerialDenseVector>&     plane,
            const vector <Epetra_SerialDenseVector>&    surfaceNodes,
            DRT::Element*                               surfaceElement,
            const bool                                  current
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
        \param plane            (in)    :   plane nodes
        */                         
        void updateRHSForRQINormal(   
            Epetra_SerialDenseVector&                   b,
            Epetra_SerialDenseVector&                   xsi,
            const vector <Epetra_SerialDenseVector>&    surfaceNodes,
            DRT::Element*                               surfaceElement,
            const vector<Epetra_SerialDenseVector>&     plane,
            const bool                                  current
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
        \brief  gathers information of a tetrahedron adjacent to an 
                interface 
                
        \param out                  (in)    : output data structure of TetGen
        \param tetraCornerIndices   (out)   : index list of tetrahedron corners
        \param order                (out)   : order the tetrahedron corner points appear
        \param index                (in)    : tetrahedron index of the tetrahedron list
        */
        void getTetrahedronInformation(   
            const tetgenio&     out, 
            vector<int>&        tetraCornerIndices,
            vector<int>&        order,
            const int           tetIndex,
            const int           index
            );          
        
        
        /*!
        \brief  gathers information of a tetrahedron adjacent to an 
                interface 
                
        \param out                  (in)    : output data structure of TetGen
        \param xfemElement          (in)    : xfem element
        \param tetraCornerIndices   (in)    : index list of tetrahedron corners
        \param tetraCornerIndices   (out)   : list of tetrahedron corners nodes in the current 
                                              configuration of the xfem element
        */
        void getTetrahedronNodes(
            const tetgenio&                         out, 
            DRT::Element*                           xfemElement,
            vector<Epetra_SerialDenseVector>&       tetraCornerNodes,
            vector<int>&                            tetraCornerIndices
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
        \brief  computes the normal to the interface edge of the
                adjacent tetrahedon. In addition the normal lies within the plane
                of the tetrahedron. 
                
        \param index1               (in)    : index 1 of an interface edge
        \param index2               (in)    : index 2 of an interface edge
        \param out                  (in)    : out
        \param tetraCornerIndices   (in)    : indices of tetrahedron corners
        \param plane                (in)    : plane nodes
        */    
        void computeIntersectionNormal(    
            const int                           index1,
            const int                           index2,
            const int                           globalHigherOrderIndex,
            const tetgenio&                     out,
            const vector<int>&                  tetraCornerIndices,
            vector<Epetra_SerialDenseVector>&   plane
            );
               
        
        /*!
        \brief  computes the normal to the interface edge of the
                adjacent tetrahedon. In addition the normal lies within the plane
                of the tetrahedron. 
                
        \param index1               (in)    : index 1 of an interface edge
        \param index2               (in)    : index 2 of an interface edge
        \param out                  (in)    : out
        \param tetraCornerIndices   (in)    : indices of tetrahedron corners
        \param tetraCornerNodes     (in)    : list of tetrahedron corner nodes in 
                                              current configuration of the xfem element
        \param plane                (out)   : plane nodes
        \param current              (in)    : true if nodes are current configuration , false 
                                              if they are in reference configuration
        */    
        void computeIntersectionNormal(    
            DRT::Element*                           xfemElement,   
            const int                               index1,
            const int                               index2,
            const int                               globalHigherOrderIndex, 
            const tetgenio&                         out,
            const vector<int>&                      tetraCornerIndices,
            const vector<int>&                      order,
            const vector<Epetra_SerialDenseVector>& tetraCornerNodes,
            vector<Epetra_SerialDenseVector>&       plane,
            const bool                              current
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
            const int                                   index,
            Epetra_SerialDenseVector&                   xsi, 
            const vector <Epetra_SerialDenseVector>&    surfaceNodes,
            const vector <Epetra_SerialDenseVector>&    plane,
            int                                         lineIndex, 
            DRT::Element*                               surfaceElement, 
            DRT::Element*                               xfemElement, 
            tetgenio&                                   out,
            const bool                                  current, 
            const bool                                  normal 
            );
          
         
        /*!
        \brief  checks if the considered tetrahedron facet lies on 
                a xfem-element faxt
                
        \param index1               (in)    : index 1 of an interface edge
        \param index2               (in)    : index 2 of an interface edge
        \param corner index         (in)    : index lying on the tetrahedron facet 
                                              opposite to the interface edge
        \param tetIndex             (in)    : tetrahedron index
        \param out                  (in)    : output data structure of tetGen
        */           
        bool facetOnExternalBoundary(
            int         index1,
            int         index2,
            int         cornerIndex,
            int         tetIndex,
            tetgenio&   out
            );
     
     
        bool computeRecoveryPlane( 
            vector <Epetra_SerialDenseVector>&          surfaceNodes,
            DRT::Element*                               surfaceElement,
            int&                                        lineIndex,
            const vector<Epetra_SerialDenseVector>&     plane,
            Epetra_SerialDenseVector&                   xsi,
            const bool                                  current
            );
            
            
        void updateAForRQIPlane(   
            Epetra_SerialDenseMatrix&                   A,
            const Epetra_SerialDenseVector&             xsi,
            const vector<Epetra_SerialDenseVector>&     plane,
            const vector<Epetra_SerialDenseVector>&     curve,
            DRT::Element*                               surfaceElement,
            DRT::Element*                               lineElement,
            const bool                                  current
            ); 
            
      
        void updateRHSForRQIPlane( 
            Epetra_SerialDenseVector&                   b,
            Epetra_SerialDenseVector&                   xsi,    
            const vector <Epetra_SerialDenseVector>&    curve,
            DRT::Element*                               lineElement,
            const vector<Epetra_SerialDenseVector>&     plane,
            const bool                                  current
            ); 
            
            
        bool findCommonFaceEdge( 
            int                                 faceIndex1, 
            int                                 faceIndex2, 
            vector<int>&                        adjacentFacesList,
            Epetra_SerialDenseVector&           edgePoint,
            Epetra_SerialDenseVector&           oppositePoint,
            tetgenio&                           out
            );
            
      
        bool findCommonSurfaceEdge(  
            int                                             faceIndex1, 
            int                                             faceIndex2,
            int&                                            surfaceIndex,
            int&                                            lineIndex
            );
          
        int findAdjacentFaceMarker(
            int         edgeIndex1, 
            int         edgeIndex2, 
            int         faceMarkerIndex,
            tetgenio&   out
            );   
            
        int findIntersectingSurfaceEdge(
            DRT::Element*                       xfemElement,
            DRT::Element*                       surfaceElement,
            Epetra_SerialDenseVector            edgeNode1,
            Epetra_SerialDenseVector            edgeNode2
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
            
            
        /*!
        \brief  computes unique coordinates for a certain region in the 
                xfem element. It is copied into the data structure of Tetgen
                   
        \param indec            (in)        : index of the xfem element
        \param out              (in)        : tetgenio output data structure
        */        
        void printTetViewOutput(
            int             index,
            tetgenio&       out
            ); 
            
            
        void printTetViewOutputPLC(
            DRT::Element*   xfemElement,
            int             index,
            tetgenio&       in
            );   
                                               
};


};


#endif  // #ifndef INTERSECTION_H
#endif  // #ifdef CCADISCRET

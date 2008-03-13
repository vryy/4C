/*!
\file intersection_service.H

\brief collection of math tools for the interface computation of two
       curved meshes

<pre>
Maintainer: Ursula Mayer
</pre>
*/

#ifdef CCADISCRET
#ifndef INTERSECTION_SERVICE_H_
#define INTERSECTION_SERVICE_H_

#include <blitz/array.h>
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_element.H"
#include "Epetra_SerialDenseMatrix.h"
#include "Teuchos_RefCountPtr.hpp"
#include "xfem.H"


namespace XFEM
{

    //! tolerance 1e-14 used for residuals in all Newton-methods
    const double                 TOL14 = 1e-14;

    //! tolerance 1e-7 
    const double                 TOL7 = 1e-7;

    
    /*!
    \brief Returns the sum of two Epetra_SerialDenseVectors

    \param v1 (in) : arbitrary Epetra_SerialDenseVector
    \param v2 (in) : arbitrary Epetra_SerialDenseVector
    \return sum of two Epetra_SerialDenseVectors
    */
    template<int length>
    inline BlitzVec addTwoVectors(   
        const BlitzVec&   v1,
        const BlitzVec&   v2)
    {   
        return BlitzVec(v1 + v2);
    }  
    
    
    /*!
    \brief Returns the difference of Epetra_SerialDenseVectors v1 - v2
    \param v1 (in) : Epetra_SerialDenseVector      
    \param v2 (in) : Epetra_SerialDenseVector
    \return difference of Epetra_SerialDenseVectors v1 - v2
    */
    template<int length>
    inline BlitzVec subtractsTwoVectors( 
        const BlitzVec& v1,
        const BlitzVec& v2)
    {   
        return BlitzVec(v1 - v2);
    }
    
    /*!
    \brief computes the cross product of two BlitzVecs a x b
    \param a (in) : arbitrary BlitzVec     
    \param b (in) : arbitrary BlitzVec
    \return cross product of two BlitzVecs a x b
    */
    BlitzVec computeCrossProduct(
        const BlitzVec& a,
        const BlitzVec& b
        );
    
     
    /*!
    \brief normalizes a BlitzVec
    \param v (in) : BlitzVec
    \return normalized BlitzVec 
    */   
    void normalizeVector(
        BlitzVec&     v
        );
    
    
    /*!
    \brief  computes the Theorem of Pythagoras
            (a^2 + b^2)^(1/2)
            (modified from NUMERICAL RECIPES)
                   
    \param a    (in)        : a
    \param b    (in)        : b
    \return result of the Theorem of Pythagoras
    */    
    double pythagoras(
        const double  a, 
        const double  b
        );
    
    
//    /*!
//    \brief checks if element is a volume element
//    
//    \param element       (in)    : element
//    \return true if it is a volume element, false otherwise
//    */      
//    bool checkIfVolumeElement(
//            const DRT::Element::DiscretizationType distype
//        );
    
    
//    /*!
//    \brief checks if element is a surface element
//    
//    \param element       (in)    : element
//    \return true if it is a surface element, false otherwise
//    */     
//    bool checkIfSurfaceElement(
//            const DRT::Element::DiscretizationType distype
//        );
    
    
//    /*!
//    \brief checks if element is a line element
//
//    \param element       (in)    : element
//    \return true if it is a line element, false otherwise
//    */     
//    bool checkIfLineElement(
//            const DRT::Element::DiscretizationType distype
//        );
    
    
    /*!
    \brief transforms a point in element coordinates to a point
           in current coordinates with respect to a given element

    \param element          (in)        : element
    \param eleCoord         (in)    	: vector of element coordinates
    \return                               vector of current coordinates
    */
    BlitzVec elementToCurrentCoordinates(   
        const DRT::Element*  element, 
        const BlitzVec&      eleCoord
        );
    
    /*!
    \brief transforms a point in element coordinates to a point
           in current coordinates with respect to a given element

    \param element          (in)        : element
    \param eleCoord         (in)        : vector of element coordinates
    \return                               vector of current coordinates
    */
    void elementToCurrentCoordinatesInPlace(   
        const DRT::Element*  element, 
        BlitzVec&      eleCoord
        );

    
    /*!
    \brief transforms a point in current coordinates to a point
           in element coordinates with respect to a given element

    \param element          (in)        : element
    \param xsi              (in/out)    : vector of current coordinates/
                                          vector of element coordinates
    */
    BlitzVec currentToElementCoordinatesExact(
            const DRT::Element*         element,
            const BlitzVec&             x
            );
    
    
    /*!
    \brief transforms a point in current coordinates to a point
           in element coordinates with respect to a given element       
           The nonlinear system of equation is solved with help of the Newton-method.

    \param element              (in)        : element 
    \param x                    (in)        : node in current coordinates (x, y, z)
    \param xsi                  (out)       : node in element coordinates
    \return bool, indicates, whether positions was found within element coordinates
    */  
    bool currentToElementCoordinates(
        const DRT::Element*                 element,
        const BlitzVec&                     x,
        BlitzVec&                           xsi
        );
    
    
    /*!
    \brief compares two points (overloaded method)

    \param point1       (in)    : first point   (double*)
    \param point2       (in)    : second point  (double*)
    \param lenght       (in)    : length of arrays
    \return true if both points equal each other, false otherwise
    */
    template<int length>
    bool comparePoints(
        const double*     point1,
        const double*     point2) 
    {   
        bool equal = true;
                 
        for(int i = 0; i < length; i++)
            if(fabs(point1[i] - point2[i]) > TOL7)
            {
                equal = false;
                break;
            }
      
        return equal;
    }
    
    
    /*!
     \brief compares two points (overloaded method)
 
     \param point1       (in)    : first point   (vector<double>)
     \param point2       (in)    : second point  (double*)
     \return true if both points equal each other, false otherwise
     */
    template <int length>
    inline bool comparePoints(   
        const vector<double>&     point1,
        const double*             point2) 
    {   
        bool equal = true;
        dsassert(point1.size() == length, "not of proper length!");
        for(unsigned int i = 0; i < point1.size() ; ++i)
            if(fabs(point1[i] - point2[i]) > TOL7)
            {
                equal = false;
                break;
            }
        
        return equal;
    }
     
     
     /*!
     \brief compares two points (overloaded method)
     
     \param point1       (in)    : first point   (vector<double>)
     \param point2       (in)    : second point  (vector<double>)
     \return true if both points equal each other, false otherwise
     */
     template <int length>
     inline bool comparePoints(   
         const vector<double>& point1,
         const vector<double>& point2) 
     {   
         bool equal = true;
         dsassert(point1.size() == length, "not of proper length!");
         dsassert(point2.size() == length, "not of proper length!");
         for(unsigned int i = 0; i < length ; ++i)
             if(fabs(point1[i] - point2[i]) > TOL7)
             {
                 equal = false;
                 break;
             }
       
         return equal;
     }
     
     
     /*!
     \brief compares two points (overloaded method)
 
     \param point1       (in)    : first point   (BlitzVec)
     \param point2       (in)    : second point  (BlitzVec)
     \return true if both points equal each other, false otherwise
     */
     template <int length>
     inline bool comparePoints(    
         const BlitzVec&     point1,
         const BlitzVec&     point2) 
     {   
         bool equal = true;
         dsassert(point1.size() == length, "not of proper length!");
         dsassert(point2.size() == length, "not of proper length!");
         for(int i = 0; i < length ; i++)
             if(fabs(point1(i) - point2(i)) > TOL7)
             {
                 equal = false;
                 break;
             }
       
         return equal;
     }
     
     
    /*!
    \brief Checks if a node is within an 
           extended axis-aligned bounding box (XAABB)

    \param node (in)        : node to be checked
    \param XAABB (in)       : extended axis-aligned bouning box
    \return true if node is within the XAABB or false otherwise
    */
    bool isPositionWithinXAABB(
            const BlitzVec&                    pos,
            const Epetra_SerialDenseMatrix&    XAABB
            );
    
    
    /*!
    \brief Checks if a straightline with 2 nodes is within an 
           extended axis-aligned bounding box (XAABB)

    \param node1 (in)       : node to be checked
    \param node2 (in)       : node to be checked
    \param XAABB (in)       : extended axis-aligned bouning box
    \return true if node is within the XAABB or false otherwise
    */
    bool isLineWithinXAABB(
            const BlitzVec&                    pos1,
            const BlitzVec&                    pos2,
            const Epetra_SerialDenseMatrix&    XAABB
            );
    
    
    /*!
    \brief checks if a position in current coordinates lies within a certain element           
           The nonlinear system of equation is solved with help of the Newton-method.
    
    \param element              (in)        : element 
    \param x                    (in)        : node in current coordinates
    */  
    bool checkPositionWithinElement(  
        const DRT::Element*                 element,
        const BlitzVec&                     x
        );   
    
    
    /*!
    \brief checks if a position in current coordinates lies within a certain element           
           The nonlinear system of equation is solved with help of the Newton-method.

    \param discretization       (in)        : mesh 
    \param x                    (in)        : node in current coordinates
    */  
    bool PositionWithinDiscretization(  
        const RCP<DRT::Discretization>      dis,
        const BlitzVec&                     x);
    
    
    /*!
    \brief checks if a position in current coordinates lies within a certain element           
           The nonlinear system of equation is solved with help of the Newton-method.
    */  
    bool PositionWithinCondition(
            const BlitzVec&                   x,                     ///< node in current coordinates
            const int                         xfem_condition_label, 
            const RCP<DRT::Discretization>    cutterdis              ///< mesh
        );
    
    
    /*!
    \brief  finds the nearest point for a given point on a surface element
            the element coordinates of this point are stored in eleCoord, the normal
            point from the given point to the nearest point on the surface,
            distance returns the absolut length of the normal. 
            If the computed point lies outside the element , the distance which is retuned equals -1
            and the normal is set to zero.

    \param surfaceElement       (in)       : surface element
    \param physCoord            (out)      : physical coordinates      
    \param eleCoord             (in)       : element coordinates
    \param normal               (in)       : normal point from a given point to the surface point
    \param distance             (in)       : length of the normal
    */         
    bool searchForNearestPointOnSurface(
            const DRT::Element*                     surfaceElement,
            const BlitzVec&                         physCoord,
            BlitzVec&                               eleCoord,
            BlitzVec&                               normal,
            double&                                 distance
            );
    
    
    /*!
       \brief checks if a position in current coordinates lies within a certain surfaceElement           
          The nonlinear system of equation is solved with help of the Newton-method.
    
       \param surfaceElement   (in)        : surfaceElement 
       \param physCoord           	(in)        : node in physical coordinates (x, y, z)
       \return									  node in element coordinates (r, s)
    */  
    BlitzVec CurrentToSurfaceElementCoordinates(
        const DRT::Element*      surfaceElement,
        const BlitzVec&          physCoord
        );
    
    
    /*!
    \brief checks if a position in element coordinates lies within a certain Element parameter space

    \return  true if point lies on the surface element, false otherwise
    */  
    bool checkPositionWithinElementParameterSpace(
            const BlitzVec&                        eleCoord, ///< position in element coordinates
            const DRT::Element::DiscretizationType distype   ///< shape of the surface element
            );
    
    
    /*!
    \brief 	updates the 3x2 Jacobian matrix for the computation of the
    		surface element coordinates for a point in physical coordinates

    \param xsi          		(in)          : element coordinates      
    \param surface element		(in)        : surface element
    
    \return Jacobian matrix
    */         
    BlitzMat updateJacobianForMap3To2(   
        const BlitzVec&	xsi,
        const DRT::Element*            	surfaceElement
        );      
    
   
    /*!
    \brief 	updates the nonlinear equations for the computation of the
    		surface element coordinates for a surface point in physical coordinates

    \param F   				(out)       : nonlinear equations
    \param xsi          	(in)      	: element coordinates    
    \param x         		(in)      	: physical coordinates      
    \param surface element	(in)       	: surface element
    */    
    BlitzVec updateFForMap3To2(   
        const BlitzVec&                 xsi,
        const BlitzVec&	                x,
        const DRT::Element*            	surfaceElement
        );
    
    
    /*!
    \brief 	updates the nonlinear equations for the computation of the
    		surface element coordinates for a surface point in physical coordinates

    \param A  				(out)       : system matrix
    \param Jacobi          	(in)      	: Jacobian matrix
    \param F         		(in)      	: nonlinear equations  
    \param xsi         		(in)      	: element coordinates      
    \param surface element	(in)       	: surface element
    */    
    void updateAForMap3To2(   
            BlitzMat&                       A,
            const BlitzMat&                 Jacobi,
            const BlitzVec&                 F,
            const BlitzVec&                 xsi,
            const DRT::Element*             surfaceElement
            );
    
    
    /*!
    \brief Computes a rough overestimating extended  
           axis-aligned bounding box for an element (XAABB)

    \param  element          (in)      : element
    \return extended axis-aligned bounding box  (XAABB) for an element
    */
    Epetra_SerialDenseMatrix computeFastXAABB( 
            const DRT::Element*   element
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
    \brief checks if a position is within an extended axis aligned bounding boxe

    \param pos     (in)         : position
    \param XAABB   (in)         : XAABB of the xfem element
    \return true if the position is within the XAABB, false otherwise
    */
    bool PositionInXAABB(
            const BlitzVec&                    pos,
            const Epetra_SerialDenseMatrix&    XAABB
            );
    
}


#endif  //INTERSECTION_MATH_H_
#endif  // #ifdef CCADISCRET



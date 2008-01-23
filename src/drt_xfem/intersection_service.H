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
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RefCountPtr.hpp"


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
    Epetra_SerialDenseVector addTwoVectors(     
        const Epetra_SerialDenseVector&   v1,
        const Epetra_SerialDenseVector&   v2
        );      
    
    
    /*!
    \brief Returns the sum of two vector<double>

    \param v1 (in) :    arbitrary vector<double>
    \param v2 (in) :    arbitrary vector<double>
    \return sum of two vector<double>
    */
    vector<double> addTwoVectors(   
        const vector<double>&   v1,
        const vector<double>&   v2
        );
    
                                          
    /*!
    \brief Returns the difference of Epetra_SerialDenseVectors v1 - v2
    \param v1 (in) : Epetra_SerialDenseVector      
    \param v2 (in) : Epetra_SerialDenseVector
    \return difference of Epetra_SerialDenseVectors v1 - v2
    */
    Epetra_SerialDenseVector subtractsTwoVectors(   
        const Epetra_SerialDenseVector&   v1,
        const Epetra_SerialDenseVector&   v2
        );       
    
                                        
    /*!
    \brief Returns the difference of two vector<double> v1 - v2
    \param v1 (in) : arbitrary vector<double>       
    \param v2 (in) : arbitrary vector<double>
    \return difference of two vector<double> v1 - v2
    */
    vector<double> subtractsTwoVectors( 
        const vector <double>&   v1,
        const vector <double>&   v2
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
    \brief normalizes a Epetra_SerialDenseVector
    \param v (in) : Epetra_SerialDenseVector
    \return normalized Epetra_SerialDenseVector 
    */   
    void normalizeVector(   
        Epetra_SerialDenseVector&     v
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

    
    
    /*!
    \brief checks if element is a volume element

    \param element       (in)    : element
    \return true if it is a volume element, false otherwise
    */      
    bool checkIfVolumeElement(
            const DRT::Element::DiscretizationType distype
        );
    
    
    /*!
    \brief checks if element is a surface element

    \param element       (in)    : element
    \return true if it is a surface element, false otherwise
    */     
    bool checkIfSurfaceElement(
            const DRT::Element::DiscretizationType distype
        );
    
    
    /*!
    \brief checks if element is a line element

    \param element       (in)    : element
    \return true if it is a line element, false otherwise
    */     
    bool checkIfLineElement(
            const DRT::Element::DiscretizationType distype
        );
    
    
    
    /*!
    \brief Checks if a node is within an 
           extended axis-aligned bounding box (XAABB)

    \param node (in)        : node to be checked
    \param XAABB (in)       : extended axis-aligned bouning box
    \return true if node is within the XAABB or false otherwise
    */
    bool isPositionWithinXAABB(  
        const Epetra_SerialDenseVector&    pos,
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
        const Epetra_SerialDenseVector&    pos1,
        const Epetra_SerialDenseVector&    pos2,
        const Epetra_SerialDenseMatrix&    XAABB
        );
        
    
    /*!
    \brief checks if a position in current coordinates lies within a certain element           
           The nonlinear system of equation is solved with help of the Newton-method.

    \param element              (in)        : element 
    \param x                    (in)        : node in current coordinates
    \param xsi                  (inout)     : node in element coordinates
    */  
    bool checkPositionWithinElement(    
        DRT::Element*                       element,
        const Epetra_SerialDenseVector&     x,
        Epetra_SerialDenseVector&           xsi
        );
    
    /*!
    \brief checks if a position in current coordinates lies within a certain element           
           The nonlinear system of equation is solved with help of the Newton-method.

    \param element              (in)        : element 
    \param x                    (in)        : node in current coordinates
    */  
    bool checkPositionWithinElement(    
        DRT::Element*                       element,
        const Epetra_SerialDenseVector&     x
        );   
    
    /*!
    \brief checks if a position in current coordinates lies within a certain element           
           The nonlinear system of equation is solved with help of the Newton-method.

    \param discretization       (in)        : mesh 
    \param x                    (in)        : node in current coordinates
    */  
    bool PositionWithinDiscretization(  
        RCP<DRT::Discretization>            dis,
        const Epetra_SerialDenseVector&     x);
    
    /*!
    \brief checks if a position in current coordinates lies within a certain element           
           The nonlinear system of equation is solved with help of the Newton-method.
    */  
    bool PositionWithinCondition(
            const blitz::Array<double,1>&     x,                     ///< node in current coordinates
            const int                         xfem_condition_label, 
            RCP<DRT::Discretization>          cutterdis              ///< mesh
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
    \brief Computes a rough overestimating extended  
           axis-aligned bounding box for an element (XAABB)

    \param  element          (in)      : element (DRT::Element*)
    \return extended axis-aligned bounding box  (XAABB) for an element
    */
    Epetra_SerialDenseMatrix computeFastXAABB( 
        DRT::Element*   element
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
        const Epetra_SerialDenseVector&    pos,
        const Epetra_SerialDenseMatrix&    XAABB
        );
    
}


#endif  //INTERSECTION_MATH_H_
#endif  // #ifdef CCADISCRET



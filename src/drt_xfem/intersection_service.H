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


#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_element.H"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RefCountPtr.hpp"


using namespace std;
using namespace DRT;
using namespace DRT::Utils;
using namespace Teuchos;


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
       
        
};


#endif  //INTERSECTION_MATH_H_
#endif  // #ifdef CCADISCRET



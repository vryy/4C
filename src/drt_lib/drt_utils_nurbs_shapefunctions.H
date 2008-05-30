/*!----------------------------------------------------------------------
\file drt_utils_nurbs_shapefunctions.H

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef DRT_UTILS_NURBS_SHAPEFUNCTION
#define DRT_UTILS_NURBS_SHAPEFUNCTION

#include<vector>
#include<iostream>

#include <blitz/array.h>
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_utils_bspline.H"


namespace DRT
{
namespace NURBS
{
namespace UTILS
{
  
/*!
\brief Evaluate basis functions of nurbs basis functions.

\verbatim

       reference element          element in knotspan
      
             v^                        eta^
              |                           | 
         -1   |  +1                 *   *   *   * 
       +1 +---+---+ +1             *+***+---+***+*
          |XXX|XXX|         Psi     *   |XXX|   * 
        --+---+---+-->    ------>   *   |XXX|   * --->
          |XXX|XXX|   u             *   |XXX|   *   xi
       -1 +---+---+ -1		   *+***+---+***+*
         -1   |  +1		    *   *   *   * 
                                        *   *
                \       (knot[0])(degree)   (knot[0])(degree+1)
                 \
                  \                       |
                   \                      |
                    \                     |
                     \                    |
                      \                   |
                       \                  |
                        \                 |
                         \                |      +----
                          \               |       \
               Phi o Psi   \              | Phi =  +  N_{i,j}(xi,eta)*f_{i,j}
              ^             \             |       /
              N_{i,j}(u,v)   \            |      +----
                              \           |
                               \          |
                                \         |
                                 \        |
                                  \       |
                                          v
                                  *              
      			         *+-------+*    
      				   \XXXXXX|     
                               	    \XXXXX|     
                                     \XXXX|         
      				    **+---+**   
      				      *   *     
      
                            'real' geometry element

\endverbatim

PART   I: EVALUATION OF  BEZIER SHAPE FUNCTIONS
              
\verbatim

           B_{i,j}(xi,eta) = N_i(xi)*M_j(eta)

\endverbatim

PART  II: PROJECTING TO NURBS SHAPE FUNCTIONS

\verbatim
      
                      w_{i,j} * B_{i,j}(xi,eta)
                
       N_{i,j} = ----------------------------------
                  +----+				
                   \				
                    +    w_{k,l} * B_{k,l}(xi,eta)
      		   /				
      		  +----+				
      		    k,l                           

\endverbatim


\param nurbs_shape_funct  (o) 
\param uv                 (i) 
\param knots              (i) 
\param weights            (i) 
\param distype            (i) 

\return void

*/

void nurbs_get_2D_funct(
    blitz::Array<double,1>                 & nurbs_shape_funct,
    const blitz::Array<double,1>           & uv               ,
    const vector<blitz::Array<double, 1> > & knots            ,
    const blitz::Array<double, 1>          & weights          ,
    const DRT::Element::DiscretizationType & distype          );


/*!
\brief Evaluate basis functions and derivatives (with 
respect to the coordinates of the reference element) 
of nurbs basis functions.

\verbatim

       reference element          element in knotspan
      
             v^                        eta^
              |                           | 
         -1   |  +1                 *   *   *   * 
       +1 +---+---+ +1             *+***+---+***+*
          |XXX|XXX|         Psi     *   |XXX|   * 
        --+---+---+-->    ------>   *   |XXX|   * --->
          |XXX|XXX|   u             *   |XXX|   *   xi
       -1 +---+---+ -1		   *+***+---+***+*
         -1   |  +1		    *   *   *   * 
                                        *   *
                \       (knot[0])(degree)   (knot[0])(degree+1)
                 \
                  \                       |
                   \                      |
                    \                     |
                     \                    |
                      \                   |
                       \                  |
                        \                 |
                         \                |      +----
                          \               |       \
               Phi o Psi   \              | Phi =  +  N_{i,j}(xi,eta)*f_{i,j}
              ^             \             |       /
              N_{i,j}(u,v)   \            |      +----
                              \           |
                               \          |
                                \         |
                                 \        |
                                  \       |
                                          v
                                  *              
      			         *+-------+*    
      				   \XXXXXX|     
                               	    \XXXXX|     
                                     \XXXX|         
      				    **+---+**   
      				      *   *     
      
                            'real' geometry element

\endverbatim

PART   I: EVALUATION OF  BEZIER SHAPE FUNCTIONS
              
\verbatim


           B_{i,j}(xi,eta) = N_i(xi) *M_j (eta)


          dB_{i,j}   
          --------(xi,eta) = N_i'(xi)*M_j (eta)
             dxi 

    
          dB_{i,j}   
          --------(xi,eta) = N_i (xi)*M_j'(eta)
            deta 

\endverbatim

PART  II: PROJECTING TO NURBS SHAPE FUNCTIONS

\verbatim
      
                      w_{i,j} * B_{i,j}(xi,eta)
                
       N_{i,j} = ----------------------------------
                  +----+				
                   \				
                    +    w_{k,l} * B_{k,l}(xi,eta)
      		   /				
      		  +----+				
      		    k,l                           

\endverbatim

Derivatives are evaluated using the chain rule


PART III: TRANSFORMING DERIVATIVES FROM xi/eta TO u/v

\verbatim

              +-  dxi   dxi  -+   +-   du       -+
              |   ---   ---   |   |   ---   0    |
              |    du    dv   |   |   2.0        |
       Psi' = |               | = |              |
              |  deta  deta   |   |         dv   |
              |  ----  ----   |   |    0   ---   |
              +-  du    dv   -+   +-       2.0  -+
    
\endverbatim

\param nurbs_shape_funct  (o) 
\param nurbs_shape_deriv  (o) 
\param uv                 (i) 
\param knots              (i) 
\param weights            (i) 
\param distype            (i) 

\return void

*/

void nurbs_get_2D_funct_deriv(
    blitz::Array<double,1>                 & nurbs_shape_funct,
    blitz::Array<double,2>                 & nurbs_shape_deriv,
    const blitz::Array<double,1>           & uv               ,
    const vector<blitz::Array<double, 1> > & knots            ,
    const blitz::Array<double, 1>          & weights          ,
    const DRT::Element::DiscretizationType & distype          );


/*!

\brief Evaluate basis functions, first and second derivatives 
       (with respect to the coordinates of the reference 
       element) of nurbs basis functions.

\verbatim

       reference element          element in knotspan
      
             v^                        eta^
              |                           | 
         -1   |  +1                 *   *   *   * 
       +1 +---+---+ +1             *+***+---+***+*
          |XXX|XXX|         Psi     *   |XXX|   * 
        --+---+---+-->    ------>   *   |XXX|   * --->
          |XXX|XXX|   u             *   |XXX|   *   xi
       -1 +---+---+ -1		   *+***+---+***+*
         -1   |  +1		    *   *   *   * 
                                        *   *
                \       (knot[0])(degree)   (knot[0])(degree+1)
                 \
                  \                       |
                   \                      |
                    \                     |
                     \                    |
                      \                   |
                       \                  |
                        \                 |
                         \                |      +----
                          \               |       \
               Phi o Psi   \              | Phi =  +  N_{i,j}(xi,eta)*f_{i,j}
              ^             \             |       /
              N_{i,j}(u,v)   \            |      +----
                              \           |
                               \          |
                                \         |
                                 \        |
                                  \       |
                                          v
                                  *              
      			         *+-------+*    
      				   \XXXXXX|     
                               	    \XXXXX|     
                                     \XXXX|         
      				    **+---+**   
      				      *   *     
      
                            'real' geometry element

\endverbatim

PART   I: EVALUATION OF  BEZIER SHAPE FUNCTIONS
              
\verbatim


           B_{i,j}(xi,eta) = N_i(xi) *M_j (eta)


          dB_{i,j}   
          --------(xi,eta) = N_i'(xi)*M_j (eta)
             dxi 

    
          dB_{i,j}   
          --------(xi,eta) = N_i (xi)*M_j'(eta)
            deta 

\endverbatim

PART  II: PROJECTING TO NURBS SHAPE FUNCTIONS

\verbatim
      
                      w_{i,j} * B_{i,j}(xi,eta)
                
       N_{i,j} = ----------------------------------
                  +----+				
                   \				
                    +    w_{k,l} * B_{k,l}(xi,eta)
      		   /				
      		  +----+				
      		    k,l                           

\endverbatim

First and second derivatives are evaluated using 
the chain rule. The second derivatives are already 
quite painful to compute...


PART III: TRANSFORMING DERIVATIVES FROM xi/eta TO u/v

\verbatim

              +-  dxi   dxi  -+   +-   du       -+
              |   ---   ---   |   |   ---   0    |
              |    du    dv   |   |   2.0        |
       Psi' = |               | = |              |
              |  deta  deta   |   |         dv   |
              |  ----  ----   |   |    0   ---   |
              +-  du    dv   -+   +-       2.0  -+
    
\endverbatim


\param nurbs_shape_funct  (o) 
\param nurbs_shape_deriv  (o) 
\param nurbs_shape_deriv2 (o) 
\param uv                 (i) 
\param knots              (i) 
\param weights            (i) 
\param distype            (i) 

\return void

*/

void nurbs_get_2D_funct_deriv_deriv2(
    blitz::Array<double,1>                 & nurbs_shape_funct ,
    blitz::Array<double,2>                 & nurbs_shape_deriv ,
    blitz::Array<double,2>                 & nurbs_shape_deriv2,
    const blitz::Array<double,1>           & uv                ,
    const vector<blitz::Array<double, 1> > & knots             ,
    const blitz::Array<double, 1>          & weights           ,
    const DRT::Element::DiscretizationType & distype           );


}// namespace UTILS

}// namespace NURBS

}// namespace DRT

#endif // DRT_UTILS_NURBS_SHAPEFUNCTION

#endif // CCADISCRET

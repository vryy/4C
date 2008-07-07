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
#include "../drt_fem_general/drt_utils_bspline.H"


namespace DRT
{
namespace NURBS
{
namespace UTILS
{

  
  /*
   
  \brief Evaluate basis functions of nurbs basis functions (line element version).
  
  \verbatim
  
          reference element          element in knotspan
                              Psi
          --+XXX+XXX+-->    ------>   ****+XXX+**** --->
           -1      +1   u                             xi
                  \            knot(degree)   knot(degree+1)D
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
                 Phi o Psi   \              | Phi =  +  N_{i}(xi)*f_{i}
                  ^           \             |       /
                  N_{i}(u,v)   \            |      +----
                                \           |
                                 \          |
                                  \         |
                                   \        |
                                    \       |
                                            v
                                         *
                                          *
                                           +
                                           X
                                           X
                                           X
                                           +
                                            *
                                             *
                                 'real' geometry element

   \endverbatim

   PART   I: EVALUATION OF  BEZIER SHAPE FUNCTIONS

   \verbatim

   B_{i}(xi) = N_i(xi)

   \endverbatim

   PART  II: PROJECTING TO NURBS SHAPE FUNCTIONS
   
   \verbatim

                        w_{i} * B_{i}(xi)
            N_{i} = --------------------------
                     +----+
                      \
                       +    w_{k} * B_{k}(xi)
                      /
                     +----+
                       k
   \endverbatim
   
   
    \param nurbs_shape_funct  (o)
    \param u                  (i)
    \param knots              (i)
    \param weights            (i)
    \param distype            (i)
   
    \return void
    
  */

  void nurbs_get_1D_funct(
    blitz::Array<double,1>                 & nurbs_shape_funct  ,
    const double                           & u                  ,
    const blitz::Array<double, 1>          & knots              ,
    const blitz::Array<double, 1>          & weights            ,
    const DRT::Element::DiscretizationType & distype            );


  /*!
    \brief Evaluate basis functions and derivatives (with
    respect to the coordinates of the reference element)
    of nurbs basis functions --- version for 1d line
    elements.
   
   \verbatim
   
           reference element          element in knotspan
                                 Psi
          --+XXX+XXX+-->    ------>   ****+XXX+**** --->
           -1      +1   u                             xi
                  \            knot(degree)   knot(degree+1)D
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
                 Phi o Psi   \              | Phi =  +  N_{i}(xi)*f_{i}
                  ^           \             |       /
                  N_{i}(u,v)   \            |      +----
                                \           |
                                 \          |
                                  \         |
                                   \        |
                                    \       |
                                            v
                                         *
                                          *
                                           +
                                           X
                                           X
                                           X
                                           +
                                            *
                                             *
                                 'real' geometry element

   \endverbatim

   PART   I: EVALUATION OF  BEZIER SHAPE FUNCTIONS
   
  \verbatim

              B_{i}(xi) = N_i(xi)


              dB_{i}
              ------(xi) = N_i'(xi)
               dxi

  \endverbatim

  PART  II: PROJECTING TO NURBS SHAPE FUNCTIONS

  \verbatim

                      w_{i} * B_{i}(xi)
          N_{i} = --------------------------
                  +----+
                   \
                    +    w_{k} * B_{k}(xi)
                   /
                  +----+
                    k
                    
  \endverbatim
  
  Derivatives are evaluated using the chain rule
  
  PART III: TRANSFORMING DERIVATIVES FROM xi TO u
  
  \verbatim
   
            +-   -+   +-   -+
            | dxi |   |  du |
     Psi' = | --- | = | --- |
            |  du |   | 2.0 |
            +-   -+   +-   -+
  
  \endverbatim

  \param nurbs_shape_funct  (o)
  \param nurbs_shape_deriv  (o)
  \param uv                 (i)
  \param knots              (i)
  \param weights            (i)
  \param distype            (i)
  
  \return void
  
  */

  void nurbs_get_1D_funct_deriv(
    blitz::Array<double,1>                 & nurbs_shape_funct,
    blitz::Array<double,1>                 & nurbs_shape_deriv,
    const double                           & u                ,
    const blitz::Array<double, 1>          & knots            ,
    const blitz::Array<double, 1>          & weights          ,
    const DRT::Element::DiscretizationType & distype          );

  
/*!
\brief Evaluate 2d basis functions of nurbs basis functions.

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
\brief Evaluate 2d basis functions and derivatives (with 
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

\brief Evaluate 2d basis functions, first and second derivatives 
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


/*!
\brief Evaluate 2d basis functions of nurbs basis functions.

\verbatim

       reference element          element in knotspan
                             
                                              nu ^		      
				                 |		      
                w		                 |   		           
                |		                 |  		      
                |		                 |  		      
              +-|--+----+                        +---------+            
             /  | /    /|                       /         /|            
            /   |/    / |                      /         / |            
           +----+----+  +        Psi          /         /  |            
          /    / 1  /| /|      ------->      /         /   |            
         /    /    / |/ |                   /         /    |                 
        +----+----+ 1+-------- v           +---------+     +-------->     
        |    |    | /| /                   |         |    /          xi   
        |    |1   |/ |/                    |         |   /                
        +----+----+  +                     |         |  /                 
        |   /|    | /                      |         | /                  
        |  / |    |/                       |         |/                   
        +-/--+----+                        +---------+                    
         /                                /|         |                      
      u /     \                          / |         |                      
               \                        /  | 	     |		        
		\		       /   |         | 			      
		 \		  eta v    |         |                     
                  \                        |         |
                   \                       |         |
                    \                      |         |
                     \                     |         |
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
                     Phi o Psi   \              | Phi =  +  N_{i,j,k}(xi,eta,nu)*f_{i,j,k}
                ^                 \             |       /
                N_{i,j,k}(u,v,w)   \            |      +----
                                    \           |
                                     \          |
                                      \         |
                                       \        |
                                        \       |
                        
                                             _____    
                                          +--     --+                               
                                         /         / \
                                        /         /   |
                                       /         /    |
                                      /         /     |
                                     / _____   /     /
                                    +--     --+     +
                                     \       /     / 
                                      |     |     /  
                                      |     |    /   
                                      |     |   /    
                                     /       \ /     
                                    +---------+      
                         
                             'real' geometry element

\endverbatim

PART   I: EVALUATION OF  BEZIER SHAPE FUNCTIONS
              
\verbatim

           B_{i,j,k}(xi,eta,nu) = N_i(xi) *M_j (eta)*L_k (nu)

\endverbatim

PART  II: PROJECTING TO NURBS SHAPE FUNCTIONS

\verbatim
      
                        w_{i,j,k} * B_{i,j,k}(xi,eta,nu)
                
       N_{i,j,k} = ---------------------------------------
                   +----+				
                    \				
                     +    w_{s,t,r} * B_{s,t,r}(xi,eta,nu)
      		    /				
      		   +----+				
      		   s,t,r                           

\endverbatim


\param nurbs_shape_funct  (o) 
\param uv                 (i) 
\param knots              (i) 
\param weights            (i) 
\param distype            (i) 

\return void

*/

void nurbs_get_3D_funct(
    blitz::Array<double,1>                 & nurbs_shape_funct,
    const blitz::Array<double,1>           & uv               ,
    const vector<blitz::Array<double, 1> > & knots            ,
    const blitz::Array<double, 1>          & weights          ,
    const DRT::Element::DiscretizationType & distype          );


/*!
\brief Evaluate 3d basis functions and derivatives (with 
respect to the coordinates of the reference element) 
of nurbs basis functions.

\verbatim

       reference element          element in knotspan
                                  
                                              nu ^		      
				                 |		      
                w		                 |   		           
                |		                 |  		      
                |		                 |  		      
              +-|--+----+                        +---------+            
             /  | /    /|                       /         /|            
            /   |/    / |                      /         / |            
           +----+----+  +        Psi          /         /  |            
          /    / 1  /| /|      ------->      /         /   |            
         /    /    / |/ |                   /         /    |                 
        +----+----+ 1+-------- v           +---------+     +-------->     
        |    |    | /| /                   |         |    /          xi   
        |    |1   |/ |/                    |         |   /                
        +----+----+  +                     |         |  /                 
        |   /|    | /                      |         | /                  
        |  / |    |/                       |         |/                   
        +-/--+----+                        +---------+                    
         /                                /|         |                      
      u /     \                          / |         |                      
               \                        /  | 	     |		        
		\		       /   |         | 			      
		 \		  eta v    |         |                     
                  \                        |         |
                   \                       |         |
                    \                      |         |
                     \                     |         |
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
                     Phi o Psi   \              | Phi =  +  N_{i,j,k}(xi,eta,nu)*f_{i,j,k}
                ^                 \             |       /
                N_{i,j,k}(u,v,w)   \            |      +----
                                    \           |
                                     \          |
                                      \         |
                                       \        |
                                        \       |
                        
                                             _____    
                                          +--     --+                               
                                         /         / \
                                        /         /   |
                                       /         /    |
                                      /         /     |
                                     / _____   /     /
                                    +--     --+     +
                                     \       /     / 
                                      |     |     /  
                                      |     |    /   
                                      |     |   /    
                                     /       \ /     
                                    +---------+      
                         
                             'real' geometry element

\endverbatim

PART   I: EVALUATION OF  BEZIER SHAPE FUNCTIONS
              
\verbatim


           B_{i,j,k}(xi,eta,nu) = N_i(xi) *M_j (eta)*L_k (nu)


          dB_{i,j,k}   
          ----------(xi,eta,nu) = N_i'(xi)*M_j (eta)*L_k (nu)
              dxi 

    
          dB_{i,j,k}   
          ----------(xi,eta,nu) = N_i (xi)*M_j'(eta)*L_k (nu)
             deta 

    
          dB_{i,j,k}   
          ----------(xi,eta,nu) = N_i (xi)*M_j (eta)*L_k'(nu)
             dnu 

\endverbatim

PART  II: PROJECTING TO NURBS SHAPE FUNCTIONS

\verbatim
      
                        w_{i,j,k} * B_{i,j,k}(xi,eta,nu)
                
       N_{i,j,k} = ---------------------------------------
                   +----+				
                    \				
                     +    w_{s,t,r} * B_{s,t,r}(xi,eta,nu)
      		    /				
      		   +----+				
      		   s,t,r                           

\endverbatim

Derivatives are evaluated using the chain rule


PART III: TRANSFORMING DERIVATIVES FROM xi/eta/nu TO u/v/w

\verbatim

              +   dxi   dxi   dxi -+   +-  du           -+
              |   ---   ---   ---  |   |  ---   0    0   |
              |    du    dv    dw  |   |  2.0            |
              |                    |   |                 |
              |  deta  deta  deta  |   |        dv       |
       Psi' = |  ----  ----  ----  | = |   0   ---   0   |
              |   du    dv    dw   |   |       2.0       |
              |                    |   |                 |
              |   dnu   dnu   dnu  |   |             dw  |
              |  ----  ----  ----  |   |   0    0   ---  |
              +-  du    dv    dw  -+   +-           2.0 -+

\endverbatim

\param nurbs_shape_funct  (o) 
\param nurbs_shape_deriv  (o) 
\param uv                 (i) 
\param knots              (i) 
\param weights            (i) 
\param distype            (i) 

\return void

*/

void nurbs_get_3D_funct_deriv(
    blitz::Array<double,1>                 & nurbs_shape_funct,
    blitz::Array<double,2>                 & nurbs_shape_deriv,
    const blitz::Array<double,1>           & uv               ,
    const vector<blitz::Array<double, 1> > & knots            ,
    const blitz::Array<double, 1>          & weights          ,
    const DRT::Element::DiscretizationType & distype          );


/*!

\brief Evaluate 3d basis functions, first and second derivatives 
       (with respect to the coordinates of the reference 
       element) of nurbs basis functions.

\verbatim

       reference element          element in knotspan
                                  
                                              nu ^		      
				                 |		      
                w		                 |   		           
                |		                 |  		      
                |		                 |  		      
              +-|--+----+                        +---------+            
             /  | /    /|                       /         /|            
            /   |/    / |                      /         / |            
           +----+----+  +        Psi          /         /  |            
          /    / 1  /| /|      ------->      /         /   |            
         /    /    / |/ |                   /         /    |                 
        +----+----+ 1+-------- v           +---------+     +-------->     
        |    |    | /| /                   |         |    /          xi   
        |    |1   |/ |/                    |         |   /                
        +----+----+  +                     |         |  /                 
        |   /|    | /                      |         | /                  
        |  / |    |/                       |         |/                   
        +-/--+----+                        +---------+                    
         /                                /|         |                      
      u /     \                          / |         |                      
               \                        /  | 	     |		        
		\		       /   |         | 			      
		 \		  eta v    |         |                     
                  \                        |         |
                   \                       |         |
                    \                      |         |
                     \                     |         |
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
                     Phi o Psi   \              | Phi =  +  N_{i,j,k}(xi,eta,nu)*f_{i,j,k}
                ^                 \             |       /
                N_{i,j,k}(u,v,w)   \            |      +----
                                    \           |
                                     \          |
                                      \         |
                                       \        |
                                        \       |
                        
                                             _____    
                                          +--     --+                               
                                         /         / \
                                        /         /   |
                                       /         /    |
                                      /         /     |
                                     / _____   /     /
                                    +--     --+     +
                                     \       /     / 
                                      |     |     /  
                                      |     |    /   
                                      |     |   /    
                                     /       \ /     
                                    +---------+                               
                             'real' geometry element

\endverbatim

PART   I: EVALUATION OF  BEZIER SHAPE FUNCTIONS
              
\verbatim



           B_{i,j,k}(xi,eta,nu) = N_i(xi) *M_j (eta)*L_k (nu)


          dB_{i,j,k}   
          ----------(xi,eta,nu) = N_i'(xi)*M_j (eta)*L_k (nu)
              dxi 

    
          dB_{i,j,k}   
          ----------(xi,eta,nu) = N_i (xi)*M_j'(eta)*L_k (nu)
             deta 

    
          dB_{i,j,k}   
          ----------(xi,eta,nu) = N_i (xi)*M_j (eta)*L_k'(nu)
             dnu 


          2
         d B_{i,j,k}   
         -----------(xi,eta,nu) = N_i"(xi)*M_j (eta)*L_k (nu)
           dxi dxi

          2
         d B_{i,j,k}   
         -----------(xi,eta,nu) = N_i (xi)*M_j"(eta)*L_k (nu)
          deta deta 

          2
         d B_{i,j,k}   
         -----------(xi,eta,nu) = N_i (xi)*M_j (eta)*L_k"(nu)
           dnu dnu 

          2
         d B_{i,j,k}   
         -----------(xi,eta,nu) = N_i'(xi)*M_j'(eta)*L_k (nu)
           dxi deta

          2
         d B_{i,j,k}   
         -----------(xi,eta,nu) = N_i'(xi)*M_j (eta)*L_k'(nu)
           dxi dnu 

          2
         d B_{i,j,k}   
         -----------(xi,eta,nu) = N_i (xi)*M_j'(eta)*L_k'(nu)
           deta dnu 


\endverbatim

PART  II: PROJECTING TO NURBS SHAPE FUNCTIONS

\verbatim

      
                        w_{i,j,k} * B_{i,j,k}(xi,eta,nu)
                
       N_{i,j,k} = ---------------------------------------
                   +----+				
                    \				
                     +    w_{s,t,r} * B_{s,t,r}(xi,eta,nu)
      		    /				
      		   +----+				
      		   s,t,r                           
      
\endverbatim

First and second derivatives are evaluated using 
the chain rule. The second derivatives are already 
quite painful to compute...


PART III: TRANSFORMING DERIVATIVES FROM xi/eta TO u/v

\verbatim

              +   dxi   dxi   dxi -+   +-  du           -+
              |   ---   ---   ---  |   |  ---   0    0   |
              |    du    dv    dw  |   |  2.0            |
              |                    |   |                 |
              |  deta  deta  deta  |   |        dv       |
       Psi' = |  ----  ----  ----  | = |   0   ---   0   |
              |   du    dv    dw   |   |       2.0       |
              |                    |   |                 |
              |   dnu   dnu   dnu  |   |             dw  |
              |  ----  ----  ----  |   |   0    0   ---  |
              +-  du    dv    dw  -+   +-           2.0 -+
    
\endverbatim
\verbatim

   bezier_shape_deriv2(.,...): 

   +-------------+---    ......   ---+-------------------+
   |  2          |                   |  2                |
   | d B_{0,0,0} |                   | d B_{N-1,M-1,L-1} |
   | ----------- |                   | ----------------- |
   |    dxi dxi  |                   |       dxi dxi     |
   |             |                   |                   |
   +-------------+---    ......   ---+-------------------+
   |             |                   |                   |
   |  2          |                   |  2                |
   | d B_{0,0,0} |                   | d B_{N-1,M-1,L-1} |
   | ----------- |                   | ----------------- |
   |  deta deta  |                   |      deta deta    |
   |             |                   |                   |
   +-------------+---    ......   ---+-------------------+
   |             |                   |                   |
   |  2          |                   |  2                |
   | d B_{0,0,0} |                   | d B_{N-1,M-1,L-1} |
   | ----------- |                   | ----------------- |
   |   dnu dnu   |                   |       dnu dnu     |
   |             |                   |                   |
   +-------------+---   ......    ---+-------------------+
   |             |                   |                   |
   |  2          |                   |  2                |
   | d B_{0,0,0} |                   | d B_{N-1,M-1,L-1} |
   | ----------- |                   | ----------------- |
   |   dxi deta  |                   |      dxi deta     |
   |             |                   |                   |
   +-------------+---   ......    ---+-------------------+
   |             |                   |                   |
   |  2          |                   |  2                |
   | d B_{0,0,0} |                   | d B_{N-1,M-1,L-1} |
   | ----------- |                   | ----------------- |
   |   dxi dnu   |                   |       dxi dnu     |
   |             |                   |                   |
   +-------------+---   ......    ---+-------------------+
   |             |                   |                   |
   |  2          |                   |  2                |
   | d B_{0,0,0} |                   | d B_{N-1,M-1,L-1} |
   | ----------- |                   | ----------------- |
   |   deta dnu  |                   |      deta dnu     |
   +-------------+---   .......   ---+-------------------+


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

void nurbs_get_3D_funct_deriv_deriv2(
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

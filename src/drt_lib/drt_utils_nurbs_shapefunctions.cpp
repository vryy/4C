/*!----------------------------------------------------------------------
\file drt_utils_nurbs_shapefunctions.cpp

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "drt_utils_nurbs_shapefunctions.H"

using namespace std;

/*

Evaluate basis functions of nurbs basis functions (line element version).

*/

void DRT::NURBS::UTILS::nurbs_get_1D_funct(
  blitz::Array<double,1>                 & nurbs_shape_funct  ,
  const double                           & u                  ,
  const blitz::Array<double, 1>          & knots              ,
  const blitz::Array<double, 1>          & weights            ,
  const DRT::Element::DiscretizationType & distype            )
{

  int degree = 0;

  switch (distype)
  {
  case DRT::Element::nurbs2:
  {
    degree=1;
    
    break;
  }
  case DRT::Element::nurbs3:
  {
    degree=2;
    
    break;
  }
  default:
  {
    dserror("Unknown distype for nurbs line element evaluation\n");
  }
  }
 
  // size is the number of control points/basis 
  // functions of this element
  const int size   = (degree+1);

  // Gausspoint is in [-1:1]
  //
  // has to be mapped into knotvector evaluation interval
  
  /*
       reference element          element in knotspan
      
                            Psi                     
        --+XXX+XXX+-->    ------>   ****+XXX+**** --->

         -1      +1   u                             xi
                                                 
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
  */
    
  // the jacobian matrix of Psi is of the very 
  // simple form
  //
  //              Ju
  //              |
  //              |
  //          +-  du -+
  //   Psi' = |  ---  |
  //          +- 2.0 -+
  //                 
  //                 
  //                 
  //

  double du = knots(degree+1)-knots(degree);

  double Ju =du/2.0;

  // get mapped coordinates
  
  double xi = knots(degree)+(u+1.0)*Ju;
  
  // -----------------------------------------------
  //  PART I: EVALUATION OF  BEZIER SHAPE FUNCTIONS
  // -----------------------------------------------
  
  blitz::Array<double,1> bezier_shape_funct(size);
    
  // allocate bspline polynomials for both direction
  DRT::NURBS::UTILS::BsplinePolynomial bspline_xi (degree,knots);
  
  // get temporary double
  double bspline_xi_value;
  
  //          
  //  Bezier basis function 
  //          |
  //          |
  //       B_{i}(xi) = N_i(xi)
  //                     |       
  //                     |       
  //            bspline_u_value  
  //
  
  // loop all basis functions (corresponding to 
  // control points)
  for(int rr=0;rr<degree+1;++rr)
  {
    // in first direction:
    
    // get bsplinevalue
    bspline_xi.EvaluateBspline(bspline_xi_value,
			       xi              ,
			       rr              );

    // add value to bezier_shape_funct
    bezier_shape_funct(rr)=bspline_xi_value;
  }

  // -----------------------------------------------
  //  PART II: PROJECTING TO NURBS SHAPE FUNCTIONS
  // -----------------------------------------------
  
  // alloc temporary doubles, initialise to zero
  double                 sum_funct_weight;
  
  sum_funct_weight=0;
  
  /*loop all control points, compute sums
      
                         +----+
                          \
       sum_funct_weight =  +    w_{i} * B_{i}(xi)
                          /
                         +----+
  */

  for(int rr=0;rr<size;++rr)
  {
    sum_funct_weight+=weights(rr)*bezier_shape_funct(rr);
  }
  
  /* Compute Nurbs basis functions
      
      
                    w_{i} * B_{i}(xi)
                
       N_{i} = ------------------------
               +----+				
                \				
                 +    w_{k} * B_{k}(xi)
      		/				
      	       +----+				
      		  k                           
  */

  // loop all basis functions
  for(int rr=0;rr<size;++rr)
  {    
    // get shape function
    nurbs_shape_funct(rr)=weights(rr)*bezier_shape_funct(rr)/sum_funct_weight;
  }

  return;
}


/*

Evaluate basis functions and derivatives of nurbs basis functions (1d).

*/
void DRT::NURBS::UTILS::nurbs_get_1D_funct_deriv(
  blitz::Array<double,1>                 & nurbs_shape_funct,
  blitz::Array<double,1>                 & nurbs_shape_deriv,
  const double                           & u                ,
  const blitz::Array<double, 1>          & knots            ,
  const blitz::Array<double, 1>          & weights          ,
  const DRT::Element::DiscretizationType & distype
  )
{
  int degree = 0;

  switch (distype)
  {
  case DRT::Element::nurbs2:
  {
    degree=1;

    break;
  }
  case DRT::Element::nurbs3:
  {
    degree=2;
    
    break;
  }
  default:
  {
    dserror("Unknown distype for nurbs line element evaluation\n");
  }
  }
 
  // size is the number of control points/basis 
  // functions of this element
  const int size   = (degree+1);

  // Gausspoint is in [-1:1]
  //
  // has to be mapped into knotvector evaluation interval
  
  /*
       reference element          element in knotspan
      
                            Psi                     
        --+XXX+XXX+-->    ------>   ****+XXX+**** --->

         -1      +1   u                             xi
                                                 
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
  */
    
  // the jacobian matrix of Psi is of the very 
  // simple form
  //
  //              Ju
  //              |
  //              |
  //          +-  du -+
  //   Psi' = |  ---  |
  //          +- 2.0 -+
  //                 
  //                 
  //                 
  //

  double du = knots(degree+1)-knots(degree);

  double Ju =du/2.0;

  // get mapped coordinates
  
  double xi = knots(degree)+(u+1.0)*Ju;
  
  // -----------------------------------------------
  //  PART I: EVALUATION OF  BEZIER SHAPE FUNCTIONS
  // -----------------------------------------------
  
  blitz::Array<double,1> bezier_shape_funct(size);
  blitz::Array<double,1> bezier_shape_deriv(size);
    
  // allocate bspline polynomials for both direction
  DRT::NURBS::UTILS::BsplinePolynomial bspline_xi (degree,knots);
  
  // get temporary doubles
  double bspline_xi_value;
  double bspline_xi_derivative;
  
  //          
  //  Bezier basis function 
  //          |
  //          |
  //       B_{i}(xi,eta) = N_i(xi)
  //                         |       
  //                         |       
  //                bspline_u_value  
  //
  
  //      dB_{i}   
  //      ------(xi,eta) = N_i'(xi)
  //       dxi 
  //
  
  // loop all basis functions (corresponding to 
  // control points)
  for(int rr=0;rr<degree+1;++rr)
  {
    // in first direction:
    
    // get bsplinevalue and derivative
    bspline_xi.EvaluateBsplineAndDeriv(bspline_xi_value,
				       bspline_xi_derivative,
				       xi,
				       rr);
      
    // add value to bezier_shape_funct
    bezier_shape_funct(rr)=bspline_xi_value;
      
    // add values to bezier_shape_deriv
    bezier_shape_deriv(rr)=bspline_xi_derivative;
  }

  // -----------------------------------------------
  //  PART II: PROJECTING TO NURBS SHAPE FUNCTIONS
  // -----------------------------------------------
  
  // alloc temporary doubles, initialise to zero
  double sum_funct_weight;
  double sum_deriv_weight;
  
  sum_funct_weight=0;
  sum_deriv_weight=0;
  
  /*loop all control points, compute sums
      
                         +----+
                          \
       sum_funct_weight =  +    w_{i} * B_{i}(xi)
                          /
                         +----+
      
      
                         +----+
                          \             dB'_{i}
       sum_deriv_weight =  +    w_{i} * -------(xi)
                          /               dxi
                         +----+
  */

  for(int rr=0;rr<size;++rr)
  {
    sum_funct_weight+=weights(rr)*bezier_shape_funct(rr);
    
    sum_deriv_weight+=weights(rr)*bezier_shape_deriv(rr);
  }
  
  /* Compute Nurbs basis functions
      
      
                  w_{i} * B_{i}(xi)
               
       N_{i} = ------------------------
               +----+				
                \				
                 +    w_{k} * B_{k}(xi)
      	        /				
      	       +----+				
      	          k                           
  */

  // Nurbs derivatives are defined by the chain rule
  //
  //                         +- +----+          -+               +- +----+             -+ 
  //                         |   \               |               |   \	         dB_{k} |
  //                  dB_{i} |    +  w_{k}*B_{k} |               |    +    w_{k}*------ |
  //            w_{i}*------*|   /               | - w_{i}*B_{i}*|   /             dxi  |
  //	  	       dxi   |  +----+           |               |  +----+              |
  // dN_{i}                  +-   k             -+               +-   k                -+
  // ------ = ----------------------------------------------------------------------------
  //  dxi                         +- +----+                  -+ 2
  //                              |   \                       |
  //                              |    +    w_{k} * B_{k}(xi) |
  //		                  |   /                       |
  //		                  |  +----+                   |
  //		                  +-   k                     -+

  // loop all basis functions
  for(int rr=0;rr<size;++rr)
  {    
    // get shape function
    nurbs_shape_funct(rr)=weights(rr)*bezier_shape_funct(rr)/sum_funct_weight;
    
    // loop directions to compute derivatives
    nurbs_shape_deriv(rr)
      =
      bezier_shape_deriv(rr)*sum_funct_weight
      -
      bezier_shape_funct(rr)*sum_deriv_weight;
      
      nurbs_shape_deriv(rr)*=weights(rr)/(sum_funct_weight*sum_funct_weight);
  }
    
  // -----------------------------------------------
  //  PART III: TRANSFORMING DERIVATIVES FROM xi
  //            TO u
  // -----------------------------------------------
   
  // we already know the jacobian matrix of psi
  //
  //          +- dxi -+   +-  du -+
  //   Psi' = |  ---  | = |  ---  |
  //          +-  du -+   +- 2.0 -+
  //
  // we will obtain the derivatives with respect to  
  // u,v just by multiplication  
  
  // loop all basis function derivatives
  for(int rr=0;rr<size;++rr)
  {
    nurbs_shape_deriv(rr)*=Ju;
  }    

  return;
}

/*

Evaluate basis functions of 2d nurbs basis functions.

*/

void DRT::NURBS::UTILS::nurbs_get_2D_funct(
    blitz::Array<double,1>                 & nurbs_shape_funct  ,
    const blitz::Array<double,1>           & uv                 ,
    const vector<blitz::Array<double, 1> > & knots              ,
    const blitz::Array<double, 1>          & weights            ,
    const DRT::Element::DiscretizationType & distype            )
{

  int degree = 1;

  switch (distype)
  {
  case DRT::Element::nurbs4:
  {
    degree=1;

    break;
  }
  case DRT::Element::nurbs9:
  {
    degree=2;

    break;
  }
  default:
  {
    dserror("Unknown distype for nurbs element evaluation\n");
  }
  }
  
  // size is the number of control points/basis 
  // functions of this element
  const int size   = (degree+1)*(degree+1);

  // Gausspoint is in [-1:1]x[-1:1]
  //
  // has to be mapped into knotvector evaluation interval
  
  /*
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
  */
    
  // the jacobian matrix of Psi is of the very 
  // simple form
  //
  //               Ju
  //               |
  //               |
  //          +-   du       -+
  //          |   ---   0    |
  //          |   2.0        |
  //   Psi' = |              |
  //          |         dv   |
  //          |    0   ---   |
  //          +-       2.0  -+
  //                    |
  //                    |
  //                    Jv
  //

  double du = (knots[0])(degree+1)-(knots[0])(degree);
  double dv = (knots[1])(degree+1)-(knots[1])(degree);

  double Ju =du/2.0;
  double Jv =dv/2.0;

  // get mapped coordinates
  
  double xi = (knots[0])(degree)+((uv(0)+1))*Ju;
  double eta= (knots[1])(degree)+((uv(1)+1))*Jv;
  
  // -----------------------------------------------
  //  PART I: EVALUATION OF  BEZIER SHAPE FUNCTIONS
  // -----------------------------------------------
  
  blitz::Array<double,1> bezier_shape_funct(  size);
    
  // allocate bspline polynomials for both direction
  DRT::NURBS::UTILS::BsplinePolynomial bspline_xi (degree,knots[0]);
  DRT::NURBS::UTILS::BsplinePolynomial bspline_eta(degree,knots[1]);
  
  // get temporary doubles
  double bspline_xi_value;
  double bspline_eta_value;
  
  // define temporary int variable to compute the 
  // number of the basis function from i,j
  int id;
  
  //          
  //  Bezier basis function 
  //          |
  //          |
  //       B_{i,j}(xi,eta) = N_i(xi)*M_j(eta)
  //                           |       |
  //                           |       |
  //                  bspline_u_value  |
  //                                   |
  //                             bspline_v_value
  //
  
  // loop all basis functions (corresponding to 
  // control points)
  for(int rr=0;rr<degree+1;++rr)
  {
    // in first direction:
    
    // get bsplinevalue and derivative
    bspline_xi.EvaluateBspline(bspline_xi_value,
			       xi,
			       rr);
    
    for(int mm=0;mm<degree+1;++mm)
    {
      // in second direction:
      
      // get bsplinevalue and derivative
      bspline_eta.EvaluateBspline(bspline_eta_value,
				  eta,
				  mm);

      // get the number of the basis function
      id = rr+mm*(degree+1);
      
      // add value to bezier_shape_funct
      bezier_shape_funct(id)=bspline_xi_value*bspline_eta_value;
    }
  }

  // -----------------------------------------------
  //  PART II: PROJECTING TO NURBS SHAPE FUNCTIONS
  // -----------------------------------------------
  
  // alloc temporary doubles, initialise to zero
  double                 sum_funct_weight;
  
  sum_funct_weight=0;
  
  /*loop all control points, compute sums
      
                         +----+
                          \
       sum_funct_weight =  +    w_{i,j} * B_{i,j}(xi,eta)
                          /
                         +----+
  */

  for(int rr=0;rr<size;++rr)
  {
    sum_funct_weight+=weights(rr)*bezier_shape_funct(rr);
  }
  
  /* Compute Nurbs basis functions
      
      
                      w_{i,j} * B_{i,j}(xi,eta)
                
       N_{i,j} = ----------------------------------
                  +----+				
                   \				
                    +    w_{k,l} * B_{k,l}(xi,eta)
      		   /				
      		  +----+				
      		    k,l                           
  */

  // loop all basis functions
  for(int rr=0;rr<size;++rr)
  {    
    // get shape function
    nurbs_shape_funct(rr)=weights(rr)*bezier_shape_funct(rr)/sum_funct_weight;
  }

  return;
}

/*

Evaluate basis functions and derivatives (with 
respect to the coordinates of the reference element) 
of 2d nurbs basis functions.

*/

void DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv(
    blitz::Array<double,1>                 & nurbs_shape_funct  ,
    blitz::Array<double,2>                 & nurbs_shape_deriv  ,
    const blitz::Array<double,1>           & uv                 ,
    const vector<blitz::Array<double, 1> > & knots              ,
    const blitz::Array<double, 1>          & weights            ,
    const DRT::Element::DiscretizationType & distype            )
{

  int degree = 1;

  switch (distype)
  {
  case DRT::Element::nurbs4:
  {
    degree=1;

    break;
  }
  case DRT::Element::nurbs9:
  {
    degree=2;

    break;
  }
  default:
  {
    dserror("Unknown distype for nurbs element evaluation\n");
  }
  }
  
  // size is the number of control points/basis 
  // functions of this element
  const int size   = (degree+1)*(degree+1);

  // Gausspoint is in [-1:1]x[-1:1]
  //
  // has to be mapped into knotvector evaluation interval
  
  /*
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
  */
    
  // the jacobian matrix of Psi is of the very 
  // simple form
  //
  //               Ju
  //               |
  //               |
  //          +-   du       -+
  //          |   ---   0    |
  //          |   2.0        |
  //   Psi' = |              |
  //          |         dv   |
  //          |    0   ---   |
  //          +-       2.0  -+
  //                    |
  //                    |
  //                    Jv
  //

  double du = (knots[0])(degree+1)-(knots[0])(degree);
  double dv = (knots[1])(degree+1)-(knots[1])(degree);

  double Ju =du/2.0;
  double Jv =dv/2.0;

  // get mapped coordinates
  
  double xi = (knots[0])(degree)+((uv(0)+1))*Ju;
  double eta= (knots[1])(degree)+((uv(1)+1))*Jv;
  
  // -----------------------------------------------
  //  PART I: EVALUATION OF  BEZIER SHAPE FUNCTIONS
  // -----------------------------------------------
  
  blitz::Array<double,1> bezier_shape_funct(  size);
  blitz::Array<double,2> bezier_shape_deriv(2,size);
    
  // allocate bspline polynomials for both direction
  DRT::NURBS::UTILS::BsplinePolynomial bspline_xi (degree,knots[0]);
  DRT::NURBS::UTILS::BsplinePolynomial bspline_eta(degree,knots[1]);
  
  // get temporary doubles
  double bspline_xi_value;
  double bspline_xi_derivative;
  double bspline_eta_value;
  double bspline_eta_derivative;
  
  // define temporary int variable to compute the 
  // number of the basis function from i,j
  int id;
  
  //          
  //  Bezier basis function 
  //          |
  //          |
  //       B_{i,j}(xi,eta) = N_i(xi)*M_j(eta)
  //                           |       |
  //                           |       |
  //                  bspline_u_value  |
  //                                   |
  //                             bspline_v_value
  //
  
  //      dB_{i,j}   
  //      --------(xi,eta) = N_i'(xi)*M_j (eta)
  //         dxi 
  //
  //      dB_{i,j}   
  //      --------(xi,eta) = N_i (xi)*M_j'(eta)
  //        deta 
  
  // loop all basis functions (corresponding to 
  // control points)
  for(int rr=0;rr<degree+1;++rr)
  {
    // in first direction:
    
    // get bsplinevalue and derivative
    bspline_xi.EvaluateBsplineAndDeriv(bspline_xi_value,
				       bspline_xi_derivative,
				       xi,
				       rr);
    
    for(int mm=0;mm<degree+1;++mm)
    {
      // in second direction:
      
      // get bsplinevalue and derivative
      bspline_eta.EvaluateBsplineAndDeriv(bspline_eta_value,
					  bspline_eta_derivative,
					  eta,
					  mm);

      // get the number of the basis function
      id = rr+mm*(degree+1);
      
      // add value to bezier_shape_funct
      bezier_shape_funct(id)=bspline_xi_value*bspline_eta_value;
      
      // add values to bezier_shape_deriv
      bezier_shape_deriv(0,id)=bspline_xi_derivative*bspline_eta_value;
      bezier_shape_deriv(1,id)=bspline_xi_value*bspline_eta_derivative;
    }
  }

  // -----------------------------------------------
  //  PART II: PROJECTING TO NURBS SHAPE FUNCTIONS
  // -----------------------------------------------
  
  // alloc temporary doubles, initialise to zero
  double                 sum_funct_weight;
  blitz::Array<double,1> sum_deriv_weight(2);
  
  sum_funct_weight=0;
  sum_deriv_weight=0;
  
  /*loop all control points, compute sums
      
                         +----+
                          \
       sum_funct_weight =  +    w_{i,j} * B_{i,j}(xi,eta)
                          /
                         +----+
      
      
                            +----+
                             \               dB'_{i,j}
       sum_deriv_weight[.] =  +    w_{i,j} * ---------(xi,eta)
                             /                  d.
                            +----+
  */

  for(int rr=0;rr<size;++rr)
  {
    sum_funct_weight+=weights(rr)*bezier_shape_funct(rr);
    
    for(int mm=0;mm<2;++mm)
    {
      sum_deriv_weight(mm)+=weights(rr)*bezier_shape_deriv(mm,rr);
    }
  }
  
  /* Compute Nurbs basis functions
      
      
                      w_{i,j} * B_{i,j}(xi,eta)
                
       N_{i,j} = ----------------------------------
                  +----+				
                   \				
                    +    w_{k,l} * B_{k,l}(xi,eta)
      		   /				
      		  +----+				
      		    k,l                           
  */

  // Nurbs derivatives are defined by the chain rule
  //
  //                             +- +----+              -+                   +- +----+                 -+ 
  //                             |   \                   |                   |   \	       dB_{k,l} |
  //                    dB_{i,j} |    +  w_{k,l}*B_{k,l} |                   |    +    w_{k,l}*-------- |
  //            w_{i,j}*--------*|   /                   | - w_{i,j}*B_{i,j}*|   /                dxi   |
  //	  	            dxi  |  +----+               |                   |  +----+                  |
  // dN_{i,j}               	 +-   k,l               -+                   +-   k,l                  -+
  // -------- = -----------------------------------------------------------------------------------------
  //   dxi                             +- +----+                          -+ 2
  //                                   |   \                               |
  //                                   |    +    w_{k,l} * B_{k,l}(xi,eta) |
  //		                       |   /                               |
  //		                       |  +----+                           |
  //		                       +-   k,l                           -+

  // loop all basis functions
  for(int rr=0;rr<size;++rr)
  {    
    // get shape function
    nurbs_shape_funct(rr)=weights(rr)*bezier_shape_funct(rr)/sum_funct_weight;
    
    // loop directions to compute derivatives
    for(int mm=0;mm<2;++mm)
    {
      nurbs_shape_deriv(mm,rr)
	=
	bezier_shape_deriv(mm,rr)*sum_funct_weight
	-
	bezier_shape_funct(rr)*sum_deriv_weight(mm);
      
      nurbs_shape_deriv(mm,rr)*=weights(rr)/(sum_funct_weight*sum_funct_weight);
    }
  }
    
  // -----------------------------------------------
  //  PART III: TRANSFORMING DERIVATIVES FROM xi/eta
  //            TO u/v
  // -----------------------------------------------
   
  // we already know the jacobian matrix of psi
  //
  //          +-   dxi   dxi  -+   +-   du       -+
  //          |    ---   ---   |   |   ---   0    |
  //          |     du    dv   |   |   2.0        |
  //   Psi' = |                | = |              |
  //          |   deta  deta   |   |         dv   |
  //          |   ----  ----   |   |    0   ---   |
  //          +-   du    dv   -+   +-       2.0  -+
  //
  // we will obtain the derivatives with respect to  
  // u,v just by multiplication  
  
  // loop all basis function derivatives
  for(int rr=0;rr<size;++rr)
  {
    nurbs_shape_deriv(0,rr)*=Ju;
    nurbs_shape_deriv(1,rr)*=Jv;
  }    

  return;
}


/*

Evaluate basis functions, first and second derivatives 
(with respect to the coordinates of the reference 
element) of 2d nurbs basis functions.

*/

void DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv_deriv2(
    blitz::Array<double,1>                 & nurbs_shape_funct  ,
    blitz::Array<double,2>                 & nurbs_shape_deriv  ,
    blitz::Array<double,2>                 & nurbs_shape_deriv2 ,
    const blitz::Array<double,1>           & uv                 ,
    const vector<blitz::Array<double, 1> > & knots              ,
    const blitz::Array<double, 1>          & weights            ,
    const DRT::Element::DiscretizationType & distype            )
{

  int degree = 1;

  switch (distype)
  {
  case DRT::Element::nurbs4:
  {
    degree=1;

    break;
  }
  case DRT::Element::nurbs9:
  {
    degree=2;

    break;
  }
  default:
  {
    dserror("Unknown distype for nurbs element evaluation\n");
  }
  }
  
  // size is the number of control points/basis 
  // functions of this element
  const int size   = (degree+1)*(degree+1);

  // Gausspoint is in [-1:1]x[-1:1]
  //
  // has to be mapped into knotvector evaluation interval
  
  /*
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
  */
    
  // the jacobian matrix of Psi is of the very 
  // simple form
  //
  //               Ju
  //               |
  //               |
  //          +-   du       -+
  //          |   ---   0    |
  //          |   2.0        |
  //   Psi' = |              |
  //          |         dv   |
  //          |    0   ---   |
  //          +-       2.0  -+
  //                    |
  //                    |
  //                    Jv
  //

  double du = (knots[0])(degree+1)-(knots[0])(degree);
  double dv = (knots[1])(degree+1)-(knots[1])(degree);

  double Ju =du/2.0;
  double Jv =dv/2.0;

  // get mapped coordinates
  
  double xi = (knots[0])(degree)+((uv(0)+1))*Ju;
  double eta= (knots[1])(degree)+((uv(1)+1))*Jv;
  
  // -----------------------------------------------
  //  PART I: EVALUATION OF  BEZIER SHAPE FUNCTIONS
  // -----------------------------------------------
  
  blitz::Array<double,1> bezier_shape_funct (  size);
  blitz::Array<double,2> bezier_shape_deriv (2,size);
  blitz::Array<double,2> bezier_shape_deriv2(3,size);
    
  // allocate bspline polynomials for both direction
  DRT::NURBS::UTILS::BsplinePolynomial bspline_xi (degree,knots[0]);
  DRT::NURBS::UTILS::BsplinePolynomial bspline_eta(degree,knots[1]);
  
  // get temporary doubles for derivatives and 
  // values of the above devlared polynomials
  double bspline_xi_value;
  double bspline_eta_value;

  double bspline_xi_derivative;
  double bspline_eta_derivative;

  double bspline_xi_deriv2;
  double bspline_eta_deriv2;
  
  // define temporary int variable to compute the 
  // number of the basis function from i,j
  int id;
  
  //          
  //  Bezier basis function 
  //          |
  //          |
  //       B_{i,j}(xi,eta) = N_i(xi)*M_j(eta)
  //                           |       |
  //                           |       |
  //                  bspline_u_value  |
  //                                   |
  //                             bspline_v_value
  //
  
  //      dB_{i,j}   
  //      --------(xi,eta) = N_i'(xi)*M_j (eta)
  //         dxi 
  //
  //      dB_{i,j}   
  //      --------(xi,eta) = N_i (xi)*M_j'(eta)
  //        deta 

  //       2
  //      d B_{i,j}   
  //      ---------(xi,eta) = N_i"(xi)*M_j (eta)
  //           2
  //        dxi 
  //
  //       2
  //      d B_{i,j}   
  //      ---------(xi,eta) = N_i (xi)*M_j"(eta)
  //           2
  //       deta 
  //
  //       2
  //      d B_{i,j}   
  //      ---------(xi,eta) = N_i'(xi)*M_j'(eta)
  //      dxi deta
  //           

  
  // loop all basis functions (corresponding to 
  // control points)
  for(int rr=0;rr<degree+1;++rr)
  {
    // in first direction:
    
    // get bsplinevalue and derivative
    bspline_xi.EvaluateBsplineFirstAndSecondDeriv(
      bspline_xi_value,
      bspline_xi_derivative,
      bspline_xi_deriv2,
      xi,
      rr);
    
    for(int mm=0;mm<degree+1;++mm)
    {
      // in second direction:
      
      // get bsplinevalue and derivative
      bspline_eta.EvaluateBsplineFirstAndSecondDeriv(
	bspline_eta_value,
	bspline_eta_derivative,
	bspline_eta_deriv2,
	eta,
	mm);

      // get the number of the basis function
      id = rr+mm*(degree+1);
      
      // set value to bezier_shape_funct
      bezier_shape_funct(id)=bspline_xi_value*bspline_eta_value;
      
      // set values to bezier_shape_deriv
      bezier_shape_deriv(0,id)=bspline_xi_derivative*bspline_eta_value;
      bezier_shape_deriv(1,id)=bspline_xi_value*bspline_eta_derivative;

      // set values to bezier_shape_deriv2
      bezier_shape_deriv2(0,id)=bspline_xi_deriv2    *bspline_eta_value     ;
      bezier_shape_deriv2(1,id)=bspline_xi_value     *bspline_eta_deriv2    ;
      bezier_shape_deriv2(2,id)=bspline_xi_derivative*bspline_eta_derivative;
    }
  }

  // -----------------------------------------------
  //  PART II: PROJECTING TO NURBS SHAPE FUNCTIONS
  // -----------------------------------------------
  
  // alloc temporary doubles, initialise to zero
  double                 sum_funct_weight    ;
  blitz::Array<double,1> sum_deriv_weight (2);
  blitz::Array<double,1> sum_deriv2_weight(3);
  
  sum_funct_weight =0;
  sum_deriv_weight =0;
  sum_deriv2_weight=0;

  /*loop all control points, compute sums
      
                         +----+
                          \
       sum_funct_weight =  +    w_{i,j} * B_{i,j}(xi,eta)
                          /
                         +----+
      
      
                            +----+
                             \               dB_{i,j}
       sum_deriv_weight[0] =  +    w_{i,j} * --------(xi,eta)
                             /                  dxi
                            +----+
      
                            +----+
                             \               dB_{i,j}
       sum_deriv_weight[1] =  +    w_{i,j} * --------(xi,eta)
                             /                 deta
                            +----+


      
                             +----+            2
                              \               d B_{i,j}
       sum_deriv2_weight[0] =  +    w_{i,j} * --------(xi,eta)
                              /                     2
                             +----+              dxi
      
                             +----+            2 
                              \               d B_{i,j}
       sum_deriv2_weight[1] =  +    w_{i,j} * ---------(xi,eta)
                              /                     2
                             +----+             deta
      
                             +----+            2
                              \               d B_{i,j}
       sum_deriv2_weight[2] =  +    w_{i,j} * ---------(xi,eta)
                              /               dxi deta
                             +----+

  */

  for(int rr=0;rr<size;++rr)
  {
    sum_funct_weight+=weights(rr)*bezier_shape_funct(rr);
    
    for(int mm=0;mm<2;++mm)
    {
      sum_deriv_weight(mm)+=weights(rr)*bezier_shape_deriv(mm,rr);
      sum_deriv2_weight(mm)+=weights(rr)*bezier_shape_deriv2(mm,rr);
    }
    sum_deriv2_weight(2)+=weights(rr)*bezier_shape_deriv2(2,rr);
  }
  
  /* Compute Nurbs basis functions
      
      
                      w_{i,j} * B_{i,j}(xi,eta)
                
       N_{i,j} = ----------------------------------
                  +----+				
                   \				
                    +    w_{k,l} * B_{k,l}(xi,eta)
      		   /				
      		  +----+				
      		    k,l                           
  */

  /* Nurbs derivatives are defined by the chain rule
    
                                 +- +----+              -+                   +- +----+                 -+ 
                                 |   \                   |                   |   \	       dB_{k,l} |
                        dB_{i,j} |    +  w_{k,l}*B_{k,l} |                   |    +    w_{k,l}*-------- |
                w_{i,j}*--------*|   /                   | - w_{i,j}*B_{i,j}*|   /                dxi   |
    	  	            dxi  |  +----+               |                   |  +----+                  |
     dN_{i,j}               	 +-   k,l               -+                   +-   k,l                  -+
     -------- = -----------------------------------------------------------------------------------------
       dxi                             +- +----+                          -+ 2
                                       |   \                               |
                                       |    +    w_{k,l} * B_{k,l}(xi,eta) |
    		                       |   /                               |
    		                       |  +----+                           |
    		                       +-   k,l                           -+
  */

  /* Nurbs second derivatives are calculated by a 
     second application of the chain rule
    
                
                
                                                         +-
                                                         |
                                                         |
                                                         |
    	  	              w_{i,j}                    |   2
     d N_{i,j}                                           |  d B_{i,j}
     --------- = ------------------------------------- * |  --------- -
     dxi  dxi    +- +----+                          -+   |  dxi  dxi
        k    l   |   \                               |   |     k    l
                 |    +    w_{n,m} * B_{n,m}(xi,eta) |   |
                 |   /                               |   |
                 |  +----+                           |   |
                 +-   n,m                           -+   +-
  
                   +- +----+                           -+
          dB_{i,j} |   \               dB_{n,m}         |
          -------- |    +    w_{n,m} * --------(xi,eta) |
           dxi     |   /                  dxi           |
              k    |  +----+                 l          |
                   +-   n,m                            -+
        - ----------------------------------------------- -
               +- +----+                          -+            
               |   \                               |             
               |    +    w_{n,m} * B_{n,m}(xi,eta) |   
               |   /                               |   
               |  +----+                           |   
               +-   n,m                           -+   


  
                   +- +----+                           -+
          dB_{i,j} |   \               dB_{n,m}         |
          -------- |    +    w_{n,m} * --------(xi,eta) |
           dxi     |   /                  dxi           |
              l    |  +----+                 k          |
                   +-   n,m                            -+
        - ----------------------------------------------- -
               +- +----+                          -+            
               |   \                               |             
               |    +    w_{n,m} * B_{n,m}(xi,eta) |   
               |   /                               |   
               |  +----+                           |   
               +-   n,m                           -+   


                     +- +----+            2               -+
                     |   \               d B_{n,m}         |
                     |    +    w_{n,m} * ---------(xi,eta) |
           B_{i,j} * |   /               dxi  dxi          |
                     |  +----+              k    l         |
                     +-   n,m                             -+
        - -------------------------------------------------- +
                 +- +----+                          -+            
                 |   \                               |             
                 |    +    w_{n,m} * B_{n,m}(xi,eta) |   
                 |   /                               |   
                 |  +----+                           |   
                 +-   n,m                           -+   


                         +- +----+                            -+   +- +----+                            -+  -+
                         |   \               dB_{n,m}          |   |   \               dB_{n,m}          |   |
                         |    +    w_{n,m} * ---------(xi,eta) |   |    +    w_{n,m} * ---------(xi,eta) |   |
           2 * B_{i,j} * |   /                 dxi             | * |   /                 dxi             |   |
                         |  +----+                k            |   |  +----+                l            |   |
                         +-   n,m                             -+   +-   n,m                             -+   |
        + ------------------------------------------------------------------------------------------------   |
                                  +- +----+                          -+ 2                                    |
                                  |   \                               |                                      |
                                  |    +    w_{n,m} * B_{n,m}(xi,eta) |                                      |
                                  |   /                               |                                      |
                                  |  +----+                           |                                      |
                                  +-   n,m                           -+                                     -+
  */


  // loop all basis functions
  for(int rr=0;rr<size;++rr)
  {    
    // get shape function
    nurbs_shape_funct(rr)=weights(rr)*bezier_shape_funct(rr)/sum_funct_weight;
    
    // loop directions to compute derivatives
    for(int mm=0;mm<2;++mm)
    {
      nurbs_shape_deriv(mm,rr)
	=
	bezier_shape_deriv(mm,rr)*sum_funct_weight
	-
	bezier_shape_funct(rr)*sum_deriv_weight(mm);
      
      nurbs_shape_deriv(mm,rr)*=weights(rr)/(sum_funct_weight*sum_funct_weight);
    }

    // we apply a Horner-type scheme for the 
    // multiplication with sum_funct_weight
    nurbs_shape_deriv2(0,rr) =2*bezier_shape_funct(rr)*sum_deriv_weight(0)*sum_deriv_weight(0);
    nurbs_shape_deriv2(0,rr)/=sum_funct_weight;
    nurbs_shape_deriv2(0,rr)-=bezier_shape_funct(rr)  *sum_deriv2_weight(0);
    nurbs_shape_deriv2(0,rr)-=bezier_shape_deriv(0,rr)*sum_deriv_weight(0);
    nurbs_shape_deriv2(0,rr)-=bezier_shape_deriv(0,rr)*sum_deriv_weight(0);
    nurbs_shape_deriv2(0,rr)/=sum_funct_weight;
    nurbs_shape_deriv2(0,rr)+=bezier_shape_deriv2(0,rr);
    nurbs_shape_deriv2(0,rr)*=weights(rr)/sum_funct_weight;

    nurbs_shape_deriv2(1,rr) =2*bezier_shape_funct(rr)*sum_deriv_weight(1)*sum_deriv_weight(1);
    nurbs_shape_deriv2(1,rr)/=sum_funct_weight;
    nurbs_shape_deriv2(1,rr)-=bezier_shape_funct(rr)  *sum_deriv2_weight(1);
    nurbs_shape_deriv2(1,rr)-=bezier_shape_deriv(1,rr)*sum_deriv_weight(1);
    nurbs_shape_deriv2(1,rr)-=bezier_shape_deriv(1,rr)*sum_deriv_weight(1);
    nurbs_shape_deriv2(1,rr)/=sum_funct_weight;
    nurbs_shape_deriv2(1,rr)+=bezier_shape_deriv2(1,rr);
    nurbs_shape_deriv2(1,rr)*=weights(rr)/sum_funct_weight;

    nurbs_shape_deriv2(2,rr) =2*bezier_shape_funct(rr)*sum_deriv_weight(0)*sum_deriv_weight(1);
    nurbs_shape_deriv2(2,rr)/=sum_funct_weight;
    nurbs_shape_deriv2(2,rr)-=bezier_shape_funct(rr)  *sum_deriv2_weight(2);
    nurbs_shape_deriv2(2,rr)-=bezier_shape_deriv(0,rr)*sum_deriv_weight(1);
    nurbs_shape_deriv2(2,rr)-=bezier_shape_deriv(1,rr)*sum_deriv_weight(0);
    nurbs_shape_deriv2(2,rr)/=sum_funct_weight;
    nurbs_shape_deriv2(2,rr)+=bezier_shape_deriv2(2,rr);
    nurbs_shape_deriv2(2,rr)*=weights(rr)/sum_funct_weight;
  }

  // -----------------------------------------------
  //  PART III: TRANSFORMING DERIVATIVES FROM xi/eta
  //            TO u/v
  // -----------------------------------------------
   
  // we already know the jacobian matrix of psi
  //
  //          +-   dxi   dxi  -+   +-   du       -+
  //          |    ---   ---   |   |   ---   0    |
  //          |     du    dv   |   |   2.0        |
  //   Psi' = |                | = |              |
  //          |   deta  deta   |   |         dv   |
  //          |   ----  ----   |   |    0   ---   |
  //          +-   du    dv   -+   +-       2.0  -+
  //
  // we will obtain the derivatives with respect to  
  // u,v just by multiplication  
  
  // loop all basis function derivatives
  for(int rr=0;rr<size;++rr)
  {
    nurbs_shape_deriv(0,rr)*=Ju;
    nurbs_shape_deriv(1,rr)*=Jv;

    nurbs_shape_deriv2(0,rr)*=Ju*Ju;
    nurbs_shape_deriv2(1,rr)*=Jv*Jv;
    nurbs_shape_deriv2(2,rr)*=Ju*Jv;
  }    

  return;
}

/*

  Evaluate 2d basis functions of nurbs basis functions.

*/

void DRT::NURBS::UTILS::nurbs_get_3D_funct(
    blitz::Array<double,1>                 & nurbs_shape_funct,
    const blitz::Array<double,1>           & uv               ,
    const vector<blitz::Array<double, 1> > & knots            ,
    const blitz::Array<double, 1>          & weights          ,
    const DRT::Element::DiscretizationType & distype          )
{

  int degree = 1;

  switch (distype)
  {
  case DRT::Element::nurbs8:
  {
    degree=1;

    break;
  }
  case DRT::Element::nurbs27:
  {
    degree=2;

    break;
  }
  default:
  {
    dserror("Unknown distype for nurbs element evaluation\n");
  }
  }
  
  // size is the number of control points/basis 
  // functions of this element
  const int size   = (degree+1)*(degree+1)*(degree+1);

  // Gausspoint is in [-1:1]x[-1:1]x[-1:1]
  //
  // has to be mapped into knotvector evaluation interval
  
  /*
       reference element                element in knotspan
                                  
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

  */
    
  // the jacobian matrix of Psi is of the very 
  // simple form
  //
  //
  //                                        Ju     Jv
  //                                        |       |
  //                                        |       |		 				                
  //           +   dxi   dxi   dxi -+   +-  du      |    -+  
  //           |   ---   ---   ---  |   |  ---   0  | 0   | 
  //           |    du    dv    dw  |   |  2.0      |     |  
  //           |                    |   |          /      |  
  //           |  deta  deta  deta  |   |        dv       |  
  //    Psi' = |  ----  ----  ----  | = |   0   ---   0   |  
  //           |   du    dv    dw   |   |       2.0       |  
  //           |                    |   |                 |  
  //           |   dnu   dnu   dnu  |   |             dw  |  
  //           |  ----  ----  ----  |   |   0    0   ---  |
  //           +-  du    dv    dw  -+   +-           2.0 -+
  //                                                  |
  //                                                  |
  //                                                  Jw

  double du = (knots[0])(degree+1)-(knots[0])(degree);
  double dv = (knots[1])(degree+1)-(knots[1])(degree);
  double dw = (knots[2])(degree+1)-(knots[2])(degree);

  double Ju =du/2.0;
  double Jv =dv/2.0;
  double Jw =dw/2.0;

  // get mapped coordinates
  
  double xi = (knots[0])(degree)+((uv(0)+1))*Ju;
  double eta= (knots[1])(degree)+((uv(1)+1))*Jv;
  double nu = (knots[2])(degree)+((uv(2)+1))*Jw;
  
  // -----------------------------------------------
  //  PART I: EVALUATION OF  BEZIER SHAPE FUNCTIONS
  // -----------------------------------------------
  
  blitz::Array<double,1> bezier_shape_funct (  size);
    
  // allocate bspline polynomials for both direction
  DRT::NURBS::UTILS::BsplinePolynomial bspline_xi (degree,knots[0]);
  DRT::NURBS::UTILS::BsplinePolynomial bspline_eta(degree,knots[1]);
  DRT::NURBS::UTILS::BsplinePolynomial bspline_nu (degree,knots[2]);

  // get temporary doubles for derivatives and 
  // values of the above devlared polynomials
  double bspline_xi_value;
  double bspline_eta_value;
  double bspline_nu_value;
  
  // define temporary int variable to compute the 
  // number of the basis function from i,j,k
  int id;
  

  /*
       Bezier basis function 
                |
                |
           B_{i,j,k}(xi,eta,nu) = N_i(xi) *M_j (eta)*L_k (nu)
                                  |        |         |
                                  |        |         | 
                           bspline_u_value |         |
                                           |         |
                                    bspline_v_value  |
                                                     |
                                               bspline_w_value
  */
    
  // loop all basis functions (corresponding to 
  // control points)
  for(int rr=0;rr<degree+1;++rr)
  {
    // in first direction:
    
    // get bsplinevalue and derivative
    bspline_xi.EvaluateBspline(
      bspline_xi_value,
      xi,
      rr);
    
    for(int mm=0;mm<degree+1;++mm)
    {
      // in second direction:
      
      // get bsplinevalue and derivative
      bspline_eta.EvaluateBspline(
	bspline_eta_value,
	eta,
	mm);

      for(int nn=0;nn<degree+1;++nn)
      {
	// in third direction:
            
	// get bsplinevalue and derivative
	bspline_nu.EvaluateBspline(
	  bspline_nu_value,
	  nu,
	  nn);

	// get the number of the basis function
	id = rr+(degree+1)*(mm+nn*(degree+1));
      
	// set value to bezier_shape_funct
	bezier_shape_funct(id)=bspline_xi_value*bspline_eta_value*bspline_nu_value;
      }
    }
  }

  // -----------------------------------------------
  //  PART II: PROJECTING TO NURBS SHAPE FUNCTIONS
  // -----------------------------------------------
  
  // alloc temporary doubles, initialise to zero
  double                 sum_funct_weight;
  
  sum_funct_weight =0;

  /*loop all control points, compute sums
      
                         +----+
                          \
       sum_funct_weight =  +    w_{i,j,k} * B_{i,j,k}(xi,eta,nu)
                          /
                         +----+

  */

  for(int rr=0;rr<size;++rr)
  {
    sum_funct_weight+=weights(rr)*bezier_shape_funct(rr);
  }
  
  // loop all basis functions
  for(int rr=0;rr<size;++rr)
  {    

    /* Compute Nurbs basis functions
      
      
                      w_{i,j,k} * B_{i,j,k}(xi,eta,nu)
                
       N_{i,j,k} = ---------------------------------------
                   +----+				
                    \				
                     +    w_{r,s,t} * B_{r,s,t}(xi,eta,nu)
      		    /				
      		   +----+				
      		    r,s,t                           
    */

    nurbs_shape_funct(rr)=weights(rr)*bezier_shape_funct(rr)/sum_funct_weight;
  }

  return;
}

/*

  Evaluate 3d basis functions and derivatives (with 
  respect to the coordinates of the reference element) 
  of nurbs basis functions.

*/

void DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv(
    blitz::Array<double,1>                 & nurbs_shape_funct,
    blitz::Array<double,2>                 & nurbs_shape_deriv,
    const blitz::Array<double,1>           & uv               ,
    const vector<blitz::Array<double, 1> > & knots            ,
    const blitz::Array<double, 1>          & weights          ,
    const DRT::Element::DiscretizationType & distype          )
{

  int degree = 1;

  switch (distype)
  {
  case DRT::Element::nurbs8:
  {
    degree=1;

    break;
  }
  case DRT::Element::nurbs27:
  {
    degree=2;

    break;
  }
  default:
  {
    dserror("Unknown distype for nurbs element evaluation\n");
  }
  }
  
  // size is the number of control points/basis 
  // functions of this element
  const int size   = (degree+1)*(degree+1)*(degree+1);

  // Gausspoint is in [-1:1]x[-1:1]x[-1:1]
  //
  // has to be mapped into knotvector evaluation interval
  
  /*
       reference element                element in knotspan
                                  
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

  */
    
  // the jacobian matrix of Psi is of the very 
  // simple form
  //
  //
  //                                        Ju     Jv
  //                                        |       |
  //                                        |       |		 				                
  //           +   dxi   dxi   dxi -+   +-  du      |    -+  
  //           |   ---   ---   ---  |   |  ---   0  | 0   | 
  //           |    du    dv    dw  |   |  2.0      |     |  
  //           |                    |   |          /      |  
  //           |  deta  deta  deta  |   |        dv       |  
  //    Psi' = |  ----  ----  ----  | = |   0   ---   0   |  
  //           |   du    dv    dw   |   |       2.0       |  
  //           |                    |   |                 |  
  //           |   dnu   dnu   dnu  |   |             dw  |  
  //           |  ----  ----  ----  |   |   0    0   ---  |
  //           +-  du    dv    dw  -+   +-           2.0 -+
  //                                                  |
  //                                                  |
  //                                                  Jw

  double du = (knots[0])(degree+1)-(knots[0])(degree);
  double dv = (knots[1])(degree+1)-(knots[1])(degree);
  double dw = (knots[2])(degree+1)-(knots[2])(degree);

  double Ju =du/2.0;
  double Jv =dv/2.0;
  double Jw =dw/2.0;

  // get mapped coordinates
  
  double xi = (knots[0])(degree)+((uv(0)+1))*Ju;
  double eta= (knots[1])(degree)+((uv(1)+1))*Jv;
  double nu = (knots[2])(degree)+((uv(2)+1))*Jw;
  
  // -----------------------------------------------
  //  PART I: EVALUATION OF  BEZIER SHAPE FUNCTIONS
  // -----------------------------------------------
  
  blitz::Array<double,1> bezier_shape_funct (  size);
  blitz::Array<double,2> bezier_shape_deriv (2,size);
    
  // allocate bspline polynomials for both direction
  DRT::NURBS::UTILS::BsplinePolynomial bspline_xi (degree,knots[0]);
  DRT::NURBS::UTILS::BsplinePolynomial bspline_eta(degree,knots[1]);
  DRT::NURBS::UTILS::BsplinePolynomial bspline_nu (degree,knots[2]);

  // get temporary doubles for derivatives and 
  // values of the above devlared polynomials
  double bspline_xi_value;
  double bspline_eta_value;
  double bspline_nu_value;

  double bspline_xi_derivative;
  double bspline_eta_derivative;
  double bspline_nu_derivative;
  
  // define temporary int variable to compute the 
  // number of the basis function from i,j,k
  int id;
  

  /*
       Bezier basis function 
                |
                |
           B_{i,j,k}(xi,eta,nu) = N_i(xi) *M_j (eta)*L_k (nu)
                                  |        |         |
                                  |        |         | 
                           bspline_u_value |         |
                                           |         |
                                    bspline_v_value  |
                                                     |
                                               bspline_w_value


          dB_{i,j,k}   
          ----------(xi,eta,nu) = N_i'(xi)*M_j (eta)*L_k (nu)
              dxi 

    
          dB_{i,j,k}   
          ----------(xi,eta,nu) = N_i (xi)*M_j'(eta)*L_k (nu)
             deta 

    
          dB_{i,j,k}   
          ----------(xi,eta,nu) = N_i (xi)*M_j (eta)*L_k'(nu)
             dnu 

  */
    
  // loop all basis functions (corresponding to 
  // control points)
  for(int rr=0;rr<degree+1;++rr)
  {
    // in first direction:
    
    // get bsplinevalue and derivative
    bspline_xi.EvaluateBsplineAndDeriv(
      bspline_xi_value,
      bspline_xi_derivative,
      xi,
      rr);
    
    for(int mm=0;mm<degree+1;++mm)
    {
      // in second direction:
      
      // get bsplinevalue and derivative
      bspline_eta.EvaluateBsplineAndDeriv(
	bspline_eta_value,
	bspline_eta_derivative,
	eta,
	mm);

      for(int nn=0;nn<degree+1;++nn)
      {
	// in third direction:
            
	// get bsplinevalue and derivative
	bspline_nu.EvaluateBsplineAndDeriv(
	  bspline_nu_value,
	  bspline_nu_derivative,
	  nu,
	  nn);

	// get the number of the basis function
	id = rr+(degree+1)*(mm+nn*(degree+1));
      
	// set value to bezier_shape_funct
	bezier_shape_funct(id)=bspline_xi_value*bspline_eta_value*bspline_nu_value;
      
	// set values to bezier_shape_deriv
	bezier_shape_deriv(0,id)=bspline_xi_derivative*bspline_eta_value*bspline_nu_value;
	bezier_shape_deriv(1,id)=bspline_xi_value*bspline_eta_derivative*bspline_nu_value;
	bezier_shape_deriv(2,id)=bspline_xi_value*bspline_eta_value*bspline_nu_derivative;
      }
    }
  }

  // -----------------------------------------------
  //  PART II: PROJECTING TO NURBS SHAPE FUNCTIONS
  // -----------------------------------------------
  
  // alloc temporary doubles, initialise to zero
  double                 sum_funct_weight    ;
  blitz::Array<double,1> sum_deriv_weight (3);
  
  sum_funct_weight =0;
  sum_deriv_weight =0;

  /*loop all control points, compute sums
      
                         +----+
                          \
       sum_funct_weight =  +    w_{i,j,k} * B_{i,j,k}(xi,eta,nu)
                          /
                         +----+
      
      
                            +----+
                             \                 dB_{i,j,k}
       sum_deriv_weight[0] =  +    w_{i,j,k} * ----------(xi,eta,nu)
                             /                     dxi
                            +----+
      
                            +----+
                             \                 dB_{i,j,k}
       sum_deriv_weight[1] =  +    w_{i,j,k} * ----------(xi,eta,nu)
                             /                     deta
                            +----+
      
                            +----+
                             \                 dB_{i,j,k}
       sum_deriv_weight[2] =  +    w_{i,j,k} * ----------(xi,eta,nu)
                             /                     dnu
                            +----+

  */

  for(int rr=0;rr<size;++rr)
  {
    sum_funct_weight+=weights(rr)*bezier_shape_funct(rr);
    
    for(int mm=0;mm<3;++mm)
    {
      sum_deriv_weight(mm)+=weights(rr)*bezier_shape_deriv(mm,rr);
    }
  }
  
  // loop all basis functions
  for(int rr=0;rr<size;++rr)
  {    

    /* Compute Nurbs basis functions
      
      
                      w_{i,j,k} * B_{i,j,k}(xi,eta,nu)
                
       N_{i,j,k} = ---------------------------------------
                   +----+				
                    \				
                     +    w_{r,s,t} * B_{r,s,t}(xi,eta,nu)
      		    /				
      		   +----+				
      		    r,s,t                           
    */
    nurbs_shape_funct(rr)=weights(rr)*bezier_shape_funct(rr)/sum_funct_weight;

    /* Nurbs derivatives are defined by the chain rule
    
                                       +- +----+                  -+                       +- +----+                     -+ 
                                       |   \                       |                       |   \	       dB_{r,s,t} |
                            dB_{i,j,k} |    +  w_{r,s,t}*B_{r,s,t} |                       |    +    w_{r,s,t}*---------- |
                  w_{i,j,k}*----------*|   /                       | - w_{i,j,k}*B_{i,j,k}*|   /                   dxi    |
    	  	                dxi    |  +----+                   |                       |  +----+                      |
     dN_{i,j,k}                        +-  r,s,t                  -+                       +-  r,s,t                     -+
     ---------- = ---------------------------------------------------------------------------------------------------------
         dxi                           +- +----+                                 -+ 2
                                       |   \                                      |
                                       |    +    w_{r,s,t} * B_{r,s,t}(xi,eta,nu) |
    		                       |   /                                      |
    		                       |  +----+                                  |
    		                       +-  r,s,t                                 -+
    */
    // loop directions to compute derivatives
    for(int mm=0;mm<3;++mm)
    {
      nurbs_shape_deriv(mm,rr)
	=
	bezier_shape_deriv(mm,rr)*sum_funct_weight
	-
	bezier_shape_funct(rr)*sum_deriv_weight(mm);
      
      nurbs_shape_deriv(mm,rr)*=weights(rr)/(sum_funct_weight*sum_funct_weight);
    }
  }


  // -----------------------------------------------
  //  PART III: TRANSFORMING DERIVATIVES FROM xi/eta
  //            TO u/v
  // -----------------------------------------------
   
  // we already know the jacobian matrix of psi
  //
  //            +   dxi   dxi   dxi -+   +-  du           -+
  //            |   ---   ---   ---  |   |  ---   0    0   |
  //            |    du    dv    dw  |   |  2.0            |
  //            |                    |   |                 |
  //            |  deta  deta  deta  |   |        dv       |
  //     Psi' = |  ----  ----  ----  | = |   0   ---   0   |
  //            |   du    dv    dw   |   |       2.0       |
  //            |                    |   |                 |
  //            |   dnu   dnu   dnu  |   |             dw  |
  //            |  ----  ----  ----  |   |   0    0   ---  |
  //            +-  du    dv    dw  -+   +-           2.0 -+
  //
  //
  // we will obtain the derivatives with respect to  
  // u,v,w just by multiplication  
  
  // loop all basis function derivatives
  for(int rr=0;rr<size;++rr)
  {
    nurbs_shape_deriv(0,rr)*=Ju;
    nurbs_shape_deriv(1,rr)*=Jv;
    nurbs_shape_deriv(2,rr)*=Jw;
  }    

  return;
}


/*

 Evaluate 3d basis functions, first and second derivatives 
 (with respect to the coordinates of the reference 
 element) of nurbs basis functions.

*/

void DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv_deriv2(
    blitz::Array<double,1>                 & nurbs_shape_funct ,
    blitz::Array<double,2>                 & nurbs_shape_deriv ,
    blitz::Array<double,2>                 & nurbs_shape_deriv2,
    const blitz::Array<double,1>           & uv                ,
    const vector<blitz::Array<double, 1> > & knots             ,
    const blitz::Array<double, 1>          & weights           ,
    const DRT::Element::DiscretizationType & distype           )
{

  int degree = 1;

  switch (distype)
  {
  case DRT::Element::nurbs8:
  {
    degree=1;

    break;
  }
  case DRT::Element::nurbs27:
  {
    degree=2;

    break;
  }
  default:
  {
    dserror("Unknown distype for nurbs element evaluation\n");
  }
  }
  
  // size is the number of control points/basis 
  // functions of this element
  const int size   = (degree+1)*(degree+1)*(degree+1);

  // Gausspoint is in [-1:1]x[-1:1]x[-1:1]
  //
  // has to be mapped into knotvector evaluation interval
  
  /*
       reference element                element in knotspan
                                  
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
  */
    
  // the jacobian matrix of Psi is of the very 
  // simple form
  //
  //
  //                                        Ju     Jv
  //                                        |       |
  //                                        |       |		 				                
  //           +   dxi   dxi   dxi -+   +-  du      |    -+  
  //           |   ---   ---   ---  |   |  ---   0  | 0   | 
  //           |    du    dv    dw  |   |  2.0      |     |  
  //           |                    |   |          /      |  
  //           |  deta  deta  deta  |   |        dv       |  
  //    Psi' = |  ----  ----  ----  | = |   0   ---   0   |  
  //           |   du    dv    dw   |   |       2.0       |  
  //           |                    |   |                 |  
  //           |   dnu   dnu   dnu  |   |             dw  |  
  //           |  ----  ----  ----  |   |   0    0   ---  |
  //           +-  du    dv    dw  -+   +-           2.0 -+
  //                                                  |
  //                                                  |
  //                                                  Jw

  double du = (knots[0])(degree+1)-(knots[0])(degree);
  double dv = (knots[1])(degree+1)-(knots[1])(degree);
  double dw = (knots[2])(degree+1)-(knots[2])(degree);

  double Ju =du/2.0;
  double Jv =dv/2.0;
  double Jw =dw/2.0;

  // get mapped coordinates
  
  double xi = (knots[0])(degree)+((uv(0)+1))*Ju;
  double eta= (knots[1])(degree)+((uv(1)+1))*Jv;
  double nu = (knots[2])(degree)+((uv(2)+1))*Jw;
  
  // -----------------------------------------------
  //  PART I: EVALUATION OF  BEZIER SHAPE FUNCTIONS
  // -----------------------------------------------
  
  blitz::Array<double,1> bezier_shape_funct (  size);
  blitz::Array<double,2> bezier_shape_deriv (2,size);
  blitz::Array<double,2> bezier_shape_deriv2(3,size);
    
  // allocate bspline polynomials for both direction
  DRT::NURBS::UTILS::BsplinePolynomial bspline_xi (degree,knots[0]);
  DRT::NURBS::UTILS::BsplinePolynomial bspline_eta(degree,knots[1]);
  DRT::NURBS::UTILS::BsplinePolynomial bspline_nu (degree,knots[2]);

  // get temporary doubles for derivatives and 
  // values of the above devlared polynomials
  double bspline_xi_value;
  double bspline_eta_value;
  double bspline_nu_value;

  double bspline_xi_derivative;
  double bspline_eta_derivative;
  double bspline_nu_derivative;

  double bspline_xi_deriv2;
  double bspline_eta_deriv2;
  double bspline_nu_deriv2;
  
  // define temporary int variable to compute the 
  // number of the basis function from i,j,k
  int id;
  

  /*
       Bezier basis function 
                |
                |
           B_{i,j,k}(xi,eta,nu) = N_i(xi) *M_j (eta)*L_k (nu)
                                  |        |         |
                                  |        |         | 
                           bspline_u_value |         |
                                           |         |
                                    bspline_v_value  |
                                                     |
                                               bspline_w_value


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

  */
    
  // loop all basis functions (corresponding to 
  // control points)
  for(int rr=0;rr<degree+1;++rr)
  {
    // in first direction:
    
    // get bsplinevalue and derivative
    bspline_xi.EvaluateBsplineFirstAndSecondDeriv(
      bspline_xi_value,
      bspline_xi_derivative,
      bspline_xi_deriv2,
      xi,
      rr);
    
    for(int mm=0;mm<degree+1;++mm)
    {
      // in second direction:
      
      // get bsplinevalue and derivative
      bspline_eta.EvaluateBsplineFirstAndSecondDeriv(
	bspline_eta_value,
	bspline_eta_derivative,
	bspline_eta_deriv2,
	eta,
	mm);

      for(int nn=0;nn<degree+1;++nn)
      {
	// in third direction:
            
	// get bsplinevalue and derivative
	bspline_nu.EvaluateBsplineFirstAndSecondDeriv(
	  bspline_nu_value,
	  bspline_nu_derivative,
	  bspline_nu_deriv2,
	  nu,
	  nn);

	// get the number of the basis function
	id = rr+(degree+1)*(mm+nn*(degree+1));
      
	// set value to bezier_shape_funct
	bezier_shape_funct(id)=bspline_xi_value*bspline_eta_value*bspline_nu_value;
      
	// set values to bezier_shape_deriv
	bezier_shape_deriv(0,id)=bspline_xi_derivative*bspline_eta_value*bspline_nu_value;
	bezier_shape_deriv(1,id)=bspline_xi_value*bspline_eta_derivative*bspline_nu_value;
	bezier_shape_deriv(2,id)=bspline_xi_value*bspline_eta_value*bspline_nu_derivative;

	// set values to bezier_shape_deriv2
	bezier_shape_deriv2(0,id)=bspline_xi_deriv2    *bspline_eta_value     *bspline_nu_value     ;
	bezier_shape_deriv2(1,id)=bspline_xi_value     *bspline_eta_deriv2    *bspline_nu_value     ;
	bezier_shape_deriv2(2,id)=bspline_xi_value     *bspline_eta_value     *bspline_nu_deriv2    ;
	bezier_shape_deriv2(3,id)=bspline_xi_derivative*bspline_eta_derivative*bspline_nu_value     ;
	bezier_shape_deriv2(4,id)=bspline_xi_derivative*bspline_eta_value     *bspline_nu_derivative;
	bezier_shape_deriv2(5,id)=bspline_xi_value     *bspline_eta_derivative*bspline_nu_derivative;

      }
    }
  }

  // -----------------------------------------------
  //  PART II: PROJECTING TO NURBS SHAPE FUNCTIONS
  // -----------------------------------------------
  
  // alloc temporary doubles, initialise to zero
  double                 sum_funct_weight    ;
  blitz::Array<double,1> sum_deriv_weight (3);
  blitz::Array<double,1> sum_deriv2_weight(6);
  
  sum_funct_weight =0;
  sum_deriv_weight =0;
  sum_deriv2_weight=0;

  /*loop all control points, compute sums
      
                         +----+
                          \
       sum_funct_weight =  +    w_{i,j,k} * B_{i,j,k}(xi,eta,nu)
                          /
                         +----+
      
      
                            +----+
                             \                 dB_{i,j,k}
       sum_deriv_weight[0] =  +    w_{i,j,k} * ----------(xi,eta,nu)
                             /                     dxi
                            +----+
      
                            +----+
                             \                 dB_{i,j,k}
       sum_deriv_weight[1] =  +    w_{i,j,k} * ----------(xi,eta,nu)
                             /                     deta
                            +----+
      
                            +----+
                             \                 dB_{i,j,k}
       sum_deriv_weight[2] =  +    w_{i,j,k} * ----------(xi,eta,nu)
                             /                     dnu
                            +----+


      
                             +----+              2
                              \                 d B_{i,j,k}
       sum_deriv2_weight[0] =  +    w_{i,j,k} * -----------(xi,eta,nu)
                              /                        2
                             +----+                 dxi
      
                             +----+              2 
                              \                 d B_{i,j,k}
       sum_deriv2_weight[1] =  +    w_{i,j,k} * -----------(xi,eta,nu)
                              /                        2
                             +----+                deta
      
                             +----+              2
                              \                 d B_{i,j,k}
       sum_deriv2_weight[2] =  +    w_{i,j,k} * -----------(xi,eta,nu)
                              /                       2
                             +----+                dnu

      
                             +----+              2
                              \                 d B_{i,j,k}
       sum_deriv2_weight[3] =  +    w_{i,j,k} * -----------(xi,eta,nu)
                              /                  dxi deta   
                             +----+              
      
                             +----+              2 
                              \                 d B_{i,j,k}
       sum_deriv2_weight[4] =  +    w_{i,j,k} * -----------(xi,eta,nu)
                              /                   dxi dnu
                             +----+             
      
                             +----+              2
                              \                 d B_{i,j,k}
       sum_deriv2_weight[5] =  +    w_{i,j,k} * -----------(xi,eta,nu)
                              /                   deta dnu
                             +----+
  */

  for(int rr=0;rr<size;++rr)
  {
    sum_funct_weight+=weights(rr)*bezier_shape_funct(rr);
    
    for(int mm=0;mm<3;++mm)
    {
      sum_deriv_weight(mm)+=weights(rr)*bezier_shape_deriv(mm,rr);
    }


    for(int mm=0;mm<6;++mm)
    {
      sum_deriv2_weight(mm)+=weights(rr)*bezier_shape_deriv2(mm,rr);
    }
  }
  
  // loop all basis functions
  for(int rr=0;rr<size;++rr)
  {    

    /* Compute Nurbs basis functions
      
      
                      w_{i,j,k} * B_{i,j,k}(xi,eta,nu)
                
       N_{i,j,k} = ---------------------------------------
                   +----+				
                    \				
                     +    w_{r,s,t} * B_{r,s,t}(xi,eta,nu)
      		    /				
      		   +----+				
      		    r,s,t                           
    */
    nurbs_shape_funct(rr)=weights(rr)*bezier_shape_funct(rr)/sum_funct_weight;

    /* Nurbs derivatives are defined by the chain rule
    
                                       +- +----+                  -+                       +- +----+                     -+ 
                                       |   \                       |                       |   \	       dB_{r,s,t} |
                            dB_{i,j,k} |    +  w_{r,s,t}*B_{r,s,t} |                       |    +    w_{r,s,t}*---------- |
                  w_{i,j,k}*----------*|   /                       | - w_{i,j,k}*B_{i,j,k}*|   /                   dxi    |
    	  	                dxi    |  +----+                   |                       |  +----+                      |
     dN_{i,j,k}                        +-  r,s,t                  -+                       +-  r,s,t                     -+
     ---------- = ---------------------------------------------------------------------------------------------------------
         dxi                           +- +----+                                 -+ 2
                                       |   \                                      |
                                       |    +    w_{r,s,t} * B_{r,s,t}(xi,eta,nu) |
    		                       |   /                                      |
    		                       |  +----+                                  |
    		                       +-  r,s,t                                 -+
    */
    // loop directions to compute derivatives
    for(int mm=0;mm<3;++mm)
    {
      nurbs_shape_deriv(mm,rr)
	=
	bezier_shape_deriv(mm,rr)*sum_funct_weight
	-
	bezier_shape_funct(rr)*sum_deriv_weight(mm);
      
      nurbs_shape_deriv(mm,rr)*=weights(rr)/(sum_funct_weight*sum_funct_weight);
    }



    /* Nurbs second derivatives are calculated by a 
       second application of the chain rule
    
                                                                  +-
                                                                  |
                                                                  |
                                                                  |
    	  	                     w_{i,j,k}                    |   2
     d N_{i,j,k}                                                  |  d B_{i,j,k}
     ----------- = -------------------------------------------- * |  ----------- -
      dxi  dxi     +- +----+                                 -+   |   dxi  dxi
         n    m    |   \                                      |   |      n    m
                   |    +    w_{r,s,t} * B_{r,s,t}(xi,eta,nu) |   |
                   |   /                                      |   |
                   |  +----+                                  |   |
                   +-  r,s,t                                 -+   +-
  
                     +- +----+                                  -+
          dB_{i,j,k} |   \                 dB_{r,s,t}            |
          ---------- |    +    w_{r,s,t} * ----------(xi,eta,nu) |
             dxi     |   /                    dxi                |
                n    |  +----+                   m               |
                     +-  r,s,t                                  -+
        - -------------------------------------------------------- -
               +- +----+                                 -+            
               |   \                                      |             
               |    +    w_{r,s,t} * B_{r,s,t}(xi,eta,nu) |   
               |   /                                      |   
               |  +----+                                  |   
               +-  r,s,t                                 -+   


  
                     +- +----+                                  -+
          dB_{i,j,k} |   \                 dB_{r,s,t}            |
          ---------- |    +    w_{r,s,t} * ----------(xi,eta,nu) |
             dxi     |   /                    dxi                |
                m    |  +----+                   n               |
                     +-  r,s,t                                  -+
        - -------------------------------------------------------- -
               +- +----+                                 -+            
               |   \                                      |             
               |    +    w_{r,s,t} * B_{r,s,t}(xi,eta,nu) |   
               |   /                                      |   
               |  +----+                                  |   
               +-  r,s,t                                 -+   


                       +- +----+              2                    -+
                       |   \                 d B_{r,s,t}            |
                       |    +    w_{r,s,t} * -----------(xi,eta,nu) |
           B_{i,j,k} * |   /                  dxi  dxi              |
                       |  +----+                 n    m             |
                       +-  r,s,t                                   -+
        - ----------------------------------------------------------- +
                 +- +----+                                 -+            
                 |   \                                      |             
                 |    +    w_{r,s,t} * B_{r,s,t}(xi,eta,nu) |   
                 |   /                                      |   
                 |  +----+                                  |   
                 +-  r,s,t                                 -+   


                           +- +----+                                  -+   +- +----+                                  -+  -+
                           |   \                 dB_{r,s,t}            |   |   \                 dB_{r,s,t}            |   |
                           |    +    w_{r,s,t} * ----------(xi,eta,nu) |   |    +    w_{r,s,t} * ----------(xi,eta,nu) |   |
           2 * B_{i,j,k} * |   /                    dxi                | * |   /                    dxi                |   |
                           |  +----+                   n               |   |  +----+                   m               |   |
                           +-  r,s,t                                  -+   +-  r,s,t                                  -+   |
        + --------------------------------------------------------------------------------------------------------------   |
                                       +- +----+                                 -+ 2                                      |
                                       |   \                                      |                                        |
                                       |    +    w_{r,s,t} * B_{r,s,t}(xi,eta,nu) |                                        |
                                       |   /                                      |                                        |
                                       |  +----+                                  |                                        |
                                       +-  r,s,t                                 -+                                       -+
    */

    // we apply a Horner-type scheme for the 
    // multiplication with sum_funct_weight
    nurbs_shape_deriv2(0,rr) =2*bezier_shape_funct(rr)*sum_deriv_weight(0)*sum_deriv_weight(0);
    nurbs_shape_deriv2(0,rr)/=sum_funct_weight;
    nurbs_shape_deriv2(0,rr)-=bezier_shape_funct(rr)  *sum_deriv2_weight(0);
    nurbs_shape_deriv2(0,rr)-=bezier_shape_deriv(0,rr)*sum_deriv_weight (0);
    nurbs_shape_deriv2(0,rr)-=bezier_shape_deriv(0,rr)*sum_deriv_weight (0);
    nurbs_shape_deriv2(0,rr)/=sum_funct_weight;
    nurbs_shape_deriv2(0,rr)+=bezier_shape_deriv2(0,rr);
    nurbs_shape_deriv2(0,rr)*=weights(rr)/sum_funct_weight;

    nurbs_shape_deriv2(1,rr) =2*bezier_shape_funct(rr)*sum_deriv_weight(1)*sum_deriv_weight(1);
    nurbs_shape_deriv2(1,rr)/=sum_funct_weight;
    nurbs_shape_deriv2(1,rr)-=bezier_shape_funct(rr)  *sum_deriv2_weight(1);
    nurbs_shape_deriv2(1,rr)-=bezier_shape_deriv(1,rr)*sum_deriv_weight (1);
    nurbs_shape_deriv2(1,rr)-=bezier_shape_deriv(1,rr)*sum_deriv_weight (1);
    nurbs_shape_deriv2(1,rr)/=sum_funct_weight;
    nurbs_shape_deriv2(1,rr)+=bezier_shape_deriv2(1,rr);
    nurbs_shape_deriv2(1,rr)*=weights(rr)/sum_funct_weight;

    nurbs_shape_deriv2(2,rr) =2*bezier_shape_funct(rr)*sum_deriv_weight(2)*sum_deriv_weight(2);
    nurbs_shape_deriv2(2,rr)/=sum_funct_weight;
    nurbs_shape_deriv2(2,rr)-=bezier_shape_funct(rr)  *sum_deriv2_weight(2);
    nurbs_shape_deriv2(2,rr)-=bezier_shape_deriv(2,rr)*sum_deriv_weight (2);
    nurbs_shape_deriv2(2,rr)-=bezier_shape_deriv(2,rr)*sum_deriv_weight (2);
    nurbs_shape_deriv2(2,rr)/=sum_funct_weight;
    nurbs_shape_deriv2(2,rr)+=bezier_shape_deriv2(2,rr);
    nurbs_shape_deriv2(2,rr)*=weights(rr)/sum_funct_weight;

    nurbs_shape_deriv2(3,rr) =2*bezier_shape_funct(rr)*sum_deriv_weight(0)*sum_deriv_weight(1);
    nurbs_shape_deriv2(3,rr)/=sum_funct_weight;
    nurbs_shape_deriv2(3,rr)-=bezier_shape_funct(rr)  *sum_deriv2_weight(3);
    nurbs_shape_deriv2(3,rr)-=bezier_shape_deriv(0,rr)*sum_deriv_weight (1);
    nurbs_shape_deriv2(3,rr)-=bezier_shape_deriv(1,rr)*sum_deriv_weight (0);
    nurbs_shape_deriv2(3,rr)/=sum_funct_weight;
    nurbs_shape_deriv2(3,rr)+=bezier_shape_deriv2(3,rr);
    nurbs_shape_deriv2(3,rr)*=weights(rr)/sum_funct_weight;

    nurbs_shape_deriv2(4,rr) =2*bezier_shape_funct(rr)*sum_deriv_weight(0)*sum_deriv_weight(2);
    nurbs_shape_deriv2(4,rr)/=sum_funct_weight;
    nurbs_shape_deriv2(4,rr)-=bezier_shape_funct(rr)  *sum_deriv2_weight(4);
    nurbs_shape_deriv2(4,rr)-=bezier_shape_deriv(0,rr)*sum_deriv_weight (2);
    nurbs_shape_deriv2(4,rr)-=bezier_shape_deriv(2,rr)*sum_deriv_weight (0);
    nurbs_shape_deriv2(4,rr)/=sum_funct_weight;
    nurbs_shape_deriv2(4,rr)+=bezier_shape_deriv2(4,rr);
    nurbs_shape_deriv2(4,rr)*=weights(rr)/sum_funct_weight;

    nurbs_shape_deriv2(5,rr) =2*bezier_shape_funct(rr)*sum_deriv_weight(1)*sum_deriv_weight(2);
    nurbs_shape_deriv2(5,rr)/=sum_funct_weight;
    nurbs_shape_deriv2(5,rr)-=bezier_shape_funct(rr)  *sum_deriv2_weight(5);
    nurbs_shape_deriv2(5,rr)-=bezier_shape_deriv(1,rr)*sum_deriv_weight (2);
    nurbs_shape_deriv2(5,rr)-=bezier_shape_deriv(2,rr)*sum_deriv_weight (1);
    nurbs_shape_deriv2(5,rr)/=sum_funct_weight;
    nurbs_shape_deriv2(5,rr)+=bezier_shape_deriv2(5,rr);
    nurbs_shape_deriv2(5,rr)*=weights(rr)/sum_funct_weight;
  }


  // -----------------------------------------------
  //  PART III: TRANSFORMING DERIVATIVES FROM xi/eta
  //            TO u/v
  // -----------------------------------------------
   
  // we already know the jacobian matrix of psi
  //
  //            +   dxi   dxi   dxi -+   +-  du           -+
  //            |   ---   ---   ---  |   |  ---   0    0   |
  //            |    du    dv    dw  |   |  2.0            |
  //            |                    |   |                 |
  //            |  deta  deta  deta  |   |        dv       |
  //     Psi' = |  ----  ----  ----  | = |   0   ---   0   |
  //            |   du    dv    dw   |   |       2.0       |
  //            |                    |   |                 |
  //            |   dnu   dnu   dnu  |   |             dw  |
  //            |  ----  ----  ----  |   |   0    0   ---  |
  //            +-  du    dv    dw  -+   +-           2.0 -+
  //
  //
  // we will obtain the derivatives with respect to  
  // u,v,w just by multiplication  
  
  // loop all basis function derivatives
  for(int rr=0;rr<size;++rr)
  {
    nurbs_shape_deriv(0,rr)*=Ju;
    nurbs_shape_deriv(1,rr)*=Jv;
    nurbs_shape_deriv(2,rr)*=Jw;

    nurbs_shape_deriv2(0,rr)*=Ju*Ju;
    nurbs_shape_deriv2(1,rr)*=Jv*Jv;
    nurbs_shape_deriv2(2,rr)*=Jw*Jw;
    nurbs_shape_deriv2(3,rr)*=Ju*Jv;
    nurbs_shape_deriv2(4,rr)*=Ju*Jw;
    nurbs_shape_deriv2(5,rr)*=Jv*Jw;
  }    

  return;
}


#endif // CCADISCRET

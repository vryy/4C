/*!----------------------------------------------------------------------
\file
\brief evaluate galerkin part of stiffness matrix

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
/*!---------------------------------------------------------------------                                         
\brief evaluate galerkin part of Kkapeps

<pre>                                                         he 12/02

In this routine the galerkin part of matrix Kkapeps is calculated:

    /
   |  u * grad(kapeps) * psi                    d_omega
  /

    /
   |  (nue+nue_t/sig) * grad(kapeps) * grad(psi)    d_omega
  /


    /
   |  factor * kapeps_old * kapeps * psi     d_omega
  /

LOW-REYNOLD's MODEL only for kappa:

    /
   | 2.0 * visc * [ 2*grad(k_old)/(4*k_old)  * grad (k) - grad(k_old)*grad(k_old)/(4*k_old^2) * k ] * psi     d_omega
  /
  


NOTE: there's only one elestif
      --> Kkapeps is stored in estif[0..(iel-1)][0..(iel-1)]
      
</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **estif        double	   (i/o)  ele stiffness matrix
\param   kapepsint    double	   (i)    kapeps at int point
\param  *velint       double	   (i)    velocity at int point
\param   eddyint      double	   (i)    eddy-viscosity at int point
\param  *kapepsderxy  double	   (i)    kapeps deriv. at int point
\param  *funct        double	   (i)    nat. shape funcs
\param **derxy        double	   (i)    global coord. deriv.
\param   fac 	    double	   (i)    weighting factor	      
\param   visc         double	   (i)    fluid viscosity	     
\param   factor        double	   (i)    facor 
\param   sig          double	   (i)    facor 
\param   iel	    int  	   (i)    number of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calkkapeps(
                FLUID_DYN_CALC  *dynvar,
		    double         **estif,   
		    double           kapepsint, 
	          double          *velint,
                double           eddyint, 
                double          *kapepsderxy, 
                double          *funct,  
		    double         **derxy,  
		    double           fac,    
		    double           visc,   
		    double           factor,
		    double           sig,
                int              iel     
              )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   icn  - column node: number of node considered for matrix column    |  
/*----------------------------------------------------------------------*/
int     irow, icol,irn,icn;
double  auxc;

#ifdef DEBUG 
dstrc_enter("f2_calkkapeps");
#endif		

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix K:
    /
   |  u * grad(kapeps) * psi   d_omega
  /
 *----------------------------------------------------------------------*/
icol=0;
for (icn=0;icn<iel;icn++) 
{
   irow=0;
   auxc = (derxy[0][icn]*velint[0]+derxy[1][icn]*velint[1])*fac;

   for (irn=0;irn<iel;irn++)
   {
      estif[irow][icol] += auxc*funct[irn];
      irow += 1;
   } /* end loop over irn */
   icol += 1;
} /* end loop over icn */

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix K:
    /
   |  (nue+nue_t/sig) * grad(kapeps) * grad(psi)   d_omega
  /
 *----------------------------------------------------------------------*/
icol=0;
 for (icn=0;icn<iel;icn++)
 {
      irow=0;  
      auxc = (visc+eddyint/sig)*fac;   
   
      for (irn=0;irn<iel;irn++)
      {
	   estif[irow][icol] += auxc*(derxy[0][icn]*derxy[0][irn]+derxy[1][icn]*derxy[1][irn]);
         irow += 1;	 
      } /* end loop over irn */
      icol += 1;
 } /* end loop over icn */

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix K:
    /
   |  factor * kapeps_old * kapeps * psi   d_omega
  /
 *----------------------------------------------------------------------*/
icol=0;
   for (icn=0;icn<iel;icn++)
   {
      irow=0;  
      auxc=factor*kapepsint*funct[icn]*fac;

      for (irn=0;irn<iel;irn++)
      {
        estif[irow][icol]   += auxc*funct[irn];
        irow += 1;
      } /* end loop over irn */
      icol += 1;
   } /* end loop over icn */


if(dynvar->kapeps_flag==0) 
{
/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix K for LOW-REYNOLD's MODEL:
    /
   | 2.0 * visc * [ 2*grad(k_old)/(4*k_old)  * grad (k) - grad(k_old)*grad(k_old)/(4*k_old^2) * k ]   * psi  d_omega
  /
*----------------------------------------------------------------------*/
icol=0;
   for (icn=0;icn<iel;icn++)
   {
      irow =0;  

      auxc  = 2/(4*kapepsint)*
              (kapepsderxy[0]*derxy[0][icn]+kapepsderxy[1]*derxy[1][icn]);
      auxc -= (pow(kapepsderxy[0],2)+pow(kapepsderxy[1],2))/
              (4*kapepsint*kapepsint) * funct[icn];
              
      for (irn=0;irn<iel;irn++)
      {
        estif[irow][icol]  += 2.0*visc*fac*auxc*funct[irn];
	  
        irow += 1;
      } /* end loop over irn */
      icol += 1;
   } /* end loop over icn */

} /* endif */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calkvv */


/*!---------------------------------------------------------------------                                         
\brief evaluate galerkin part of Mkapeps

<pre>                                                         he  12/02

In this routine the galerkin part of matrix Mvv is calculated:

    /
   |  psi * kapeps    d_omega
  /
  
NOTE: there's only one elemass  			      
      --> Mkapeps is stored in emass[0..(iel-1)][0..(iel-1)]  
      
</pre>
\param **emass     double	   (i/o)  ele mass matrix
\param  *funct     double	   (i)    nat. shape funcs
\param   fac       double	   (i)    weighting factor
\param   iel       int   	   (i)    number of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calmkapeps(
		double         **emass,  
		double          *funct, 
		double           fac,   
		int              iel    
              )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   icn  - column node: number of node considered for matrix column    |  
/*----------------------------------------------------------------------*/
int     irow, icol,irn,icn;  
double  auxc;

#ifdef DEBUG 
dstrc_enter("f2_calmkapeps");
#endif		


/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Mkapeps:
    /
   |  psi * kapeps    d_omega
  /
 *----------------------------------------------------------------------*/
icol=0;
for(icn=0;icn<iel;icn++)
{
   irow = 0;
   auxc = funct[icn]*fac;
   
   for(irn=0;irn<iel;irn++)
   {
      emass[irow][icol] += auxc*funct[irn];
      irow += 1;
   } /* end loop over irn */
 
  icol += 1;

} /* end loop over icn */

/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calmkapeps */


#endif

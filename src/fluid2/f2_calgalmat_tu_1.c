/*!----------------------------------------------------------------------
\file
\brief evaluate galerkin part of stiffness matrix

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
/*!---------------------------------------------------------------------                                         
\brief evaluate galerkin part of Kkapome

<pre>                                                         he 02/03

In this routine the galerkin part of matrix Kkapeps is calculated:

    /
   |  u * grad(kapome) * psi                    d_omega
  /

    /
   |  (nue+nue_t*sig) * grad(kapome) * grad(psi)    d_omega
  /


    /
   |  factor * kapome_old * kapome * psi     d_omega
  /
  


NOTE: there's only one elestif
      --> Kkapome is stored in estif[0..(iel-1)][0..(iel-1)]
      
</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **estif        double	   (i/o)  ele stiffness matrix
\param   kapomeint    double	   (i)    kapome at int point
\param  *velint       double	   (i)    velocity at int point
\param   eddyint      double	   (i)    eddy-viscosity at int point
\param  *funct        double	   (i)    nat. shape funcs
\param **derxy        double	   (i)    global coord. deriv.
\param   fac 	    double	   (i)    weighting factor	      
\param   visc         double	   (i)    fluid viscosity	     
\param   factor       double	   (i)    factor 
\param   sig          double	   (i)    factor 
\param   iel	    int  	   (i)    number of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calkkapome(
                FLUID_DYN_CALC  *dynvar,
		    double         **estif,   
		    double           kapomeint, 
	          double          *velint,
                double           eddyint, 
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
dstrc_enter("f2_calkkapome");
#endif		

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix K:
    /
   |  u * grad(kapome) * psi   d_omega
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
   |  (nue+nue_t*sig) * grad(kapome) * grad(psi)   d_omega
  /
 *----------------------------------------------------------------------*/
icol=0;
 for (icn=0;icn<iel;icn++)
 {
      irow=0;  
      auxc = (visc+eddyint*sig)*fac;   
   
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
   |  factor * kapome_old * kapome * psi   d_omega
  /
 *----------------------------------------------------------------------*/
icol=0;
   for (icn=0;icn<iel;icn++)
   {
      irow=0;  
      auxc=factor*kapomeint*funct[icn]*fac;

      for (irn=0;irn<iel;irn++)
      {
        estif[irow][icol]   += auxc*funct[irn];
        irow += 1;
      } /* end loop over irn */
      icol += 1;
   } /* end loop over icn */


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calkvv */


/*!---------------------------------------------------------------------                                         
\brief evaluate galerkin part of Mkapome

<pre>                                                         he  02/03

In this routine the galerkin part of matrix Mvv is calculated:

    /
   |  psi * kapome    d_omega
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
void f2_calmkapome(
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
dstrc_enter("f2_calmkapome");
#endif		


/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Mkapome:
    /
   |  psi * kapome    d_omega
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
} /* end of f2_calmkapome */


#endif

/*!----------------------------------------------------------------------
\file
\brief evaluate galerkin part of stiffness matrix

<pre>
Maintainer: Thomas Hettich
            hettich@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hettich/
            0771 - 685-6575
</pre>

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
\param **estif        DOUBLE	   (i/o)  ele stiffness matrix
\param   kapomeint    DOUBLE	   (i)    kapome at INT point
\param  *velint       DOUBLE	   (i)    velocity at INT point
\param   eddyint      DOUBLE	   (i)    eddy-viscosity at INT point
\param  *funct        DOUBLE	   (i)    nat. shape funcs
\param **derxy        DOUBLE	   (i)    global coord. deriv.
\param   fac 	    DOUBLE	   (i)    weighting factor	      
\param   visc         DOUBLE	   (i)    fluid viscosity	     
\param   factor       DOUBLE	   (i)    factor 
\param   sig          DOUBLE	   (i)    factor 
\param   iel	    INT  	   (i)    number of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calkkapome(
                FLUID_DYN_CALC  *dynvar,
		    DOUBLE         **estif,   
		    DOUBLE           kapomeint, 
	          DOUBLE          *velint,
                DOUBLE           eddyint, 
                DOUBLE          *funct,  
		    DOUBLE         **derxy,  
		    DOUBLE           fac,    
		    DOUBLE           visc,   
		    DOUBLE           factor,
		    DOUBLE           sig,
                INT              iel     
              )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   icn  - column node: number of node considered for matrix column    |  
/*----------------------------------------------------------------------*/
INT     irow, icol,irn,icn;
DOUBLE  auxc;

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
\param **emass     DOUBLE	   (i/o)  ele mass matrix
\param  *funct     DOUBLE	   (i)    nat. shape funcs
\param   fac       DOUBLE	   (i)    weighting factor
\param   iel       INT   	   (i)    number of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calmkapome(
		DOUBLE         **emass,  
		DOUBLE          *funct, 
		DOUBLE           fac,   
		INT              iel    
              )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   icn  - column node: number of node considered for matrix column    |  
/*----------------------------------------------------------------------*/
INT     irow, icol,irn,icn;  
DOUBLE  auxc;

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

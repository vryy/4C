/*!----------------------------------------------------------------------
\file
\brief stabilisation part of element stiffness matrix for fluid2

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
#include "fluid2_tu.h"
/*!---------------------------------------------------------------------
\brief evaluate stabilisaton part of Kkapome

<pre>                                                         he  02/03

In this routine the stabilisation part of matrix Kvv is calculated:

    /
   |  tau_tu * u * grad(kapome) * grad(psi) * u   d_omega + D.C. 
  /

    /
   |  -tau_tu * div[(nue+nue_t*sig)*grad(kapome)] * grad(psi) * u   d_omega + D.C. 
  /
  
    /
   |  tau_tu * factor * kapome_old * kapome * grad(psi) * u    d_omega + D.C. 
  /  



NOTE: there's only one elestif
      --> Kkapome is stored in estif[0..(iel-1)][0..(iel-1)]
      
</pre>
\param  *ele	     ELEMENT	   (i)	   actual element
\param  *elev	     ELEMENT	   (i)	   actual element for vel.
\param  *dynvar        FLUID_DYN_CALC  (i)
\param **estif          double	   (i/o)   ele stiffness matrix
\param   kapomeint      double	   (i)     kapome at integr. point
\param  *velint         double	   (i)     vel. at integr. point
\param  *velint_dc      double	   (i)     vel. at integr. point for D.C.
\param   eddyint        double	   (i)     eddy at integr. point
\param  *funct          double	   (i)     nat. shape functions
\param **derxy          double	   (i)     global derivatives
\param **derxy2         double	   (i)     2nd global derivatives
\param   fac	      double	   (i)     weighting factor	   
\param   visc	      double	   (i)     fluid viscosity
\param   factor	      double	   (i)     factor
\param   sig	      double	   (i)     factor
\param   iel	      int		   (i)     num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabkkapome(			      
                ELEMENT         *ele,    
		    ELEMENT         *elev,    
                FLUID_DYN_CALC  *dynvar,
		    double         **estif,  
		    double           kapomeint, 
		    double          *velint, 
		    double          *velint_dc, 
                double           eddyint, 
                double          *funct,  
		    double         **derxy,  
		    double         **derxy2, 
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
 |   ird  - row dim.: number of spatial dimension at row node           |  
/*----------------------------------------------------------------------*/
int    irow,icol,irn,icn,ird;
double taumu,taumu_dc;
double c,c_dc;
double auxc,aux;

#ifdef DEBUG 
dstrc_enter("f2_calstabkkapome");
#endif

/*---------------------------------------- set stabilisation parameter */
taumu     = dynvar->tau_tu;
taumu_dc  = dynvar->tau_tu_dc;
c    = fac*taumu;
c_dc = fac*taumu_dc;
/*----------------------------------------------------------------------*
    /
   | tau_tu * u * grad(kapome) * grad(psi) * u  d_omega
  /
 *----------------------------------------------------------------------*/
 icol=0;
   for(icn=0;icn<iel;icn++)
   {
     irow=0;
     auxc = velint[0]*derxy[0][icn]+velint[1]*derxy[1][icn];

      for(irn=0;irn<iel;irn++)
      {
         estif[irow][icol] += auxc*c   *(velint[0]   *derxy[0][irn]+velint[1]   *derxy[1][irn]);
         estif[irow][icol] += auxc*c_dc*(velint_dc[0]*derxy[0][irn]+velint_dc[1]*derxy[1][irn]); 
         irow += 1;
      } /* end of loop over irn */
      icol += 1;
   } /* end of loop over icn */

/*----------------------------------------------------------------------*
    /
   |  -tau_tu * div[(nue+nue_t*sig)*grad(kapome)] * grad(psi) * u   d_omega =
  /
 
    /
   |  -tau_tu *  (nue+nue_t*sig) *  div grad(kapome) * grad(psi) * u   d_omega
  /                                    
 
 *----------------------------------------------------------------------*/
   icol=0;
      for (icn=0;icn<iel;icn++)
      {
 	 irow = 0;
       auxc = (visc+eddyint*sig)*(derxy2[0][icn]+derxy2[1][icn]);

	 for (irn=0;irn<iel;irn++)
	 {          
          estif[irow][icol] -= auxc*c   *(derxy[0][irn]*velint[0]   +derxy[1][irn]*velint[1]);
          estif[irow][icol] -= auxc*c_dc*(derxy[0][irn]*velint_dc[0]+derxy[1][irn]*velint_dc[1]);  
	    irow += 1;
	 } /* end of loop over irn */
	 icol += 1;
      } /* end of loop over icn */

/*----------------------------------------------------------------------*
    /
   |  tau_tu * factor * kapome_old * kapome * grad(psi) * u   d_omega
  /
 *----------------------------------------------------------------------*/
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
	 irow = 0;
       auxc = factor*kapomeint*funct[icn];

	 for (irn=0;irn<iel;irn++)
	 {
          estif[irow][icol]     += auxc*c   *(velint[0]   *derxy[0][irn] + velint[1]   *derxy[1][irn]);
          estif[irow][icol]     += auxc*c_dc*(velint_dc[0]*derxy[0][irn] + velint_dc[1]*derxy[1][irn]);  
          irow += 1;
	 } /* end of loop over irn */
	 icol += 1;
      } /* end of loop over icn */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calstabkvv */


/*!--------------------------------------------------------------------- 
\brief evaluate stabilisaton part of Mkapome

<pre>                                                         he  02/03

In this routine the stabilisation part of matrix Mkapome is calculated:

    /
   |   tau_tu * grad(psi) * u * kapome  d_omega + D.C.
  /
  
  
NOTE: there's only one elestif  			    
      --> Mkapeps is stored in emass[0..(iel-1)][0..(iel-1)]
      
</pre>
\param  *ele	 ELEMENT	     (i)	   actual element
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **emass     double	   (i/o)   ele mass matrix
\param  *velint    double	   (i)     vel. at integr. point
\param  *velint_dc double	   (i)     vel. at integr. point for D.C.
\param  *funct     double	   (i)     nat. shape functions
\param **derxy     double	   (i)     global derivatives
\param   fac	   double	   (i)     weighting factor
\param   iel	   int  	   (i)	   num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabmkapome(
                    ELEMENT         *ele,     
		        FLUID_DYN_CALC  *dynvar,
		        double         **emass,  
    		        double          *velint, 
    		        double          *velint_dc, 
                    double          *funct,  
		        double         **derxy,  
		        double           fac,    
		        int              iel    
                     )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |   
/*----------------------------------------------------------------------*/
int    irow,icol,irn,icn;
double taumu,taumu_dc;
double c,c_dc;
double auxc;

#ifdef DEBUG 
dstrc_enter("f2_calstabmkapome");
#endif

/*---------------------------------------- set stabilisation parameter */
taumu      = dynvar->tau_tu;
taumu_dc   = dynvar->tau_tu_dc;

c    = fac * taumu;
c_dc = fac * taumu_dc;

/*----------------------------------------------------------------------*
   Calculate convection stabilisation part:
    /
   |   tau_tu * grad(psi) * u * kapome  d_omega
  /
 *----------------------------------------------------------------------*/
   icol=0;
   for (icn=0;icn<iel;icn++)
   {
      irow = 0;
      auxc = funct[icn];

      for (irn=0;irn<iel;irn++)
      {
	 emass[irow][icol] += auxc*c   *(velint[0]   *derxy[0][irn]+velint[1]   *derxy[1][irn]);
 	 emass[irow][icol] += auxc*c_dc*(velint_dc[0]*derxy[0][irn]+velint_dc[1]*derxy[1][irn]); 
	 irow += 1;
      } /* end loop over irn */
      icol += 1;      
   } /* end loop over icn */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calstabmvv */

#endif

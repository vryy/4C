/*!----------------------------------------------------------------------
\file
\brief evaluate galerkin part of stiffness matrix

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
/*!---------------------------------------------------------------------                                         
\brief evaluate galerkin part of Kvv

<pre>                                                         genk 04/02

In this routine the galerkin part of matrix Kvv is calculated:

    /
   |  2 * nue * eps(v) : eps(u)   d_omega
  /

    /
   |  v * u_old * grad(u)     d_omega
  /

    /
   |  v * u * grad(u_old)     d_omega
  /

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif
      --> Kvv is stored in estif[0..(2*iel-1)][0..(2*iel-1)]
      
</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **estif     double	   (i/o)  ele stiffness matrix
\param  *velint    double	   (i)    vel at int point
\param **vderxy    double	   (i)    global vel derivatives
\param  *funct     double	   (i)    nat. shape funcs
\param **derxy     double	   (i)    global coord. deriv.
\param   fac 	   double	   (i)    weighting factor	      
\param   visc      double	   (i)    fluid viscosity	     
\param   iel	   int  	   (i)	  number of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calkvv(
                FLUID_DYN_CALC  *dynvar,
		double         **estif,   
		double          *velint, 
		double         **vderxy, 
		double          *funct,  
		double         **derxy,  
		double           fac,    
		double           visc,   
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
double  c,aux;

#ifdef DEBUG 
dstrc_enter("f2_calkvv");
#endif		

c=fac*visc;

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix K:
    /
   |  2 * nue * eps(v) : eps(u)   d_omega
  /
 *----------------------------------------------------------------------*/
icol=0;
for (icn=0;icn<iel;icn++) 
{
   irow=0;
   for (irn=0;irn<iel;irn++)
   {
      estif[irow][icol]     += c*(TWO*derxy[0][irn]*derxy[0][icn] \
                                    + derxy[1][irn]*derxy[1][icn]); 	     
      estif[irow+1][icol]   += c*(    derxy[0][irn]*derxy[1][icn]);
      estif[irow+1][icol+1] += c*(TWO*derxy[1][irn]*derxy[1][icn] \
                                    + derxy[0][irn]*derxy[0][icn]);
      estif[irow][icol+1]   += c*(    derxy[1][irn]*derxy[0][icn]);
      irow += 2;
   } /* end loop over irn */
   icol += 2;
} /* end loop over icn */

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Nc(u):
    /
   |  v * u_old * grad(u)     d_omega
  /
 *----------------------------------------------------------------------*/
if(dynvar->nic != 0) /* evaluate for Newton- and fixed-point-like-iteration */
{
   icol=0;
   for (icn=0;icn<iel;icn++)
   {
      irow=0;  
      for (irn=0;irn<iel;irn++)
      {
         aux = (velint[0]*derxy[0][icn] \
               +velint[1]*derxy[1][icn])*funct[irn]*fac;
	 estif[irow][icol]     += aux;
         estif[irow+1][icol+1] += aux;
         irow += 2;	 
      } /* end loop over irn */
      icol += 2;
   } /* end loop over icn */
} /* endif (dynvar->nic != 0) */

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Nr(u):
    /
   |  v * u * grad(u_old)     d_omega
  /
 *----------------------------------------------------------------------*/
if (dynvar->nir != 0) /* evaluate for Newton iteraton */
{
   icol=0;
   for (icn=0;icn<iel;icn++)
   {
      irow=0;  
      for (irn=0;irn<iel;irn++)
      {
         aux = funct[irn]*funct[icn]*fac;
	 estif[irow][icol]     += aux*vderxy[0][0];
	 estif[irow+1][icol]   += aux*vderxy[1][0];
	 estif[irow+1][icol+1] += aux*vderxy[1][1];
	 estif[irow][icol+1]   += aux*vderxy[0][1];
	 irow += 2;
      } /* end loop over irn */
      icol += 2;
   } /* end loop over icn */
} /* endif (dynvar->nir != 0) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calkvv */

/*!---------------------------------------------------------------------                                         
\brief evaluate galerkin part of Kvp

<pre>                                                         genk 04/02

In this routine the galerkin part of matrix Kvp is calculated:

    /
   |  - div(v) * p     d_omega
  /

    /
   | - q * div(u)      d_omega
  / 

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  				   
      --> Kvp is stored in estif[(0..(2*iel-1)][(2*iel)..(3*iel-1)]
      --> Kpv is stored in estif[((2*iel)..(3*iel-1)][0..(2*iel-1)]
      
</pre>
\param **estif     double	   (i/o)  ele stiffness matrix
\param  *funct     double	   (i)    nat. shape funcs
\param **derxy     double	   (i)    global coord. deriv.
\param   fac       double	   (i)    weighting factor
\param   iel       int  	   (i)    number of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calkvp(
		double         **estif,   
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
 |   posc - since there's only one full element stiffness matrix the    |
 |          column number has to be changed!                            |
/*----------------------------------------------------------------------*/
int     irow, icol,irn,ird;  
int     posc;
double  aux;

#ifdef DEBUG 
dstrc_enter("f2_calkvp");
#endif		

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Kvp:
    /
   |  - div(v) * p     d_omega
  /
  
   and matrix Kpv: 
    /
   | - q * div(u)      d_omega
  /      
 *----------------------------------------------------------------------*/
for (icol=0;icol<iel;icol++)
{
  irow=-1;
  posc = icol + 2*iel;
  for (irn=0;irn<iel;irn++)
  {
     for(ird=0;ird<2;ird++)
     {      
	aux = funct[icol]*derxy[ird][irn]*fac;
	irow++;
	estif[irow][posc] -= aux;
	estif[posc][irow] -= aux;		
     } /* end loop over ird */
  } /* end loop over irn */
} /* end loop over icol */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calkvp */

/*!---------------------------------------------------------------------                                         
\brief evaluate galerkin part of Mvv

<pre>                                                         genk 04/02

In this routine the galerkin part of matrix Mvv is calculated:

    /
   |  v * u    d_omega
  /

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elemass  			      
      --> Mvv is stored in emass[0..(2*iel-1)][0..(2*iel-1)]  
      
</pre>
\param **emass     double	   (i/o)  ele mass matrix
\param  *funct     double	   (i)    nat. shape funcs
\param   fac       double	   (i)    weighting factor
\param   iel       int   	   (i)    number of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calmvv(
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
int     nvdfe;             /* number of velocity dofs of actual element */
double  aux;

#ifdef DEBUG 
dstrc_enter("f2_calmvv");
#endif		

nvdfe = NUM_F2_VELDOF*iel;

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Mvv:
    /
   |  v * u    d_omega
  /
 *----------------------------------------------------------------------*/
icn=-1;
for(icol=0;icol<nvdfe;icol+=2)
{
   icn++;
   irn=-1;
   for(irow=0;irow<nvdfe;irow+=2)
   {
      irn++;
      aux = funct[icn]*funct[irn]*fac;
      emass[irow][icol]     += aux;
      emass[irow+1][icol+1] += aux;
   } /* end loop over irow */
} /* end loop over icol */

/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calmvv */



#endif

/*!----------------------------------------------------------------------
\file
\brief stabilisation part of element stiffness matrix for fluid2

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
/*!---------------------------------------------------------------------
\brief evaluate stabilisaton part of Kvv

<pre>                                                         genk 04/02

In this routine the stabilisation part of matrix Kvv is calculated:

    /
   |  tau_c * div(v) * div(u)   d_omega
  /

    /
   |  tau_mu * u_old * grad(v) * u_old * grad(u)   d_omega
  /
  
    /
   |  tau_mu * u_old * grad(v) * u * grad(u_old)   d_omega
  /  

    /
   |  -tau_mu * 2 * nue * u_old * grad(v) * div(eps(u))   d_omega
  /

    /
   |  +/- tau_mp  * 4 * nue**2 * div(eps(v)) d_omega
  /  

    /
   |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u_old * grad(u) d_omega
  /

    /
   |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u * grad(u_old) d_omega
  /

  
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif
      --> Kvv is stored in estif[0..(2*iel-1)][0..(2*iel-1)]
      
</pre>
\param  *ele	   ELEMENT	   (i)	   actual element
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **estif     double	   (i/o)   ele stiffness matrix
\param  *velint    double	   (i)     vel. at integr. point
\param **vderxy    double	   (i)     global vel. deriv.
\param  *funct     double	   (i)     nat. shape functions
\param **derxy     double	   (i)     global derivatives
\param **derxy2    double	   (i)     2nd global derivatives
\param   fac	   double	   (i)	   weighting factor	   
\param   visc	   double	   (i)	   fluid viscosity
\param   iel	   int		   (i)	   num. of nodes in ele
\param   ihoel	   int		   (i)	   flag for higer ord. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabkvv(			      
                    ELEMENT         *ele,    
		    FLUID_DYN_CALC  *dynvar,
		    double         **estif,  
		    double          *velint, 
		    double         **vderxy, 
		    double          *funct,  
		    double         **derxy,  
		    double         **derxy2, 
		    double           fac,    
		    double           visc,   
		    int              iel,    
                    int              ihoel   
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
double taumu;
double taump;
double tauc;
double c,cc;
double aux,auxr,auxc;
double sign;

#ifdef DEBUG 
dstrc_enter("f2_calstabkvv");
#endif

/*---------------------------------------- set stabilisation parameter */
taumu = dynvar->tau[0];
taump = dynvar->tau[1];
tauc  = dynvar->tau[2];

/*----------------------------------------------------------------------*
   Calculate continuity stabilisation part:
    /
   |  tau_c * div(v) * div(u)   d_omega
  /
 *----------------------------------------------------------------------*/
if (ele->e.f2->icont!=0)
{
   c = fac*tauc;

   icol=0;
   for(icn=0;icn<iel;icn++)
   {
      irow=0;
      for(irn=0;irn<iel;irn++)
      {
         estif[irow][icol]     += derxy[0][icn]*derxy[0][irn]*c;
         estif[irow+1][icol]   += derxy[0][icn]*derxy[1][irn]*c;
         estif[irow][icol+1]   += derxy[1][icn]*derxy[0][irn]*c;
         estif[irow+1][icol+1] += derxy[1][icn]*derxy[1][irn]*c;
         irow += 2;
      } /* end of loop over irn */
      icol += 2;
   } /* end of loop over icn */
} /* endif (ele->e.f2->icont!=0) */
c = fac*taumu;
cc = c;

/*------------------------------ calculate advection stabilisation part */
if (ele->e.f2->iadvec!=0)
{
/*----------------------------------------------------------------------*
   Calculate advection stabilisation part Nc(u):
    /
   |  tau_mu * u_old * grad(v) * u_old * grad(u)   d_omega
  /
 *----------------------------------------------------------------------*/
   if (dynvar->nic!=0) /* evaluate for Newton- and fixed-point-like-iteration */
   {
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
         auxc = (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn])*cc;
	 irow=0;
	 for (irn=0;irn<iel;irn++)
	 {
	    aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn])*auxc;
	    estif[irow][icol]     += aux;
	    estif[irow+1][icol+1] += aux;
	    irow += 2;
	 } /* end of loop over irn */
	 icol += 2;
      } /* end of loop over icn */
   } /* endif (dynvar->nic!=0) */

/*----------------------------------------------------------------------*
   Calculate advection stabilisation part Nr(u):
    /
   |  tau_mu * u_old * grad(v) * u * grad(u_old)   d_omega
  /
 *----------------------------------------------------------------------*/
   if (dynvar->nir!=0)  /* evaluate for Newton iteraton */
   {
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
         auxc = funct[icn]*cc;
	 irow=0;
	 for (irn=0;irn<iel;irn++)
	 {
	    aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn])*auxc;
	    estif[irow][icol]     += aux*vderxy[0][0];
	    estif[irow+1][icol]   += aux*vderxy[1][0];
	    estif[irow][icol+1]   += aux*vderxy[0][1];
	    estif[irow+1][icol+1] += aux*vderxy[1][1];
	    irow += 2;
	 } /* end of loop over irn */
	 icol += 2;
      } /* end of loop over icn */
   } /* endif (dynvar->nir!=0) */

/*----------------------------------------------------------------------*
   Calculate advection stabilisation part for higher order elements:
    /
   |  -tau_mu * 2 * nue * u_old * grad(v) * div(eps(u))   d_omega
  /
 *----------------------------------------------------------------------*/
   if (ihoel!=0)
   {
      cc = c*visc;
      
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
         irow=0;
	 auxc = derxy2[0][icn] + derxy2[1][icn];
	 for (irn=0;irn<iel;irn++)
	 {
	    aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn])*cc;
	    estif[irow][icol]     -= aux*(derxy2[0][icn] + auxc);
	    estif[irow+1][icol]   -= aux* derxy2[2][icn];
	    estif[irow+1][icol+1] -= aux*(derxy2[1][icn] + auxc);
	    estif[irow][icol+1]   -= aux* derxy2[2][icn];
	    irow += 2;
	 } /* end of loop over irn */
         icol += 2;
      } /* end of loop over icn */
   } /* end of stabilisation for higher order elements */
} /* end of advection stabilisation */

/*-------------------------------- calculate viscous stabilisation part */
if (ihoel!=0 && ele->e.f2->ivisc!=0)
{   
   switch (ele->e.f2->ivisc) /* choose stabilisation type --> sign */
   {
   case 1: /* GLS- */
      sign = ONE;
   break;
   case 2: /* GLS+ */
      sign = -ONE;
   break;
   default:
      dserror("viscous stabilisation parameter unknown: IVISC");
   } /* end switch(ele->e.f2->ivisc) */

/*----------------------------------------------------------------------*
   Calculate viscous stabilisation part for higher order elements:
    /
   |  +/- tau_mp  * 4 * nue**2 * div(eps(v)) d_omega
  /
 *----------------------------------------------------------------------*/   
   cc = fac * taump * visc*visc * sign;
   
   icol=0;
   for (icn=0;icn<iel;icn++)
   {
      irow=0;
      auxc = derxy2[0][icn] + derxy2[1][icn];
      for (irn=0;irn<iel;irn++)
      {
         auxr = derxy2[0][irn] + derxy2[1][irn];
	 aux  = auxc*auxr;
	 estif[irow][icol]     +=  (derxy2[0][icn]*(derxy2[0][irn]+auxr) \
	                          + derxy2[2][icn]* derxy2[2][irn]       \
			          + derxy2[0][irn]* auxc + aux)*cc       ;
	 estif[irow+1][icol]   +=  (derxy2[0][icn]* derxy2[2][irn]       \
	                          + derxy2[2][icn]*(derxy2[1][irn]+auxr) \
			          + derxy2[2][irn]* auxc)*cc             ;  
	 estif[irow+1][icol+1] +=  (derxy2[2][icn]* derxy2[2][irn]	 \
	                          + derxy2[1][icn]*(derxy2[1][irn]+auxr) \
			          + derxy2[1][irn]* auxc + aux)*cc       ;
	 estif[irow][icol+1]   +=  (derxy2[2][icn]*(derxy2[0][irn]+auxr) \
	                          + derxy2[1][icn]* derxy2[2][irn]	 \
			          + derxy2[2][irn]* auxc)*cc             ;
	 irow += 2;		      		      
      } /* end of loop over irn */
      icol += 2;
   } /* end of loop over icn */

/*----------------------------------------------------------------------*
   Calculate viscous stabilisation part Nc(u) for higher order elements:
    /
   |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u_old * grad(u) d_omega
  /
 *----------------------------------------------------------------------*/
   cc = fac * taump * visc * sign;
   
   if (dynvar->nic!=0) /* evaluate for Newton- and fixed-point-like-iteration */
   {
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
         irow=0;
	 aux = velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn];
	 for (irn=0;irn<iel;irn++)
	 {
	    auxr = derxy2[0][irn] + derxy2[1][irn];
	    estif[irow][icol]     -= (derxy2[0][irn]+auxr)*aux*cc;
	    estif[irow+1][icol]   -=  derxy2[2][irn]*aux*cc;
	    estif[irow][icol+1]   -=  derxy2[2][irn]*aux*cc;
	    estif[irow+1][icol+1] -= (derxy2[1][irn]+auxr)*aux*cc;
	    irow += 2;
	 } /* end of loop over irn */
	 icol += 2;
      } /* end of loop over icn */
   } /* endif (dynvar->nic!=0) */
   
/*----------------------------------------------------------------------*
   Calculate viscous stabilisation part Nr(u) for higher order elements:
    /
   |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u * grad(u_old) d_omega
  /
 *----------------------------------------------------------------------*/   
   if (dynvar->nir!=0)  /* evaluate for Newton iteraton */
   {
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
         irow=0;
	 aux = funct[icn]*cc;
	 for (irn=0;irn<iel;irn++)
	 {
	    auxr = derxy2[0][irn]*derxy2[1][irn];
	    estif[irow][icol]     -=  (derxy2[0][irn]*vderxy[0][0]   \
	                             + derxy2[2][irn]*vderxy[1][0]   \
				     + auxr*vderxy[0][0])*aux        ;
	    estif[irow+1][icol]   -=  (derxy2[2][irn]*vderxy[0][0]   \
	                             + derxy2[1][irn]*vderxy[1][0]   \
				     + auxr*vderxy[1][0])*aux        ;
	    estif[irow][icol+1]   -=  (derxy2[0][irn]*vderxy[0][1]   \
	                             + derxy2[2][irn]*vderxy[1][1]   \
				     + auxr*vderxy[0][1])*aux        ;
	    estif[irow+1][icol+1] -=  (derxy2[2][irn]*vderxy[0][1]   \
	                             + derxy2[1][irn]*vderxy[1][1]   \
				     + auxr*vderxy[1][1])*aux        ;
	    irow += 2;
	 } /* end of loop over irn */
	 icol += 2;
      } /* end of loop over icn */
   } /* endif (dynvar->nir!=0) */
} /* end of viscous stabilisation for higher order elments */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calstabkvv */

/*!---------------------------------------------------------------------
\brief evaluate stabilisaton part of Kvp

<pre>                                                         genk 04/02

In this routine the stabilisation part of matrix Kvv is calculated:

    /
   |  tau_mu * u_old * grad(v) * grad(p)   d_omega
  /

    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * grad(p)  d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
 NOTE: there's only one elestif 				    
       --> Kvp is stored in estif[(0..(2*iel-1)][(2*iel)..(3*iel-1)]
      
</pre>
\param  *ele	   ELEMENT	   (i)	   actual element
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **estif     double	   (i/o)   ele stiffness matrix
\param  *velint    double	   (i)     vel. at integr. point
\param  *funct     double	   (i)     nat. shape functions
\param **derxy     double	   (i)     global derivatives
\param **derxy2    double	   (i)     2nd global derivatives
\param   fac	   double	   (i)     weighting factor
\param   visc      double	   (i)     fluid viscosity
\param   iel	   int  	   (i)	   num. of nodes in ele
\param   ihoel     int  	   (i)	   flag for higer ord. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabkvp(
                    ELEMENT         *ele,    
		    FLUID_DYN_CALC  *dynvar,
		    double         **estif, 
		    double          *velint,
		    double          *funct, 
		    double         **derxy, 
		    double         **derxy2,
		    double           fac,   
		    double           visc,  
		    int              iel,   
		    int              ihoel   	    
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
int    irow,icol,irn,ird,posc;
double taumu;
double taump;
double c;
double aux;
double sign;

#ifdef DEBUG 
dstrc_enter("f2_calstabkvp");
#endif

/*---------------------------------------- set stabilisation parameter */
taumu = dynvar->tau[0];
taump = dynvar->tau[1];

c = fac * taumu;

/*------------------------------ calculate advection stabilisation part */
if (ele->e.f2->iadvec!=0)
{
/*----------------------------------------------------------------------*
   Calculate advection stabilisation:
    /
   |  tau_mu * u_old * grad(v) * grad(p)   d_omega
  /
 *----------------------------------------------------------------------*/
   for (icol=0;icol<iel;icol++)
   {
      irow=0;
      posc=2*iel+icol;
      for (irn=0;irn<iel;irn++)
      {
         aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn])*c;
	 estif[irow][posc]   += derxy[0][icol]*aux;
	 estif[irow+1][posc] += derxy[1][icol]*aux;
	 irow += 2;
      } /* end loop over irn */
   } /* end loop over icol */
} /* end of advection stabilisation */

/*-------------------------------- calculate viscous stabilisation part */
if (ele->e.f2->ivisc!=0 && ihoel!=0)
{
/*----------------------------------------------------------------------*
   Calculate viscous stabilisation parts for higher order elements:
    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * grad(p)  d_omega
  /
 *----------------------------------------------------------------------*/
   switch (ele->e.f2->ivisc) /* choose stabilisation type --> sign */
   {
   case 1: /* GLS- */
      sign = ONE;
   break;
   case 2: /* GLS+ */
      sign = -ONE;
   break;
   default:
      dserror("viscous stabilisation parameter unknown: IVISC");
   } /* end switch(ele->e.f2->ivisc) */
   c = fac * taump * visc * sign;
   
   for (icol=0;icol<iel;icol++)
   {
      irow=0;
      posc=2*iel+icol;
      for (irn=0;irn<iel;irn++)
      {
         aux = derxy2[0][irn] + derxy2[1][irn];
	 estif[irow][posc]   -=  (derxy2[0][irn]*derxy[0][icol] \
	                        + derxy2[2][irn]*derxy[1][icol] \
			        + aux*derxy[0][icol])*c         ;
	 estif[irow+1][posc] -=  (derxy2[2][irn]*derxy[0][icol] \
	                        + derxy2[1][irn]*derxy[1][icol] \
			        + aux*derxy[1][icol])*c         ;
	 irow += 2;				
      } /* end loop over irn */
   } /* end loop over icol */
} /*------------ end of viscous stabilisation for higher order elements */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calstabkvp */

/*!--------------------------------------------------------------------- 
\brief evaluate stabilisaton part of Mvv

<pre>                                                         genk 04/02

In this routine the stabilisation part of matrix Mvv is calculated:

    /
   |  -/+ tau_mu * u_old * grad(v) * u d_omega
  /
  
    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * u  d_omega
  /  
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  			    
      --> Mvv is stored in emass[0..(2*iel-1)][0..(2*iel-1)]
      
</pre>
\param  *ele	   ELEMENT	   (i)	   actual element
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **emass     double	   (i/o)   ele mass matrix
\param  *velint    double	   (i)     vel. at integr. point
\param  *funct     double	   (i)     nat. shape functions
\param **derxy     double	   (i)     global derivatives
\param **derxy2    double	   (i)     2nd global derivatives
\param   fac	   double	   (i)     weighting factor
\param   visc      double	   (i)     fluid viscosity
\param   iel	   int  	   (i)	   num. of nodes in ele
\param   ihoel     int  	   (i)	   flag for higer ord. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabmvv(
                    ELEMENT         *ele,     
		    FLUID_DYN_CALC  *dynvar,
		    double         **emass,  
		    double          *velint, 
    		    double          *funct,  
		    double         **derxy,  
		    double         **derxy2, 
		    double           fac,    
		    double           visc,   
		    int              iel,    
		    int              ihoel           
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
double taumu;
double taump;
double c,cc;
double aux,auxc;
double sign;

#ifdef DEBUG 
dstrc_enter("f2_calstabmvv");
#endif

/*---------------------------------------- set stabilisation parameter */
taumu = dynvar->tau[0];
taump = dynvar->tau[1];

c = fac * taumu;
cc = c;

/*------------------------------ calculate advection stabilisation part */
if (ele->e.f2->iadvec!=0)
{
/*----------------------------------------------------------------------*
   Calculate convection stabilisation part:
    /
   |  -/+ tau_mu * u_old * grad(v) * u d_omega
  /
 *----------------------------------------------------------------------*/
   icol=0;
   for (icn=0;icn<iel;icn++)
   {
      auxc = funct[icn]*cc;
      irow=0;
      for (irn=0;irn<iel;irn++)
      {
         aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn])*auxc;
	 emass[irow][icol]     += aux;
	 emass[irow+1][icol+1] += aux;
	 irow += 2;
      } /* end loop over irn */
      icol += 2;      
   } /* end loop over icn */
} /* end of advection stabilisation */

/*-------------------------------- calculate viscous stabilisation part */
if (ele->e.f2->ivisc!=0 && ihoel!=0)
{
/*----------------------------------------------------------------------*
   Calculate viscous stabilisation parts for higher order elements:
    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * u  d_omega
  /
 *----------------------------------------------------------------------*/
   switch (ele->e.f2->ivisc) /* choose stabilisation type --> sign */
   {
   case 1: /* GLS- */
      sign = ONE;
   break;
   case 2: /* GLS+ */
      sign = -ONE;
   break;
   default:
      dserror("viscous stabilisation parameter unknown: IVISC");
   } /* end switch(ele->e.f2->ivisc) */
   c = fac * taump * visc * sign;
   
   icol=0;
   for (icn=0;icn<iel;icn++)
   {      
      aux = funct[icn]*c;
      irow=0;
      for(irn=0;irn<iel;irn++)
      {
         emass[irow][icol]     -= (TWO*derxy2[0][irn] + derxy2[1][irn])*aux;
	 emass[irow+1][icol]   -=      derxy2[2][irn]*aux;
	 emass[irow+1][icol+1] -= (TWO*derxy2[1][irn] + derxy2[0][irn])*aux;
	 emass[irow][icol+1]   -=      derxy2[2][irn]*aux;
	 irow += 2;
      } /* end loop over irn */
      icol += 2;
   } /* end loop over icn */
} /* end of viscous stabilisation for higher order elements */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calstabmvv */

/*!--------------------------------------------------------------------- 
\brief evaluate stabilisaton part of Kpv

<pre>                                                         genk 04/02

In this routine the stabilisation part of matrix Kpv is calculated:

    /
   |  - tau_mp * grad(q) * u_old * grad(u) d_omega
  /
  
    /
   |  - tau_mp * grad(q) * u * grad(u_old) d_omega
  /  

    /
   |  tau_mp * 2 * nue *grad(q) * div(eps(u)) d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  				    
      --> Kpv is stored in estif[((2*iel)..(3*iel-1)][0..(2*iel-1)] 
      
</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)      
\param **estif     double	   (i/o)   ele stiffness matrix
\param  *velint    double	   (i)     vel. at integr. point
\param **vderxy    double	   (i)     global vel. deriv.
\param  *funct     double	   (i)     nat. shape functions
\param **derxy     double	   (i)     global derivatives
\param **derxy2    double	   (i)     2nd global derivatives
\param   fac	   double	   (i)     weighting factor
\param   visc      double	   (i)     fluid viscosity
\param   iel	   int  	   (i)	   num. of nodes in ele
\param   ihoel     int  	   (i)	   flag for higer ord. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabkpv(
		    FLUID_DYN_CALC  *dynvar,
		    double         **estif,   
		    double          *velint, 
		    double         **vderxy, 
		    double          *funct,  
		    double         **derxy,  
		    double         **derxy2, 
		    double           fac,    
		    double           visc,   
		    int              iel,    
		    int              ihoel          
                   )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |
 |   posr - since there's only one full element stiffness matrix the    |
 |          row number has to be changed!                               |
/*----------------------------------------------------------------------*/
int    irow,icol,irn,ird,icn,posr;
double c;
double aux;
double sign;
double taump;

#ifdef DEBUG 
dstrc_enter("f2_calstabkpv");
#endif

/*---------------------------------------- set stabilisation parameter */
taump = dynvar->tau[1];

c = fac * taump;

/*----------------------------------------------------------------------*
   Calculate stabilisation part Nc(u):
    /
   |  - tau_mp * grad(q) * u_old * grad(u) d_omega
  /
 *----------------------------------------------------------------------*/
if (dynvar->nic!=0) /* evaluate for Newton- and fixed-point-like-iteration */
{
   icol=0;
   for (icn=0;icn<iel;icn++)
   {
      aux = (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn])*c;
      for (irow=0;irow<iel;irow++)
      {
         posr = irow + 2*iel;
	 estif[posr][icol]   -= derxy[0][irow]*aux;
	 estif[posr][icol+1] -= derxy[1][irow]*aux;
      } /* end loop over irow */
      icol += 2;
   } /* end loop over icn */
} /* endif (dynvar->nic!=0) */ 

/*----------------------------------------------------------------------*
   Calculate stabilisation part Nr(u):
    /
   |  - tau_mp * grad(q) * u * grad(u_old) d_omega
  /
 *----------------------------------------------------------------------*/
if (dynvar->nir!=0) /* evaluate for Newton iteration */
{
   icol=0;
   for (icn=0;icn<iel;icn++)
   {
      aux = funct[icn]*c;
      for (irow=0;irow<iel;irow++)
      {
         posr = irow + 2*iel;
	 estif[posr][icol]   -= aux*(derxy[0][irow]*vderxy[0][0]   \
	                           + derxy[1][irow]*vderxy[1][0])  ;
	 estif[posr][icol+1] -= aux*(derxy[0][irow]*vderxy[0][1]   \
	                           + derxy[1][irow]*vderxy[1][1])  ;
      } /* end loop over irow */
      icol += 2;
   } /* end loop over icn */
} /* endif (dynvar->nir!=0) */

/*----------------------------------------------------------------------*
   Calculate stabilisation part for higher order elements:
    /
   |  tau_mp * 2 * nue *grad(q) * div(eps(u)) d_omega
  /
 *----------------------------------------------------------------------*/
if (ihoel!=0)
{
   c = c * visc;
   
   icol=0;
   for (icn=0;icn<iel;icn++)
   {
      aux = derxy2[0][icn] + derxy2[1][icn];
      for (irow=0;irow<iel;irow++)
      {
         posr = irow + 2*iel;
         estif[posr][icol]   += ((derxy2[0][icn]+aux)*derxy[0][irow]   \
	                        + derxy2[2][icn]     *derxy[1][irow])*c;
         estif[posr][icol+1] +=  (derxy2[2][icn]     *derxy[0][irow]   \
	                        +(derxy2[1][icn]+aux)*derxy[1][irow])*c;
      } /* end loop over irow */
      icol += 2;
   } /* end loop over icn */
} /* end of stabilisation for higher order elements */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calstabkpv */

/*!--------------------------------------------------------------------- 
\brief evaluate stabilisaton part of Kpp

<pre>                                                         genk 04/02

In this routine the stabilisation part of matrix Kpp is calculated:

    /
   |  - tau_mp * grad(q) *grad(p) d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  			     
      --> Kpp is stored in				     
	      estif[((2*iel)..(3*iel-1)][((2*iel)..(3*iel-1)] 
      
</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **estif     double	   (i/o)   ele stiffness matrix
\param **derxy     double	   (i)     global derivatives
\param   fac	   double	   (i)     weighting factor
\param   iel	   int  	   (i)     num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabkpp(
		    FLUID_DYN_CALC  *dynvar,
		    double         **estif,   
		    double         **derxy,  
		    double           fac,    
		    int              iel             
                   )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   posr - since there's only one full element stiffness matrix the    |
 |          row number has to be changed!                               |
 |   posc - since there's only one full element stiffness matrix the    |
 |          column number has to be changed!                            |
/*----------------------------------------------------------------------*/
int    irow,icol,posc,posr;
double c;
double taump;

#ifdef DEBUG 
dstrc_enter("f2_calstabkpp");
#endif

/*---------------------------------------- set stabilisation parameter */
taump = dynvar->tau[1];

c = fac * taump;

/*----------------------------------------------------------------------*
   Calculate stabilisation part for matrix Kpp:
    /
   |  - tau_mp * grad(q) *grad(p) d_omega
  /
 *----------------------------------------------------------------------*/
for (icol=0;icol<iel;icol++)
{
   posc = icol + 2*iel;
   for (irow=0;irow<iel;irow++)
   {
      posr = irow + 2*iel;
      estif[posr][posc] -= (derxy[0][irow]*derxy[0][icol]    \
                           +derxy[1][irow]*derxy[1][icol])*c ;
   } /* end of loop over irow */
} /* end of loop over icol */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calstabkpp */

/*!---------------------------------------------------------------------
\brief evaluate stabilisaton part of Mpv

<pre>                                                         genk 04/02

In this routine the stabilisation part of matrix Mpv is calculated:

    /
   |  - tau_mp * grad(q) * u d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elemass  				    
      --> Mpv is stored in emass[((2*iel)..(3*iel-1)][0..(2*iel-1)] 
      
</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **emass     double	   (i/o)   ele mass matrix
\param  *funct     double	   (i)     nat. shape functions
\param **derxy     double	   (i)     global derivatives
\param   fac	   double	   (i)     weighting factor
\param   iel	   int		   (i)	   num. of nodes in ele

\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabmpv(
		    FLUID_DYN_CALC  *dynvar,
		    double         **emass,   
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
 |   posr - since there's only one full element stiffness matrix the    |
 |          row number has to be changed!                               |
/*----------------------------------------------------------------------*/
int    irow,icol,icn,posr;
double c;
double taump;
double auxc;

#ifdef DEBUG 
dstrc_enter("f2_calstabmpv");
#endif

/*---------------------------------------- set stabilisation parameter */
taump = dynvar->tau[1];

c = fac * taump;

/*----------------------------------------------------------------------*
   Calculate stabilisation part for matrix Mpv:
    /
   |  - tau_mp * grad(q) * u d_omega
  /
 *----------------------------------------------------------------------*/
icol=0;
for (icn=0;icn<iel;icn++)
{
   auxc = funct[icn]*c;
   for (irow=0;irow<iel;irow++)
   {
      posr = irow + 2*iel;
      emass[posr][icol]   -= derxy[0][irow]*auxc;
      emass[posr][icol+1] -= derxy[1][irow]*auxc;
   } /* end loop over irow */
   icol +=2;
} /* end loop over icn */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calstabmpv */



#endif

/*!----------------------------------------------------------------------
\file
\brief stabilisation part of element stiffness matrix for fluid3

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID3 
#include "../headers/standardtypes.h"
#include "fluid3_prototypes.h"
#include "fluid3.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;   
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

static FLUID_DYNAMIC   *fdyn;
/*!---------------------------------------------------------------------
\brief evaluate stabilisaton part of Kvv

<pre>                                                         genk 05/02

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
      --> Kvv is stored in estif[0..(3*iel-1)][0..(3*iel-1)]
      
</pre>
\param  *ele	   ELEMENT	   (i)	   actual element
\param **estif     DOUBLE	   (i/o)   ele stiffness matrix
\param  *velint    DOUBLE	   (i)     vel. at integr. point
\param **vderxy    DOUBLE	   (i)     global vel. deriv.
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param **derxy2    DOUBLE	   (i)     2nd global derivatives
\param   fac	   DOUBLE	   (i)	   weighting factor	   
\param   visc	   DOUBLE	   (i)	   fluid viscosity
\param   iel	   INT		   (i)	   num. of nodes in ele
\param   ihoel	   INT		   (i)	   flag for higer ord. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabkvv(
                    ELEMENT         *ele,
		    DOUBLE         **estif,
		    DOUBLE          *velint,
		    DOUBLE         **vderxy,
		    DOUBLE          *funct,
		    DOUBLE         **derxy,
		    DOUBLE         **derxy2,
		    DOUBLE           fac,
		    DOUBLE           visc,
		    INT              iel,	 
                    INT              ihoel
                   )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |  
 *----------------------------------------------------------------------*/
INT    irow,icol,irn,icn;
DOUBLE taumu;
DOUBLE taump;
DOUBLE tauc;
DOUBLE c,cc;
DOUBLE aux,auxr,auxc;
DOUBLE sign;
STAB_PAR_GLS *gls;	/* pointer to GLS stabilisation parameters	*/

#ifdef DEBUG 
dstrc_enter("f3_calstabkvv");
#endif

/*---------------------------------------------------------- initialise */
fdyn    = alldyn[genprob.numff].fdyn;
gls    = ele->e.f3->stabi.gls;

if (ele->e.f3->stab_type != stab_gls) 
   dserror("routine with no or wrong stabilisation called");

/*---------------------------------------- set stabilisation parameter */
taumu = fdyn->tau[0];
taump = fdyn->tau[1];
tauc  = fdyn->tau[2];

/*----------------------------------------------------------------------*
   Calculate continuity stabilisation part:
    /
   |  tau_c * div(v) * div(u)   d_omega
  /
 *----------------------------------------------------------------------*/
if (gls->icont!=0)
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
         estif[irow+2][icol]   += derxy[0][icn]*derxy[2][irn]*c;

         estif[irow][icol+1]   += derxy[1][icn]*derxy[0][irn]*c;
         estif[irow+1][icol+1] += derxy[1][icn]*derxy[1][irn]*c;
         estif[irow+2][icol+1] += derxy[1][icn]*derxy[2][irn]*c;

         estif[irow][icol+2]   += derxy[2][icn]*derxy[0][irn]*c;
         estif[irow+1][icol+2] += derxy[2][icn]*derxy[1][irn]*c;
         estif[irow+2][icol+2] += derxy[2][icn]*derxy[2][irn]*c;
         irow += 3;
      } /* end of loop over irn */
      icol += 3;
   } /* end of loop over icn */
} /* endif (ele->e.f3->icont!=0) */
c = fac*taumu;
cc = c;
/*------------------------------ calculate advection stabilisation part */
if (gls->iadvec!=0)
{

/*----------------------------------------------------------------------*
   Calculate advection stabilisation part Nc(u):
    /
   |  tau_mu * u_old * grad(v) * u_old * grad(u)   d_omega
  /
 *----------------------------------------------------------------------*/
   if (fdyn->nic!=0) /* evaluate for Newton- and fixed-point-like-iteration */
   {
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
         auxc =  (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn] \
	        + velint[2]*derxy[2][icn])*cc;
	 irow=0;
	 for (irn=0;irn<iel;irn++)
	 {
	    aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn] \
	         + velint[2]*derxy[2][irn])*auxc;
	    estif[irow][icol]     += aux;
	    estif[irow+1][icol+1] += aux;
	    estif[irow+2][icol+2] += aux;
	    irow += 3;
	 } /* end of loop over irn */
	 icol += 3;
      } /* end of loop over icn */
   } /* endif (fdyn->nic!=0) */

/*----------------------------------------------------------------------*
   Calculate advection stabilisation part Nr(u):
    /
   |  tau_mu * u_old * grad(v) * u * grad(u_old)   d_omega
  /
 *----------------------------------------------------------------------*/
   if (fdyn->nir!=0)  /* evaluate for Newton iteraton */
   {
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
         auxc = funct[icn]*cc;
	 irow=0;
	 for (irn=0;irn<iel;irn++)
	 {
	    aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn] \
	         + velint[2]*derxy[2][irn])*auxc;

	    estif[irow][icol]     += aux*vderxy[0][0];
	    estif[irow+1][icol]   += aux*vderxy[1][0];
	    estif[irow+2][icol]   += aux*vderxy[2][0];
	    
	    estif[irow][icol+1]   += aux*vderxy[0][1];
	    estif[irow+1][icol+1] += aux*vderxy[1][1];
	    estif[irow+2][icol+1] += aux*vderxy[2][1];	    

	    estif[irow][icol+2]   += aux*vderxy[0][2];
	    estif[irow+1][icol+2] += aux*vderxy[1][2];
	    estif[irow+2][icol+2] += aux*vderxy[2][2];	 
	    irow += 3;
	 } /* end of loop over irn */
	 icol += 3;
      } /* end of loop over icn */
   } /* endif (fdyn->nir!=0) */

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
	 auxc = derxy2[0][icn] + derxy2[1][icn] + derxy2[2][icn];
	 for (irn=0;irn<iel;irn++)
	 {
	    aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn] \
	         + velint[2]*derxy[2][irn])*cc;

	    estif[irow][icol]     -= aux*(derxy2[0][icn] + auxc);
	    estif[irow+1][icol]   -= aux* derxy2[3][icn];
	    estif[irow+2][icol]   -= aux* derxy2[4][icn];

	    estif[irow][icol+1]   -= aux* derxy2[3][icn];
	    estif[irow+1][icol+1] -= aux*(derxy2[1][icn] + auxc);
	    estif[irow+2][icol+1] -= aux* derxy2[5][icn];

	    estif[irow][icol+2]   -= aux* derxy2[4][icn];
	    estif[irow+1][icol+2] -= aux* derxy2[5][icn];
	    estif[irow+2][icol+2] -= aux*(derxy2[2][icn] + auxc);
	    irow += 3;
	 } /* end of loop over irn */
         icol += 3;
      } /* end of loop over icn */
   } /* end of stabilisation for higher order elements */
} /* end of advection stabilisation */

/*-------------------------------- calculate viscous stabilisation part */
if (ihoel!=0 && gls->ivisc!=0)
{   
   switch (gls->ivisc) /* choose stabilisation type --> sign */
   {
   case 1: /* GLS- */
      sign = ONE;
   break;
   case 2: /* GLS+ */
      sign = -ONE;
   break;
   default:
      dserror("viscous stabilisation parameter unknown: IVISC");
   } /* end switch (ele->e.f3->ivisc) */

/*----------------------------------------------------------------------*
   Calculate viscous stabilisation part for higher order elements:
    /
   |  +/- tau_mp  * 4 * nue**2 * div(eps(v)) * div(eps(u))d_omega
  /
 *----------------------------------------------------------------------*/   
   cc = fac * taump * visc*visc * sign;

   icol=0;
   for (icn=0;icn<iel;icn++)
   {
      irow=0;
      auxc = derxy2[0][icn] + derxy2[1][icn] + derxy2[2][icn] ;
      for (irn=0;irn<iel;irn++)
      {
         auxr = derxy2[0][irn] + derxy2[1][irn] + derxy2[2][irn];

	 estif[irow][icol]     += ( (auxc + derxy2[0][icn])*(auxr + derxy2[0][irn])  
	                          + derxy2[3][icn]*derxy2[3][irn]      
				  + derxy2[4][icn]*derxy2[4][irn] )*cc;	 
	 estif[irow+1][icol]   += ( derxy2[3][icn]*(auxr + derxy2[0][irn])            
	                          + (auxc + derxy2[1][icn])*derxy2[3][irn]  
				  + derxy2[5][icn]*derxy2[4][irn] )*cc;
	 estif[irow+2][icol]   += ( derxy2[4][icn]*(auxr + derxy2[0][irn])  
	                          + derxy2[5][icn]*derxy2[3][irn]   
				  + (auxc + derxy2[2][icn])*derxy2[4][irn] )*cc; 

	 estif[irow][icol+1]   += ( (auxc + derxy2[0][icn])*derxy2[3][irn] 
	                          + derxy2[3][icn]*(auxr + derxy2[1][irn]) 
				  + derxy2[4][icn]*derxy2[5][irn] )*cc;	 
	 estif[irow+1][icol+1] += ( derxy2[3][icn]*derxy2[3][irn] 
	                          + (auxc + derxy2[1][icn])*(auxr + derxy2[1][irn]) 
				  + derxy2[5][icn]*derxy2[5][irn] )*cc;	
	 estif[irow+2][icol+1] += ( derxy2[4][icn]*derxy2[3][irn] 
	                          + derxy2[5][icn]*(auxr + derxy2[1][irn]) 
				  + (auxc + derxy2[2][icn])*derxy2[5][irn] )*cc;
			          				  
	 estif[irow][icol+2]   += ( (auxc + derxy2[0][icn])*derxy2[4][irn] 
	                          + derxy2[3][icn]*derxy2[5][irn] 
				  + derxy2[4][icn]*(auxr + derxy2[2][irn]) )*cc;
	 estif[irow+1][icol+2] += ( derxy2[3][icn]*derxy2[4][irn] 
	                          + (auxc + derxy2[1][icn])*derxy2[5][irn] 
				  + derxy2[5][icn]*(auxr + derxy2[2][irn]) )*cc; 

	 estif[irow+2][icol+2] += ( derxy2[4][icn]*derxy2[4][irn] 
	                          + derxy2[5][icn]*derxy2[5][irn] 
				  + (auxc + derxy2[2][irn])*(auxr + derxy2[2][irn]) )*cc;	          
	 irow += 3;		      		      
      } /* end of loop over irn */
      icol += 3;
   } /* end of loop over icn */

/*----------------------------------------------------------------------*
   Calculate viscous stabilisation part Nc(u) for higher order elements:
    /
   |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u_old * grad(u) d_omega
  /
 *----------------------------------------------------------------------*/
   cc = fac * taump * visc * sign;

   if (fdyn->nic!=0) /* evaluate for Newton- and fixed-point-like-iteration */
   {
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
         irow=0;
	 aux = (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn] \
	      + velint[2]*derxy[2][icn])*cc;
	 for (irn=0;irn<iel;irn++)
	 {
	    auxr = derxy2[0][irn] + derxy2[1][irn] + derxy2[2][irn];
	    estif[irow][icol]     -= (derxy2[0][irn]+auxr)*aux;
	    estif[irow+1][icol]   -=  derxy2[3][irn]*aux;
	    estif[irow+2][icol]   -=  derxy2[4][irn]*aux;

	    estif[irow][icol+1]   -=  derxy2[3][irn]*aux;
	    estif[irow+1][icol+1] -= (derxy2[1][irn]+auxr)*aux;
	    estif[irow+2][icol+1] -= (derxy2[5][irn])*aux;

	    estif[irow][icol+2]   -=  derxy2[4][irn]*aux;
	    estif[irow+1][icol+2] -= (derxy2[5][irn])*aux;
	    estif[irow+2][icol+2] -= (derxy2[2][irn]+auxr)*aux;
	    irow += 3;
	 } /* end of loop over irn */
	 icol += 3;
      } /* end of loop over icn */
   } /* endif (fdyn->nic!=0) */
/*----------------------------------------------------------------------*
   Calculate viscous stabilisation part Nr(u) for higher order elements:
    /
   |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u * grad(u_old) d_omega
  /
 *----------------------------------------------------------------------*/   
   if (fdyn->nir!=0)  /* evaluate for Newton iteraton */
   {
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
         irow=0;
	 aux = funct[icn]*cc;
	 for (irn=0;irn<iel;irn++)
	 {
	    auxr = derxy2[0][irn]*derxy2[1][irn]*derxy2[2][irn];
	    estif[irow][icol]     -=  (derxy2[0][irn]*vderxy[0][0]   \
	                             + derxy2[3][irn]*vderxy[1][0]   \
	                             + derxy2[4][irn]*vderxy[2][0]   \
				     + auxr*vderxy[0][0])*aux        ;
	    estif[irow+1][icol]   -=  (derxy2[3][irn]*vderxy[0][0]   \
	                             + derxy2[1][irn]*vderxy[1][0]   \
	                             + derxy2[5][irn]*vderxy[2][0]   \
				     + auxr*vderxy[1][0])*aux        ;
	    estif[irow+2][icol]   -=  (derxy2[4][irn]*vderxy[0][0]   \
	                             + derxy2[5][irn]*vderxy[1][0]   \
	                             + derxy2[2][irn]*vderxy[2][0]   \
				     + auxr*vderxy[2][0])*aux        ;

	    estif[irow][icol+1]   -=  (derxy2[0][irn]*vderxy[0][1]   \
	                             + derxy2[3][irn]*vderxy[1][1]   \
	                             + derxy2[4][irn]*vderxy[2][1]   \
				     + auxr*vderxy[0][1])*aux        ;
	    estif[irow+1][icol+1] -=  (derxy2[3][irn]*vderxy[0][1]   \
	                             + derxy2[1][irn]*vderxy[1][1]   \
	                             + derxy2[5][irn]*vderxy[2][1]   \
				     + auxr*vderxy[1][1])*aux        ;
	    estif[irow+2][icol+1] -=  (derxy2[4][irn]*vderxy[0][1]   \
	                             + derxy2[5][irn]*vderxy[1][1]   \
	                             + derxy2[2][irn]*vderxy[2][1]   \
				     + auxr*vderxy[2][1])*aux        ;

	    estif[irow][icol+2]   -=  (derxy2[0][irn]*vderxy[0][2]   \
	                             + derxy2[3][irn]*vderxy[1][2]   \
	                             + derxy2[4][irn]*vderxy[2][2]   \
				     + auxr*vderxy[0][2])*aux        ;
	    estif[irow+1][icol+2] -=  (derxy2[3][irn]*vderxy[0][2]   \
	                             + derxy2[1][irn]*vderxy[1][2]   \
	                             + derxy2[5][irn]*vderxy[2][2]   \
				     + auxr*vderxy[1][2])*aux        ;
	    estif[irow+2][icol+2] -=  (derxy2[4][irn]*vderxy[0][2]   \
	                             + derxy2[5][irn]*vderxy[1][2]   \
	                             + derxy2[2][irn]*vderxy[2][2]   \
				     + auxr*vderxy[2][2])*aux        ;
	    irow += 3;
	 } /* end of loop over irn */
	 icol += 3;
      } /* end of loop over icn */
   } /* endif (fdyn->nir!=0) */
} /* end of viscous stabilisation for higher order elments */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calstabkvv */ 

/*!---------------------------------------------------------------------
\brief evaluate stabilisaton part of Kvp

<pre>                                                         genk 05/02

In this routine the stabilisation part of matrix Kvv is calculated:

    /
   |  tau_mu * u_old * grad(v) * grad(p)   d_omega
  /

    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * grad(p)  d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  				   
      --> Kvp is stored in estif[(0..(3*iel-1)][(3*iel)..(4*iel-1)]
      
</pre>
\param  *ele	   ELEMENT	   (i)	   actual element
\param **estif     DOUBLE	   (i/o)   ele stiffness matrix
\param  *velint    DOUBLE	   (i)     vel. at integr. point
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param **derxy2    DOUBLE	   (i)     2nd global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   visc      DOUBLE	   (i)     fluid viscosity
\param   iel	   INT  	   (i)	   num. of nodes in ele
\param   ihoel     INT  	   (i)	   flag for higer ord. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabkvp(
                    ELEMENT         *ele,
		    DOUBLE         **estif,
		    DOUBLE          *velint,
		    DOUBLE          *funct,
		    DOUBLE         **derxy,
		    DOUBLE         **derxy2,
		    DOUBLE           fac,
		    DOUBLE           visc,
		    INT              iel,
		    INT              ihoel		 
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
 *----------------------------------------------------------------------*/
INT    irow,icol,irn,posc;
DOUBLE taumu;
DOUBLE taump;
DOUBLE c;
DOUBLE aux;
DOUBLE sign;
STAB_PAR_GLS *gls;	/* pointer to GLS stabilisation parameters	*/

#ifdef DEBUG 
dstrc_enter("f3_calstabkvp");
#endif	

/*---------------------------------------------------------- initialise */
fdyn   = alldyn[genprob.numff].fdyn;
gls    = ele->e.f3->stabi.gls;

if (ele->e.f3->stab_type != stab_gls) 
   dserror("routine with no or wrong stabilisation called");

/*---------------------------------------- set stabilisation parameter */
taumu = fdyn->tau[0];
taump = fdyn->tau[1];

c = fac * taumu;

/*------------------------------ calculate advection stabilisation part */
if (gls->iadvec!=0)
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
      posc=3*iel+icol;
      for (irn=0;irn<iel;irn++)
      {
         aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn] \
	      + velint[2]*derxy[2][irn])*c;
	 estif[irow][posc]   += derxy[0][icol]*aux;
	 estif[irow+1][posc] += derxy[1][icol]*aux;
	 estif[irow+2][posc] += derxy[2][icol]*aux;
	 irow += 3;
      } /* end of loop over irn */
   } /* end of loop over icol */
} /* end of advection stabilisation */

/*-------------------------------- calculate viscous stabilisation part */
if (gls->ivisc!=0 && ihoel!=0)
{
/*----------------------------------------------------------------------*
   Calculate viscous stabilisation parts for higher order elements:
    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * grad(p)  d_omega
  /
 *----------------------------------------------------------------------*/
   switch (gls->ivisc) /* choose stabilisation type --> sign */
   {
   case 1: /* GLS- */
      sign = ONE;
   break;
   case 2: /* GLS+ */
      sign = -ONE;
   break;
   default:
      dserror("viscous stabilisation parameter unknown: IVISC");
   } /* end switch (ele->e.f3->ivisc) */
   c = fac * taump * visc * sign;
 
   for (icol=0;icol<iel;icol++)
   {
      irow=0;
      posc=3*iel+icol;
      for (irn=0;irn<iel;irn++)
      {
         aux = derxy2[0][irn] + derxy2[1][irn] + derxy2[2][irn];
	 estif[irow][posc]   -=  (derxy2[0][irn]*derxy[0][icol] \
	                        + derxy2[3][irn]*derxy[1][icol] \
	                        + derxy2[4][irn]*derxy[2][icol] \
			        + aux*derxy[0][icol])*c         ;
	 estif[irow+1][posc] -=  (derxy2[3][irn]*derxy[0][icol] \
	                        + derxy2[1][irn]*derxy[1][icol] \
	                        + derxy2[5][irn]*derxy[2][icol] \
			        + aux*derxy[1][icol])*c         ;
	 estif[irow+2][posc] -=  (derxy2[4][irn]*derxy[0][icol] \
	                        + derxy2[5][irn]*derxy[1][icol] \
	                        + derxy2[2][irn]*derxy[2][icol] \
			        + aux*derxy[2][icol])*c         ;
	 irow += 3;				
      } /* end of loop over irn */
   } /* end of loop over icol */
} /*------------ end of viscous stabilisation for higher order elements */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calstabkvp */

/*!--------------------------------------------------------------------- 
\brief evaluate stabilisaton part of Mvv

<pre>                                                         genk 05/02

In this routine the stabilisation part of matrix Mvv is calculated:

    /
   |  -/+ tau_mu * u_old * grad(v) * u d_omega
  /
  
    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * u  d_omega
  /  
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  			    
      --> Mvv is stored in estif[0..(3*iel-1)][0..(3*iel-1)]
      
</pre>
\param  *ele	   ELEMENT	   (i)	   actual element
\param **emass     DOUBLE	   (i/o)   ele mass matrix
\param  *velint    DOUBLE	   (i)     vel. at integr. point
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param **derxy2    DOUBLE	   (i)     2nd global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   visc      DOUBLE	   (i)     fluid viscosity
\param   iel	   INT  	   (i)	   num. of nodes in ele
\param   ihoel     INT  	   (i)	   flag for higer ord. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabmvv(
                    ELEMENT         *ele,
		    DOUBLE         **estif,
		    DOUBLE          *velint,
    		    DOUBLE          *funct,
		    DOUBLE         **derxy,
		    DOUBLE         **derxy2,
		    DOUBLE           fac,
		    DOUBLE           visc,
		    INT              iel,
		    INT              ihoel		 
                   )
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |   
 *----------------------------------------------------------------------*/
INT    irow,icol,irn,icn;
DOUBLE taumu;
DOUBLE taump;
DOUBLE c,cc;
DOUBLE aux,auxc,auxr;
DOUBLE sign;
STAB_PAR_GLS *gls;	/* pointer to GLS stabilisation parameters	*/

#ifdef DEBUG 
dstrc_enter("f3_calstabmvv");
#endif

/*---------------------------------------------------------- initialise */
fdyn   = alldyn[genprob.numff].fdyn;
gls    = ele->e.f3->stabi.gls;

if (ele->e.f3->stab_type != stab_gls) 
   dserror("routine with no or wrong stabilisation called");

/*---------------------------------------- set stabilisation parameter */
taumu = fdyn->tau[0];
taump = fdyn->tau[1];

c = fac * taumu;
cc = c;

/*------------------------------ calculate advection stabilisation part */
if (gls->iadvec!=0)
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
         aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn] \
	      + velint[2]*derxy[2][irn] )*auxc;
	 estif[irow][icol]     += aux;
	 estif[irow+1][icol+1] += aux;
	 estif[irow+2][icol+2] += aux;
	 irow += 3;
      } /* end of loop over irn */
      icol += 3;      
   } /* end of loop over icn */
} /* end of advection stabilisation */

/*-------------------------------- calculate viscous stabilisation part */
if (gls->ivisc!=0 && ihoel!=0)
{
/*----------------------------------------------------------------------*
   Calculate viscous stabilisation parts for higher order elements:
    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * u  d_omega
  /
 *----------------------------------------------------------------------*/
   switch (gls->ivisc) /* choose stabilisation type --> sign */
   {
   case 1: /* GLS- */
      sign = ONE;
   break;
   case 2: /* GLS+ */
      sign = -ONE;
   break;
   default:
      dserror("viscous stabilisation parameter unknown: IVISC");
   } /* end switch (ele->e.f3->ivisc) */
   c = fac * taump * visc * sign;
   
   icol=0;
   for (icn=0;icn<iel;icn++)
   {      
      aux = funct[icn]*c;
      irow=0;
      for(irn=0;irn<iel;irn++)
      {
         auxr = derxy2[0][irn] + derxy2[1][irn] + derxy2[2][irn];
	 estif[irow][icol]     -= (derxy2[0][irn] + auxr)*aux;
	 estif[irow+1][icol]   -=  derxy2[3][irn]*aux;
	 estif[irow+2][icol]   -=  derxy2[4][irn]*aux;

	 estif[irow][icol+1]   -=  derxy2[3][irn]*aux;
	 estif[irow+1][icol+1] -= (derxy2[1][irn] + auxr)*aux;
	 estif[irow+2][icol+1] -=  derxy2[5][irn]*aux;

	 estif[irow][icol+2]   -=  derxy2[4][irn]*aux;
	 estif[irow+1][icol+2] -=  derxy2[5][irn]*aux;
	 estif[irow+2][icol+2] -= (derxy2[2][irn] + auxr)*aux;
	 irow += 3;
      } /* end of loop over irn */
      icol += 3;
   } /* end of loop over icn */
} /* end of viscous stabilisation for higher order elements */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calstabmvv */

/*!--------------------------------------------------------------------- 
\brief evaluate stabilisaton part of Kpv

<pre>                                                         genk 05/02

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
      --> Kpv is stored in estif[((3*iel)..(4*iel-1)][0..(3*iel-1)] 
      
</pre>
\param **estif     DOUBLE	   (i/o)   ele stiffness matrix
\param  *velint    DOUBLE	   (i)     vel. at integr. point
\param **vderxy    DOUBLE	   (i)     global vel. deriv.
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param **derxy2    DOUBLE	   (i)     2nd global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   visc      DOUBLE	   (i)     fluid viscosity
\param   iel	   INT  	   (i)	   num. of nodes in ele
\param   ihoel     INT  	   (i)	   flag for higer ord. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabkpv(
		    DOUBLE         **estif,
		    DOUBLE          *velint,
		    DOUBLE         **vderxy,
		    DOUBLE          *funct,
		    DOUBLE         **derxy,
		    DOUBLE         **derxy2,
		    DOUBLE           fac,
		    DOUBLE           visc,
		    INT              iel,
		    INT              ihoel		 
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
 *----------------------------------------------------------------------*/
INT    irow,icol,icn,posr;
DOUBLE c;
DOUBLE aux;
DOUBLE taump;

#ifdef DEBUG 
dstrc_enter("f3_calstabkpv");
#endif

/*---------------------------------------- set stabilisation parameter */
fdyn  = alldyn[genprob.numff].fdyn;
taump = fdyn->tau[1];

c = fac * taump;

/*----------------------------------------------------------------------*
   Calculate stabilisation part Nc(u):
    /
   |  - tau_mp * grad(q) * u_old * grad(u) d_omega
  /
 *----------------------------------------------------------------------*/
if (fdyn->nic!=0) /* evaluate for Newton- and fixed-point-like-iteration */
{
   icol=0;
   for (icn=0;icn<iel;icn++)
   {
      aux = (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn] \
           + velint[2]*derxy[2][icn])*c;
      for (irow=0;irow<iel;irow++)
      {
         posr = irow + 3*iel;
	 estif[posr][icol]   -= derxy[0][irow]*aux;
	 estif[posr][icol+1] -= derxy[1][irow]*aux;
	 estif[posr][icol+2] -= derxy[2][irow]*aux;
      } /* end of loop over irow */
      icol += 3;
   } /* end of loop over icn */
} /* endif (fdyn->nic!=0) */

/*----------------------------------------------------------------------*
   Calculate stabilisation part Nr(u):
    /
   |  - tau_mp * grad(q) * u * grad(u_old) d_omega
  /
 *----------------------------------------------------------------------*/
if (fdyn->nir!=0) /* evaluate for Newton iteration */
{
   icol=0;
   for (icn=0;icn<iel;icn++)
   {
      aux = funct[icn]*c;
      for (irow=0;irow<iel;irow++)
      {
         posr = irow + 3*iel;
	 estif[posr][icol]   -= aux*(derxy[0][irow]*vderxy[0][0]   \
                                   + derxy[1][irow]*vderxy[1][0]   \
	                           + derxy[2][irow]*vderxy[2][0])  ;
	 estif[posr][icol+1] -= aux*(derxy[0][irow]*vderxy[0][1]   \
                                   + derxy[1][irow]*vderxy[1][1]   \
	                           + derxy[2][irow]*vderxy[2][1])  ;
	 estif[posr][icol+2] -= aux*(derxy[0][irow]*vderxy[0][2]   \
                                   + derxy[1][irow]*vderxy[1][2]   \
	                           + derxy[2][irow]*vderxy[2][2])  ;
      } /* end of loop over irow */
      icol += 3;
   } /* end of loop over icn */
} /* endif (fdyn->nir!=0) */

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
      aux = derxy2[0][icn] + derxy2[1][icn] + derxy2[2][icn];
      for (irow=0;irow<iel;irow++)
      {
         posr = irow + 3*iel;
         estif[posr][icol]   += ((derxy2[0][icn]+aux)*derxy[0][irow]   \
	                        + derxy2[3][icn]     *derxy[1][irow]   \
	                        + derxy2[4][icn]     *derxy[2][irow])*c;

         estif[posr][icol+1] += ( derxy2[3][icn]     *derxy[0][irow]   \
	                        +(derxy2[1][icn]+aux)*derxy[1][irow]   \
	                        + derxy2[5][icn]     *derxy[2][irow])*c;

         estif[posr][icol+2] += ( derxy2[4][icn]     *derxy[0][irow]   \
	                        + derxy2[5][icn]     *derxy[1][irow]   \
	                        +(derxy2[2][icn]+aux)*derxy[2][irow])*c;
      } /* end of loop over irow */
      icol += 3;
   } /* end of loop over icn */
} /* end of stabilisation for higher order elements */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calstabkpv */

/*!--------------------------------------------------------------------- 
\brief evaluate stabilisaton part of Kpp

<pre>                                                         genk 05/02

In this routine the stabilisation part of matrix Kpp is calculated:

    /
   |  - tau_mp * grad(q) *grad(p) d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  			     
      --> Kpp is stored in				     
	      estif[((3*iel)..(4*iel-1)][((3*iel)..(4*iel-1)] 
      
</pre>
\param **estif     DOUBLE	   (i/o)   ele stiffness matrix
\param **derxy     DOUBLE	   (i)     global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   iel	   INT  	   (i)     num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabkpp(
		    DOUBLE         **estif,
		    DOUBLE         **derxy,
		    DOUBLE           fac,
		    INT              iel		 
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
 *----------------------------------------------------------------------*/
INT    irow,icol,posc,posr;
DOUBLE c;
DOUBLE taump;

#ifdef DEBUG 
dstrc_enter("f3_calstabkpp");
#endif

/*---------------------------------------- set stabilisation parameter */
fdyn  = alldyn[genprob.numff].fdyn;
taump = fdyn->tau[1];

c = fac * taump;

/*----------------------------------------------------------------------*
   Calculate stabilisation part for matrix Kpp:
    /
   |  - tau_mp * grad(q) *grad(p) d_omega
  /
 *----------------------------------------------------------------------*/
for (icol=0;icol<iel;icol++)
{
   posc = icol + 3*iel;
   for (irow=0;irow<iel;irow++)
   {
      posr = irow + 3*iel;
      estif[posr][posc] -= (derxy[0][irow]*derxy[0][icol]    \
                           +derxy[1][irow]*derxy[1][icol]    \
                           +derxy[2][irow]*derxy[2][icol])*c ;
   } /* end of loop over irow */
} /* end of loop over icol */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calstabkpp */

/*!---------------------------------------------------------------------
\brief evaluate stabilisaton part of Mpv

<pre>                                                         genk 05/02

In this routine the stabilisation part of matrix Mpv is calculated:

    /
   |  - tau_mp * grad(q) * u d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  				   
      --> Mpv is stored in estif[((3*iel)..(4*iel-1)][0..(3*iel-1)]
      
</pre>
\param **emass     DOUBLE	   (i/o)   ele mass matrix
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   iel	   INT		   (i)	   num. of nodes in ele

\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabmpv(
		    DOUBLE         **estif,
		    DOUBLE          *funct,
		    DOUBLE         **derxy,
		    DOUBLE           fac,
		    INT              iel		 
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
 *----------------------------------------------------------------------*/
INT    irow,icol,icn,posr;
DOUBLE c;
DOUBLE taump;
DOUBLE auxc;

#ifdef DEBUG 
dstrc_enter("f3_calstabmpv");
#endif

/*---------------------------------------- set stabilisation parameter */
fdyn  = alldyn[genprob.numff].fdyn;
taump = fdyn->tau[1];

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
      posr = irow + 3*iel;
      estif[posr][icol]   -= derxy[0][irow]*auxc;
      estif[posr][icol+1] -= derxy[1][irow]*auxc;
      estif[posr][icol+2] -= derxy[2][irow]*auxc;
   } /* end of loop over irow */
   icol += 3;
} /* end of loop over icn */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calstabmpv */



#endif

/*!----------------------------------------------------------------------
\file
\brief evaluate stabilization part of submesh matrices for fluid3

<pre>
Maintainer: Volker Gravemeier
            vgravem@stanford.edu
            
            
</pre>

------------------------------------------------------------------------*/
#ifdef FLUID3_ML 
#include "../headers/standardtypes.h"
#include "fluid3ml_prototypes.h"
#include "../fluid3/fluid3.h"
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

static FLUID_DYNAMIC *fdyn;
/*!---------------------------------------------------------------------                                         
\brief evaluate stabilization part of sm stiffness matrix SMK for fluid3

<pre>                                                       gravem 07/03

In this routine, the stabilization part of the submesh stiffness matrix 
SMK is calculated.

</pre>
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smestif   DOUBLE	   (i/o)  submesh ele stiffness matrix
\param  *velint    DOUBLE	   (i)    velocity at int. point
\param **vderxy    DOUBLE	   (i)    global vel. deriv. at int. p.
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun. 
\param **smderxy2  DOUBLE	   (i)    sm 2nd global deriv. of sh. fun. 
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc      DOUBLE	   (i)    fluid viscosity	     
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\param   ihoelsm   INT  	   (i)	  flag for higher-order elements
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabsmk(FLUID_DYN_ML    *mlvar, 
		   DOUBLE         **smestif,  
		   DOUBLE          *velint, 
		   DOUBLE         **vderxy, 
		   DOUBLE          *smfunct,  
		   DOUBLE         **smderxy,  
		   DOUBLE         **smderxy2, 
		   DOUBLE           fac,    
		   DOUBLE           visc,   
		   INT              smiel,    
                   INT              ihoelsm)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
/*----------------------------------------------------------------------*/
INT    irow,icol;
DOUBLE tau;
DOUBLE con,ccon,beta,divv,cb,cbb,ccb;
DOUBLE aux,auxc;
DOUBLE sign;

#ifdef DEBUG 
dstrc_enter("f3_calstabsmk");
#endif

fdyn = alldyn[genprob.numff].fdyn;
/*---------------------------------------- set stabilization parameter */
tau = mlvar->smtau;
con = fac*tau;

switch (fdyn->conte) /* determine parameter beta of convective term */
{
case 0: 
  beta = ZERO;
break;
case 1: 
  beta = ONE;
  divv= vderxy[0][0]+vderxy[1][1]+vderxy[2][2]; 
break;
case 2: 
  beta = ONE/TWO;
  divv= vderxy[0][0]+vderxy[1][1]+vderxy[2][2]; 
break;
default:
   dserror("unknown form of convective term");
} /* end switch (fdyn->conte) */

if (mlvar->smstado<0) sign = -ONE;  /* USFEM */
else                  sign = ONE;   /* GLS- */

/*----------------------------------------------------------------------*
   Calculate convective stabilization part:
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
    /
   |  tau * u_old * grad(w) * u_old * grad(bub)   d_omega
  /
 *----------------------------------------------------------------------*/
for (icol=0;icol<smiel;icol++)
{
  auxc = (velint[0]*smderxy[0][icol]+velint[1]*smderxy[1][icol]\
         +velint[2]*smderxy[2][icol])*con;
  for (irow=0;irow<smiel;irow++)
  {
    aux = (velint[0]*smderxy[0][irow]+velint[1]*smderxy[1][irow]\
          +velint[2]*smderxy[2][icol])*auxc;
    smestif[irow][icol] += aux;
  } /* end of loop over irow */
} /* end of loop over icol */

if (fdyn->conte!=0)  
{
  cb  = con*beta;
  cbb = cb*beta*sign;
/*----------------------------------------------------------------------*
    /
   |  beta * tau * u_old * grad(w) * bub * div(u_old)   d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<smiel;icol++)
  {
    auxc = smfunct[icol]*divv*cb;
    for (irow=0;irow<smiel;irow++)
    {
      aux = (velint[0]*smderxy[0][irow]+velint[1]*smderxy[1][irow]\
            +velint[2]*smderxy[2][icol])*auxc;
      smestif[irow][icol] += aux;
    } /* end of loop over irow */
  } /* end of loop over icol */
/*----------------------------------------------------------------------*
    /
   |  +/- beta * tau * w * div(u_old) * u_old * grad(bub)   d_omega
  /
 *----------------------------------------------------------------------*/
  cb = cb*sign;
  for (icol=0;icol<smiel;icol++)
  {
    aux = (velint[0]*smderxy[0][icol]+velint[1]*smderxy[1][icol]\
          +velint[2]*smderxy[2][icol])*divv*cb;
    for (irow=0;irow<smiel;irow++)
    {
      smestif[irow][icol] += aux*smfunct[irow];
    } /* end of loop over irow */
  } /* end of loop over icol */
/*----------------------------------------------------------------------*
    /
   |  +/- beta * beta * tau * w * div(u_old) * bub * div(u_old)   d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<smiel;icol++)
  {
    aux = smfunct[icol]*divv*divv*cbb;
    for (irow=0;irow<smiel;irow++)
    {
      smestif[irow][icol] += aux*smfunct[irow];
    } /* end of loop over irow */
  } /* end of loop over icol */
}

/*----------------------------------------------------------------------*
   Calculate convective stabilization part for higher-order elements:
 *----------------------------------------------------------------------*/
if (ihoelsm!=0)
{
  ccon = con*visc;
/*----------------------------------------------------------------------*
    /
   |  - tau * u_old * grad(w) * nue * delta(bub)   d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<smiel;icol++)
  {
    auxc = (smderxy2[0][icol]+smderxy2[1][icol]+smderxy2[2][icol])*ccon;
    for (irow=0;irow<smiel;irow++)
    {
      aux = (velint[0]*smderxy[0][irow]+velint[1]*smderxy[1][irow]\
            +velint[2]*smderxy[2][irow])*auxc;
      smestif[irow][icol] -= aux;
    } /* end of loop over irow */
  } /* end of loop over icol */
  
  if (fdyn->conte!=0)  
  {
    ccb  = ccon*beta*sign;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau_mu * nue * w * div(u_old) * delta(bub)   d_omega
  /
 *----------------------------------------------------------------------*/
    for (icol=0;icol<smiel;icol++)
    {
      aux = (smderxy2[0][icol]+smderxy2[1][icol]+smderxy2[2][icol])*divv*ccb;
      for (irow=0;irow<smiel;irow++)
      {
  	smestif[irow][icol] -= smfunct[irow]*aux;
      } /* end of loop over irow */
    } /* end of loop over icol */
  } 

/*----------------------------------------------------------------------*
   Calculate viscous stabilization part for higher-order elements:
 *----------------------------------------------------------------------*/   
  ccon = con * visc * visc * sign;
/*----------------------------------------------------------------------*
    /
   |  +/- tau  * nue * delta(w) * nue * delta(bub) d_omega
  /
 *----------------------------------------------------------------------*/   
  for (icol=0;icol<smiel;icol++)
  {
    auxc = (smderxy2[0][icol]+smderxy2[1][icol]+smderxy2[2][icol])*ccon;
    for (irow=0;irow<smiel;irow++)
    {
      aux = (smderxy2[0][irow]+smderxy2[1][irow]+smderxy2[2][irow])*auxc;
      smestif[irow][icol] += aux;
    } /* end of loop over irow */
  } /* end of loop over icol */

/*----------------------------------------------------------------------*
    /
   |  -/+ tau  * nue * delta(w) * u_old * grad(bub) d_omega
  /
 *----------------------------------------------------------------------*/
  ccon = con * visc * sign;
  
  for (icol=0;icol<smiel;icol++)
  {
    auxc = (velint[0]*smderxy[0][icol]+velint[1]*smderxy[1][icol]\
           +velint[2]*smderxy[2][icol])*ccon;
    for (irow=0;irow<smiel;irow++)
    {
      aux = (smderxy2[0][irow]+smderxy2[1][irow]+smderxy2[2][irow])*auxc;
      smestif[irow][icol] -= aux;
    } /* end of loop over irow */
  } /* end of loop over icol */

  if (fdyn->conte!=0)  
  {
    ccb  = ccon*beta;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau * nue * delta(w) * bub * div(u_old) d_omega
  /
 *----------------------------------------------------------------------*/
    for (icol=0;icol<smiel;icol++)
    {
      auxc = smfunct[icol]*divv*ccb;
      for (irow=0;irow<smiel;irow++)
      {
    	aux = (smderxy2[0][irow]+smderxy2[1][irow]+smderxy2[2][irow])*auxc;
    	smestif[irow][icol] -= aux;
      } /* end of loop over irow */
    } /* end of loop over icol */
  }
} /* end of stabilization for higher order elements */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calstabsmk */

/*!---------------------------------------------------------------------                                         
\brief evaluate stabilization part of submesh mass matrix SMM for fluid3

<pre>                                                       gravem 07/03

In this routine, the stabilization part of the submesh mass matrix SMM 
is calculated.

</pre>
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smemass   DOUBLE	   (i/o)  submesh element mass matrix
\param  *velint    DOUBLE	   (i)    velocity at int. point
\param **vderxy    DOUBLE	   (i)    global vel. deriv. at int. p.
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun. 
\param **smderxy2  DOUBLE	   (i)    sm 2nd global deriv. of sh. fun. 
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc      DOUBLE	   (i)    fluid viscosity	     
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\param   ihoelsm   INT  	   (i)	  flag for higher-order elements
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabsmm(FLUID_DYN_ML    *mlvar, 
		   DOUBLE         **smemass,  
		   DOUBLE          *velint, 
		   DOUBLE         **vderxy, 
		   DOUBLE          *smfunct,  
		   DOUBLE         **smderxy,  
		   DOUBLE         **smderxy2, 
		   DOUBLE           fac,    
		   DOUBLE           visc,   
		   INT              smiel,    
                   INT              ihoelsm)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
/*----------------------------------------------------------------------*/
INT    irow,icol;
DOUBLE tau;
DOUBLE con,ccon,beta,divv,cb;
DOUBLE aux,auxc;
DOUBLE sign;

#ifdef DEBUG 
dstrc_enter("f3_calstabsmm");
#endif

fdyn = alldyn[genprob.numff].fdyn;
/*---------------------------------------- set stabilization parameter */
tau = mlvar->smtau;

switch (fdyn->conte) /* determine parameter beta of convective term */
{
case 0: 
  beta = ZERO;
break;
case 1: 
  beta = ONE;
  divv= vderxy[0][0]+vderxy[1][1]+vderxy[2][2]; 
break;
case 2: 
  beta = ONE/TWO;
  divv= vderxy[0][0]+vderxy[1][1]+vderxy[2][2]; 
break;
default:
   dserror("unknown form of convective term");
} /* end switch (fdyn->conte) */

if (mlvar->smstado<0) sign = -ONE;  /* USFEM */
else                  sign = ONE;   /* GLS- */

/*----------------------------------------------------------------------*
   Calculate temporal stabilization part:
 *----------------------------------------------------------------------*/
if (ABS(mlvar->smstado)<3)
{ 
  if (mlvar->smstado==-1) con = fac*tau;
  else                    con = fac*tau*sign; 
/*---------------------------------------------------------------------*
    /
   |  +/- tau * (1/(theta*dt)) * w * bub d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<smiel;icol++)
  {
    auxc = smfunct[icol]*con*(ONE/fdyn->thsl);
    for (irow=0;irow<smiel;irow++)
    {
      aux = smfunct[irow]*auxc;
      smemass[irow][icol] += aux;
    } /* end loop over irow */
  } /* end loop over icol */

/*---------------------------------------------------------------------*
    /
   |  +/- tau * w * u_old * grad(bub) d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<smiel;icol++)
  {
    auxc = (velint[0]*smderxy[0][icol]+velint[1]*smderxy[1][icol]\
           +velint[2]*smderxy[2][icol])*con;
    for (irow=0;irow<smiel;irow++)
    {
      aux = smfunct[irow]*auxc;
      smemass[irow][icol] += aux;
    } /* end loop over irow */
  } /* end loop over icol */

  if (fdyn->conte!=0)  
  {
    cb = con*beta;
/*---------------------------------------------------------------------*
    /
   |  +/- tau * w * bub * div(u_old) d_omega
  /
 *----------------------------------------------------------------------*/
    for (icol=0;icol<smiel;icol++)
    {
      auxc = smfunct[icol]*divv*cb;
      for (irow=0;irow<smiel;irow++)
      {
        aux = smfunct[irow]*auxc;
        smemass[irow][icol] += aux;
      } /* end loop over irow */
    } /* end loop over icol */
  }  

/*----------------------------------------------------------------------*
   Calculate temporal stabilization part for higher-order elements:
 *----------------------------------------------------------------------*/
  if (ihoelsm!=0)
  {
    ccon = con * visc;
/*----------------------------------------------------------------------*
    /
   |  -/+ tau * w * nue * delta(bub)  d_omega
  /
 *----------------------------------------------------------------------*/
    for (icol=0;icol<smiel;icol++)
    {	 
      aux = (smderxy2[0][icol]+smderxy2[1][icol]+smderxy2[2][icol])*ccon;
      for(irow=0;irow<smiel;irow++)
      {
        smemass[irow][icol] -= smfunct[irow]*aux;
      } /* end loop over irow */
    } /* end loop over icol */
  } /* end of temporal stabilization for higher-order elements */
} /* end of temporal stabilization */   

/*----------------------------------------------------------------------*
   Calculate convective stabilization part:
 *----------------------------------------------------------------------*/
con = fac*tau;
/*----------------------------------------------------------------------*
    /
   |  tau * u_old * grad(w) * bub d_omega
  /
 *----------------------------------------------------------------------*/
for (icol=0;icol<smiel;icol++)
{
  auxc = smfunct[icol]*con;
  for (irow=0;irow<smiel;irow++)
  {
    aux = (velint[0]*smderxy[0][irow]+velint[1]*smderxy[1][irow]\
          +velint[2]*smderxy[2][irow])*auxc;
    smemass[irow][icol]	+= aux;
  } /* end loop over irow */
} /* end loop over icol */
    
if (fdyn->conte!=0)  
{
  cb = con*beta*sign;
/*----------------------------------------------------------------------*
    /
   |  +/- beta * tau * w * div(u_old) * bub  d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<smiel;icol++)
  {
    auxc = smfunct[icol]*divv*cb;
    for (irow=0;irow<smiel;irow++)
    {
      aux = smfunct[irow]*auxc;
      smemass[irow][icol] += aux;
    } /* end of loop over irow */
  } /* end of loop over icol */
}  

/*----------------------------------------------------------------------*
   Calculate viscous stabilization part for higher-order elements:
 *----------------------------------------------------------------------*/
if (ihoelsm!=0)
{
  ccon = con * visc * sign;
/*----------------------------------------------------------------------*
    /
   |  -/+ tau * nue * delta(w) * bub  d_omega
  /
 *----------------------------------------------------------------------*/
  for (icol=0;icol<smiel;icol++)
  {	 
    auxc = smfunct[icol]*ccon;
    for(irow=0;irow<smiel;irow++)
    {
      aux = (smderxy2[0][irow]+smderxy2[1][irow]+smderxy2[2][irow])*auxc;
      smemass[irow][icol] -= aux;
    } /* end loop over irow */
  } /* end loop over icol */
} /* end of viscous stabilization for higher-order elements */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calstabsmm */

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
void f3_lscalstabkvv(ELEMENT         *ele,
		   DOUBLE	  **estif,
		   DOUBLE	   *velint,
		   DOUBLE	  **vderxy,
		   DOUBLE	   *funct,
		   DOUBLE	  **derxy,
		   DOUBLE	  **derxy2,
		   DOUBLE	    fac,
		   DOUBLE	    visc,
		   INT  	    iel,	
                   INT  	    ihoel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |  
/*----------------------------------------------------------------------*/
INT    irow,icol,irn,icn,ird;
DOUBLE taumu;
DOUBLE taump;
DOUBLE tauc;
DOUBLE con,ccon,beta,divv,cb,cbb,ccb;
DOUBLE aux,auxr,auxc;
DOUBLE signvi,signre;
STAB_PAR_GLS *gls;	/* pointer to GLS stabilisation parameters	*/

#ifdef DEBUG 
dstrc_enter("f3_lscalstabkvv");
#endif
/*---------------------------------------------------------- initialise */
fdyn = alldyn[genprob.numff].fdyn;

gls    = ele->e.f3->stabi.gls;

if (ele->e.f3->stab_type != stab_gls) 
   dserror("routine with no or wrong stabilisation called");
   
/*---------------------------------------- set stabilisation parameter */
taumu = fdyn->tau[0];
taump = fdyn->tau[1];
tauc  = fdyn->tau[2];

switch (fdyn->conte) /* determine parameter beta of convective term */
{
case 0: 
  beta = ZERO;
break;
case 1: 
  beta = ONE;
  divv= vderxy[0][0]+vderxy[1][1]+vderxy[2][2]; 
break;
case 2: 
  beta = ONE/TWO;
  divv= vderxy[0][0]+vderxy[1][1]+vderxy[2][2]; 
break;
default:
   dserror("unknown form of convective term");
} /* end switch (fdyn->conte) */

switch (gls->ivisc) /* determine type of stabilization */
{
case -1: /* USFEM (but no viscous stabilization) */
   signre = -ONE;
break;
case 0: /* GLS (but no viscous stabilization) */
   signre = ZERO;
break;
case 1: /* GLS- */
   signre = ONE;
   signvi = ONE;
break;
case 2: /* GLS+ */
   signre = ONE;
   signvi = -ONE;
break;
case 3: /* USFEM */
   signre = -ONE;
   signvi = -ONE;
break;
default:
   dserror("unknown type of stabilization");
} /* end switch (ele->e.f3->ivisc) */

/*----------------------------------------------------------------------*
   Calculate continuity stabilisation part:
    /
   |  tau_c * div(v) * div(u)   d_omega
  /
 *----------------------------------------------------------------------*/
if (gls->icont!=0)
{
   con = fac*tauc;

   icol=0;
   for(icn=0;icn<iel;icn++)
   {
      irow=0;
      for(irn=0;irn<iel;irn++)
      {
         estif[irow][icol]     += derxy[0][icn]*derxy[0][irn]*con;
         estif[irow][icol+1]   += derxy[1][icn]*derxy[0][irn]*con;
         estif[irow][icol+2]   += derxy[2][icn]*derxy[0][irn]*con;

         estif[irow+1][icol]   += derxy[0][icn]*derxy[1][irn]*con;
         estif[irow+1][icol+1] += derxy[1][icn]*derxy[1][irn]*con;
         estif[irow+1][icol+2] += derxy[2][icn]*derxy[1][irn]*con;
         estif[irow+2][icol]   += derxy[0][icn]*derxy[2][irn]*con;
         estif[irow+2][icol+1] += derxy[1][icn]*derxy[2][irn]*con;
         estif[irow+2][icol+2] += derxy[2][icn]*derxy[2][irn]*con;
         irow += 3;
      } /* end of loop over irn */
      icol += 3;
   } /* end of loop over icn */
} /* endif (ele->e.f3->icont!=0) */

con = fac*taumu;

if (gls->iadvec!=0)
{
/*----------------------------------------------------------------------*
   Calculate convection stabilisation part Nc(u):
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
    /
   |  tau_mu * u_old * grad(v) * u_old * grad(u)   d_omega
  /
 *----------------------------------------------------------------------*/
    icol=0;
    for (icn=0;icn<iel;icn++)
    {
      auxc =  (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn] \
             + velint[2]*derxy[2][icn])*con;
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

    if (fdyn->conte!=0)  
    {
      cb  = con*beta;
      cbb = cb*beta*signre;
/*----------------------------------------------------------------------*
    /
   |  beta * tau_mu * u_old * grad(v) * u * div(u_old)   d_omega
  /
 *----------------------------------------------------------------------*/
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        auxc = funct[icn]*divv*cb;
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
/*----------------------------------------------------------------------*
    /
   |  +/- beta * tau_mu * v * div(u_old) * u_old * grad(u)   d_omega
  /
 *----------------------------------------------------------------------*/
      cb = cb*signre;
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        aux = (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn] \
	     + velint[2]*derxy[2][irn])*divv*cb;
        irow=0;
        for (irn=0;irn<iel;irn++)
        {
          estif[irow][icol]     += aux*funct[irn];
          estif[irow+1][icol+1] += aux*funct[irn];
          estif[irow+2][icol+2] += aux*funct[irn];
          irow += 3;
        } /* end of loop over irn */
        icol += 3;
      } /* end of loop over icn */
/*----------------------------------------------------------------------*
    /
   |  +/- beta * beta * tau_mu * v * div(u_old) * u * div(u_old)   d_omega
  /
 *----------------------------------------------------------------------*/
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        aux = funct[icn]*divv*divv*cbb;
        irow=0;
        for (irn=0;irn<iel;irn++)
        {
          estif[irow][icol]     += aux*funct[irn];
          estif[irow+1][icol+1] += aux*funct[irn];
          estif[irow+2][icol+2] += aux*funct[irn];
          irow += 3;
        } /* end of loop over irn */
        icol += 3;
      } /* end of loop over icn */
    }

/*----------------------------------------------------------------------*
   Calculate convection stabilisation part Nr(u):
 *----------------------------------------------------------------------*/
  if (fdyn->nir!=0)  /* evaluate for Newton iteraton */
  {
/*----------------------------------------------------------------------*
    /
   |  tau_mu * u_old * grad(v) * u * grad(u_old)   d_omega
  /
 *----------------------------------------------------------------------*/
    icol=0;
    for (icn=0;icn<iel;icn++)
    {
      auxc = funct[icn]*con;
      irow=0;
      for (irn=0;irn<iel;irn++)
      {
        aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn] \
             + velint[2]*derxy[2][irn])*auxc;

        estif[irow][icol]     += aux*vderxy[0][0];
        estif[irow][icol+1]   += aux*vderxy[0][1];
        estif[irow][icol+2]   += aux*vderxy[0][2];
        
        estif[irow+1][icol]   += aux*vderxy[1][0];
        estif[irow+1][icol+1] += aux*vderxy[1][1];
        estif[irow+1][icol+2] += aux*vderxy[1][2];

        estif[irow+2][icol]   += aux*vderxy[2][0];
        estif[irow+2][icol+1] += aux*vderxy[2][1];	
        estif[irow+2][icol+2] += aux*vderxy[2][2];   
        irow += 3;
      } /* end of loop over irn */
      icol += 3;
    } /* end of loop over icn */

    if (fdyn->conte!=0)  
    {
      cb  = con*beta;
      cbb = cb*beta*signre;
/*----------------------------------------------------------------------*
    /
   |  beta * tau_mu * u_old * grad(v) * u_old * div(u)   d_omega
  /
 *----------------------------------------------------------------------*/
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        irow=0;
        for (irn=0;irn<iel;irn++)
        {
          aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn] \
               + velint[2]*derxy[2][irn])*cb;
          estif[irow][icol]     += aux*velint[0]*derxy[0][icn];
          estif[irow][icol+1]   += aux*velint[0]*derxy[1][icn];
          estif[irow][icol+2]   += aux*velint[0]*derxy[2][icn];

          estif[irow+1][icol]   += aux*velint[1]*derxy[0][icn];
          estif[irow+1][icol+1] += aux*velint[1]*derxy[1][icn];
          estif[irow+1][icol+2] += aux*velint[1]*derxy[2][icn];

          estif[irow+2][icol]   += aux*velint[2]*derxy[0][icn];
          estif[irow+2][icol+1] += aux*velint[2]*derxy[1][icn];
          estif[irow+2][icol+2] += aux*velint[2]*derxy[2][icn];
          irow += 3;
        } /* end of loop over irn */
        icol += 3;
      } /* end of loop over icn */
/*----------------------------------------------------------------------*
    /
   |  +/- beta * tau_mu * v * div(u_old) * u * grad(u_old)   d_omega
  /
 *----------------------------------------------------------------------*/
      cb = cb*signre;
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        auxc = funct[icn]*divv*cb;
        irow=0;
        for (irn=0;irn<iel;irn++)
        {
          aux = funct[irn]*auxc;
          estif[irow][icol]     += aux*vderxy[0][0];
          estif[irow][icol+1]   += aux*vderxy[0][1];
          estif[irow][icol+2]   += aux*vderxy[0][2];

          estif[irow+1][icol]   += aux*vderxy[1][0];
          estif[irow+1][icol+1] += aux*vderxy[1][1];
          estif[irow+1][icol+2] += aux*vderxy[1][2];

          estif[irow+2][icol]   += aux*vderxy[2][0];
          estif[irow+2][icol+1] += aux*vderxy[2][1];
          estif[irow+2][icol+2] += aux*vderxy[2][2];
          irow += 3;
        } /* end of loop over irn */
        icol += 3;
      } /* end of loop over icn */
/*----------------------------------------------------------------------*
    /
   |  +/- beta * beta * tau_mu * v * div(u_old) * u_old * div(u)   d_omega
  /
 *----------------------------------------------------------------------*/
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        irow=0;
        for (irn=0;irn<iel;irn++)
        {
          aux = funct[irn]*divv*cbb;
          estif[irow][icol]     += aux*velint[0]*derxy[0][icn];
          estif[irow][icol+1]   += aux*velint[0]*derxy[1][icn];
          estif[irow][icol+2]   += aux*velint[0]*derxy[2][icn];

          estif[irow+1][icol]   += aux*velint[1]*derxy[0][icn];
          estif[irow+1][icol+1] += aux*velint[1]*derxy[1][icn];
          estif[irow+1][icol+2] += aux*velint[1]*derxy[2][icn];

          estif[irow+2][icol]   += aux*velint[2]*derxy[0][icn];
          estif[irow+2][icol+1] += aux*velint[2]*derxy[1][icn];
          estif[irow+2][icol+2] += aux*velint[2]*derxy[2][icn];
          irow += 3;
        } /* end of loop over irn */
        icol += 3;
      } /* end of loop over icn */
    }
  } /* endif (fdyn->nir!=0) */

/*----------------------------------------------------------------------*
   Calculate convection stabilisation part for higher order elements:
 *----------------------------------------------------------------------*/
  if (ihoel!=0)
  {
    ccon = con*visc;
    
    if (fdyn->vite==0)
    {
/*----------------------------------------------------------------------*
    /
   |  - tau_mu * nue * u_old * grad(v) * delta(u)   d_omega
  /
 *----------------------------------------------------------------------*/
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        irow=0;
        auxc = (derxy2[0][icn]+derxy2[1][icn]+derxy2[2][icn])*ccon;
        for (irn=0;irn<iel;irn++)
        {
          aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn] \
               + velint[2]*derxy[2][irn])*auxc;
          estif[irow][icol]     -= aux;
          estif[irow+1][icol+1] -= aux;
          estif[irow+2][icol+2] -= aux;
          irow += 3;
        } /* end of loop over irn */
        icol += 3;
      } /* end of loop over icn */
      
      if (fdyn->conte!=0)  
      {
        ccb  = ccon*beta*signre;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau_mu * nue * v * div(u_old) * delta(u)   d_omega
  /
 *----------------------------------------------------------------------*/
        icol=0;
        for (icn=0;icn<iel;icn++)
        {
          auxc = (derxy2[0][icn]+derxy2[1][icn]+derxy2[2][icn])*divv*ccb;
          irow=0;
          for (irn=0;irn<iel;irn++)
          {
            aux = funct[irn]*auxc;
            estif[irow][icol]     -= aux;
            estif[irow+1][icol+1] -= aux;
            estif[irow+2][icol+2] -= aux;
            irow += 3;
          } /* end of loop over irn */
          icol += 3;
        } /* end of loop over icn */
      }	
    }
    else
    { 
/*----------------------------------------------------------------------*
    /
   |  - tau_mu * 2 * nue * u_old * grad(v) * div(eps(u))   d_omega
  /
 *----------------------------------------------------------------------*/
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        irow=0;
	auxc = derxy2[0][icn] + derxy2[1][icn] + derxy2[2][icn];
	for (irn=0;irn<iel;irn++)
	{
	  aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn] \
	       + velint[2]*derxy[2][irn])*ccon;

	  estif[irow][icol]	-= aux*(derxy2[0][icn] + auxc);
	  estif[irow][icol+1]	-= aux* derxy2[3][icn];
	  estif[irow][icol+2]	-= aux* derxy2[4][icn];

	  estif[irow+1][icol]	-= aux* derxy2[3][icn];
	  estif[irow+1][icol+1] -= aux*(derxy2[1][icn] + auxc);
	  estif[irow+1][icol+2] -= aux* derxy2[5][icn];

	  estif[irow+2][icol]	-= aux* derxy2[4][icn];
	  estif[irow+2][icol+1] -= aux* derxy2[5][icn];
	  estif[irow+2][icol+2] -= aux*(derxy2[2][icn] + auxc);
	  irow += 3;
	} /* end of loop over irn */
        icol += 3;
      } /* end of loop over icn */
      
      if (fdyn->conte!=0)  
      {
        ccb  = ccon*beta*signre;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau_mu * 2 * nue * v * div(u_old) * div(eps(u))   d_omega
  /
 *----------------------------------------------------------------------*/
        icol=0;
        for (icn=0;icn<iel;icn++)
        {
          auxc = derxy2[0][icn] + derxy2[1][icn] + derxy2[2][icn];
          irow=0;
          for (irn=0;irn<iel;irn++)
          {
            aux = funct[irn]*divv*ccb;
            estif[irow][icol]     -= aux*(derxy2[0][icn] + auxc);
            estif[irow][icol+1]   -= aux* derxy2[3][icn];
            estif[irow][icol+2]   -= aux* derxy2[4][icn];

            estif[irow+1][icol]   -= aux* derxy2[3][icn];
            estif[irow+1][icol+1] -= aux*(derxy2[1][icn] + auxc);
            estif[irow+1][icol+2] -= aux* derxy2[5][icn];

            estif[irow+2][icol]   -= aux* derxy2[4][icn];
            estif[irow+2][icol+1] -= aux* derxy2[5][icn];
            estif[irow+2][icol+2] -= aux*(derxy2[2][icn] + auxc);
            irow += 3;
          } /* end of loop over irn */
          icol += 3;
        } /* end of loop over icn */
      }	
    }      
  } /* end of stabilisation for higher order elements */
} /* end of convection stabilisation */

/*----------------------------------------------------------------------*
   Calculate viscous stabilisation part for higher order elements:
 *----------------------------------------------------------------------*/   
if (ihoel!=0 && gls->ivisc>0)
{   
  ccon = fac * taump * visc * visc * signvi;
  
  if (fdyn->vite==0)
  {
/*----------------------------------------------------------------------*
    /
   |  +/- tau_mp  * nue**2 * delta(v) * delta(u) d_omega
  /
 *----------------------------------------------------------------------*/   
    icol=0;
    for (icn=0;icn<iel;icn++)
    {
      irow=0;
      auxc = (derxy2[0][icn] + derxy2[1][icn] + derxy2[2][icn])*ccon;
      for (irn=0;irn<iel;irn++)
      {
        aux = (derxy2[0][irn] + derxy2[1][irn] + derxy2[2][irn])*auxc;
	estif[irow][icol]     += aux;
	estif[irow+1][icol+1] += aux;
	estif[irow+2][icol+2] += aux;
	irow += 3;				     
      } /* end of loop over irn */
      icol += 3;
    } /* end of loop over icn */

/*----------------------------------------------------------------------*
   Calculate viscous stabilisation part Nc(u) for higher order elements:
    /
   |  -/+ tau_mp  * nue * delta(v) * u_old * grad(u) d_omega
  /
 *----------------------------------------------------------------------*/
    ccon = fac * taump * visc * signvi;
   
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        irow=0;
	auxc = (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn] \
	      + velint[2]*derxy[2][icn])*ccon;
	for (irn=0;irn<iel;irn++)
	{
	  aux = (derxy2[0][irn] + derxy2[1][irn] + derxy2[2][irn])*auxc;
	  estif[irow][icol]	-= aux;
	  estif[irow+1][icol+1] -= aux;
	  estif[irow+2][icol+2] -= aux;
	  irow += 3;
	} /* end of loop over irn */
	icol += 3;
      } /* end of loop over icn */

      if (fdyn->conte!=0)  
      {
        ccb  = ccon*beta;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau_mp  * nue * delta(v) * u * div(u_old) d_omega
  /
 *----------------------------------------------------------------------*/
        icol=0;
        for (icn=0;icn<iel;icn++)
        {
          irow=0;
	  auxc = funct[icn]*divv*ccb;
	  for (irn=0;irn<iel;irn++)
	  {
	    aux = (derxy2[0][irn] + derxy2[1][irn] + derxy2[2][irn])*auxc;
	    estif[irow][icol]	  -= aux;
	    estif[irow+1][icol+1] -= aux;
	    estif[irow+2][icol+2] -= aux;
	    irow += 3;
	  } /* end of loop over irn */
	  icol += 3;
        } /* end of loop over icn */
      }
   
/*----------------------------------------------------------------------*
   Calculate viscous stabilisation part Nr(u) for higher order elements:
    /
   |  -/+ tau_mp  * nue * delta(v) * u * grad(u_old) d_omega
  /
 *----------------------------------------------------------------------*/   
    if (fdyn->nir!=0)  /* evaluate for Newton iteraton */
    {
      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        irow=0;
	auxc = funct[icn]*ccon;
	for (irn=0;irn<iel;irn++)
	{
	  aux = (derxy2[0][irn] + derxy2[1][irn] + derxy2[2][irn])*auxc;
	  estif[irow][icol]	-= vderxy[0][0]*aux;
	  estif[irow][icol+1]	-= vderxy[0][1]*aux;
	  estif[irow][icol+2]	-= vderxy[0][2]*aux;

	  estif[irow+1][icol]	-= vderxy[1][0]*aux;
	  estif[irow+1][icol+1] -= vderxy[1][1]*aux;
	  estif[irow+1][icol+2] -= vderxy[1][2]*aux;

	  estif[irow+2][icol]	-= vderxy[2][0]*aux;
	  estif[irow+2][icol+1] -= vderxy[2][1]*aux;
	  estif[irow+2][icol+2] -= vderxy[2][2]*aux;
	  irow += 3;
	} /* end of loop over irn */
	icol += 3;
      } /* end of loop over icn */

      if (fdyn->conte!=0)  
      {
        ccb  = ccon*beta;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau_mp  * nue * delta(v) * u_old * div(u) d_omega
  /
 *----------------------------------------------------------------------*/
        icol=0;
        for (icn=0;icn<iel;icn++)
        {
          irow=0;
	  for (irn=0;irn<iel;irn++)
	  {
	    aux = (derxy2[0][irn] + derxy2[1][irn] + derxy2[2][irn])*ccb;
	    estif[irow][icol]	  -= velint[0]*derxy[0][icn]*aux;
	    estif[irow][icol+1]	  -= velint[0]*derxy[1][icn]*aux;
	    estif[irow][icol+2]	  -= velint[0]*derxy[2][icn]*aux;

	    estif[irow+1][icol]	  -= velint[1]*derxy[0][icn]*aux;
	    estif[irow+1][icol+1] -= velint[1]*derxy[1][icn]*aux;
	    estif[irow+1][icol+2] -= velint[1]*derxy[2][icn]*aux;

	    estif[irow+2][icol]	  -= velint[2]*derxy[0][icn]*aux;
	    estif[irow+2][icol+1] -= velint[2]*derxy[1][icn]*aux;
	    estif[irow+2][icol+2] -= velint[2]*derxy[2][icn]*aux;
	    irow += 3;
	  } /* end of loop over irn */
	  icol += 3;
        } /* end of loop over icn */
      }
    } /* endif (fdyn->nir!=0) */
  }
  else
  {
/*----------------------------------------------------------------------*
    /
   |  +/- tau_mp  * 4 * nue**2 * div(eps(v)) * div(eps(u)) d_omega
  /
 *----------------------------------------------------------------------*/   
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
	        		 + derxy2[4][icn]*derxy2[4][irn] )*ccon;  
	estif[irow][icol+1]   += ( (auxc + derxy2[0][icn])*derxy2[3][irn]  
				 + derxy2[3][icn]*(auxr + derxy2[1][irn])  
	        		 + derxy2[4][icn]*derxy2[5][irn] )*ccon;  
	estif[irow][icol+2]   += ( (auxc + derxy2[0][icn])*derxy2[4][irn]  
				 + derxy2[3][icn]*derxy2[5][irn]  
	        		 + derxy2[4][icn]*(auxr + derxy2[2][irn]) )*ccon;

	estif[irow+1][icol]   += ( derxy2[3][icn]*(auxr + derxy2[0][irn]) 	     
				 + (auxc + derxy2[1][icn])*derxy2[3][irn]   
	        		 + derxy2[5][icn]*derxy2[4][irn] )*ccon;
	estif[irow+1][icol+1] += ( derxy2[3][icn]*derxy2[3][irn]  
				 + (auxc + derxy2[1][icn])*(auxr + derxy2[1][irn])  
	        		 + derxy2[5][icn]*derxy2[5][irn] )*ccon; 
	estif[irow+1][icol+2] += ( derxy2[3][icn]*derxy2[4][irn]  
				 + (auxc + derxy2[1][icn])*derxy2[5][irn]  
	        		 + derxy2[5][icn]*(auxr + derxy2[2][irn]) )*ccon; 

	estif[irow+2][icol]   += ( derxy2[4][icn]*(auxr + derxy2[0][irn])   
				 + derxy2[5][icn]*derxy2[3][irn]    
	        		 + (auxc + derxy2[2][icn])*derxy2[4][irn] )*ccon; 
	estif[irow+2][icol+1] += ( derxy2[4][icn]*derxy2[3][irn]  
				 + derxy2[5][icn]*(auxr + derxy2[1][irn])  
	        		 + (auxc + derxy2[2][icn])*derxy2[5][irn] )*ccon;
	estif[irow+2][icol+2] += ( derxy2[4][icn]*derxy2[4][irn]  
				 + derxy2[5][icn]*derxy2[5][irn]  
	        		 + (auxc + derxy2[2][irn])*(auxr + derxy2[2][irn]) )*ccon;		 
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
    ccon = fac * taump * visc * signvi;

      icol=0;
      for (icn=0;icn<iel;icn++)
      {
        irow=0;
	auxc = (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn] \
	      + velint[2]*derxy[2][icn])*ccon;
	for (irn=0;irn<iel;irn++)
	{
	  auxr = derxy2[0][irn] + derxy2[1][irn] + derxy2[2][irn];
	  estif[irow][icol]	-= (derxy2[0][irn]+auxr)*auxc;
	  estif[irow+1][icol]	-=  derxy2[3][irn]*auxc;
	  estif[irow+2][icol]	-=  derxy2[4][irn]*auxc;

	  estif[irow][icol+1]	-=  derxy2[3][irn]*auxc;
	  estif[irow+1][icol+1] -= (derxy2[1][irn]+auxr)*auxc;
	  estif[irow+2][icol+1] -=  derxy2[5][irn]*auxc;

	  estif[irow][icol+2]	-=  derxy2[4][irn]*auxc;
	  estif[irow+1][icol+2] -=  derxy2[5][irn]*auxc;
	  estif[irow+2][icol+2] -= (derxy2[2][irn]+auxr)*auxc;
	  irow += 3;
	} /* end of loop over irn */
	icol += 3;
      } /* end of loop over icn */

      if (fdyn->conte!=0)  
      {
        ccb  = ccon*beta;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau_mp  * 2 * nue * div(eps(v)) * u * div(u_old) d_omega
  /
 *----------------------------------------------------------------------*/
        icol=0;
        for (icn=0;icn<iel;icn++)
        {
          irow=0;
	  auxc = funct[icn]*divv*ccb;
	  for (irn=0;irn<iel;irn++)
	  {
	    auxr = derxy2[0][irn] + derxy2[1][irn] + derxy2[2][irn];
	    estif[irow][icol]	  -= (derxy2[0][irn]+auxr)*auxc;
	    estif[irow][icol+1]   -=  derxy2[3][irn]*auxc;
	    estif[irow][icol+2]   -=  derxy2[4][irn]*auxc;

	    estif[irow+1][icol]   -=  derxy2[3][irn]*auxc;
	    estif[irow+1][icol+1] -= (derxy2[1][irn]+auxr)*auxc;
	    estif[irow+1][icol+2] -=  derxy2[5][irn]*auxc;

	    estif[irow+2][icol]   -=  derxy2[4][irn]*auxc;
	    estif[irow+2][icol+1] -=  derxy2[5][irn]*auxc;
	    estif[irow+2][icol+2] -= (derxy2[2][irn]+auxr)*auxc;
	    irow += 3;
	  } /* end of loop over irn */
	  icol += 3;
        } /* end of loop over icn */
      }
    
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
	auxc = funct[icn]*ccon;
	for (irn=0;irn<iel;irn++)
	{
	  auxr = derxy2[0][irn]+derxy2[1][irn]+derxy2[2][irn];
	  estif[irow][icol]	-=  (derxy2[0][irn]*vderxy[0][0]   
	        		   + derxy2[3][irn]*vderxy[1][0]   
	        		   + derxy2[4][irn]*vderxy[2][0]   
	  			   + auxr*vderxy[0][0])*auxc	   ;
	  estif[irow][icol+1]	-=  (derxy2[0][irn]*vderxy[0][1]   
	        		   + derxy2[3][irn]*vderxy[1][1]   
	        		   + derxy2[4][irn]*vderxy[2][1]   
	  			   + auxr*vderxy[0][1])*auxc	   ;
	  estif[irow][icol+2]	-=  (derxy2[0][irn]*vderxy[0][2]   
	        		   + derxy2[3][irn]*vderxy[1][2]   
	        		   + derxy2[4][irn]*vderxy[2][2]   
	  			   + auxr*vderxy[0][2])*auxc	   ;

	  estif[irow+1][icol]	-=  (derxy2[3][irn]*vderxy[0][0]   
	        		   + derxy2[1][irn]*vderxy[1][0]   
	        		   + derxy2[5][irn]*vderxy[2][0]   
	  			   + auxr*vderxy[1][0])*auxc	   ;
	  estif[irow+1][icol+1] -=  (derxy2[3][irn]*vderxy[0][1]   
	        		   + derxy2[1][irn]*vderxy[1][1]   
	        		   + derxy2[5][irn]*vderxy[2][1]   
	  			   + auxr*vderxy[1][1])*auxc	   ;
	  estif[irow+1][icol+2] -=  (derxy2[3][irn]*vderxy[0][2]   
	        		   + derxy2[1][irn]*vderxy[1][2]   
	        		   + derxy2[5][irn]*vderxy[2][2]   
	  			   + auxr*vderxy[1][2])*auxc	   ;

	  estif[irow+2][icol]	-=  (derxy2[4][irn]*vderxy[0][0]   
	        		   + derxy2[5][irn]*vderxy[1][0]   
	        		   + derxy2[2][irn]*vderxy[2][0]   
	  			   + auxr*vderxy[2][0])*auxc	   ;
	  estif[irow+2][icol+1] -=  (derxy2[4][irn]*vderxy[0][1]   
	        		   + derxy2[5][irn]*vderxy[1][1]   
	        		   + derxy2[2][irn]*vderxy[2][1]   
	  			   + auxr*vderxy[2][1])*auxc	   ;
	  estif[irow+2][icol+2] -=  (derxy2[4][irn]*vderxy[0][2]   
	        		   + derxy2[5][irn]*vderxy[1][2]   
	        		   + derxy2[2][irn]*vderxy[2][2]   
	  			   + auxr*vderxy[2][2])*auxc	   ;
	  irow += 3;
	} /* end of loop over irn */
	icol += 3;
      } /* end of loop over icn */

      if (fdyn->conte!=0)  
      {
        ccb  = ccon*beta;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau_mp  * 2 * nue * div(eps(v)) * u_old * div(u) d_omega
  /
 *----------------------------------------------------------------------*/
        icol=0;
        for (icn=0;icn<iel;icn++)
        {
          irow=0;
	  for (irn=0;irn<iel;irn++)
	  {
	    aux = derxy2[0][irn] + derxy2[1][irn] + derxy2[2][irn];
	    estif[irow][icol]	  -=  ((derxy2[0][irn]+aux)*velint[0] 
	        	  	      + derxy2[3][irn]*velint[1]       
	        	  	      + derxy2[4][irn]*velint[2])      
				      *derxy[0][icn]*ccb              ;
	    estif[irow][icol+1]	  -=  ((derxy2[0][irn]+aux)*velint[0] 
	        	  	      + derxy2[3][irn]*velint[1]       
	        	  	      + derxy2[4][irn]*velint[2])      
				      *derxy[1][icn]*ccb              ;
	    estif[irow][icol+2]	  -=  ((derxy2[0][irn]+aux)*velint[0] 
	        	  	      + derxy2[3][irn]*velint[1]       
	        	  	      + derxy2[4][irn]*velint[2])      
				      *derxy[2][icn]*ccb              ;

	    estif[irow+1][icol]	  -=  ((derxy2[1][irn]+aux)*velint[1] 
	        	  	      + derxy2[3][irn]*velint[0]       
	        	  	      + derxy2[5][irn]*velint[2])      
				      *derxy[0][icn]*ccb              ;
	    estif[irow+1][icol+1] -=  ((derxy2[1][irn]+aux)*velint[1] 
	        	  	      + derxy2[3][irn]*velint[0]       
	        	  	      + derxy2[5][irn]*velint[2])      
				      *derxy[1][icn]*ccb              ;
	    estif[irow+1][icol+2] -=  ((derxy2[1][irn]+aux)*velint[1] 
	        	  	      + derxy2[3][irn]*velint[0]       
	        	  	      + derxy2[5][irn]*velint[2])      
				      *derxy[2][icn]*ccb              ;

	    estif[irow+2][icol]	  -=  ((derxy2[2][irn]+aux)*velint[2] 
	        	  	      + derxy2[4][irn]*velint[0]       
	        	  	      + derxy2[5][irn]*velint[1])      
				      *derxy[0][icn]*ccb              ;
	    estif[irow+2][icol+1] -=  ((derxy2[2][irn]+aux)*velint[2] 
	        	  	      + derxy2[4][irn]*velint[0]       
	        	  	      + derxy2[5][irn]*velint[1])      
				      *derxy[1][icn]*ccb              ;
	    estif[irow+2][icol+2] -=  ((derxy2[2][irn]+aux)*velint[2] 
	        	  	      + derxy2[4][irn]*velint[0]       
	        	  	      + derxy2[5][irn]*velint[1])      
				      *derxy[2][icn]*ccb              ;
	    irow += 3;
	  } /* end of loop over irn */
	  icol += 3;
        } /* end of loop over icn */
      }
    } /* endif (fdyn->nir!=0) */
  }    
} /* end of viscous stabilisation for higher order elments */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_lscalstabkvv */ 

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
void f3_lscalstabkvp(ELEMENT         *ele,
		   DOUBLE	  **estif,
		   DOUBLE	   *velint,
		   DOUBLE	  **vderxy,
		   DOUBLE	   *funct,
		   DOUBLE	  **derxy,
		   DOUBLE	  **derxy2,
		   DOUBLE	    fac,
		   DOUBLE	    visc,
		   INT  	    iel,
		   INT  	    ihoel)
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
INT    irow,icol,irn,ird,posc;
DOUBLE taumu;
DOUBLE taump;
DOUBLE con,beta,divv,cb;
DOUBLE aux;
DOUBLE signvi,signre;
STAB_PAR_GLS *gls;	/* pointer to GLS stabilisation parameters	*/

#ifdef DEBUG 
dstrc_enter("f3_lscalstabkvp");
#endif

/*---------------------------------------------------------- initialise */
fdyn = alldyn[genprob.numff].fdyn;

gls    = ele->e.f3->stabi.gls;

if (ele->e.f3->stab_type != stab_gls) 
   dserror("routine with no or wrong stabilisation called");

/*---------------------------------------- set stabilisation parameter */
taumu = fdyn->tau[0];
taump = fdyn->tau[1];

con = fac * taumu;

switch (fdyn->conte) /* determine parameter beta of convective term */
{
case 0: 
  beta = ZERO;
break;
case 1: 
  beta = ONE;
  divv= vderxy[0][0]+vderxy[1][1]+vderxy[2][2]; 
break;
case 2: 
  beta = ONE/TWO;
  divv= vderxy[0][0]+vderxy[1][1]+vderxy[2][2]; 
break;
default:
   dserror("unknown form of convective term");
} /* end switch (fdyn->conte) */

switch (gls->ivisc) /* determine type of stabilization */
{
case -1: /* USFEM (but no viscous stabilization) */
   signre = -ONE;
break;
case 0: /* GLS (but no viscous stabilization) */
   signre = ZERO;
break;
case 1: /* GLS- */
   signre = ONE;
   signvi = ONE;
break;
case 2: /* GLS+ */
   signre = ONE;
   signvi = -ONE;
break;
case 3: /* USFEM */
   signre = -ONE;
   signvi = -ONE;
break;
default:
   dserror("unknown type of stabilization");
} /* end switch (ele->e.f3->ivisc) */

/*----------------------------------------------------------------------*
   Calculate convection stabilisation:
 *----------------------------------------------------------------------*/
if (gls->iadvec!=0)
{
/*----------------------------------------------------------------------*
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
           + velint[2]*derxy[2][irn])*con;
      estif[irow][posc]   += derxy[0][icol]*aux;
      estif[irow+1][posc] += derxy[1][icol]*aux;
      estif[irow+2][posc] += derxy[2][icol]*aux;
      irow += 3;
    } /* end of loop over irn */
  } /* end of loop over icol */
    
  if (fdyn->conte!=0)  
  {
    cb = con*beta*signre;
/*----------------------------------------------------------------------*
    /
   |  +/- beta * tau_mu * v * div(u_old) * grad(p)  d_omega
  /
 *----------------------------------------------------------------------*/
    for (icol=0;icol<iel;icol++)
    {
      irow=0;
      posc=3*iel+icol;
      for (irn=0;irn<iel;irn++)
      {
    	aux = funct[irn]*divv*cb;
        estif[irow][posc]   += derxy[0][icol]*aux;
        estif[irow+1][posc] += derxy[1][icol]*aux;
        estif[irow+2][posc] += derxy[2][icol]*aux;
    	irow += 3;
      } /* end of loop over irn */
      icol += 3;
    } /* end of loop over icn */
  }  
} /* end of convection stabilisation */

/*----------------------------------------------------------------------*
   Calculate viscous stabilisation parts for higher order elements:
 *----------------------------------------------------------------------*/
if (gls->ivisc>0 && ihoel!=0)
{
  con = fac * taump * visc * signvi;
  
  if (fdyn->vite==0)
  {
/*----------------------------------------------------------------------*
    /
   |  -/+ tau_mp * nue * delta(v) * grad(p)  d_omega
  /
 *----------------------------------------------------------------------*/
    for (icol=0;icol<iel;icol++)
    {
      irow=0;
      posc=3*iel+icol;
      for (irn=0;irn<iel;irn++)
      {
        aux = (derxy2[0][irn] + derxy2[1][irn] + derxy2[2][irn])*con;
	estif[irow][posc]   -= aux*derxy[0][icol];
	estif[irow+1][posc] -= aux*derxy[1][icol];
	estif[irow+2][posc] -= aux*derxy[2][icol];
	irow += 3;			       
      } /* end loop over irn */
    } /* end loop over icol */
  }
  else
  {
/*----------------------------------------------------------------------*
    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * grad(p)  d_omega
  /
 *----------------------------------------------------------------------*/
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
	        	       + aux*derxy[0][icol])*con       ;
	estif[irow+1][posc] -=  (derxy2[3][irn]*derxy[0][icol] \
			       + derxy2[1][irn]*derxy[1][icol] \
			       + derxy2[5][irn]*derxy[2][icol] \
	        	       + aux*derxy[1][icol])*con       ;
	estif[irow+2][posc] -=  (derxy2[4][irn]*derxy[0][icol] \
			       + derxy2[5][irn]*derxy[1][icol] \
			       + derxy2[2][irn]*derxy[2][icol] \
	        	       + aux*derxy[2][icol])*con       ;
	irow += 3;			       
      } /* end of loop over irn */
    } /* end of loop over icol */
  }
} /*------------ end of viscous stabilisation for higher order elements */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_lscalstabkvp */

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
void f3_lscalstabmvv(ELEMENT         *ele,
		   DOUBLE	  **emass,
		   DOUBLE	   *velint,
		   DOUBLE	  **vderxy,
    		   DOUBLE	   *funct,
		   DOUBLE	  **derxy,
		   DOUBLE	  **derxy2,
		   DOUBLE	    fac,
		   DOUBLE	    visc,
		   INT  	    iel,
		   INT  	    ihoel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
 |   irn  - row node: number of node considered for matrix-row          |
 |   ird  - row dim.: number of spatial dimension at row node           |   
/*----------------------------------------------------------------------*/
INT    irow,icol,irn,icn,ird;
DOUBLE taumu;
DOUBLE taump;
DOUBLE con,beta,divv,cb;
DOUBLE aux,auxc,auxr;
DOUBLE signvi,signre;
STAB_PAR_GLS *gls;	/* pointer to GLS stabilisation parameters	*/

#ifdef DEBUG 
dstrc_enter("f3_lscalstabmvv");
#endif

/*---------------------------------------------------------- initialise */
fdyn = alldyn[genprob.numff].fdyn;

gls    = ele->e.f3->stabi.gls;

if (ele->e.f3->stab_type != stab_gls) 
   dserror("routine with no or wrong stabilisation called");

/*---------------------------------------- set stabilisation parameter */
taumu = fdyn->tau[0];
taump = fdyn->tau[1];

con = fac * taumu;

switch (fdyn->conte) /* determine parameter beta of convective term */
{
case 0: 
  beta = ZERO;
break;
case 1: 
  beta = ONE;
  divv= vderxy[0][0]+vderxy[1][1]+vderxy[2][2]; 
break;
case 2: 
  beta = ONE/TWO;
  divv= vderxy[0][0]+vderxy[1][1]+vderxy[2][2]; 
break;
default:
   dserror("unknown form of convective term");
} /* end switch (fdyn->conte) */

switch (gls->ivisc) /* determine type of stabilization */
{
case -1: /* USFEM (but no viscous stabilization) */
   signre = -ONE;
break;
case 0: /* GLS (but no viscous stabilization) */
   signre = ZERO;
break;
case 1: /* GLS- */
   signre = ONE;
   signvi = ONE;
break;
case 2: /* GLS+ */
   signre = ONE;
   signvi = -ONE;
break;
case 3: /* USFEM */
   signre = -ONE;
   signvi = -ONE;
break;
default:
   dserror("unknown type of stabilization");
} /* end switch (ele->e.f3->ivisc) */

/*----------------------------------------------------------------------*
   Calculate convection stabilisation part:
 *----------------------------------------------------------------------*/
if (gls->iadvec!=0)
{
/*----------------------------------------------------------------------*
    /
   |  -/+ tau_mu * u_old * grad(v) * u d_omega
  /
 *----------------------------------------------------------------------*/
  icol=0;
  for (icn=0;icn<iel;icn++)
  {
    auxc = funct[icn]*con;
    irow=0;
    for (irn=0;irn<iel;irn++)
    {
      aux = (velint[0]*derxy[0][irn] + velint[1]*derxy[1][irn] \
           + velint[2]*derxy[2][irn] )*auxc;
      emass[irow][icol]     += aux;
      emass[irow+1][icol+1] += aux;
      emass[irow+2][icol+2] += aux;
      irow += 3;
    } /* end of loop over irn */
    icol += 3;      
  } /* end of loop over icn */
    
  if (fdyn->conte!=0)  
  {
    cb = con*beta*signre;
/*----------------------------------------------------------------------*
    /
   |  +/- beta * tau_mu * v * div(u_old) * u  d_omega
  /
 *----------------------------------------------------------------------*/
    icol=0;
    for (icn=0;icn<iel;icn++)
    {
      irow=0;
      auxc = funct[icn]*divv*cb;
      for (irn=0;irn<iel;irn++)
      {
    	aux = funct[irn]*auxc;
        emass[irow][icol]     += aux;
        emass[irow+1][icol+1] += aux;
        emass[irow+2][icol+2] += aux;
    	irow += 3;
      } /* end of loop over irn */
      icol += 3;
    } /* end of loop over icn */
  }  
} /* end of convection stabilisation */

/*----------------------------------------------------------------------*
   Calculate viscous stabilisation parts for higher order elements:
 *----------------------------------------------------------------------*/
if (gls->ivisc>0 && ihoel!=0)
{
  con = fac * taump * visc * signvi;
  
  if (fdyn->vite==0)
  {
/*----------------------------------------------------------------------*
    /
   |  -/+ tau_mp * nue * delta(v) * u  d_omega
  /
 *----------------------------------------------------------------------*/
    icol=0;
    for (icn=0;icn<iel;icn++)
    {      
      auxc = funct[icn]*con;
      irow=0;
      for(irn=0;irn<iel;irn++)
      {
        aux = (derxy2[0][irn] + derxy2[1][irn] + derxy2[2][irn])*auxc;
	emass[irow][icol]     -= aux;
	emass[irow+1][icol+1] -= aux;
	emass[irow+2][icol+2] -= aux;
	irow += 3;
      } /* end loop over irn */
      icol += 3;
    } /* end loop over icn */
  }
  else
  {
/*----------------------------------------------------------------------*
    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * u  d_omega
  /
 *----------------------------------------------------------------------*/
    icol=0;
    for (icn=0;icn<iel;icn++)
    {      
      aux = funct[icn]*con;
      irow=0;
      for(irn=0;irn<iel;irn++)
      {
        auxr = derxy2[0][irn] + derxy2[1][irn] + derxy2[2][irn];
	emass[irow][icol]     -= (derxy2[0][irn] + auxr)*aux;
	emass[irow][icol+1]   -=  derxy2[3][irn]*aux;
	emass[irow][icol+2]   -=  derxy2[4][irn]*aux;

	emass[irow+1][icol]   -=  derxy2[3][irn]*aux;
	emass[irow+1][icol+1] -= (derxy2[1][irn] + auxr)*aux;
	emass[irow+1][icol+2] -=  derxy2[5][irn]*aux;

	emass[irow+2][icol]   -=  derxy2[4][irn]*aux;
	emass[irow+2][icol+1] -=  derxy2[5][irn]*aux;
	emass[irow+2][icol+2] -= (derxy2[2][irn] + auxr)*aux;
	irow += 3;
      } /* end of loop over irn */
      icol += 3;
    } /* end of loop over icn */
  }
} /* end of viscous stabilisation for higher order elements */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_lscalstabmvv */

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
void f3_lscalstabkpv(DOUBLE	  **estif,
		   DOUBLE	   *velint,
		   DOUBLE	  **vderxy,
		   DOUBLE	   *funct,
		   DOUBLE	  **derxy,
		   DOUBLE	  **derxy2,
		   DOUBLE	    fac,
		   DOUBLE	    visc,
		   INT  	    iel,
		   INT  	    ihoel)
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
INT    irow,icol,irn,ird,icn,posr;
DOUBLE con,beta,divv,cb;
DOUBLE aux;
DOUBLE taump;

#ifdef DEBUG 
dstrc_enter("f3_lscalstabkpv");
#endif

fdyn = alldyn[genprob.numff].fdyn;
/*---------------------------------------- set stabilisation parameter */
taump = fdyn->tau[1];

con = fac * taump;

switch (fdyn->conte) /* determine parameter beta of convective term */
{
case 0: 
  beta = ZERO;
break;
case 1: 
  beta = ONE;
  divv= vderxy[0][0]+vderxy[1][1]+vderxy[2][2]; 
break;
case 2: 
  beta = ONE/TWO;
  divv= vderxy[0][0]+vderxy[1][1]+vderxy[2][2]; 
break;
default:
   dserror("unknown form of convective term");
} /* end switch (fdyn->conte) */

/*----------------------------------------------------------------------*
   Calculate stabilisation part Nc(u):
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
    /
   |  - tau_mp * grad(q) * u_old * grad(u) d_omega
  /
 *----------------------------------------------------------------------*/
  icol=0;
  for (icn=0;icn<iel;icn++)
  {
    aux = (velint[0]*derxy[0][icn] + velint[1]*derxy[1][icn] \
         + velint[2]*derxy[2][icn])*con;
    for (irow=0;irow<iel;irow++)
    {
      posr = irow + 3*iel;
      estif[posr][icol]   -= derxy[0][irow]*aux;
      estif[posr][icol+1] -= derxy[1][irow]*aux;
      estif[posr][icol+2] -= derxy[2][irow]*aux;
    } /* end of loop over irow */
    icol += 3;
  } /* end of loop over icn */
  
  if (fdyn->conte!=0)  
  {
    cb = con*beta;
/*----------------------------------------------------------------------*
    /
   |  - beta * tau_mp * grad(q) * u * div(u_old) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=0;
    for (icn=0;icn<iel;icn++)
    {
      aux = funct[icn]*divv*cb;
      for (irow=0;irow<iel;irow++)
      {
        posr = irow + 3*iel;
        estif[posr][icol]   -= derxy[0][irow]*aux;
        estif[posr][icol+1] -= derxy[1][irow]*aux;
        estif[posr][icol+2] -= derxy[2][irow]*aux;
      } /* end loop over irow */
      icol += 3;
    } /* end loop over icn */
  }

/*----------------------------------------------------------------------*
   Calculate stabilisation part Nr(u):
 *----------------------------------------------------------------------*/
if (fdyn->nir!=0) /* evaluate for Newton iteration */
{
/*----------------------------------------------------------------------*
    /
   |  - tau_mp * grad(q) * u * grad(u_old) d_omega
  /
 *----------------------------------------------------------------------*/
  icol=0;
  for (icn=0;icn<iel;icn++)
  {
    aux = funct[icn]*con;
    for (irow=0;irow<iel;irow++)
    {
      posr = irow + 3*iel;
      estif[posr][icol]   -= aux*(derxy[0][irow]*vderxy[0][0]	\
        			+ derxy[1][irow]*vderxy[1][0]	\
        			+ derxy[2][irow]*vderxy[2][0])  ;
      estif[posr][icol+1] -= aux*(derxy[0][irow]*vderxy[0][1]	\
        			+ derxy[1][irow]*vderxy[1][1]	\
        			+ derxy[2][irow]*vderxy[2][1])  ;
      estif[posr][icol+2] -= aux*(derxy[0][irow]*vderxy[0][2]	\
        			+ derxy[1][irow]*vderxy[1][2]	\
        			+ derxy[2][irow]*vderxy[2][2])  ;
    } /* end of loop over irow */
    icol += 3;
  } /* end of loop over icn */
  
  if (fdyn->conte!=0)  
  {
    cb = con*beta;
/*----------------------------------------------------------------------*
    /
   |  - beta * tau_mp * grad(q) * u_old * div(u) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=0;
    for (icn=0;icn<iel;icn++)
    {
      for (irow=0;irow<iel;irow++)
      {
        posr = irow + 3*iel;
        aux = (velint[0]*derxy[0][irow] + velint[1]*derxy[1][irow] \
	     + velint[2]*derxy[2][irow])*cb;
	estif[posr][icol]   -= aux*derxy[0][icn];
	estif[posr][icol+1] -= aux*derxy[1][icn];
	estif[posr][icol+2] -= aux*derxy[2][icn];
      } /* end loop over irow */
      icol += 3;
    } /* end loop over icn */
  }
} /* endif (fdyn->nir!=0) */

/*----------------------------------------------------------------------*
   Calculate viscous stabilisation parts for higher order elements:
 *----------------------------------------------------------------------*/
if (ihoel!=0)
{
  con = con * visc;
  
  if (fdyn->vite==0)
  {
/*----------------------------------------------------------------------*
    /
   |  tau_mp * nue *grad(q) * delta(u) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=0;
    for (icn=0;icn<iel;icn++)
    {
      aux = (derxy2[0][icn] + derxy2[1][icn] + derxy2[2][icn])*con;
      for (irow=0;irow<iel;irow++)
      {
        posr = irow + 3*iel;
        estif[posr][icol]   += aux*derxy[0][irow];
        estif[posr][icol+1] += aux*derxy[1][irow];
        estif[posr][icol+1] += aux*derxy[2][irow];
      } /* end loop over irow */
      icol += 3;
    } /* end loop over icn */
  }
  else
  {
/*----------------------------------------------------------------------*
    /
   |  tau_mp * 2 * nue *grad(q) * div(eps(u)) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=0;
    for (icn=0;icn<iel;icn++)
    {
      aux = derxy2[0][icn] + derxy2[1][icn] + derxy2[2][icn];
      for (irow=0;irow<iel;irow++)
      {
        posr = irow + 3*iel;
        estif[posr][icol]   += ((derxy2[0][icn]+aux)*derxy[0][irow]   \
			       + derxy2[3][icn]     *derxy[1][irow]   \
			       + derxy2[4][icn]     *derxy[2][irow])*con;

        estif[posr][icol+1] += ( derxy2[3][icn]     *derxy[0][irow]   \
			       +(derxy2[1][icn]+aux)*derxy[1][irow]   \
			       + derxy2[5][icn]     *derxy[2][irow])*con;

        estif[posr][icol+2] += ( derxy2[4][icn]     *derxy[0][irow]   \
			       + derxy2[5][icn]     *derxy[1][irow]   \
			       +(derxy2[2][icn]+aux)*derxy[2][irow])*con;
      } /* end of loop over irow */
      icol += 3;
    } /* end of loop over icn */
  }
} /* end of stabilisation for higher order elements */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_lscalstabkpv */

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
void f3_lscalstabkpp(DOUBLE	  **estif,
		   DOUBLE	  **derxy,
		   DOUBLE	    fac,
		   INT  	    iel)
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
INT    irow,icol,posc,posr;
DOUBLE con;
DOUBLE taump;

#ifdef DEBUG 
dstrc_enter("f3_lscalstabkpp");
#endif

fdyn = alldyn[genprob.numff].fdyn;
/*---------------------------------------- set stabilisation parameter */
taump = fdyn->tau[1];

con = fac * taump;

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
                           +derxy[2][irow]*derxy[2][icol])*con ;
   } /* end of loop over irow */
} /* end of loop over icol */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_lscalstabkpp */

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
void f3_lscalstabmpv(DOUBLE	  **emass,
		   DOUBLE	   *funct,
		   DOUBLE	  **derxy,
		   DOUBLE	    fac,
		   INT  	    iel)
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
INT    irow,icol,icn,posr;
DOUBLE con;
DOUBLE taump;
DOUBLE auxc;

#ifdef DEBUG 
dstrc_enter("f3_lscalstabmpv");
#endif

fdyn = alldyn[genprob.numff].fdyn;
/*---------------------------------------- set stabilisation parameter */
taump = fdyn->tau[1];

con = fac * taump;

/*----------------------------------------------------------------------*
   Calculate stabilisation part for matrix Mpv:
    /
   |  - tau_mp * grad(q) * u d_omega
  /
 *----------------------------------------------------------------------*/
icol=0;
for (icn=0;icn<iel;icn++)
{
   auxc = funct[icn]*con;
   for (irow=0;irow<iel;irow++)
   {
      posr = irow + 3*iel;
      emass[posr][icol]   -= derxy[0][irow]*auxc;
      emass[posr][icol+1] -= derxy[1][irow]*auxc;
      emass[posr][icol+2] -= derxy[2][irow]*auxc;
   } /* end of loop over irow */
   icol += 3;
} /* end of loop over icn */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_lscalstabmpv */


#endif

/*!----------------------------------------------------------------------
\file
\brief evaluate "Time" force vectors for submesh element for fluid3

<pre>
Maintainer:  name
             e-mail
             homepage
             telephone number


</pre>

------------------------------------------------------------------------*/
#ifdef FLUID3_ML 
#include "../headers/standardtypes.h"
#include "fluid3ml_prototypes.h"
#include "../fluid3/fluid3.h"
/*!---------------------------------------------------------------------                                         
\brief evaluate Galerkin part of submesh "Time" force vector for fluid3

<pre>                                                       gravem 07/03

In this routine, the Galerkin part of the submesh "Time" force vector
on the rhs is calculated.

</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smetfor   DOUBLE	   (i/o)  submesh element time force vec.
\param  *velintn   DOUBLE          (i)    velocity at int. point at n
\param  *velintnt  DOUBLE          (i)    'temporal' vel. at int. p. at n
\param  *velintnc  DOUBLE          (i)    'convective' vel. at int. p. at n
\param **vderxyn   DOUBLE          (i)    global vel. derivatives at n
\param **vderxync  DOUBLE          (i)    global 'convective' vel. der. at n
\param **vderxynv  DOUBLE          (i)    global 'viscous' vel. der. at n
\param **vderxy2n  DOUBLE          (i)    2nd global vel. derivatives at n
\param  *pderxyn   DOUBLE          (i)    global pres. derivatives at n
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun.
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc 	   DOUBLE	   (i)    physical viscosity	      
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\param   iel	   INT  	   (i)	  number of nodes of l-s element
\param   ihoelsm   INT  	   (i)	  flag for higher-order sm ele.
\param   ihoel	   INT  	   (i)	  flag for higher-order l-s ele.
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calsmft(FLUID_DYN_CALC  *dynvar,
	        FLUID_DYN_ML	*mlvar, 
		DOUBLE         **smetfor,  
		DOUBLE  	*velintn, 
		DOUBLE  	*velintnt, 
		DOUBLE  	*velintnc, 
		DOUBLE         **vderxyn, 
		DOUBLE         **vderxync, 
		DOUBLE         **vderxynv, 
		DOUBLE         **vderxy2n, 
		DOUBLE          *pderxyn, 
		DOUBLE  	*smfunct,  
		DOUBLE         **smderxy,  
		DOUBLE  	 fac,	 
		DOUBLE  	 visc,   
		INT		 smiel,    
		INT		 iel,    
                INT		 ihoel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
/*----------------------------------------------------------------------*/
INT    irow,icol,icn,isd;
DOUBLE con,ccon,beta,divvc,cb,facsr,facpr;
DOUBLE aux;
DOUBLE sign;

#ifdef DEBUG 
dstrc_enter("f3_calsmft");
#endif

facsr = fac*dynvar->thsr;
facpr = fac*dynvar->thpr;
con   = facsr*visc;

switch (dynvar->conte) /* determine parameter beta of convective term */
{
case 0: 
  beta = ZERO;
break;
case 1: 
  beta = ONE;
  divvc= vderxync[0][0]+vderxync[1][1]+vderxync[2][2]; 
break;
case 2: 
  beta = ONE/TWO;
  divvc= vderxync[0][0]+vderxync[1][1]+vderxync[2][2]; 
break;
default:
   dserror("unknown form of convective term");
} /* end switch (dynvar->conte) */

/*----------------------------------------------------------------------*
   Calculate temporal forces (u_transient = u_n or ls_u_n or sm_u_n):
 *----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*
    /
   |  w * u_transient d_omega
  /
 *----------------------------------------------------------------------*/
icol=mlvar->nelbub-3;
for (icn=0;icn<3;icn++)
{
  aux = velintnt[icn]*fac;
  for (irow=0;irow<smiel;irow++)
  {
    smetfor[irow][icol] += smfunct[irow]*aux;
  } /* end loop over irow */
  icol++;
} /* end loop over icol */

if (dynvar->thsr!=ZERO)
{
/*----------------------------------------------------------------------*
   Calculate convective forces (u_convective = u_n or ls_u_n):
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
    /
   |  - ((1-theta)*dt) * w * u_convective * grad(u_n) d_omega
  /
 *----------------------------------------------------------------------*/
  icol=mlvar->nelbub-3;
  for (icn=0;icn<3;icn++)
  {
    aux = (velintnc[0]*vderxyn[icn][0]+velintnc[1]*vderxyn[icn][1]\
          +velintnc[2]*vderxyn[icn][2])*facsr;
    for (irow=0;irow<smiel;irow++)
    {
      smetfor[irow][icol] -= smfunct[irow]*aux;
    } /* end loop over irow */
    icol++;
  } /* end loop over icol */

  if (dynvar->conte!=0)  
  {
    cb = beta*divvc*facsr;
/*----------------------------------------------------------------------*
    /
   |  - ((1-theta)*dt) * beta * v * u_n * div(u_convective)  d_omega
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-3;
    for (icn=0;icn<3;icn++)
    {
      aux = velintn[icn]*cb;
      for (irow=0;irow<smiel;irow++)
      {
        smetfor[irow][icol] -= smfunct[irow]*aux;
      } /* end loop over irow */
      icol++;
    } /* end loop over icol */
  }  

/*----------------------------------------------------------------------*
   Calculate viscous forces (separate large- and small-scale part of u_n):
 *----------------------------------------------------------------------*/
  if (mlvar->quastabub==0)
  {
    ccon = facsr*(visc + mlvar->smsgvisc);
/*----------------------------------------------------------------------*
    /
   |  - ((1-theta)*dt) * grad(w) * (nue+nueT) * grad(sm_u_n) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-3;
    for (icn=0;icn<3;icn++)
    {
      for (irow=0;irow<smiel;irow++)
      {
        smetfor[irow][icol] -= (smderxy[0][irow]*vderxynv[0][icn]\
                               +smderxy[1][irow]*vderxynv[1][icn]\
                               +smderxy[2][irow]*vderxynv[2][icn]);
      } /* end loop over irow */
      icol++;
    } /* end loop over icol */
  }  

  if (ihoel!=0)
  {
/*----------------------------------------------------------------------*
    /
   |  ((1-theta)*dt) * w * nue * delta(ls_u_n)  d_omega
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-3;
    for (icn=0;icn<3;icn++)
    {
      aux = (vderxy2n[icn][0]+vderxy2n[icn][1]+vderxy2n[icn][2])*con;
      for (irow=0;irow<smiel;irow++)
      {
        smetfor[irow][icol] += smfunct[irow]*aux;
      } /* end loop over irow */
      icol++;
    } /* end loop over icol */
  } /* end of viscous stabilisation for higher order elements */
}  

/*----------------------------------------------------------------------*
   Calculate pressure forces:
 *----------------------------------------------------------------------*/
if (dynvar->thpr!=ZERO && dynvar->iprerhs>0)
{
/*----------------------------------------------------------------------*
    /
   |  - ((1-theta)*dt) * w * grad(p_n) d_omega
  /
 *----------------------------------------------------------------------*/
  icol=mlvar->nelbub-3;
  for (icn=0;icn<3;icn++)
  {
    aux = pderxyn[icn]*facpr;
    for (irow=0;irow<smiel;irow++)
    {
      smetfor[irow][icol] -= smfunct[irow]*aux;
    } /* end loop over irow */
    icol++;
  } /* end loop over icol */
}
    
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calsmft */

/*!---------------------------------------------------------------------                                         
\brief evaluate stabilization part of sm "Time" force vector for fluid3

<pre>                                                       gravem 07/03

In this routine, the stabilization part of the submesh "Time" force 
vector on the rhs is calculated.

</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smetfor   DOUBLE	   (i/o)  submesh element time force vec.
\param  *velintn   DOUBLE          (i)    velocity at int. point at n
\param  *velintnt  DOUBLE          (i)    'temporal' vel. at int. p. at n
\param  *velintnc  DOUBLE          (i)    'convective' vel. at int. p. at n
\param **vderxyn   DOUBLE          (i)    global vel. derivatives at n
\param **vderxync  DOUBLE          (i)    global 'convective' vel. der. at n
\param **vderxynv  DOUBLE          (i)    global 'viscous' vel. der. at n
\param **vderxy2n  DOUBLE          (i)    2nd global vel. derivatives at n
\param  *pderxyn   DOUBLE          (i)    global pres. derivatives at n
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun.
\param **smderxy2  DOUBLE	   (i)    sm 2nd global deriv. of shape fun.
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc 	   DOUBLE	   (i)    physical viscosity	      
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\param   iel	   INT  	   (i)	  number of nodes of l-s element
\param   ihoelsm   INT  	   (i)	  flag for higher-order sm ele.
\param   ihoel	   INT  	   (i)	  flag for higher-order l-s ele.
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabsmft(FLUID_DYN_CALC  *dynvar,
	            FLUID_DYN_ML    *mlvar, 
		    DOUBLE         **smetfor,  
		    DOUBLE  	    *velintn, 
		    DOUBLE  	    *velintnt, 
		    DOUBLE  	    *velintnc, 
		    DOUBLE         **vderxyn, 
		    DOUBLE         **vderxync, 
		    DOUBLE         **vderxynv, 
		    DOUBLE         **vderxy2n, 
		    DOUBLE          *pderxyn, 
	 	    DOUBLE  	    *smfunct,  
		    DOUBLE         **smderxy,  
		    DOUBLE         **smderxy2, 
		    DOUBLE  	     fac,	 
		    DOUBLE  	     visc,   
		    INT		     smiel,    
		    INT		     iel,    
                    INT		     ihoelsm,
                    INT		     ihoel)
{
/*----------------------------------------------------------------------*
 | NOTATION:                                                            |
 |   irow - row number in element matrix                                |
 |   icol - column number in element matrix                             |
/*----------------------------------------------------------------------*/
INT    irow,icol,icn;
DOUBLE con,ccon,beta,divvc,cb,ccb,cbb;
DOUBLE aux,auxc;
DOUBLE sign,tau;

#ifdef DEBUG 
dstrc_enter("f3_calstabsmft");
#endif

/*---------------------------------------- set stabilization parameter */
tau = mlvar->smtau;

switch (dynvar->conte) /* determine parameter beta of convective term */
{
case 0: 
  beta = ZERO;
break;
case 1: 
  beta = ONE;
  divvc= vderxync[0][0]+vderxync[1][1]+vderxync[2][2]; 
break;
case 2: 
  beta = ONE/TWO;
  divvc= vderxync[0][0]+vderxync[1][1]+vderxync[2][2]; 
break;
default:
   dserror("unknown form of convective term");
} /* end switch (dynvar->conte) */

if (mlvar->smstado<0) sign = -ONE;  /* USFEM */
else                  sign = ONE;   /* GLS- */

/*----------------------------------------------------------------------*
   Calculate temporal stabilization:
 *----------------------------------------------------------------------*/
if (ABS(mlvar->smstado)<3)
{ 
  if (mlvar->smstado==-1) con = fac*tau;
  else  		  con = fac*tau*sign; 
/*---------------------------------------------------------------------*
    /
   |  +/-  tau * (1/(theta*dt)) * w * u_transient d_omega
  /
 *----------------------------------------------------------------------*/
  icol=mlvar->nelbub-3;
  for (icn=0;icn<3;icn++)
  {
    aux = velintnt[icn]*con*(ONE/dynvar->thsl);
    for (irow=0;irow<smiel;irow++)
    {
      smetfor[irow][icol] += smfunct[irow]*aux;
    } /* end loop over irow */
    icol++;
  } /* end loop over icol */
      
  ccon = con*(dynvar->thsr/dynvar->thsl);
/*---------------------------------------------------------------------*
    /
   |  -/+  tau * ((1-theta)*dt/(theta*dt)) * w * u_convective * grad(u_n)
  /
 *----------------------------------------------------------------------*/
  icol=mlvar->nelbub-3;
  for (icn=0;icn<3;icn++)
  {
    aux = (velintnc[0]*vderxyn[icn][0]+velintnc[1]*vderxyn[icn][1]\
          +velintnc[2]*vderxyn[icn][2])*ccon;
    for (irow=0;irow<smiel;irow++)
    {
      smetfor[irow][icol] -= smfunct[irow]*aux;
    } /* end loop over irow */
    icol++;
  } /* end loop over icol */

  if (dynvar->conte!=0)  
  {
    cb = beta*divvc*ccon;
/*----------------------------------------------------------------------*
    /
   |  -/+ tau * ((1-theta)*dt/(theta*dt)) * beta * w * u_n * div(u_conv.)
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-3;
    for (icn=0;icn<3;icn++)
    {
      aux = velintn[icn]*cb;
      for (irow=0;irow<smiel;irow++)
      {
        smetfor[irow][icol] -= smfunct[irow]*aux;
      } /* end loop over irow */
      icol++;
    } /* end loop over icol */
  }  

  if (ihoel!=0)  
  {
    ccon = ccon*visc;
/*---------------------------------------------------------------------*
    /
   |  +/- tau * ((1-theta)*dt/(theta*dt)) * w * nue * delta(u_n) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-3;
    for (icn=0;icn<3;icn++)
    {
      aux = (vderxy2n[icn][0]+vderxy2n[icn][1]+vderxy2n[icn][2])*ccon;
      for (irow=0;irow<smiel;irow++)
      {
        smetfor[irow][icol] += smfunct[irow]*aux;
      } /* end loop over irow */
      icol++;
    } /* end loop over icol */
  }
  
  if (dynvar->thpr!=ZERO && dynvar->iprerhs>0)
  {
    ccon = con*(dynvar->thpr/dynvar->thpl);
/*----------------------------------------------------------------------*
    /
   |  -/+ tau * ((1-theta)*dt/(theta*dt)) * w * grad(p_n) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-3;
    for (icn=0;icn<3;icn++)
    {
      aux = pderxyn[icn]*ccon;
      for (irow=0;irow<smiel;irow++)
      {
        smetfor[irow][icol] -= smfunct[irow]*aux;
      } /* end loop over irow */
      icol++;
    }  
  }
}  

/*----------------------------------------------------------------------*
   Calculate convective stabilization:
 *----------------------------------------------------------------------*/
con = fac*tau;
/*----------------------------------------------------------------------*
    /
   |  tau * u_convective * grad(w) * u_transient d_omega
  /
 *----------------------------------------------------------------------*/
icol=mlvar->nelbub-3;
for (icn=0;icn<3;icn++)
{
  auxc = velintnt[icn]*con;
  for (irow=0;irow<smiel;irow++)
  {
    aux = (velintnc[0]*smderxy[0][irow]+velintnc[1]*smderxy[1][irow]\
          +velintnc[2]*smderxy[2][irow])*auxc;
    smetfor[irow][icol] += aux;
  } /* end loop over irow */
  icol++;
} /* end loop over icol */
      
if (dynvar->conte!=0)  
{
  cb = con*beta;
/*----------------------------------------------------------------------*
    /
   |  +/- tau * beta * w * div(u_convective) * u_transient  d_omega
  /
 *----------------------------------------------------------------------*/
  icol=mlvar->nelbub-3;
  for (icn=0;icn<3;icn++)
  {
    aux = velintnt[icn]*divvc*cb;
    for (irow=0;irow<smiel;irow++)
    {
      smetfor[irow][icol] += smfunct[irow]*aux;
    } /* end loop over irow */
    icol++;
  } /* end loop over icol */
}  

if (dynvar->thsr!=ZERO)
{
  ccon = con*dynvar->thsr;
/*----------------------------------------------------------------------*
    /
   |  - tau * ((1-theta)*dt) * u_con. * grad(w) * u_con. * grad(u_n) d_omega
  /
 *----------------------------------------------------------------------*/
  icol=mlvar->nelbub-3;
  for (icn=0;icn<3;icn++)
  {
    auxc = (velintnc[0]*vderxyn[icn][0]+velintnc[1]*vderxyn[icn][1]\
           +velintnc[2]*vderxyn[icn][2])*ccon;
     for (irow=0;irow<smiel;irow++)
    {
      aux = (velintnc[0]*smderxy[0][irow]+velintnc[1]*smderxy[1][irow]\
            +velintnc[1]*smderxy[2][irow])*auxc;
      smetfor[irow][icol] -= aux;
    } /* end loop over irow */
    icol++;
  } /* end loop over icol */
      
  if (dynvar->conte!=0)  
  {
    cb  = ccon*beta;
    cbb = cb*beta*sign;
/*----------------------------------------------------------------------*
    /
   |  - beta * tau * ((1-th)*dt) * u_con. * grad(w) * u_n * div(u_con.)
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-3;
    for (icn=0;icn<3;icn++)
    {
      auxc = velintn[icn]*divvc*cb;
       for (irow=0;irow<smiel;irow++)
      {
        aux = (velintnc[0]*smderxy[0][irow]+velintnc[1]*smderxy[1][irow]\
	      +velintnc[2]*smderxy[2][irow])*auxc;
        smetfor[irow][icol] -= aux;
      } /* end loop over irow */
      icol++;
    } /* end loop over icol */
  
    cb = cb*sign;
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * tau * ((1-th)*dt) * w * div(u_con.) * u_con. * grad(u_n)
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-3;
    for (icn=0;icn<3;icn++)
    {
      aux = (velintnc[0]*vderxyn[icn][0]+velintnc[1]*vderxyn[icn][1]\
            +velintnc[2]*vderxyn[icn][2])*divvc*cb;
      for (irow=0;irow<smiel;irow++)
      {
        smetfor[irow][icol] -= aux*smfunct[irow];
      } /* end loop over irow */
      icol++;
    } /* end loop over icol */
  
/*----------------------------------------------------------------------*
    /
   |  -/+ beta * beta * tau * ((1-th)*dt) * v * div(u_con.) * u_n * div(u_con.)
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-3;
    for (icn=0;icn<3;icn++)
    {
      aux = velintn[icn]*divvc*divvc*cbb;
      for (irow=0;irow<smiel;irow++)
      {
        smetfor[irow][icol] -= aux*smfunct[irow];
      } /* end loop over irow */
      icol++;
    } /* end loop over icol */
  }

  if (ihoel!=0)
  {
    ccon = con*dynvar->thsr*visc;
/*----------------------------------------------------------------------*
    /
   |  tau * ((1-theta)*dt) * u_con. * grad(w) * nue * delta(u_n) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-3;
    for (icn=0;icn<3;icn++)
    {
      auxc = (vderxy2n[icn][0]+vderxy2n[icn][1]+vderxy2n[icn][2])*ccon;
      for (irow=0;irow<smiel;irow++)
      {
        aux = (velintnc[0]*smderxy[0][irow]+velintnc[1]*smderxy[1][irow]\
	      +velintnc[2]*smderxy[2][irow])*auxc;
        smetfor[irow][icol] += aux;
      } /* end loop over irow */
      icol++;
    } /* end loop over icol */
      
    if (dynvar->conte!=0)  
    {
      cb = ccon*beta*sign;
/*----------------------------------------------------------------------*
    /
   |  +/- beta * tau * ((1-th)*dt) * w * div(u_con.) * nue * delta(u_n)
  /
 *----------------------------------------------------------------------*/
      icol=mlvar->nelbub-3;
      for (icn=0;icn<3;icn++)
      {
        aux = (vderxy2n[icn][0]+vderxy2n[icn][1]+vderxy2n[icn][2])*divvc*cb;
        for (irow=0;irow<smiel;irow++)
        {
          smetfor[irow][icol] += aux*smfunct[irow];
        } /* end loop over irow */
        icol++;
      } /* end loop over icol */
    }
  }
}      
  
if (dynvar->thpr!=ZERO && dynvar->iprerhs>0)
{
  ccon = con*dynvar->thpr;
/*----------------------------------------------------------------------*
    /
   |  - tau * ((1-theta)*dt) * u_convective * grad(w) * grad(p_n) d_omega
  /
 *----------------------------------------------------------------------*/
  icol=mlvar->nelbub-3;
  for (icn=0;icn<3;icn++)
  {
    auxc = pderxyn[icn]*ccon;
    for (irow=0;irow<smiel;irow++)
    {
      aux = (velintnc[0]*smderxy[0][irow]+velintnc[1]*smderxy[1][irow]\
            +velintnc[2]*smderxy[2][irow])*auxc;
      smetfor[irow][icol] -= aux;
    } /* end loop over irow */
    icol++;
  }  

  if (dynvar->conte!=0)  
  {
    cb = ccon*beta*sign;
/*----------------------------------------------------------------------*
    /
   |  -/+ tau * ((1-theta)*dt) * beta * w * div(u_old) * grad(p_n) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-3;
    for (icn=0;icn<3;icn++)
    {
      aux = pderxyn[icn]*divvc*cb;
      for (irow=0;irow<smiel;irow++)
      {
        smetfor[irow][icol] -= smfunct[irow]*aux;
      } /* end loop over irow */
      icol++;
    }
  }    
}
    
if (ihoelsm!=0)
{
/*----------------------------------------------------------------------*
   Calculate viscous stabilization for higher-order elements:
 *----------------------------------------------------------------------*/
  ccon = con*visc*sign;
/*----------------------------------------------------------------------*
    /
   |  -/+ tau * nue * delta(w) * u_transient d_omega
  /
 *----------------------------------------------------------------------*/
  icol=mlvar->nelbub-3;
  for (icn=0;icn<3;icn++)
  {
    auxc = velintnt[icn]*ccon;
    for (irow=0;irow<smiel;irow++)
    {
      aux = (smderxy2[0][irow]+smderxy2[1][irow]+smderxy2[2][irow])*auxc;
      smetfor[irow][icol] -= aux;
    } /* end loop over irow */
    icol++;
  } /* end loop over icol */

  if (dynvar->thsr!=ZERO)
  {
    ccon = con*dynvar->thsr*visc*sign;
/*---------------------------------------------------------------------*
    /
   |  +/- tau * ((1-theta)*dt) * nue * delta(w) * u_convective * grad(u_n)
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-3;
    for (icn=0;icn<3;icn++)
    {
      auxc = (velintnc[0]*vderxyn[icn][0]+velintnc[1]*vderxyn[icn][1]\
             +velintnc[2]*vderxyn[icn][2])*ccon;
      for (irow=0;irow<smiel;irow++)
      {
        aux = (smderxy2[0][irow] + smderxy2[1][irow])*auxc;
        smetfor[irow][icol] += aux;
      }  /* end loop over irow */
      icol++;
    } /* end loop over icol */

    if (dynvar->conte!=0)  
    {
      cb = beta*ccon;
/*----------------------------------------------------------------------*
    /
   |  +/- tau * ((1-theta)*dt) * beta * nue * delta(w) *  * u_n * div(u_conv.)
  /
 *----------------------------------------------------------------------*/
      icol=mlvar->nelbub-3;
      for (icn=0;icn<3;icn++)
      {
        auxc = velintn[icn]*divvc*cb;
        for (irow=0;irow<smiel;irow++)
        {
         aux = (smderxy2[0][irow]+smderxy2[1][irow]+smderxy2[2][irow])*auxc;
         smetfor[irow][icol] += aux;
        } /* end loop over irow */
        icol++;
      } /* end loop over icol */
    }  
    
    if (ihoel!=0)
    {
      ccon = ccon*visc;
/*---------------------------------------------------------------------*
    /
   |  -/+ tau * ((1-theta)*dt) * nue * delta(w) * nue * delta(u_n) d_omega
  /
 *----------------------------------------------------------------------*/
      icol=mlvar->nelbub-3;
      for (icn=0;icn<3;icn++)
      {
        auxc = (vderxy2n[icn][0]+vderxy2n[icn][1]+vderxy2n[icn][2])*ccon;
        for (irow=0;irow<smiel;irow++)
        {
          aux = (smderxy2[0][irow]+smderxy2[1][irow]+smderxy2[2][irow])*auxc;
          smetfor[irow][icol] -= aux;
        } /* end loop over irow */
        icol++;
      } /* end loop over icol */
    }  
  }  
 
  if (dynvar->thpr!=ZERO && dynvar->iprerhs>0)
  {
    ccon = con*dynvar->thpr*visc*sign;
/*----------------------------------------------------------------------*
    /
   |  +/- tau * ((1-theta)*dt) * nue * delta(w) * grad(p_n) d_omega
  /
 *----------------------------------------------------------------------*/
    icol=mlvar->nelbub-3;
    for (icn=0;icn<3;icn++)
    {
      auxc = pderxyn[icn]*ccon;
      for (irow=0;irow<smiel;irow++)
      {
        aux = (smderxy2[0][irow]+smderxy2[1][irow]+smderxy2[2][irow])*auxc;
        smetfor[irow][icol] += aux;
      } /* end loop over irow */
      icol++;
    }  
  }
}  

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calstabsmft */

#endif

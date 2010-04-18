/*!----------------------------------------------------------------------
\file
\brief submesh stability parameter for fluid2 element

<pre>
Maintainer: Volker Gravemeier
            vgravem@stanford.edu


</pre>

------------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef FLUID2_ML
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "fluid2ml_prototypes.h"
#include "../fluid2/fluid2.h"
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

static DOUBLE Q13  = ONE/THREE;
static DOUBLE Q112 = ONE/TWELVE;

/*!---------------------------------------------------------------------
\brief routine to calculate submesh stability parameter for fluid2

<pre>                                                       gravem 07/03

In this routine, the stability parameter for the actual submesh element
is calculated.

</pre>

\param   *ele         ELEMENT	      (i)    actual element
\param   *mlvar       FLUID_DYN_ML    (i/o)
\param   *velint      DOUBLE	      (i)    l-s velocity at ele. center
\param    visc        DOUBLE	      (i)    viscosity
\param    iel         INT	      (i)    number of nodes
\param	  typ         DIS_TYP	      (i)    element type
\return void

------------------------------------------------------------------------*/
void f2_smstabpar(ELEMENT         *ele,
		  FLUID_DYN_ML    *mlvar,
		  DOUBLE	  *velint,
		  DOUBLE	   visc,
		  INT		   iel,
		  DIS_TYP          typ)
{
DOUBLE hdiv=ONE;   /* element length quotient for higher-order elements */
DOUBLE velno;      /* velocity norm                                     */
DOUBLE c_mk;       /* parameter mk                                      */
DOUBLE dt;         /* time step size                                    */
DOUBLE rc;         /* reactive coefficient                              */
DOUBLE re,re1,re2; /* various element Reynolds numbers                  */
DOUBLE hk;         /* characteristic element length                     */
DOUBLE aux1;       /* auxiliary parameter                               */

#ifdef DEBUG
dstrc_enter("f2_smstabpar");
#endif

fdyn = alldyn[genprob.numff].fdyn;
/*----------------------- higher order element diameter modifications ? */
switch(mlvar->smstamk)
{
case -1:
   c_mk = Q13;
   if (typ==quad8 || typ == quad9)
   {
      if (iel<10)
         hdiv = TWO;
      else
         hdiv = THREE;
   }
   else if (typ==tri6)
   {
      if (iel==6)
         hdiv = TWO;
      else
         hdiv = THREE;
   }
break;
case 0:
   if (iel>=6)
      c_mk=Q112;
   else
      c_mk=Q13;
break;
default:
   dserror("mk > 0 not implemented yet!\n");
} /* end switch (mlvar->smstamk) */

/*---------------------------------- choose stability-parameter version */
switch(mlvar->smstapa)
{
case 35: /*-------------------------- version diss. Wall - instationary */
  velno = sqrt(velint[0]*velint[0] + velint[1]*velint[1]); /*norm of vel */
  dt = fdyn->dta;	/* check if dta or dt has to be chosen!!!!!!!! */
  hk = ele->e.f2->smcml/hdiv;
  if (velno>EPS15)
  {
    aux1 = DMIN(hk/TWO/velno,c_mk*hk*hk/FOUR/visc);
    mlvar->smtau = DMIN(dt,aux1);
  }
  else
    mlvar->smtau = DMIN(dt,c_mk*hk*hk/FOUR/visc);
break;

case 36: /*---------------------------- version diss. Wall - stationary */
  velno = sqrt(velint[0]*velint[0] + velint[1]*velint[1]); /*norm of vel*/
  aux1= velno*c_mk/FOUR/visc;
  hk = ele->e.f2->smcml/hdiv;
  re = aux1*hk;
  if (re<ONE) mlvar->smtau = c_mk*hk*hk/FOUR/visc;
  else        mlvar->smtau = hk/TWO/velno;
break;

case 37: /*------------------------------ version Franca/Valentin (2000) */
  velno = sqrt(velint[0]*velint[0] + velint[1]*velint[1]); /*norm of vel */
  rc = fdyn->thsl;	 /* reaction coefficient = theta times dt */
  hk = ele->e.f2->smcml/hdiv;
  re1 = TWO*visc*rc/c_mk/hk/hk; /* first element reynolds number */
  if (re1<ONE) re1 = ONE;
  re2 = c_mk*velno*hk/visc; /* second element reynolds number */
  if (re2<ONE) re2 = ONE;
  mlvar->smtau = hk*hk/(hk*hk*re1/rc+TWO*visc*re2/c_mk);
break;

default:
   dserror("submesh stability parameter version unknown!\n");
} /* end switch (mlvar->smstapa) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_smstabpar*/

/*!---------------------------------------------------------------------
\brief routine to calculate submesh subgrid viscosity for fluid2

<pre>                                                       gravem 07/03

In this routine, the subgrid viscosity for the actual submesh element
is calculated.

</pre>

\param   *ele         ELEMENT	      (i)    actual element
\param   *mlvar       FLUID_DYN_ML    (i/o)
\param   *velint      DOUBLE	      (i)    l-s velocity at ele. center
\param  **vderxy      DOUBLE	      (i)    l-s vel. der. at ele. center
\param    visc        DOUBLE	      (i)    viscosity
\param    iel         INT	      (i)    number of nodes
\param	  typ         DIS_TYP	      (i)    element type
\return void

------------------------------------------------------------------------*/
void f2_smsgvisc(ELEMENT         *ele,
                 FLUID_DYN_ML    *mlvar,
		 DOUBLE 	 *velint,
		 DOUBLE 	**vderxy,
		 DOUBLE 	  visc,
		 INT		  iel,
		 DIS_TYP          typ)
{
DOUBLE hdiv=ONE;   /* element length quotient for higher-order elements */
DOUBLE velno;      /* velocity norm                                     */
DOUBLE norovt;     /* norm of rate-of-velocity tensor                   */
DOUBLE pe;         /* element Peclet number                             */
DOUBLE hk;         /* characteristic element length                     */
DOUBLE auxfu;      /* auxiliary function                                */

#ifdef DEBUG
dstrc_enter("f2_smsgvisc");
#endif


/*----------------------- higher order element diameter modifications ? */
if (typ==quad8 || typ==quad9)
{
   if (iel<10)
      hdiv = TWO;
   else
      hdiv = THREE;
}
else if (typ==tri6)
{
   if (iel==6)
      hdiv = TWO;
   else
      hdiv = THREE;
}

/*-------------------------------------- characteristic element length */
hk = ele->e.f2->smcml/hdiv;

/*----------------------------------- choose type of subgrid viscosity */
switch(mlvar->smsgvi)
{
case 1: /*--------- element-Peclet-number version Brooks/Hughes (1982) */
/*--------------------------------------------------- norm of velocity */
   velno = sqrt(velint[0]*velint[0] + velint[1]*velint[1]);
/*---------------------------------------------- element Peclet number */
   pe = velno*hk/TWO/visc;
/*------------------------------------------------- auxiliary function */
   if (pe>EPS15) auxfu=ONE/((exp(pe)-exp(-pe))/(exp(pe)+exp(-pe)))-ONE/pe;
   else auxfu = ZERO;
/*-------------------------------------------------- subgrid viscosity */
   mlvar->smsgvisc = pe*visc*auxfu;
break;

case 2: /*--------------------------------- Smagorinsky version (1963) */
/*------------------------------------ norm of rate-of-velocity tensor */
   norovt = sqrt(TWO*(vderxy[0][0]*vderxy[0][0]+vderxy[1][1]*vderxy[1][1]\
                     +vderxy[0][1]*vderxy[1][0])+vderxy[0][1]*vderxy[0][1]\
                     +vderxy[1][0]*vderxy[1][0]);
/*-------------------------------------------------- subgrid viscosity */
   mlvar->smsgvisc = mlvar->smsmagcon*mlvar->smsmagcon*hk*hk*norovt;
break;

default:
   dserror("subgrid viscosity version unknown!\n");
} /* end switch (mlvar->smsgvi) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_smsgvisc*/

/*!---------------------------------------------------------------------
\brief routine to calculate stability parameter for large-scale element

<pre>                                                       gravem 05/05

   ele->e.f2->stabi.gls->iadvec: advection stab.
      0 = no
      1 = yes
   ele->e.f2->stabi.gls->ipres: pressure stab.
      0 = no
      1 = yes
   ele->e.f2->stabi.gls->ivisc: diffusion stab.
      0 = no
      1 = GLS-
      2 = GLS+
   ele->e.f2->stabi.gls->icont: continuity stab.
      0 = no
      1 = yes
   ele->e.f2->stabi.gls->istapa: version of stab. parameter
      35 = diss wall instationary
      36 = diss wall stationanary
   ele->e.f2->stabi.gls->norm_P: p-norm
      p = 1<=p<=oo
      0 = max.-norm (p=oo)
   ele->e.f2->stabi.gls->mk: higher order elements control flag
      0 = mk fixed (--> (bi)linear: 1/3, biquadr.: 1/12)
      1 = min(1/3,2*C)
     -1 = mk=1/3  (--> element order via approx. nodal-distance)
   ele->e.f2->stabi.gls->ihele[]:
      x/y/z = length-def. for velocity/pressure/continuity stab
      0 = don't compute
      1 = sqrt(area)
      2 = area equivalent diameter
      3 = diameter/sqrt(2)
      4 = sqrt(2)*area/diagonal (rectangle) 4*area/s (triangle)
      5 = streamlength (element length in flow direction
   ele->e.f2->stabi.gls->ninths: number of integration points for streamlength
      1 = at center of element
      2 = at every INT pt used for element.-stab.-matrices
   ele->e.f2->stabi.gls->istapc: flag for stabilisation parameter calculation
      1 = at center of element
      2 = at every integration point
   ele->e.f2->stabi.gls->clamb \
   ele->e.f2->c1               |_>> stabilisation constants (input)
   ele->e.f2->c2               |
   ele->e.f2->c3              /
   ele->e.f2->stabi.gls->istrle: has streamlength to be computed
   ele->e.f2->stabi.gls->iarea: calculation of area length
   ele->e.f2->stabi.gls->iduring: calculation during INT.-pt.loop
   ele->e.f2-->stabi.gls>itau[0]: flag for tau_mu calc. (-1: before, 1:during)
   ele->e.f2->stabi.gls->itau[1]: flag for tau_mp calc. (-1: before, 1:during)
   ele->e.f2->stabi.gls->itau[2]: flag for tau_c calc. (-1: before, 1:during)
   ele->e.f2->hk[i]: element sizes (vel / pre / cont)
   ele->e.f2->stabi.gls->idiaxy: has diagonals to be computed
   fdyn->tau[0]: stability parameter momentum / velocity (tau_mu)
   fdyn->tau[1]: stability parameter momentum / pressure (tau_mp)
   fdyn->tau[2]: stability parameter continuity (tau_c)

</pre>

\param   *ele,        ELEMENT	      (i)    actual element
\param   *velint,     DOUBLE	      (i)    vel at center
\param    visc,       DOUBLE	      (i)    viscosity
\param    iel,        INT	      (i)    number of nodes
\param	  typ,        DIS_TYP	      (i)    element type
\param	  iflag       INT	      (i)    flag for evaluation
\return void

------------------------------------------------------------------------*/
void f2_lsstabpar(ELEMENT       *ele,
		  DOUBLE	*velint,
		  DOUBLE         visc,
		  INT  	         iel,
		  DIS_TYP        typ,
		  INT  	         iflag)
{
INT    isp;
DOUBLE hdiv=ONE;
DOUBLE velno;
DOUBLE c_mk;
DOUBLE dt,rc;
DOUBLE re,re1,re2;
DOUBLE hk;
DOUBLE aux1;
STAB_PAR_GLS *gls;	/* pointer to GLS stabilisation parameters	*/

#ifdef DEBUG
dstrc_enter("f2_lsstabpar");
#endif

/*---------------------------------------------------------- initialise */
fdyn = alldyn[genprob.numff].fdyn;
gls    = ele->e.f2->stabi.gls;

/*----------------------- higher order element diameter modifications ? */
switch(gls->mk)
{
case -1:
   c_mk = Q13;
   if (typ == quad8 || typ == quad9)
   {
      if (iel<10)
         hdiv = TWO;
      else
         hdiv = THREE;
   }
   else if (typ==tri6)
   {
      if (iel==6)
         hdiv = TWO;
      else
         hdiv = THREE;
   }
break;
case 0:
   if (iel>=6)
      c_mk=Q112;
   else
      c_mk=Q13;
break;
default:
   dserror("mk > 0 not implemented yet!\n");
} /* end switch (ele->e.f2->mk) */
/*---------------------------------- choose stability-parameter version */
switch(gls->istapa)
{
case 35: /*-------------------------- version diss. Wall - instationary */
   velno = sqrt(velint[0]*velint[0] + velint[1]*velint[1]); /*norm of vel */
   dt = fdyn->dta;     /* check if dta or dt has to be chosen!!!!!!!! */
   for (isp=0;isp<3;isp++)
   {
      if (gls->itau[isp]!=iflag)
         continue;
      hk = ele->e.f2->hk[isp]/hdiv;
      switch(isp)
      {
      case 2:/* continiuty stabilisation */
         re = c_mk*hk*velno/TWO/visc;  /* element reynolds number */
	 fdyn->tau[isp] = (gls->clamb)*velno*hk/TWO*DMIN(ONE,re);
      break;
      default: /* velocity / pressure stabilisation */
         if (velno>EPS15)
	 {
	    aux1 = DMIN(hk/TWO/velno , c_mk*hk*hk/FOUR/visc);
            fdyn->tau[isp] = DMIN(dt , aux1);
         }
	 else
            fdyn->tau[isp] = DMIN(dt , c_mk*hk*hk/FOUR/visc);
       break;
      } /* end switch (isp) */
   } /* end of loop over isp */
break;

case 36: /*---------------------------- version diss. Wall - stationary */
   velno = sqrt(velint[0]*velint[0] + velint[1]*velint[1]); /*norm of vel*/
   aux1= velno*c_mk/FOUR/visc;
   for (isp=0;isp<3;isp++)
   {
      if (gls->itau[isp]!=iflag)
         continue;
      hk = ele->e.f2->hk[isp]/hdiv;
      re = aux1*hk;
      switch(isp)
      {
      case 2: /* continuity stabilisation */
         fdyn->tau[isp] = (gls->clamb)*velno*hk/TWO*DMIN(ONE,re);
/*        fdyn->tau[isp] = velno*hk/TWO*DMIN(ONE,re); */
      break;
      default: /* velocity / pressure stabilisation */
         if (re<ONE)
	    fdyn->tau[isp] = c_mk*hk*hk/FOUR/visc;
	 else
	    fdyn->tau[isp] = hk/TWO/velno;
      break;
      }
   }
break;

case 37: /*------------------------------- version Franca/Valentin (2000) */
   velno = sqrt(velint[0]*velint[0] + velint[1]*velint[1]); /*norm of vel */
   rc = fdyn->thsl;     /* reaction coefficient = theta times dt */
   for (isp=0;isp<3;isp++)
   {
     if (gls->itau[isp]!=iflag)
       continue;
     hk = ele->e.f2->hk[isp]/hdiv;
     switch(isp)
     {
     case 2:/* continuity stabilisation */
       re = c_mk*hk*velno/TWO/visc;  /* element reynolds number */
       fdyn->tau[isp] = (gls->clamb)*velno*hk/TWO*DMIN(ONE,re);
     break;
     default: /* velocity / pressure stabilisation */
       re1 = TWO*visc*rc/c_mk/hk/hk; /* first element reynolds number */
       if (re1<ONE)
	 re1 = ONE;

       re2 = c_mk*velno*hk/visc; /* second element reynolds number */
       if (re2<ONE)
	 re2 = ONE;

       fdyn->tau[isp] = hk*hk/(hk*hk*re1/rc+TWO*visc*re2/c_mk);
     break;
     } /* end switch (isp) */
   } /* end of loop over isp */
break;

default:
   dserror("stability parameter version ISTAP unknown!\n");
} /* end switch (ele->e.f2->istapa) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_lsstabpar*/

/*!---------------------------------------------------------------------
\brief routine to calculate subgrid viscosity for large-scale element

<pre>                                                      gravem 07/03

In this routine, two versions of a subgrid viscosity (based on the
element Peclet number or the Smagorinsky version) may be calculated.

</pre>

\param   *ele,        ELEMENT	      (i)    actual element
\param   *velint      DOUBLE	      (i)    vel. at center
\param  **vderxy      DOUBLE	      (i)    vel. der. at center
\param    visc        DOUBLE	      (i)    viscosity
\param    iel,        int	      (i)    number of nodes
\param	  typ         DIS_TYP	      (i)    element type
\return void

------------------------------------------------------------------------*/
void f2_lssgvisc(ELEMENT         *ele,
		 DOUBLE          *velint,
		 DOUBLE         **vderxy,
		 DOUBLE           visc,
		 INT              iel,
		 DIS_TYP          typ)
{
DOUBLE hdiv=ONE;
DOUBLE velno,norovt;
DOUBLE pek;
DOUBLE hk;
DOUBLE auxfu;

#ifdef DEBUG
dstrc_enter("f2_lssgvisc");
#endif


fdyn = alldyn[genprob.numff].fdyn;

/*----------------------- higher order element diameter modifications ? */
if (typ==quad8 || typ==quad9)
{
   if (iel<10)
      hdiv = TWO;
   else
      hdiv = THREE;
}
else if (typ==tri6)
{
   if (iel==6)
      hdiv = TWO;
   else
      hdiv = THREE;
}
hk = ele->e.f2->hk[2]/hdiv;

/*---------------------------------- choose type of subgrid viscosity */
switch(fdyn->sgvisc)
{
case 1: /*-------- element Peclet number version Brooks/Hughes (1982) */
/*-------------------------------------------------- norm of velocity */
   velno = sqrt(velint[0]*velint[0] + velint[1]*velint[1]);
/*------------------------------------------- element Peclet number */
   pek = velno*hk/TWO/visc;
/*------------------------------------------------- auxiliary function */
   if (pek>EPS15) auxfu=ONE/((exp(pek)-exp(-pek))/(exp(pek)+exp(-pek)))-ONE/pek;
   else auxfu = ZERO;
/*-------------------------------------------------- subgrid viscosity */
   fdyn->sugrvisc = pek*visc*auxfu;
break;

case 2: /*--------------------------------- Smagorinsky version (1963) */
/*----------------------------------- norm of rate-of-velocity tensor */
   norovt = sqrt(TWO*(vderxy[0][0]*vderxy[0][0]+vderxy[1][1]*vderxy[1][1]\
                     +vderxy[0][1]*vderxy[1][0])+vderxy[0][1]*vderxy[0][1]\
                     +vderxy[1][0]*vderxy[1][0]);
/*-------------------------------------------------- subgrid viscosity */
   fdyn->sugrvisc = fdyn->smagcon*fdyn->smagcon*hk*hk*norovt;
break;

default:
   dserror("subgrid viscosity version unknown!\n");
} /* end switch (fdyn->sgvisc) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_lssgvisc*/

#endif
#endif

/*!----------------------------------------------------------------------
\file
\brief stability parameter for fluid2 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID2
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
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

static DOUBLE Q13  = ONE/THREE;
static DOUBLE Q112 = ONE/TWELVE;
static FLUID_DYNAMIC *fdyn;
/*!---------------------------------------------------------------------
\brief routine to calculate stability parameter

<pre>                                                         genk 04/02

   ele->e.f2->iadvec: advection stab.
      0 = no
      1 = yes
   ele->e.f2->ipres: pressure stab.
      0 = no
      1 = yes
   ele->e.f2->ivisc: diffusion stab.
      0 = no
      1 = GLS-
      2 = GLS+
   ele->e.f2->icont: continuity stab.
      0 = no
      1 = yes
   ele->e.f2->istapa: version of stab. parameter
      35 = diss wall instationary
      36 = diss wall stationanary
   ele->e.f2->norm_P: p-norm
      p = 1<=p<=oo
      0 = max.-norm (p=oo)
   ele->e.f2->mk: higher order elements control flag
      0 = mk fixed (--> (bi)linear: 1/3, biquadr.: 1/12)
      1 = min(1/3,2*C)
     -1 = mk=1/3  (--> element order via approx. nodal-distance)
   ele->e.f2->ihele[]:
      x/y/z = length-def. for velocity/pressure/continuity stab
      0 = don't compute
      1 = sqrt(area)
      2 = area equivalent diameter
      3 = diameter/sqrt(2)
      4 = sqrt(2)*area/diagonal (rectangle) 4*area/s (triangle)
      5 = streamlength (element length in flow direction
   ele->e.f2->ninths: number of integration points for streamlength
      1 = at center of element
      2 = at every INT pt used for element.-stab.-matrices
   ele->e.f2->istapc: flag for stabilisation parameter calculation
      1 = at center of element
      2 = at every integration point
   ele->e.f2->clamb \
   ele->e.f2->c1     |_>> stabilisation constants (input)
   ele->e.f2->c2     |
   ele->e.f2->c3    /
   ele->e.f2->istrle: has streamlength to be computed
   ele->e.f2->iarea: calculation of area length
   ele->e.f2->iduring: calculation during INT.-pt.loop
   ele->e.f2->itau[0]: flag for tau_mu calc. (-1: before, 1:during)
   ele->e.f2->itau[1]: flag for tau_mp calc. (-1: before, 1:during)
   ele->e.f2->itau[2]: flag for tau_c calc. (-1: before, 1:during)
   ele->e.f2->hk[i]: element sizes (vel / pre / cont)
   ele->e.f2->idiaxy: has diagonals to be computed
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
void f2_calstabpar(
	            ELEMENT         *ele,
		    DOUBLE          *velint,
		    DOUBLE           visc,
		    INT              iel,
		    DIS_TYP          typ,
		    INT              iflag
                  )
{
INT    isp;
DOUBLE hdiv=ONE;
DOUBLE velno;
DOUBLE c_mk;
DOUBLE dt;
DOUBLE re;
DOUBLE hk;
DOUBLE aux1;
STAB_PAR_GLS *gls;	/* pointer to GLS stabilisation parameters	*/

#ifdef DEBUG
dstrc_enter("f2_calstabpar");
#endif

/*---------------------------------------------------------- initialise */
gls  = ele->e.f2->stabi.gls;
fdyn = alldyn[genprob.numff].fdyn;

dsassert(ele->e.f2->stab_type == stab_gls,
   "routine with no or wrong stabilisation called");

/*----------------------- higher order element diameter modifications ? */
switch(gls->mk)
{
case -1:
   c_mk = Q13;
   if (typ==quad8 || typ==quad9)
   {
      dsassert(iel<10,"number of nodes per element not valid!\n");
/*      if (iel<10)
         hdiv = TWO;
      else */
      hdiv = THREE;
   }
   else if (typ==tri6)
   {
      hdiv = TWO;
/*      if (iel==6)
         hdiv = TWO;
      else
         hdiv = THREE; */
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
} /* end swtich (ele->e.f2->mk) */
/*---------------------------------- choose stability-parameter version */
switch(gls->istapa)
{
case 35: /*-------------------------- version diss. Wall - instationary */
   velno = sqrt(velint[0]*velint[0] + velint[1]*velint[1]); /*norm of vel */
   dt = fdyn->dta;     /* check if dta or dt has to be chosen!!!!!!!! */
   for (isp=0;isp<3;isp++)
   {
      if (gls->itau[isp]!=iflag) continue;
      hk = ele->e.f2->hk[isp]/hdiv;
      if(isp==2)
      {
         re = c_mk*hk*velno/TWO/visc;  /* element reynolds number */
         fdyn->tau[isp] = (gls->clamb)*velno*hk/TWO*DMIN(ONE,re);

      }
      else
      {
         if (velno>EPS15)
	 {
	    aux1 = DMIN(hk/TWO/velno , c_mk*hk*hk/FOUR/visc);
            fdyn->tau[isp] = DMIN(dt , aux1);
         }
	 else
            fdyn->tau[isp] = DMIN(dt , c_mk*hk*hk/FOUR/visc);
      } /* end */
   } /* end of loop over isp */
break;

case 36: /*---------------------------- version diss. Wall - stationary */
   dserror("stationary stabilisation not checked yet!!!\n");
   velno = sqrt(velint[0]*velint[0] + velint[1]*velint[1]); /*norm of vel*/
   aux1= velno*c_mk/FOUR/visc;
   for (isp=0;isp<3;isp++)
   {
      if (gls->itau[isp]!=iflag) continue;
      hk = ele->e.f2->hk[isp]/hdiv;
      re = aux1*hk;
      if (isp==2)
         fdyn->tau[isp] = (gls->clamb)*velno*hk/TWO*DMIN(ONE,re);
      else
      {
         if (re<ONE)
	    fdyn->tau[isp] = c_mk*hk*hk/FOUR/visc;
	 else
	    fdyn->tau[isp] = hk/TWO/velno;
      }
   }
break;

default:
   dserror("stability parameter version ISTAP unknown!\n");
} /* end switch (ele->e.f2->istapa) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calstabpar*/

/*!---------------------------------------------------------------------
\brief routine to calculate stability parameter for height function

<pre>                                                         genk 06/03
</pre>

\return void

------------------------------------------------------------------------*/
void f2_calstabpar_hf(
                        ELEMENT          *ele,
                        INT               ngnode,
                        INT              *iedgnod,
                        DOUBLE          **xyze,
                        DOUBLE            det,
                        DOUBLE           *velint,
                        DOUBLE            phiintng,
                        DOUBLE            phiintn,
                        DOUBLE            phiderx
		     )
{
DOUBLE phidot,hs,aux,dx,dy;
DOUBLE velno,rphi;
static DOUBLE C_HF=ONE; /* constant for height function, usually ONE
                           see GUELER, BEHR, TEZDUYAR (1999)             */

#ifdef DEBUG
dstrc_enter("f2_calstabpar_hf");
#endif

dsassert(ngnode==2,"Stab-Parameter at free surface only implemented for lin. elements!\n");
fdyn = alldyn[genprob.numff].fdyn;

/*--------------------------------------- compute hs according to Onate: */
dx = xyze[0][iedgnod[0]]-xyze[0][iedgnod[1]];
dy = xyze[1][iedgnod[0]]-xyze[1][iedgnod[1]];
hs = FABS(dx*velint[0]+dy*velint[1]);


#if 1
/*------------------------------------------------------ compute tau_DC */
if(0)
{
   dserror("heightfun stabilisation");
   velno=sqrt(velint[0]*velint[0] + velint[1]*velint[1]);
#if 0
   switch (ele->e.f2->ihfs[0])
   {
   case 1: /* GUELER, BEHR, TEZDUYAR (1999) */
      hs=TWO*det;
      phidot = (phiintng - phiintn)/fdyn->dta;
      aux  = phidot + velint[0]*phiderx - velint[1];
      rphi = FABS(aux)/fdyn->dphimax;
      fdyn->tau[3] = DMIN(ONE,C_HF*hs*hs*rphi);
   break;
   case 2: /* GUELER, BEHR, TEZDUYAR (1999) */
      velno=sqrt(velint[0]*velint[0] + velint[1]*velint[1]);
      C_HF=TWO*det*velno;
   break;
   default:
      dserror("parameter ihfs[0] out of range!\n");
   }
#endif
}

/*---------------------------------------------------- compute tau_SUPG */
if (0)
{
   dserror("heightfun stabilisation");
#if 0
   switch (ele->e.f2->ihfs[1])
   {
   case 1:
     fdyn->tau[4]= hs/TWO;
     fdyn->tau[4]= 0.01;
   break;
   default:
      dserror("parameter ihfs[0] out of range!\n");
   }
#endif
}

/*------------------------------------------------- print to error file */
#if 1
   fprintf(allfiles.out_err,"tau_dc = %12.10f  tau_supg = %12.10f\n",
          fdyn->tau[3],fdyn->tau[4]);
   fflush(allfiles.out_err);
#endif

#endif

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calstabpar_hf*/

#endif
/*! @} (documentation module close)*/

/*!----------------------------------------------------------------------
\file
\brief submesh stability parameter for fluid2 element

<pre>
Maintainer:  name
             e-mail
             homepage
             telephone number


</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "fluid2ml_prototypes.h"
#include "../fluid2/fluid2.h"

static DOUBLE Q13  = ONE/THREE;
static DOUBLE Q112 = ONE/TWELVE;

/*!--------------------------------------------------------------------- 
\brief routine to calculate submesh stability parameter for fluid2                

<pre>                                                       gravem 07/03   
  									 
In this routine, the stability parameter for the actual submesh element
is calculated.

</pre>

\param   *ele         ELEMENT	      (i)    actual element
\param   *dynvar      FLUID_DYN_CALC  (i)
\param   *mlvar       FLUID_DYN_ML    (i/o)
\param   *velint      DOUBLE	      (i)    l-s velocity at ele. center
\param    visc        DOUBLE	      (i)    viscosity
\param    iel         INT	      (i)    number of nodes	     
\param	  ntyp        INT	      (i)    element type
\return void                                                                       

------------------------------------------------------------------------*/ 
void f2_smstabpar(ELEMENT         *ele,      
		  FLUID_DYN_CALC  *dynvar,
		  FLUID_DYN_ML    *mlvar,
		  DOUBLE	  *velint,  
		  DOUBLE	   visc,    
		  INT		   iel,     
		  INT		   ntyp)
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


/*----------------------- higher order element diameter modifications ? */
switch(mlvar->smstamk)
{
case -1:
   c_mk = Q13;
   if (ntyp==1 && iel>4)
   {
      if (iel<10)
         hdiv = TWO;
      else
         hdiv = THREE;               
   }
   else if (ntyp==2 && iel>3)
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
  dt = dynvar->dta;	/* check if dta or dt has to be chosen!!!!!!!! */
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
  rc = dynvar->thsl;	 /* reaction coefficient = theta times dt */
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
\param	  ntyp        INT	      (i)    element type
\return void                                                                       

------------------------------------------------------------------------*/ 
void f2_smsgvisc(ELEMENT         *ele,
                 FLUID_DYN_ML    *mlvar,
		 DOUBLE 	 *velint,  
		 DOUBLE 	**vderxy,  
		 DOUBLE 	  visc,    
		 INT		  iel,     
		 INT		  ntyp)
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
if (ntyp==1 && iel>4)
{
   if (iel<10)
      hdiv = TWO;
   else
      hdiv = THREE;		  
}
else if (ntyp==2 && iel>3)
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

#endif

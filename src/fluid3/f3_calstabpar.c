/*!----------------------------------------------------------------------
\file
\brief stability parameter for fluid3 element

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

static DOUBLE Q13  = ONE/THREE;
static DOUBLE Q112 = ONE/TWELVE;
static FLUID_DYNAMIC *fdyn;
/*!--------------------------------------------------------------------- 
\brief routine to calculate stability parameter                

<pre>                                                         genk 05/02   
  									 
   ele->e.f3->stabi.gls->iadvec: advection stab.					
      0 = no								
      1 = yes								
   ele->e.f3->stabi.gls->ipres: pressure stab.					
      0 = no								
      1 = yes								
   ele->e.f3->stabi.gls->ivisc: diffusion stab.					
      0 = no								
      1 = GLS-  							
      2 = GLS+  							
   ele->e.f3->stabi.gls->icont: continuity stab.					
      0 = no								
      1 = yes								
   ele->e.f3->stabi.gls->istapa: version of stab. parameter			
      35 = diss wall instationary					
      36 = diss wall stationanary					
   ele->e.f3->stabi.gls->norm_P: p-norm						
      p = 1<=p<=oo							
      0 = max.-norm (p=oo)						
   ele->e.f3->stabi.gls->mk: higher order elements control flag			
      0 = mk fixed (--> (bi)linear: 1/3, biquadr.: 1/12)		
      1 = min(1/3,2*C)  						
     -1 = mk=1/3  (--> element order via approx. nodal-distance)	
   ele->e.f3->stabi.gls->ihele[]:  						
      x/y/z = length-def. for velocity/pressure/continuity stab 	
      0 = don't compute 						
      1 = sqrt(area)							
      2 = area equivalent diameter					
      3 = diameter/sqrt(2)						
      4 = sqrt(2)*area/diagonal (rectangle) 4*area/s (triangle) 	
      5 = streamlength (element length in flow direction		
   ele->e.f3->stabi.gls->ninths: number of integration points for streamlength	
      1 = at center of element  					
      2 = at every INT pt used for element.-stab.-matrices		
   ele->e.f3->stabi.gls->istapc: flag for stabilisation parameter calculation	
      1 = at center of element  					
      2 = at every integration point					
   ele->e.f3->stabi.gls->clamb \							
   ele->e.f3->c1               |_>> stabilisation constants (input)		
   ele->e.f3->c2               |  						
   ele->e.f3->c3              /							
   ele->e.f3->stabi.gls->istrle: has streamlength to be computed			
   ele->e.f3->stabi.gls->iarea: calculation of area length 			
   ele->e.f3->stabi.gls->iduring: calculation during INT.-pt.loop  		
   ele->e.f3-->stabi.gls>itau[0]: flag for tau_mu calc. (-1: before, 1:during)	
   ele->e.f3->stabi.gls->itau[1]: flag for tau_mp calc. (-1: before, 1:during)	
   ele->e.f3->stabi.gls->itau[2]: flag for tau_c calc. (-1: before, 1:during)	
   ele->e.f3->hk[i]: element sizes (vel / pre / cont)			
   ele->e.f3->stabi.gls->idiaxy: has diagonals to be computed			
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
void f3_calstabpar(
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
dstrc_enter("f3_calstabpar");
#endif		

/*---------------------------------------------------------- initialise */
gls   = ele->e.f3->stabi.gls;
fdyn  = alldyn[genprob.numff].fdyn;

dsassert(ele->e.f3->stab_type == stab_gls, 
         "routine with no or wrong stabilisation called");

/*------------------------ higher order element diameter modifications ? */
switch(gls->mk)
{
case -1:
   c_mk = Q13;
   if (typ==hex20 || typ==hex27)
   {
      hdiv=TWO;
/*      if (iel<32)
         hdiv = TWO;
      else
         hdiv = THREE; */ 
   }
   else if (typ==tet10)
   {
      hdiv=TWO;
/*      if (iel==10)
         hdiv = TWO;
      else
         hdiv = THREE; */
   }
break;
case 0:
   if (iel>=8)
      c_mk=Q112;
   else
      c_mk=Q13;
break;   
default:
   c_mk = 0.0;
   dserror("mk > 0 not implemented yet!");
} /* end swtich (ele->e.f3->mk) */
/*---------------------------------- choose stability-parameter version */
switch(gls->istapa)
{
case 35: /*-------------------------- version diss. Wall - instationary */
   velno = sqrt(velint[0]*velint[0] \
               +velint[1]*velint[1] \
	       +velint[2]*velint[2]); /*norm of vel */
   dt = fdyn->dta;     /* check if dta or dt has to be chosen!!!!!!!! */
   for (isp=0;isp<3;isp++)
   {
      if (gls->itau[isp]!=iflag)
         continue;
      hk = ele->e.f3->hk[isp]/hdiv;
      if (isp== 2)/* continiuty stabilisation */
      {   
         re = c_mk*hk*velno/TWO/visc;  /* element reynolds number */
	 fdyn->tau[isp] = (gls->clamb)*velno*hk/TWO*DMIN(ONE,re);        
      }
      else /* velocity / pressure stabilisation */
      {
         if (velno>EPS15)
	 { 
	    aux1 = DMIN(hk/TWO/velno , c_mk*hk*hk/FOUR/visc);
            fdyn->tau[isp] = DMIN(dt , aux1);
         }
	 else
            fdyn->tau[isp] = DMIN(dt , c_mk*hk*hk/FOUR/visc);
      } /* endif(isp) */
   } /* end of loop over isp */
break;
   
case 36: /*---------------------------- version diss. Wall - stationary */
   dserror("stationary stabilisation not checked yet!!!");
   velno = sqrt(velint[0]*velint[0] \
               +velint[1]*velint[1] \
	       +velint[2]*velint[2]); /*norm of vel */
   aux1= velno*c_mk/FOUR/visc;
   for (isp=0;isp<3;isp++)
   {
      if (gls->itau[isp]!=iflag)
         continue;
      hk = ele->e.f3->hk[isp]/hdiv;
      re = aux1*hk;
      if (isp==2) /* continiuty stabilisation ### TWO VERSIONS ??? ###*/
      {
         fdyn->tau[isp] = (gls->clamb)*velno*hk/TWO*DMIN(ONE,re);
/*         dynvar->tau[isp] = velno*hk/TWO*DMIN(ONE,re);*/
      }
      else
      {
         if (re<ONE)
	    fdyn->tau[isp] = c_mk*hk*hk/FOUR/visc;
	 else
	    fdyn->tau[isp] = hk/TWO/velno;
      }  /* end switch (isp) */
   } /* end loop over isp */   
break;

default:
   dserror("stability parameter version ISTAP unknown!");   
} /* end switch (ele->e.f3->istapa) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calstabpar*/	



#endif

/*!----------------------------------------------------------------------
\file
\brief stability parameter for fluid2 element

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
static double Q13  = ONE/THREE;
static double Q112 = ONE/TWELVE;
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
      2 = at every int pt used for element.-stab.-matrices		
   ele->e.f2->istapc: flag for stabilisation parameter calculation	
      1 = at center of element  					
      2 = at every integration point					
   ele->e.f2->clamb \							
   ele->e.f2->c1     |_>> stabilisation constants (input)		
   ele->e.f2->c2     |  						
   ele->e.f2->c3    /							
   ele->e.f2->istrle: has streamlength to be computed			
   ele->e.f2->iarea: calculation of area length 			
   ele->e.f2->iduring: calculation during int.-pt.loop  		
   ele->e.f2->itau[0]: flag for tau_mu calc. (-1: before, 1:during)	
   ele->e.f2->itau[1]: flag for tau_mp calc. (-1: before, 1:during)	
   ele->e.f2->itau[2]: flag for tau_c calc. (-1: before, 1:during)	
   ele->e.f2->hk[i]: element sizes (vel / pre / cont)			
   ele->e.f2->idiaxy: has diagonals to be computed			
   dynvar->tau[0]: stability parameter momentum / velocity (tau_mu)	
   dynvar->tau[1]: stability parameter momentum / pressure (tau_mp)	
   dynvar->tau[2]: stability parameter continuity (tau_c)		

</pre>

\param   *ele,        ELEMENT	      (i)    actual element
\param   *dynvar,     FLUID_DYN_CALC  (i/o)
\param   *velint,     double	      (i)    vel at center
\param    visc,       double	      (i)    viscosity
\param    iel,        int	      (i)    number of nodes	     
\param	  ntyp,       int	      (i)    element type
\param	  iflag       int	      (i)    flag for evaluation
\return void                                                                       

------------------------------------------------------------------------*/ 
void f2_calstabpar(
	            ELEMENT         *ele,      
		    FLUID_DYN_CALC  *dynvar,
		    double          *velint,  
		    double           visc,    
		    int              iel,     
		    int              ntyp,    
		    int              iflag    
                  )
{
int    isp;
double hdiv=ONE; 
double velno;
double c_mk;
double dt;
double re;
double hk;
double aux1;

#ifdef DEBUG 
dstrc_enter("f2_calstabpar");
#endif


/*----------------------- higher order element diameter modifications ? */
switch(ele->e.f2->mk)
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
} /* end swtich (ele->e.f2->mk) */
/*---------------------------------- choose stability-parameter version */
switch(ele->e.f2->istapa)
{
case 35: /*-------------------------- version diss. Wall - instationary */
   velno = sqrt(velint[0]*velint[0] + velint[1]*velint[1]); /*norm of vel */
   dt = dynvar->dta;     /* check if dta or dt has to be chosen!!!!!!!! */
   for (isp=0;isp<3;isp++)
   {
      if (ele->e.f2->itau[isp]!=iflag)
         continue;
      hk = ele->e.f2->hk[isp]/hdiv;
      switch(isp)
      {
      case 2:/* continiuty stabilisation */
         re = c_mk*hk*velno/TWO/visc;  /* element reynolds number */
	 dynvar->tau[isp] = (ele->e.f2->clamb)*velno*hk/TWO*DMIN(ONE,re);         
      break;
      default: /* velocity / pressure stabilisation */
         if (velno>EPS15)
	 { 
	    aux1 = DMIN(hk/TWO/velno , c_mk*hk*hk/FOUR/visc);
            dynvar->tau[isp] = DMIN(dt , aux1);
         }
	 else
            dynvar->tau[isp] = DMIN(dt , c_mk*hk*hk/FOUR/visc);
       break;
      } /* end switch (isp) */
   } /* end of loop over isp */
break;
   
case 36: /*---------------------------- version diss. Wall - stationary */
   dserror("stationary stabilisation not checked yet!!!\n");
   velno = sqrt(velint[0]*velint[0] + velint[1]*velint[1]); /*norm of vel*/
   aux1= velno*c_mk/FOUR/visc;
   for (isp=0;isp<3;isp++)
   {
      if (ele->e.f2->itau[isp]!=iflag)
         continue;
      hk = ele->e.f2->hk[isp]/hdiv;
      re = aux1*hk;
      switch(isp)
      {
      case 2: /* continiuty stabilisation ### TWO VERSIONS ??? ###*/
         dynvar->tau[isp] = (ele->e.f2->clamb)*velno*hk/TWO*DMIN(ONE,re);
/*        dynvar->tau[isp] = velno*hk/TWO*DMIN(ONE,re); */
      break;
      default: /* velocity / pressure stabilisation */
         if (re<ONE)
	    dynvar->tau[isp] = c_mk*hk*hk/FOUR/visc;
	 else
	    dynvar->tau[isp] = hk/TWO/velno;
      break;
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

#endif

/*!----------------------------------------------------------------------
\file
\brief stability parameter for fluid2 element

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
#include "fluid2_tu.h"
/*!--------------------------------------------------------------------- 
\brief routine to calculate stability parameter                

<pre>                                                           he  02/03   
  									 

</pre>

\param   *ele,        ELEMENT	      (i)    actual element
\param   *elev,       ELEMENT	      (i)    element with velocity-infos
\param   *dynvar,     FLUID_DYN_CALC(i/o)  
\param   *eddyint,    double	      (i)    eddyint at center
\param   *velint,     double	      (i)    vel at center
\param   *velint_dc,  double	      (i)    vel at center for D.C.
\param    visc,       double	      (i)    viscosity
\return void                                                                       

------------------------------------------------------------------------*/ 
void f2_calstabpar_tu_1(
	            ELEMENT         *ele,      
		      ELEMENT         *elev,
                  FLUID_DYN_CALC  *dynvar,
		      double           eddyint, 
                  double          *velint, 
                  double          *velint_dc, 
                  double           visc    
                  )
{ 
double peclet,xi_Pe; 
double vel_abs;
double kappa; 

#ifdef DEBUG 
dstrc_enter("f2_calstabpar_tu_1");
#endif


/*----------------------------------- calculate peclet number           */  
ele->e.f2_tu->strom = elev->e.f2->hk[0];

vel_abs = sqrt(pow(velint[0],2) + pow(velint[1],2));      /*norm of vel */

if (dynvar->kapomega_flag == 0) kappa = visc+eddyint*0.5;
if (dynvar->kapomega_flag == 1) kappa = visc+eddyint*0.5;

peclet  = ele->e.f2_tu->strom * vel_abs / (2*kappa);

xi_Pe = DMIN(peclet/3, 1);

/*----------------------------------- calculate stabilisation parameter */  
if(vel_abs!=0.0)
dynvar->tau_tu      = (ele->e.f2_tu->strom*xi_Pe) / (2*vel_abs); 
else
dynvar->tau_tu      = dynvar->dta; 
  
dynvar->tau_tu      = DMIN(dynvar->tau_tu,dynvar->dta);   

/*------------------- calculate stabilisation parameter for DISC. CAPT. */
vel_abs = sqrt(pow(velint_dc[0],2) + pow(velint_dc[1],2));      

peclet  = ele->e.f2_tu->strom_dc * vel_abs / (2*kappa);

xi_Pe  = DMIN(peclet/3, 1);

if(vel_abs != 0.0)
dynvar->tau_tu_dc  = (ele->e.f2_tu->strom_dc*xi_Pe) / (2*vel_abs); 
else
dynvar->tau_tu_dc  = dynvar->dta; 
  
/*------------------------------------- VERSION: DC1 -------------------*/                  
dynvar->tau_tu_dc  = DMIN(dynvar->tau_tu_dc,dynvar->dta);  

/*------------------------------------- VERSION: DC2 -------------------*/                  
/*dynvar->tau_tu_dc  = DMAX(0,dynvar->tau_tu_dc-dynvar->tau_tu);*/  

/*----------------------------------------------------------------------*/

#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calstabpar_tu_1*/	

#endif

/*!----------------------------------------------------------------------
\file
\brief stability parameter for fluid2 element

<pre>
Maintainer: Thomas Hettich
            hettich@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hettich/
            0771 - 685-6575
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID2TU 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
#include "fluid2_tu.h"
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
\brief routine to calculate stability parameter                

<pre>                                                           he  12/02   
  									 

</pre>

\param   *ele,        ELEMENT	      (i)    actual element
\param   *elev,       ELEMENT	      (i)    element with velocity-infos
\param   *eddyint,    DOUBLE	      (i)    eddyint at center
\param   *velint,     DOUBLE	      (i)    vel at center
\param   *velint_dc,  DOUBLE	      (i)    vel at center for D.C.
\param    visc,       DOUBLE	      (i)    viscosity
\return void                                                                       

------------------------------------------------------------------------*/ 
void f2_calstabpar_tu(
	            ELEMENT         *ele,      
		      ELEMENT         *elev,
		      DOUBLE           eddyint, 
                  DOUBLE          *velint, 
                  DOUBLE          *velint_dc, 
                  DOUBLE           visc    
                  )
{ 
DOUBLE peclet,xi_Pe; 
DOUBLE vel_abs;
DOUBLE kappa; 

#ifdef DEBUG 
dstrc_enter("f2_calstabpar_tu");
#endif


/*--------------------------------------------- calculate peclet number */
fdyn   = alldyn[genprob.numff].fdyn;  
ele->e.f2_tu->strom = elev->e.f2->hk[0];

vel_abs = sqrt(pow(velint[0],2) + pow(velint[1],2));      

if (fdyn->kapeps_flag == 0) kappa = visc+eddyint;
if (fdyn->kapeps_flag == 1) kappa = visc+eddyint/1.3;

peclet  = ele->e.f2_tu->strom * vel_abs / (2*kappa);

xi_Pe  = DMIN(peclet/3, 1);

/*----------------------------------- calculate stabilisation parameter */  
if(vel_abs!=0)
fdyn->tau_tu      = (ele->e.f2_tu->strom*xi_Pe) / (2*vel_abs); 
else
fdyn->tau_tu      = fdyn->dta; 
  
fdyn->tau_tu      = DMIN(fdyn->tau_tu,fdyn->dta);   


/*------------------- calculate stabilisation parameter for DISC. CAPT. */  
vel_abs = sqrt(pow(velint_dc[0],2) + pow(velint_dc[1],2));      

peclet  = ele->e.f2_tu->strom_dc * vel_abs / (2*kappa);

xi_Pe  = DMIN(peclet/3, 1);

if(vel_abs != 0.0)
fdyn->tau_tu_dc   = (ele->e.f2_tu->strom_dc*xi_Pe) / (2*vel_abs); 
else
fdyn->tau_tu_dc   = fdyn->dta; 
  
/*------------------------------------- VERSION: DC1 -------------------*/                  
fdyn->tau_tu_dc   = DMIN(fdyn->tau_tu_dc,fdyn->dta);  

/*------------------------------------- VERSION: DC2 -------------------*/                  
/*fdyn->tau_tu_dc   = DMAX(0,fdyn->tau_tu_dc-fdyn->tau_tu);*/  

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calstabpar_tu*/	

#endif

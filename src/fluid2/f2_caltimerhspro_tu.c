/*!----------------------------------------------------------------------
\file
\brief time RHS for fluid2 element

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2_tu.h"
/*!--------------------------------------------------------------------- 
\brief galerkin part of time forces for kapeps dof

<pre>                                                         he  12/02

In this routine the galerkin part of the time forces for kapeps dof
is calculated:

PRODUKTIONSTERM:
                     /
   (+)    THETA*dt  |  0.5 * factor1* eddyint * ((grad(u) + [grad(u)]^T)^2 * psi d_omega
                   /		  		      


</pre>
\param   *dynvar      FLUID_DYN_CALC  (i)
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param    eddynint     DOUBLE	      (i)    eddy-visc. at integr. point
\param   *funct       DOUBLE	      (i)    nat. shape functions      
\param    visc	    DOUBLE	      (i)    fluid viscosity
\param    fac	    DOUBLE	      (i)    weighting factor
\param    factor1	    DOUBLE	      (i)    factor
\param    production  DOUBLE	      (i)    factor
\param    iel	    INT	      (i)    num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calgalprofkapeps(
                  FLUID_DYN_CALC  *dynvar, 
                  DOUBLE          *eforce,    
		      DOUBLE           eddynint,
                  DOUBLE          *funct,    
                  DOUBLE           visc,     
		      DOUBLE           fac,      
                  DOUBLE           factor1,  
                  DOUBLE           production,  
                  INT              iel       
                  )  
{
INT    j,irow,isd,inode;  
DOUBLE c;
DOUBLE aux;
DOUBLE facsl;

#ifdef DEBUG 
dstrc_enter("f2_calgaltfkapeps");
#endif


/*--------------------------------------------------- set some factors */
facsl = fac * dynvar->thsl;

/*----------------------------------------------------------------------*
   Calculate forces of time force vector:
                 /
   (+) THETA*dt |   0.5 * nue_t * factor1 * (grad(u)+[grad(u)^T])^2 * psi   d_omega
               /   
*----------------------------------------------------------------------*/
   irow = 0;
   for (inode=0;inode<iel;inode++)
   {
      eforce[irow]   += eddynint * factor1 * production * funct[inode]  * facsl;
      irow ++;
   } /* end loop over inode */


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_calele */

/*!--------------------------------------------------------------------- 
\brief stabilisation part of time forces for kapeps dof

<pre>                                                         he  12/02

In this routine the stabilisation part of the time forces for kapeps dofs
is calculated:
		 		                   	  		      

PRODUKTIONSTERM:

              /
(+) THETA*dt |  tau_tu * 0.5 * nue_t * factor1 * (grad(u)+[grad(u)^T])^2 * grad(psi) * u   d_omega + D. C.
            /   
              


       
</pre>
\param   *dynvar      FLUID_DYN_CALC(i)
\param   *ele         ELEMENT	      (i)    actual element
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param    eddynint    DOUBLE	      (i)    eddy-visc. at integr. point
\param   *funct       DOUBLE	      (i)    nat. shape functions      
\param    visc	    DOUBLE	      (i)    fluid viscosity
\param    fac	    DOUBLE	      (i)    weighting factor
\param    factor1	    DOUBLE	      (i)    factor
\param    production  DOUBLE	      (i)    factor
\param   *velint      DOUBLE	      (i)    vel. at integr. point
\param   *velint_dc   DOUBLE	      (i)    vel. at integr. point for D.C.
\param   **derxy      DOUBLE	      (i)    global deriv.
\param    iel	    INT	      (i)    num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabprofkapeps(
                   FLUID_DYN_CALC  *dynvar,
                   ELEMENT         *ele,      
	            DOUBLE          *eforce,  
		      DOUBLE           eddynint, 
                  DOUBLE          *funct,    
                  DOUBLE           visc,     
                  DOUBLE           fac,     
                  DOUBLE           factor1,
                  DOUBLE           production,  
                  DOUBLE           *velint,  
                  DOUBLE           *velint_dc,  
                  DOUBLE          **derxy,  
                  INT              iel      
                  )
{
INT    j,irow,isd,inode;
DOUBLE aux;
DOUBLE taumu,taumu_dc;
DOUBLE facsl,facsl_dc;

#ifdef DEBUG 
dstrc_enter("f2_calstabtfkapeps");
#endif

/*--------------------------------------------------- set some factors */
taumu    = dynvar->tau_tu;
taumu_dc = dynvar->tau_tu_dc;
facsl    = fac * dynvar->thsl * taumu;
facsl_dc = fac * dynvar->thsl * taumu_dc;

/*----------------------------------------------------------------------*
   Calculate forces of time force vector:
                 /
   (+) THETA*dt |  tau_tu * 0.5 * nue_t * factor1 * (grad(u)+[grad(u)^T])^2 * grad(psi) * u   d_omega
               /   
*----------------------------------------------------------------------*/
   irow = 0;
   for (inode=0;inode<iel;inode++)
   {
    aux = eddynint*factor1*production;
    
      for (isd=0;isd<2;isd++)
      {
        eforce[irow]   += aux*facsl   *velint[isd]   *derxy[isd][inode];
        eforce[irow]   += aux*facsl_dc*velint_dc[isd]*derxy[isd][inode]; 
      } /* end loop over isd */
    irow ++;
   } /* end loop over inode */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_calele */

#endif

/*!----------------------------------------------------------------------
\file
\brief time RHS for fluid2 element

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

static FLUID_DYNAMIC *fdyn;
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
fdyn  = alldyn[genprob.numff].fdyn;
facsl = fac * fdyn->thsl;

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
fdyn     = alldyn[genprob.numff].fdyn;

taumu    = fdyn->tau_tu;
taumu_dc = fdyn->tau_tu_dc;
facsl    = fac * fdyn->thsl * taumu;
facsl_dc = fac * fdyn->thsl * taumu_dc;

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

/*!----------------------------------------------------------------------
\file
\brief time RHS for fluid2 element

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2_tu.h"
/*!--------------------------------------------------------------------- 
\brief galerkin part of time forces for kapome dof

<pre>                                                         he  02/03

In this routine the galerkin part of the time forces for kapome dof
is calculated:

PRODUKTIONSTERM:
                     /
   (+)    THETA*dt  |  0.5 * factor1* eddyint * ((grad(u) + [grad(u)]^T)^2 * psi d_omega
                   /		  		      


</pre>
\param   *dynvar      FLUID_DYN_CALC  (i)
\param   *eforce      double	      (i/o)  element force vector
\param    eddynint     double	      (i)    eddy-visc. at integr. point
\param   *funct       double	      (i)    nat. shape functions      
\param    fac	    double	      (i)    weighting factor
\param    factor1	    double	      (i)    factor
\param    production  double	      (i)    factor
\param    iel	    int	      (i)    num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calgalprofkapome(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,    
		      double           eddynint,
                  double          *funct,    
		      double           fac,      
                  double           factor1,  
                  double           production,  
                  int              iel       
                  )  
{
int    j,irow,isd,inode;  
double c;
double aux;
double facsl;

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
\brief stabilisation part of time forces for kapome dof

<pre>                                                         he  02/03

In this routine the stabilisation part of the time forces for kapome dofs
is calculated:
		 		                   	  		      

PRODUKTIONSTERM:

              /
(+) THETA*dt |  tau_tu * 0.5 * nue_t * factor1 * (grad(u)+[grad(u)^T])^2 * grad(psi) * u   d_omega + D. C.
            /   
              


       
</pre>
\param   *dynvar      FLUID_DYN_CALC(i)
\param   *ele         ELEMENT	      (i)    actual element
\param   *eforce      double	      (i/o)  element force vector
\param    eddynint    double	      (i)    eddy-visc. at integr. point
\param   *funct       double	      (i)    nat. shape functions      
\param    fac	    double	      (i)    weighting factor
\param    factor1	    double	      (i)    factor
\param    production  double	      (i)    factor
\param   *velint      double	      (i)    vel. at integr. point
\param   *velint_dc   double	      (i)    vel. at integr. point for D.C.
\param   **derxy      double	      (i)    global deriv.
\param    iel	    int	      (i)    num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabprofkapome(
                   FLUID_DYN_CALC  *dynvar,
                   ELEMENT         *ele,      
	            double          *eforce,  
		      double           eddynint, 
                  double          *funct,    
                  double           fac,     
                  double           factor1,
                  double           production,  
                  double           *velint,  
                  double           *velint_dc,  
                  double          **derxy,  
                  int              iel      
                  )
{
int    j,irow,isd,inode;
double aux;
double taumu,taumu_dc;
double facsl,facsl_dc;

#ifdef DEBUG 
dstrc_enter("f2_calstabtfkapome");
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

/*!----------------------------------------------------------------------
\file
\brief iteration RHS for fluid2 element

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2_tu.h"
/*!--------------------------------------------------------------------- 
\brief galerkin part of iteration forces for kapome dofs

<pre>                                                        he  03/02

In this routine the galerkin part of the iteration forces for kapome dofs
is calculated:

                    / 
          THETA*dt | factor * (kapome_old)^2 * psi  d_omega
                  /  
               

</pre>
\param  *dynvar      FLUID_DYN_CALC  (i)
\param  *eforce      double	    (i/o)   element force vector
\param   kapomeint   double	     (i)    kapome at integr. point
\param  *funct       double	     (i)    nat. shape funcs
\param   fac 	   double	     (i)    weighting factor
\param   factor2 	   double	     (i)    factor
\param   iel	   int           (i)	num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calgalifkapome(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,
                  double           kapomeint,  
                  double          *funct,   
		      double           fac,     
                  double           factor2,    
                  int              iel      
                 )  
{
int    inode,isd;
int    irow;  
double facsl;

#ifdef DEBUG 
dstrc_enter("f2_calgalifkapome");
#endif

/*----------------------------------------------- set some factors */
facsl = fac * dynvar->thsl;

/*------------------------------------------------------------------*
   Calculate forces of iteration force vector:
             /
   THETA*dt |  factor * (kapome_old)^2 * psi  d_omega
           /  
 *------------------------------------------------------------------*/ 
irow = -1;
for (inode=0;inode<iel;inode++)
{
      irow++;
      eforce[irow] += factor2 * pow(kapomeint,2) * funct[inode] * facsl;
} /* end loop over inode */


/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calgalifv */

/*!---------------------------------------------------------------------  
\brief stabilisation part of iteration forces for kapome dofs

<pre>                                                         he  03/02

In this routine the stabilisation part of the iteration forces for kapome


           /
 THETA*dt | tau_tu * factor * (kapome_old)^2 * grad(psi) * u  d_omega + D.C.
         /



</pre>
\param   *dynvar   FLUID_DYN_CALC  (i)
\param   *ele      ELEMENT	   (i)    actual element
\param   *eforce   double	   (i/o)  element force vector
\param    kapomeint double	   (i)    kapome at integr. point
\param   *velint   double	   (i)    vel at integr. point
\param   *velint_dc  double	   (i)    vel at integr. point for D.C.
\param   *funct    double	   (i)    nat. shape funcs
\param  **derxy    double	   (i)    global derivative
\param    fac 	 double	   (i)    weighting factor
\param    factor2  double	   (i)    factor
\param    iel	   int	   (i)    num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabifkapome(
                  FLUID_DYN_CALC  *dynvar, 
                  ELEMENT         *ele,      
	            double          *eforce,  
		      double           kapomeint,  
                  double          *velint,  
                  double          *velint_dc,  
                  double          *funct,   
		      double         **derxy,   
		      double           fac,     
                  double           factor2,    
		      int              iel      
                  )  
{
int    inode,isd;
int    irow;  
double facsl,facsl_dc;
double aux,aux_dc;
double taumu,taumu_dc;

#ifdef DEBUG 
dstrc_enter("f2_calstabifkapome");
#endif

/*--------------------------------------------------- set some factors */
taumu    = dynvar->tau_tu;
taumu_dc = dynvar->tau_tu_dc;
facsl    = fac * dynvar->thsl * taumu;
facsl_dc = fac * dynvar->thsl * taumu_dc;
/*---------------------------------------------------------------------- 
   Calculate  stabilastion of iteration force vector:

           /
 THETA*dt |tau_tu * factor * (kapome_old)^2 * grad(psi) * u  d_omega
         /
   
*----------------------------------------------------------------------*/ 
 irow=0;
   for (inode=0;inode<iel;inode++)
   {
      aux    = derxy[0][inode]*velint[0]    + derxy[1][inode]*velint[1];
      aux_dc = derxy[0][inode]*velint_dc[0] + derxy[1][inode]*velint_dc[1];
     
      eforce[irow] += factor2*pow(kapomeint,2)*aux   *facsl;
      eforce[irow] += factor2*pow(kapomeint,2)*aux_dc*facsl_dc;  
     
      irow ++;		      
   } /* end loop over inode */


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calgalifkapeps */

#endif

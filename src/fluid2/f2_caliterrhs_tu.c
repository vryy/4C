/*!----------------------------------------------------------------------
\file
\brief iteration RHS for fluid2 element

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2_tu.h"
/*!--------------------------------------------------------------------- 
\brief galerkin part of iteration forces for kapeps dofs

<pre>                                                        he  12/02

In this routine the galerkin part of the iteration forces for kapeps dofs
is calculated:

                    / 
          THETA*dt | factor * (kapeps_old)^2 * psi  d_omega
                  /  
               

LOW-REYNOLD's MODEL only for epsilon:
    
                   / 
   (+)   THETA*dt | 2.0*visc*nue_t* (vderxy2_12)^2 * psi  d_omega
                 /  


</pre>
\param  *dynvar      FLUID_DYN_CALC  (i)
\param  *eforce      double	    (i/o)   element force vector
\param   eddyint     double	     (i)    eddy-visc at integr. point
\param   kapepsint   double	     (i)    kapeps at integr. point
\param  *funct       double	     (i)    nat. shape funcs
\param   fac 	   double	     (i)    weighting factor
\param   factor2 	   double	     (i)    factor
\param   vderxy_12   double	     (i)    factor
\param   visc 	   double	     (i)    viscosity
\param   iel	   int           (i)	num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calgalifkapeps(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,
		      double           eddyint,  
                  double           kapepsint,  
                  double          *funct,   
		      double           fac,     
                  double           factor2,    
                  double           vderxy_12,    
                  double           visc,    
                  int              iel      
                 )  
{
int    inode,isd;
int    irow;  
double facsl;

#ifdef DEBUG 
dstrc_enter("f2_calgalifkapeps");
#endif

/*----------------------------------------------- set some factors */
facsl = fac * dynvar->thsl;

/*------------------------------------------------------------------*
   Calculate forces of iteration force vector:
             /
   THETA*dt |  factor * (kapeps_old)^2 * psi  d_omega
           /  
 *------------------------------------------------------------------*/ 
irow = -1;
for (inode=0;inode<iel;inode++)
{
      irow++;
      eforce[irow] += factor2 * pow(kapepsint,2) * funct[inode] * facsl;
} /* end loop over inode */


if(dynvar->kapeps_flag==1) 
{
/*------------------------------------------------------------------*
   Calculate forces of iteration force vector:

LOW-REYNOLD's MODEL:

                    / 
          THETA*dt | 2.0*visc* nue_t* (vderxy_12)^2 * psi   d_omega
                  /  
*------------------------------------------------------------------*/ 
irow = -1;
for (inode=0;inode<iel;inode++)
{
    irow++;
    eforce[irow] += 2.0*visc*eddyint* vderxy_12 * funct[inode] * facsl;
} /* end loop over inode */

} /* endif */

/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calgalifv */

/*!---------------------------------------------------------------------  
\brief stabilisation part of iteration forces for kapeps dofs

<pre>                                                         he  12/02

In this routine the stabilisation part of the iteration forces for kapeps


           /
 THETA*dt | tau_tu * factor * (kapeps_old)^2 * grad(psi) * u  d_omega + D. C.
         /


LOW-REYNOLD's MODEL only for epsilon:


              / 
(+) THETA*dt | tau_tu * 2.0*visc*nue_t* (vderxy2_12)^2 *  grad(psi) * u  d_omega + D. C.
            /  


</pre>
\param   *dynvar   FLUID_DYN_CALC  (i)
\param   *ele      ELEMENT	   (i)    actual element
\param   *eforce   double	   (i/o)  element force vector
\param    kapepsint double	   (i)    kapeps at integr. point
\param   *velint   double	   (i)    vel at integr. point
\param   *velint_dc  double	   (i)    vel at integr. point for D.C.
\param    eddyint  double	   (i)    eddy-visc. at integr. point
\param   *funct    double	   (i)    nat. shape funcs
\param  **derxy    double	   (i)    global derivative
\param    fac 	 double	   (i)    weighting factor
\param    factor2  double	   (i)    factor
\param    vderxy_12 double	   (i)    factor
\param    visc     double	   (i)    fluid viscosity
\param    iel	   int	   (i)    num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabifkapeps(
                  FLUID_DYN_CALC  *dynvar, 
                  ELEMENT         *ele,      
	            double          *eforce,  
		      double           kapepsint,  
                  double          *velint,  
                  double          *velint_dc,  
		      double           eddyint,  
                  double          *funct,   
		      double         **derxy,   
		      double           fac,     
                  double           factor2,    
                  double           vderxy_12,    
                  double           visc,    
		      int              iel      
                  )  
{
int    inode,isd;
int    irow;  
double facsl,facsl_dc;
double aux,aux_dc;
double taumu,taumu_dc;

#ifdef DEBUG 
dstrc_enter("f2_calstabifkapeps");
#endif

/*--------------------------------------------------- set some factors */
taumu    = dynvar->tau_tu;
taumu_dc = dynvar->tau_tu_dc;
facsl    = fac * dynvar->thsl * taumu;
facsl_dc = fac * dynvar->thsl * taumu_dc;
/*---------------------------------------------------------------------- 
   Calculate  stabilastion of iteration force vector:

           /
 THETA*dt |tau_tu * factor * (kapeps_old)^2 * grad(psi) * u  d_omega
         /
   
*----------------------------------------------------------------------*/ 
 irow=0;
   for (inode=0;inode<iel;inode++)
   {
      aux    = derxy[0][inode]*velint[0]    + derxy[1][inode]*velint[1];
      aux_dc = derxy[0][inode]*velint_dc[0] + derxy[1][inode]*velint_dc[1];
     
      eforce[irow] += factor2*pow(kapepsint,2)*aux   *facsl;
      eforce[irow] += factor2*pow(kapepsint,2)*aux_dc*facsl_dc; 
     
      irow ++;		      
   } /* end loop over inode */


if(dynvar->kapeps_flag==1) 
{
/*---------------------------------------------------------------------- 
   Calculate  stabilastion of iteration force vector:

LOW-REYNOLD's MODEL:

             / 
   THETA*dt | tau_tu * 2.0*visc* nue_t* (vderxy2_12)^2 * grad(psi) * u  d_omega
           /  

*----------------------------------------------------------------------*/ 
 irow=0;
   for (inode=0;inode<iel;inode++)
   {
     aux = 2.0*visc*eddyint*vderxy_12;

      for (isd=0;isd<2;isd++)
      {
        eforce[irow]   += aux*facsl   *velint[isd]   *derxy[isd][inode];
        eforce[irow]   += aux*facsl_dc*velint_dc[isd]*derxy[isd][inode]; 
      } /* end loop over isd */

     irow ++;
   } /* end loop over inode */

} /* endif */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calgalifkapeps */

#endif

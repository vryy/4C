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
\param  *eforce      DOUBLE	    (i/o)   element force vector
\param   eddyint     DOUBLE	     (i)    eddy-visc at integr. point
\param   kapepsint   DOUBLE	     (i)    kapeps at integr. point
\param  *funct       DOUBLE	     (i)    nat. shape funcs
\param   fac 	   DOUBLE	     (i)    weighting factor
\param   factor2 	   DOUBLE	     (i)    factor
\param   vderxy_12   DOUBLE	     (i)    factor
\param   visc 	   DOUBLE	     (i)    viscosity
\param   iel	   INT           (i)	num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calgalifkapeps(
                  FLUID_DYN_CALC  *dynvar, 
                  DOUBLE          *eforce,
		      DOUBLE           eddyint,  
                  DOUBLE           kapepsint,  
                  DOUBLE          *funct,   
		      DOUBLE           fac,     
                  DOUBLE           factor2,    
                  DOUBLE           vderxy_12,    
                  DOUBLE           visc,    
                  INT              iel      
                 )  
{
INT    inode,isd;
INT    irow;  
DOUBLE facsl;

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
\param   *eforce   DOUBLE	   (i/o)  element force vector
\param    kapepsint DOUBLE	   (i)    kapeps at integr. point
\param   *velint   DOUBLE	   (i)    vel at integr. point
\param   *velint_dc  DOUBLE	   (i)    vel at integr. point for D.C.
\param    eddyint  DOUBLE	   (i)    eddy-visc. at integr. point
\param   *funct    DOUBLE	   (i)    nat. shape funcs
\param  **derxy    DOUBLE	   (i)    global derivative
\param    fac 	 DOUBLE	   (i)    weighting factor
\param    factor2  DOUBLE	   (i)    factor
\param    vderxy_12 DOUBLE	   (i)    factor
\param    visc     DOUBLE	   (i)    fluid viscosity
\param    iel	   INT	   (i)    num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabifkapeps(
                  FLUID_DYN_CALC  *dynvar, 
                  ELEMENT         *ele,      
	            DOUBLE          *eforce,  
		      DOUBLE           kapepsint,  
                  DOUBLE          *velint,  
                  DOUBLE          *velint_dc,  
		      DOUBLE           eddyint,  
                  DOUBLE          *funct,   
		      DOUBLE         **derxy,   
		      DOUBLE           fac,     
                  DOUBLE           factor2,    
                  DOUBLE           vderxy_12,    
                  DOUBLE           visc,    
		      INT              iel      
                  )  
{
INT    inode,isd;
INT    irow;  
DOUBLE facsl,facsl_dc;
DOUBLE aux,aux_dc;
DOUBLE taumu,taumu_dc;

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

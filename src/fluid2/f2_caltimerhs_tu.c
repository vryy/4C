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
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2_tu.h"
/*!--------------------------------------------------------------------- 
\brief galerkin part of time forces for kapeps dof

<pre>                                                         he  12/02

In this routine the galerkin part of the time forces for kapeps dof
is calculated:

        /
   (+) |  kapeps * psi     d_omega
      /  
               

                      /
   (-) (1-THETA)*dt  |  u * grad(kapeps) * psi     d_omega
                    /
		  
                      /
   (-) (1-THETA)*dt  |  (nue+nue_t/sig) *grad(kapeps) * grad(psi)  d_omega
                    /  
		  
                      /
   (-) (1-THETA)*dt  |  factor * kapeps^2 * psi d_omega
                    /		  		      

LOW-REYNOLD's MODEL only for kappa:

                       /
   (-) (1-THETA)*dt   | 2.0 * visc * grad(k_old)*grad(k_old)/(4*k_old)  * psi  d_omega
                     /
 

                      /
   (+) (1-THETA)*dt  |  0.5 * factor1* eddyint * ((grad(u) + [grad(u)]^T)^2 * psi d_omega
                    /		  		      

                     / 
   (+) (1-THETA)*dt | factor2 * (kapeps)^2 * psi  d_omega
                   /  

LOW-REYNOLD's MODEL only for epsilon:

                     / 
   (+) (1-THETA)*dt | 2.0* visc * nue_t * (vderxy2_12)^2 *  psi  d_omega
                   /

</pre>
\param   *dynvar      FLUID_DYN_CALC  (i)
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param    kapepsint   DOUBLE	      (i)    vel. at integr. point
\param   *velint      DOUBLE	      (i)    vel. at integr. point
\param    eddyint     DOUBLE	      (i)    eddy-visc. at integr. point
\param   *funct       DOUBLE	      (i)    nat. shape functions      
\param  **derxy       DOUBLE	      (i)    global derivatives
\param  **vderxy      DOUBLE	      (i)    global vel. deriv.
\param   *kapepsderxy DOUBLE	      (i)    global kapeps deriv.
\param    visc	    DOUBLE	      (i)    fluid viscosity
\param    fac	    DOUBLE	      (i)    weighting factor
\param    factor	    DOUBLE	      (i)    factor
\param    factor1	    DOUBLE	      (i)    factor
\param    factor2	    DOUBLE	      (i)    factor
\param    sig	    DOUBLE	      (i)    factor
\param    vderxy_12   DOUBLE	      (i)    factor
\param    production  DOUBLE	      (i)    factor
\param    iel	    INT	      (i)    num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calgaltfkapeps(
                  FLUID_DYN_CALC  *dynvar, 
                  DOUBLE          *eforce,    
		      DOUBLE           kapepsint,  
                  DOUBLE          *velint,   
		      DOUBLE           eddyint,
                  DOUBLE          *funct,    
		      DOUBLE         **derxy,    
		      DOUBLE         **vderxy,   
		      DOUBLE          *kapepsderxy,   
                  DOUBLE           visc,     
		      DOUBLE           fac,      
                  DOUBLE           factor,  
                  DOUBLE           factor1,  
                  DOUBLE           factor2,  
                  DOUBLE           sig,  
                  DOUBLE           vderxy_12,  
                  DOUBLE           production,  
                  INT              iel       
                  )  
{
INT    j,irow,isd,inode;  
DOUBLE c;
DOUBLE aux;
DOUBLE facsr;

#ifdef DEBUG 
dstrc_enter("f2_calgaltfkapeps");
#endif


/*--------------------------------------------------- set some factors */
facsr = fac * dynvar->thsr;

/*----------------------------------------------------------------------*
   Calculate intertia forces of time force vector:
      /
     |  kapeps * psi     d_omega
    /  
 *----------------------------------------------------------------------*/ 
irow=-1;
for (inode=0;inode<iel;inode++)
{
      irow++;
      eforce[irow] += funct[inode]*kapepsint*fac;
} /* end of loop over inode */

/*----------------------------------------------------------------------*
   Calculate  forces of time force vector:
    
                      /
   (-) (1-THETA)*dt  |  u * grad(kapeps) * psi     d_omega
                    /
*----------------------------------------------------------------------*/
irow=0;
for (inode=0;inode<iel;inode++)
{
 aux = funct[inode]*facsr;

   for (isd=0;isd<2;isd++)
   {
      eforce[irow] -= aux*kapepsderxy[isd]*velint[isd];
   } /* end of loop over isd */
 irow++;
} /* end of loop over irwo */

/*----------------------------------------------------------------------*
   Calculate  forces of time force vector:
                      /
   (-) (1-THETA)*dt  |  (nue+nue_t/sig) *grad(kapeps) * grad(psi)  d_omega
                    /  
 *----------------------------------------------------------------------*/ 
irow=0;
for (inode=0;inode<iel;inode++)
{
 aux = (visc+eddyint/sig)*facsr;

   for (isd=0;isd<2;isd++)
   {
      eforce[irow] -= aux*derxy[isd][inode]*kapepsderxy[isd];
   } /* end of loop over isd */
 irow++;
} /* end of loop over inode */

/*----------------------------------------------------------------------*
   Calculate  forces of time force vector:
                      /
   (-) (1-THETA)*dt  |  factor * kapeps^2 * psi d_omega
                    /		  		      
 *----------------------------------------------------------------------*/ 
  irow=-1;
   for (inode=0;inode<iel;inode++)
   {
      irow++;
      eforce[irow] -= factor*funct[inode]*pow(kapepsint,2)*facsr;
   } /* end of loop over inode */

if(dynvar->kapeps_flag==0) 
{
/*----------------------------------------------------------------------*
   Calculate  forces of time force vector:

LOW-REYNOLD's MODEL:

                      /
  (-) (1-THETA)*dt   | 2.0 * visc * grad(k_old)*grad(k_old)/(4*k_old)  * psi  d_omega
                    /
*----------------------------------------------------------------------*/ 
  irow=-1;
   for (inode=0;inode<iel;inode++)
   {
      irow++;
      eforce[irow] -= 1/(4*kapepsint)*2.0*visc*
                     (pow(kapepsderxy[0],2)+pow(kapepsderxy[1],2))*funct[inode]*facsr;
   } /* end of loop over inode */

} /* endif */
/*----------------------------------------------------------------------*
   Calculate  forces of time force vector:
                      /
   (+) (1-THETA)*dt  |  0.5 * factor1 * eddyint* ((grad(u) + [grad(u)]^T)^2 * psi d_omega
                    /		  		      
 *----------------------------------------------------------------------*/ 
irow=-1;
for (inode=0;inode<iel;inode++)
{
      irow++;
      eforce[irow] += factor1*eddyint*production*funct[inode]*facsr;
} /* end of loop over inode */

/*----------------------------------------------------------------------*
   Calculate  forces of time force vector:
                      /
   (+) (1-THETA)*dt  |  factor2 * kapeps^2 * psi d_omega
                    /		  		      
*----------------------------------------------------------------------*/ 
  irow=-1;
   for (inode=0;inode<iel;inode++)
   {
         irow++;
         eforce[irow] += factor2*funct[inode]*pow(kapepsint,2)*facsr;
   } /* end of loop over inode */

if(dynvar->kapeps_flag==1) 
{
/*------------------------------------------------------------------*
   Calculate forces of iteration force vector:

LOW-REYNOLD's MODEL:

                    / 
 (+) (1-THETA)*dt  | 2.0*visc* nue_t* (vderxy_12)^2 * psi   d_omega
                  /  
*------------------------------------------------------------------*/ 
irow = -1;
for (inode=0;inode<iel;inode++)
{
      irow++;
      eforce[irow] += 2.0*visc*eddyint*vderxy_12*funct[inode]*facsr;
} /* end loop over inode */

} /* endif */

 
/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calgaltfv */


/*!--------------------------------------------------------------------- 
\brief stabilisation part of time forces for kapeps dof

<pre>                                                         he  12/02

In this routine the stabilisation part of the time forces for kapeps dofs
is calculated:
		 		                   	  		      
           
                     /
                (+) | tau_tu * u * kapeps * grad(psi)  d_omega + D. C.
                   /  


                     /
   (-) (1-THETA)*dt |  tau_tu * u * grad(kapeps) * u * grad(psi) d_omega + D. C.
                   /          
    
                     /
   (+) (1-THETA)*dt |  tau_tu * div((nue+nue_t/sig)*grad(kapeps)) * u * grad(psi)  d_omega + D. C.
                   /

                     /
   (-) (1-THETA)*dt |  tau_tu * factor * kapeps^2 * grad(psi) * u  d_omega + D. C.
                   /

LOW-REYNOLD's MODEL only for kappa:

                       /
   (-) (1-THETA)*dt   |  tau_tu * grad(k_old)*grad(k_old)/(4*k_old) * 2.0 * visc* grad(psi) * u  d_omega + D. C.
                     /
  


                     /
   (+) (1-THETA)*dt |  tau_tu * 0.5 * nue_t * factor1 * (grad(u)+[grad(u)^T])^2 * grad(psi) * u   d_omega + D. C.
                   /   
      	 
                     /
   (+) (1-THETA)*dt |  tau_tu * factor2 * kapeps^2 * grad(psi) * u  d_omega + D. C.
                   /
		 	   		
LOW-REYNOLD's MODEL only for epsilon:

                     / 
   (+) (1-THETA)*dt | tau_tu * 2.0*visc*nue_t* (vderxy2_12)^2 *  grad(psi) * u  d_omega + D. C.
                   /  
       
</pre>
\param   *dynvar        FLUID_DYN_CALC    (i)
\param   *ele           ELEMENT	      (i)    actual element
\param   *eforce        DOUBLE	      (i/o)  element force vector
\param    kapepsint     DOUBLE	      (i)    kapeps at integr. point
\param   *velint        DOUBLE	      (i)    vel. at integr. point
\param   *velint_dc     DOUBLE	      (i)    vel. at integr. point for D.C.
\param    eddyint       DOUBLE	      (i)    eddy-visc. at integr. point
\param  **derxy         DOUBLE	      (i)    global derivatives
\param   *kapepsderxy2  DOUBLE	      (i)    kapeps 2nd derivatives
\param  **vderxy        DOUBLE	      (i)    global vel. deriv.
\param   *kapepsderxy   DOUBLE	      (i)    global kapeps deriv.
\param    fac	      DOUBLE	      (i)    weighting factor
\param    visc	      DOUBLE	      (i)    fluid viscosity
\param    factor	      DOUBLE	      (i)    factor
\param    factor1	      DOUBLE	      (i)    factor
\param    factor2	      DOUBLE	      (i)    factor
\param    sig	      DOUBLE	      (i)    factor
\param    vderxy_12     DOUBLE	      (i)    factor
\param    production    DOUBLE	      (i)    factor
\param    iel	      INT	            (i)    num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabtfkapeps(
                   FLUID_DYN_CALC  *dynvar,
                   ELEMENT         *ele,      
	            DOUBLE          *eforce,  
	 	      DOUBLE           kapepsint,  
       	      DOUBLE          *velint,  
       	      DOUBLE          *velint_dc,  
		      DOUBLE           eddyint, 
                  DOUBLE         **derxy,   
		      DOUBLE          *kapepsderxy2,   
                  DOUBLE         **vderxy,  
		      DOUBLE          *kapepsderxy,
                  DOUBLE           fac,     
		      DOUBLE           visc,    
                  DOUBLE           factor,
                  DOUBLE           factor1,
                  DOUBLE           factor2,
                  DOUBLE           sig,
                  DOUBLE           vderxy_12,  
                  DOUBLE           production,  
                  INT              iel      
                  )
{
INT    j,irow,isd,inode;
DOUBLE aux,aux_dc;
DOUBLE taumu,taumu_dc;
DOUBLE facsr,facsr_dc;

#ifdef DEBUG 
dstrc_enter("f2_calstabtfkapeps");
#endif

/*--------------------------------------------------- set some factors */
taumu    = dynvar->tau_tu;
taumu_dc = dynvar->tau_tu_dc;
facsr    = fac * dynvar->thsr * taumu;
facsr_dc = fac * dynvar->thsr * taumu_dc;

/*----------------------------------------------------------------------*
   Calculate forces of time force vector:
                     /
                (+) | tau_tu * u * kapeps * grad(psi)     d_omega
                   /  
*----------------------------------------------------------------------*/
   irow=0;
   for (inode=0;inode<iel;inode++)
   {
     aux = kapepsint;
     
      for (isd=0;isd<2;isd++)
      {
	 eforce[irow] += aux*fac*taumu   *derxy[isd][inode]*velint[isd];
	 eforce[irow] += aux*fac*taumu_dc*derxy[isd][inode]*velint_dc[isd]; 
      } /* end loop over isd */
     irow++;
   } /* end loop over inode */

/*----------------------------------------------------------------------*
   Calculate forces of time force vector:
                     /
   (-) (1-THETA)*dt |  tau_tu * u * grad(kapeps) * u * grad(psi) d_omega
                   /          
*----------------------------------------------------------------------*/
   irow=0;
   for (inode=0;inode<iel;inode++)
   {
       aux    =(velint[0]*velint[0]*kapepsderxy[0]+
                velint[0]*velint[1]*kapepsderxy[1])*derxy[0][inode]+
               (velint[1]*velint[0]*kapepsderxy[0]+
                velint[1]*velint[1]*kapepsderxy[1])*derxy[1][inode]; 

       aux_dc =(velint_dc[0]*velint_dc[0]*kapepsderxy[0]+
                velint_dc[0]*velint_dc[1]*kapepsderxy[1])*derxy[0][inode]+
               (velint_dc[1]*velint_dc[0]*kapepsderxy[0]+
                velint_dc[1]*velint_dc[1]*kapepsderxy[1])*derxy[1][inode]; 

	 eforce[irow] -= aux   *facsr;
	 eforce[irow] -= aux_dc*facsr_dc; 
    irow++;
   } /* end loop over inode */
 
/*----------------------------------------------------------------------*
   Calculate forces of time force vector:
                     /
   (+) (1-THETA)*dt |  tau_tu * div((nue+nue_t/sig)*grad(kapeps)) * u * grad(psi)  d_omega=
                   /

                     /
   (+) (1-THETA)*dt |  tau_tu *(nue+nue_t/sig) * div grad(kapeps)] * u * grad(psi)  d_omega 
                   /   
                  
*----------------------------------------------------------------------*/
  irow = 0;
   for (inode=0;inode<iel;inode++)
   {
    aux = (visc+eddyint/sig)*(kapepsderxy2[0]+kapepsderxy2[1]);

    for (isd=0;isd<2;isd++)
    {
        eforce[irow] += aux*facsr   *derxy[isd][inode]*velint[isd];			 
        eforce[irow] += aux*facsr_dc*derxy[isd][inode]*velint_dc[isd];			 
    } /* end loop over irn */
    irow ++;
   } /* end loop over inode*/

/*----------------------------------------------------------------------*
   Calculate forces of time force vector:
                     /
   (-) (1-THETA)*dt |  tau_tu * factor * kapeps^2 * grad(psi) * u  d_omega
                   /
*----------------------------------------------------------------------*/ 
  irow=0;
   for (inode=0;inode<iel;inode++)
   {
    aux = factor*pow(kapepsint,2);
   
      for (isd=0;isd<2;isd++)
      {
         eforce[irow] -= aux*facsr   *velint[isd]   *derxy[isd][inode];
         eforce[irow] -= aux*facsr_dc*velint_dc[isd]*derxy[isd][inode]; 		
      } /* end loop over isd */
   irow++;
   } /* end loop ove inode */

if(dynvar->kapeps_flag==0) 
{
/*----------------------------------------------------------------------*
   Calculate forces of time force vector:

LOW-REYNOLD's MODEL:

                       /
   (-) (1-THETA)*dt   |  tau_tu * grad(k_old)*grad(k_old)/(4*k_old) * 2.0 * visc * grad(psi) * u  d_omega
                     /

*----------------------------------------------------------------------*/ 
  irow=0;
   for (inode=0;inode<iel;inode++)
   {
    aux = 1/(4*kapepsint)*2.0*visc*
          (pow(kapepsderxy[0],2)+pow(kapepsderxy[1],2));
    
      for (isd=0;isd<2;isd++)
      {
         eforce[irow] -= aux*facsr   *velint[isd]   *derxy[isd][inode];
         eforce[irow] -= aux*facsr_dc*velint_dc[isd]*derxy[isd][inode];  
      } /* end loop over isd */
   irow++;
   } /* end loop ove inode */
} /* endif */


/*----------------------------------------------------------------------*
   Calculate forces of time force vector:
                     /
   (+) (1-THETA)*dt |  tau_tu * 0.5 * nue_t * factor1 * (grad(u)+[grad(u)^T])^2 * grad(psi) * u   d_omega
                   /   
 *----------------------------------------------------------------------*/
   irow = 0;
   for (inode=0;inode<iel;inode++)
   {
    aux = eddyint*factor1*production;

      for (isd=0;isd<2;isd++)
      {
        eforce[irow] += aux*facsr   *velint[isd]   *derxy[isd][inode];
        eforce[irow] += aux*facsr_dc*velint_dc[isd]*derxy[isd][inode]; 
      } /* end loop over isd */
  irow ++;
} /* end loop over inode */

/*----------------------------------------------------------------------*
   Calculate forces of time force vector:
                     /
   (+) (1-THETA)*dt |  tau_tu * factor2 * kapeps^2 * grad(psi) * u  d_omega
                   /
*----------------------------------------------------------------------*/ 
  irow=0;
   for (inode=0;inode<iel;inode++)
   {
    aux = factor2*pow(kapepsint,2);
    
      for (isd=0;isd<2;isd++)
      {
         eforce[irow] += aux*facsr   *velint[isd]   *derxy[isd][inode];
         eforce[irow] += aux*facsr_dc*velint_dc[isd]*derxy[isd][inode]; 
      } /* end loop over isd */
   irow++;
   } /* end loop ove inode */


if(dynvar->kapeps_flag==1) 
{
/*---------------------------------------------------------------------- 
   Calculate  stabilastion of iteration force vector:

LOW-REYNOLD's MODEL:

                   / 
 (+) (1-THETA)*dt | tau_tu * 2.0*visc* nue_t* (vderxy2_12)^2 * grad(psi) * u  d_omega
                 /  

*----------------------------------------------------------------------*/ 
 irow=0; 
   for (inode=0;inode<iel;inode++)
   {
     aux = 2.0*visc*eddyint*vderxy_12;
     
      for (isd=0;isd<2;isd++)
      {
         eforce[irow] += aux*facsr   *velint[isd]   *derxy[isd][inode];
         eforce[irow] += aux*facsr_dc*velint_dc[isd]*derxy[isd][inode]; 
      } /* end loop over isd */
     irow ++;		      
   } /* end loop over inode */

} /* endif */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calstabtfkapeps */




#endif

/*!----------------------------------------------------------------------
\file
\brief external RHS for fluid2 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2 
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
/*!--------------------------------------------------------------------- 
\brief galerkin part of external forces for vel dofs

<pre>                                                         genk 09/02

In this routine the galerkin part of the time forces for vel dofs
is calculated:

                   /
   + (1-THETA)*dt |  v * b_old   d_omega
                 /  

               /
   + THETA*dt |  v * b   d_omega
             /      	  		      

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
     
      
</pre>
\param   *dynvar      FLUID_DYN_CALC  (i)
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param	 *funct       DOUBLE	      (i)    nat. shape functions      
\param   *edeadn      DOUBLE          (i)    ele dead load at n
\param   *edeadng     DOUBLE          (i)    ele dead load at n+1
\param	  fac	      DOUBLE	      (i)    weighting factor
\param	  iel	      INT	      (i)    num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calgalexfv(
                  FLUID_DYN_CALC  *dynvar, 
                  DOUBLE          *eforce,     
		  DOUBLE          *funct,       
                  DOUBLE          *edeadn,
		  DOUBLE          *edeadng,
		  DOUBLE           fac,      
		  INT              iel       
              ) 
{
DOUBLE  facsl, facsr;
INT     inode,irow,isd;

#ifdef DEBUG 
dstrc_enter("f2_calgalexfv");
#endif		

/*--------------------------------------------------- set some factors */
facsl = fac*dynvar->thsl;
facsr = fac*dynvar->thsr;

/*----------------------------------------------------------------------*
   Calculate galerkin part of external forces:
   
                   /
   + (1-THETA)*dt |  v * b_old   d_omega
                 /  

               /
   + THETA*dt |  v * b   d_omega
             /
   
 *----------------------------------------------------------------------*/
irow=-1;
for (inode=0;inode<iel;inode++)
{
   for (isd=0;isd<2;isd++)
   {
      irow++;
      eforce[irow] += funct[inode]*(edeadn[isd]*facsr+edeadng[isd]*facsl);
   } /* end of loop over isd */
} /* end of loop over inode */ 
 
 
/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calgalexfv */ 

/*!--------------------------------------------------------------------- 
\brief stabilisation part of external forces for vel dofs

<pre>                                                         genk 09/02

In this routine the stabilisation part of the time forces for vel dofs
is calculated:

EULER:
                     /
   + thetas(l,r)*dt |  taum_mu * u * grad(v) * b^   d_omega
                   /  

ALE:
                     /
   + thetas(l,r)*dt |  taum_mu * c * grad(v) * b^   d_omega
                   /

EULER/ALE:
                       /
   -/+ thetas(l,r)*dt |  tau_mp * 2*nue * div( eps(v) ) * b^  d_omega
                     /      	  		      

This routine is called twice with different values:
1. values are added to Iteration RHS (evaluation at n+1):
    thetas(l,r) = THETA*dt
    b^ = b = deadng
2. values are added to Time RHS (evaluation at n):
   thetas(l,r) = (1-THETA)*dt
   b^ = b_old = deadn    

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE: for EULER
      velint = U

NOTE: for ALE
      velint = C (ale-convective velocity)           
      
</pre>
\param   *dynvar      FLUID_DYN_CALC  (i)
\param   *ele         ELEMENT	      (i)    actual element
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param  **derxy,      DOUBLE	      (i)    global derivatives
\param  **derxy2,     DOUBLE	      (i)    2nd. global derivatives    
\param   *edead       DOUBLE          (i)    ele dead load at n or n+1
\param   *velint      DOUBLE	      (i)    vel. at integr. point
\param	  fac	      DOUBLE	      (i)    weighting factor
\param	  visc	      DOUBLE	      (i)    fluid viscosity
\param	  iel	      INT	      (i)    num. of nodes in ele
\param	  ihoel       INT	      (i)    flag for higer ord. ele
\param	  flag	      INT	      (i)    flag for n or n+1
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabexfv(
                    FLUID_DYN_CALC  *dynvar, 
                    ELEMENT         *ele,  
                    DOUBLE          *eforce,     
		    DOUBLE         **derxy,
		    DOUBLE         **derxy2,      
                    DOUBLE          *edead,
	 	    DOUBLE          *velint,  
		    DOUBLE           fac,      
                    DOUBLE           visc,
		    INT              iel,
		    INT              ihoel,
		    INT              flag      
                   ) 
{

INT    irow,inode,isd;
DOUBLE sign,aux;
DOUBLE taumu,taump;
DOUBLE c,fvts; 
DOUBLE fact[2];
STAB_PAR_GLS *gls;	/* pointer to GLS stabilisation parameters	*/

#ifdef DEBUG 
dstrc_enter("f2_calstabexfv");
#endif		

/*---------------------------------------------------------- initialise */
gls    = ele->e.f2->stabi.gls;

if (ele->e.f2->stab_type != stab_gls) 
   dserror("routine with no or wrong stabilisation called");

/*--------------------------------------------------- set some factors */
taumu = dynvar->tau[0];
taump = dynvar->tau[1];
switch (flag)
{
case 0: /* evaluation at n */
   c = dynvar->thsr;
break;
case 1: /* evaluation at n+1 */
   c = dynvar->thsl;
break;
default:
   dserror("value of flag not valid!!!\n");
}

/*----------------------------------------------------------------------*
   Calculate external/convective stab-forces of time force vector:
EULER:
                     /
   + thetas(l,r)*dt |  taum_mu * u * grad(v) * b^   d_omega
                   /  
ALE:
                     /
   + thetas(l,r)*dt |  taum_mu * c * grad(v) * b^   d_omega
                   /
 *----------------------------------------------------------------------*/
if (gls->iadvec!=0)
{
   fact[0] = edead[0]*fac*taumu*c;
   fact[1] = edead[1]*fac*taumu*c;
   irow=-1;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1];
      for (isd=0;isd<2;isd++)
      {
         irow++;
	 eforce[irow] += aux*fact[isd];
      } /* end loop over isd */
   } /* end loop over inode */
} /* endif (ele->e.f2->iadvec!=0) */

/*----------------------------------------------------------------------*
   Calculate external/viscous stab-forces of time force vector:
                       /
   -/+ thetas(l,r)*dt |  tau_mp * 2*nue * div( eps(v) ) * b  d_omega
                     /
 *----------------------------------------------------------------------*/
if (gls->ivisc!=0 && ihoel!=0)
{
   switch (gls->ivisc) /* choose stabilisation type --> sign */
   {
   case 1: /* GLS- */
      sign = ONE;
   break;
   case 2: /* GLS+ */
      sign = -ONE;
   break;
   default:
      dserror("viscous stabilisation parameter unknown: IVISC");
   } /* end switch (ele->e.f2->ivisc) */

   fvts = fac*visc*taump*sign*c;
   irow = 0;
   for (inode=0;inode<iel;inode++)
   {
      eforce[irow]   -= ((TWO*derxy2[0][inode] + derxy2[1][inode])*edead[0] \
                            + derxy2[2][inode]*edead[1])*fvts;
      eforce[irow+1] -= ((TWO*derxy2[1][inode] + derxy2[0][inode])*edead[1] \
                            + derxy2[2][inode]*edead[0])*fvts;
      irow += 2;			 
   } /* end loop over inode */
} /* endif (ele->e.f2->ivisc!=0 && ihoel!=0) */

/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_stabexfv */ 

/*!--------------------------------------------------------------------- 
\brief stabilisation part of external forces for pre dofs

<pre>                                                         genk 09/02

In this routine the stabilisation part of the time forces for pre dofs
is calculated:

                      /
   - thetas(l,r)*dt  |  tau_mp * grad(q) * b^  d_omega
                    /	  		      

This routine is called twice with different values:
1. values are added to Iteration RHS (evaluation at n+1):
    thetas(l,r) = THETA*dt
    b^ = b = deadng
2. values are added to Time RHS (evaluation at n):
   thetas(l,r) = (1-THETA)*dt
   b^ = b_old = deadn    

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE:							
    there's only one full element force vector  	
    for pre-dofs the pointer eforce points to the entry 
    eforce[2*iel]     
      
</pre>
\param   *dynvar      FLUID_DYN_CALC  (i)
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param  **derxy,      DOUBLE	      (i)    global derivatives    
\param   *edead       DOUBLE          (i)    ele dead load at n or n+1
\param   *velint      DOUBLE	      (i)    vel. at integr. point
\param	  fac	      DOUBLE	      (i)    weighting factor
\param	  iel	      INT	      (i)    num. of nodes in ele
\param	  flag	      INT	      (i)    flag for n or n+1
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabexfp(
                    FLUID_DYN_CALC  *dynvar, 
                    DOUBLE          *eforce,     
		    DOUBLE         **derxy,       
                    DOUBLE          *edead,  
		    DOUBLE           fac,      
		    INT              iel,
		    INT              flag      
                   ) 
{
INT    inode;
DOUBLE c;
DOUBLE taump;
DOUBLE fact[2];

/*--------------------------------------------------- set some factors */
taump = dynvar->tau[1];
switch (flag)
{
case 0: /* evaluation at n */
   c = dynvar->thpr;
break;
case 1: /* evaluation at n+1 */
   c = dynvar->thpl;
break;
default:
   dserror("value of flag not valid!!!\n");
}

/*----------------------------------------------------------------------*
   Calculate inertia/pressure stab forces of time force vector:
                      /
   - thetas(l,r)*dt  |  tau_mp * grad(q) * b  d_omega
                    /	  		      
 *----------------------------------------------------------------------*/
fact[0] = edead[0]*taump*fac*c;
fact[1] = edead[1]*taump*fac*c;
for (inode=0;inode<iel;inode++)
{
   eforce[inode] -= (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1]);
}

/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_stabexfp */ 

#endif
/*! @} (documentation module close)*/

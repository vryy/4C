/*!----------------------------------------------------------------------
\file
\brief external RHS for fluid2 element

------------------------------------------------------------------------*/
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
\param   *eforce      double	      (i/o)  element force vector
\param	 *funct       double	      (i)    nat. shape functions      
\param   *edeadn      double          (i)    ele dead load at n
\param   *edeadng     double          (i)    ele dead load at n+1
\param	  fac	      double	      (i)    weighting factor
\param	  iel	      int	      (i)    num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calgalexfv(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,     
		  double          *funct,       
                  double          *edeadn,
		  double          *edeadng,
		  double           fac,      
		  int              iel       
              ) 
{
double  facsl, facsr;
int     inode,irow,isd;

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

                     /
   + thetas(l,r)*dt |  taum_mu * u * grad(v) * b^   d_omega
                   /  

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
     
      
</pre>
\param   *dynvar      FLUID_DYN_CALC  (i)
\param   *ele         ELEMENT	      (i)    actual element
\param   *eforce      double	      (i/o)  element force vector
\param  **derxy,      double	      (i)    global derivatives
\param  **derxy2,     double	      (i)    2nd. global derivatives    
\param   *edead       double          (i)    ele dead load at n or n+1
\param   *velint      double	      (i)    vel. at integr. point
\param	  fac	      double	      (i)    weighting factor
\param	  visc	      double	      (i)    fluid viscosity
\param	  iel	      int	      (i)    num. of nodes in ele
\param	  ihoel       int	      (i)    flag for higer ord. ele
\param	  flag	      int	      (i)    flag for n or n+1
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabexfv(
                    FLUID_DYN_CALC  *dynvar, 
                    ELEMENT         *ele,  
                    double          *eforce,     
		    double         **derxy,
		    double         **derxy2,      
                    double          *edead,
	 	    double          *velint,  
		    double           fac,      
                    double           visc,
		    int              iel,
		    int              ihoel,
		    int              flag      
                   ) 
{

int    irow,inode,isd;
double sign,aux;
double taumu,taump;
double c,fvts; 
double fact[2];
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
                     /
   + thetas(l,r)*dt |  tau_mu * u * grad(v) * b  d_omega
                   /
 *----------------------------------------------------------------------*/
if (ele->e.f2->iadvec!=0)
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
if (ele->e.f2->ivisc!=0 && ihoel!=0)
{
   switch (ele->e.f2->ivisc) /* choose stabilisation type --> sign */
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
\param   *eforce      double	      (i/o)  element force vector
\param  **derxy,      double	      (i)    global derivatives    
\param   *edead       double          (i)    ele dead load at n or n+1
\param   *velint      double	      (i)    vel. at integr. point
\param	  fac	      double	      (i)    weighting factor
\param	  iel	      int	      (i)    num. of nodes in ele
\param	  flag	      int	      (i)    flag for n or n+1
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabexfp(
                    FLUID_DYN_CALC  *dynvar, 
                    double          *eforce,     
		    double         **derxy,       
                    double          *edead,  
		    double           fac,      
		    int              iel,
		    int              flag      
                   ) 
{
int    inode;
double c;
double taump;
double fact[2];

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

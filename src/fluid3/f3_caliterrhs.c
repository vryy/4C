/*!----------------------------------------------------------------------
\file
\brief iteration RHS for fluid2 element

------------------------------------------------------------------------*/
#ifdef D_FLUID3 
#include "../headers/standardtypes.h"
#include "fluid3_prototypes.h"
#include "fluid3.h"

/*!--------------------------------------------------------------------- 
\brief galerkin part of iteration forces for vel dofs

<pre>                                                         genk 05/02

In this routine the galerkin part of the iteration forces for vel dofs
is calculated:

                   /
   (+/-) THETA*dt |  v * u * grad(u)  d_omega
                 /  


see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
      
</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *eforce    double	   (i/o)  element force vector
\param  *covint    double	   (i)	  conv. vels at int. point
\param  *funct     double	   (i)    nat. shape funcs
\param   fac 	   double	   (i)    weighting factor
\param	 iel	   int		   (i)	  num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calgalifv(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,    
		  double          *covint,
		  double          *funct,
		  double           fac,
		  int              iel
                 )  
{
int    inode,isd;
int    irow;  
double facsl;
 
#ifdef DEBUG 
dstrc_enter("f3_calgalifv");
#endif

/*--------------------------------------------------- set some factors */
facsl = fac * dynvar->thsl * dynvar->sigma;

/*----------------------------------------------------------------------*
   Calculate convective forces of iteration force vector:
                   /
   (+/-) THETA*dt |  v * u * grad(u)  d_omega
    |            /  
    |-> signs due to nonlin. iteration scheme (dynvar->sigma) 
 *----------------------------------------------------------------------*/ 
irow = -1;
for (inode=0;inode<iel;inode++)
{
   for (isd=0;isd<3;isd++)
   {
      irow++;
      eforce[irow] += funct[inode]*covint[isd]*facsl;
   } /* end of loop over isd */
} /* end of loop over inode */

/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calgalifv */

/*!---------------------------------------------------------------------  
\brief stabilisation part of iteration forces for vel dofs

<pre>                                                         genk 05/02

In this routine the stabilisation part of the iteration forces for vel dofs
is calculated:

                   /
   (+/-) THETA*dt |  tau_mu * u * grad(v) * u * grad(u)  d_omega
                 /    

                       /
   (-/+) -/+ THETA*dt |  tau_mp * 2*nue * div( eps(v) ) *  * grad(u)  d_omega
                     /



see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
      
</pre>
\param   *dynvar   FLUID_DYN_CALC  (i)
\param   *ele      ELEMENT	   (i)    actual element
\param   *eforce   double	   (i/o)  element force vector
\param   *covint   double	   (i)    conv. vels at int. point
\param   *velint   double	   (i)    vel at integr. point
\param	 *funct    double	   (i)	  nat. shape funcs
\param  **derxy    double	   (i)	  global derivative
\param  **derxy2   double	   (i)    2nd global derivative
\param    fac 	   double	   (i)    weighting factor
\param    visc     double	   (i)    fluid viscosity
\param    ihoel    int		   (i)    flag for higher order ele
\param	  iel	   int		   (i)    num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabifv(
                  FLUID_DYN_CALC  *dynvar, 
                  ELEMENT         *ele,
		  double          *eforce,   /* element force vector */
		  double          *covint,
		  double          *velint,
		  double          *funct,
		  double         **derxy,
		  double         **derxy2,
		  double           fac,
		  double           visc,
		  int              ihoel,
		  int              iel
                 )  
{
int    inode,isd;
int    irow;  
double facsl,cc;
double aux, sign;
double taumu,taump;
double fact[3];
 
#ifdef DEBUG 
dstrc_enter("f3_calgalifv");
#endif

/*--------------------------------------------------- set some factors */
facsl = fac * dynvar->thsl * dynvar->sigma;
taumu = dynvar->tau[0];
taump = dynvar->tau[1];

/*----------------------------------------------------------------------*
   Calculate convective/convective stabilastion of iteration force vector:
                   /
   (+/-) THETA*dt |  tau_mu * u * grad(v) * u * grad(u)  d_omega
    |            /  
    |           
    |-> signs due to nonlin. iteration scheme (dynvar->sigma)
 *----------------------------------------------------------------------*/ 
if (ele->e.f3->iadvec!=0)
{
   fact[0] = taumu*covint[0]*facsl;
   fact[1] = taumu*covint[1]*facsl;
   fact[2] = taumu*covint[2]*facsl;
   irow=-1;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1] \
          + derxy[2][inode]*velint[2];
      for (isd=0;isd<3;isd++)
      {
         irow++;
	 eforce[irow] += aux*fact[isd];
      } /* end of loop over isd */
   } /* end of loop over inode */
} /* endif (ele->e.f3->iadvec!=0) */

/*----------------------------------------------------------------------*
   Calculate convective/viscous stabilastion of iteration force vector:
                       /
   (-/+) -/+ THETA*dt |  tau_mp * 2*nue * div( eps(v) ) *  * grad(u)  d_omega
    |                /  
    |           
    |-> signs due to nonlin. iteration scheme (dynvar->sigma)
 *----------------------------------------------------------------------*/ 
if (ele->e.f3->ivisc!=0 && ihoel!=0)
{
   switch (ele->e.f3->ivisc) /* choose stabilisation type --> sign */
   {
   case 1: /* GLS- */
      sign = ONE;
   break;
   case 2: /* GLS+ */
      sign = -ONE;
   break;
   default:
      dserror("viscous stabilisation parameter unknown: IVISC");
   } /* end switch (ele->e.f3->ivisc) */
   
   cc = facsl*visc*taump*sign;
   irow=0;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy2[0][inode] + derxy2[1][inode] + derxy2[2][inode];
      eforce[irow]   -= ((derxy2[0][inode]+aux)*covint[0]    \
                        + derxy2[3][inode]     *covint[1]    \
                        + derxy2[4][inode]     *covint[2])*cc;
      eforce[irow+1] -= ((derxy2[1][inode]+aux)*covint[1]    \
                        + derxy2[3][inode]     *covint[0]    \
                        + derxy2[5][inode]     *covint[2])*cc;
      eforce[irow+2] -= ((derxy2[2][inode]+aux)*covint[2]    \
                        + derxy2[4][inode]     *covint[0]    \
                        + derxy2[5][inode]     *covint[1])*cc;
      irow += 3;		      
   } /* end of loop over inode */
} /* endif (ele->e.f3->ivisc!=0 && ihoel!=0) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calgalifv */

/*!--------------------------------------------------------------------- 
\brief stabilisation part of iteration forces for pre dofs

<pre>                                                         genk 05/02

In this routine the stabilisation part of the iteration forces for pre 
dofs is calculated:

                   /
   (-/+) THETA*dt |  tau_mp * grad(q) * u * grad(u)  d_omega
                 /  

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE:							
    there's only one full element force vector  	
    for pre-dofs the pointer eforce points to the entry  
    eforce[3*iel]					
      
</pre>
\param   *dynvar   FLUID_DYN_CALC  (i)
\param   *eforce   double	   (i/o)  element force vector
\param   *covint   double	   (i)    conv. vels at int. point
\param  **derxy    double	   (i)    global derivative
\param    fac      double	   (i)    weighting factor
\param	  iel	   int  	   (i)	  num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabifp(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,    
		  double          *covint,
		  double          **derxy,
		  double           fac,
		  int              iel
                 )  
{
int    inode; 
double facsl;
double taump;
double fact[3];
 
#ifdef DEBUG 
dstrc_enter("f3_calgstabifp");
#endif

/*--------------------------------------------------- set some factors */
facsl = fac * dynvar->thsl * dynvar->sigma;
taump = dynvar->tau[1];

/*----------------------------------------------------------------------*
   Calculate convective pressure stabilisation iteration force vector:
                   /
   (-/+) THETA*dt |  tau_mp * grad(q) * u * grad(u)  d_omega
    |            /  
    |-> signs due to nonlin. iteration scheme (dynvar->sigma) 
 *----------------------------------------------------------------------*/ 
fact[0] = covint[0]*taump*facsl;
fact[1] = covint[1]*taump*facsl;
fact[2] = covint[2]*taump*facsl;
for (inode=0;inode<iel;inode++)
{
   eforce[inode] -= (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1] \
                   + derxy[2][inode]*fact[2] );
} /* end of loop over inode */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calstabifp */

#endif

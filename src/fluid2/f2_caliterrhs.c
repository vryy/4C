/*!----------------------------------------------------------------------
\file
\brief iteration RHS for fluid2 element

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2 
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
/*!--------------------------------------------------------------------- 
\brief galerkin part of iteration forces for vel dofs

<pre>                                                         genk 04/02

In this routine the galerkin part of the iteration forces for vel dofs
is calculated:

EULER:

                   /
   (+/-) THETA*dt |  v * u * grad(u)  d_omega
                 /  

ALE:
                   /
   (+/-) THETA*dt |  v * c * grad(u)  d_omega
                 /
NOTE:
  EULER:  covint = u*grad(u)
  ALE:    covint = c*grad(u)

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
void f2_calgalifv(
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
dstrc_enter("f2_calgalifv");
#endif

/*--------------------------------------------------- set some factors */
facsl = fac * dynvar->thsl * dynvar->sigma;

/*----------------------------------------------------------------------*
   Calculate convective forces of iteration force vector:

EULER:
                   /
   (+/-) THETA*dt |  v * u * grad(u)  d_omega
    |            /  
    |-> signs due to nonlin. iteration scheme (dynvar->sigma) 

ALE:
                   /
   (+/-) THETA*dt |  v * c * grad(u)  d_omega
                 /
 *----------------------------------------------------------------------*/ 
irow = -1;
for (inode=0;inode<iel;inode++)
{
   for (isd=0;isd<2;isd++)
   {
      irow++;
      eforce[irow] += funct[inode]*covint[isd]*facsl;
   } /* end loop over isd */
} /* end loop over inode */

/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calgalifv */

/*!---------------------------------------------------------------------  
\brief stabilisation part of iteration forces for vel dofs

<pre>                                                         genk 04/02
                                            modified for ALE  genk 01/03

In this routine the stabilisation part of the iteration forces for vel dofs
is calculated:

EULER:
                   /
   (+/-) THETA*dt |  tau_mu * u * grad(v) * u * grad(u)  d_omega
                 /    

                       /
   (-/+) -/+ THETA*dt |  tau_mp * 2*nue * div( eps(v) ) * u * grad(u)  d_omega
                     /

ALE:
                   /
   (+/-) THETA*dt |  tau_mu * c * grad(v) * c * grad(u)  d_omega
                 /    

                       /
   (-/+) -/+ THETA*dt |  tau_mp * 2*nue * div( eps(v) ) * c * grad(u)  d_omega
                     /

NOTE:
  EULER:  covint = u*grad(u)   velint = u
  ALE:    covint = c*grad(u)   velint = alecovint (c)

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
void f2_calstabifv(
                  FLUID_DYN_CALC  *dynvar, 
                  ELEMENT         *ele,      
		  double          *eforce,  
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
double fact[2];
 
#ifdef DEBUG 
dstrc_enter("f2_calgalifv");
#endif

/*--------------------------------------------------- set some factors */
facsl = fac * dynvar->thsl * dynvar->sigma;
taumu = dynvar->tau[0];
taump = dynvar->tau[1];

/*----------------------------------------------------------------------*
   Calculate convective/convective stabilastion of iteration force vector:
EULER:
                   /
   (+/-) THETA*dt |  tau_mu * u * grad(v) * u * grad(u)  d_omega
    |            /  
    |           
    |-> signs due to nonlin. iteration scheme (dynvar->sigma)

ALE:
                   /
   (+/-) THETA*dt |  tau_mu * c * grad(v) * c * grad(u)  d_omega
                 / 
 *----------------------------------------------------------------------*/ 
if (ele->e.f2->iadvec!=0)
{
   fact[0] = taumu*covint[0]*facsl;
   fact[1] = taumu*covint[1]*facsl;
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
   Calculate convective/viscous stabilastion of iteration force vector:
EULER:
                       /
   (-/+) -/+ THETA*dt |  tau_mp * 2*nue * div( eps(v) ) * u * grad(u)  d_omega
    |                /  
    |           
    |-> signs due to nonlin. iteration scheme (dynvar->sigma)

ALE:
                       /
   (-/+) -/+ THETA*dt |  tau_mp * 2*nue * div( eps(v) ) * c * grad(u)  d_omega
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
   } /* end switch(ele->e.f2->ivisc) */
   cc = facsl*visc*taump*sign;
   irow=0;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy2[0][inode] + derxy2[1][inode];
      eforce[irow]   -= ((derxy2[0][inode]+aux)*covint[0]    \
                        + derxy2[2][inode]     *covint[1])*cc;
      eforce[irow+1] -= ((derxy2[1][inode]+aux)*covint[1]    \
                        + derxy2[2][inode]     *covint[0])*cc;
      irow += 2;		      
   } /* end loop over inode */
} /* endif (ele->e.f2->ivisc!=0 && ihoel!=0) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calgalifv */

/*!--------------------------------------------------------------------- 
\brief stabilisation part of iteration forces for pre dofs

<pre>                                                         genk 04/02

In this routine the stabilisation part of the iteration forces for pre 
dofs is calculated:

EULER:
                   /
   (-/+) THETA*dt |  tau_mp * grad(q) * u * grad(u)  d_omega
                 /  

ALE:
                   /
   (-/+) THETA*dt |  tau_mp * grad(q) * c * grad(u)  d_omega
                 / 
		 
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE:
  EULER:  covint = u*grad(u)  
  ALE:    covint = c*grad(u)  

NOTE:							
    there's only one full element force vector  	
    for pre-dofs the pointer eforce points to the entry 
    eforce[2*iel]					
      
</pre>
\param   *dynvar   FLUID_DYN_CALC  (i)
\param   *eforce   double	   (i/o)  element force vector
\param   *covint   double	   (i)    conv. vels at int. point
\param  **derxy    double	   (i)    global derivative
\param    fac      double	   (i)    weighting factor
\param	  iel	   int  	   (i)	  num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabifp(
                   FLUID_DYN_CALC  *dynvar, 
                   double          *eforce,   
		   double          *covint,  
		   double         **derxy,   
		   double           fac,     
		   int              iel        
                 )  
{
int    inode; 
double facsl;
double taump;
double fact[2];
 
#ifdef DEBUG 
dstrc_enter("f2_calgstabifp");
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
for (inode=0;inode<iel;inode++)
{
   eforce[inode] -= (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1]);
} /* end loop over inode */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calstabifp */



#endif
/*! @} (documentation module close)*/

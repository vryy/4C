#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#define ONE (1.0)
#define TWO (2.0)
/*----------------------------------------------------------------------*
 | routine to calculate galerkin part of iteration forces for vel dofs  |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_calgalifv(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,   /* element force vector */
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
                   /
   (+/-) THETA*dt |  v * u * grad(u)  d_omega
    |            /  
    |-> signs due to nonlin. iteration scheme (dynvar->sigma) 
 *----------------------------------------------------------------------*/ 
irow = -1;
for (inode=0;inode<iel;inode++)
{
   for (isd=0;isd<2;isd++)
   {
      irow++;
      eforce[irow] += funct[inode]*covint[isd]*facsl;
   }
}


#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_calgalifv */

/*----------------------------------------------------------------------*
 | routine to calculate stabilisation part of iteration forces          |
 | for vel dofs                                             genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_calstabifv(
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
                   /
   (+/-) THETA*dt |  tau_mu * u * grad(v) * u * grad(u)  d_omega
    |            /  
    |           
    |-> signs due to nonlin. iteration scheme (dynvar->sigma)
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
      } 
   }
}

/*----------------------------------------------------------------------*
   Calculate convective/viscous stabilastion of iteration force vector:
                       /
   (-/+) -/+ THETA*dt |  tau_mp * 2*nue * div( eps(v) ) *  * grad(u)  d_omega
    |                /  
    |           
    |-> signs due to nonlin. iteration scheme (dynvar->sigma)
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
   }
   cc = facsl*visc*taump*sign;
   irow=0;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy2[0][inode] + derxy2[1][inode];
      eforce[irow]   -= ((derxy2[0][inode]+aux)*covint[0] \
                        + derxy2[2][inode]*covint[1])*cc;
      eforce[irow+1] -= ((derxy2[1][inode]+aux)*covint[1] \
                        + derxy2[2][inode]*covint[0])*cc;
      irow += 2;		      
   }   
}   


#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_calgalifv */

/*----------------------------------------------------------------------*
 | routine to calculate stabilisation part of iteration forces          |
 | for pre dofs                                                         |
 | NOTE:                                                                |
 |     there's only one full element force vector                       |  
 |     for pre-dofs the pointer eforce points to the entry              |
 |     eforce[2*iel]                                                    |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_calstabifp(
                  FLUID_DYN_CALC  *dynvar, 
                  double          *eforce,   /* element force vector */
		  double          *covint,
		  double          **derxy,
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
}


#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_calstabifp */



#endif

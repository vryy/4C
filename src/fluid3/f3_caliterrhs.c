/*!----------------------------------------------------------------------
\file
\brief iteration RHS for fluid2 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID3 
#include "../headers/standardtypes.h"
#include "fluid3_prototypes.h"
#include "fluid3.h"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;   
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

static FLUID_DYNAMIC   *fdyn;
/*!--------------------------------------------------------------------- 
\brief galerkin part of iteration forces for vel dofs

<pre>                                                         genk 05/02

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
\param  *eforce    DOUBLE	   (i/o)  element force vector
\param  *covint    DOUBLE	   (i)	  conv. vels at INT. point
\param  *funct     DOUBLE	   (i)    nat. shape funcs
\param   fac 	   DOUBLE	   (i)    weighting factor
\param	 iel	   INT		   (i)	  num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calgalifv(
                  DOUBLE          *eforce,    
                  DOUBLE          *covint,
                  DOUBLE          *funct,
                  DOUBLE           fac,
                  INT              iel
                 )  
{
INT    inode,isd;
INT    irow;  
DOUBLE facsl;
 
#ifdef DEBUG 
dstrc_enter("f3_calgalifv");
#endif

/*--------------------------------------------------- set some factors */   
fdyn  = alldyn[genprob.numff].fdyn;
facsl = fac * fdyn->thsl * fdyn->sigma;

/*----------------------------------------------------------------------*
   Calculate convective forces of iteration force vector:

EULER:
               /
   + THETA*dt |  v * u * grad(u)  d_omega
             /  

ALE:
               /
   + THETA*dt |  v * c * grad(u)  d_omega
             /
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
\param   *gls      STAB_PAR_GLS    (i)    stabilisation
\param   *ele      ELEMENT	   (i)    actual element
\param   *eforce   DOUBLE	   (i/o)  element force vector
\param   *covint   DOUBLE	   (i)    conv. vels at INT. point
\param   *velint   DOUBLE	   (i)    vel at integr. point
\param	 *funct    DOUBLE	   (i)	  nat. shape funcs
\param  **derxy    DOUBLE	   (i)	  global derivative
\param  **derxy2   DOUBLE	   (i)    2nd global derivative
\param    fac 	   DOUBLE	   (i)    weighting factor
\param    visc     DOUBLE	   (i)    fluid viscosity
\param    ihoel    INT		   (i)    flag for higher order ele
\param	  iel	   INT		   (i)    num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabifv(
                     STAB_PAR_GLS    *gls,  
                     ELEMENT         *ele,
                     DOUBLE          *eforce,   
                     DOUBLE          *covint,
                     DOUBLE          *velint,
                     DOUBLE          *funct,
                     DOUBLE         **derxy,
                     DOUBLE         **derxy2,
                     DOUBLE           fac,
                     DOUBLE           visc,
                     INT              ihoel,
                     INT              iel
                 )  
{
INT    inode,isd;
INT    irow;  
DOUBLE facsl,cc;
DOUBLE aux;
DOUBLE sign=ONE;
DOUBLE taumu,taump;
DOUBLE fact[3];

#ifdef DEBUG 
dstrc_enter("f3_calgalifv");
#endif	

/*---------------------------------------------------------- initialise */
fdyn  = alldyn[genprob.numff].fdyn;
gls   = ele->e.f3->stabi.gls;

dsassert(ele->e.f3->stab_type == stab_gls, 
        "routine with no or wrong stabilisation called");
 
/*--------------------------------------------------- set some factors */
facsl = fac * fdyn->thsl * fdyn->sigma;
taumu = fdyn->tau[0];
taump = fdyn->tau[1];

/*----------------------------------------------------------------------*
   Calculate convective/convective stabilastion of iteration force vector:
EULER:
               /
   + THETA*dt |  tau_mu * u * grad(v) * u * grad(u)  d_omega
             /  

ALE:
               /
   + THETA*dt |  tau_mu * c * grad(v) * c * grad(u)  d_omega
             / 
 *----------------------------------------------------------------------*/ 
if (gls->iadvec!=0)
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
EULER:
                  /
    -/+ THETA*dt |  tau_mp * 2*nue * div( eps(v) ) * u * grad(u)  d_omega
                /  

ALE:
                 /
   -/+ THETA*dt |  tau_mp * 2*nue * div( eps(v) ) * c * grad(u)  d_omega
               /

 *----------------------------------------------------------------------*/

if (gls->ivisc!=0 && ihoel!=0)
{
   if (gls->ivisc==2) sign*=-ONE; /* GLS+ stabilisation */
   
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

EULER:
               /
   - THETA*dt |  tau_mp * grad(q) * u * grad(u)  d_omega
             /  

ALE:
               /
   - THETA*dt |  tau_mp * grad(q) * c * grad(u)  d_omega
             / 
		 
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE:
  EULER:  covint = u*grad(u)  
  ALE:    covint = c*grad(u)  

NOTE:							
    there's only one full element force vector  	
    for pre-dofs the pointer eforce points to the entry  
    eforce[3*iel]					
      
</pre>
\param   *gls      STAB_PAR_GLS    (i)    stabilisation
\param   *eforce   DOUBLE	   (i/o)  element force vector
\param   *covint   DOUBLE	   (i)    conv. vels at INT. point
\param  **derxy    DOUBLE	   (i)    global derivative
\param    fac      DOUBLE	   (i)    weighting factor
\param	  iel	   INT  	   (i)	  num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabifp(
                     STAB_PAR_GLS    *gls,  
                     DOUBLE          *eforce,    
                     DOUBLE          *covint,
                     DOUBLE          **derxy,
                     DOUBLE           fac,
                     INT              iel
                 ) 
{
INT    inode; 
DOUBLE facsl;
DOUBLE taump;
DOUBLE fact[3];
 
#ifdef DEBUG 
dstrc_enter("f3_calgstabifp");
#endif

if (gls->ipres==0) goto end; /* no pressure stabilisation */

/*--------------------------------------------------- set some factors */
fdyn  = alldyn[genprob.numff].fdyn;
facsl = fac * fdyn->thsl * fdyn->sigma;
taump = fdyn->tau[1];

/*----------------------------------------------------------------------*
   Calculate convective pressure stabilisation iteration force vector:
EULER:
               /
   - THETA*dt |  tau_mp * grad(q) * u * grad(u)  d_omega
             /  
                   
ALE:
               /
   - THETA*dt |  tau_mp * grad(q) * c * grad(u)  d_omega
             / 
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
end:
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calstabifp */

#endif

/*!----------------------------------------------------------------------
\file
\brief time RHS for fluid2_pro element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2_PRO 
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2_PRO 
#include "../headers/standardtypes.h"
#include "fluid2pro_prototypes.h"
#include "fluid2pro.h"
/*!--------------------------------------------------------------------- 
\brief galerkin part of time forces for vel dofs

<pre>                                                        basol 11/02

In this routine the galerkin part of the time forces for vel dofs
is calculated:

      /
   + |  v * u     d_omega   = MUn
    /  
    
         /
   (-)  |  v * u * grad(u)  d_omega  = - N(Un)Un
       /
		  

        /
    +  |  div(v) * p  d_omega  =  CPn
      /		  		      

etforce = MUn+dt*(fn-N(Un)Un) 
eiforce = -dt*CPn     
</pre>
\param   *dynvar      FLUID_DYN_CALC  (i)
\param   *etforce     DOUBLE	      (i/o)  element force vector part1
\param   *eiforce     DOUBLE	      (i/o)  element force vector part2
\param   *velint      DOUBLE	      (i)    vel. at integr. point
\param   *covint      DOUBLE	      (i)    conv. vel. at integr. p.
\param   **vderxy     DOUBLE	      (i)    velocity derivatives 
\param	 *funct       DOUBLE	      (i)    nat. shape functions      
\param	**derxy       DOUBLE	      (i)    global derivatives
\param	  preint      DOUBLE	      (i)    pres. at integr. point
\param	  visc	      DOUBLE	      (i)    fluid viscosity
\param	  fac	      DOUBLE	      (i)    weighting factor
\param	  dt	      DOUBLE	      (i)    incremental time step
\param	  iel	      INT	      (i)    num. of nodes in ele
\return void                                                                       
------------------------------------------------------------------------*/
void f2pro_calgaltfv(
                  FLUID_DYN_CALC  *dynvar, 
                  DOUBLE          *etforce,
		  DOUBLE          *eiforce,
		  DOUBLE          *velint,   
		  DOUBLE          *covint,
		  DOUBLE         **vderxy,   
		  DOUBLE          *funct,    
		  DOUBLE         **derxy,		      
		  DOUBLE           preint,   
		  DOUBLE           visc,     
		  DOUBLE           fac,
		  DOUBLE           dt,      
		  INT               iel       
              )  
{
INT    j,irow,isd,inode;  
DOUBLE aux,c;
DOUBLE fact[2];
DOUBLE dta=dynvar->dta;

#ifdef DEBUG 
dstrc_enter("f2pro_calgaltfv");
#endif

/*----------------------------------------------------------------------*
Calculate intertia forces of time force vector:
           /
MUn =   + |  v * u     d_omega
         /  
 *----------------------------------------------------------------------*/ 
fact[0] = velint[0]*fac;
fact[1] = velint[1]*fac;
c = fac*visc;
irow=-1;
for (inode=0;inode<iel;inode++)
{
   for (isd=0;isd<2;isd++)
   {
      irow++;
      eiforce[irow] += funct[inode]*fact[isd];
   } /* end of loop over isd */
} /* end of loop over inode */

/*----------------------------------------------------------------------*
   Calculate viscous forces of time force vector:
                    /
   - (1-THETA)*dt  |  2*nue * eps(v) : eps(u)  d_omega
                  /  
   It is needed when the viscous term is integrated via trapezoidal rule
   dynvar->theta = 0.5 		  
 *----------------------------------------------------------------------*/ 
irow=-1;
aux = dta*(ONE-dynvar->theta);
for (inode=0;inode<iel;inode++)
{
   for (isd=0;isd<2;isd++)
   {
      irow++;
      for (j=0;j<2;j++)
      {
         eiforce[irow] -= aux*(derxy[j][inode]*(vderxy[isd][j]
	                  +vderxy[j][isd])*c);
      } /* end of loop over j */
   } /* end of loop over isd */
} /* end of loop over inode */

/*----------------------------------------------------------------------*
   Calculate convective forces of time force vector:
                     /
   N(Un)Un  = - dt  |  v * u * grad(u)     d_omega
                   / 
 *----------------------------------------------------------------------*/
irow=-1;
aux = fac*dta;
for (inode=0;inode<iel;inode++)
{
   for (isd=0;isd<2;isd++)
   {
      irow++;
      eiforce[irow] -= funct[inode]*covint[isd]*aux;
   } /* end of loop over isd */
} /* end of loop over irwo */

/*----------------------------------------------------------------------*
   Calculate pressure forces of time force vector:
                      /
    CPn  =    -dt *  |  div(v) * p  d_omega
                    /
REMARK:: C term is defined as (-) another (-) comes from taking 
to the other side and the below term becomes (+)
*------------------------------------------------------------------------*/ 
aux = preint*fac*dta ;
irow=-1;
for (inode=0;inode<iel;inode++)
{
   for (isd=0;isd<2;isd++)
   {
      irow++;
      etforce[irow] += derxy[isd][inode]*aux;
   } /* end of loop over isd */
} /* end of loop over inode */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2pro_calgaltfv */

#endif
/*! @} (documentation module close)*/

/*!----------------------------------------------------------------------
\file
\brief external RHS for fluid3 element

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
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param	 *funct       DOUBLE	      (i)    nat. shape functions
\param   *edeadn      DOUBLE          (i)    ele dead load at n
\param   *edeadng     DOUBLE          (i)    ele dead load at n+1
\param	  fac	      DOUBLE	      (i)    weighting factor
\param	  iel	      INT	      (i)    num. of nodes in ele
\return void

------------------------------------------------------------------------*/
void f3_calgalexfv(
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

fdyn    = alldyn[genprob.numff].fdyn;
/*--------------------------------------------------- set some factors */
facsl = fac*fdyn->thsl;
facsr = fac*fdyn->thsr;
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
   for (isd=0;isd<3;isd++)
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
} /* end of f3_calgalexfv */

/*!---------------------------------------------------------------------
\brief stabilisation part of external forces for vel dofs

<pre>                                                         genk 09/02

In this routine the stabilisation part of the time forces for vel dofs
is calculated:

EULER:
                     /
   + thetas(l,r)*dt |  tau_mu * u * grad(v) * b^   d_omega
                   /

ALE:
                     /
   + thetas(l,r)*dt |  tau_mu * c * grad(v) * b^   d_omega
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


</pre>
\param   *gls         STAB_PAR_GLS    (i)    stabilisation
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
void f3_calstabexfv(
                     STAB_PAR_GLS    *gls,
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
DOUBLE sign=ONE;
DOUBLE aux;
DOUBLE taumu,taump;
DOUBLE c,fvts;
DOUBLE fact[3];

#ifdef DEBUG
dstrc_enter("f3_calstabexfv");
#endif

/*---------------------------------------------------------- initialise */
fdyn   = alldyn[genprob.numff].fdyn;

dsassert(ele->e.f3->stab_type == stab_gls,
         "routine with no or wrong stabilisation called");

/*--------------------------------------------------- set some factors */
taumu = fdyn->tau[0];
taump = fdyn->tau[1];
switch (flag)
{
case 0: /* evaluation at n */
   c = fdyn->thsr;
break;
case 1: /* evaluation at n+1 */
   c = fdyn->thsl;
break;
default:
   c = 0.0;
   dserror("value of flag not valid!!!\n");
}

/*----------------------------------------------------------------------*
   Calculate external/convective stab-forces of time force vector:
EULER:
                     /
   + thetas(l,r)*dt |  tau_mu * u * grad(v) * b^   d_omega
                   /
ALE:
                     /
   + thetas(l,r)*dt |  tau_mu * c * grad(v) * b^   d_omega
                   /
 *----------------------------------------------------------------------*/
if (gls->iadvec!=0)
{
   fact[0] = edead[0]*fac*taumu*c;
   fact[1] = edead[1]*fac*taumu*c;
   fact[2] = edead[2]*fac*taumu*c;
   irow=-1;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1]
          + derxy[2][inode]*velint[2];
      for (isd=0;isd<3;isd++)
      {
         irow++;
         eforce[irow] += aux*fact[isd];
      } /* end loop over isd */
   } /* end loop over inode */
} /* endif (ele->e.f3->iadvec!=0) */

/*----------------------------------------------------------------------*
   Calculate external/viscous stab-forces of time force vector:
                       /
   -/+ thetas(l,r)*dt |  tau_mp * 2*nue * div( eps(v) ) * b  d_omega
                     /
 *----------------------------------------------------------------------*/
if (gls->ivisc!=0 && ihoel!=0)
{
   if (gls->ivisc==2) sign*=-ONE; /* GLS+ stabilisation */

   fvts = fac*visc*taump*sign*c;
   irow = 0;
   for (inode=0;inode<iel;inode++)
   {
      aux = derxy2[0][inode] + derxy2[1][inode] + derxy2[2][inode];
      eforce[irow]   -= ( derxy2[0][inode]*edead[0]  \
                        + derxy2[3][inode]*edead[1]  \
                        + derxy2[4][inode]*edead[2] + aux*edead[0])*fvts;
      eforce[irow+1] -= ( derxy2[3][inode]*edead[0]  \
                        + derxy2[1][inode]*edead[1]  \
                        + derxy2[5][inode]*edead[2] + aux*edead[1])*fvts;
      eforce[irow+2] -= ( derxy2[4][inode]*edead[0]  \
                        + derxy2[5][inode]*edead[1]  \
                        + derxy2[2][inode]*edead[2] + aux*edead[2])*fvts;
      irow += 3;
   } /* end loop over inode */
} /* endif (ele->e.f3->ivisc!=0 && ihoel!=0) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_stabexfv */

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
\param   *gls      STAB_PAR_GLS    (i)    stabilisation
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param  **derxy,      DOUBLE	      (i)    global derivatives
\param   *edead       DOUBLE          (i)    ele dead load at n or n+1
\param   *velint      DOUBLE	      (i)    vel. at integr. point
\param	  fac	      DOUBLE	      (i)    weighting factor
\param	  iel	      INT	      (i)    num. of nodes in ele
\param	  flag	      INT	      (i)    flag for n or n+1
\return void

------------------------------------------------------------------------*/
void f3_calstabexfp(
                     STAB_PAR_GLS    *gls,
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
DOUBLE fact[3];

if (gls->ipres==0) goto end; /* no pressure stabilisation */
fdyn    = alldyn[genprob.numff].fdyn;

/*--------------------------------------------------- set some factors */
taump = fdyn->tau[1];
switch (flag)
{
case 0: /* evaluation at n */
   c = fdyn->thpr;
break;
case 1: /* evaluation at n+1 */
   c = fdyn->thpl;
break;
default:
   c = 0.0;
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
fact[2] = edead[2]*taump*fac*c;
for (inode=0;inode<iel;inode++)
{
   eforce[inode] -= (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1]
                   + derxy[2][inode]*fact[2]);
}

/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_stabexfp */

#endif

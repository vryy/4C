/*!----------------------------------------------------------------------
\file
\brief external RHS for fluid2_xfem element

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_XFEM
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "../fluid2/fluid2.h"
#include "xfem_prototypes.h"
/*!
\addtogroup XFEM
*//*! @{ (documentation module open)*/



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA            *alldyn;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;



static FLUID_DYNAMIC      *fdyn;



/*!---------------------------------------------------------------------
\brief galerkin part of external forces for vel dofs

<pre>                                                            irhan 05/04

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
\param    index	      INT  	      (i)    index for local assembly
\param    DENS	      DOUBLE  	      (i)    fluid density
\return void

------------------------------------------------------------------------*/
void xfem_f2_calgalexfv(
  DOUBLE             *eforce,
  DOUBLE             *funct,
  DOUBLE             *edeadn,
  DOUBLE             *edeadng,
  DOUBLE              fac,
  INT                 iel,
  INT                *index,
  DOUBLE              DENS
  )
{
  DOUBLE  facsl,facsr;
  INT     inode,irow,isd;

#ifdef DEBUG
  dstrc_enter("xfem_f2_calgalexfv");
#endif
/*----------------------------------------------------------------------*/

  /* initialize */
  fdyn = alldyn[genprob.numff].fdyn;

  /* set some factors */
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
  for (inode=0; inode<TWO*iel; inode++)
  {
    irow = index[inode];
    for (isd=0;isd<2;isd++)
    {
      eforce[irow] += funct[inode]*(edeadn[isd]*facsr+edeadng[isd]*facsl)*DENS;
      irow++;
    }
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calgalexfv */



/*!---------------------------------------------------------------------
\brief stabilisation part of external forces for vel dofs

<pre>                                                            irhan 05/04

In this routine the stabilisation part of the time forces for vel dofs
is calculated:

EULER:
                     /
   + thetas(l,r)*dt |  taum_mu * u * grad(v) * b^   d_omega
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
\param    index	      INT  	      (i)    index for local assembly
\param    DENS	      DOUBLE  	      (i)    fluid density
\return void

------------------------------------------------------------------------*/
void xfem_f2_calstabexfv(
  ELEMENT            *ele,
  DOUBLE             *eforce,
  DOUBLE           **derxy,
  DOUBLE           **derxy2,
  DOUBLE            *edead,
  DOUBLE            *velint,
  DOUBLE             fac,
  DOUBLE             visc,
  INT                iel,
  INT                ihoel,
  INT                flag,
  INT               *index,
  DOUBLE             DENS
  )
{
  INT         irow,inode,isd;
  DOUBLE      sign,aux;
  DOUBLE      taumu,taump;
  DOUBLE      c,fvts;
  DOUBLE      fact[2];

  STAB_PAR_GLS  *gls;	/* pointer to GLS stabilisation parameters	*/

#ifdef DEBUG
  dstrc_enter("xfem_f2_calstabexfv");
#endif
/*---------------------------------------------------------------------*/

  /* initialize */
  fdyn = alldyn[genprob.numff].fdyn;
  gls = ele->e.f2->stabi.gls;

  /* set some factors */
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
        dserror("value of flag not valid!!!\n");
  }

/*----------------------------------------------------------------------*
   Calculate external/convective stab-forces of time force vector:
EULER:
                     /
   + thetas(l,r)*dt |  taum_mu * u * grad(v) * b^   d_omega
                   /
 *----------------------------------------------------------------------*/
  if (gls->iadvec!=0)
  {
    fact[0] = edead[0]*fac*taumu*c*DENS*DENS;
    fact[1] = edead[1]*fac*taumu*c*DENS*DENS;
    for (inode=0;inode<TWO*iel;inode++)
    {
      irow = index[inode];
      aux = derxy[0][inode]*velint[0] + derxy[1][inode]*velint[1];
      for (isd=0;isd<2;isd++)
      {
        eforce[irow] += aux*fact[isd];
        irow++;
      }
    }
  }

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
    }

    fvts = fac*visc*taump*sign*c*DENS;
    for (inode=0;inode<TWO*iel;inode++)
    {
      irow = index[inode];
      eforce[irow]   -= ((TWO*derxy2[0][inode] + derxy2[1][inode])*edead[0]
                         + derxy2[2][inode]*edead[1])*fvts;
      eforce[irow+1] -= ((TWO*derxy2[1][inode] + derxy2[0][inode])*edead[1]
                         + derxy2[2][inode]*edead[0])*fvts;
    }
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of xfem_f2_calstabexfv */



/*!---------------------------------------------------------------------
\brief stabilisation part of external forces for pre dofs

<pre>                                                            irhan 05/04

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
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param  **derxy,      DOUBLE	      (i)    global derivatives
\param   *edead       DOUBLE          (i)    ele dead load at n or n+1
\param   *velint      DOUBLE	      (i)    vel. at integr. point
\param	  fac	      DOUBLE	      (i)    weighting factor
\param	  iel	      INT	      (i)    num. of nodes in ele
\param	  flag	      INT	      (i)    flag for n or n+1
\param    DENS	      DOUBLE  	      (i)    fluid density
\return void

------------------------------------------------------------------------*/
void xfem_f2_calstabexfp(
  DOUBLE          *eforce,
  DOUBLE         **derxy,
  DOUBLE          *edead,
  DOUBLE           fac,
  INT              iel,
  INT              flag,
  DOUBLE           DENS
  )
{
  INT        inode;
  DOUBLE     c;
  DOUBLE     taump;
  DOUBLE     fact[2];

#ifdef DEBUG
  dstrc_enter("xfem_f2_calstabexfp");
#endif
/*---------------------------------------------------------------------*/

  /* initialize */
  fdyn = alldyn[genprob.numff].fdyn;

  /* set some factors */
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
        dserror("value of flag not valid!!!\n");
  }

/*----------------------------------------------------------------------*
   Calculate inertia/pressure stab forces of time force vector:
                      /
   - thetas(l,r)*dt  |  tau_mp * grad(q) * b  d_omega
                    /
 *----------------------------------------------------------------------*/
  fact[0] = edead[0]*taump*fac*c*DENS;
  fact[1] = edead[1]*taump*fac*c*DENS;
  for (inode=0;inode<iel;inode++)
  {
    eforce[inode] -= (derxy[0][inode]*fact[0] + derxy[1][inode]*fact[1]);
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of xfem_f2_stabexfp */
/*! @} (documentation module close)*/
#endif

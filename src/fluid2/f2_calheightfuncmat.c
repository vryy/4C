/*!----------------------------------------------------------------------
\file
\brief stiffness matrix for free surface via height function

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID2
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
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

/*!---------------------------------------------------------------------
\brief evaluate galerkin part for vertical height function

<pre>                                                         genk 05/03

In this routine the height function matrices are calculated:

    /
   |  psi * phi    d_gamma_FS
  /

    /
   |  psi * (Ux)_old * phi,x    d_gamma_FS
  /

##################################################
    /
   |  psi * Ux * (phi,x)_old    d_gamma_FS
  /
##################################################

     /
  - |  psi * Uy    d_gamma_FS
   /

    /
   |  tau_SUPG * Ux * psi,x * phi    d_gamma_FS
  /

   /
  |   tau_SUPG * Ux * psi,x * Ux * phi,x  d_gamma_FS
 /

    /
-  |   tau_SUPG * Ux * psi,x * Uy   d_gamma_FS
  /


NOTE: there's only one elemass
      --> Mpsipsi is stored in emass[(3*iel)..(4*iel-1)][(3*iel)..(4*iel-1)]
NOTE: there's only one elestif
      --> Kpsipsi is stored in estif[(3*iel)..(4*iel-1)][(3*iel)..(4*iel-1)]
NOTE: there's only one elestif
      --> Kpsiv is stored in estif[(3*iel)..(4*iel-1)][(0)..(2*iel-1)]

</pre>
\param **emass     DOUBLE     (i/o)  ele mass matrix
\param **estif     DOUBLE     (i/o)  ele stiffnes matrix
\param  *funct     DOUBLE     (i)    nat. shape funcs
\param **derxy     DOUBLE     (i)    global derivatives of shape funcs
\param  *velint    DOUBLE     (i)    velocity at int. point
\param   fac       DOUBLE     (i)    weighting factor
\param   iel       INT        (i)    number of nodes of act. ele
\param   ngnode    IN         (i)    number of nodes of act. edge
\param  *iedgnod   INT        (i)    edge node numbers
\return void

------------------------------------------------------------------------*/
void f2_calmat_vhf(
                        DOUBLE         **emass,
                        DOUBLE         **estif,
                        DOUBLE          *funct,
                        DOUBLE         **derxy,
                        DOUBLE          *velint,
                        DOUBLE           phiderxng,
                        DOUBLE           fac,
                        INT              iel,
                        INT              ngnode,
                        INT             *iedgnod
                         )
{

INT     irow, icol,irn,icn;
INT     nd;
DOUBLE  c,tau_supg,tau_dc;
FLUID_DYNAMIC *fdyn;

#ifdef DEBUG
dstrc_enter("f2_calmat_vhf");
#endif

fdyn    = alldyn[genprob.numff].fdyn;
nd = NUMDOF_FLUID2*iel;

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Mpsipsi:
    /
   |  psi * phi    d_gamma_FS
  /
 *----------------------------------------------------------------------*/
for(icn=0;icn<ngnode;icn++)
{
   icol = nd+iedgnod[icn];
   for(irn=0;irn<ngnode;irn++)
   {
      irow = nd+iedgnod[irn];
      emass[irow][icol]     += funct[irn]*funct[icn]*fac;
   } /* end loop over irn */
} /* end loop over icn */


/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Kpsipsi:
    /
   |  psi * (Ux)_old * phi,x    d_gamma_FS
  /
 *----------------------------------------------------------------------*/
for(icn=0;icn<ngnode;icn++)
{
   icol = nd+iedgnod[icn];
   for(irn=0;irn<ngnode;irn++)
   {
      irow = nd+iedgnod[irn];
      estif[irow][icol]     += velint[0]*funct[irn]*derxy[0][icn]*fac;
   } /* end loop over irn */
} /* end loop over icn */

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Kpsiv:
    /
   |  psi * Ux * (phi,x)_old    d_gamma_FS
  /
 *----------------------------------------------------------------------*/
for(icn=0;icn<ngnode;icn++)
{
   icol = NUM_F2_VELDOF*iedgnod[icn];
   for(irn=0;irn<ngnode;irn++)
   {
      irow = nd+iedgnod[irn];
      estif[irow][icol]     += funct[irn]*funct[icn]*phiderxng*fac;
   } /* end loop over irn */
} /* end loop over icn */

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Kpsiv:
     /
  - |  psi * Uy    d_gamma_FS
   /
 *----------------------------------------------------------------------*/
for(icn=0;icn<ngnode;icn++)
{
   icol = NUM_F2_VELDOF*iedgnod[icn]+1;
   for(irn=0;irn<ngnode;irn++)
   {
      irow = nd+iedgnod[irn];
      estif[irow][icol]     -= funct[irn]*funct[icn]*fac;
   } /* end loop over irn */
} /* end loop over icn */


/*------------------------------------------------------- STABILISATION */
if (fdyn->hf_stab>0)
{
   tau_supg = fdyn->tau[3];
   tau_dc   = fdyn->tau[4];

/*----------------------------------------------------------------------*
   Calculate full stabilisation part of matrix Mpsipsi:
    /
   |  tau_SUPG * Ux * psi,x * phi    d_gamma_FS
  /
 *----------------------------------------------------------------------*/
   c=fac*tau_supg*velint[0];
   for(icn=0;icn<ngnode;icn++)
   {
      icol = nd+iedgnod[icn];
      for(irn=0;irn<ngnode;irn++)
      {
         irow = nd+iedgnod[irn];
         emass[irow][icol]     += derxy[0][irn]*funct[icn]*c;
      } /* end loop over irn */
   } /* end loop over icn */

/*----------------------------------------------------------------------*
   Calculate full stabilisation part of matrix Kpsipsi:
   /
  |   tau_SUPG * Ux * psi,x * Ux_old * phi,x  d_gamma_FS
 /
 *----------------------------------------------------------------------*/
   c=fac*tau_supg*velint[0]*velint[0];
   for(icn=0;icn<ngnode;icn++)
   {
      icol = nd+iedgnod[icn];
      for(irn=0;irn<ngnode;irn++)
      {
         irow = nd+iedgnod[irn];
         estif[irow][icol]     += derxy[0][irn]*derxy[0][icn]*c;
      } /* end loop over irn */
   } /* end loop over icn */

/*----------------------------------------------------------------------*
   Calculate full stabilisation part of matrix Kpsipsi:
   /
  |   tau_SUPG * Ux * psi,x * Ux * (phi,x)_old  d_gamma_FS
 /
 *----------------------------------------------------------------------*/
   c=fac*tau_supg*velint[0]*phiderxng;
   for(icn=0;icn<ngnode;icn++)
   {
      icol = NUM_F2_VELDOF*iedgnod[icn];
      for(irn=0;irn<ngnode;irn++)
      {
         irow = nd+iedgnod[irn];
         estif[irow][icol]     += derxy[0][irn]*funct[icn]*c;
      } /* end loop over irn */
   } /* end loop over icn */


/*----------------------------------------------------------------------*
   Calculate full stabilisation part of matrix Kpsiv:
    /
-  |   tau_SUPG * Ux * psi,x * Uy   d_gamma_FS
  /
 *----------------------------------------------------------------------*/
   c=fac*tau_supg*velint[0];
   for(icn=0;icn<ngnode;icn++)
   {
      icol = NUM_F2_VELDOF*iedgnod[icn]+1;
      for(irn=0;irn<ngnode;irn++)
      {
         irow = nd+iedgnod[irn];
         estif[irow][icol]     -= derxy[0][irn]*funct[icn]*c;
      } /* end loop over irn */
   } /* end loop over icn */
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calgalmat_vhf */


/*!---------------------------------------------------------------------
\brief evaluate stabilisation part for vertical height function

<pre>                                                         genk 05/03

In this routine the galerkin part of matrix Mpsipsi is calculated:

    /
   |  tau_SUPG * Ux * psi,x * phi    d_gamma_FS
  /

   /
  |   tau_SUPG * Ux * psi,x * Ux * phi,x  d_gamma_FS
 /

    /
-  |   tau_SUPG * Ux * psi,x * Uy   d_gamma_FS
  /

    /
   |  tau_DC * psi,x * phi,x    d_gamma_FS
  /

NOTE: there's only one elemass
      --> Mpsipsi is stored in emass[(3*iel)..(4*iel-1)][(3*iel)..(4*iel-1)]
NOTE: there's only one elestif
      --> Kpsipsi is stored in estif[(3*iel)..(4*iel-1)][(3*iel)..(4*iel-1)]
NOTE: there's only one elestif
      --> Kpsiv is stored in estif[(3*iel)..(4*iel-1)][(0)..(2*iel-1)]

</pre>
\param **emass     DOUBLE           (i/o)  ele mass matrix
\param **estif     DOUBLE           (i/o)  ele stiffnes matrix
\param  *funct     DOUBLE           (i)    nat. shape funcs
\param **derxy     DOUBLE           (i)    global derivatives of shape funcs
\param  *velint    DOUBLE           (i)    velocity at int. point
\param   fac       DOUBLE           (i)    weighting factor
\param   iel       INT              (i)    number of nodes of act. ele
\param   ngnode    INT              (i)    number of nodes of act. edge
\param  *iedgnod   INT              (i)    edge node numbers
\return void

------------------------------------------------------------------------*/
void f2_calstabmat_vhf(
                           ELEMENT         *ele,
		           DOUBLE         **emass,
                           DOUBLE         **estif,
		           DOUBLE          *funct,
                           DOUBLE         **derxy,
                           DOUBLE          *velint,
		           DOUBLE           fac,
		           INT              iel,
	  	           INT              ngnode,
                           INT             *iedgnod
                         )
{

INT     irow, icol,irn,icn;
INT     nd;
DOUBLE  c;
DOUBLE  tau_dc, tau_supg;
FLUID_DYNAMIC *fdyn;

#ifdef DEBUG
dstrc_enter("f2_calgalmat_vhf");
#endif

fdyn    = alldyn[genprob.numff].fdyn;
nd = NUMDOF_FLUID2*iel;

#if 1
if (iel==4)
{
   dserror("heightfunc stabilisation");
   tau_supg = fdyn->tau[4];
/*----------------------------------------------------------------------*
   Calculate full stabilisation part of matrix Mpsipsi:
    /
   |  tau_SUPG * Ux * psi,x * phi    d_gamma_FS
  /
 *----------------------------------------------------------------------*/
   c=fac*tau_supg*velint[0];
   for(icn=0;icn<ngnode;icn++)
   {
      icol = nd+iedgnod[icn];
      for(irn=0;irn<ngnode;irn++)
      {
         irow = nd+iedgnod[irn];
         emass[irow][icol]     += funct[icn]*derxy[0][irn]*c;
      } /* end loop over irn */
   } /* end loop over icn */

/*----------------------------------------------------------------------*
   Calculate full stabilisation part of matrix Kpsipsi:
   /
  |   tau_SUPG * Ux * psi,x * Ux * phi,x  d_gamma_FS
 /
 *----------------------------------------------------------------------*/
   c=fac*tau_supg*velint[0]*velint[0];
   for(icn=0;icn<ngnode;icn++)
   {
      icol = nd+iedgnod[icn];
      for(irn=0;irn<ngnode;irn++)
      {
         irow = nd+iedgnod[irn];
         estif[irow][icol]     += derxy[0][irn]*derxy[0][icn]*c;
      } /* end loop over irn */
   } /* end loop over icn */

/*----------------------------------------------------------------------*
   Calculate full stabilisation part of matrix Kpsiv:
    /
-  |   tau_SUPG * Ux * psi,x * Uy   d_gamma_FS
  /
 *----------------------------------------------------------------------*/
   c=fac*tau_supg*velint[0];
   for(icn=0;icn<ngnode;icn++)
   {
      icol = NUM_F2_VELDOF*iedgnod[icn]+1;
      for(irn=0;irn<ngnode;irn++)
      {
         irow = nd+iedgnod[irn];
         estif[irow][icol]     -= derxy[0][irn]*funct[icn]*c;
      } /* end loop over irn */
   } /* end loop over icn */
} /* endif (ele->e.f2->ifssu!=0) */

if (iel==4)
{
   dserror("heightfunc stabilisation");
   tau_dc   = fdyn->tau[3];
/*----------------------------------------------------------------------*
   Calculate full stabilisation part of matrix Kpsipsi:
    /
   |  tau_DC * psi,x * phi,x    d_gamma_FS
  /
 *----------------------------------------------------------------------*/
   c=fac*tau_dc;
   for(icn=0;icn<ngnode;icn++)
   {
      icol = nd+iedgnod[icn];
      for(irn=0;irn<ngnode;irn++)
      {
         irow = nd+iedgnod[irn];
         estif[irow][icol]     += derxy[0][irn]*derxy[0][icn]*c;
      } /* end loop over irn */
   } /* end loop over icn */
}/*endif (ele->e.f2->ifsdc!=0) */

#endif

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calgalmat_vhf */


#endif
/*! @} (documentation module close)*/

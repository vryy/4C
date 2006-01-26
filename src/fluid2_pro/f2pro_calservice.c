/*!----------------------------------------------------------------------
\file
\brief service routines for fluid2 element

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
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
\brief set all arrays for element calculation

<pre>                                                         basol 11/02

get the element velocities and the pressure at different times

</pre>
\param   *elevel      ELEMENT	   (i)    actual element for velocity
\param   *elepre      ELEMENT	   (i)    actual element for pressure
\param   **xyze       DOUBLE	   (o)    coordinates of the element
\param   **eveln      DOUBLE	   (o)    ele velocities at time n
\param   *epren       DOUBLE	   (o)    ele pressures at time n
\param   *ipos                     (i)    node array positions
\return void

------------------------------------------------------------------------*/
void f2pro_calset(
  ELEMENT         *ele,
  DOUBLE         **xyze,
  DOUBLE         **eveln,
  DOUBLE         **evelng,
  DOUBLE         **evhist,
  DOUBLE          *epren,
  ARRAY_POSITION  *ipos
  )
{
  INT i, veln, velnp, hist;
  INT numpdof;
  NODE *actnode;                /* actual node for element */

#ifdef DEBUG
  dstrc_enter("f2pro_calset");
#endif

  veln  = ipos->veln;
  velnp = ipos->velnp;
  hist  = ipos->hist;

  /*--------------------------------------------- set element coordinates */
  for (i=0;i<ele->numnp;i++)
  {
    xyze[0][i]=ele->node[i]->x[0];
    xyze[1][i]=ele->node[i]->x[1];
  }

  /*---------------------------------------------------------------------*
   | position of the different solutions:                                |
   | node->sol_incement: solution history used for calculations          |
   |       sol_increment[ipos][i]: solution at some time level           |
   |       ipos->velnp .. solution at time n+1                           |
   |       ipos->veln  .. solution at time n                             |
   |       ipos->velnm .. solution at time n-1                           |
   *---------------------------------------------------------------------*/

  /* -> computation of time forces -------------------
     -> velocities and pressure at (n) are needed ----*/

  for(i=0;i<ele->numnp;i++)
  {
    actnode=ele->node[i];

    /*------------------------------------- set element velocities at (n) */
    eveln[0][i]=actnode->sol_increment.a.da[veln][0];
    eveln[1][i]=actnode->sol_increment.a.da[veln][1];

    /*----------------------------------- set element velocities (n+gamma) */
    evelng[0][i]=actnode->sol_increment.a.da[velnp][0];
    evelng[1][i]=actnode->sol_increment.a.da[velnp][1];

    /*--------------------------------------- set vel. histories at (n) ---*/
    evhist[0][i] = actnode->sol_increment.a.da[hist][0];
    evhist[1][i] = actnode->sol_increment.a.da[hist][1];
  }

  switch (ele->e.f2pro->dm)
  {
  case dm_q2pm1:
    numpdof = 3;
    break;
  case dm_q1p0:
    numpdof = 1;
    break;
  default:
    dserror("unsupported discretization mode %d", ele->e.f2pro->dm);
  }

  for (i=0; i<numpdof; ++i)
  {
    /*---------------------------------------------- set pressures (n+1) ---*/
    epren[i]   = ele->e.f2pro->press[i];
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of f2pro_calset */


/*----------------------------------------------------------------------*/
/*!
  \brief calculate shape functions and derivatives of the pressure

  \param pfunct       (o) shape function array
  \param pderiv       (o) shape function derivative array
  \param r            (i) natural coordinate
  \param s            (i) natural coordinate
  \param dm           (i) discontinuous mode
  \param numpdof      (i) number of pressure dofs

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void f2pro_prec(DOUBLE* pfunct, DOUBLE** pderiv, DOUBLE r, DOUBLE s, DISMODE dm, INT* numpdof)
{
#ifdef DEBUG
  dstrc_enter("f2pro_prec");
#endif

  switch (dm)
  {
  case dm_q2pm1:
    /* discontinuous pressure on quadratic velocity */
    /* 9/3 */

    *numpdof = 3;

    pfunct[0]=1;
    pfunct[1]=r;
    pfunct[2]=s;

    pderiv[0][0]= 0;
    pderiv[1][0]= 0;

    pderiv[0][1]= 1;
    pderiv[1][1]= 0;

    pderiv[0][2]= 0;
    pderiv[1][2]= 1;
    break;
  case dm_q2q1:
  {
    DOUBLE rp;
    DOUBLE rm;
    DOUBLE sp;
    DOUBLE sm;

    rp=1.+r;
    rm=1.-r;
    sp=1.+s;
    sm=1.-s;

    /* Taylor-Hood */

    *numpdof = 4;

    pfunct[0]=.25*rp*sp;
    pfunct[1]=.25*rm*sp;
    pfunct[2]=.25*rm*sm;
    pfunct[3]=.25*rp*sm;

    pderiv[0][0]= .25*sp;
    pderiv[1][0]= .25*rp;

    pderiv[0][1]=-.25*sp;
    pderiv[1][1]= .25*rm;

    pderiv[0][2]=-.25*sm;
    pderiv[1][2]=-.25*rm;

    pderiv[0][3]= .25*sm;
    pderiv[1][3]=-.25*rp;
    break;
  }
  case dm_q1p0:
    /* constant pressure */
    *numpdof = 1;

    pfunct[0]=1;

    pderiv[0][0]= 0;
    pderiv[1][0]= 0;

    break;
  default:
    dserror("unsupported discretisation mode %d", dm);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


#endif
/*! @} (documentation module close)*/

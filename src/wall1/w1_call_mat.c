/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_call_mat' which selects the proper
       material law for a wall element
       contains the routine 'w1_getdensity' which gets the density
       out of the material law

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | select proper material law                               al 01/02    |
 *----------------------------------------------------------------------*/
void w1_call_mat(ELEMENT   *ele,
                 MATERIAL  *mat,
                 WALL_TYPE wtype,
                 DOUBLE **bop,
                 DOUBLE  *gop,
                 DOUBLE  *alpha,
                 INT ip,
                 DOUBLE *stress,
                 DOUBLE **d,
                 INT istore,/* controls storing of new stresses to wa */
                 INT newval)/* controls evaluation of new stresses    */
{
INT i;
INT j;
DOUBLE disd[4];
/*----------------------------------------------------------------------*/
DOUBLE strain[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1_call_mat");
#endif
/*------------------------------------------------ call material law ---*/
  switch(mat->mattyp)
  {
  case m_stvenant:/*------------------------------- linear elastic ---*/
    w1_mat_linel(mat->m.stvenant->youngs,
                 mat->m.stvenant->possionratio,
                 wtype,
                 d);
    w1_disd(ele,bop,gop,alpha,wtype,disd);
    w1_eps (disd,ele->e.w1->wtype,strain);
    for (i=0; i<4; i++) stress[i] = 0.0;
    for (i=0; i<4; i++) for (j=0; j<4; j++) stress[i] += d[i][j]*strain[j];
  break;
  case m_stvenpor:/*------------------------ porous linear elastic ---*/
    w1_mat_stvpor(mat, ele->e.w1->elewa->matdata, wtype, d);
    w1_disd(ele,bop,gop,alpha,wtype,disd);
    w1_eps (disd,ele->e.w1->wtype,strain);
    for (i=0; i<4; i++) stress[i] = 0.0;
    for (i=0; i<4; i++) for (j=0; j<4; j++) stress[i] += d[i][j]*strain[j];
  break;
  case m_pl_mises:/*--------------------- von mises material law ---*/
    w1_mat_plast_mises(mat->m.pl_mises->youngs,
                       mat->m.pl_mises->possionratio,
                       mat->m.pl_mises->Sigy,
                       mat->m.pl_mises->Hard,
                       mat->m.pl_mises->GF,
                       mat->m.pl_mises->betah,
                       ele,
                       wtype,
                       bop,
                       gop,
                       alpha,
                       ip,
                       stress,
                       d,
                       istore,
                       newval);
  break;
#ifdef D_MAT
  case m_pl_mises_3D:/*--------------------- von mises material law -> 3D---*/
    w1_mat_plast_mises_3D(mat->m.pl_mises->youngs,
                          mat->m.pl_mises->possionratio,
                          mat->m.pl_mises->Sigy,
                          mat->m.pl_mises->Hard,
                          mat->m.pl_mises->GF,
                          mat->m.pl_mises->betah,
                          ele,
                          wtype,
                          bop,
                          gop,
                          alpha,
                          ip,
                          stress,
                          d,
                          istore,
                          newval);
  break;
#endif
  case m_pl_dp:/*------------------- drucker prager material law ---*/
    w1_mat_plast_dp(   mat->m.pl_dp->youngs,
                       mat->m.pl_dp->possionratio,
                       mat->m.pl_dp->Sigy,
                       mat->m.pl_dp->Hard,
                       mat->m.pl_dp->PHI,
                       ele,
                       wtype,
                       bop,
                       gop,
                       alpha,
                       ip,
                       stress,
                       d,
                       istore,
                       newval);
  break;
  case m_pl_epc:/*---------- elastoplastic concrete material law ---*/
    w1_mat_plast_epc(  mat->m.pl_epc->youngs      ,
                       mat->m.pl_epc->possionratio,
                       mat->m.pl_epc->ftm         ,
                       mat->m.pl_epc->fcm         ,
                       mat->m.pl_epc->gt          ,
                       mat->m.pl_epc->gc          ,
                       mat->m.pl_epc->gamma1      ,
                       mat->m.pl_epc->gamma2      ,
                       mat->m.pl_epc->nstiff      ,
                       mat->m.pl_epc->maxreb      ,
                       mat->m.pl_epc->reb_area    ,
                       mat->m.pl_epc->reb_ang     ,
                       mat->m.pl_epc->reb_so      ,
                       mat->m.pl_epc->reb_ds      ,
                       mat->m.pl_epc->reb_rgamma  ,
                       ele,
                       wtype,
                       bop,
                       gop,
                       alpha,
                       ip,
                       stress,
                       d,
                       istore,
                       newval);
  break;
#ifdef D_MAT
  case m_pl_epc3D:/*------ elastoplastic concrete material law 3D --- sh 12/03 -*/
    w1_mat_plast_epc3D(mat->m.pl_epc->youngs      ,
                       mat->m.pl_epc->possionratio,
                       mat->m.pl_epc->ftm         ,
                       mat->m.pl_epc->fcm         ,
                       mat->m.pl_epc->gt          ,
                       mat->m.pl_epc->gc          ,
                       mat->m.pl_epc->gamma1      ,
                       mat->m.pl_epc->gamma2      ,
                       mat->m.pl_epc->gamma3      ,
                       mat->m.pl_epc->gamma4      ,
                       ele,
                       wtype,
                       bop,
                       gop,
                       alpha,
                       ip,
                       stress,
                       d,
                       istore,
                       newval);
  break;
#endif
  case m_dam_mp:/*--------------------- isotropic damage ---*/
    w1_mat_dam_mp(mat->m.dam_mp->youngs,
                   mat->m.dam_mp->nue,
                   mat->m.dam_mp->kappa_0,
                   mat->m.dam_mp->alpha,
                   mat->m.dam_mp->beta,
                   ele,
                   wtype,
                   bop,
                   gop,
                   alpha,
                   ip,
                   stress,
                   d,
                   istore,
                   newval);
  break;
  case m_damage:   /*------------------------------------- CDM law ---*/
    w1_mat_damage(mat->m.damage->youngs,
                       mat->m.damage->possionratio,
                       mat->m.damage->Equival,
                       mat->m.damage->Damtyp,
                       mat->m.damage->Kappa_0,
                       mat->m.damage->Kappa_m,
                       mat->m.damage->Alpha,
                       mat->m.damage->Beta,
                       mat->m.damage->k_fac,
                       ele,
                       wtype,
                       bop,
                       gop,
                       alpha,
                       ip,
                       stress,
                       d,
                       istore,
                       newval);
  break;
  default:
    dserror(" unknown type of material law");
  break;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_call_mat */
/*----------------------------------------------------------------------*/

#if 1
/*----------------------------------------------------------------------*/
/* get density out of material law                        ah 06/02      */
/*----------------------------------------------------------------------*/
void w1_getdensity(MATERIAL   *mat, DOUBLE *density)
{
#ifdef DEBUG
dstrc_enter("w1_getdensity");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------ switch material type */
switch(mat->mattyp)
{
case m_stvenant:/*-------------------------------------- linear elastic */
   *density = mat->m.stvenant->density;
break;
case m_neohooke:/*------------------------------ kompressible neo-hooke */
   *density = mat->m.neohooke->density;
break;
case m_stvenpor:/*------------------------ porous linear elastic ---*/
   *density = mat->m.stvenpor->density;
break;
case m_pl_mises:/*--------------------------- von mises material law ---*/
  dserror("Ilegal typ of material for this element");
break;
case m_pl_mises_3D:/*-------------Stefan's von mises 3D material law ---*/
  dserror("Ilegal typ of material for this element");
break;
case m_pl_dp:/*------------------------- drucker prager material law ---*/
   dserror("Ilegal typ of material for this element");
break;
default:
   dserror("Ilegal typ of material for this element");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_getdensity */
#endif
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/

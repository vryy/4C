/*----------------------------------------------------------------------------*/
/*! \file
\brief shell8

\level 1


*/
/*---------------------------------------------------------------------------*/
#include "../../headers/standardtypes.h"
#include "shell8.h"

#include "../../drt_lib/drt_dserror.H"

#if 0
/*----------------------------------------------------------------------*
 | call material laws                                     m.gee 6/01    |
 *----------------------------------------------------------------------*/
void s8_tmat(ELEMENT    *ele,
             MATERIAL   *mat,
             DOUBLE     *stress,
             DOUBLE     *strain,
             DOUBLE    **C,
             DOUBLE    **gmkovc,
             DOUBLE    **gmkonc,
             DOUBLE    **gmkovr,
             DOUBLE    **gmkonr,
             DOUBLE    **gkovc,
             DOUBLE    **gkonc,
             DOUBLE    **gkovr,
             DOUBLE    **gkonr,
             DOUBLE      detc,
             DOUBLE      detr,
             DOUBLE      e3,
             INT         option,
             INT         ngauss)
{
DOUBLE C4[3][3][3][3];


#ifdef DEBUG
dstrc_enter("s8_tmat");
#endif


/*----------------------------------------------------------------------*/
/*------------------------------------------------ switch material type */
switch(mat->mattyp)
{
case m_stvenant:/*------------------------ st.venant-kirchhoff-material */
   s8_eps(strain,gmkovc,gmkovr);
   /*=============================== make material in curvilinear bases */
   s8_mat_linel(mat->m.stvenant,gmkonr,C);
   s8_mat_stress1(stress,strain,C);
   /*===================================================================*/
break;
case m_neohooke:/*------------------------------ kompressible neo-hooke */
   s8_eps(strain,gmkovc,gmkovr);
   s8_mat_neohooke(mat->m.neohooke,stress,C,gmkonr,gmkonc,detr,detc);
break;
case m_compogden:/*--------------------------------- kompressible ogden */
   /*-------------------------------------- make Green-Lagrange strains */
   s8_eps(strain,gmkovc,gmkovr);
   /*----------------------------- call compressible ogden material law */
   /*         Ogden hyperelasticity without deviatoric-volumetric split */
/*
   s8_mat_ogden_coupled(mat->m.compogden,stress,C4,gkonr,gmkovc);
*/
   /*            Ogden hyperelasticity with deviatoric-volumetric split */
   s8_mat_ogden_uncoupled2(mat->m.compogden,stress,C4,gkonr,gmkovc);
   /* PK2 stresses are cartesian ->  return stresses to curvilinear bases */
   s8_kon_cacu(stress,gkonr);
   /*---------------- C4 is cartesian -> return C4 to curvilinear bases */
   s8_4kon_cacu(C4,gkonr);
   /*---------------------- sort material tangent from tensor to matrix */
   s8_c4_to_C2(C4,C);
   /*-------------------------------------------------------------------*/
break;
case m_viscohyper:/*-------------------------viscous kompressible ogden */
   /*-------------------------------------- make Green-Lagrange strains */
   s8_eps(strain,gmkovc,gmkovr);
   /*----------------------------- call compressible ogden material law */
   /*    viscous Ogden hyperelasticity with deviatoric-volumetric split */
   s8_mat_ogden_viscous(ele,mat->m.viscohyper,stress,C4,gkonr,gmkovc,ngauss);
   /* PK2 stresses are cartesian ->  return stresses to curvilinear bases */
   s8_kon_cacu(stress,gkonr);
   /*---------------- C4 is cartesian -> return C4 to curvilinear bases */
   s8_4kon_cacu(C4,gkonr);
   /*---------------------- sort material tangent from tensor to matrix */
   s8_c4_to_C2(C4,C);
   /*-------------------------------------------------------------------*/
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
} /* end of s8_tmat */
#endif /* if 0 */

/*----------------------------------------------------------------------*
 | get density out of material law                        m.gee 2/02    |
 *----------------------------------------------------------------------*/
void s8_getdensity(MATERIAL *mat, DOUBLE *density)
{
  /*----------------------------------------------------------------------*/
  /*------------------------------------------------ switch material type */
  switch (mat->mattyp)
  {
    case m_stvenant: /*-------------------------- ST.VENANT-KIRCHHOFF-MATERIAL */
      *density = mat->m.stvenant->density;
      break;
    case m_neohooke: /*------------------------------ kompressible neo-hooke */
      *density = mat->m.neohooke->density;
      break;
    case m_compogden: /*-------------------------------- kompressible ogden */
      *density = mat->m.compogden->density;
      break;
    case m_viscohyper: /*---------------------------viscos compressible ogden */
      *density = mat->m.viscohyper->density;
      break;
    default:
      dserror("Ilegal typ of material for this element");
      break;
  }
  /*----------------------------------------------------------------------*/
  return;
} /* end of s8_getdensity */

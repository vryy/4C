/*!----------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0771 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
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
             INT         option)
{
INT    i;
DOUBLE C4[3][3][3][3];
DOUBLE stress2[6],strain2[6];
DOUBLE E;
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
   
   
   
   /*----------------------------------- transform strains to cartesian */
   /*s8_kov_cuca(strain,gkonr);  */
   
   /*=========================== in here is material in cartesian bases */
   /*-------------------------------------------------- do material law */
   /*s8_mat_linel_cart(mat->m.stvenant,C4,C,strain); */
   /*----------------------------------------------- cartesian stresses */
   /*s8_mat_stress1(stress,strain,C); */
   /*===================================================================*/
   
   /*------------------------------ return strains to curvilinear bases */
   /*s8_kov_cacu(strain,gkovr); */
   /*----------------------------------- return C4 to curvilinear bases */
   /*s8_4kon_cacu(C4,gkonr); */
   /*-------------------- return material tangent from tensor to matrix */
   /*s8_c4_to_C2(C4,C); */
   /*----------------------------- return stresses to curvilinear bases */
   /*s8_kon_cacu(stress,gkonr);*/



break;
case m_neohooke:/*------------------------------ kompressible neo-hooke */
   s8_eps(strain,gmkovc,gmkovr); 
   s8_mat_neohooke(mat->m.neohooke,stress,C,gmkonr,gmkonc,detr,detc);
break;
case m_compogden:/*--------------------------------- kompressible ogden */
   /*-------------------------------------- make Green-Lagrange strains */
   s8_eps(strain,gmkovc,gmkovr); 
   /*================================== temporary calculating st.venant */
   /*E = 0.0;
   for (i=0; i<3; i++) E += (mat->m.compogden->alfap[i])*(mat->m.compogden->mup[i]);
   E *= (1.0+(mat->m.compogden->nue)); */


   /* This is venant in curvilinear coodinates 
   s8_mat_lineltmp(E,mat->m.compogden->nue,gmkonr,C);
   s8_mat_stress1(stress,strain,C);
   */


   /*----------------------------------- transform strains to cartesian */
   /*s8_kov_cuca(strain,gkonr);  */
   /*-------------------------------------------------- do material law */
   /*s8_mat_linel_carttmp(E,mat->m.compogden->nue,C4); */
   /*-------------------- return material tangent from tensor to matrix */
   /*s8_c4_to_C2(C4,C); */
   /*----------------------------------------------- cartesian stresses */
   /*s8_mat_stress1(stress,strain,C); */
   /*------------------------------ return strains to curvilinear bases */
   /*s8_kov_cacu(strain,gkovr); */
   /*----------------------------------- return C4 to curvilinear bases */
   /*s8_4kon_cacu(C4,gkonr); */
   /*-------------------- return material tangent from tensor to matrix */
   /*s8_c4_to_C2(C4,C); */
   /*----------------------------- return stresses to curvilinear bases */
   /*s8_kon_cacu(stress,gkonr); */
   
   /*===================================================================*/
   
   /*----------------------------- call compressible ogden material law */
   /* strain2 was only used for testing purposes and unused at the moment */
   s8_mat_ogden_coupled(mat->m.compogden,stress,strain2,C4,
                        gkovr,gkonr,gkovc,gkonc,gmkovr,gmkonr,gmkovc,gmkonc);
   /*s8_mat_ogden_uncoupled(mat->m.compogden,stress,C4,
                            gkovr,gkonr,gkovc,gkonc,gmkovr,gmkonr,gmkovc,gmkonc);*/
   /*------------------------------ return strains to curvilinear bases */
   /* strain2 was only used for testing purposes and unused at the moment */
   /*s8_kov_cacu(strain2,gkovr);*/
   /* PK2 stresses are cartesian ->  return stresses to curvilinear bases */
   s8_kon_cacu(stress,gkonr);
   /*---------------- C4 is cartesian -> return C4 to curvilinear bases */
   s8_4kon_cacu(C4,gkonr);
   /*---------------------- sort material tangent from tensor to matrix */
   s8_c4_to_C2(C4,C);
   
   
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



/*----------------------------------------------------------------------*
 | get density out of material law                        m.gee 2/02    |
 *----------------------------------------------------------------------*/
void s8_getdensity(MATERIAL   *mat, DOUBLE *density)
{
#ifdef DEBUG 
dstrc_enter("s8_getdensity");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------ switch material type */
switch(mat->mattyp)
{
case m_stvenant:/*-------------------------- ST.VENANT-KIRCHHOFF-MATERIAL */
   *density = mat->m.stvenant->density;
break;
case m_neohooke:/*------------------------------ kompressible neo-hooke */
   *density = mat->m.neohooke->density;
break;
case m_compogden:/*-------------------------------- kompressible ogden */
   *density = mat->m.compogden->density;
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
} /* end of s8_getdensity */
#endif

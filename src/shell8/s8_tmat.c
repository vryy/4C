#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 | call material laws                                     m.gee 6/01    |
 *----------------------------------------------------------------------*/
void s8_tmat(ELEMENT    *ele,
             MATERIAL   *mat,
             double     *stress,
             double     *strain,
             double    **C,
             double    **gmkovc,
             double    **gmkonc,
             double    **gmkovr,
             double    **gmkonr,
             double    **gkovc,
             double    **gkonc,
             double    **gkovr,
             double    **gkonr,
             double      detc,
             double      detr,
             double      e3,
             int         option)
{
int i,j;
#ifdef DEBUG 
dstrc_enter("s8_tmat");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------ switch material type */
switch(mat->mattyp)
{
case m_stvenant:/*------------------------ st.venant-kirchhoff-material */
   s8_eps(strain,gmkovc,gmkovr); 
   s8_mat_linel(mat->m.stvenant,gmkonr,C);
   s8_mat_stress1(stress,strain,C);
break;
case m_neohooke:/*------------------------------ kompressible neo-hooke */
   s8_eps(strain,gmkovc,gmkovr); 
   s8_mat_neohooke(mat->m.neohooke,stress,C,gmkonr,gmkonc,detr,detc);
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
void s8_getdensity(MATERIAL   *mat, double *density)
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

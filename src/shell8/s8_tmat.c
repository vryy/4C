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
#ifdef DEBUG 
dstrc_enter("s8_tmat");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------ switch material type */
switch(mat->mattyp)
{
case m_stvenant:/*-------------------------- ST.VENANT-KIRCHHOFF-MATERIAL */
   s8_eps(strain,gmkovc,gmkovr); 
   s8_mat_linel(mat->m.stvenant,gmkonr,C);
   s8_mat_stress1(stress,strain,C);
break;
case m_neohooke:/*------------------------------ kompressible neo-hooke */
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

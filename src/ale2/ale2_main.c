#include "../headers/standardtypes.h"
#include "ale2.h"
/*----------------------------------------------------------------------*
 |                                                         mn 06/02     |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*----------------------------------------------------------------------*
 | main ale2    control routine                            mn 06/02     |
 |                                                                      |
 *----------------------------------------------------------------------*/
void ale2(     PARTITION   *actpart,
               INTRA       *actintra,
               ELEMENT     *ele,
               ARRAY       *estif_global,
               CALC_ACTION *action)
{
#ifdef D_ALE
int  i;
ALE2_DATA     actdata;
MATERIAL    *actmat;

#ifdef DEBUG 
dstrc_enter("ale2");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------- switch to do option */
switch (*action)
{
/*------------------------------------------- init the element routines */
case calc_ale_init:
   ale2_static_ke(NULL,NULL,NULL,NULL,1);
break;
/*----------------------------------- calculate linear stiffness matrix */
case calc_ale_stiff:
   actmat = &(mat[ele->mat-1]);
   ale2_static_ke(ele,&actdata,actmat,estif_global,0);
break;
/*----------------------------------------------------------- do nothig */
case calc_ale_rhs:
break;
/*-------------------------------------------------------------- defaul */
default:
   dserror("action unknown");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif
return; 
} /* end of ale3 */

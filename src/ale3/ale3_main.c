#include "../headers/standardtypes.h"
#include "ale3.h"
/*----------------------------------------------------------------------*
 |                                                         mn 06/02     |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*----------------------------------------------------------------------*
 | main ale3    control routine                            mn 06/02     |
 |                                                                      |
 *----------------------------------------------------------------------*/
void ale3(     PARTITION   *actpart,
               INTRA       *actintra,
               ELEMENT     *ele,
               ARRAY       *estif_global,
               CALC_ACTION *action)
{
#ifdef D_ALE
int  i;
ALE3_DATA     actdata;
MATERIAL    *actmat;

#ifdef DEBUG 
dstrc_enter("ale1");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------- switch to do option */
switch (*action)
{
/*------------------------------------------- init the element routines */
case calc_ale_init:
   ale_static_ke(NULL,NULL,NULL,NULL,1);
break;
/*----------------------------------- calculate linear stiffness matrix */
case calc_ale_stiff:
   actmat = &(mat[ele->mat-1]);
   ale_static_ke(ele,&actdata,actmat,estif_global,0);
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

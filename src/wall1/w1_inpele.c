#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*----------------------------------------------------------------------*
 | read wall element                                         al 9/01    |
 *----------------------------------------------------------------------*/
void w1inp(ELEMENT *ele)
{
int  i;
int  ierr=0;
#ifdef DEBUG 
dstrc_enter("winp");
#endif
/*------------------------------------------------ allocate the element */      
ele->e.w1 = (WALL1*)CALLOC(1,sizeof(WALL1));
if (ele->e.w1==NULL) dserror("Allocation of element failed");
/*---------------------------------------------- read elements topology */
frchk("QUAD4",&ierr);
if (ierr==1) 
{
   ele->distyp = quad4;
   ele->numnp=4;
   ele->lm = (int*)CALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("QUAD4",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("QUAD8",&ierr);
if (ierr==1) 
{
   ele->distyp = quad8;
   ele->numnp=8;
   ele->lm = (int*)calloc(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("QUAD8",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("QUAD9",&ierr);
if (ierr==1) 
{
   ele->distyp = quad9;
   ele->numnp=9;
   ele->lm = (int*)calloc(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("QUAD9",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
/*------------------------------------------ reduce node numbers by one */
for (i=0; i<ele->numnp; i++) (ele->lm[i])--;
/*-------------------------------------------- read the material number */
frint("MAT",&(ele->mat),&ierr);
if (ierr!=1) dserror("Reading of WALL1 element failed");
/*---------------------------------------------- read element thickness */
frdouble("THICK",&(ele->e.w1->thick),&ierr);
if (ierr!=1) dserror("Reading of WALL1 element failed");
/*-------------------------- read gaussian points for wall elements */
frint_n("GP",&(ele->e.w1->nGP[0]),2,&ierr);
if (ierr!=1) dserror("Reading of WALL1 element failed");
/*---------------------------------------------- read 2D problem type */
ele->e.w1->wtype = plane_stress; /* default */
frchk("PLANE_STRESS",&ierr);
if (ierr==1) ele->e.w1->wtype = plane_stress;
frchk("PLANE_STRAIN",&ierr);
if (ierr==1) ele->e.w1->wtype = plane_strain;
frchk("ROTATIONAL_SYMMETRY",&ierr);
if (ierr==1) ele->e.w1->wtype = rotat_symmet;
    
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of w1inp */

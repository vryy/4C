#include "../headers/standardtypes.h"
#include "brick1.h"
/*----------------------------------------------------------------------*
 | read shell8 element                                    m.gee 8/00    |
 *----------------------------------------------------------------------*/
void b1inp(ELEMENT *ele)
{
int  i;
int  ierr=0;
int  quad;
int  counter;
long int  topology[100];
char *colpointer;
char buffer[50];
#ifdef DEBUG 
dstrc_enter("b1inp");
#endif
/*------------------------------------------------ allocate the element */      
ele->e.b1 = (BRICK1*)calloc(1,sizeof(BRICK1));
if (ele->e.b1==NULL) dserror("Allocation of element failed");
/*---------------------------------------------- read elements topology */
frchk("HEX8",&ierr);
if (ierr==1) 
{
   ele->numnp=8;
   ele->lm = (int*)calloc(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("HEX8",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("HEX20",&ierr);
if (ierr==1) 
{
   ele->numnp=20;
   ele->lm = (int*)calloc(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("HEX20",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("HEX27",&ierr);
if (ierr==1) 
{
   ele->numnp=27;
   ele->lm = (int*)calloc(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("HEX27",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("TET4",&ierr);
if (ierr==1) 
{
   ele->numnp=4;
   ele->lm = (int*)calloc(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("TET4",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("TET10",&ierr);
if (ierr==1) 
{
   ele->numnp=10;
   ele->lm = (int*)calloc(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("TET10",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
/*------------------------------------------ reduce node numbers by one */
for (i=0; i<ele->numnp; i++) (ele->lm[i])--;
/*-------------------------------------------- read the material number */
frint("MAT",&(ele->mat),&ierr);
if (ierr!=1) dserror("Reading of BRICK1 element failed");
/*-------------------------------------------- read the gaussian points */
frint_n("GP",&(ele->e.b1->nGP[0]),3,&ierr);
if (ierr!=1) dserror("Reading of BRICK1 element failed");
/*-------------------------- read gaussian points for triangle elements */
frint("GP_TET",&(ele->e.b1->nGP_tet),&ierr);
if (ierr!=1) dserror("Reading of BRICK1 element failed");
    
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of b1inp */

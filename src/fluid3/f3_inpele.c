#include "../headers/standardtypes.h"
#include "fluid3.h"
/*----------------------------------------------------------------------*
 | read fluid3 element                                    m.gee 8/00    |
 *----------------------------------------------------------------------*/
void f3inp(ELEMENT *ele)
{
int  i;
int  ierr=0;
int  quad;
int  counter;
long int  topology[100];
char *colpointer;
char buffer[50];
#ifdef DEBUG 
dstrc_enter("f3inp");
#endif
/*------------------------------------------------ allocate the element */      
ele->e.f3 = (FLUID3*)CALLOC(1,sizeof(FLUID3));
if (ele->e.f3==NULL) dserror("Allocation of element FLUID3 failed");
/*---------------------------------------------- read the element nodes */
frchk("HEX8",&ierr);
if (ierr==1)
{
   ele->numnp=8;
   ele->lm = (int*)CALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("HEX8",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("HEX20",&ierr);
if (ierr==1)
{
   ele->numnp=20;
   ele->lm = (int*)CALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("HEX20",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("TET4",&ierr);
if (ierr==1)
{
   ele->numnp=4;
   ele->lm = (int*)CALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("TET4",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("TET10",&ierr);
if (ierr==1)
{
   ele->numnp=10;
   ele->lm = (int*)CALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("TET10",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
/*------------------------------------------ reduce node numbers by one */
for (i=0; i<ele->numnp; i++) (ele->lm[i])--;
/*-------------------------------------------- read the material number */
frint("MAT",&(ele->mat),&ierr);
if (ierr!=1) dserror("Reading of FLUID3 element failed");
/*-------------------------------------------- read the gaussian points */
frint_n("GP",&(ele->e.f3->nGP[0]),2,&ierr);
if (ierr!=1) dserror("Reading of FLUID3 element failed");
/*-------------------------- read gaussian points for triangle elements */
frint("GP_TRI",&(ele->e.f3->nGP_tri),&ierr);
if (ierr!=1) dserror("Reading of FLUID3 element failed");
/*-------------------------------------------------- read is ale or not */
frint("ALE",&(ele->e.f3->is_ale),&ierr);
if (ierr!=1) dserror("Reading of FLUID3 element failed");
/*----------------------------- her now read all the other stuff needed */
/*----------------------------- her now read all the other stuff needed */
/*----------------------------- her now read all the other stuff needed */
/*----------------------------- her now read all the other stuff needed */
/*----------------------------- her now read all the other stuff needed */
    
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3inp */

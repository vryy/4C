#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "../fluid3/fluid3.h"
#include "ale2.h"
/*----------------------------------------------------------------------*
 | read ale element                                       m.gee 8/00    |
 *----------------------------------------------------------------------*/
void ale2inp(ELEMENT *ele)
{
int  i;
int  ierr=0;
int  quad;
int  counter;
long int  topology[100];
char *colpointer;
char buffer[50];
#ifdef DEBUG 
dstrc_enter("ale2inp");
#endif
/*------------------------------------------------ allocate the element */      
ele->e.ale2 = (ALE2*)CALLOC(1,sizeof(ALE2));
if (ele->e.ale2==NULL) dserror("Allocation of element ALE failed");
/*----------------------------------- read stuff needed for ALE element */
/*---------------------------------------------- read the element nodes */
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
   ele->lm = (int*)CALLOC(ele->numnp,sizeof(int));
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
if (ierr!=1) dserror("Reading of ALE2 element failed");
/*-------------------------------------------- read the gaussian points */
frint_n("GP",&(ele->e.ale2->nGP[0]),2,&ierr);
if (ierr!=1) dserror("Reading of ALE2 element failed");
/*-------------------------- read gaussian points for triangle elements */
frint("JAC",&(ele->e.ale2->jacobi),&ierr);
if (ierr!=1) dserror("Reading of ALE element failed");

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ale2inp */
#endif

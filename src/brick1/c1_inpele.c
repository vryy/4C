/*!----------------------------------------------------------------------
\file
\brief contains the routine 'c1inp' which reads data for a 3D hex element

*----------------------------------------------------------------------*/
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

/*! 
\addtogroup BRICK1 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief input of data for a 3D hex element

<pre>                                                              al 06/02
This routine reads data for an 3D-hex-element.

</pre>
\param *ele     ELEMENT    (o)   the element

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: inp_struct_field()

*----------------------------------------------------------------------*/
void c1inp(ELEMENT *ele)
{
int  i;
int  ierr=0;
int  quad;
int  counter;
int  eletopo[100];
long int  topology[100];
char *colpointer;
char buffer[50];
#ifdef DEBUG 
dstrc_enter("c1inp");
#endif
/*------------------------------------------------ allocate the element */      
ele->e.c1 = (BRICK1*)CCACALLOC(1,sizeof(BRICK1));
if (ele->e.c1==NULL) dserror("Allocation of element failed");
/*---------------------------------------------- read elements topology */
frchk("HEX8",&ierr);
if (ierr==1) 
{
   ele->distyp = hex8;
   ele->numnp=8;
   ele->lm = (int*)CCACALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("HEX8",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("HEX20",&ierr);
if (ierr==1) 
{
   ele->distyp = hex20;
   ele->numnp=20;
   ele->lm = (int*)CCACALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("HEX20",&(ele->lm[0]),ele->numnp,&ierr);
   /*
   ele->lm[ 0] = eletopo[ 0];
   ele->lm[ 1] = eletopo[ 1];
   ele->lm[ 2] = eletopo[ 2];
   ele->lm[ 3] = eletopo[ 3];
   ele->lm[ 4] = eletopo[ 4];
   ele->lm[ 5] = eletopo[ 5];
   ele->lm[ 6] = eletopo[ 6];
   ele->lm[ 7] = eletopo[ 7];
   ele->lm[ 8] = eletopo[ 8];
   ele->lm[ 9] = eletopo[ 9];
   ele->lm[10] = eletopo[10];
   ele->lm[11] = eletopo[11];
   ele->lm[12] = eletopo[16];
   ele->lm[13] = eletopo[17];
   ele->lm[14] = eletopo[18];
   ele->lm[15] = eletopo[19];
   ele->lm[16] = eletopo[12];
   ele->lm[17] = eletopo[13];
   ele->lm[18] = eletopo[14];
   ele->lm[19] = eletopo[15];
   */
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("HEX27",&ierr);
if (ierr==1) 
{
   ele->numnp=27;
   ele->distyp = hex27;
   ele->lm = (int*)CCACALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("HEX27",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("TET4",&ierr);
if (ierr==1) 
{
   ele->numnp=4;
   ele->lm = (int*)CCACALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("TET4",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("TET10",&ierr);
if (ierr==1) 
{
   ele->numnp=10;
   ele->lm = (int*)CCACALLOC(ele->numnp,sizeof(int));
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
frint_n("GP",&(ele->e.c1->nGP[0]),3,&ierr);
if (ierr!=1) dserror("Reading of BRICK1 element failed");
/*-------------------------- read eas-information */
frint("HYB",&(ele->e.c1->nhyb),&ierr);
/*------------------------ read element formulation: linear, T.L., U.L. */
frint("FORM",&(ele->e.c1->form),&ierr);
    
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of c1inp */
#endif
/*! @} (documentation module close)*/


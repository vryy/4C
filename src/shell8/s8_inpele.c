#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 | read shell8 element                                    m.gee 8/00    |
 *----------------------------------------------------------------------*/
void s8inp(ELEMENT *ele)
{
int  i;
int  ierr=0;
int  quad;
int  counter;
long int  topology[100];
char *colpointer;
char buffer[50];
int  nhyb=0;
#ifdef DEBUG 
dstrc_enter("s8inp");
#endif
/*------------------------------------------------ allocate the element */      
ele->e.s8 = (SHELL8*)CALLOC(1,sizeof(SHELL8));
if (ele->e.s8==NULL) dserror("Allocation of element failed");
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
   ele->lm = (int*)CALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("QUAD9",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("TRI3",&ierr);
if (ierr==1) 
{
   ele->distyp = tri3;
   ele->numnp=3;
   ele->lm = (int*)CALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("TRI3",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("TRI6",&ierr);
if (ierr==1) 
{
   ele->distyp = tri6;
   ele->numnp=6;
   ele->lm = (int*)CALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("TRI6",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
/*---------------------------------- allocate array for internal forces */
amdef("intforce",&(ele->e.s8->intforce),NUMDOF_SHELL8*ele->numnp,1,"DV");
/*------------------------------------------ reduce node numbers by one */
for (i=0; i<ele->numnp; i++) (ele->lm[i])--;
/*-------------------------------------------- read the material number */
frint("MAT",&(ele->mat),&ierr);
if (ierr!=1) dserror("Reading of SHELL8 element failed");
/*-------------------------------------------- read the shell thickness */
frdouble("THICK",&(ele->e.s8->thick),&ierr);
if (ierr!=1) dserror("Reading of SHELL8 element failed");
/*-------------------------------------------- read the gaussian points */
frint_n("GP",&(ele->e.s8->nGP[0]),3,&ierr);
if (ierr!=1) dserror("Reading of SHELL8 element failed");
/*-------------------------- read gaussian points for triangle elements */
frint("GP_TRI",&(ele->e.s8->nGP_tri),&ierr);
if (ierr!=1) dserror("Reading of SHELL8 element failed");
/*--------------------------------------- read local or global stresses */
frchar("FORCES",buffer,&ierr);
if (ierr)
{
   if (strncmp(buffer,"XYZ",3)==0)       ele->e.s8->forcetyp = s8_xyz;
   if (strncmp(buffer,"RST",3)==0)       ele->e.s8->forcetyp = s8_rst;
   if (strncmp(buffer,"RST_ortho",9)==0) ele->e.s8->forcetyp = s8_rst_ortho;
}
/*------------------------------------------------------------ read EAS */    
colpointer = strstr(allfiles.actplace,"EAS");
colpointer+=3;

colpointer = strpbrk(colpointer,"Nn");
ierr = sscanf(colpointer," %s ",buffer);    
if (ierr!=1) dserror("Reading of shell8 eas failed");
if (strncmp(buffer,"none",4)==0)  ele->e.s8->eas[0]=0;
if (strncmp(buffer,"N4_1",4)==0)  ele->e.s8->eas[0]=1;
if (strncmp(buffer,"N4_2",4)==0)  ele->e.s8->eas[0]=2;
if (strncmp(buffer,"N4_3",4)==0)  ele->e.s8->eas[0]=3;
if (strncmp(buffer,"N4_4",4)==0)  ele->e.s8->eas[0]=4;
if (strncmp(buffer,"N4_5",4)==0)  ele->e.s8->eas[0]=5;
if (strncmp(buffer,"N4_7",4)==0)  ele->e.s8->eas[0]=7;
if (strncmp(buffer,"N9_7",4)==0)  ele->e.s8->eas[0]=7;
if (strncmp(buffer,"N9_9",4)==0)  ele->e.s8->eas[0]=9;
if (strncmp(buffer,"N9_11",4)==0) ele->e.s8->eas[0]=11;
colpointer += strlen(buffer);    
    
colpointer = strpbrk(colpointer,"Nn");
ierr = sscanf(colpointer," %s ",buffer);    
if (ierr!=1) dserror("Reading of shell8 eas failed");
if (strncmp(buffer,"none",4)==0)  ele->e.s8->eas[1]=0;
if (strncmp(buffer,"N4_4",4)==0)  ele->e.s8->eas[1]=4;
if (strncmp(buffer,"N4_5",4)==0)  ele->e.s8->eas[1]=5;
if (strncmp(buffer,"N4_6",4)==0)  ele->e.s8->eas[1]=6;
if (strncmp(buffer,"N4_7",4)==0)  ele->e.s8->eas[1]=7;
if (strncmp(buffer,"N9_9",4)==0)  ele->e.s8->eas[1]=9;
if (strncmp(buffer,"N9_11",4)==0) ele->e.s8->eas[1]=11;
colpointer += strlen(buffer);    
    
colpointer = strpbrk(colpointer,"Nn");
ierr = sscanf(colpointer," %s ",buffer);    
if (ierr!=1) dserror("Reading of shell8 eas failed");
if (strncmp(buffer,"none",4)==0)  ele->e.s8->eas[2]=0;
if (strncmp(buffer,"N_1",4)==0)   ele->e.s8->eas[2]=1;
if (strncmp(buffer,"N_3",4)==0)   ele->e.s8->eas[2]=3;
if (strncmp(buffer,"N_4",4)==0)   ele->e.s8->eas[2]=4;
if (strncmp(buffer,"N_6",4)==0)   ele->e.s8->eas[2]=6;
if (strncmp(buffer,"N_8",4)==0)   ele->e.s8->eas[2]=8;
if (strncmp(buffer,"N_9",4)==0)   ele->e.s8->eas[2]=9;
colpointer += strlen(buffer);    

colpointer = strpbrk(colpointer,"Nn");
ierr = sscanf(colpointer," %s ",buffer);    
if (ierr!=1) dserror("Reading of shell8 eas failed");
if (strncmp(buffer,"none",4)==0)  ele->e.s8->eas[3]=0;
if (strncmp(buffer,"N4_2",4)==0)  ele->e.s8->eas[3]=2;
if (strncmp(buffer,"N4_4",4)==0)  ele->e.s8->eas[3]=4;
if (strncmp(buffer,"N9_2",4)==0)  ele->e.s8->eas[3]=2;
if (strncmp(buffer,"N9_4",4)==0)  ele->e.s8->eas[3]=4;
if (strncmp(buffer,"N9_6",4)==0)  ele->e.s8->eas[3]=6;
colpointer += strlen(buffer);    

colpointer = strpbrk(colpointer,"Nn");
ierr = sscanf(colpointer," %s ",buffer);    
if (ierr!=1) dserror("Reading of shell8 eas failed");
if (strncmp(buffer,"none",4)==0)  ele->e.s8->eas[4]=0;
if (strncmp(buffer,"N4_2",4)==0)  ele->e.s8->eas[4]=2;
if (strncmp(buffer,"N4_4",4)==0)  ele->e.s8->eas[4]=4;
if (strncmp(buffer,"N9_2",4)==0)  ele->e.s8->eas[4]=2;
if (strncmp(buffer,"N9_4",4)==0)  ele->e.s8->eas[4]=4;
if (strncmp(buffer,"N9_6",4)==0)  ele->e.s8->eas[4]=6;
/*--------------------- count nhyb and allocate storage for eas strains */
for (i=0; i<5; i++) nhyb+=ele->e.s8->eas[i];
ele->e.s8->nhyb=nhyb;
if (nhyb>0)
{
   amdef("alfa",&(ele->e.s8->alfa),1,nhyb,"DA");
   amzero(&(ele->e.s8->alfa));

   amdef("Dtildinv",&(ele->e.s8->Dtildinv),nhyb,nhyb,"DA");
   amzero(&(ele->e.s8->Dtildinv));
   
   amdef("Lt",&(ele->e.s8->Lt),nhyb,ele->numnp*NUMDOF_SHELL8,"DA");
   amzero(&(ele->e.s8->Lt));
   
   amdef("Rtilde",&(ele->e.s8->Rtilde),nhyb,1,"DV");
   amzero(&(ele->e.s8->Rtilde));
}
/*------------------------------------------------------------ read ANS */
frchar("ANS",buffer,&ierr);
if (ierr!=1) dserror("reading of shell8 ans failed");
if (strncmp(buffer,"none",4)==0)  ele->e.s8->ans=0;
if (strncmp(buffer,"Q",4)==0)     ele->e.s8->ans=1;
if (strncmp(buffer,"T",4)==0)     ele->e.s8->ans=2;
if (strncmp(buffer,"QT",4)==0)    ele->e.s8->ans=3;
if (strncmp(buffer,"TQ",4)==0)    ele->e.s8->ans=3;
/*------------------------------------------------------------ read sdc */
frdouble("SDC",&(ele->e.s8->sdc),&ierr);
if (ierr!=1) dserror("Reading of shell8 sdc failed");
    
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of s8inp */
#endif

/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9inp: which reads the input data of a shell9 element. 
         NOTE: The material has to be read before this element, 
               as the section data is defined within the material 
               (number of kinematic and material layers,...)

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/

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
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


/*!----------------------------------------------------------------------
\brief read shell9 element                                       

<pre>                     m.gee 6/01              modified by    sh 10/02
This routine reads all input data for a shell9 element
</pre>
\param  ELEMENT *ele   (i/o) element arrays which have to be filled

\warning The material has to be read before this element, as the section
         data is defined within the material (num_klay, ...)
\return void                                               
\sa calling: ---; called by: inp_struct_field()   [input_mesh.c]

*----------------------------------------------------------------------*/
void s9inp(ELEMENT *ele)
{
int          i;
int          ierr=0;
int          quad;
int          counter;
long int     topology[100];
char        *colpointer;
char         buffer[50];
int          nhyb=0;
int          numdof_shell9;
int          num_klay;
KINLAY      *actkinlay;
int          dummy_int;
MATERIAL    *actmat;

#ifdef DEBUG 
dstrc_enter("s9inp");
#endif
/*------------------------------------------------ allocate the element */      
ele->e.s9 = (SHELL9*)CCACALLOC(1,sizeof(SHELL9));
if (ele->e.s9==NULL) dserror("Allocation of element failed");
/*---------------------------------------------- read elements topology */
frchk("QUAD4",&ierr);
if (ierr==1) 
{
   ele->distyp = quad4;
   ele->numnp=4;
   ele->lm = (int*)CCACALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("QUAD4",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("QUAD8",&ierr);
if (ierr==1) 
{
   ele->distyp = quad8;
   ele->numnp=8;
   ele->lm = (int*)CCACALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("QUAD8",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("QUAD9",&ierr);
if (ierr==1) 
{
   ele->distyp = quad9;
   ele->numnp=9;
   ele->lm = (int*)CCACALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("QUAD9",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("TRI3",&ierr);
if (ierr==1) 
{
   ele->distyp = tri3;
   ele->numnp=3;
   ele->lm = (int*)CCACALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("TRI3",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("TRI6",&ierr);
if (ierr==1) 
{
   ele->distyp = tri6;
   ele->numnp=6;
   ele->lm = (int*)CCACALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("TRI6",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
/*-------------------------------------------- read the material number */
frint("MAT",&(ele->mat),&ierr);
if (ierr!=1) dserror("Reading of SHELL9 element failed");
/*----------------------------------- MultilayerMaterial to this element*/
actmat = &(mat[ele->mat-1]);
if (actmat->mattyp != m_multi_layer) dserror("Wrong mattyp to SHELL9 element -> has to be MultiLayer");

/*-------- read cross sectional data from MAT_Multilayer to this element*/
  /*------------------------ number of kinematic layers to this element */
  num_klay = actmat->m.multi_layer->num_klay;
  ele->e.s9->num_klay = num_klay;
  /*---------------- calculate numdof of actual element an set it to it */
  numdof_shell9 = 3 + 3*ele->e.s9->num_klay;
  ele->e.s9->numdf = numdof_shell9;
  /*-------------- allocate vector for layerhgt of kinematic layers ------*/
  ele->e.s9->klayhgt = (double*)CCACALLOC(ele->e.s9->num_klay,sizeof(double));
  if (ele->e.s9->klayhgt==NULL) dserror("Allocation of SHELL9 element failed -> klayhgt");
  ele->e.s9->klayhgt = actmat->m.multi_layer->klayhgt;
  /*--------------------- allocate a kinlay for each kinematic layer ------*/
  ele->e.s9->kinlay = (KINLAY*)CCACALLOC(ele->e.s9->num_klay,sizeof(KINLAY));
  if (ele->e.s9->kinlay==NULL) dserror("Allocation of SHELL9 element failed -> kinlay");
  /*-------------------- put information for each kinematic layer --------*/
  for (i=0; i<ele->e.s9->num_klay; i++)
  {
     ele->e.s9->kinlay[i] = actmat->m.multi_layer->kinlay[i];
  }
 
/*---------------------------------- allocate array for internal forces */
amdef("intforce",&(ele->e.s9->intforce),numdof_shell9*ele->numnp,1,"DV");
/*------------------------------------------ reduce node numbers by one */
for (i=0; i<ele->numnp; i++) (ele->lm[i])--;
/*-------------------------------------------- read the shell thickness */
frdouble("THICK",&(ele->e.s9->thick),&ierr);
if (ierr!=1) dserror("Reading of SHELL9 element failed");
/*-------------------------------------------- read the gaussian points */
frint_n("GP",&(ele->e.s9->nGP[0]),3,&ierr);
if (ierr!=1) dserror("Reading of SHELL9 element failed");
if (ele->e.s9->nGP[2] != 2) dserror("nGP[2] != 2 -> only 2 GP in Thickness direction implemented: Reading of SHELL9 element failed");
/*write a warning if a 8/9-noded Element is calculated with less than 3 GPs*/
if (ele->distyp == quad9 && ele->e.s9->nGP[0] < 3) printf("WARNING in s9_inpele.c: QUAD9 but GP[0]<3 ; this could lead to ZEMs \n");
if (ele->distyp == quad8 && ele->e.s9->nGP[0] < 3) printf("WARNING in s9_inpele.c: QUAD8 but GP[0]<3 ; this could lead to ZEMs \n");
/*-------------------------- read gaussian points for triangle elements */
frint("GP_TRI",&(ele->e.s9->nGP_tri),&ierr);
if (ierr!=1) dserror("Reading of SHELL9 element failed");
/*--------------------------------------- read local or global stresses */
frchar("FORCES",buffer,&ierr);
if (ierr)
{
   if (strncmp(buffer,"XYZ",3)==0)       ele->e.s9->forcetyp = s9_xyz;
   if (strncmp(buffer,"RST",3)==0)       ele->e.s9->forcetyp = s9_rst;
   if (strncmp(buffer,"RST_ortho",9)==0) ele->e.s9->forcetyp = s9_rst_ortho;
}
/*------------------------------------------------------------ read EAS */    
colpointer = strstr(allfiles.actplace,"EAS");
colpointer+=3;

colpointer = strpbrk(colpointer,"Nn");
ierr = sscanf(colpointer," %s ",buffer);    
if (ierr!=1) dserror("Reading of shell9 eas failed");
if (strncmp(buffer,"none",4)==0)  ele->e.s9->eas[0]=0;
if (strncmp(buffer,"N4_1",4)==0)  ele->e.s9->eas[0]=1;
if (strncmp(buffer,"N4_2",4)==0)  ele->e.s9->eas[0]=2;
if (strncmp(buffer,"N4_3",4)==0)  ele->e.s9->eas[0]=3;
if (strncmp(buffer,"N4_4",4)==0)  ele->e.s9->eas[0]=4;
if (strncmp(buffer,"N4_5",4)==0)  ele->e.s9->eas[0]=5;
if (strncmp(buffer,"N4_7",4)==0)  ele->e.s9->eas[0]=7;
if (strncmp(buffer,"N9_7",4)==0)  ele->e.s9->eas[0]=7;
if (strncmp(buffer,"N9_9",4)==0)  ele->e.s9->eas[0]=9;
if (strncmp(buffer,"N9_11",4)==0) ele->e.s9->eas[0]=11;
colpointer += strlen(buffer);    
    
colpointer = strpbrk(colpointer,"Nn");
ierr = sscanf(colpointer," %s ",buffer);    
if (ierr!=1) dserror("Reading of shell9 eas failed");
if (strncmp(buffer,"none",4)==0)  ele->e.s9->eas[1]=0;
if (strncmp(buffer,"N4_4",4)==0)  ele->e.s9->eas[1]=4;
if (strncmp(buffer,"N4_5",4)==0)  ele->e.s9->eas[1]=5;
if (strncmp(buffer,"N4_6",4)==0)  ele->e.s9->eas[1]=6;
if (strncmp(buffer,"N4_7",4)==0)  ele->e.s9->eas[1]=7;
if (strncmp(buffer,"N9_9",4)==0)  ele->e.s9->eas[1]=9;
if (strncmp(buffer,"N9_11",4)==0) ele->e.s9->eas[1]=11;
colpointer += strlen(buffer);    
    
colpointer = strpbrk(colpointer,"Nn");
ierr = sscanf(colpointer," %s ",buffer);    
if (ierr!=1) dserror("Reading of shell9 eas failed");
if (strncmp(buffer,"none",4)==0)  ele->e.s9->eas[2]=0;
if (strncmp(buffer,"N_1",4)==0)   ele->e.s9->eas[2]=1;
if (strncmp(buffer,"N_3",4)==0)   ele->e.s9->eas[2]=3;
if (strncmp(buffer,"N_4",4)==0)   ele->e.s9->eas[2]=4;
if (strncmp(buffer,"N_6",4)==0)   ele->e.s9->eas[2]=6;
if (strncmp(buffer,"N_8",4)==0)   ele->e.s9->eas[2]=8;
if (strncmp(buffer,"N_9",4)==0)   ele->e.s9->eas[2]=9;
colpointer += strlen(buffer);    

colpointer = strpbrk(colpointer,"Nn");
ierr = sscanf(colpointer," %s ",buffer);    
if (ierr!=1) dserror("Reading of shell9 eas failed");
if (strncmp(buffer,"none",4)==0)  ele->e.s9->eas[3]=0;
if (strncmp(buffer,"N4_2",4)==0)  ele->e.s9->eas[3]=2;
if (strncmp(buffer,"N4_4",4)==0)  ele->e.s9->eas[3]=4;
if (strncmp(buffer,"N9_2",4)==0)  ele->e.s9->eas[3]=2;
if (strncmp(buffer,"N9_4",4)==0)  ele->e.s9->eas[3]=4;
if (strncmp(buffer,"N9_6",4)==0)  ele->e.s9->eas[3]=6;
colpointer += strlen(buffer);    

colpointer = strpbrk(colpointer,"Nn");
ierr = sscanf(colpointer," %s ",buffer);    
if (ierr!=1) dserror("Reading of shell9 eas failed");
if (strncmp(buffer,"none",4)==0)  ele->e.s9->eas[4]=0;
if (strncmp(buffer,"N4_2",4)==0)  ele->e.s9->eas[4]=2;
if (strncmp(buffer,"N4_4",4)==0)  ele->e.s9->eas[4]=4;
if (strncmp(buffer,"N9_2",4)==0)  ele->e.s9->eas[4]=2;
if (strncmp(buffer,"N9_4",4)==0)  ele->e.s9->eas[4]=4;
if (strncmp(buffer,"N9_6",4)==0)  ele->e.s9->eas[4]=6;
/*--------------------- count nhyb and allocate storage for eas strains */
for (i=0; i<5; i++) nhyb+=ele->e.s9->eas[i];
ele->e.s9->nhyb=nhyb;
if (nhyb>0)
{
   amdef("alfa",&(ele->e.s9->alfa),num_klay,nhyb,"DA");
   amzero(&(ele->e.s9->alfa));

   for (i=0; i<num_klay; i++) amdef("Dtildinv",&(ele->e.s9->Dtildinv[i]),nhyb,nhyb,"DA");
   for (i=0; i<num_klay; i++) amzero(&(ele->e.s9->Dtildinv[i]));

   for (i=0; i<num_klay; i++) amdef("Lt",&(ele->e.s9->Lt[i]),nhyb,ele->numnp*numdof_shell9,"DA");
   for (i=0; i<num_klay; i++) amzero(&(ele->e.s9->Lt[i]));
   
   for (i=0; i<num_klay; i++) amdef("Rtilde",&(ele->e.s9->Rtilde[i]),nhyb,1,"DV");
   for (i=0; i<num_klay; i++) amzero(&(ele->e.s9->Rtilde[i]));
}
/*------------------------------------------------------------ read ANS */
frchar("ANS",buffer,&ierr);
if (ierr!=1) dserror("reading of shell9 ans failed");
if (strncmp(buffer,"none",4)==0)  ele->e.s9->ans=0;
if (strncmp(buffer,"Q",4)==0)     ele->e.s9->ans=1;
/*if (strncmp(buffer,"T",4)==0)     ele->e.s9->ans=2;
if (strncmp(buffer,"QT",4)==0)    ele->e.s9->ans=3;
if (strncmp(buffer,"TQ",4)==0)    ele->e.s9->ans=3;
/************************************************************************/
/************************************************************************/
/*------------------------------------------------------------ read sdc */
frdouble("SDC",&(ele->e.s9->sdc),&ierr);
if (ierr!=1) dserror("Reading of shell9 sdc failed");

/*----------------------------------------------------------------------*/    
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of s9inp */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/

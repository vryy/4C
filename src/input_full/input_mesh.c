#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 | Global variables for this file                        m.gee 11/00    |
 *----------------------------------------------------------------------*/
static ARRAY tmpnodes1;
static ARRAY tmpnodes2;



/*----------------------------------------------------------------------*
 | input of fields                                        m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inpfield()
{
int  i;
int  ierr;
#ifdef DEBUG 
dstrc_enter("inpfield");
#endif
/*--------------------------------------read node coordinates from file */
inpnodes();
/*--------------------------------------------------------- read field  */
/*-------- dependent on the problem typ need different number of fields */
/*----------------------------------------------- FSI 3D typ of problem */
if (genprob.probtyp == prb_fsi)
{
   if (genprob.numfld!=3) dserror("numfld != 3 for FSI");
   
   field = (FIELD*)calloc(genprob.numfld,sizeof(FIELD));
   if (field==NULL) dserror("Allocation of fields failed");

   field[0].fieldtyp = structure;
   inp_struct_field(&(field[0]));
   
   field[1].fieldtyp = fluid;
   inp_fluid_field (&(field[1]));
   
   field[2].fieldtyp = ale;
   inp_ale_field  (&(field[2]));

/*--- ale and fluid field are supposed to be compatible, so inherit info */
/*
   fluid_to_ale(&(field[1]),&(field[2]));
*/
}
/*------------------------------------------- structure type of problem */
if (genprob.probtyp==prb_structure)
{
   if (genprob.numfld!=1) dserror("numfld != 1 for structural problem");
   field = (FIELD*)calloc(genprob.numfld,sizeof(FIELD));
   if (field==NULL) dserror("Allocation of fields failed");

   field[0].fieldtyp = structure;
   inp_struct_field(&(field[0]));
}
/*----------------------------------------------- fluid type of problem */
if (genprob.probtyp==prb_fluid)
{
   if (genprob.numfld!=1) dserror("numfld != 1 for fluid problem");
   field = (FIELD*)calloc(genprob.numfld,sizeof(FIELD));
   if (field==NULL) dserror("Allocation of fields failed");
   
   field[0].fieldtyp = fluid;
   inp_fluid_field (&(field[0]));
}
/*---------------------------------------- Optimisation type of problem */
if (genprob.probtyp == prb_opt)
{
/* not yet implemented */
}
/*-------------------------------------- assign the nodes to the fields */
for (i=0; i<genprob.numfld; i++)
{
   inp_assign_nodes(&(field[i]));
}
amdel(&tmpnodes1);
amdel(&tmpnodes2);
/*---------------------------------- make element-node-element topology */
for (i=0; i<genprob.numfld; i++)
{
   inp_topology(&(field[i]));
}
/*----------------------------------------------- FSI 3D typ of problem */
if (genprob.probtyp==prb_fsi)
{
/*--- ale and fluid field are supposed to be compatible, so inherit info */
   fluid_to_ale(&(field[1]),&(field[2]));
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpfield */





/*----------------------------------------------------------------------*
 | sort nodes to the fields                               m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_assign_nodes(FIELD *field)
{
int  i,j,k;
int  ierr;
int  node_Id;
int counter;
int minusone=-1;
ARRAY nodeflag;
ARRAY coords;
ELEMENT *actele;
#ifdef DEBUG 
dstrc_enter("inp_assign_nodes");
#endif
amdef("nodeflag",&nodeflag,genprob.nnode,1,"IV");
aminit(&nodeflag,&minusone);
/*----------------  set a flag to the node_id for each node in the field */
for (i=0; i<field->numele; i++)
{
   actele = &(field->element[i]);
   for (j=0; j<actele->numnp; j++)
   {
      node_Id = actele->lm[j];
      nodeflag.a.iv[node_Id]=node_Id;
   }
}
/*----------------------------------------------------- count the flags */
counter=0;
for (i=0; i<genprob.nnode; i++)
{
   if (nodeflag.a.iv[i]!=-1) counter++;
}
field->numnp=counter;
/*-------------------------------------- Allocate the nodes to the field */
field->node = (NODE*)calloc(field->numnp,sizeof(NODE));
if (field->node==NULL) dserror("Allocation of nodes failed");
/*---------------- assign the node Ids and coords to the NODE structure */
counter=0;
for (i=0; i<genprob.nnode; i++)
{
   if (nodeflag.a.iv[i]!=-1)
   {
      field->node[counter].Id = nodeflag.a.iv[i];
      for (j=0; j<3; j++)
      {
         field->node[counter].x[j] = tmpnodes1.a.da[i][j];
      }
      counter++;
   }
}
/*----------------------------------------------------------------------*/
amdel(&nodeflag);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_assign_nodes */







/*----------------------------------------------------------------------*
 | input of node coords                                   m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inpnodes()
{
int  ierr=0;
int  counter;
#ifdef DEBUG 
dstrc_enter("inpnodes");
#endif
/*--------------- allocate temporary array for coordinates of all nodes */
amdef("tempnod1",&tmpnodes1,(genprob.nnode),3,"DA");
amdef("tempnod2",&tmpnodes2,(genprob.nnode),1,"IV");
/*--------------------------------------------------------- rewind file */
frrewind();
/*---------------------------------------------------------- read nodes */
frfind("--NODE COORDS");
frread();
counter=0;
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   frint("NODE",&(tmpnodes2.a.iv[counter]),&ierr);
   if (ierr!=1) dserror("reading of nodes failed");
   (tmpnodes2.a.iv[counter])--;

   frdouble_n("COORD",&(tmpnodes1.a.da[counter][0]),3,&ierr);
   if (ierr!=1) dserror("reading of nodes failed");

   counter++;

   frread();
}
frrewind();
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpnodes */





/*----------------------------------------------------------------------*
 | input of structure field                               m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_struct_field(FIELD *structfield)
{
int  ierr;
int  counter=0;
int  elenumber;
int  isquad;
char *colpointer;
#ifdef DEBUG 
dstrc_enter("inp_struct_field");
#endif
/*-------------------------------------------- count number of elements */
frrewind();
frfind("--STRUCTURE ELEMENTS");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   counter++;
   frread();
}
frrewind();
structfield->numele = counter;
/*--------------------------------------------------- allocate elements */
structfield->element=(ELEMENT*)calloc(structfield->numele,sizeof(ELEMENT));
if (structfield->element==NULL) dserror("Allocation of ELEMENT failed");
/*------------------------------------------------------- read elements */
frrewind();
frfind("--STRUCTURE ELEMENTS");
frread();
counter=0;
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   colpointer = allfiles.actplace;
   elenumber  = strtol(colpointer,&colpointer,10);
   structfield->element[counter].Id = --elenumber;
/*---------- read the typ of element and call element readning function */
/*------------------------------------------------ elementtyp is SHELL8 */
   frchk("SHELL8",&ierr);
   if (ierr==1)
   {
#ifndef D_SHELL8 
      dserror("SHELL8 needed but not defined in Makefile");
#endif
   }
#ifdef D_SHELL8 
   if (ierr==1) 
   {
      structfield->element[counter].eltyp=el_shell8;
      s8inp(&(structfield->element[counter]));
   }
#endif
/*------------------------------------------------ elementtyp is BRICK1 */
   frchk("BRICK1",&ierr);
   if (ierr==1)
   {
#ifndef D_BRICK1 
      dserror("BRICK1 needed but not defined in Makefile");
#endif
   }
#ifdef D_BRICK1 
   if (ierr==1) 
   {
      structfield->element[counter].eltyp=el_brick1;
      b1inp(&(structfield->element[counter]));
   }
#endif
/*------------------------------------------------ elementtyp is WALL  */
   frchk("WALL",&ierr);
   if (ierr==1)
   {
#ifndef D_WALL1 
      dserror("WALL1 needed but not defined in Makefile");
#endif
   }
#ifdef D_WALL1 
   if (ierr==1) 
   {
      structfield->element[counter].eltyp=el_wall1;
      w1inp(&(structfield->element[counter]));
   }
#endif
/*--------------------------------------------other structural elements */
   counter++;
   frread();
}
frrewind();
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_struct_field */








/*----------------------------------------------------------------------*
 | input of fluid field                                   m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_fluid_field(FIELD *fluidfield)
{
int  ierr;
int  counter=0;
int  elenumber;
int  isquad;
char *colpointer;
#ifdef DEBUG 
dstrc_enter("inp_fluid_field");
#endif
/*-------------------------------------------- count number of elements */
frrewind();
frfind("--FLUID ELEMENTS");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   counter++;
   frread();
}
frrewind();
fluidfield->numele = counter;
/*--------------------------------------------------- allocate elements */
fluidfield->element=(ELEMENT*)calloc(fluidfield->numele,sizeof(ELEMENT));
if (fluidfield->element==NULL) dserror("Allocation of ELEMENT failed");
/*------------------------------------------------------- read elements */
frrewind();
frfind("--FLUID ELEMENTS");
frread();
counter=0;
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   colpointer = allfiles.actplace;
   elenumber  = strtol(colpointer,&colpointer,10);
   fluidfield->element[counter].Id = --elenumber;
/*---------- read the typ of element and call element readning function */
/*------------------------------------------------ elementtyp is FLUID3 */
   frchk("FLUID3",&ierr);
   if (ierr==1)
   {
#ifndef D_FLUID3 
      dserror("FLUID3 needed but not defined in Makefile");
#endif
   }
#ifdef D_FLUID3 
   if (ierr==1) 
   {
      fluidfield->element[counter].eltyp=el_fluid3;
      f3inp(&(fluidfield->element[counter]));
   }
#endif
/*------------------------------------------------ elementtyp is FLUID1 */
   frchk("FLUID1",&ierr);
   if (ierr==1)
   {
#ifndef D_FLUID1 
      dserror("FLUID3 needed but not defined in Makefile");
#endif
   }
#ifdef D_FLUID1 
   if (ierr==1) 
   {
      fluidfield->element[counter].eltyp=el_fluid1;
      /* not implemented reading of one fluid1 line */
   }
#endif
   counter++;
   frread();
}
frrewind();
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_fluid_field */








/*----------------------------------------------------------------------*
 | input of ale field                                     m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_ale_field(FIELD *alefield)
{
int  ierr;
int  counter=0;
int  elenumber;
int  isquad;
char *colpointer;
#ifdef DEBUG 
dstrc_enter("inp_ale_field");
#endif

/*-------------------------------------------- count number of elements */
frrewind();
frfind("--ALE ELEMENTS");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   counter++;
   frread();
}
frrewind();
alefield->numele = counter;
/*--------------------------------------------------- allocate elements */
alefield->element=(ELEMENT*)calloc(alefield->numele,sizeof(ELEMENT));
if (alefield->element==NULL) dserror("Allocation of ELEMENT failed");
/*------------------------------------------------------- read elements */
frrewind();
frfind("--ALE ELEMENTS");
frread();
counter=0;
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   colpointer = allfiles.actplace;
   elenumber  = strtol(colpointer,&colpointer,10);
   alefield->element[counter].Id = --elenumber;
/*---------- read the typ of element and call element reading function */
/*------------------------------------------------ elementtyp is ALE */
   frchk("ALE3",&ierr);
   if (ierr==1)
   {
#ifndef D_ALE 
      dserror("ALE3 needed but not defined in Makefile");
#endif
   }
#ifdef D_ALE 
   if (ierr==1) 
   {
      alefield->element[counter].eltyp=el_ale;
      ale3inp(&(alefield->element[counter]));
   }
#endif
   counter++;
   frread();
}
frrewind();
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_ale_field */

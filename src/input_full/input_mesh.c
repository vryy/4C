/*!----------------------------------------------------------------------
\file
\brief input of mesh data

------------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../shell8/shell8.h"
#include "../shell9/shell9.h"
#include "../wall1/wall1.h"
#include "../brick1/brick1.h"
#include "../fluid3/fluid3.h"
#include "../fluid3/fluid3_prototypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid2/fluid2_prototypes.h"
#include "../ale3/ale3.h"
#include "../ale2/ale2.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
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
 | Global variables for this file                        m.gee 11/00    |
 *----------------------------------------------------------------------*/
static ARRAY tmpnodes1;
static ARRAY tmpnodes2;

/*----------------------------------------------------------------------*
 | input of fields                                        m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inpfield()
{
int  i,j;
int  ierr;
#ifdef DEBUG 
dstrc_enter("inpfield");
#endif

/*--------------------------------------read node coordinates from file */
inpnodes();
/*--------------------------------------------------------- read field  */
/*----------------------------------------------- FSI 3D typ of problem */
if (genprob.probtyp == prb_fsi)
{
   if (genprob.numfld!=3) dserror("numfld != 3 for FSI");
   
   field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));

   field[genprob.numsf].fieldtyp = structure;   
   if (genprob.multidis>0) inpdis(&(field[genprob.numsf]));
   else field[genprob.numsf].ndis=1;
   inp_struct_field(&(field[genprob.numsf]));
   
   field[genprob.numff].fieldtyp = fluid;
   if (genprob.multidis>0) inpdis(&(field[genprob.numff]));
   else field[genprob.numff].ndis=1;
   inp_fluid_field (&(field[genprob.numff]));
   
   field[genprob.numaf].fieldtyp = ale;
   if (genprob.multidis>0) inpdis(&(field[genprob.numaf]));
   else field[genprob.numaf].ndis=1;
   inp_ale_field  (&(field[genprob.numaf]));

/*--- ale and fluid field are supposed to be compatible, so inherit info */
#ifdef D_ALE
/*
   fluid_to_ale(&(field[1]),&(field[2]));
*/
#endif
}
/*------------------------------------------- structure type of problem */
if (genprob.probtyp==prb_structure)
{
   if (genprob.numfld!=1) dserror("numfld != 1 for structural problem");
   field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
   
   field[genprob.numsf].fieldtyp = structure;
   if (genprob.multidis>0) inpdis(&(field[genprob.numsf]));
   else field[genprob.numsf].ndis=1;
   inp_struct_field(&(field[genprob.numsf]));
}
/*---------------------------------------- Optimisation type of problem */
if (genprob.probtyp == prb_opt)
{  /*-- structure type of problem */
   if (genprob.numfld!=1) dserror("numfld != 1 for structural problem");
   field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));

   field[0].fieldtyp = structure;
   if (genprob.multidis>0) inpdis(&(field[0]));
   else field[0].ndis=1;
   inp_struct_field(&(field[0]));
}
/*----------------------------------------------- fluid type of problem */
if (genprob.probtyp==prb_fluid)
{
   if (genprob.numfld==1) /* single field fluid problem                 */
   {
      field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
   
      field[genprob.numff].fieldtyp = fluid;
      if (genprob.multidis>0) inpdis(&(field[genprob.numff]));
      else field[genprob.numff].ndis=1;
      inp_fluid_field (&(field[genprob.numff]));
   }
   else if (genprob.numfld==2) /* two field fluid problem (fluid+ale)       */
   {
      field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
      
      field[genprob.numff].fieldtyp = fluid;
      if (genprob.multidis>0) inpdis(&(field[genprob.numff]));
      else field[genprob.numff].ndis=1;
      inp_fluid_field (&(field[genprob.numff]));

      field[genprob.numaf].fieldtyp = ale;
      if (genprob.multidis>0) inpdis(&(field[genprob.numaf]));
      else field[genprob.numaf].ndis=1;
      inp_ale_field  (&(field[genprob.numaf]));      
   }
   else dserror("NUMFLD>2 not allowed for Problemtype FLUID\n");   
}
/*------------------------------------------------- ale type of problem */
if (genprob.probtyp==prb_ale)
{
   if (genprob.numfld!=1) dserror("numfld != 1 for ale problem");
   field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
   
   field[genprob.numaf].fieldtyp = ale;
   if (genprob.multidis>0) inpdis(&(field[genprob.numaf]));
   else field[genprob.numaf].ndis=1;
   inp_ale_field(&(field[genprob.numaf]));
}

/*---------------------------------------- Optimisation type of problem */
if (genprob.probtyp == prb_opt)
{
/* not yet implemented */
}
/*-------------------------------------- assign the nodes to the fields */
for (i=0; i<genprob.numfld; i++)
for (j=0;j<field[i].ndis;j++)
   inp_assign_nodes(&(field[i].dis[j]));

amdel(&tmpnodes1);
amdel(&tmpnodes2);
/*---------------------------------- make element-node-element topology */
for (i=0; i<genprob.numfld; i++)
for (j=0;j<field[i].ndis;j++)
   inp_topology(&(field[i].dis[j]));

/*----------------------------------------------- FSI 3D typ of problem */
if (genprob.probtyp==prb_fsi)
{
/*--- ale and fluid field are supposed to be compatible, so inherit info */
#ifdef D_ALE
   fluid_to_ale(&(field[1]),&(field[2]));
#endif
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
void inp_assign_nodes(DISCRET *actdis)
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
for (i=0; i<actdis->numele; i++)
{
   actele = &(actdis->element[i]);
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
actdis->numnp=counter;
/*-------------------------------------- Allocate the nodes to the field */
actdis->node = (NODE*)CCACALLOC(counter,sizeof(NODE));
/*---------------- assign the node Ids and coords to the NODE structure */
counter=0;
for (i=0; i<genprob.nnode; i++)
{
   if (nodeflag.a.iv[i]!=-1)
   {
      actdis->node[counter].Id = nodeflag.a.iv[i];
      for (j=0; j<3; j++)
      {
         actdis->node[counter].x[j] = tmpnodes1.a.da[i][j];
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



/*!---------------------------------------------------------------------                                         
\brief input of disretisation data

<pre>                                                         genk 08/02		     
</pre>

\return void                                                                             

------------------------------------------------------------------------*/
void inpdis(FIELD *actfield)
{
int  ierr=0;
char buffer[50];
#ifdef DEBUG 
dstrc_enter("inpdis");
#endif

/*------------------------------------------------- read discretisation */
if (frfind("--DISCRETISATION")==0) goto end;
frread();
switch(actfield->fieldtyp)
{
case fluid: 
   while(strncmp(allfiles.actplace,"------",6)!=0)
   {
      frint("NUMFLUIDDIS", &(actfield->ndis),&ierr);
      frread();
   }
break;
case structure:
   while(strncmp(allfiles.actplace,"------",6)!=0)
   {
      frint("NUMSTRUCDIS", &(actfield->ndis),&ierr);
      frread();
   }
break;
case ale:
   while(strncmp(allfiles.actplace,"------",6)!=0)
   {
      frint("NUMALEDIS", &(actfield->ndis),&ierr);
      frread();
   }
break;
}   
frrewind();
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdis */



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
/*---------------------------------------------------------- read nodes */
if (frfind("--NODE COORDS")==0) dserror("frfind: NODE COORDS is not in input file");
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

end:
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
/*----------------------------------------- allocate one discretization */
/*structfield->ndis=1;*/
if (structfield->ndis>1)
   dserror("different discretisations not implemented yet for structural elements\n");
structfield->dis = (DISCRET*)CCACALLOC(structfield->ndis,sizeof(DISCRET));
/*-------------------------------------------- count number of elements */
if (frfind("--STRUCTURE ELEMENTS")==1)
{
  frread();
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    counter++;
    frread();
  }
}
structfield->dis[0].numele = counter;
/*--------------------------------------------------- allocate elements */
structfield->dis[0].element=(ELEMENT*)CCACALLOC(structfield->dis[0].numele,sizeof(ELEMENT));
/*------------------------------------------------------- read elements */
if (frfind("--STRUCTURE ELEMENTS")==0) goto end;
frread();
counter=0;
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   colpointer = allfiles.actplace;
   elenumber  = strtol(colpointer,&colpointer,10);
   structfield->dis[0].element[counter].Id = --elenumber;
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
      structfield->dis[0].element[counter].eltyp=el_shell8;
      s8inp(&(structfield->dis[0].element[counter]));
   }
#endif
/*------------------------------------------------ elementtyp is SHELL9 */
   frchk("SHELL9",&ierr);
   if (ierr==1)
   {
#ifndef D_SHELL9 
      dserror("SHELL9 needed but not defined in Makefile");
#endif
   }
#ifdef D_SHELL9 
   if (ierr==1) 
   {
      structfield->dis[0].element[counter].eltyp=el_shell9;
      s9inp(&(structfield->dis[0].element[counter]));
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
      structfield->dis[0].element[counter].eltyp=el_brick1;
      c1inp(&(structfield->dis[0].element[counter]));
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
      structfield->dis[0].element[counter].eltyp=el_wall1;
      w1inp(&(structfield->dis[0].element[counter]));
   }
#endif
/*--------------------------------------------other structural elements */
   counter++;
   frread();
}
frrewind();
/*----------------------------------------------------------------------*/

end:
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
int  cpro=0;
int  elenumber;
int  isquad;
char *colpointer;
#ifdef DEBUG 
dstrc_enter("inp_fluid_field");
#endif
/*-------------------------------------------- allocate discretizations */
/*fluidfield->ndis=1; */
fluidfield->dis = (DISCRET*)CCACALLOC(fluidfield->ndis,sizeof(DISCRET));

/*
remarks about different discretisations:
we asume to read in one "global" discretisation from the input file.
from this discretisation all the other ones can be directly derived!!!  */
/*-------------------------------------------- count number of elements */
if (frfind("--FLUID ELEMENTS")==1)
{
  frread();
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    counter++;
    frread();
  }
}
fluidfield->dis[0].numele = counter;
/*--------------------------------------------------- allocate elements */
fluidfield->dis[0].element=(ELEMENT*)CCACALLOC(fluidfield->dis[0].numele,sizeof(ELEMENT));
/*------------------------------------------------------- read elements */
if (frfind("--FLUID ELEMENTS")==0) goto end;
frread();
counter=0;
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   colpointer = allfiles.actplace;
   elenumber  = strtol(colpointer,&colpointer,10);
   fluidfield->dis[0].element[counter].Id = --elenumber;
/*---------- read the typ of element and call element readning function */
/*-------------------------------------------- elementtyp is FLUID2_PRO */
   frchk("FLUID2_PRO",&ierr);
   if (ierr==1)
   {
#ifndef D_FLUID2_PRO 
      dserror("FLUID2_PRO needed but not defined in Makefile");
#endif
   }
#ifdef D_FLUID2_PRO 
   if (ierr==1) 
   {
      /*-------------------------- allocate elements of second discretisation */
      if (cpro==0)
      {
         if(fluidfield->ndis<2) 
            dserror("NUMFLUIDDIS has to be g.t. 1 for FLUID2_PRO Elements!\n");
         fluidfield->dis[1].numele = fluidfield->dis[0].numele;
	 fluidfield->dis[1].element=(ELEMENT*)CCACALLOC(fluidfield->dis[1].numele,sizeof(ELEMENT));
         cpro++;
      } /* endif (cpro==0) */      
      fluidfield->dis[0].element[counter].eltyp=el_fluid2_pro;
      fluidfield->dis[1].element[counter].eltyp=el_fluid2_pro;
      f2pro_inp(&(fluidfield->dis[0].element[counter]));
      f2pro_dis(&(fluidfield->dis[0].element[counter]),&(fluidfield->dis[1].element[counter]));       
      goto read;
   }
#endif
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
      fluidfield->dis[0].element[counter].eltyp=el_fluid3;
      f3inp(&(fluidfield->dis[0].element[counter]),counter);
   }
#endif   
/*------------------------------------------------ elementtyp is FLUID2 */
   frchk("FLUID2",&ierr);
   if (ierr==1)
   {
#ifndef D_FLUID2 
      dserror("FLUID2 needed but not defined in Makefile");
#endif
   }
#ifdef D_FLUID2 
   if (ierr==1) 
   {
      fluidfield->dis[0].element[counter].eltyp=el_fluid2;
      f2_inp(&(fluidfield->dis[0].element[counter]),counter);

      if (fluidfield->dis[0].element[counter].e.f2->turbu == 2 || 
          fluidfield->dis[0].element[counter].e.f2->turbu == 3 )
      {
       if(cpro==0)
       {
        fluidfield->dis[1].numele = fluidfield->dis[0].numele;
	  fluidfield->dis[1].element=(ELEMENT*)CCACALLOC(fluidfield->dis[1].numele,sizeof(ELEMENT));
         cpro++;
       } /* endif (cpro==0) */      
       fluidfield->dis[1].element[counter].eltyp=el_fluid2_tu;
       f2tu_dis(&(fluidfield->dis[0].element[counter]),&(fluidfield->dis[1].element[counter]));       
      } /* endif (e.f2->turbu == 2 || 3) */      
   }
#endif
/*----------------------------------------------------------------------*/
read:
   counter++;
   frread();
}
frrewind();
/*----------------------------------------------------------------------*/

end:
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

/*----------------------------------------- allocate one discretization */
/*alefield->ndis=1;*/
if (alefield->ndis>1) 
   dserror("different discretisations not implemented yet for structural elements\n");
alefield->dis = (DISCRET*)CCACALLOC(alefield->ndis,sizeof(DISCRET));
/*-------------------------------------------- count number of elements */
if (frfind("--ALE ELEMENTS")==1)
{
  frread();
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    counter++;
    frread();
  }
}
frrewind();
alefield->dis[0].numele = counter;
/*--------------------------------------------------- allocate elements */
alefield->dis[0].element=(ELEMENT*)CCACALLOC(alefield->dis[0].numele,sizeof(ELEMENT));
/*------------------------------------------------------- read elements */
frrewind();
if (frfind("--ALE ELEMENTS")==0) goto end;
frread();
counter=0;
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   colpointer = allfiles.actplace;
   elenumber  = strtol(colpointer,&colpointer,10);
   alefield->dis[0].element[counter].Id = --elenumber;
/*---------- read the typ of element and call element reading function */
/*------------------------------------------------ elementtyp is ALE3 */
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
      alefield->dis[0].element[counter].eltyp=el_ale3;
      ale3inp(&(alefield->dis[0].element[counter]));
   }
#endif
/*------------------------------------------------ elementtyp is ALE2 */
   frchk("ALE2",&ierr);
   if (ierr==1)
   {
#ifndef D_ALE 
      dserror("ALE2 needed but not defined in Makefile");
#endif
   }
#ifdef D_ALE 
   if (ierr==1) 
   {
      alefield->dis[0].element[counter].eltyp=el_ale2;
      ale2inp(&(alefield->dis[0].element[counter]));
   }
#endif
/*------------------------------------------------------- elementtyp is */
   counter++;
   frread();
}
frrewind();
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_ale_field */

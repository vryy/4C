#include "../headers/standardtypes.h"
#include "../shell8/shell8.h"
#include "../ale/ale.h"
/*----------------------------------------------------------------------*
 |  routine to write solution to GID                     m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_sol_init(FIELD *actfield)
{
int           i,j,k,l,n;
FILE *out     = allfiles.gidres;
#ifdef DEBUG 
dstrc_enter("out_gid_sol_init");
#endif
/*----------------------------------------------------------------------*/
fprintf(out,"Gid Post Results File 1.0\n");
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"# P_CARAT postprocessing output to GID\n");
fprintf(out,"#-------------------------------------------------------------------------------\n");
/*----------------------------------------------------------------------*/
fflush(out);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_gid_sol_init */




/*----------------------------------------------------------------------*
 |  routine to write a result range table to solution file  m.gee 12/01 |
 *----------------------------------------------------------------------*/
void out_gid_sol_setrangetable(FIELD *actfield, 
                               char string[],
                               double from,
                               double to,
                               int    steps)
{
int           i,j,k,l,n;

double        interval;

FILE *out     = allfiles.gidres;
int           strlenght;
char          sign='"';
#ifdef DEBUG 
dstrc_enter("out_gid_sol_setrangetable");
#endif
/*----------------------------------------------------------------------*/
strlenght = strlen(string);
interval  = (to-from)/(double)steps;
/*----------------------------------------------------------------------*/
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"# GiD postprocessing Range Table %s\n",string);
fprintf(out,"#-------------------------------------------------------------------------------\n");
/*----------------------------------------- result type is displacement */
fprintf(out,"RESULTRANGESTABLE %c%s%c\n",sign,string,sign);
fprintf(out,"                    - %-18.5#f : %crange0%c\n",from,sign,sign);
for (i=0; i<steps; i++)
{
fprintf(out," %-18.5#f - %-18.5#f : %crange%d%c\n",from+i*interval,
                                                 from+(i+1)*interval,
                                                 sign,i+1,sign); 
}
fprintf(out," %-18.5#f -                    : %crange%d%c\n",to,sign,i+1,sign); 
fprintf(out,"END RESULTRANGESTABLE\n");
/*----------------------------------------------------------------------*/
fflush(out);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_gid_sol_setrangetable */



/*----------------------------------------------------------------------*
 |  routine to define a gausspoint set for Gid postproc     m.gee 12/01 |
 *----------------------------------------------------------------------*/
void out_gid_sol_gausspointsets(FIELD *actfield)
{
int           i,j;
ELEMENT      *actele;
ELEMENT      *ele_quad4=NULL;
ELEMENT      *ele_quad8=NULL;
ELEMENT      *ele_quad9=NULL;
ELEMENT      *ele_tri3 =NULL;
ELEMENT      *ele_tri6 =NULL;
ELEMENT      *ele_hex8 =NULL;
ELEMENT      *ele_hex20=NULL;
ELEMENT      *ele_hex27=NULL;
ELEMENT      *ele_tet4 =NULL;
ELEMENT      *ele_tet10=NULL;


FILE *out     = allfiles.gidres;

int           nummesh;
int          *discrets;
int           ndiscrets=10;
DIS_TYP       distyp;
int          *nGP;
int           numGP;
FIELDTYP     fieldtyp;

char          sign='"';
char         *meshname  = "           ";

#ifdef DEBUG 
dstrc_enter("out_gid_sol_gausspointsets");
#endif
/*----------------------------------------------------------------------*/
fieldtyp = actfield->fieldtyp;
/*----------------------------------------------------------------------*/
discrets = (int*)calloc(ndiscrets,sizeof(int));
if (!discrets) dserror("Allocation of memory failed");
/*------------ check how many types of discrets there are in this field */
for (j=0; j<actfield->numele; j++)
{
   actele = &(actfield->element[j]);
   switch(actele->distyp)
   {
   case quad4: discrets[0]=1; ele_quad4 = actele; break;
   case quad8: discrets[1]=1; ele_quad8 = actele; break;
   case quad9: discrets[2]=1; ele_quad9 = actele; break;
   case tri3 : discrets[3]=1; ele_tri3  = actele; break;
   case tri6 : discrets[4]=1; ele_tri6  = actele; break;
   case hex8 : discrets[5]=1; ele_hex8  = actele; break;
   case hex20: discrets[6]=1; ele_hex20 = actele; break;
   case hex27: discrets[7]=1; ele_hex27 = actele; break;
   case tet4 : discrets[8]=1; ele_tet4  = actele; break;
   case tet10: discrets[9]=1; ele_tet10 = actele; break;
   default:
      dserror("Unknown type of element");
   break;
   }
}
nummesh=0;
for (j=0; j<ndiscrets; j++) 
if (discrets[j]) nummesh++;
/*----------------------------------------------------------------------*/
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"# GiD postprocessing GaussPoint Sets\n");
fprintf(out,"#-------------------------------------------------------------------------------\n");
/*--------------------------------------- loop i as meshes in one field */
for (i=0; i<nummesh; i++)
{
   fprintf(out,"# GaussPoint set %d\n",i);
   sprintf(meshname,"           ");
   switch (fieldtyp)
   {
   case structure:
   fprintf(out,"GAUSSPOINTS %cstructure%d%c ",sign,i,sign);
   sprintf(meshname,"structure%d",i);
   break;
   case fluid:
   fprintf(out,"GAUSSPOINTS %cfluid%d%c ",sign,i,sign);
   sprintf(meshname,"fluid%d",i);
   break;
   case ale:
   fprintf(out,"GAUSSPOINTS %cale%d%c ",sign,i,sign);
   sprintf(meshname,"ale%d",i);
   break;
   default:
      dserror("Unknown type of field");
   break;
   }
/*----------------------------------------------------- check for quad4 */ 
   if (discrets[0])
   {
      distyp=quad4;
      fprintf(out,"ELEMTYPE Quadrilateral %s \n",meshname);
      numGP=4;
      /*---------------------------- print number of Gaussian points */
      fprintf(out,"NUMBER OF GAUSS POINTS: %d\n",numGP);
      fprintf(out,"NATURAL COORDINATES: Internal\n");
      fprintf(out,"END GAUSSPOINTS\n");
      /*------------------ set flag, that this mesh has been printed */
      discrets[0]=0;
      goto nextmesh;
   }
/*----------------------------------------------------- check for quad8 */ 
   if (discrets[1])
   {
      distyp=quad8;
      fprintf(out,"ELEMTYPE Quadrilateral %s \n",meshname);
      numGP=9;
      /*---------------------------- print number of Gaussian points */
      fprintf(out,"NUMBER OF GAUSS POINTS: %d\n",numGP);
      fprintf(out,"NATURAL COORDINATES: Internal\n");
      fprintf(out,"END GAUSSPOINTS\n");
      /*------------------ set flag, that this mesh has been printed */
      discrets[1]=0;
      goto nextmesh;
   }
/*----------------------------------------------------- check for quad9 */ 
   if (discrets[2])
   {
      distyp=quad9;
      fprintf(out,"ELEMTYPE Quadrilateral %s \n",meshname);
      numGP=9;
      /*---------------------------- print number of Gaussian points */
      fprintf(out,"NUMBER OF GAUSS POINTS: %d\n",numGP);
      fprintf(out,"NATURAL COORDINATES: Internal\n");
      fprintf(out,"END GAUSSPOINTS\n");
      /*------------------ set flag, that this mesh has been printed */
      discrets[2]=0;
      goto nextmesh;
   }
/*------------------------------------------------------ check for tri3 */ 
   if (discrets[3])
   {
      distyp=tri3;
      fprintf(out,"ELEMTYPE Quadrilateral %s \n",meshname);
      numGP=1;
      /*---------------------------- print number of Gaussian points */
      fprintf(out,"NUMBER OF GAUSS POINTS: %d\n",numGP);
      fprintf(out,"NATURAL COORDINATES: Internal\n");
      fprintf(out,"END GAUSSPOINTS\n");
      /*------------------ set flag, that this mesh has been printed */
      discrets[3]=0;
      goto nextmesh;
   }
/*------------------------------------------------------ check for tri6 */ 
   if (discrets[4])
   {
      distyp=tri6;
      fprintf(out,"ELEMTYPE Quadrilateral %s \n",meshname);
      numGP=3;
      /*---------------------------- print number of Gaussian points */
      fprintf(out,"NUMBER OF GAUSS POINTS: %d\n",numGP);
      fprintf(out,"NATURAL COORDINATES: Internal\n");
      fprintf(out,"END GAUSSPOINTS\n");
      /*------------------ set flag, that this mesh has been printed */
      discrets[4]=0;
      goto nextmesh;
   }
/*------------------------------------------------------ check for hex8 */ 
   if (discrets[5])
   {
      distyp=hex8;
      fprintf(out,"ELEMTYPE Quadrilateral %s \n",meshname);
      numGP=8;
      /*---------------------------- print number of Gaussian points */
      fprintf(out,"NUMBER OF GAUSS POINTS: %d\n",numGP);
      fprintf(out,"NATURAL COORDINATES: Internal\n");
      fprintf(out,"END GAUSSPOINTS\n");
      /*------------------ set flag, that this mesh has been printed */
      discrets[5]=0;
      goto nextmesh;
   }
/*----------------------------------------------------- check for hex20 */ 
   if (discrets[6])
   {
      distyp=hex20;
      fprintf(out,"ELEMTYPE Quadrilateral %s \n",meshname);
      numGP=27;
      /*---------------------------- print number of Gaussian points */
      fprintf(out,"NUMBER OF GAUSS POINTS: %d\n",numGP);
      fprintf(out,"NATURAL COORDINATES: Internal\n");
      fprintf(out,"END GAUSSPOINTS\n");
      /*------------------ set flag, that this mesh has been printed */
      discrets[6]=0;
      goto nextmesh;
   }
/*----------------------------------------------------- check for hex27 */ 
   if (discrets[7])
   {
      distyp=hex27;
      fprintf(out,"ELEMTYPE Quadrilateral %s \n",meshname);
      numGP=27;
      /*---------------------------- print number of Gaussian points */
      fprintf(out,"NUMBER OF GAUSS POINTS: %d\n",numGP);
      fprintf(out,"NATURAL COORDINATES: Internal\n");
      fprintf(out,"END GAUSSPOINTS\n");
      /*------------------ set flag, that this mesh has been printed */
      discrets[7]=0;
      goto nextmesh;
   }
/*------------------------------------------------------ check for tet4 */ 
   if (discrets[8])
   {
      distyp=tet4;
      fprintf(out,"ELEMTYPE Quadrilateral %s \n",meshname);
      numGP=4;
      /*---------------------------- print number of Gaussian points */
      fprintf(out,"NUMBER OF GAUSS POINTS: %d\n",numGP);
      fprintf(out,"NATURAL COORDINATES: Internal\n");
      fprintf(out,"END GAUSSPOINTS\n");
      /*------------------ set flag, that this mesh has been printed */
      discrets[8]=0;
      goto nextmesh;
   }
/*----------------------------------------------------- check for tet10 */ 
   if (discrets[9])
   {
      distyp=tet10;
      fprintf(out,"ELEMTYPE Quadrilateral %s \n",meshname);
      numGP=10;
      /*---------------------------- print number of Gaussian points */
      fprintf(out,"NUMBER OF GAUSS POINTS: %d\n",numGP);
      fprintf(out,"NATURAL COORDINATES: Internal\n");
      fprintf(out,"END GAUSSPOINTS\n");
      /*------------------ set flag, that this mesh has been printed */
      discrets[9]=0;
      goto nextmesh;
   }
/*----------------------------------------------------------------------*/
nextmesh:;
}
/*----------------------------------------------------------------------*/
free(discrets);
fflush(out);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_gid_sol_gausspointsets */



/*----------------------------------------------------------------------*
 |  routine to write a result range table to solution file  m.gee 12/01 |
 *----------------------------------------------------------------------*/
void out_gid_domains(FIELD *actfield)
{
int           i,j,k,l,n;

ELEMENT      *actele;
ELEMENT      *ele_quad4=NULL;
ELEMENT      *ele_quad8=NULL;
ELEMENT      *ele_quad9=NULL;
ELEMENT      *ele_tri3 =NULL;
ELEMENT      *ele_tri6 =NULL;
ELEMENT      *ele_hex8 =NULL;
ELEMENT      *ele_hex20=NULL;
ELEMENT      *ele_hex27=NULL;
ELEMENT      *ele_tet4 =NULL;
ELEMENT      *ele_tet10=NULL;

int           nummesh;
int          *discrets;
int           ndiscrets=10;
DIS_TYP       distyp;
int          *nGP;
int           numGP;
FIELDTYP     fieldtyp;

char          sign='"';
char         *meshname  = "           ";

FILE *out     = allfiles.gidres;

#ifdef DEBUG 
dstrc_enter("out_gid_domains");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
fieldtyp = actfield->fieldtyp;
/*----------------------------------------------------------------------*/
discrets = (int*)calloc(ndiscrets,sizeof(int));
if (!discrets) dserror("Allocation of memory failed");
/*------------ check how many types of discrets there are in this field */
for (j=0; j<actfield->numele; j++)
{
   actele = &(actfield->element[j]);
   switch(actele->distyp)
   {
   case quad4: discrets[0]=1; ele_quad4 = actele; break;
   case quad8: discrets[1]=1; ele_quad8 = actele; break;
   case quad9: discrets[2]=1; ele_quad9 = actele; break;
   case tri3 : discrets[3]=1; ele_tri3  = actele; break;
   case tri6 : discrets[4]=1; ele_tri6  = actele; break;
   case hex8 : discrets[5]=1; ele_hex8  = actele; break;
   case hex20: discrets[6]=1; ele_hex20 = actele; break;
   case hex27: discrets[7]=1; ele_hex27 = actele; break;
   case tet4 : discrets[8]=1; ele_tet4  = actele; break;
   case tet10: discrets[9]=1; ele_tet10 = actele; break;
   default:
      dserror("Unknown type of element");
   break;
   }
}
nummesh=0;
for (j=0; j<ndiscrets; j++) 
if (discrets[j]) nummesh++;
/*----------------------------------------------------------------------*/
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"# Domains\n");
fprintf(out,"#-------------------------------------------------------------------------------\n");
/*--------------------------------------- loop i as meshes in one field */
for (i=0; i<nummesh; i++)
{
   sprintf(meshname,"           ");
   switch (fieldtyp)
   {
   case structure:
   fprintf(out,"RESULT %cDomains%c %cp_carat%c 0 Scalar OnGaussPoints %cstructure%d%c \n",sign,
                                                                                        sign,
                                                                                        sign,
                                                                                        sign,
                                                                                        sign,
                                                                                        i,
                                                                                        sign);
   sprintf(meshname,"structure%d",i);
   break;
   case fluid:
   fprintf(out,"RESULT %cDomains%c %cp_carat%c 0 Scalar OnGaussPoints %cfluid%d%c \n",sign,
                                                                                        sign,
                                                                                        sign,
                                                                                        sign,
                                                                                        sign,
                                                                                        i,
                                                                                        sign);
   sprintf(meshname,"fluid%d",i);
   break;
   case ale:
   fprintf(out,"RESULT %cDomains%c %cp_carat%c 0 Scalar OnGaussPoints %cale%d%c \n",sign,
                                                                                        sign,
                                                                                        sign,
                                                                                        sign,
                                                                                        sign,
                                                                                        i,
                                                                                        sign);
   sprintf(meshname,"ale%d",i);
   break;
   default:
      dserror("Unknown type of field");
   break;
   }
/*----------------------------------------------------- check for quad4 */ 
   if (discrets[0])
   {
      distyp=quad4;
      numGP=4;
      fprintf(out,"VALUES\n");
      /*---------------------------------------- print owner of element */
      for (j=0; j<actfield->numele; j++)
      {
         actele = &(actfield->element[j]);
         if (actele->distyp != distyp) continue;
         fprintf(out," %6d  %18.5#E \n",actele->Id+1,(double)(actele->proc));
         for (k=0; k<numGP-1; k++)
         fprintf(out,"         %-18.5#f \n",(double)(actele->proc));
      }
      /*------------------ set flag, that this mesh has been printed */
      fprintf(out,"END VALUES\n");
      discrets[0]=0;
      goto nextmesh;
   }
/*----------------------------------------------------- check for quad8 */ 
   if (discrets[1])
   {
      distyp=quad8;
      numGP=9;
      fprintf(out,"VALUES\n");
      /*---------------------------------------- print owner of element */
      for (j=0; j<actfield->numele; j++)
      {
         actele = &(actfield->element[j]);
         if (actele->distyp != distyp) continue;
         fprintf(out," %6d  %-18.5#f \n",actele->Id+1,(double)(actele->proc));
         for (k=0; k<numGP-1; k++)
         fprintf(out,"         %-18.5#f \n",(double)(actele->proc));
      }
      /*------------------ set flag, that this mesh has been printed */
      fprintf(out,"END VALUES\n");
      discrets[1]=0;
      goto nextmesh;
   }
/*----------------------------------------------------- check for quad9 */ 
   if (discrets[2])
   {
      distyp=quad9;
      numGP=9;
      fprintf(out,"VALUES\n");
      /*---------------------------------------- print owner of element */
      for (j=0; j<actfield->numele; j++)
      {
         actele = &(actfield->element[j]);
         if (actele->distyp != distyp) continue;
         fprintf(out," %6d  %-18.5#f \n",actele->Id+1,(double)(actele->proc));
         for (k=0; k<numGP-1; k++)
         fprintf(out,"         %-18.5#f \n",(double)(actele->proc));
      }
      /*------------------ set flag, that this mesh has been printed */
      fprintf(out,"END VALUES\n");
      discrets[2]=0;
      goto nextmesh;
   }
/*------------------------------------------------------ check for tri3 */ 
   if (discrets[3])
   {
      distyp=tri3;
      numGP=1;
      fprintf(out,"VALUES\n");
      /*---------------------------------------- print owner of element */
      for (j=0; j<actfield->numele; j++)
      {
         actele = &(actfield->element[j]);
         if (actele->distyp != distyp) continue;
         fprintf(out," %6d  %-18.5#f \n",actele->Id+1,(double)(actele->proc));
         for (k=0; k<numGP-1; k++)
         fprintf(out,"         %-18.5#f \n",(double)(actele->proc));
      }
      /*------------------ set flag, that this mesh has been printed */
      fprintf(out,"END VALUES\n");
      discrets[3]=0;
      goto nextmesh;
   }
/*------------------------------------------------------ check for tri6 */ 
   if (discrets[4])
   {
      distyp=tri6;
      numGP=3;
      fprintf(out,"VALUES\n");
      /*---------------------------------------- print owner of element */
      for (j=0; j<actfield->numele; j++)
      {
         actele = &(actfield->element[j]);
         if (actele->distyp != distyp) continue;
         fprintf(out," %6d  %-18.5#f \n",actele->Id+1,(double)(actele->proc));
         for (k=0; k<numGP-1; k++)
         fprintf(out,"         %-18.5#f \n",(double)(actele->proc));
      }
      /*------------------ set flag, that this mesh has been printed */
      fprintf(out,"END VALUES\n");
      discrets[4]=0;
      goto nextmesh;
   }
/*------------------------------------------------------ check for hex8 */ 
   if (discrets[5])
   {
      distyp=hex8;
      numGP=8;
      fprintf(out,"VALUES\n");
      /*---------------------------------------- print owner of element */
      for (j=0; j<actfield->numele; j++)
      {
         actele = &(actfield->element[j]);
         if (actele->distyp != distyp) continue;
         fprintf(out," %6d  %-18.5#f \n",actele->Id+1,(double)(actele->proc));
         for (k=0; k<numGP-1; k++)
         fprintf(out,"         %-18.5#f \n",(double)(actele->proc));
      }
      /*------------------ set flag, that this mesh has been printed */
      fprintf(out,"END VALUES\n");
      discrets[5]=0;
      goto nextmesh;
   }
/*----------------------------------------------------- check for hex20 */ 
   if (discrets[6])
   {
      distyp=hex20;
      numGP=27;
      fprintf(out,"VALUES\n");
      /*---------------------------------------- print owner of element */
      for (j=0; j<actfield->numele; j++)
      {
         actele = &(actfield->element[j]);
         if (actele->distyp != distyp) continue;
         fprintf(out," %6d  %-18.5#f \n",actele->Id+1,(double)(actele->proc));
         for (k=0; k<numGP-1; k++)
         fprintf(out,"         %-18.5#f \n",(double)(actele->proc));
      }
      /*------------------ set flag, that this mesh has been printed */
      fprintf(out,"END VALUES\n");
      discrets[6]=0;
      goto nextmesh;
   }
/*----------------------------------------------------- check for hex27 */ 
   if (discrets[7])
   {
      distyp=hex27;
      numGP=27;
      fprintf(out,"VALUES\n");
      /*---------------------------------------- print owner of element */
      for (j=0; j<actfield->numele; j++)
      {
         actele = &(actfield->element[j]);
         if (actele->distyp != distyp) continue;
         fprintf(out," %6d  %-18.5#f \n",actele->Id+1,(double)(actele->proc));
         for (k=0; k<numGP-1; k++)
         fprintf(out,"         %-18.5#f \n",(double)(actele->proc));
      }
      /*------------------ set flag, that this mesh has been printed */
      fprintf(out,"END VALUES\n");
      discrets[7]=0;
      goto nextmesh;
   }
/*------------------------------------------------------ check for tet4 */ 
   if (discrets[8])
   {
      distyp=tet4;
      numGP=4;
      fprintf(out,"VALUES\n");
      /*---------------------------------------- print owner of element */
      for (j=0; j<actfield->numele; j++)
      {
         actele = &(actfield->element[j]);
         if (actele->distyp != distyp) continue;
         fprintf(out," %6d  %-18.5#f \n",actele->Id+1,(double)(actele->proc));
         for (k=0; k<numGP-1; k++)
         fprintf(out,"         %-18.5#f \n",(double)(actele->proc));
      }
      /*------------------ set flag, that this mesh has been printed */
      fprintf(out,"END VALUES\n");
      discrets[8]=0;
      goto nextmesh;
   }
/*----------------------------------------------------- check for tet10 */ 
   if (discrets[9])
   {
      distyp=tet10;
      numGP=10;
      fprintf(out,"VALUES\n");
      /*---------------------------------------- print owner of element */
      for (j=0; j<actfield->numele; j++)
      {
         actele = &(actfield->element[j]);
         if (actele->distyp != distyp) continue;
         fprintf(out," %6d  %-18.5#f \n",actele->Id+1,(double)(actele->proc));
         for (k=0; k<numGP-1; k++)
         fprintf(out,"         %-18.5#f \n",(double)(actele->proc));
      }
      /*------------------ set flag, that this mesh has been printed */
      fprintf(out,"END VALUES\n");
      discrets[9]=0;
      goto nextmesh;
   }
/*----------------------------------------------------------------------*/
nextmesh:;
}
/*----------------------------------------------------------------------*/
free(discrets);
fflush(out);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_gid_domains */






/*----------------------------------------------------------------------*
 |  routine to write solution of a step to GID           m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_sol(char string[], FIELD *actfield, INTRA  *actintra, int step)
{
int           i,j;
int           counter;

FILE         *out     = allfiles.gidres;

int           strlenght;
int           ncomponents;

char          sign='"';
char         *vector  ="VECTOR";
char         *scalar  ="SCALAR";
char         *onnodes ="ONNODES";
char         *ongausspoints ="ONGAUSSPOINTS";

char         *my_result_type;
char         *my_result_location;
char         *my_location_name   ;
char         *my_component_name1 ;
char         *my_component_name2 ;
char         *my_component_name3 ;
char         *my_component_name4 ;
char         *my_component_name5 ;
char         *my_component_name6 ;
char         *my_component_name7 ;
char         *my_component_name8 ;
char         *my_component_name9 ;
char         *my_component_name10;
char         *my_component_name11;
char         *my_component_name12;
char         *my_component_name13;
char         *my_component_name14;
char         *my_component_name15;
char         *my_component_name16;
char         *my_component_name17;
char         *my_component_name18;

NODE         *actnode;
int           nodeId;
double        disp[3];

#ifdef DEBUG 
dstrc_enter("out_gid_sol");
#endif
/*----------------------------------------------------------------------*/
strlenght = strlen(string);
/*==========================================result type is displacement */
if (strncmp(string,"displacement",strlenght)==0)
{
   /*--------------------------- set type of result */
   my_result_type     = vector;

   /*-------------------------- set result location */
   my_result_location = onnodes;

   /*--set location name (only for Gaussian points) */
   my_location_name = "";
   /*------ number of components is ncomponentnames */
   ncomponents=3;
   my_component_name1 = "X-Disp";
   my_component_name2 = "Y-Disp";
   my_component_name3 = "Z-Disp";
}/*====================================================================*/
/*================================================result type is stress */
if (strncmp(string,"stress______",strlenght)==0)
{
   /*--------------------------- set type of result */
   my_result_type     = scalar;

   /*-------------------------- set result location */
   my_result_location = ongausspoints;

   /*--set location name (only for Gaussian points) */
   my_location_name = "";
   /*------ number of components is ncomponentnames */
   ncomponents=3;
   my_component_name1 = "X-Disp";
   my_component_name2 = "Y-Disp";
   my_component_name3 = "Z-Disp";
}/*====================================================================*/
/*--------------------------------------------------now start printing */
/*---------------------------------------------- print the result head */
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"RESULT %c%s%c %cp_carat%c %d %s %s %c%s%c\n",
                                           sign,string,sign,
                                           sign,sign,
                                           step,
                                           my_result_type,
                                           my_result_location,
                                           sign,my_location_name,sign);
/*---------------------------------- print the result ranges table name */
fprintf(out,"# the result ranges table used for this result:\n");
fprintf(out,"RESULTRANGESTABLE %c%s%c\n",sign,string,sign);
/*------------------------------------------- print the component names */
fprintf(out,"# the component names used for this result:\n");
switch(ncomponents)
{
case 1:
break;
case 2:
break;
case 3:
fprintf(out,"COMPONENTNAMES %c%s%c, %c%s%c, %c%s%c\n",sign,my_component_name1,sign,
                                                      sign,my_component_name2,sign,
                                                      sign,my_component_name3,sign);
break;
case 4:
break;
case 5:
break;
case 6:
break;
case 7:
break;
case 8:
break;
case 9:
break;
case 10:
break;
case 11:
break;
case 12:
break;
default:
   dserror("don't know how many components to print to GiD");
break;
}
/*----------------------------------------------------- get the results */
/*----------------------------------------- result type is displacement */
if (strncmp(string,"displacement",strlenght)==0)
{
   fprintf(out,"# the displacements values\n");
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->numnp; i++)
   {
      actnode = &(actfield->node[i]);
      nodeId  = actnode->Id+1;
      disp[0] = actnode->sol.a.da[step][0];
      disp[1] = actnode->sol.a.da[step][1];
      disp[2] = actnode->sol.a.da[step][2];
      fprintf(out," %6d %-18.5#f %-18.5#f %-18.5#f \n",nodeId,disp[0],disp[1],disp[2]);
   }
   fprintf(out,"END VALUES\n");  
}
/*----------------------------------------------------------------------*/
fflush(out);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_gid_sol */

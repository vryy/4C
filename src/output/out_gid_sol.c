#include "../headers/standardtypes.h"
#include "gid.h"
#include "../shell8/shell8.h"
#include "../wall1/wall1.h"
#include "../brick1/brick1.h"
#include "../ale2/ale2.h"
#include "../ale3/ale3.h"
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
 |  routine to write solution to GID                     m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_sol_init()
{
int               i,j;
FIELD            *actfield;
GIDSET           *actgid;
ELEMENT          *actele;
FILE             *out       = allfiles.gidres;
char              sign='"';
char             *charptr;

#ifdef DEBUG 
dstrc_enter("out_gid_sol_init");
#endif
/*----------------------------------------------------------------------*/
fprintf(out,"Gid Post Results File 1.0\n");
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"# P_CARAT postprocessing output to GID\n");
fprintf(out,"#-------------------------------------------------------------------------------\n");
/*----------------------------------------------------------------------*/
/*------------------------- check number of fields and allocate storage */
gid = (GIDSET*)CCACALLOC(genprob.numfld,sizeof(GIDSET));
if (!gid) dserror("Allocation of memory failed");
/*--------------------------------------- loop all fields and init data */
for (i=0; i<genprob.numfld; i++)
{
   actfield         = &(field[i]);
   actgid           = &(gid[i]);
   actgid->fieldtyp = actfield->fieldtyp;
   /*------------------------------------------------ set the fieldname */
   switch (actgid->fieldtyp)
   {
   case structure:
      actgid->fieldnamelenght = 9;
      actgid->fieldname       = "structure";
   break;
   case fluid:
      actgid->fieldnamelenght = 5;
      actgid->fieldname       = "fluid";
   break;
   case ale:
      actgid->fieldnamelenght = 3;
      actgid->fieldname       = "ale";
   break;
   default:
      dserror("Unknown type of field");
   break;
   }
   /*------------------------------------------------- set range tables */
   strncpy(actgid->standardrangetable,"standard_         ",18);
   charptr = actgid->standardrangetable + 9;
   strncpy(charptr,actgid->fieldname,actgid->fieldnamelenght);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RANGETABLES %s\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULTRANGESTABLE %c%s%c\n",sign,actgid->standardrangetable,sign);
   fprintf(out,"            - -1000000.0 : %cvery small%c\n",sign,sign);
   fprintf(out," -1000000.0 -  1000000.0 : %cnormal%c\n",sign,sign);
   fprintf(out,"  1000000.0 -            : %cvery large%c\n",sign,sign);
   fprintf(out,"END RESULTRANGESTABLE\n");
   /*--------------------------- find and set meshes and gausspointsets */
   actgid->is_shell8_22 = 0;
   actgid->is_shell8_33 = 0;
   actgid->is_brick1_222= 0;
   actgid->is_brick1_333= 0;
   actgid->is_fluid2_22 = 0;
   actgid->is_fluid2_33 = 0;
   actgid->is_fluid3_222= 0;
   actgid->is_fluid3_333= 0;
   actgid->is_ale_11    = 0;
   actgid->is_ale_22    = 0;
   actgid->is_ale_tri_1 = 0;
   actgid->is_ale_tri_3 = 0;
   actgid->is_ale_111   = 0;
   actgid->is_ale_222   = 0;
   actgid->is_ale_tet_1 = 0;
   actgid->is_ale_tet_4 = 0;
   actgid->is_wall1_22  = 0;
   actgid->is_wall1_33  = 0;
   /*---------------------------- check for different types of elements */
   for (j=0; j<actfield->dis[0].numele; j++)
   {
      actele = &(actfield->dis[0].element[j]);
      switch(actele->eltyp)
      {
      case el_shell8: 
         if (actele->numnp==4)  
         {
            actgid->is_shell8_22   = 1;
            actgid->shell8_22_name = "shell8_22";
         }
         if (actele->numnp==8 || actele->numnp==9)  
         {
            actgid->is_shell8_33   = 1;
            actgid->shell8_33_name = "shell8_33";
         }
      break;
      case el_brick1: 
         if (actele->numnp==8)  
         {
            actgid->is_brick1_222   = 1;
            actgid->brick1_222_name = "brick1_222";
         }
         if (actele->numnp==20 || actele->numnp==27) 
         {
            actgid->is_brick1_333   = 1;
            actgid->brick1_333_name = "brick1_333"; 
         }
      break;
      case el_fluid2: 
         if (actele->numnp==4)  
         {
            actgid->is_fluid2_22    = 1;
            actgid->fluid2_22_name  = "fluid2_22";
         }
         if (actele->numnp==8  || actele->numnp==9) 
         {
            actgid->is_fluid2_33    = 1;
            actgid->fluid2_33_name  = "fluid2_33";
         }
      break;
      case el_fluid3: 
         if (actele->numnp==8)  
         {
            actgid->is_fluid3_222   = 1;
            actgid->fluid3_222_name = "fluid3_222";
         }
         if (actele->numnp==20 || actele->numnp==27) 
         {
            actgid->is_fluid3_333   = 1;
            actgid->fluid3_333_name = "fluid3_333";
         }
      break;
#ifdef D_ALE
      case el_ale2:    
         if (actele->numnp==4)  
         {
             if (actele->e.ale2->nGP[0]==1 && actele->e.ale2->nGP[1]==1 )
             {
                 actgid->is_ale_11    = 1;
                 actgid->ale_11_name  = "ale_11";
             }
             else if (actele->e.ale2->nGP[0]==2 && actele->e.ale2->nGP[1]==2 )
             {
                 actgid->is_ale_22    = 1;
                 actgid->ale_22_name  = "ale_22";
             }
         }
         if (actele->numnp==3)  
         {
             if ( actele->e.ale2->nGP[0] == 1)
             {
                 actgid->is_ale_tri_1    = 1;
                 actgid->ale_tri_1_name  = "ale_tri_1";
             }
             else if ( actele->e.ale2->nGP[0] == 3)
             {
                 actgid->is_ale_tri_3    = 1;
                 actgid->ale_tri_3_name  = "ale_tri_3";
             }
         }
      break;
      case el_ale3:    
         if (actele->numnp==8)  
         {
             if (actele->e.ale3->nGP[0]==1 && actele->e.ale3->nGP[1]==1 && actele->e.ale3->nGP[2]==1 )
             {
                 actgid->is_ale_111   = 1;
                 actgid->ale_111_name = "ale_111";
             }
             else if (actele->e.ale3->nGP[0]==2 && actele->e.ale3->nGP[1]==2 && actele->e.ale3->nGP[2]==2 )
             {
                 actgid->is_ale_222   = 1;
                 actgid->ale_222_name = "ale_222";
             }
         }
         if (actele->numnp==4)  
         {
             if ( actele->e.ale3->nGP[0] == 1)
             {
                 actgid->is_ale_tet_1    = 1;
                 actgid->ale_tet_1_name  = "ale_tet_1";
             }
             else if ( actele->e.ale3->nGP[0] == 4)
             {
                 actgid->is_ale_tet_4    = 1;
                 actgid->ale_tet_4_name  = "ale_tet_4";
             }
         }
      break;
#endif
/*---------------------------------------------------------fh 06/02----*/
      case el_wall1:    
         if (actele->numnp==4)  
         {
            actgid->is_wall1_22    = 1;
            actgid->wall1_22_name  = "wall1_22";
         }
         if (actele->numnp==8 | actele->numnp==9)  
         {
            actgid->is_wall1_33   = 1;
            actgid->wall1_33_name = "wall1_33";
         }
      break;
      default:
         dserror("Unknown type of element");
      break;
      }
   } /* end of (j=0; j<actfield->numele; j++) */
   /*----------------------------- now we can write the gausspoint sets */
#if 0/* this is the shell visualization using Quadrilateral */
   if (actgid->is_shell8_22)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s SHELL8 2x2 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral %c%s%c\n",
                                                                   sign,actgid->shell8_22_name,sign,
                                                                   sign,actgid->shell8_22_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 4\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
#endif
   /* this is the shell visualization using Hexahedra */
   if (actgid->is_shell8_22)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s SHELL8 2x2x2 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Hexahedra %c%s%c\n",
                                                                   sign,actgid->shell8_22_name,sign,
                                                                   sign,actgid->shell8_22_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 8\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   /*-------------------------------------------------*/
   if (actgid->is_shell8_33)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s SHELL8 3x3 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral %c%s%c\n",
                                                                    sign,actgid->shell8_33_name,sign,
                                                                    sign,actgid->shell8_33_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 9\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_brick1_222)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s BRICK1 2x2x2 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Hexahedra %c%s%c\n",
                                                                sign,actgid->brick1_222_name,sign,
                                                                sign,actgid->brick1_222_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 8\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_brick1_333)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s BRICK1 3x3x3 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Hexahedra %c%s%c\n",
                                                                sign,actgid->brick1_333_name,sign,
                                                                sign,actgid->brick1_333_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 27\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_fluid2_22)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s FLUID2 2x2 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral %c%s%c\n",
                                                                sign,actgid->fluid2_22_name,sign,
                                                                sign,actgid->fluid2_22_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 4\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_fluid2_33)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s FLUID2 3x3 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral %c%s%c\n",
                                                                sign,actgid->fluid2_33_name,sign,
                                                                sign,actgid->fluid2_33_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 9\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_fluid3_222)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s FLUID3 2x2x2 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Hexahedra %c%s%c\n",
                                                                sign,actgid->fluid3_222_name,sign,
                                                                sign,actgid->fluid3_222_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 8\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_fluid3_333)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s FLUID3 3x3x3 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Hexahedra %c%s%c\n",
                                                                sign,actgid->fluid3_333_name,sign,
                                                                sign,actgid->fluid3_333_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 27\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_ale_11)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s ALE 1x1 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral %c%s%c\n",
                                                                    sign,actgid->ale_11_name,sign,
                                                                    sign,actgid->ale_11_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 1\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_ale_22)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s ALE 2x2 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral %c%s%c\n",
                                                                    sign,actgid->ale_22_name,sign,
                                                                    sign,actgid->ale_22_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 4\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_ale_tri_1)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s ALE 1 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Triangle %c%s%c\n",
                                                                    sign,actgid->ale_tri_1_name,sign,
                                                                    sign,actgid->ale_tri_1_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 1\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_ale_tri_3)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s ALE 3 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Triangle %c%s%c\n",
                                                                    sign,actgid->ale_tri_3_name,sign,
                                                                    sign,actgid->ale_tri_3_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 3\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_ale_111)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s ALE 1x1x1 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%cELEMTYPE Hexahedra %c%s%c\n",
                                                               sign,actgid->ale_111_name,sign,
                                                               sign,actgid->ale_111_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 1\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_ale_222)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s ALE 2x2x2 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%cELEMTYPE Hexahedra %c%s%c\n",
                                                               sign,actgid->ale_222_name,sign,
                                                               sign,actgid->ale_222_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 8\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_ale_tet_1)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s ALE 1 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Tetrahedra %c%s%c\n",
                                                                    sign,actgid->ale_tet_1_name,sign,
                                                                    sign,actgid->ale_tet_1_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 1\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_ale_tet_4)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s ALE 4 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Tetrahedra %c%s%c\n",
                                                                    sign,actgid->ale_tet_4_name,sign,
                                                                    sign,actgid->ale_tet_4_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 4\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
/*---------------------------------------------------------fh 06/02----*/
   if (actgid->is_wall1_22)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s WALL1 2x2 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral %c%s%c\n",
                                                                   sign,actgid->wall1_22_name,sign,
                                                                   sign,actgid->wall1_22_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 4\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_wall1_33)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s WALL1 3x3 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral %c%s%c\n",
                                                                    sign,actgid->wall1_33_name,sign,
                                                                    sign,actgid->wall1_33_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 9\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
/*----------------------------------------------------------------------*/
} /* end of (i=0; i<genprob.numfld; i++) */
/*----------------------------------------------------------------------*/
fflush(out);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_gid_sol_init */



/*----------------------------------------------------------------------*
 |  write the domain numbers to GiD as gausspointvalues     m.gee 12/01 |
 *----------------------------------------------------------------------*/
void out_gid_domains(FIELD *actfield)
{
int           i,j,k;

FILE         *out = allfiles.gidres;
ELEMENT      *actele;
GIDSET       *actgid = NULL;
char          sign='"';

#ifdef DEBUG 
dstrc_enter("out_gid_domains");
#endif
/*----------------------------------------------------------------------*/
/*-------------------------------------- find the correct gid structure */
for (i=0; i<genprob.numfld; i++)
{
   if (gid[i].fieldtyp == actfield->fieldtyp)
   {
      actgid = &(gid[i]);
      break;
   }
}
if (!actgid) dserror("Cannot find correct field");
/*----------------------------------------------------------------------*/
if (actgid->is_shell8_22) 
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->shell8_22_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->shell8_22_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_shell8 || actele->numnp != 4) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(double)actele->proc);
#if 0
      for (j=1; j<4; j++)/* quadrilateral version */
#endif
      for (j=1; j<8; j++)/* hexahedra version */
      fprintf(out,"            %18.5E\n",(double)actele->proc); 
   }
   fprintf(out,"END VALUES\n");
}

if (actgid->is_shell8_33)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->shell8_33_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->shell8_33_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_shell8 || actele->numnp != 9) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(double)actele->proc);
      for (j=1; j<9; j++)
      fprintf(out,"            %18.5E\n",(double)actele->proc); 
   }
   fprintf(out,"END VALUES\n");
}

if (actgid->is_brick1_222)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->brick1_222_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->brick1_222_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_brick1 || actele->numnp != 8) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(double)actele->proc);
      for (j=1; j<8; j++)
      fprintf(out,"            %18.5E\n",(double)actele->proc); 
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_brick1_333)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->brick1_333_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->brick1_333_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_brick1 || (actele->numnp != 20 || actele->numnp != 27)) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(double)actele->proc);
      for (j=1; j<27; j++)
      fprintf(out,"            %18.5E\n",(double)actele->proc); 
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_fluid3_222)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->fluid3_222_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->fluid3_222_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_fluid3 || actele->numnp != 8) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(double)actele->proc);
      for (j=1; j<8; j++)
      fprintf(out,"            %18.5E\n",(double)actele->proc); 
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_fluid3_333)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->fluid3_333_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->fluid3_333_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_fluid3 || (actele->numnp != 20 || actele->numnp != 27)) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(double)actele->proc);
      for (j=1; j<27; j++)
      fprintf(out,"            %18.5E\n",(double)actele->proc); 
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_ale_11)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_11_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                                         sign,actgid->ale_11_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_ale2 || actele->numnp != 4) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(double)actele->proc);
      for (j=1; j<4; j++)
      fprintf(out,"            %18.5E\n",(double)actele->proc); 
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_ale_22)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_22_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                                         sign,actgid->ale_22_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_ale2 || actele->numnp != 4) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(double)actele->proc);
      for (j=1; j<4; j++)
      fprintf(out,"            %18.5E\n",(double)actele->proc); 
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_ale_tri_1)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_tri_1_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                                   sign,actgid->ale_tri_1_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_ale2 || actele->numnp != 3) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(double)actele->proc);
      for (j=1; j<3; j++)
      fprintf(out,"            %18.5E\n",(double)actele->proc); 
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_ale_tri_3)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_tri_3_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                                   sign,actgid->ale_tri_3_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_ale2 || actele->numnp != 3) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(double)actele->proc);
      for (j=1; j<3; j++)
      fprintf(out,"            %18.5E\n",(double)actele->proc); 
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_ale_111)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_111_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->ale_111_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_ale3 || actele->numnp != 8) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(double)actele->proc);
      for (j=1; j<8; j++)
      fprintf(out,"            %18.5E\n",(double)actele->proc); 
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_ale_222)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_222_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->ale_222_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_ale3 || actele->numnp != 8) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(double)actele->proc);
      for (j=1; j<8; j++)
      fprintf(out,"            %18.5E\n",(double)actele->proc); 
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_ale_tet_1)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_tet_1_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->ale_tet_1_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_ale3 || actele->numnp != 4) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(double)actele->proc);
      for (j=1; j<4; j++)
      fprintf(out,"            %18.5E\n",(double)actele->proc); 
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_ale_tet_4)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_tet_4_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->ale_tet_4_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_ale3 || actele->numnp != 4) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(double)actele->proc);
      for (j=1; j<4; j++)
      fprintf(out,"            %18.5E\n",(double)actele->proc); 
   }
   fprintf(out,"END VALUES\n");
}


/*---------------------------------------------------------fh 06/02----*/
if (actgid->is_wall1_22)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->wall1_22_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->wall1_22_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_wall1 || actele->numnp != 4) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(double)actele->proc);
      for (j=1; j<4; j++)
      fprintf(out,"            %18.5E\n",(double)actele->proc); 
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_wall1_33)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->wall1_33_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->wall1_33_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_wall1 || actele->numnp < 8) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(double)actele->proc);
      for (j=1; j<9; j++)
      fprintf(out,"            %18.5E\n",(double)actele->proc); 
   }
   fprintf(out,"END VALUES\n");
}
/*----------------------------------------------------------------------*/
fflush(out);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_gid_domains */






/*----------------------------------------------------------------------*
 |  routine to write solution of a step to GID           m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_sol(char string[], FIELD *actfield, INTRA  *actintra, int step,
                 int place)
{
int           i,j;

FILE         *out = allfiles.gidres;

NODE         *actnode;
ELEMENT      *actele;
GIDSET       *actgid = NULL;

char         *resulttype;
char         *resultplace;
char         *gpset;
char         *rangetable;
int           ncomponent;
char         *componentnames[18];

char          sign='"';

int           stringlenght;

int           ngauss;
double      **forces;
double      **stress;

double        a1,a2,a3,thick,scal,sdc;
int           tot_numnp;
/* 
   gausspoint permutation :
   On the Gausspoint number i in Gid, the results of Carats GP number gausspermn[i]
   have to be written
*/   

int           gaussperm4[4] = {3,1,0,2};
int           gaussperm8[8] = {0,4,2,6,1,5,3,7};
int           gaussperm9[9] = {8,2,0,6,5,1,3,7,4};
int           gaussperm27[27] = {0,9,18,3,12,21,6,15,24,1,10,19,4,13,22,7,16,25,2,11,20,5,14,23,8,17,26};

#ifdef DEBUG 
dstrc_enter("out_gid_sol");
#endif
/*----------------------------------------------------------------------*/
/*-------------------------------------- find the correct gid structure */
for (i=0; i<genprob.numfld; i++)
{
   if (gid[i].fieldtyp == actfield->fieldtyp)
   {
      actgid = &(gid[i]);
      break;
   }
}
if (!actgid) dserror("Cannot find correct field");
/*----------------------------------------------------------------------*/
stringlenght = strlen(string);
/*========================================= result type is displacement */
if (strncmp(string,"displacement",stringlenght)==0)
{
   resulttype        = "VECTOR";
   resultplace       = "ONNODES";
   gpset             = ""; 
   rangetable        = actgid->standardrangetable;
   ncomponent        = 3;
   componentnames[0] = "x-displ";
   componentnames[1] = "y-displ";
   componentnames[2] = "z-displ";
   /*-------------------------------------------------------------------*/
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT %s on FIELD %s\n",string,actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %c%s%c %cpcarat%c %d %s %s\n",
                                                             sign,string,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace
                                                             );
   fprintf(out,"RESULTRANGESTABLE %c%s%c\n",
                                            sign,actgid->standardrangetable,sign
                                            );
   switch (genprob.ndim)
   {
     case 3:
       fprintf(out,"COMPONENTNAMES %c%s%c,%c%s%c,%c%s%c\n",
                                                       sign,componentnames[0],sign,
                                                       sign,componentnames[1],sign,
                                                       sign,componentnames[2],sign
                                                       );
     break;
     case 2:
       fprintf(out,"COMPONENTNAMES %c%s%c,%c%s%c\n",
                                                       sign,componentnames[0],sign,
                                                       sign,componentnames[1],sign
                                                       );
     break;
     default:
       dserror("Unknown numer of dimensions");
     break;
   }
   fprintf(out,"VALUES\n");
#ifdef D_SHELL8
#if 1 /* this is hexahedra output */
   if (actfield->dis[0].element[0].eltyp == el_shell8 && actfield->dis[0].element[0].distyp == quad4)
   {
      tot_numnp = genprob.nnode;
      scal = 100.0;
      sdc  = actfield->dis[0].element[0].e.s8->sdc;
      for (i=0; i<actfield->dis[0].numnp; i++)
      {
         actnode = &(actfield->dis[0].node[i]);
         /* the lower surface */
         fprintf(out," %6d %23.15E %23.15E %23.15E\n",
                                                    actnode->Id+1,
                                                    actnode->sol.a.da[place][0]-actnode->sol.a.da[place][3]*scal/sdc,
                                                    actnode->sol.a.da[place][1]-actnode->sol.a.da[place][4]*scal/sdc,
                                                    actnode->sol.a.da[place][2]-actnode->sol.a.da[place][5]*scal/sdc
                                                    );
         /* the upper surface */
         fprintf(out," %6d %23.15E %23.15E %23.15E\n",
                                                    actnode->Id+1+tot_numnp,
                                                    actnode->sol.a.da[place][0]+actnode->sol.a.da[place][3]*scal/sdc,
                                                    actnode->sol.a.da[place][1]+actnode->sol.a.da[place][4]*scal/sdc,
                                                    actnode->sol.a.da[place][2]+actnode->sol.a.da[place][5]*scal/sdc
                                                    );
      }
   }
   else
#endif
#endif   
   for (i=0; i<actfield->dis[0].numnp; i++)
   {
      actnode = &(actfield->dis[0].node[i]);
      switch (genprob.ndim)
      {
	case 3:
        fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                                                   actnode->Id+1,
                                                   actnode->sol.a.da[place][0],
                                                   actnode->sol.a.da[place][1],
                                                   actnode->sol.a.da[place][2]
                                                   );
	break;
	case 2:
        fprintf(out," %6d %18.5E %18.5E \n",
                                                   actnode->Id+1,
                                                   actnode->sol.a.da[place][0],
                                                   actnode->sol.a.da[place][1]
                                                   );
        break;
	default:
	  dserror("Unknown number of dimensions");
        break;
      }
   }
   fprintf(out,"END VALUES\n");
} /* end of (strncmp(string,"displacement",stringlenght)==0) */
/*========================================= result type is eigenmodes */
if (strncmp(string,"eigenmodes",stringlenght)==0)
{
   resulttype        = "VECTOR";
   resultplace       = "ONNODES";
   gpset             = ""; 
   rangetable        = actgid->standardrangetable;
   ncomponent        = 3;
   componentnames[0] = "x-mode";
   componentnames[1] = "y-mode";
   componentnames[2] = "z-mode";
   /*-------------------------------------------------------------------*/
for (j=0; j<place; j++)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT %s on FIELD %s\n",string,actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %c%s%c %cpcarat%c %d %s %s\n",
                                                             sign,string,sign,
                                                             sign,sign,
                                                             j,
                                                             resulttype,
                                                             resultplace
                                                             );
   fprintf(out,"RESULTRANGESTABLE %c%s%c\n",
                                            sign,actgid->standardrangetable,sign
                                            );
   fprintf(out,"COMPONENTNAMES %c%s%c,%c%s%c,%c%s%c\n",
                                                       sign,componentnames[0],sign,
                                                       sign,componentnames[1],sign,
                                                       sign,componentnames[2],sign
                                                       );
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numnp; i++)
   {
      actnode = &(actfield->dis[0].node[i]);
      fprintf(out," %6d %18.5e %18.5e %18.5e\n",
                                                   actnode->Id+1,
                                                   actnode->sol.a.da[j][0],
                                                   actnode->sol.a.da[j][1],
                                                   actnode->sol.a.da[j][2]
                                                   );
   }
   fprintf(out,"END VALUES\n");
}
} /* end of (strncmp(string,"eigenmodes",stringlenght)==0) */
/*======================================================================*/
/*=============================================== result type is stress */
if (strncmp(string,"stress",stringlenght)==0)
{
   /*--------- ---------now go through the meshes and print the results */
   /*===============================shell8 element with 2x2 gausspoints */
   /* these shells have 18 stresses, do 2 x 3D Matrix */
   /* for shell8 stresses permutation: */
   /* ii[18] = {0,2,8,1,3,16,4,17,9,5,7,14,6,10,12,11,13,15};*/
#if 0
#ifdef D_SHELL8
   if (actgid->is_shell8_22)
   {
      ngauss=4;
      resulttype        = "MATRIX";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->shell8_22_name;
      rangetable        = actgid->standardrangetable;
      /*--- print the constant-in-thickness direction forces */
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT shell8_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cshell8_forces%c %cpcarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      fprintf(out,"COMPONENTNAMES %cForce-11%c,%cForce-22%c,%cForce-33%c,%cForce-12%c,%cForce-23%c,%cForce-13%c\n",
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign);
      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_shell8 || actele->numnp !=4) continue;
         forces = actele->e.s8->forces.a.d3[place];
         fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                             actele->Id+1,
                             forces[0][gaussperm4[0]],
                             forces[1][gaussperm4[0]],
                             forces[9][gaussperm4[0]],
                             forces[2][gaussperm4[0]],
                             forces[4][gaussperm4[0]],
                             forces[3][gaussperm4[0]]
                             );
         for (j=1; j<ngauss; j++)
         fprintf(out,"        %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                             forces[0][gaussperm4[j]],
                             forces[1][gaussperm4[j]],
                             forces[9][gaussperm4[j]],
                             forces[2][gaussperm4[j]],
                             forces[4][gaussperm4[j]],
                             forces[3][gaussperm4[j]]
                             );
 
      }
      fprintf(out,"END VALUES\n");
      /*--- print the linear-in-thickness direction forces */
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT shell8_moments on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cshell8_moments%c %cpcarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      fprintf(out,"COMPONENTNAMES %cMoment-11%c,%cMoment-22%c,%cMoment-33%c,%cMoment-12%c,%cMoment-23%c,%cMoment-13%c\n",
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign);
      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_shell8 || actele->numnp !=4) continue;
         forces = actele->e.s8->forces.a.d3[place];
         fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                             actele->Id+1,
                             forces[5] [gaussperm4[0]],
                             forces[6] [gaussperm4[0]],
                             forces[15][gaussperm4[0]],
                             forces[7] [gaussperm4[0]],
                             forces[11][gaussperm4[0]],
                             forces[10][gaussperm4[0]]
                             );
         for (j=1; j<ngauss; j++)
         fprintf(out,"        %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                             forces[5] [gaussperm4[j]],
                             forces[6] [gaussperm4[j]],
                             forces[15][gaussperm4[j]],
                             forces[7] [gaussperm4[j]],
                             forces[11][gaussperm4[j]],
                             forces[10][gaussperm4[j]]
                             );
 
      }
      fprintf(out,"END VALUES\n");
   }
   /*===============================shell8 element with 3x3 gausspoints */
   /* these shells have 18 stresses, do 2 x 3D Matrix */
   /* for shell8 stresses permutation: */
   /* ii[18] = {0,2,8,1,3,16,4,17,9,5,7,14,6,10,12,11,13,15};*/
   if (actgid->is_shell8_33)
   {
      ngauss=9;
      resulttype        = "MATRIX";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->shell8_33_name;
      rangetable        = actgid->standardrangetable;
      /*--- print the constant-in-thickness direction forces */
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT shell8_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cshell8_forces%c %cpcarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      fprintf(out,"COMPONENTNAMES %cForce-11%c,%cForce-22%c,%cForce-33%c,%cForce-12%c,%cForce-23%c,%cForce-13%c\n",
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign);
      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_shell8 || actele->numnp !=9) continue;
         forces = actele->e.s8->forces.a.d3[place];
         fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                             actele->Id+1,
                             forces[0][gaussperm9[0]],
                             forces[1][gaussperm9[0]],
                             forces[9][gaussperm9[0]],
                             forces[2][gaussperm9[0]],
                             forces[4][gaussperm9[0]],
                             forces[3][gaussperm9[0]]
                             );
         for (j=1; j<ngauss; j++)
         fprintf(out,"        %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                             forces[0][gaussperm9[j]],
                             forces[1][gaussperm9[j]],
                             forces[9][gaussperm9[j]],
                             forces[2][gaussperm9[j]],
                             forces[4][gaussperm9[j]],
                             forces[3][gaussperm9[j]]
                             );
 
      }
      fprintf(out,"END VALUES\n");
      /*--- print the linear-in-thickness direction forces */
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT shell8_moments on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cshell8_moments%c %cpcarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      fprintf(out,"COMPONENTNAMES %cMoment-11%c,%cMoment-22%c,%cMoment-33%c,%cMoment-12%c,%cMoment-23%c,%cMoment-13%c\n",
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign);
      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_shell8 || actele->numnp !=9) continue;
         forces = actele->e.s8->forces.a.d3[place];
         fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                             actele->Id+1,
                             forces[5] [gaussperm9[0]],
                             forces[6] [gaussperm9[0]],
                             forces[15][gaussperm9[0]],
                             forces[7] [gaussperm9[0]],
                             forces[11][gaussperm9[0]],
                             forces[10][gaussperm9[0]]
                             );
         for (j=1; j<ngauss; j++)
         fprintf(out,"        %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                             forces[5] [gaussperm9[j]],
                             forces[6] [gaussperm9[j]],
                             forces[15][gaussperm9[j]],
                             forces[7] [gaussperm9[j]],
                             forces[11][gaussperm9[j]],
                             forces[10][gaussperm9[j]]
                             );
 
      }
      fprintf(out,"END VALUES\n");
   }
#endif
#endif
#ifdef D_WALL1
/*---------------------------------------------------------fh 06/02----*/
   /*--------- ---------now go through the meshes and print the results */
   /*===============================wall1 element with 2x2 gausspoints */
   /* these walls have 4 stresses, do 1D Matrix */
   if (actgid->is_wall1_22)
   {
      ngauss=4;
      resulttype        = "MATRIX";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->wall1_22_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT wall1_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cwall1_forces%c %cpcarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-xy%c,%cStress-zz%c,%cdummy%c,%cdummy%c\n",
             sign,sign,
	     sign,sign,
	     sign,sign,
	     sign,sign,
	     sign,sign,
	     sign,sign);
      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_wall1 || actele->numnp !=4) continue;
         stress=actele->e.w1->stress_GP.a.d3[place];
	 fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             actele->Id+1,
			     stress[0][gaussperm4[0]],
			     stress[1][gaussperm4[0]],
			     stress[2][gaussperm4[0]],
			     stress[3][gaussperm4[0]]
			     );
         for (j=1; j<ngauss; j++)
         fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             stress[0][gaussperm4[j]],
			     stress[1][gaussperm4[j]],
			     stress[2][gaussperm4[j]],
			     stress[3][gaussperm4[j]]
                             );
      }
      fprintf(out,"END VALUES\n");
      
   }
   /*--------- ---------now go through the meshes and print the results */
   /*===============================wall1 element with 3x3 gausspoints */
   /* these walls have 4 stresses, do 1D Matrix */
   if (actgid->is_wall1_33)
   {
      ngauss=9;
      resulttype        = "MATRIX";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->wall1_33_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT wall1_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cwall1_forces%c %cpcarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-xy%c,%cStress-zz%c,%cdummy%c,%cdummy%c\n",
             sign,sign,
	     sign,sign,
	     sign,sign,
	     sign,sign,
	     sign,sign,
	     sign,sign);
      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_wall1 || actele->numnp < 8) continue;
         stress=actele->e.w1->stress_GP.a.d3[place];
	 fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             actele->Id+1,
			     stress[0][gaussperm9[0]],
			     stress[1][gaussperm9[0]],
			     stress[2][gaussperm9[0]],
			     stress[3][gaussperm9[0]]
			     );
         for (j=1; j<ngauss; j++)
	 fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             stress[0][gaussperm9[j]],
			     stress[1][gaussperm9[j]],
			     stress[2][gaussperm9[j]],
			     stress[3][gaussperm9[j]]
                             );
      }
      fprintf(out,"END VALUES\n");
/*----------------------------------------------------------------------*/      
   }
#endif   
#ifdef D_BRICK1   
   /* bricks have 6 stress - use 3D matrix */
   if (actgid->is_brick1_222)
   {  ngauss=8;
      resulttype        = "MATRIX";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->brick1_222_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT Brick1_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cbrick1_forces%c %cpcarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-zz%c,%cStress-xy%c,%cStress-yz%c,%cStress-xz%c\n",
              sign,sign,
              sign,sign);
      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_brick1 || actele->numnp !=8) continue;
         stress=actele->e.c1->stress_GP.a.d3[place];
	 fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             actele->Id+1,
			     stress[0][gaussperm8[0]],
			     stress[1][gaussperm8[0]],
			     stress[2][gaussperm8[0]],
			     stress[3][gaussperm8[0]],
			     stress[4][gaussperm8[0]],
			     stress[5][gaussperm8[0]]
			     );
         for (j=1; j<ngauss; j++)
         fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             stress[0][gaussperm8[j]],
			     stress[1][gaussperm8[j]],
			     stress[2][gaussperm8[j]],
			     stress[3][gaussperm8[j]],
                             stress[4][gaussperm8[j]],
			     stress[5][gaussperm8[j]]
			     );
      }
      fprintf(out,"END VALUES\n");
   }
       
   if (actgid->is_brick1_333)
   {  ngauss=27;
      resulttype        = "MATRIX";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->brick1_333_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT Brick1_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cbrick1_forces%c %cpcarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-zz%c,%cStress-xy%c,%cStress-yz%c,%cStress-xz%c\n",
              sign,sign,
              sign,sign);
      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_brick1 || actele->numnp !=20) continue;
         stress=actele->e.c1->stress_GP.a.d3[place];
	 fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             actele->Id+1,
			     stress[0][gaussperm27[0]],
			     stress[1][gaussperm27[0]],
			     stress[2][gaussperm27[0]],
			     stress[3][gaussperm27[0]],
			     stress[4][gaussperm27[0]],
			     stress[5][gaussperm27[0]]
			     );
         for (j=1; j<ngauss; j++)
         fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             stress[0][gaussperm27[j]],
			     stress[1][gaussperm27[j]],
			     stress[2][gaussperm27[j]],
			     stress[3][gaussperm27[j]],
                             stress[4][gaussperm27[j]],
			     stress[5][gaussperm27[j]]
			     );
      }
      fprintf(out,"END VALUES\n");                       
   }
#endif
#ifdef D_ALE
   if (actgid->is_ale_11)
   {
    dserror("stress output for 4-noded ale not yet impl.");
   }
   if (actgid->is_ale_22)
   {
    dserror("stress output for 4-noded ale not yet impl.");
   }
   if (actgid->is_ale_tri_1)
   {
    dserror("stress output for 3-noded ale not yet impl.");
   }
   if (actgid->is_ale_tri_3)
   {
    dserror("stress output for 3-noded ale not yet impl.");
   }
   if (actgid->is_ale_111)
   {
    dserror("stress output for 8-noded ale not yet impl.");
   }
   if (actgid->is_ale_222)
   {
    dserror("stress output for 8-noded ale not yet impl.");
   }
   if (actgid->is_ale_tet_1)
   {
    dserror("stress output for 4-noded ale not yet impl.");
   }
   if (actgid->is_ale_tet_4)
   {
    dserror("stress output for 4-noded ale not yet impl.");
   }
#endif
} /* end of (strncmp(string,"stress",stringlenght)==0) */
/*========================================= result type is velocity */
if (strncmp(string,"velocity",stringlenght)==0)
{
   resulttype        = "VECTOR";
   resultplace       = "ONNODES";
   gpset             = ""; 
   rangetable        = actgid->standardrangetable;
   ncomponent        = 3;
   componentnames[0] = "x-vel";
   componentnames[1] = "y-vel";
   componentnames[2] = "z-vel";
   /*-------------------------------------------------------------------*/
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT %s on FIELD %s\n",string,actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %c%s%c %cpcarat%c %d %s %s\n",
                                                             sign,string,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace
                                                             );
   fprintf(out,"RESULTRANGESTABLE %c%s%c\n",
                                            sign,actgid->standardrangetable,sign
                                            );
   switch (genprob.ndim)
   {
     case 3:
       fprintf(out,"COMPONENTNAMES %c%s%c,%c%s%c,%c%s%c\n",
                                                       sign,componentnames[0],sign,
                                                       sign,componentnames[1],sign,
                                                       sign,componentnames[2],sign
                                                       );
     break;
     case 2:
       fprintf(out,"COMPONENTNAMES %c%s%c,%c%s%c\n",
                                                       sign,componentnames[0],sign,
                                                       sign,componentnames[1],sign
                                                       );
     break;
     default:
       dserror("Unknown numer of dimensions");
     break;
   }
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numnp; i++)
   {
      actnode = &(actfield->dis[0].node[i]);
      switch (genprob.ndim)
      {
	case 3:
        fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                                                   actnode->Id+1,
                                                   actnode->sol.a.da[place][0],
                                                   actnode->sol.a.da[place][1],
                                                   actnode->sol.a.da[place][2]
                                                   );
	break;
	case 2:
        fprintf(out," %6d %18.5E %18.5E \n",
                                                   actnode->Id+1,
                                                   actnode->sol.a.da[place][0],
                                                   actnode->sol.a.da[place][1]
                                                   );
        break;
	default:
	  dserror("Unknown number of dimensions");
        break;
      }
   }
   fprintf(out,"END VALUES\n");
} /* end of (strncmp(string,"velocity",stringlenght)==0) */
/*========================================= result type is pressure */
if (strncmp(string,"pressure",stringlenght)==0)
{
   resulttype        = "SCALAR";
   resultplace       = "ONNODES";
   gpset             = ""; 
   rangetable        = actgid->standardrangetable;
   ncomponent        = 3;
   componentnames[0] = "pressure";
   /*-------------------------------------------------------------------*/
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT %s on FIELD %s\n",string,actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %c%s%c %cpcarat%c %d %s %s\n",
                                                             sign,string,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace
                                                             );
   fprintf(out,"RESULTRANGESTABLE %c%s%c\n",
                                            sign,actgid->standardrangetable,sign
                                            );
   switch (genprob.ndim)
   {
     case 3:
       fprintf(out,"COMPONENTNAMES %c%s%c\n",
                                                       sign,componentnames[0],sign
                                                       );
     break;
     case 2:
       fprintf(out,"COMPONENTNAMES %c%s%c\n",
                                                       sign,componentnames[0],sign
                                                       );
     break;
     default:
       dserror("Unknown numer of dimensions");
     break;
   }
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numnp; i++)
   {
      actnode = &(actfield->dis[0].node[i]);
      switch (genprob.ndim)
      {
	case 3:
        fprintf(out," %6d %18.5E \n",
                                                   actnode->Id+1,
                                                   actnode->sol.a.da[place][3]
                                                   );
	break;
	case 2:
        fprintf(out," %6d %18.5E \n",
                                                   actnode->Id+1,
                                                   actnode->sol.a.da[place][2]
                                                   );
        break;
	default:
	  dserror("Unknown number of dimensions");
        break;
      }
   }
   fprintf(out,"END VALUES\n");
} /* end of (strncmp(string,"velocity",stringlenght)==0) */
/*----------------------------------------------------------------------*/
fflush(out);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_gid_sol */




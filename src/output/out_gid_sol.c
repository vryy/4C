/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

------------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#include "gid.h"
#include "../shell8/shell8.h"
#include "../shell9/shell9.h"
#include "../wall1/wall1.h"
#include "../brick1/brick1.h"
#include "../beam3/beam3.h"
#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3.h"
#include "../fluid3/fluid3_prototypes.h"
#include "../ale2/ale2.h"
#include "../ale3/ale3.h"
#include "../axishell/axishell.h"
#include "../interf/interf.h"
#include "../wallge/wallge.h"
#include "../struct2_ml/s2ml.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       ah 03/04       |
 | vector of numfld submesh-FIELDs, defined in global_control.c         |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *sm_field;
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
#ifndef NO_TEXT_OUTPUT
INT               i,j;
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
fprintf(out,"# CCARAT postprocessing output to GID\n");
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
   actgid->is_shell8_4_22 = 0;
   actgid->is_shell8_4_33 = 0;
   actgid->is_shell8_8_22 = 0;
   actgid->is_shell8_8_33 = 0;
   actgid->is_shell8_9_22 = 0;
   actgid->is_shell8_9_33 = 0;

   actgid->is_shell9_4_22 = 0;
   actgid->is_shell9_4_33 = 0;
   actgid->is_shell9_8_22 = 0;
   actgid->is_shell9_8_33 = 0;
   actgid->is_shell9_9_22 = 0;
   actgid->is_shell9_9_33 = 0;
   actgid->is_brick1_222= 0;
   actgid->is_brick1_333= 0;
   actgid->is_fluid2_22 = 0;
   actgid->is_fluid2_33 = 0;
   actgid->is_fluid2_3 = 0;
   actgid->is_fluid2_6 = 0;
   actgid->is_fluid2_pro_22 = 0;
   actgid->is_fluid2_pro_33 = 0;
   actgid->is_fluid3_222= 0;
   actgid->is_fluid3_333= 0;

   actgid->is_f3f_8_222   = 0;
   actgid->is_f3f_8_333   = 0;
   actgid->is_f3f_20_222   = 0;
   actgid->is_f3f_20_333   = 0;

   actgid->is_ale_11    = 0;
   actgid->is_ale_22    = 0;
   actgid->is_ale_tri_1 = 0;
   actgid->is_ale_tri_3 = 0;

   actgid->is_ale_8_111   = 0;
   actgid->is_ale_8_222   = 0;
   actgid->is_ale_8_333   = 0;
   actgid->is_ale_20_111   = 0;
   actgid->is_ale_20_222   = 0;
   actgid->is_ale_20_333   = 0;

   actgid->is_wall1_11  = 0;
   actgid->is_wall1_22  = 0;
   actgid->is_wall1_33  = 0;
   actgid->is_beam3_21  = 0;
   actgid->is_beam3_22  = 0;
   actgid->is_beam3_32  = 0;
   actgid->is_beam3_33  = 0;
   actgid->is_ale_tet_1 = 0;
   actgid->is_ale_tet_4 = 0;
   actgid->is_axishell  = 0;
   actgid->is_interf_22  = 0;
   actgid->is_interf_33  = 0;
   actgid->is_wallge_22  = 0;
   actgid->is_wallge_33  = 0;




   /* check for different types of elements */
   for (j=0; j<actfield->dis[0].numele; j++)
   {
      actele = &(actfield->dis[0].element[j]);
      switch(actele->eltyp)
      {


#ifdef D_SHELL8
      case el_shell8:
         if (actele->numnp==4)
         {
           if (actele->e.s8->nGP[0]==2 && actele->e.s8->nGP[1]==2)
           {
             actgid->is_shell8_4_22   = 1;
             actgid->shell8_4_22_name = "shell8_4_22";
           }
           if (actele->e.s8->nGP[0]==3 && actele->e.s8->nGP[1]==3)
           {
             actgid->is_shell8_4_33   = 1;
             actgid->shell8_4_33_name = "shell8_4_33";
           }
         }

         if (actele->numnp==8)
         {
           if (actele->e.s8->nGP[0]==2 && actele->e.s8->nGP[1]==2)
           {
             actgid->is_shell8_8_22   = 1;
             actgid->shell8_8_22_name = "shell8_8_22";
           }
           if (actele->e.s8->nGP[0]==3 && actele->e.s8->nGP[1]==3)
           {
             actgid->is_shell8_8_33   = 1;
             actgid->shell8_8_33_name = "shell8_8_33";
           }
         }

         if (actele->numnp==9)
         {
           if (actele->e.s8->nGP[0]==2 && actele->e.s8->nGP[1]==2)
           {
             actgid->is_shell8_9_22   = 1;
             actgid->shell8_9_22_name = "shell8_9_22";
           }
           if (actele->e.s8->nGP[0]==3 && actele->e.s8->nGP[1]==3)
           {
             actgid->is_shell8_9_33   = 1;
             actgid->shell8_9_33_name = "shell8_9_33";
           }
         }
      break;
#endif /*D_SHELL8*/



#ifdef D_SHELL9
      case el_shell9:
         if (actele->numnp==4)
         {
             if (actele->e.s9->nGP[0]==2 && actele->e.s9->nGP[1]==2)
             {
                 actgid->is_shell9_4_22   = 1;
                 actgid->shell9_4_22_name = "shell9_4_22";
             }
             if (actele->e.s9->nGP[0]==3 && actele->e.s9->nGP[1]==3)
             {
                 actgid->is_shell9_4_33   = 1;
                 actgid->shell9_4_33_name = "shell9_4_33";
             }
         }
         if (actele->numnp==8)
         {
             if (actele->e.s9->nGP[0]==2 && actele->e.s9->nGP[1]==2)
             {
                 actgid->is_shell9_8_22   = 1;
                 actgid->shell9_8_22_name = "shell9_8_22";
             }
             if (actele->e.s9->nGP[0]==3 && actele->e.s9->nGP[1]==3)
             {
                 actgid->is_shell9_8_33   = 1;
                 actgid->shell9_8_33_name = "shell9_8_33";
             }
         }
         if (actele->numnp==9)
         {
             if (actele->e.s9->nGP[0]==2 && actele->e.s9->nGP[1]==2)
             {
                 actgid->is_shell9_9_22   = 1;
                 actgid->shell9_9_22_name = "shell9_9_22";
             }
             if (actele->e.s9->nGP[0]==3 && actele->e.s9->nGP[1]==3)
             {
                 actgid->is_shell9_9_33   = 1;
                 actgid->shell9_9_33_name = "shell9_9_33";
             }
         }
      break;
#endif /*D_SHELL9*/



#ifdef D_BRICK1
      case el_brick1:
         /*--- initialize stress output for hex ---*/
         c1_out_gid_sol_str(NULL, actfield, 0, 1);

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
#endif /*D_BRICK1*/



#ifdef D_FLUID2
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
#endif



#ifdef D_FLUID2_PRO
      case el_fluid2_pro:
         if (actele->numnp==4)
         {
            actgid->is_fluid2_pro_22    = 1;
            actgid->fluid2_pro_22_name  = "fluid2_pro_22";
         }
         if (actele->numnp==8  || actele->numnp==9)
         {
            actgid->is_fluid2_pro_33    = 1;
            actgid->fluid2_pro_33_name  = "fluid2_pro_33";
         }
      break;
#endif



#ifdef D_FLUID3
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
#endif



#ifdef D_FLUID3_F
      case el_fluid3_fast:
         if (actele->numnp==8)
         {
           if (actele->e.f3->nGP[0]==2 &&
               actele->e.f3->nGP[1]==2 &&
               actele->e.f3->nGP[2]==2 )
           {
             actgid->is_f3f_8_222   = 1;
             actgid->f3f_8_222_name = "f3f_8_222";
           }
           if (actele->e.f3->nGP[0]==3 &&
               actele->e.f3->nGP[1]==3 &&
               actele->e.f3->nGP[2]==3 )
           {
             actgid->is_f3f_8_333   = 1;
             actgid->f3f_8_333_name = "f3f_8_333";
           }
         }
         if (actele->numnp==20)
         {
           if (actele->e.f3->nGP[0]==2 &&
               actele->e.f3->nGP[1]==2 &&
               actele->e.f3->nGP[2]==2 )
           {
            actgid->is_f3f_20_222   = 1;
            actgid->f3f_20_222_name = "f3f_20_222";
           }
           if (actele->e.f3->nGP[0]==3 &&
               actele->e.f3->nGP[1]==3 &&
               actele->e.f3->nGP[2]==3 )
           {
            actgid->is_f3f_20_333   = 1;
            actgid->f3f_20_333_name = "f3f_20_333";
           }
         }
         if (actele->numnp==27)
           dserror("Fluid3_fast elements with 27 nodes not possible with GiD output!!");
      break;
#endif



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
        if (actele->e.ale3->nGP[0]==1 &&
            actele->e.ale3->nGP[1]==1 &&
            actele->e.ale3->nGP[2]==1 )
        {
          actgid->is_ale_8_111   = 1;
          actgid->ale_8_111_name = "ale_8_111";
        }
        if (actele->e.ale3->nGP[0]==2 &&
            actele->e.ale3->nGP[1]==2 &&
            actele->e.ale3->nGP[2]==2 )
        {
          actgid->is_ale_8_222   = 1;
          actgid->ale_8_222_name = "ale_8_222";
        }
        if (actele->e.ale3->nGP[0]==3 &&
            actele->e.ale3->nGP[1]==3 &&
            actele->e.ale3->nGP[2]==3 )
        {
          actgid->is_ale_8_333   = 1;
          actgid->ale_8_333_name = "ale_8_333";
        }
      }


      if (actele->numnp==20)
      {
        if (actele->e.ale3->nGP[0]==1 &&
            actele->e.ale3->nGP[1]==1 &&
            actele->e.ale3->nGP[2]==1 )
        {
          actgid->is_ale_20_111   = 1;
          actgid->ale_20_111_name = "ale_20_111";
        }
        if (actele->e.ale3->nGP[0]==2 &&
            actele->e.ale3->nGP[1]==2 &&
            actele->e.ale3->nGP[2]==2 )
        {
          actgid->is_ale_20_222   = 1;
          actgid->ale_20_222_name = "ale_20_222";
        }
        if (actele->e.ale3->nGP[0]==3 &&
            actele->e.ale3->nGP[1]==3 &&
            actele->e.ale3->nGP[2]==3 )
        {
          actgid->is_ale_20_333   = 1;
          actgid->ale_20_333_name = "ale_20_333";
        }
      }


      if (actele->numnp==27)
        dserror("Ale3 elements with 27 nodes not possible with GiD output");



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
#endif /*D_ALE*/



/*---------------------------------------------------------fh 06/02----*/
#ifdef D_WALL1
      case el_wall1:
         if (actele->e.w1->nGP[0]==1)
         {
            actgid->is_wall1_11    = 1;
            actgid->wall1_11_name  = "wall1_11";
         }
         if (actele->e.w1->nGP[0]==2)
         {
            actgid->is_wall1_22    = 1;
            actgid->wall1_22_name  = "wall1_22";
         }
         if (actele->e.w1->nGP[0]==3)
         {
            actgid->is_wall1_33   = 1;
            actgid->wall1_33_name = "wall1_33";
         }
      break;
#endif /*D_WALL1*/



#ifdef D_BEAM3
/*---------------------------------------------------------fh 10/02----*/
      case el_beam3:
         if (actele->numnp==2)
	 {
	    if (actele->e.b3->nGP[0]==1)
	    {
	       actgid->is_beam3_21 = 1;
	       actgid->beam3_21_name  = "beam3_21";
	    }
	    if (actele->e.b3->nGP[0]==2)
	    {
	       actgid->is_beam3_22 = 1;
	       actgid->beam3_22_name  = "beam3_22";
	    }
	 }
	 if (actele->numnp==3)
	 {
	    if (actele->e.b3->nGP[0]==2)
	    {
	       actgid->is_beam3_32 = 1;
	       actgid->beam3_32_name  = "beam3_32";
	    }
	    if (actele->e.b3->nGP[0]==3)
	    {
	       actgid->is_beam3_33 = 1;
	       actgid->beam3_33_name  = "beam3_33";
	    }
	 }
      break;
#endif /*D_BEAM3*/



#ifdef D_AXISHELL
      case el_axishell:
            actgid->is_axishell    = 1;
            actgid->axishell_name  = "axishell";
      break;
#endif /*D_AXISHELL*/



#ifdef D_INTERF
      case el_interf:
         if (actele->e.interf->nGP==2)
         {
            actgid->is_interf_22    = 1;
            actgid->interf_22_name  = "interf_22";
         }
         if (actele->e.interf->nGP==3)
         {
            actgid->is_interf_33    = 1;
            actgid->interf_33_name  = "interf_33";
         }
      break;
#endif /*D_INTERF*/



#ifdef D_WALLGE
      case el_wallge:
         if (actele->e.wallge->nGP[0]==2)
         {
            actgid->is_wallge_22    = 1;
            actgid->wallge_22_name  = "wallge_22";
         }
         if (actele->e.wallge->nGP[0]==3)
         {
            actgid->is_wallge_33    = 1;
            actgid->wallge_33_name  = "wallge_33";
         }
     break;
#endif /*D_WALLGE*/



      default:
         dserror("Unknown type of element");
      break;
      }
   } /* end of (j=0; j<actfield->numele; j++) */



   /*----------------------------- now we can write the gausspoint sets */


#ifdef D_SHELL8
#ifdef S8_HEX8  /* output of shell8 as hexahedra version, only hex8 */
   if (actgid->is_shell8_4_22 || actgid->is_shell8_4_33 ||
       actgid->is_shell8_8_22 || actgid->is_shell8_8_33 ||
       actgid->is_shell8_9_22 || actgid->is_shell8_9_33)
   {
      /* this is the shell visualization using Hexahedra */
      if (actgid->is_shell8_4_22 || actgid->is_shell8_4_33 )
      {
          if (actgid->is_shell8_4_22) /*2x2 gp*/
          {
             fprintf(out,"#-------------------------------------------------------------------------------\n");
             fprintf(out,"# GAUSSPOINTSET FOR FIELD %s SHELL8 2x2x2 GP\n",actgid->fieldname);
             fprintf(out,"#-------------------------------------------------------------------------------\n");
             fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Hexahedra %c%s%c\n",
                                                                             sign,actgid->shell8_4_22_name,sign,
                                                                             sign,actgid->shell8_4_22_name,sign);
             fprintf(out,"NUMBER OF GAUSS POINTS: 8\n");
             fprintf(out,"NATURAL COORDINATES: Internal\n");
             fprintf(out,"END GAUSSPOINTS\n");
          }
          else if (actgid->is_shell8_4_33) /*3x3 gp*/
          {
             fprintf(out,"#-------------------------------------------------------------------------------\n");
             fprintf(out,"# GAUSSPOINTSET FOR FIELD %s SHELL8 3x3 GP\n",actgid->fieldname);
             fprintf(out,"#-------------------------------------------------------------------------------\n");
             fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Hexahedra %c%s%c\n",
                                                                              sign,actgid->shell8_4_33_name,sign,
                                                                              sign,actgid->shell8_4_33_name,sign);
             fprintf(out,"NUMBER OF GAUSS POINTS: 18\n");
             fprintf(out,"NATURAL COORDINATES: Internal\n");
             fprintf(out,"END GAUSSPOINTS\n");
          }
         goto end_s8_init;
      }
      else dserror("hexahedra output for shell8 only for Quad4 !!");
   }
#endif /*S8_HEX8*/
#endif /*D_SHELL8*/

#ifdef D_SHELL8
   /* this is the shell visualization using Quadrilateral */
   if (actgid->is_shell8_4_22)
   {
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# GAUSSPOINTSET FOR FIELD %s SHELL8 2x2 GP\n",actgid->fieldname);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral %c%s%c\n",
         sign,actgid->shell8_4_22_name,sign,
         sign,actgid->shell8_4_22_name,sign);
     fprintf(out,"NUMBER OF GAUSS POINTS: 4\n");
     fprintf(out,"NATURAL COORDINATES: Internal\n");
     fprintf(out,"END GAUSSPOINTS\n");
   }

   if (actgid->is_shell8_4_33)
   {
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# GAUSSPOINTSET FOR FIELD %s SHELL8 3x3 GP\n",actgid->fieldname);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral %c%s%c\n",
         sign,actgid->shell8_4_33_name,sign,
         sign,actgid->shell8_4_33_name,sign);
     fprintf(out,"NUMBER OF GAUSS POINTS: 9\n");
     fprintf(out,"NATURAL COORDINATES: Internal\n");
     fprintf(out,"END GAUSSPOINTS\n");
   }

   if (actgid->is_shell8_8_22)
   {
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# GAUSSPOINTSET FOR FIELD %s SHELL8 2x2 GP\n",actgid->fieldname);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral %c%s%c\n",
         sign,actgid->shell8_8_22_name,sign,
         sign,actgid->shell8_8_22_name,sign);
     fprintf(out,"NUMBER OF GAUSS POINTS: 4\n");
     fprintf(out,"NATURAL COORDINATES: Internal\n");
     fprintf(out,"END GAUSSPOINTS\n");
   }

   if (actgid->is_shell8_8_33)
   {
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# GAUSSPOINTSET FOR FIELD %s SHELL8 3x3 GP\n",actgid->fieldname);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral %c%s%c\n",
         sign,actgid->shell8_8_33_name,sign,
         sign,actgid->shell8_8_33_name,sign);
     fprintf(out,"NUMBER OF GAUSS POINTS: 9\n");
     fprintf(out,"NATURAL COORDINATES: Internal\n");
     fprintf(out,"END GAUSSPOINTS\n");
   }

   if (actgid->is_shell8_9_22)
   {
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# GAUSSPOINTSET FOR FIELD %s SHELL8 2x2 GP\n",actgid->fieldname);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral %c%s%c\n",
         sign,actgid->shell8_9_22_name,sign,
         sign,actgid->shell8_9_22_name,sign);
     fprintf(out,"NUMBER OF GAUSS POINTS: 4\n");
     fprintf(out,"NATURAL COORDINATES: Internal\n");
     fprintf(out,"END GAUSSPOINTS\n");
   }

   if (actgid->is_shell8_9_33)
   {
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# GAUSSPOINTSET FOR FIELD %s SHELL8 3x3 GP\n",actgid->fieldname);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral %c%s%c\n",
         sign,actgid->shell8_9_33_name,sign,
         sign,actgid->shell8_9_33_name,sign);
     fprintf(out,"NUMBER OF GAUSS POINTS: 9\n");
     fprintf(out,"NATURAL COORDINATES: Internal\n");
     fprintf(out,"END GAUSSPOINTS\n");
   }
#endif /*D_SHELL8*/

end_s8_init:  /* end of shell8 elements */


#ifdef D_SHELL9
   /* this is the shell9 visualization using Hexahedra */
   if (actgid->is_shell9_4_22)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s SHELL9 2x2x2 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Hexahedra %c%s%c\n",
                                                                   sign,actgid->shell9_4_22_name,sign,
                                                                   sign,actgid->shell9_4_22_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 8\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_shell9_4_33)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s SHELL9 3x3x3 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Hexahedra %c%s%c\n",
                                                                   sign,actgid->shell9_4_33_name,sign,
                                                                   sign,actgid->shell9_4_33_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 27\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_shell9_8_22)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s SHELL9 2x2x2 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Hexahedra %c%s%c\n",
                                                                   sign,actgid->shell9_8_22_name,sign,
                                                                   sign,actgid->shell9_8_22_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 8\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_shell9_8_33)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s SHELL9 3x3x3 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Hexahedra %c%s%c\n",
                                                                   sign,actgid->shell9_8_33_name,sign,
                                                                   sign,actgid->shell9_8_33_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 27\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_shell9_9_22)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s SHELL9 2x2x2 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Hexahedra %c%s%c\n",
                                                                   sign,actgid->shell9_9_22_name,sign,
                                                                   sign,actgid->shell9_9_22_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 8\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_shell9_9_33)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s SHELL9 3x3x3 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Hexahedra %c%s%c\n",
                                                                   sign,actgid->shell9_9_33_name,sign,
                                                                   sign,actgid->shell9_9_33_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 27\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
#endif /*D_SHELL9*/

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
   if (actgid->is_fluid2_pro_22)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s FLUID2_PRO 2x2 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral %c%s%c\n",
                                                                sign,actgid->fluid2_pro_22_name,sign,
                                                                sign,actgid->fluid2_pro_22_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 4\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_fluid2_pro_33)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s FLUID2_PRO 3x3 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral %c%s%c\n",
                                                                sign,actgid->fluid2_pro_33_name,sign,
                                                                sign,actgid->fluid2_pro_33_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 9\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }

#ifdef D_FLUID3
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
#endif

#ifdef D_FLUID3_F
   if (actgid->is_f3f_8_222)
   {
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# GAUSSPOINTSET FOR FIELD %s FLUID3_FAST 2x2x2 GP\n",actgid->fieldname);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Hexahedra %c%s%c\n",
         sign,actgid->f3f_8_222_name,sign,
         sign,actgid->f3f_8_222_name,sign);
     fprintf(out,"NUMBER OF GAUSS POINTS: 8\n");
     fprintf(out,"NATURAL COORDINATES: Internal\n");
     fprintf(out,"END GAUSSPOINTS\n");
   }

   if (actgid->is_f3f_8_333)
   {
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# GAUSSPOINTSET FOR FIELD %s FLUID3_FAST 3x3x3 GP\n",actgid->fieldname);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Hexahedra %c%s%c\n",
         sign,actgid->f3f_8_333_name,sign,
         sign,actgid->f3f_8_333_name,sign);
     fprintf(out,"NUMBER OF GAUSS POINTS: 27\n");
     fprintf(out,"NATURAL COORDINATES: Internal\n");
     fprintf(out,"END GAUSSPOINTS\n");
   }

   if (actgid->is_f3f_20_222)
   {
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# GAUSSPOINTSET FOR FIELD %s FLUID3_FAST 2x2x2 GP\n",actgid->fieldname);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Hexahedra %c%s%c\n",
         sign,actgid->f3f_20_222_name,sign,
         sign,actgid->f3f_20_222_name,sign);
     fprintf(out,"NUMBER OF GAUSS POINTS: 8\n");
     fprintf(out,"NATURAL COORDINATES: Internal\n");
     fprintf(out,"END GAUSSPOINTS\n");
   }

   if (actgid->is_f3f_20_333)
   {
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# GAUSSPOINTSET FOR FIELD %s FLUID3_FAST 3x3x3 GP\n",actgid->fieldname);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Hexahedra %c%s%c\n",
         sign,actgid->f3f_20_333_name,sign,
         sign,actgid->f3f_20_333_name,sign);
     fprintf(out,"NUMBER OF GAUSS POINTS: 27\n");
     fprintf(out,"NATURAL COORDINATES: Internal\n");
     fprintf(out,"END GAUSSPOINTS\n");
   }
#endif


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


   if (actgid->is_ale_8_111)
   {
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# GAUSSPOINTSET FOR FIELD %s ALE 1x1x1 GP\n",actgid->fieldname);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS %c%s%cELEMTYPE Hexahedra %c%s%c\n",
         sign,actgid->ale_8_111_name,sign,
         sign,actgid->ale_8_111_name,sign);
     fprintf(out,"NUMBER OF GAUSS POINTS: 1\n");
     fprintf(out,"NATURAL COORDINATES: Internal\n");
     fprintf(out,"END GAUSSPOINTS\n");
   }

   if (actgid->is_ale_8_222)
   {
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# GAUSSPOINTSET FOR FIELD %s ALE 2x2x2 GP\n",actgid->fieldname);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS %c%s%cELEMTYPE Hexahedra %c%s%c\n",
         sign,actgid->ale_8_222_name,sign,
         sign,actgid->ale_8_222_name,sign);
     fprintf(out,"NUMBER OF GAUSS POINTS: 8\n");
     fprintf(out,"NATURAL COORDINATES: Internal\n");
     fprintf(out,"END GAUSSPOINTS\n");
   }

   if (actgid->is_ale_8_333)
   {
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# GAUSSPOINTSET FOR FIELD %s ALE 3x3x3 GP\n",actgid->fieldname);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS %c%s%cELEMTYPE Hexahedra %c%s%c\n",
         sign,actgid->ale_8_333_name,sign,
         sign,actgid->ale_8_333_name,sign);
     fprintf(out,"NUMBER OF GAUSS POINTS: 27\n");
     fprintf(out,"NATURAL COORDINATES: Internal\n");
     fprintf(out,"END GAUSSPOINTS\n");
   }

   if (actgid->is_ale_20_111)
   {
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# GAUSSPOINTSET FOR FIELD %s ALE 1x1x1 GP\n",actgid->fieldname);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS %c%s%cELEMTYPE Hexahedra %c%s%c\n",
         sign,actgid->ale_20_111_name,sign,
         sign,actgid->ale_20_111_name,sign);
     fprintf(out,"NUMBER OF GAUSS POINTS: 1\n");
     fprintf(out,"NATURAL COORDINATES: Internal\n");
     fprintf(out,"END GAUSSPOINTS\n");
   }

   if (actgid->is_ale_20_222)
   {
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# GAUSSPOINTSET FOR FIELD %s ALE 2x2x2 GP\n",actgid->fieldname);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS %c%s%cELEMTYPE Hexahedra %c%s%c\n",
         sign,actgid->ale_20_222_name,sign,
         sign,actgid->ale_20_222_name,sign);
     fprintf(out,"NUMBER OF GAUSS POINTS: 8\n");
     fprintf(out,"NATURAL COORDINATES: Internal\n");
     fprintf(out,"END GAUSSPOINTS\n");
   }

   if (actgid->is_ale_20_333)
   {
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# GAUSSPOINTSET FOR FIELD %s ALE 3x3x3 GP\n",actgid->fieldname);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"GAUSSPOINTS %c%s%cELEMTYPE Hexahedra %c%s%c\n",
         sign,actgid->ale_20_333_name,sign,
         sign,actgid->ale_20_333_name,sign);
     fprintf(out,"NUMBER OF GAUSS POINTS: 27\n");
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
   if (actgid->is_wall1_11)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s WALL1 1x1/ GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Triangle %c%s%c\n",
                                                                   sign,actgid->wall1_11_name,sign,
                                                                   sign,actgid->wall1_11_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 1\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
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
#ifdef D_BEAM3
   if (actgid->is_beam3_21)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s BEAM3 1 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Linear %c%s%c\n",
                                                             sign,actgid->beam3_21_name,sign,
							     sign,actgid->beam3_21_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 1\n");
   fprintf(out,"NODES NOT INCLUDED\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_beam3_22)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s BEAM3 2 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Linear %c%s%c\n",
                                                             sign,actgid->beam3_22_name,sign,
							     sign,actgid->beam3_22_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 2\n");
   fprintf(out,"NODES NOT INCLUDED\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_beam3_32)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s BEAM3 2 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Linear %c%s%c\n",
                                                             sign,actgid->beam3_32_name,sign,
							     sign,actgid->beam3_32_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 2\n");
   fprintf(out,"NODES NOT INCLUDED\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   if (actgid->is_beam3_33)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s BEAM3 3 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Linear %c%s%c\n",
                                                             sign,actgid->beam3_33_name,sign,
							     sign,actgid->beam3_33_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 3\n");
   fprintf(out,"NODES NOT INCLUDED\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
#endif /*D_BEAM3*/
   /* ---------------------------------------------------------------- */
   if (actgid->is_axishell)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s AXISHELL 1 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Linear\n", sign,actgid->axishell_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 1\n");
   fprintf(out,"Nodes not included\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   /* ---------------------------------------------------------------- */
   if (actgid->is_interf_22)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s INTERFACE 2(+2) GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral\n", sign,actgid->interf_22_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 4\n");
   fprintf(out,"Nodes not included\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   /* ---------------------------------------------------------------- */
   if (actgid->is_interf_33)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s INTERFACE 3(+2x3) GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral\n", sign,actgid->interf_33_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 9\n");
   fprintf(out,"Nodes not included\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   /* ---------------------------------------------------------------- */
   if (actgid->is_wallge_22)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s WALLGE 2*2 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral\n", sign,actgid->wallge_22_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 4\n");
   fprintf(out,"Nodes not included\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   /* ---------------------------------------------------------------- */
   if (actgid->is_wallge_33)
   {
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# GAUSSPOINTSET FOR FIELD %s WALLGE 3*3 GP\n",actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral\n", sign,actgid->wallge_33_name,sign);
   fprintf(out,"NUMBER OF GAUSS POINTS: 9\n");
   fprintf(out,"Nodes not included\n");
   fprintf(out,"NATURAL COORDINATES: Internal\n");
   fprintf(out,"END GAUSSPOINTS\n");
   }
   /* ---------------------------------------------------------------- */
/*----------------------------------------------------------------------*/
} /* end of (i=0; i<genprob.numfld; i++) */
/*----------------------------------------------------------------------*/
fflush(out);
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of out_gid_sol_init */



/*----------------------------------------------------------------------*
 |  write the domain numbers to GiD as gausspointvalues     m.gee 12/01 |
 *----------------------------------------------------------------------*/
void out_gid_domains(FIELD *actfield)
{
#ifndef NO_TEXT_OUTPUT
INT           i,j;

FILE         *out = allfiles.gidres;
ELEMENT      *actele;
GIDSET       *actgid = NULL;
char          sign='"';

#ifdef DEBUG
dstrc_enter("out_gid_domains");
#endif
#ifdef PERF
  perf_begin(73);
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


#ifdef D_SHELL8
#ifdef S8_HEX8  /* output of shell8 as hexahedra version, only hex8 */
    /* 4-noded shell8 element -> Hex8 */
   if (actgid->is_shell8_4_22 || actgid->is_shell8_4_33 )
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      if (actgid->is_shell8_4_22) /*2x2 gp*/
      {
         fprintf(out,"# RESULT Domains on MESH %s\n",actgid->shell8_4_22_name);
         fprintf(out,"#-------------------------------------------------------------------------------\n");
         fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                                     sign,actgid->shell8_4_22_name,sign);
      }
      else if (actgid->is_shell8_4_33) /*3x3 gp*/
      {
         fprintf(out,"# RESULT Domains on MESH %s\n",actgid->shell8_4_33_name);
         fprintf(out,"#-------------------------------------------------------------------------------\n");
         fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                                     sign,actgid->shell8_4_33_name,sign);
      }
      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_shell8 || actele->numnp != 4) continue;
         fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
         for (j=1; j<8; j++)/* hexahedra version */
         {
            fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
         }
      }
      fprintf(out,"END VALUES\n");
      goto end_8_domains;
   }
#endif /*S8_HEX8*/
#endif /*D_SHELL8*/

#ifdef D_SHELL8
if (actgid->is_shell8_4_22)
{
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"# RESULT Domains on MESH %s\n",actgid->shell8_4_22_name);
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
      sign,actgid->shell8_4_22_name,sign);
  fprintf(out,"VALUES\n");
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_shell8 || actele->numnp != 4) continue;
    fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
    for (j=1; j<4; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
  }
  fprintf(out,"END VALUES\n");
}

if (actgid->is_shell8_4_33)
{
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"# RESULT Domains on MESH %s\n",actgid->shell8_4_33_name);
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
      sign,actgid->shell8_4_33_name,sign);
  fprintf(out,"VALUES\n");
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_shell8 || actele->numnp != 4) continue;
    fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
    for (j=1; j<9; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
  }
  fprintf(out,"END VALUES\n");
}

if (actgid->is_shell8_8_22)
{
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"# RESULT Domains on MESH %s\n",actgid->shell8_8_22_name);
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
      sign,actgid->shell8_8_22_name,sign);
  fprintf(out,"VALUES\n");
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_shell8 || actele->numnp != 8) continue;
    fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
    for (j=1; j<4; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
  }
  fprintf(out,"END VALUES\n");
}

if (actgid->is_shell8_8_33)
{
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"# RESULT Domains on MESH %s\n",actgid->shell8_8_33_name);
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
      sign,actgid->shell8_8_33_name,sign);
  fprintf(out,"VALUES\n");
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_shell8 || actele->numnp != 8) continue;
    fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
    for (j=1; j<9; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
  }
  fprintf(out,"END VALUES\n");
}

if (actgid->is_shell8_9_22)
{
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"# RESULT Domains on MESH %s\n",actgid->shell8_9_22_name);
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
      sign,actgid->shell8_9_22_name,sign);
  fprintf(out,"VALUES\n");
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_shell8 || actele->numnp != 9) continue;
    fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
    for (j=1; j<4; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
  }
  fprintf(out,"END VALUES\n");
}

if (actgid->is_shell8_9_33)
{
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"# RESULT Domains on MESH %s\n",actgid->shell8_9_33_name);
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
      sign,actgid->shell8_9_33_name,sign);
  fprintf(out,"VALUES\n");
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_shell8 || actele->numnp != 9) continue;
    fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
    for (j=1; j<9; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
  }
  fprintf(out,"END VALUES\n");
}
#endif /*D_SHELL8*/

end_8_domains:  /* end of shell8 elements */

#ifdef D_SHELL9
/*4-noded*/
if (actgid->is_shell9_4_22)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->shell9_4_22_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->shell9_4_22_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_shell9 || actele->numnp != 4) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
   /* for (j=1; j<4; j++) */    /* quadrilateral version */
      for (j=1; j<8; j++)       /* hexahedra version */
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}
if (actgid->is_shell9_4_33)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->shell9_4_33_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->shell9_4_33_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_shell9 || actele->numnp != 4) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
   /* for (j=1; j<9; j++) */     /* quadrilateral version */
      for (j=1; j<27; j++)       /* hexahedra version */
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}
/*8-noded*/
if (actgid->is_shell9_8_22)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->shell9_8_22_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->shell9_8_22_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_shell9 || actele->numnp != 8) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
   /* for (j=1; j<4; j++) */     /* quadrilateral version */
      for (j=1; j<8; j++)        /* hexahedra version */
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}
if (actgid->is_shell9_8_33)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->shell9_8_33_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->shell9_8_33_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_shell9 || actele->numnp != 8) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
   /* for (j=1; j<9; j++) */    /* quadrilateral version */
      for (j=1; j<27; j++)      /* hexahedra version */
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}
/*9-noded*/
if (actgid->is_shell9_9_22)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->shell9_9_22_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->shell9_9_22_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_shell9 || actele->numnp != 9) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
   /* for (j=1; j<4; j++) */     /* quadrilateral version */
      for (j=1; j<8; j++)        /* hexahedra version */
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}
if (actgid->is_shell9_9_33)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->shell9_9_33_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->shell9_9_33_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_shell9 || actele->numnp != 9) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
   /* for (j=1; j<9; j++) */      /* quadrilateral version */
      for (j=1; j<27; j++)      /* hexahedra version */
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}
#endif /*D_SHELL9*/




if (actgid->is_brick1_222)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->brick1_222_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->brick1_222_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_brick1 || actele->numnp != 8) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
      for (j=1; j<8; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_brick1_333)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->brick1_333_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->brick1_333_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_brick1 || (actele->numnp != 20 || actele->numnp != 27)) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
      for (j=1; j<27; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_fluid3_222)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->fluid3_222_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->fluid3_222_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_fluid3 || actele->numnp != 8) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
      for (j=1; j<8; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_fluid3_333)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->fluid3_333_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->fluid3_333_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_fluid3 || (actele->numnp != 20 || actele->numnp != 27)) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
      for (j=1; j<27; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}




#ifdef D_FLUID3_F
if (actgid->is_f3f_8_222)
{
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"# RESULT Domains on MESH %s\n",actgid->f3f_8_222_name);
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
      sign,actgid->f3f_8_222_name,sign);
  fprintf(out,"VALUES\n");
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_fluid3_fast || actele->numnp != 8) continue;
    fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
    for (j=1; j<8; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
  }
  fprintf(out,"END VALUES\n");
}

if (actgid->is_f3f_8_333)
{
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"# RESULT Domains on MESH %s\n",actgid->f3f_8_333_name);
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
      sign,actgid->f3f_8_333_name,sign);
  fprintf(out,"VALUES\n");
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_fluid3_fast || actele->numnp != 8) continue;
    fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
    for (j=1; j<27; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
  }
  fprintf(out,"END VALUES\n");
}

if (actgid->is_f3f_20_222)
{
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"# RESULT Domains on MESH %s\n",actgid->f3f_20_222_name);
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
      sign,actgid->f3f_20_222_name,sign);
  fprintf(out,"VALUES\n");
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_fluid3_fast || actele->numnp != 20) continue;
    fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
    for (j=1; j<8; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
  }
  fprintf(out,"END VALUES\n");
}

if (actgid->is_f3f_20_333)
{
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"# RESULT Domains on MESH %s\n",actgid->f3f_20_333_name);
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
      sign,actgid->f3f_20_333_name,sign);
  fprintf(out,"VALUES\n");
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_fluid3_fast || actele->numnp != 20) continue;
    fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
    for (j=1; j<27; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
  }
  fprintf(out,"END VALUES\n");
}
#endif


if (actgid->is_ale_11)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_11_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                                         sign,actgid->ale_11_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_ale2 || actele->numnp != 4) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
      for (j=1; j<4; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_ale_22)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_22_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                                         sign,actgid->ale_22_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_ale2 || actele->numnp != 4) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
      for (j=1; j<4; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_ale_tri_1)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_tri_1_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                                   sign,actgid->ale_tri_1_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_ale2 || actele->numnp != 3) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
      for (j=1; j<3; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_ale_tri_3)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_tri_3_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                                   sign,actgid->ale_tri_3_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_ale2 || actele->numnp != 3) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
      for (j=1; j<3; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_ale_8_111)
{
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_8_111_name);
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
      sign,actgid->ale_8_111_name,sign);
  fprintf(out,"VALUES\n");
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_ale3 || actele->numnp != 8) continue;
    fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
  }
  fprintf(out,"END VALUES\n");
}

if (actgid->is_ale_8_222)
{
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_8_222_name);
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
      sign,actgid->ale_8_222_name,sign);
  fprintf(out,"VALUES\n");
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_ale3 || actele->numnp != 8) continue;
    fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
    for (j=1; j<8; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
  }
  fprintf(out,"END VALUES\n");
}

if (actgid->is_ale_8_333)
{
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_8_333_name);
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
      sign,actgid->ale_8_333_name,sign);
  fprintf(out,"VALUES\n");
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_ale3 || actele->numnp != 8) continue;
    fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
    for (j=1; j<27; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
  }
  fprintf(out,"END VALUES\n");
}

if (actgid->is_ale_20_111)
{
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_20_111_name);
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
      sign,actgid->ale_20_111_name,sign);
  fprintf(out,"VALUES\n");
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_ale3 || actele->numnp != 20) continue;
    fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
  }
  fprintf(out,"END VALUES\n");
}

if (actgid->is_ale_20_222)
{
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_20_222_name);
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
      sign,actgid->ale_20_222_name,sign);
  fprintf(out,"VALUES\n");
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_ale3 || actele->numnp != 20) continue;
    fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
    for (j=1; j<8; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
  }
  fprintf(out,"END VALUES\n");
}

if (actgid->is_ale_20_333)
{
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_20_333_name);
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
      sign,actgid->ale_20_333_name,sign);
  fprintf(out,"VALUES\n");
  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_ale3 || actele->numnp != 20) continue;
    fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
    for (j=1; j<27; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
  }
  fprintf(out,"END VALUES\n");
}


if (actgid->is_ale_tet_1)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_tet_1_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->ale_tet_1_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_ale3 || actele->numnp != 4) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
      for (j=1; j<4; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_ale_tet_4)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->ale_tet_4_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->ale_tet_4_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_ale3 || actele->numnp != 4) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
      for (j=1; j<4; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}


/*---------------------------------------------------------fh 06/02----*/
if (actgid->is_wall1_11)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->wall1_11_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->wall1_11_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_wall1 || actele->numnp != 3) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(double)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}

if (actgid->is_wall1_22)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->wall1_22_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->wall1_22_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_wall1) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
      for (j=1; j<4; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_wall1_33)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->wall1_33_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->wall1_33_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_wall1) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
      for (j=1; j<9; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}

#ifdef D_BEAM3
if (actgid->is_beam3_21)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->beam3_21_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->beam3_21_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_beam3) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
      for (j=1; j<1; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_beam3_22)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->beam3_22_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->beam3_22_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_beam3) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
      for (j=1; j<2; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}
if (actgid->is_beam3_32)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->beam3_32_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->beam3_32_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_beam3) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
      for (j=1; j<2; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}


if (actgid->is_beam3_33)
{
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT Domains on MESH %s\n",actgid->beam3_33_name);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cDomains%c %cpcarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,
                                                                               sign,actgid->beam3_33_name,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actfield->dis[0].numele; i++)
   {
      actele = &(actfield->dis[0].element[i]);
      if (actele->eltyp != el_beam3) continue;
      fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
      for (j=1; j<3; j++)
      fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
   }
   fprintf(out,"END VALUES\n");
}
#endif
/*----------------------------------------------------------------------*/
fflush(out);
#ifdef PERF
  perf_end(73);
#endif
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of out_gid_domains */






/*----------------------------------------------------------------------*
 |  write the domain numbers to GiD as gausspointvalues     m.gee 12/01 |
 *----------------------------------------------------------------------*/
void out_gid_domains_ssi(FIELD *actfield, INT numaf)
{
#ifndef NO_TEXT_OUTPUT
INT           i,j;

FILE         *out = allfiles.gidres;
ELEMENT      *actele;
GIDSET       *actgid = NULL;
char          sign='"';

#ifdef DEBUG
dstrc_enter("out_gid_domains");
#endif
/*----------------------------------------------------------------------*/
/*-------------------------------------- find the correct gid structure */
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
if (actgid->is_wall1_22)
{
  /*---------------------------------------------------------fh 06/02----*/
  if (numaf == 0) /*print header only for the master field */
  {
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"# RESULT Domains on MESH %s\n",actgid->wall1_22_name);
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"RESULT %cDomains%c %cccarat%c 0 SCALAR ONGAUSSPOINTS %c%s%c\n",sign,sign,sign,sign,sign,actgid->wall1_22_name,sign);
    fprintf(out,"VALUES\n");
  }

  for (i=0; i<actfield->dis[0].numele; i++)
  {
    actele = &(actfield->dis[0].element[i]);
    if (actele->eltyp != el_wall1) continue;
    fprintf(out,"    %6d  %18.5E\n",actele->Id+1,(DOUBLE)actele->proc);
    for (j=1; j<4; j++)
    fprintf(out,"            %18.5E\n",(DOUBLE)actele->proc);
  }
  if (numaf == 1) /* print end values only for the slave field */
  {
    fprintf(out,"END VALUES\n");
  }
}
else
{
  dserror("SSI output only for wall1_22 elements");
}
/*----------------------------------------------------------------------*/
fflush(out);
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of out_gid_domains */






/*----------------------------------------------------------------------*
 |  routine to write solution of a step to GID           m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_sol(char string[], FIELD *actfield, INTRA  *actintra, INT step,
                 INT place, DOUBLE time)
{
#ifndef NO_TEXT_OUTPUT
INT           i,j;

FILE         *out = allfiles.gidres;

NODE         *actnode;
ELEMENT      *actele;
GIDSET       *actgid = NULL;

char         *resulttype;
char         *resultplace;
char         *gpset;
char         *rangetable;
INT           ncomponent;
char         *componentnames[18];

char          sign='"';

INT           stringlenght;

INT           ngauss;
DOUBLE      **stress;

#ifdef D_SHELL8
DOUBLE      **forces;/* for shell8 */
DOUBLE        scal,sdc; /* for shell8 */
INT           tot_numnp;
#endif

#ifdef D_AXISHELL
DOUBLE        pv,ph,pw,px;
DOUBLE        thick;
#endif

/*
   gausspoint permutation :
   On the Gausspoint number i in Gid, the results of Carats GP number gausspermn[i]
   have to be written
*/

INT           gaussperm4[4] = {3,1,0,2};
/*INT           gaussperm8[8] = {0,4,2,6,1,5,3,7};*/
INT           gaussperm9[9] = {8,2,0,6,5,1,3,7,4};
/*INT           gaussperm27[27] = {0,9,18,3,12,21,6,15,24,1,10,19,4,13,22,7,16,25,2,11,20,5,14,23,8,17,26};*/

#ifdef DEBUG
dstrc_enter("out_gid_sol");
#endif
#ifdef PERF
  perf_begin(73);
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
   fprintf(out,"RESULT %c%s%c %cccarat%c %d %s %s\n",
                                                             sign,string,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace
                                                             );
   fprintf(out,"RESULTRANGESTABLE %c%s%c\n",
                                            sign,actgid->standardrangetable,sign
                                            );
#if 0
   /* This is a sin! Never change such a crucial global variable!
    * Everybody depends on it!
    * (It breaks the new io module in a non-obvious way!) */
#ifdef D_AXISHELL
   genprob.ndim = 3;
#endif
#endif



#ifdef D_SHELL9
   /*if shell9, the displacement have to be written in 3 components */
   /*NOTE: only first element is tested, as it is assumed that there
           is no coupling of different element types            sh 12/02 */
   if (actfield->dis[0].element[0].eltyp == el_shell9)
   {
       fprintf(out,"COMPONENTNAMES %c%s%c,%c%s%c,%c%s%c\n",
                                                       sign,componentnames[0],sign,
                                                       sign,componentnames[1],sign,
                                                       sign,componentnames[2],sign
                                                       );

       if (ioflags.struct_stress_smo == 0 ) s9_out_gid_sol_dis_unsmo(out,actfield,place);    /*unsmoothed stresses to gid*/
       if (ioflags.struct_stress_smo == 1 ) s9_out_gid_sol_dis(out,actfield,place);          /*  smoothed stresses to gid*/

    goto next;
   }
#endif /*D_SHELL9*/

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
#ifdef S8_HEX8  /* output of shell8 as hexahedra version, only hex8 */
   if (actfield->dis[0].element[0].eltyp == el_shell8 && actfield->dis[0].element[0].distyp == quad4)
   {
      tot_numnp = genprob.nnode;
      scal = 1.0;
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
#endif /*S8_HEX8*/
#endif /*D_SHELL8*/



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



#ifdef SPLIT_HEX20
   /* write solutions for additional nodes */
   for (j=0; j<actfield->dis[0].numele; j++)
   {
     INT l;
     INT eleid,nodebase;
     DOUBLE x[3];
     actele = &(actfield->dis[0].element[j]);


     if ( !(actele->eltyp == el_ale3 || actele->eltyp == el_fluid3_fast)
         || actele->numnp !=20) continue;


     eleid = actele->Id+1;
     nodebase = eleid * NODESHIFT_HEX20;

     /* node 01 */
     for (l=0;l<3;l++)
       x[l] = 0.25*(actele->node[ 8]->sol.a.da[place][l] + actele->node[ 9]->sol.a.da[place][l]
                  + actele->node[10]->sol.a.da[place][l] + actele->node[11]->sol.a.da[place][l] );
     fprintf(out," %6d %18.5e %18.5e %18.5e\n",
         nodebase + 1,
         x[0],
         x[1],
         x[2]);

     /* node 02 */
     for (l=0;l<3;l++)
       x[l] = 0.25*(actele->node[16]->sol.a.da[place][l] + actele->node[17]->sol.a.da[place][l]
                  + actele->node[18]->sol.a.da[place][l] + actele->node[19]->sol.a.da[place][l] );
     fprintf(out," %6d %18.5e %18.5e %18.5e\n",
         nodebase + 2,
         x[0],
         x[1],
         x[2]);

     /* node 03 */
     for (l=0;l<3;l++)
       x[l] = 0.25*(actele->node[ 8]->sol.a.da[place][l] + actele->node[13]->sol.a.da[place][l]
                  + actele->node[16]->sol.a.da[place][l] + actele->node[12]->sol.a.da[place][l] );
     fprintf(out," %6d %18.5e %18.5e %18.5e\n",
         nodebase + 3,
         x[0],
         x[1],
         x[2]);

     /* node 04 */
     for (l=0;l<3;l++)
       x[l] = 0.25*(actele->node[ 9]->sol.a.da[place][l] + actele->node[14]->sol.a.da[place][l]
                  + actele->node[17]->sol.a.da[place][l] + actele->node[13]->sol.a.da[place][l] );
     fprintf(out," %6d %18.5e %18.5e %18.5e\n",
         nodebase + 4,
         x[0],
         x[1],
         x[2]);

     /* node 05 */
     for (l=0;l<3;l++)
       x[l] = 0.25*(actele->node[10]->sol.a.da[place][l] + actele->node[14]->sol.a.da[place][l]
                  + actele->node[18]->sol.a.da[place][l] + actele->node[15]->sol.a.da[place][l] );
     fprintf(out," %6d %18.5e %18.5e %18.5e\n",
         nodebase + 5,
         x[0],
         x[1],
         x[2]);

     /* node 06 */
     for (l=0;l<3;l++)
       x[l] = 0.25*(actele->node[11]->sol.a.da[place][l] + actele->node[15]->sol.a.da[place][l]
                  + actele->node[19]->sol.a.da[place][l] + actele->node[12]->sol.a.da[place][l] );
     fprintf(out," %6d %18.5e %18.5e %18.5e\n",
         nodebase + 6,
         x[0],
         x[1],
         x[2]);

     /* node 07 */
     for (l=0;l<3;l++)
       x[l] = 0.0833333333
           *(actele->node[ 8]->sol.a.da[place][l] +actele->node[ 9]->sol.a.da[place][l]
           + actele->node[10]->sol.a.da[place][l] +actele->node[11]->sol.a.da[place][l]
           + actele->node[12]->sol.a.da[place][l] +actele->node[13]->sol.a.da[place][l]
           + actele->node[14]->sol.a.da[place][l] +actele->node[15]->sol.a.da[place][l]
           + actele->node[16]->sol.a.da[place][l] +actele->node[17]->sol.a.da[place][l]
           + actele->node[18]->sol.a.da[place][l] +actele->node[19]->sol.a.da[place][l] );
     fprintf(out," %6d %18.5e %18.5e %18.5e\n",
         nodebase + 7,
         x[0],
         x[1],
         x[2]);

   }  /* for (j=0; j<actfield->dis[0].numele; j++) */
#endif


   fprintf(out,"END VALUES\n");
} /* end of (strncmp(string,"displacement",stringlenght)==0) */




/*-> for shell9*/
#ifdef D_SHELL9
   next:
#endif /*D_SHELL9*/



/*============================================ result type is contact */
#ifdef D_CONTACT
if (strncmp(string,"contact",stringlenght)==0)
{
   resulttype        = "VECTOR";
   resultplace       = "ONNODES";
   gpset             = "";
   rangetable        = actgid->standardrangetable;
   ncomponent        = 3;
   componentnames[0] = "x-con";
   componentnames[1] = "y-con";
   componentnames[2] = "z-con";
   /*-------------------------------------------------------------------*/
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT %s on FIELD %s\n",string,actgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %c%s%c %cccarat%c %d %s %s\n",
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
#ifdef S8_HEX8  /* output of shell8 as hexahedra version, only hex8 */
   if (actfield->dis[0].element[0].eltyp == el_shell8 && actfield->dis[0].element[0].distyp == quad4)
   {
      tot_numnp = genprob.nnode;
      scal = 1.0;
      sdc  = actfield->dis[0].element[0].e.s8->sdc;
      for (i=0; i<actfield->dis[0].numnp; i++)
      {
         actnode = &(actfield->dis[0].node[i]);
         /* the lower surface */
         fprintf(out," %6d %23.15E %23.15E %23.15E\n",
                                                    actnode->Id+1,
                                                    -actnode->sol.a.da[place][0]+actnode->sol.a.da[place][3]*scal/sdc,
                                                    -actnode->sol.a.da[place][1]+actnode->sol.a.da[place][4]*scal/sdc,
                                                    -actnode->sol.a.da[place][2]+actnode->sol.a.da[place][5]*scal/sdc
                                                    );
         /* the upper surface */
         fprintf(out," %6d %23.15E %23.15E %23.15E\n",
                                                    actnode->Id+1+tot_numnp,
                                                    -actnode->sol.a.da[place][0]-actnode->sol.a.da[place][3]*scal/sdc,
                                                    -actnode->sol.a.da[place][1]-actnode->sol.a.da[place][4]*scal/sdc,
                                                    -actnode->sol.a.da[place][2]-actnode->sol.a.da[place][5]*scal/sdc
                                                    );
      }
   }
   else
#endif /*S8_HEX8*/ 
#endif /*D_SHELL8*/

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
} /* end of (strncmp(string,"contact",stringlenght)==0) */
#endif /*D_CONTACT*/




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
   fprintf(out,"RESULT %c%s%c %cccarat%c %d %s %s\n",
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
      fprintf(out,"RESULT %cshell8_forces%c %cccarat%c %d %s %s %c%s%c\n",
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
      fprintf(out,"RESULT %cshell8_moments%c %cccarat%c %d %s %s %c%s%c\n",
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
      fprintf(out,"RESULT %cshell8_forces%c %cccarat%c %d %s %s %c%s%c\n",
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
      fprintf(out,"RESULT %cshell8_moments%c %cccarat%c %d %s %s %c%s%c\n",
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
#endif /*D_SHELL8*/
#endif /* 0 */

#ifdef D_SHELL9
   if (actgid->is_shell9_4_22 || actgid->is_shell9_8_22 || actgid->is_shell9_9_22 ||
       actgid->is_shell9_4_33 || actgid->is_shell9_8_33 || actgid->is_shell9_9_33)
   {
      /*===============================shell9 element with 2x2 gausspoints */
      if (actgid->is_shell9_4_22 || actgid->is_shell9_8_22 || actgid->is_shell9_9_22) ngauss = 8;
      /*===============================shell9 element with 3x3 gausspoints */
      if (actgid->is_shell9_4_33 || actgid->is_shell9_8_33 || actgid->is_shell9_9_33) ngauss = 27;

      resulttype                             = "MATRIX";
      resultplace                            = "ONNODES";                 /*extrapolated to nodal points!*/
      if      (actgid->is_shell9_4_22) gpset = actgid->shell9_4_22_name;
      else if (actgid->is_shell9_8_22) gpset = actgid->shell9_8_22_name;
      else if (actgid->is_shell9_9_22) gpset = actgid->shell9_9_22_name;
      else if (actgid->is_shell9_4_33) gpset = actgid->shell9_4_33_name;
      else if (actgid->is_shell9_8_33) gpset = actgid->shell9_8_33_name;
      else if (actgid->is_shell9_9_33) gpset = actgid->shell9_9_33_name;
      rangetable        = actgid->standardrangetable;
      /*--- print the physical stresses */
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT shell9_stresses on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cshell9_stresses%c %cccarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      /*check only first element and assume, that the others are the same*/
      actele = &(actfield->dis[0].element[0]);
      switch(actele->e.s9->forcetyp)
      {
      case s9_xyz:
      fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-xy%c,%cStress-yy%c,%cStress-xz%c,%cStress-yz%c,%cStress-zz%c\n",
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign);
      break;
      case s9_rst:
      fprintf(out,"COMPONENTNAMES %cStress-rr%c,%cStress-rs%c,%cStress-ss%c,%cStress-rt%c,%cStress-st%c,%cStress-tt%c\n",
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign);
      break;
      case s9_rst_ortho:
      fprintf(out,"COMPONENTNAMES %cStress-rr%c,%cStress-rs%c,%cStress-ss%c,%cStress-rt%c,%cStress-st%c,%cStress-tt%c\n",
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign,
              sign,sign);
      break;
      }
      fprintf(out,"VALUES\n");

  /*NOTE: it is assumed that there is no coupling of different element types      sh 01/03 */

      /*write stresses at gp*/

      if (ioflags.struct_stress_smo == 0 ) s9_out_gid_sol_str_unsmo(out,actfield,place);    /*unsmoothed stresses to gid*/
      if (ioflags.struct_stress_smo == 1 ) s9_out_gid_sol_str(out,actfield,place);          /*  smoothed stresses to gid*/

      fprintf(out,"END VALUES\n");
   }
#endif /*D_SHELL9*/






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
#ifdef D_SSI
      /* if we have an ssi coupling condition on this gsurf */
      if(actfield->dis->element[0].g.gsurf->ssicouple)
      if(actfield->dis->element[0].g.gsurf->ssicouple->ssi_couptyp == ssi_master)
#endif /* D_SSI*/
      {
        fprintf(out,"#-------------------------------------------------------------------------------\n");
        fprintf(out,"# RESULT wall1_forces on FIELD %s\n",actgid->fieldname);
        fprintf(out,"#-------------------------------------------------------------------------------\n");
        fprintf(out,"RESULT %cwall1_forces%c %cccarat%c %d %s %s %c%s%c\n",
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
      }

      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
   /*---------------|| actele->numnp !=4--- */
         if (actele->eltyp != el_wall1 ) continue;
         stress =actele->e.w1->stress_GP.a.d3[place];
         if (actele->e.w1->elewa != NULL)
         {
	   fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                               actele->Id+1,
	          	     stress[0][gaussperm4[0]],
	          	     stress[1][gaussperm4[0]],
	          	     stress[2][gaussperm4[0]],
	          	     stress[3][gaussperm4[0]],
	          	     actele->e.w1->elewa[0].ipwa[gaussperm4[0]].damage,
	          	     actele->e.w1->elewa[0].ipwa[gaussperm4[0]].aequistrain
	          	     );
           for (j=1; j<ngauss; j++)
           fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                               stress[0][gaussperm4[j]],
	          	     stress[1][gaussperm4[j]],
	          	     stress[2][gaussperm4[j]],
	          	     stress[3][gaussperm4[j]],
                               actele->e.w1->elewa[0].ipwa[gaussperm4[j]].damage,
                               actele->e.w1->elewa[0].ipwa[gaussperm4[j]].aequistrain
                               );
         }
         else
         {
	   fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                               actele->Id+1,
	          	     stress[0][gaussperm4[0]],
	          	     stress[1][gaussperm4[0]],
	          	     stress[2][gaussperm4[0]],
	          	     stress[3][gaussperm4[0]],
	          	     0.0,
	          	     0.0
	          	     );
           for (j=1; j<ngauss; j++)
           fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                               stress[0][gaussperm4[j]],
	          	     stress[1][gaussperm4[j]],
	          	     stress[2][gaussperm4[j]],
	          	     stress[3][gaussperm4[j]],
                               0.0,
                               0.0
                               );
         }
      }
      fprintf(out,"END VALUES\n");

   }
/*---------------------------------------------------------he  04/03----*/
   /*------------------ now go through the meshes and print the results */
   /*================== wall1 triangle-element with 1x1 gausspoint  */
   /* these walls have 4 stresses, do 1D Matrix */
   if (actgid->is_wall1_11)
   {
      ngauss=1;
      resulttype        = "MATRIX";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->wall1_11_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT wall1_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cwall1_forces%c %cccarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-xy%c,%cStress-zz%c,%cdamage%c,%c||u||%c\n",
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
         if (actele->eltyp != el_wall1 ) continue;
         stress=actele->e.w1->stress_GP.a.d3[place];
	   fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                       actele->Id+1,
			     stress[0][0],
			     stress[1][0],
			     stress[2][0],
			     stress[3][0],
                       0.0,
                       0.0
                       );
      }
#ifdef D_SSI
      if(actfield->dis->element[0].g.gsurf->ssicouple->ssi_couptyp == ssi_slave)
#endif
        fprintf(out,"END VALUES\n");

   }
   /*-------------------now go through the meshes and print the results */
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
      fprintf(out,"RESULT %cwall1_forces%c %cccarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-xy%c,%cStress-zz%c,%cdamage%c,%caequiv_strain%c\n",
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
/*---------------|| actele->numnp < 8--- */
         if (actele->eltyp != el_wall1 ) continue;
         stress=actele->e.w1->stress_GP.a.d3[place];
         if (actele->e.w1->elewa != NULL)
         {
	   fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                               actele->Id+1,
	          	     stress[0][gaussperm9[0]],
	          	     stress[1][gaussperm9[0]],
	          	     stress[2][gaussperm9[0]],
	          	     stress[3][gaussperm9[0]],
	          	     actele->e.w1->elewa[0].ipwa[gaussperm9[0]].damage,
	          	     actele->e.w1->elewa[0].ipwa[gaussperm9[0]].aequistrain
	          	     );
           for (j=1; j<ngauss; j++)
	   fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                               stress[0][gaussperm9[j]],
	          	     stress[1][gaussperm9[j]],
	          	     stress[2][gaussperm9[j]],
	          	     stress[3][gaussperm9[j]],
	          	     actele->e.w1->elewa[0].ipwa[gaussperm9[j]].damage,
	          	     actele->e.w1->elewa[0].ipwa[gaussperm9[j]].aequistrain
                               );
         }
         else
         {
	   fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                               actele->Id+1,
	          	     stress[0][gaussperm9[0]],
	          	     stress[1][gaussperm9[0]],
	          	     stress[2][gaussperm9[0]],
	          	     stress[3][gaussperm9[0]],
	          	     0.0,
	          	     0.0
	          	     );
           for (j=1; j<ngauss; j++)
	   fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                               stress[0][gaussperm9[j]],
	          	     stress[1][gaussperm9[j]],
	          	     stress[2][gaussperm9[j]],
	          	     stress[3][gaussperm9[j]],
	          	     0.0,
	          	     0.0
                               );
         }
      }
      fprintf(out,"END VALUES\n");
/*----------------------------------------------------------------------*/
   }
#endif
#ifdef D_BRICK1
   /* bricks have 6 stress - use 3D matrix */
   if (actgid->is_brick1_222)
   {
      /*check only first element and assume, that the others are the same*/
     actele = &(actfield->dis[0].element[0]);
     switch(actele->e.c1->stresstyp)
     {
     case c1_npeqs:
      resulttype                             = "SCALAR";
      resultplace                            = "ONNODES";                 /*extrapolated to nodal points!*/
      gpset             = actgid->brick1_222_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT Brick1_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cbrick1_forces%c %cpcarat%c %d %s %s %c%s%c\n",
                   sign,sign, sign,sign, step, resulttype, resultplace, sign,gpset,sign );
        fprintf(out,"COMPONENTNAMES %cequivStress%c\n", sign,sign);
        fprintf(out,"VALUES\n");
          c1_out_gid_sol_str(out,actfield,place,0); /*extrapolated to nodal points!*/
        fprintf(out,"END VALUES\n");
     break;
     case c1_npxyz:
      resulttype                             = "MATRIX";
      resultplace                            = "ONNODES";                 /*extrapolated to nodal points!*/
      gpset             = actgid->brick1_222_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT Brick1_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cbrick1_forces%c %cpcarat%c %d %s %s %c%s%c\n",
                   sign,sign, sign,sign, step, resulttype, resultplace, sign,gpset,sign );
        fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-zz%c,%cStress-xy%c,%cStress-xz%c,%cStress-yz%c\n",
               sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);
        fprintf(out,"VALUES\n");
          c1_out_gid_sol_str(out,actfield,place,0); /*extrapolated to nodal points!*/
        fprintf(out,"END VALUES\n");
     break;
     case c1_nprst:
      resulttype                             = "MATRIX";
      resultplace                            = "ONNODES";                 /*extrapolated to nodal points!*/
      gpset             = actgid->brick1_222_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT Brick1_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cbrick1_forces%c %cpcarat%c %d %s %s %c%s%c\n",
                   sign,sign, sign,sign, step, resulttype, resultplace, sign,gpset,sign );
        fprintf(out,"COMPONENTNAMES %cStress-rr%c,%cStress-ss%c,%cStress-tt%c,%cStress-rs%c,%cStress-st%c,%cStress-tr%c\n",
               sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);
        fprintf(out,"VALUES\n");
          c1_out_gid_sol_str(out,actfield,place,0); /*extrapolated to nodal points!*/
        fprintf(out,"END VALUES\n");
     break;
     /*-----------------------------------*/
     default:
       fprintf(out,"no stresses available\n");
     break;
     }
   }
   if (actgid->is_brick1_333)
   {
      /*check only first element and assume, that the others are the same*/
     actele = &(actfield->dis[0].element[0]);
     switch(actele->e.c1->stresstyp)
     {
     case c1_npeqs:
      resulttype                             = "SCALAR";
      resultplace                            = "ONNODES";                 /*extrapolated to nodal points!*/
      gpset             = actgid->brick1_333_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT Brick1_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cbrick1_forces%c %cpcarat%c %d %s %s %c%s%c\n",
                   sign,sign, sign,sign, step, resulttype, resultplace, sign,gpset,sign );
        fprintf(out,"COMPONENTNAMES %cequivStress%c\n", sign,sign);
        fprintf(out,"VALUES\n");
          c1_out_gid_sol_str(out,actfield,place,0); /*extrapolated to nodal points!*/
        fprintf(out,"END VALUES\n");
     break;
     case c1_npxyz:
      resulttype                             = "MATRIX";
      resultplace                            = "ONNODES";                 /*extrapolated to nodal points!*/
      gpset             = actgid->brick1_333_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT Brick1_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cbrick1_forces%c %cpcarat%c %d %s %s %c%s%c\n",
                   sign,sign, sign,sign, step, resulttype, resultplace, sign,gpset,sign );
        fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-zz%c,%cStress-xy%c,%cStress-xz%c,%cStress-yz%c\n",
               sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);
        fprintf(out,"VALUES\n");
          c1_out_gid_sol_str(out,actfield,place,0); /*extrapolated to nodal points!*/
        fprintf(out,"END VALUES\n");
     break;
     case c1_nprst:
      resulttype                             = "MATRIX";
      resultplace                            = "ONNODES";                 /*extrapolated to nodal points!*/
      gpset             = actgid->brick1_333_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT Brick1_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cbrick1_forces%c %cpcarat%c %d %s %s %c%s%c\n",
                   sign,sign, sign,sign, step, resulttype, resultplace, sign,gpset,sign );
        fprintf(out,"COMPONENTNAMES %cStress-rr%c,%cStress-ss%c,%cStress-tt%c,%cStress-rs%c,%cStress-st%c,%cStress-tr%c\n",
               sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);
        fprintf(out,"VALUES\n");
          c1_out_gid_sol_str(out,actfield,place,0); /*extrapolated to nodal points!*/
        fprintf(out,"END VALUES\n");
     break;
     /*-----------------------------------*/
     default:
       fprintf(out,"no stresses available\n");
     break;
     }
   }
#endif

#ifdef D_FLUID3
   /* bricks have 6 stress - use 3D matrix */
   if (actgid->is_fluid3_222)
   {
     resulttype                             = "MATRIX";
     resultplace                            = "ONNODES";                 /*extrapolated to nodal points!*/
     gpset             = actgid->fluid3_222_name;
     rangetable        = actgid->standardrangetable;

     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# RESULT %s on FIELD %s\n",string,actgid->fieldname);
     fprintf(out,"# TIME %18.5E \n",time);
     fprintf(out,"# STEP %6d    \n",step);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"RESULT %cfluid_stresses%c %cccarat%c %d %s %s %c%s%c\n",
         sign,sign, sign,sign, step, resulttype, resultplace, sign,gpset,sign );
     fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-zz%c,%cStress-xy%c,%cStress-yz%c,%cStress-zx%c\n",
         sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);
     fprintf(out,"VALUES\n");
     f3_out_gid_sol_str(out,actfield,0,0); /*extrapolated to nodal points!*/
     fprintf(out,"END VALUES\n");

   }

   if (actgid->is_fluid3_333)
   {
     resulttype                             = "MATRIX";
     resultplace                            = "ONNODES";                 /*extrapolated to nodal points!*/
     gpset             = actgid->fluid3_333_name;
     rangetable        = actgid->standardrangetable;


     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# RESULT %s on FIELD %s\n",string,actgid->fieldname);
     fprintf(out,"# TIME %18.5E \n",time);
     fprintf(out,"# STEP %6d    \n",step);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"RESULT %cbrick1_forces%c %cccarat%c %d %s %s %c%s%c\n",
         sign,sign, sign,sign, step, resulttype, resultplace, sign,gpset,sign );
     fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-zz%c,%cStress-xy%c,%cStress-yz%c,%cStress-zx%c\n",
         sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);
     fprintf(out,"VALUES\n");
     f3_out_gid_sol_str(out,actfield,0,0); /*extrapolated to nodal points!*/
     fprintf(out,"END VALUES\n");

   }
#endif

#ifdef D_FLUID3_F
   /* bricks have 6 stress - use 3D matrix */
   if (actgid->is_f3f_8_222)
   {
     resulttype       = "MATRIX";
     resultplace      = "ONNODES";                 /*extrapolated to nodal points!*/
     gpset            = actgid->f3f_8_222_name;
     rangetable       = actgid->standardrangetable;

     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# RESULT %s on FIELD %s\n",string,actgid->fieldname);
     fprintf(out,"# TIME %18.5E \n",time);
     fprintf(out,"# STEP %6d    \n",step);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"RESULT %cfluid_stresses%c %cccarat%c %d %s %s %c%s%c\n",
         sign,sign, sign,sign, step, resulttype, resultplace, sign,gpset,sign );
     fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-zz%c,%cStress-xy%c,%cStress-yz%c,%cStress-zx%c\n",
         sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);
     fprintf(out,"VALUES\n");
     f3f_out_gid_sol_str(out,actfield,0); /*extrapolated to nodal points!*/
     fprintf(out,"END VALUES\n");
   }

   if (actgid->is_f3f_8_333)
   {
     resulttype       = "MATRIX";
     resultplace      = "ONNODES";                 /*extrapolated to nodal points!*/
     gpset            = actgid->f3f_8_333_name;
     rangetable       = actgid->standardrangetable;

     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# RESULT %s on FIELD %s\n",string,actgid->fieldname);
     fprintf(out,"# TIME %18.5E \n",time);
     fprintf(out,"# STEP %6d    \n",step);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"RESULT %cfluid_stresses%c %cccarat%c %d %s %s %c%s%c\n",
         sign,sign, sign,sign, step, resulttype, resultplace, sign,gpset,sign );
     fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-zz%c,%cStress-xy%c,%cStress-yz%c,%cStress-zx%c\n",
         sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);
     fprintf(out,"VALUES\n");
     f3f_out_gid_sol_str(out,actfield,0); /*extrapolated to nodal points!*/
     fprintf(out,"END VALUES\n");
   }

   if (actgid->is_f3f_20_222)
   {
     resulttype       = "MATRIX";
     resultplace      = "ONNODES";                 /*extrapolated to nodal points!*/
     gpset            = actgid->f3f_20_222_name;
     rangetable       = actgid->standardrangetable;

     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# RESULT %s on FIELD %s\n",string,actgid->fieldname);
     fprintf(out,"# TIME %18.5E \n",time);
     fprintf(out,"# STEP %6d    \n",step);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"RESULT %cfluid_stresses%c %cccarat%c %d %s %s %c%s%c\n",
         sign,sign, sign,sign, step, resulttype, resultplace, sign,gpset,sign );
     fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-zz%c,%cStress-xy%c,%cStress-yz%c,%cStress-zx%c\n",
         sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);
     fprintf(out,"VALUES\n");
     f3f_out_gid_sol_str(out,actfield,0); /*extrapolated to nodal points!*/
     fprintf(out,"END VALUES\n");
   }

   if (actgid->is_f3f_20_333)
   {
     resulttype       = "MATRIX";
     resultplace      = "ONNODES";                 /*extrapolated to nodal points!*/
     gpset            = actgid->f3f_20_333_name;
     rangetable       = actgid->standardrangetable;

     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# RESULT %s on FIELD %s\n",string,actgid->fieldname);
     fprintf(out,"# TIME %18.5E \n",time);
     fprintf(out,"# STEP %6d    \n",step);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"RESULT %cfluid_stresses%c %cccarat%c %d %s %s %c%s%c\n",
         sign,sign, sign,sign, step, resulttype, resultplace, sign,gpset,sign );
     fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-zz%c,%cStress-xy%c,%cStress-yz%c,%cStress-zx%c\n",
         sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);
     fprintf(out,"VALUES\n");
     f3f_out_gid_sol_str(out,actfield,0); /*extrapolated to nodal points!*/
     fprintf(out,"END VALUES\n");
   }




#endif

#ifdef D_AXISHELL
   /* STRESS output for AXISHELL */
   if (actgid->is_axishell)
   {
     /* -------------------------------- write lokale diaplcements first */
     resulttype        = "VECTOR";
     resultplace       = "ONNODES";
     gpset             = "";
     rangetable        = actgid->standardrangetable;
     ncomponent        = 3;
     componentnames[0] = "u";
     componentnames[1] = "w";
     componentnames[2] = "beta";
     /*-------------------------------------------------------------------*/
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"# RESULT %s on FIELD %s\n","verschiebungen",actgid->fieldname);
     fprintf(out,"#-------------------------------------------------------------------------------\n");
     fprintf(out,"RESULT %c%s%c %cccarat%c %d %s %s\n",
         sign,"verschiebungen",sign,
         sign,sign,
         step,
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
       fprintf(out," %6d %18.5E %18.5E %18.5E\n",
           actnode->Id+1,
           actnode->sol.a.da[1][0],
           actnode->sol.a.da[1][1],
           actnode->sol.a.da[1][2]
           );
     }
     fprintf(out,"END VALUES\n");

     /* --------------------------------------------- now write stresses */
      ngauss=1;
      resulttype        = "MATRIX";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->axishell_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT AXISHELL on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %caxishell_forces%c %cccarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );

      fprintf(out,"COMPONENTNAMES %cn_s%c,%cn_theta%c,%cm_s%c,%cm_theta%c,%cq_s%c,%cunused%c\n",
              sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);

      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_axishell || actele->numnp !=2) continue;

         /* ------------------------------------------------ */
         stress=actele->e.saxi->stress_GP.a.d3[place];
	 fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             actele->Id+1,
			     stress[0][0],
			     stress[1][0],
			     stress[2][0],
			     stress[3][0],
			     stress[4][0],
                             0.0
			     );
         /* ------------------------------------------------- */
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

   if (actgid->is_ale_8_111 ||
       actgid->is_ale_8_222 ||
       actgid->is_ale_8_333)
   {
    dserror("stress output for 8-noded ale not yet impl.");
   }
   if (actgid->is_ale_20_111 ||
       actgid->is_ale_20_222 ||
       actgid->is_ale_20_333)
   {
    dserror("stress output for 20-noded ale not yet impl.");
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
#ifdef D_BEAM3
   if (actgid->is_beam3_21)
   {
      ngauss=1;
      resulttype        = "MATRIX";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->beam3_21_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT beam3_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cbeam3_forces%c %cpcarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      fprintf(out,"COMPONENTNAMES %cN-x%c,%cV-y%c,%cV-z%c,%cM-x%c,%cM-y%c,%cM-z%c\n",
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
         if (actele->eltyp != el_beam3 ) continue;
         stress=actele->e.b3->force_GP.a.d3[place];
	 fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             actele->Id+1,
			     stress[0][0],
			     stress[1][0],
			     stress[2][0],
			     stress[3][0],
			     stress[4][0],
			     stress[5][0]
			     );
      }
      fprintf(out,"END VALUES\n");
   }

   if (actgid->is_beam3_22)
   {
      ngauss=2;
      resulttype        = "MATRIX";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->beam3_22_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT beam3_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cbeam3_forces%c %cpcarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      fprintf(out,"COMPONENTNAMES %cN-x%c,%cV-y%c,%cV-z%c,%cM-x%c,%cM-y%c,%cM-z%c\n",
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
         if (actele->eltyp != el_beam3 ) continue;
         stress=actele->e.b3->force_GP.a.d3[place];

	 fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             actele->Id+1,
			     stress[0][0],
			     stress[1][0],
			     stress[2][0],
			     stress[3][0],
			     stress[4][0],
			     stress[5][0]
			     );
	 fprintf(out,"     %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             actele->Id+1,
			     stress[0][1],
			     stress[1][1],
			     stress[2][1],
			     stress[3][1],
			     stress[4][1],
			     stress[5][1]
			     );
      }
      fprintf(out,"END VALUES\n");

   }
   if (actgid->is_beam3_32)
   {
      ngauss=2;
      resulttype        = "MATRIX";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->beam3_32_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT beam3_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cbeam3_forces%c %cpcarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      fprintf(out,"COMPONENTNAMES %cN-x%c,%cV-y%c,%cV-z%c,%cM-x%c,%cM-y%c,%cM-z%c\n",
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
         if (actele->eltyp != el_beam3 ) continue;
         stress=actele->e.b3->force_GP.a.d3[place];

	 fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             actele->Id+1,
			     stress[0][0],
			     stress[1][0],
			     stress[2][0],
			     stress[3][0],
			     stress[4][0],
			     stress[5][0]
			     );
	 fprintf(out,"     %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             actele->Id+1,
			     stress[0][1],
			     stress[1][1],
			     stress[2][1],
			     stress[3][1],
			     stress[4][1],
			     stress[5][1]
			     );
      }
      fprintf(out,"END VALUES\n");
   }

   if (actgid->is_beam3_33)
   {
      ngauss=2;
      resulttype        = "MATRIX";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->beam3_33_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT beam3_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cbeam3_forces%c %cpcarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      fprintf(out,"COMPONENTNAMES %cN-x%c,%cV-y%c,%cV-z%c,%cM-x%c,%cM-y%c,%cM-z%c\n",
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
         if (actele->eltyp != el_beam3 ) continue;
         stress=actele->e.b3->force_GP.a.d3[place];

	 fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             actele->Id+1,
			     stress[0][0],
			     stress[1][0],
			     stress[2][0],
			     stress[3][0],
			     stress[4][0],
			     stress[5][0]
			     );
	 fprintf(out,"     %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             actele->Id+1,
			     stress[0][1],
			     stress[1][1],
			     stress[2][1],
			     stress[3][1],
			     stress[4][1],
			     stress[5][1]
			     );
	 fprintf(out,"     %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             actele->Id+1,
			     stress[0][2],
			     stress[1][2],
			     stress[2][2],
			     stress[3][2],
			     stress[4][2],
			     stress[5][2]
			     );
      }
      fprintf(out,"END VALUES\n");

   }
#endif

#ifdef D_INTERF
   /* STRESS output for INTERFACE */
   if (actgid->is_interf_22)
   {
     /* ------------------------------------------------ write stresses */
      ngauss=4;
      resulttype        = "MATRIX";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->interf_22_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT INTERFACE on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cinterface_stresses%c %cccarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_interf) continue;
         if (actele->eltyp == el_interf) break;
      }
      switch(actele->e.interf->stresstyp)
      {
      /*---------------------------------- local stresses are asked for ---*/
        case if_tn:
         fprintf(out,"COMPONENTNAMES %cstress-tang%c,%cstress-normal%c,%cdummy%c,%cD-nomal%c,%cD-tang%c,%cdummy%c\n",
                 sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);
        break;
        /*--------------------------- transformation into global stresses ---*/
        case if_xy:
         fprintf(out,"COMPONENTNAMES %cstress-sxx%c,%cstress-syy%c,%cstress-sxy%c,%cD-nomal%c,%cD-tang%c,%cdummy%c\n",
              sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);
        break;
        default:
        break;
      }
      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_interf) continue;
         /* ------------------------- (|| actele->numnp !=4)  */
         /* ----------------------- gid's 1.and 4.GP get values of my 1.GP-- */
         /* ----------------------- gid's 2.and 3.GP get values of my 2.GP-- */
         stress = actele->e.interf->stress_GP.a.d3[place];
	 fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             actele->Id+1,
			     stress[0][0],
			     stress[1][0],
			     stress[2][0],
                       actele->e.interf->elewa[0].ipwa[0].dn,
                       actele->e.interf->elewa[0].ipwa[0].dt,
                             0.0
			     );
	 fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
			     stress[0][1],
			     stress[1][1],
			     stress[2][1],
                       actele->e.interf->elewa[0].ipwa[1].dn,
                       actele->e.interf->elewa[0].ipwa[1].dt,
                             0.0
			     );
	 fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
			     stress[0][1],
			     stress[1][1],
			     stress[2][1],
                       actele->e.interf->elewa[0].ipwa[1].dn,
                       actele->e.interf->elewa[0].ipwa[1].dt,
                             0.0
			     );
	 fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
			     stress[0][0],
			     stress[1][0],
			     stress[2][0],
                       actele->e.interf->elewa[0].ipwa[0].dn,
                       actele->e.interf->elewa[0].ipwa[0].dt,
                             0.0
			     );
         /* ------------------------------------------------- */
      }
      fprintf(out,"END VALUES\n");
   }
   if (actgid->is_interf_33)
   {
     /* ------------------------------------------------ write stresses */
      ngauss=9;
      resulttype        = "MATRIX";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->interf_33_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT INTERFACE on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cinterface_stresses%c %cccarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_interf) continue;
         if (actele->eltyp == el_interf) break;
      }
      switch(actele->e.interf->stresstyp)
      {
      /*---------------------------------- local stresses are asked for ---*/
        case if_tn:
         fprintf(out,"COMPONENTNAMES %cstress-tang%c,%cstress-normal%c,%cdummy%c,%cD-normal%c,%cD-tang%c,%cdummy%c\n",
                 sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);
        break;
        /*--------------------------- transformation into global stresses ---*/
        case if_xy:
         fprintf(out,"COMPONENTNAMES %cstress-sxx%c,%cstress-syy%c,%cstress-sxy%c,%cD-normal%c,%cD-tang%c,%cdummy%c\n",
              sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);
        break;
        default:
        break;
      }
      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_interf) continue;
         /* -------------------------(|| actele->numnp !=8) */
         /* ----------------------- gid's 1.and 4.and 8. GP get values of my 1.GP-- */
         /* ----------------------- gid's 5.and 7.and 9. GP get values of my 2.GP-- */
         /* ----------------------- gid's 2.and 3.and 6. GP get values of my 3.GP-- */
         stress=actele->e.interf->stress_GP.a.d3[place];
	 fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             actele->Id+1,
			     stress[0][0],
			     stress[1][0],
			     stress[2][0],
                       actele->e.interf->elewa[0].ipwa[0].dn,
                       actele->e.interf->elewa[0].ipwa[0].dt,
                             0.0
			     );
	 fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
			     stress[0][2],
			     stress[1][2],
			     stress[2][2],
                       actele->e.interf->elewa[0].ipwa[2].dn,
                       actele->e.interf->elewa[0].ipwa[2].dt,
                             0.0
			     );
	 fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
			     stress[0][2],
			     stress[1][2],
			     stress[2][2],
                       actele->e.interf->elewa[0].ipwa[2].dn,
                       actele->e.interf->elewa[0].ipwa[2].dt,
                             0.0
			     );
	 fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
			     stress[0][0],
			     stress[1][0],
			     stress[2][0],
                       actele->e.interf->elewa[0].ipwa[0].dn,
                       actele->e.interf->elewa[0].ipwa[0].dt,
                             0.0
			     );
	 fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
			     stress[0][1],
			     stress[1][1],
			     stress[2][1],
                       actele->e.interf->elewa[0].ipwa[1].dn,
                       actele->e.interf->elewa[0].ipwa[1].dt,
                             0.0
			     );
	 fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
			     stress[0][2],
			     stress[1][2],
			     stress[2][2],
                       actele->e.interf->elewa[0].ipwa[2].dn,
                       actele->e.interf->elewa[0].ipwa[2].dt,
                             0.0
			     );
	 fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
			     stress[0][1],
			     stress[1][1],
			     stress[2][1],
                       actele->e.interf->elewa[0].ipwa[1].dn,
                       actele->e.interf->elewa[0].ipwa[1].dt,
                             0.0
			     );
	 fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
			     stress[0][0],
			     stress[1][0],
			     stress[2][0],
                       actele->e.interf->elewa[0].ipwa[0].dn,
                       actele->e.interf->elewa[0].ipwa[0].dt,
                             0.0
			     );
	 fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
			     stress[0][1],
			     stress[1][1],
			     stress[2][1],
                       actele->e.interf->elewa[0].ipwa[1].dn,
                       actele->e.interf->elewa[0].ipwa[1].dt,
                             0.0
			     );
         /* ------------------------------------------------- */
      }
      fprintf(out,"END VALUES\n");
   }

#endif

#ifdef D_WALLGE
   if (actgid->is_wallge_22)
   {
      ngauss=4;
      resulttype        = "MATRIX";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->wallge_22_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT wallge_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cwallge_forces%c %cccarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-xy%c,%cdamage%c,%cloc_aequiv_strain%c,%cnonloc_aequiv_strain%c\n",
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
   /*---------------|| actele->numnp !=4--- */
         if (actele->eltyp != el_wallge ) continue;
         stress =actele->e.wallge->stress_GP.a.d3[place];
	 fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             actele->Id+1,
			     actele->e.wallge->elwa[0].iptwa[gaussperm4[0]].sig[0],
			     actele->e.wallge->elwa[0].iptwa[gaussperm4[0]].sig[1],
			     actele->e.wallge->elwa[0].iptwa[gaussperm4[0]].sig[2],
			     actele->e.wallge->elwa[0].iptwa[gaussperm4[0]].damage,
			     actele->e.wallge->elwa[0].iptwa[gaussperm4[0]].aequistrain,
			     actele->e.wallge->elwa[0].iptwa[gaussperm4[0]].aequistrain_nl
			     );
         for (j=1; j<ngauss; j++)
         fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
			     actele->e.wallge->elwa[0].iptwa[gaussperm4[j]].sig[0],
			     actele->e.wallge->elwa[0].iptwa[gaussperm4[j]].sig[1],
			     actele->e.wallge->elwa[0].iptwa[gaussperm4[j]].sig[2],
                       actele->e.wallge->elwa[0].iptwa[gaussperm4[j]].damage,
                       actele->e.wallge->elwa[0].iptwa[gaussperm4[j]].aequistrain,
			     actele->e.wallge->elwa[0].iptwa[gaussperm4[j]].aequistrain_nl
                             );
      }
      fprintf(out,"END VALUES\n");
   }
   if (actgid->is_wallge_33)
   {
      ngauss=9;
      resulttype        = "MATRIX";
      resultplace       = "ONGAUSSPOINTS";
      gpset             = actgid->wallge_33_name;
      rangetable        = actgid->standardrangetable;
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT wallge_forces on FIELD %s\n",actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cwallge_forces%c %cccarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-xy%c,%cdamage%c,%cloc_aequiv_strain%c,%cnonloc_aequiv_strain%c\n",
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
   /*---------------|| actele->numnp !=4--- */
         if (actele->eltyp != el_wallge ) continue;
         stress =actele->e.wallge->stress_GP.a.d3[place];
	 fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             actele->Id+1,
			     stress[0][gaussperm9[0]],
			     stress[1][gaussperm9[0]],
			     stress[2][gaussperm9[0]],
			     actele->e.wallge->elwa[0].iptwa[gaussperm9[0]].damage,
			     actele->e.wallge->elwa[0].iptwa[gaussperm9[0]].aequistrain,
			     actele->e.wallge->elwa[0].iptwa[gaussperm9[0]].aequistrain_nl
			     );
         for (j=1; j<ngauss; j++)
         fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             stress[0][gaussperm9[j]],
			     stress[1][gaussperm9[j]],
			     stress[2][gaussperm9[j]],
                       actele->e.wallge->elwa[0].iptwa[gaussperm9[j]].damage,
                       actele->e.wallge->elwa[0].iptwa[gaussperm9[j]].aequistrain,
			     actele->e.wallge->elwa[0].iptwa[gaussperm9[j]].aequistrain_nl
                             );
      }
      fprintf(out,"END VALUES\n");
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
   fprintf(out,"# TIME %18.5E \n",time);
   fprintf(out,"# STEP %6d    \n",step);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %c%s%c %cccarat%c %d %s %s\n",
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
        fprintf(out," %6d %22.9E %22.9E %22.9E\n",
                                                   actnode->Id+1,
                                                   actnode->sol.a.da[place][0],
                                                   actnode->sol.a.da[place][1],
                                                   actnode->sol.a.da[place][2]
                                                   );
	break;
	case 2:
        fprintf(out," %6d %22.9E %22.9E \n",
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


#ifdef SPLIT_HEX20
   /* write solutions for additional nodes */
   for (j=0; j<actfield->dis[0].numele; j++)
   {
     INT l;
     INT eleid,nodebase;
     DOUBLE x[3];
     actele = &(actfield->dis[0].element[j]);


     if ( !(actele->eltyp == el_ale3 || actele->eltyp == el_fluid3_fast)
         || actele->numnp !=20) continue;


     eleid = actele->Id+1;
     nodebase = eleid * NODESHIFT_HEX20;

     /* node 01 */
     for (l=0;l<3;l++)
       x[l] = 0.25*(actele->node[ 8]->sol.a.da[place][l] + actele->node[ 9]->sol.a.da[place][l]
                  + actele->node[10]->sol.a.da[place][l] + actele->node[11]->sol.a.da[place][l] );
     fprintf(out," %6d %18.5e %18.5e %18.5e\n",
         nodebase + 1,
         x[0],
         x[1],
         x[2]);

     /* node 02 */
     for (l=0;l<3;l++)
       x[l] = 0.25*(actele->node[16]->sol.a.da[place][l] + actele->node[17]->sol.a.da[place][l]
                  + actele->node[18]->sol.a.da[place][l] + actele->node[19]->sol.a.da[place][l] );
     fprintf(out," %6d %18.5e %18.5e %18.5e\n",
         nodebase + 2,
         x[0],
         x[1],
         x[2]);

     /* node 03 */
     for (l=0;l<3;l++)
       x[l] = 0.25*(actele->node[ 8]->sol.a.da[place][l] + actele->node[13]->sol.a.da[place][l]
                  + actele->node[16]->sol.a.da[place][l] + actele->node[12]->sol.a.da[place][l] );
     fprintf(out," %6d %18.5e %18.5e %18.5e\n",
         nodebase + 3,
         x[0],
         x[1],
         x[2]);

     /* node 04 */
     for (l=0;l<3;l++)
       x[l] = 0.25*(actele->node[ 9]->sol.a.da[place][l] + actele->node[14]->sol.a.da[place][l]
                  + actele->node[17]->sol.a.da[place][l] + actele->node[13]->sol.a.da[place][l] );
     fprintf(out," %6d %18.5e %18.5e %18.5e\n",
         nodebase + 4,
         x[0],
         x[1],
         x[2]);

     /* node 05 */
     for (l=0;l<3;l++)
       x[l] = 0.25*(actele->node[10]->sol.a.da[place][l] + actele->node[14]->sol.a.da[place][l]
                  + actele->node[18]->sol.a.da[place][l] + actele->node[15]->sol.a.da[place][l] );
     fprintf(out," %6d %18.5e %18.5e %18.5e\n",
         nodebase + 5,
         x[0],
         x[1],
         x[2]);

     /* node 06 */
     for (l=0;l<3;l++)
       x[l] = 0.25*(actele->node[11]->sol.a.da[place][l] + actele->node[15]->sol.a.da[place][l]
                  + actele->node[19]->sol.a.da[place][l] + actele->node[12]->sol.a.da[place][l] );
     fprintf(out," %6d %18.5e %18.5e %18.5e\n",
         nodebase + 6,
         x[0],
         x[1],
         x[2]);

     /* node 07 */
     for (l=0;l<3;l++)
       x[l] = 0.0833333333
           *(actele->node[ 8]->sol.a.da[place][l] +actele->node[ 9]->sol.a.da[place][l]
           + actele->node[10]->sol.a.da[place][l] +actele->node[11]->sol.a.da[place][l]
           + actele->node[12]->sol.a.da[place][l] +actele->node[13]->sol.a.da[place][l]
           + actele->node[14]->sol.a.da[place][l] +actele->node[15]->sol.a.da[place][l]
           + actele->node[16]->sol.a.da[place][l] +actele->node[17]->sol.a.da[place][l]
           + actele->node[18]->sol.a.da[place][l] +actele->node[19]->sol.a.da[place][l] );
     fprintf(out," %6d %18.5e %18.5e %18.5e\n",
         nodebase + 7,
         x[0],
         x[1],
         x[2]);

   }  /* for (j=0; j<actfield->dis[0].numele; j++) */
#endif


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
   fprintf(out,"# TIME %18.5E \n",time);
   fprintf(out,"# STEP %6d    \n",step);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %c%s%c %cccarat%c %d %s %s\n",
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
        fprintf(out," %6d %22.9E \n",
                                                   actnode->Id+1,
                                                   actnode->sol.a.da[place][3]
                                                   );
	break;
	case 2:
        fprintf(out," %6d %22.9E \n",
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
} /* end of (strncmp(string,"pressure",stringlenght)==0) */




/*========================================= result type is thickness */
#ifdef D_AXISHELL
if (strncmp(string,"thickness",stringlenght)==0)
{
  if (actgid->is_axishell)
  {
    resulttype        = "SCALAR";
    resultplace       = "ONNODES";
    gpset             = "";
    rangetable        = actgid->standardrangetable;
    ncomponent        = 1;
    componentnames[0] = "thickness";
    /*-------------------------------------------------------------------*/
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"# RESULT %s on FIELD %s\n",string,actgid->fieldname);
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"RESULT %c%s%c %cccarat%c %d %s %s\n",
        sign,string,sign,
        sign,sign,
        step,
        resulttype,
        resultplace
        );

    fprintf(out,"RESULTRANGESTABLE %c%s%c\n", sign,actgid->standardrangetable,sign);

    fprintf(out,"COMPONENTNAMES %c%s%c\n",sign,componentnames[0],sign);

    fprintf(out,"VALUES\n");

    for (i=0; i<actfield->dis[0].numnp; i++)
    {
      actnode = &(actfield->dis[0].node[i]);
      actele = actnode->element[0];

      if (actele==NULL)
        dserror("Error in output of thickness");

      if (actele->node[0] == actnode)
        thick = actele->e.saxi->thick[0];
      else
        thick = actele->e.saxi->thick[1];

      fprintf(out," %6d %18.5E \n", actnode->Id+1, thick);
    }
    fprintf(out,"END VALUES\n");
  }
} /* end of (strncmp(string,"thickness",stringlenght)==0) */
#endif




/*========================================= result type is axi-loads */
#ifdef D_AXISHELL
if (strncmp(string,"axi_loads",stringlenght)==0)
{
  if (actgid->is_axishell)
  {
    resulttype        = "MATRIX";
    resultplace       = "ONNODES";
    gpset             = "";
    rangetable        = actgid->standardrangetable;
    ncomponent        = 4;
    componentnames[0] = "pv";
    componentnames[1] = "ph";
    componentnames[2] = "px";
    componentnames[3] = "pw";
    componentnames[4] = "not_used";
    componentnames[5] = "not_used";
    /*-------------------------------------------------------------------*/
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"# RESULT %s on FIELD %s\n",string,actgid->fieldname);
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"RESULT %c%s%c %cccarat%c %d %s %s\n",
        sign,string,sign,
        sign,sign,
        step,
        resulttype,
        resultplace
        );

    fprintf(out,"RESULTRANGESTABLE %c%s%c\n", sign,actgid->standardrangetable,sign);

    fprintf(out,"COMPONENTNAMES %c%s%c %c%s%c %c%s%c %c%s%c %c%s%c %c%s%c\n",
        sign,componentnames[0],sign,
        sign,componentnames[1],sign,
        sign,componentnames[2],sign,
        sign,componentnames[3],sign,
        sign,componentnames[4],sign,
        sign,componentnames[5],sign);

    fprintf(out,"VALUES\n");

    for (i=0; i<actfield->dis[0].numnp; i++)
    {
      actnode = &(actfield->dis[0].node[i]);
      actele = actnode->element[0];

      if (actele==NULL)
        dserror("Error in output of axishell loads");

      if (actele->node[0] == actnode)
      {
        pv = actele->e.saxi->pv[0];
        ph = actele->e.saxi->ph[0];
        px = actele->e.saxi->px[0];
        pw = actele->e.saxi->pw[0];
      }
      else
      {
        pv = actele->e.saxi->pv[1];
        ph = actele->e.saxi->ph[1];
        px = actele->e.saxi->px[1];
        pw = actele->e.saxi->pw[1];
      }

      fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
          actnode->Id+1,pv,ph,px,pw,0.0,0.0);
    }
    fprintf(out,"END VALUES\n");
  }
} /* end of (strncmp(string,"axi_loads",stringlenght)==0) */
#endif
/*----------------------------------------------------------------------*/
fflush(out);
#ifdef PERF
  perf_end(73);
#endif
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of out_gid_sol */







#ifdef D_MLSTRUCT
/*----------------------------------------------------------------------*
 |  init output of multiscale-submesh-results for GID       ah 06/04    |
 *----------------------------------------------------------------------*/
void out_gid_smsol_init()
{
#ifndef NO_TEXT_OUTPUT
INT               j;
FIELD            *actsmfield;
GIDSET           *actsmgid;
ELEMENT          *actsmele;
FILE             *out       = allfiles.gidmicrores;
char             *charptr;
char              sign='"';

#ifdef DEBUG
dstrc_enter("out_gid_smsol_init");
#endif
/*----------------------------------------------------------------------*/
fprintf(out,"Gid Post Results File 1.0\n");
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"# CCARAT postprocessing output to GID for submesh-solutions\n");
fprintf(out,"# For postprocessing with gid, this file has to be copied to 'newprojectname.flavia.res'\n");
fprintf(out,"#-------------------------------------------------------------------------------\n");
/*----------------------------------------------------------------------*/
/*------------------------- check number of fields and allocate storage */
sm_gid = (GIDSET*)CCACALLOC(1,sizeof(GIDSET));
if (!sm_gid) dserror("Allocation of memory failed");
/*------------------------------------------------------- init the data */
actsmfield         = &(sm_field[0]);
actsmgid           = &(sm_gid[0]);
actsmgid->fieldtyp = structure;
/*--------------------------------------------------- set the fieldname */
actsmgid->fieldnamelenght = 9;
actsmgid->fieldname       = "structure";
/*---------------------------------------------------- set range tables */
strncpy(actsmgid->standardrangetable,"standard_         ",18);
charptr = actsmgid->standardrangetable + 9;
strncpy(charptr,actsmgid->fieldname,actsmgid->fieldnamelenght);
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"# RANGETABLES %s\n",actsmgid->fieldname);
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"RESULTRANGESTABLE %c%s%c\n",sign,actsmgid->standardrangetable,sign);
fprintf(out,"            - -1000000.0 : %cvery small%c\n",sign,sign);
fprintf(out," -1000000.0 -  1000000.0 : %cnormal%c\n",sign,sign);
fprintf(out,"  1000000.0 -            : %cvery large%c\n",sign,sign);
fprintf(out,"END RESULTRANGESTABLE\n");
/*--------------------------- find and set meshes and gausspointsets */
actsmgid->is_shell8_4_22 = 0;
actsmgid->is_shell8_4_33 = 0;
actsmgid->is_shell8_8_22 = 0;
actsmgid->is_shell8_8_33 = 0;
actsmgid->is_shell8_9_22 = 0;
actsmgid->is_shell8_9_33 = 0;

actsmgid->is_shell9_4_22 = 0;
actsmgid->is_shell9_4_33 = 0;
actsmgid->is_shell9_8_22 = 0;
actsmgid->is_shell9_8_33 = 0;
actsmgid->is_shell9_9_22 = 0;
actsmgid->is_shell9_9_33 = 0;
actsmgid->is_brick1_222  = 0;
actsmgid->is_brick1_333  = 0;
actsmgid->is_fluid2_22   = 0;
actsmgid->is_fluid2_33   = 0;
actsmgid->is_fluid2_3    = 0;
actsmgid->is_fluid2_6    = 0;
actsmgid->is_fluid2_pro_22 = 0;
actsmgid->is_fluid2_pro_33 = 0;
actsmgid->is_fluid3_222  = 0;
actsmgid->is_fluid3_333  = 0;

actsmgid->is_f3f_8_222   = 0;
actsmgid->is_f3f_8_333   = 0;
actsmgid->is_f3f_20_222  = 0;
actsmgid->is_f3f_20_333  = 0;

actsmgid->is_ale_11      = 0;
actsmgid->is_ale_22      = 0;
actsmgid->is_ale_tri_1   = 0;
actsmgid->is_ale_tri_3   = 0;

actsmgid->is_ale_8_111   = 0;
actsmgid->is_ale_8_222   = 0;
actsmgid->is_ale_8_333   = 0;
actsmgid->is_ale_20_111  = 0;
actsmgid->is_ale_20_222  = 0;
actsmgid->is_ale_20_333  = 0;

actsmgid->is_wall1_22  = 0;
actsmgid->is_wall1_33  = 0;
actsmgid->is_beam3_21  = 0;
actsmgid->is_beam3_22  = 0;
actsmgid->is_beam3_32  = 0;
actsmgid->is_beam3_33  = 0;
actsmgid->is_ale_tet_1 = 0;
actsmgid->is_ale_tet_4 = 0;
actsmgid->is_axishell  = 0;
actsmgid->is_interf_22  = 0;
actsmgid->is_interf_33  = 0;
actsmgid->is_wallge_22  = 0;
actsmgid->is_wallge_33  = 0;
/*---------------------------- check for different types of elements */
for (j=0; j<actsmfield->dis[0].numele; j++)
{
   actsmele = &(actsmfield->dis[0].element[j]);
   switch(actsmele->eltyp)
   {
#ifdef D_WALL1
   case el_wall1:
      if (actsmele->e.w1->nGP[0]==2)
      {
         actsmgid->is_wall1_22    = 1;
         actsmgid->wall1_22_name  = "wall1_22";
      }
      if (actsmele->e.w1->nGP[0]==3)
      {
         actsmgid->is_wall1_33   = 1;
         actsmgid->wall1_33_name = "wall1_33";
      }
   break;
#endif /*D_WALL1*/
#ifdef D_INTERF
   case el_interf:
      if (actsmele->e.interf->nGP==2)
      {
         actsmgid->is_interf_22    = 1;
         actsmgid->interf_22_name  = "interf_22";
      }
      if (actsmele->e.interf->nGP==3)
      {
         actsmgid->is_interf_33    = 1;
         actsmgid->interf_33_name  = "interf_33";
      }
   break;
#endif /*D_INTERF*/
   default:
      dserror("Unknown type of element");
   break;
   }
} /* end of (j=0; j<actamfield->numele; j++) */
/*----------------------------- now we can write the gausspoint sets */
if (actsmgid->is_wall1_22)
{
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"# GAUSSPOINTSET FOR FIELD %s WALL1 2x2 GP\n",actsmgid->fieldname);
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral %c%s%c\n",
                                                                sign,actsmgid->wall1_22_name,sign,
                                                                sign,actsmgid->wall1_22_name,sign);
fprintf(out,"NUMBER OF GAUSS POINTS: 4\n");
fprintf(out,"NATURAL COORDINATES: Internal\n");
fprintf(out,"END GAUSSPOINTS\n");
}
if (actsmgid->is_wall1_33)
{
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"# GAUSSPOINTSET FOR FIELD %s WALL1 3x3 GP\n",actsmgid->fieldname);
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral %c%s%c\n",
                                                                 sign,actsmgid->wall1_33_name,sign,
                                                                 sign,actsmgid->wall1_33_name,sign);
fprintf(out,"NUMBER OF GAUSS POINTS: 9\n");
fprintf(out,"NATURAL COORDINATES: Internal\n");
fprintf(out,"END GAUSSPOINTS\n");
}
/* ---------------------------------------------------------------- */
if (actsmgid->is_interf_22)
{
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"# GAUSSPOINTSET FOR FIELD %s INTERFACE 2(+2) GP\n",actsmgid->fieldname);
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral\n", sign,actsmgid->interf_22_name,sign);
fprintf(out,"NUMBER OF GAUSS POINTS: 4\n");
fprintf(out,"Nodes not included\n");
fprintf(out,"NATURAL COORDINATES: Internal\n");
fprintf(out,"END GAUSSPOINTS\n");
}
/* ---------------------------------------------------------------- */
if (actsmgid->is_interf_33)
{
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"# GAUSSPOINTSET FOR FIELD %s INTERFACE 3(+2x3) GP\n",actsmgid->fieldname);
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"GAUSSPOINTS %c%s%c ELEMTYPE Quadrilateral\n", sign,actsmgid->interf_33_name,sign);
fprintf(out,"NUMBER OF GAUSS POINTS: 9\n");
fprintf(out,"Nodes not included\n");
fprintf(out,"NATURAL COORDINATES: Internal\n");
fprintf(out,"END GAUSSPOINTS\n");
}
/*----------------------------------------------------------------------*/
fflush(out);
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of out_gid_smsol_init */
#endif /* D_MLSTRUCT */



#ifdef D_MLSTRUCT
/*----------------------------------------------------------------------*
 |  routine to write submesh small-scale displacements to GID  ah 8/04  |
 *----------------------------------------------------------------------*/
void out_gid_smdisp(char string[], INT step)
{
#ifndef NO_TEXT_OUTPUT
INT           i,j;
FILE         *out = allfiles.gidmicrores;
FIELD        *actsmfield;     /* the submesh-field */
FIELD        *actmafield;     /* the "macro"-field */
GIDSET       *actsmgid;
ELEMENT      *actmaele;       /* the actual macro-element */
NODE         *actsmnode;      /* the actual submesh-node */
INT           GlobalID;       /* pseudo global element number of submesh elements */

char          sign='"';
char         *resulttype;
char         *resultplace;
char         *gpset;
char         *rangetable;
INT           ncomponent;
char         *componentnames[2];

#ifdef DEBUG
dstrc_enter("out_gid_smdisp");
#endif
/*----------------------------------------------------------------------*/
actsmfield = &(sm_field[0]);
actmafield = &(field[0]);
actsmgid   = &(sm_gid[0]);
/*----------------------------------------------------------------------*/
resulttype        = "VECTOR";
resultplace       = "ONNODES";
gpset             = "";
rangetable        = actsmgid->standardrangetable;
ncomponent        = 2;
componentnames[0] = "x-displ";
componentnames[1] = "y-displ";
/*-----------------------------------------------------------------------*/
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"# RESULT %s on FIELD %s\n",string,actsmgid->fieldname);
fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"RESULT %c%s%c %cccarat%c %d %s %s\n",
                                                          sign,string,sign,
                                                          sign,sign,
                                                          step,
                                                          resulttype,
                                                          resultplace
                                                          );
fprintf(out,"RESULTRANGESTABLE %c%s%c\n",
                                         sign,actsmgid->standardrangetable,sign
                                         );
fprintf(out,"COMPONENTNAMES %c%s%c,%c%s%c\n",
                                                sign,componentnames[0],sign,
                                                sign,componentnames[1],sign
                                                );
fprintf(out,"VALUES\n");
/*-------------------------- print the submesh fine scale displacements */
/*--------------------------------------------- loop over macroelements */
for (i=0; i<actmafield->dis[0].numele; i++)
{
  actmaele = &(actmafield->dis[0].element[i]);
  if(actmaele->e.w1->isinomegaprime == 1)
  {
  /*-------------------------- write NodeID and Nodecoordinates to file */
    for (j=0; j<actsmfield->dis[0].numnp; j++)
    {
       actsmnode = &(actsmfield->dis[0].node[j]);
       GlobalID  = i * actsmfield->dis[0].numnp + (actsmnode->Id + 1);
       fprintf(out," %6d %18.5E %18.5E \n",
                                           GlobalID,
                                           actmaele->e.w1->sm_nodaldata[actsmnode->Id].displ_mi[0],
                                           actmaele->e.w1->sm_nodaldata[actsmnode->Id].displ_mi[1]
                                            );
    }
  }
}
fprintf(out,"END VALUES\n");
/*----------------------------------------------------------------------*/
fflush(out);
/*-------------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of out_gid_smdisp */
#endif /* D_MLSTRUCT */


#ifdef D_MLSTRUCT
/*----------------------------------------------------------------------*
 |  routine to write submesh stressses to GID                  ah 8/04  |
 *----------------------------------------------------------------------*/
void out_gid_smstress(char string[], INT step)
{
#ifndef NO_TEXT_OUTPUT
INT           i,j,a,b;

FILE         *out = allfiles.gidmicrores;
FIELD        *actsmfield;     /* the submesh-field */
FIELD        *actmafield;     /* the "macro"-field */
GIDSET       *actsmgid;
ELEMENT      *actsmele;       /* the actual submesh-element */
ELEMENT      *actmaele;       /* the actual macro-element */
INT           GlobalID;       /* pseudo global element number of submesh elements */
DOUBLE        stress[9][4];
DOUBLE        damage[9];

char          sign='"';
char         *resulttype;
char         *resultplace;
char         *gpset;
char         *rangetable;
INT           ngauss;
/*
   gausspoint permutation :
   On the Gausspoint number i in Gid, the results of Carats GP number gausspermn[i]
   have to be written
*/
INT           gaussperm4[4] = {3,1,0,2};
INT           gaussperm9[9] = {8,2,0,6,5,1,3,7,4};

#ifdef DEBUG
dstrc_enter("out_gid_smstress");
#endif
/*----------------------------------------------------------------------*/
actsmfield = &(sm_field[0]);
actmafield = &(field[0]);
actsmgid   = &(sm_gid[0]);
/*----------------------------------------------------------------------*/
/*------------------------------------------ print the submesh stresses */
if (actsmgid->is_wall1_22)
{
   ngauss=4;
   resulttype        = "MATRIX";
   resultplace       = "ONGAUSSPOINTS";
   gpset             = actsmgid->wall1_22_name;
   rangetable        = actsmgid->standardrangetable;
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT wall1_stresses in the submesh on FIELD %s\n",actsmgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cwall1_forces%c %cccarat%c %d %s %s %c%s%c\n",
                                                          sign,sign,
                                                          sign,sign,
                                                          step,
                                                          resulttype,
                                                          resultplace,
                                                          sign,gpset,sign
                                                          );
   fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-xy%c,%cStress-zz%c,%cdamage%c,%cdummy%c\n",
          sign,sign,
        sign,sign,
        sign,sign,
        sign,sign,
        sign,sign,
        sign,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actmafield->dis[0].numele; i++)
   {
     actmaele = &(actmafield->dis[0].element[i]);
     if(actmaele->e.w1->isinomegaprime == 1)
     {
       for (j=0; j<actsmfield->dis[0].numele; j++)
       {
         actsmele = &(actsmfield->dis[0].element[j]);
         if (actsmele->eltyp != el_wall1 ) continue;
         GlobalID = i * actsmfield->dis[0].numele + (actsmele->Id + 1);

         for (a=0; a<ngauss; a++)
         {
           damage[a] = actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[a].damage;
           for (b=0; b<4; b++)
              stress[a][b] = actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[a].stress[b];
          }
	    fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                   GlobalID,
	             stress[gaussperm4[0]][0],
	             stress[gaussperm4[0]][1],
	             stress[gaussperm4[0]][2],
	             stress[gaussperm4[0]][3],
	             damage[gaussperm4[0]],
	             0.0
	             );
           for (a=1; a<ngauss; a++)
             fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                     stress[gaussperm4[a]][0],
	               stress[gaussperm4[a]][1],
	               stress[gaussperm4[a]][2],
	               stress[gaussperm4[a]][3],
                     damage[gaussperm4[a]],
                     0.0
                     );
       }
     }
   }/* end: for (i=0; i<actmafield->dis[0].numele; i++) */
   fprintf(out,"END VALUES\n");
} /* end: if (actsmgid->is_wall1_22) */
/*----------------------------------------------------------------------*/
if (actsmgid->is_wall1_33)
{
   ngauss=9;
   resulttype        = "MATRIX";
   resultplace       = "ONGAUSSPOINTS";
   gpset             = actsmgid->wall1_33_name;
   rangetable        = actsmgid->standardrangetable;
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT wall1_forces on FIELD %s\n",actsmgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cwall1_forces%c %cccarat%c %d %s %s %c%s%c\n",
                                                          sign,sign,
                                                          sign,sign,
                                                          step,
                                                          resulttype,
                                                          resultplace,
                                                          sign,gpset,sign
                                                          );
   fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-xy%c,%cStress-zz%c,%cdamage%c,%cdummy%c\n",
          sign,sign,
          sign,sign,
          sign,sign,
          sign,sign,
          sign,sign,
          sign,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actmafield->dis[0].numele; i++)
   {
     actmaele = &(actmafield->dis[0].element[i]);
     if(actmaele->e.w1->isinomegaprime == 1)
     {
       for (j=0; j<actsmfield->dis[0].numele; j++)
       {
         actsmele = &(actsmfield->dis[0].element[j]);
         if (actsmele->eltyp != el_wall1 ) continue;
         GlobalID = i * actsmfield->dis[0].numele + (actsmele->Id + 1);

         for (a=0; a<ngauss; a++)
         {
           damage[a] = actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[a].damage;
           for (b=0; b<4; b++)
              stress[a][b] = actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[a].stress[b];
          }
          fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                     GlobalID,
                     stress[gaussperm9[0]][0],
                     stress[gaussperm9[0]][1],
                     stress[gaussperm9[0]][2],
                     stress[gaussperm9[0]][3],
                     damage[gaussperm9[0]],
                     0.0
                     );
          for (a=1; a<ngauss; a++)
             fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                     stress[gaussperm9[a]][0],
                     stress[gaussperm9[a]][1],
                     stress[gaussperm9[a]][2],
                     stress[gaussperm9[a]][3],
                     damage[gaussperm9[a]],
                     0.0
                     );
       }
     }
   }/* end: for (i=0; i<actmafield->dis[0].numele; i++) */
   fprintf(out,"END VALUES\n");
} /* end: if (actsmgid->is_wall1_33) */
/*----------------------------------------------------------------------*/
if (actsmgid->is_interf_22)
{
   resulttype        = "MATRIX";
   resultplace       = "ONGAUSSPOINTS";
   gpset             = actsmgid->interf_22_name;
   rangetable        = actsmgid->standardrangetable;
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT INTERFACE on FIELD %s\n",actsmgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cinterface_stresses%c %cccarat%c %d %s %s %c%s%c\n",
                                                          sign,sign,
                                                          sign,sign,
                                                          step,
                                                          resulttype,
                                                          resultplace,
                                                          sign,gpset,sign
                                                          );
   fprintf(out,"COMPONENTNAMES %cstress-tang%c,%cstress-normal%c,%cdummy%c,%cD-nomal%c,%cD-tang%c,%cdummy%c\n",
           sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actmafield->dis[0].numele; i++)
   {
     actmaele = &(actmafield->dis[0].element[i]);
     if(actmaele->e.w1->isinomegaprime == 1)
     {
       for (j=0; j<actsmfield->dis[0].numele; j++)
       {
         actsmele = &(actsmfield->dis[0].element[j]);
         if (actsmele->eltyp != el_interf) continue;
         GlobalID = i * actsmfield->dis[0].numele + (actsmele->Id + 1);
         /*  gid's 1.and 4.GP get values of my 1.GP */
         /*  gid's 2.and 3.GP get values of my 2.GP */
         fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                         GlobalID,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].T[0],
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].T[1],
                         0.0,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].dn,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].dt,
                               0.0
                         );
         fprintf(out,"      %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].T[0],
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].T[1],
                         0.0,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].dn,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].dt,
                               0.0
                         );
         fprintf(out,"       %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].T[0],
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].T[1],
                         0.0,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].dn,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].dt,
                               0.0
                         );
         fprintf(out,"      %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].T[0],
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].T[1],
                         0.0,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].dn,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].dt,
                               0.0
                         );
       }
     }
   }/* end: for (i=0; i<actmafield->dis[0].numele; i++) */
   fprintf(out,"END VALUES\n");
}/* end: if (actsmgid->is_interf_22) */
/*----------------------------------------------------------------------*/
if (actsmgid->is_interf_33)
{
   resulttype        = "MATRIX";
   resultplace       = "ONGAUSSPOINTS";
   gpset             = actsmgid->interf_33_name;
   rangetable        = actsmgid->standardrangetable;
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT INTERFACE on FIELD %s\n",actsmgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cinterface_stresses%c %cccarat%c %d %s %s %c%s%c\n",
                                                          sign,sign,
                                                          sign,sign,
                                                          step,
                                                          resulttype,
                                                          resultplace,
                                                          sign,gpset,sign
                                                          );
   fprintf(out,"COMPONENTNAMES %cstress-tang%c,%cstress-normal%c,%cdummy%c,%cD-normal%c,%cD-tang%c,%cdummy%c\n",
           sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);
   fprintf(out,"VALUES\n");
   for (i=0; i<actmafield->dis[0].numele; i++)
   {
     actmaele = &(actmafield->dis[0].element[i]);
     if(actmaele->e.w1->isinomegaprime == 1)
     {
       for (j=0; j<actsmfield->dis[0].numele; j++)
       {
         actsmele = &(actsmfield->dis[0].element[j]);
         if (actsmele->eltyp != el_interf) continue;
         GlobalID = i * actsmfield->dis[0].numele + (actsmele->Id + 1);
        /*  gid's 1.and 4.and 8. GP get values of my 1.GP */
        /*  gid's 5.and 7.and 9. GP get values of my 2.GP */
        /*  gid's 2.and 3.and 6. GP get values of my 3.GP */
         fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                         GlobalID,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].T[0],
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].T[1],
                         0.0,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].dn,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].dt,
                               0.0
                         );
         fprintf(out,"      %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[2].T[0],
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[2].T[1],
                         0.0,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[2].dn,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[2].dt,
                               0.0
                         );
         fprintf(out,"      %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[2].T[0],
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[2].T[1],
                         0.0,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[2].dn,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[2].dt,
                               0.0
                         );
         fprintf(out,"      %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].T[0],
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].T[1],
                         0.0,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].dn,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].dt,
                               0.0
                         );
         fprintf(out,"      %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].T[0],
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].T[1],
                         0.0,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].dn,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].dt,
                               0.0
                         );
         fprintf(out,"      %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[2].T[0],
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[2].T[1],
                         0.0,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[2].dn,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[2].dt,
                               0.0
                         );
         fprintf(out,"      %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].T[0],
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].T[1],
                         0.0,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].dn,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].dt,
                               0.0
                         );
         fprintf(out,"      %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].T[0],
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].T[1],
                         0.0,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].dn,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[0].dt,
                               0.0
                         );
         fprintf(out,"      %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].T[0],
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].T[1],
                         0.0,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].dn,
                         actmaele->e.w1->sm_eledata[actsmele->Id].sm_GPdata[1].dt,
                               0.0
                         );
       }
     }
   }/*end: for (i=0; i<actmafield->dis[0].numele; i++) */
   fprintf(out,"END VALUES\n");
} /*end: if (actgid->is_interf_33) */
/*----------------------------------------------------------------------*/
fflush(out);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of out_gid_smstress */

#endif /* D_MLSTRUCT */

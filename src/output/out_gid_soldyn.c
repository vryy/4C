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
#include "../ale2/ale2.h"
#include "../ale3/ale3.h"
#include "../interf/interf.h"
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
 |  routine to write solution of a step to GID           m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_soldyn(char string[], FIELD *actfield, INTRA  *actintra, INT step,
                   INT place, DOUBLE totaltime)
{
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
DOUBLE      **forces;
DOUBLE      **stress;

DOUBLE        a1,a2,a3,thick,scal,sdc;
INT           tot_numnp;

INT           numele;   /*for shell9*/
/* 
   gausspoint permutation :
   On the Gausspoint number i in Gid, the results of Carats GP number gausspermn[i]
   have to be written
*/   

INT           gaussperm4[4] = {3,1,0,2};
INT           gaussperm8[8] = {0,4,2,6,1,5,3,7};
INT           gaussperm9[9] = {8,2,0,6,5,1,3,7,4};
INT           gaussperm27[27] = {0,9,18,3,12,21,6,15,24,1,10,19,4,13,22,7,16,25,2,11,20,5,14,23,8,17,26};

#ifdef DEBUG 
dstrc_enter("out_gid_soldyn");
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
/*======================================================================*/
/*========================================= result type is displacement */
/*======================================================================*/
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
   fprintf(out,"# TIME %20.10f\n",totaltime);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %c%s%c %cpcarat%c %d %s %s\n",sign,string,sign,sign,sign,step,resulttype,resultplace);
   fprintf(out,"RESULTRANGESTABLE %c%s%c\n",sign,actgid->standardrangetable,sign);
   switch (genprob.ndim)
   {
     case 3:
       fprintf(out,"COMPONENTNAMES %c%s%c,%c%s%c,%c%s%c\n",sign,componentnames[0],sign,sign,componentnames[1],sign,sign,componentnames[2],sign);
     break;
     case 2:
       fprintf(out,"COMPONENTNAMES %c%s%c,%c%s%c\n",sign,componentnames[0],sign,sign,componentnames[1],sign);
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
      scal = 1.0;
      sdc  = actfield->dis[0].element[0].e.s8->sdc;
      for (i=0; i<actfield->dis[0].numnp; i++)
      {
         actnode = &(actfield->dis[0].node[i]);
         /* the lower surface */
         fprintf(out," %6d %23.15E %23.15E %23.15E\n",actnode->Id+1,
                                                      actnode->sol.a.da[place][0]-actnode->sol.a.da[place][3]*scal/sdc,
                                                      actnode->sol.a.da[place][1]-actnode->sol.a.da[place][4]*scal/sdc,
                                                      actnode->sol.a.da[place][2]-actnode->sol.a.da[place][5]*scal/sdc);
         /* the upper surface */
         fprintf(out," %6d %23.15E %23.15E %23.15E\n",actnode->Id+1+tot_numnp,
                                                      actnode->sol.a.da[place][0]+actnode->sol.a.da[place][3]*scal/sdc,
                                                      actnode->sol.a.da[place][1]+actnode->sol.a.da[place][4]*scal/sdc,
                                                      actnode->sol.a.da[place][2]+actnode->sol.a.da[place][5]*scal/sdc);
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
        fprintf(out," %6d %18.5E %18.5E %18.5E\n",actnode->Id+1,actnode->sol.a.da[place][0],actnode->sol.a.da[place][1],actnode->sol.a.da[place][2]);
	break;
	case 2:
        fprintf(out," %6d %18.5E %18.5E \n",actnode->Id+1,actnode->sol.a.da[place][0],actnode->sol.a.da[place][1]);
        break;
	default:
	  dserror("Unknown number of dimensions");
        break;
      }
   }
   fprintf(out,"END VALUES\n");
} /* end of (strncmp(string,"displacement",stringlenght)==0) */

/*======================================================================*/
/*============================================ result type is contact */
/*======================================================================*/
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
   fprintf(out,"# TIME %20.10f\n",totaltime);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %c%s%c %cpcarat%c %d %s %s\n",sign,string,sign,sign,sign,step,resulttype,resultplace);
   fprintf(out,"RESULTRANGESTABLE %c%s%c\n",sign,actgid->standardrangetable,sign);
   switch (genprob.ndim)
   {
     case 3:
       fprintf(out,"COMPONENTNAMES %c%s%c,%c%s%c,%c%s%c\n",sign,componentnames[0],sign,sign,componentnames[1],sign,sign,componentnames[2],sign);
     break;
     case 2:
       fprintf(out,"COMPONENTNAMES %c%s%c,%c%s%c\n",sign,componentnames[0],sign,sign,componentnames[1],sign);
     break;
     default:
       dserror("Unknown numer of dimensions");
     break;
   }
   fprintf(out,"VALUES\n");
#if 1 /* this is hexahedra output for shell8 element */
#ifdef D_SHELL8
   if (actfield->dis[0].element[0].eltyp == el_shell8 && actfield->dis[0].element[0].distyp == quad4)
   {
      tot_numnp = genprob.nnode;
      scal = 1.0;
      sdc  = actfield->dis[0].element[0].e.s8->sdc;
      for (i=0; i<actfield->dis[0].numnp; i++)
      {
         actnode = &(actfield->dis[0].node[i]);
         /* the lower surface */
         fprintf(out," %6d %23.15E %23.15E %23.15E\n",actnode->Id+1,
                                                    -actnode->sol.a.da[place][0]+actnode->sol.a.da[place][3]*scal/sdc,
                                                    -actnode->sol.a.da[place][1]+actnode->sol.a.da[place][4]*scal/sdc,
                                                    -actnode->sol.a.da[place][2]+actnode->sol.a.da[place][5]*scal/sdc);
         /* the upper surface */
         fprintf(out," %6d %23.15E %23.15E %23.15E\n",actnode->Id+1+tot_numnp,
                                                    -actnode->sol.a.da[place][0]-actnode->sol.a.da[place][3]*scal/sdc,
                                                    -actnode->sol.a.da[place][1]-actnode->sol.a.da[place][4]*scal/sdc,
                                                    -actnode->sol.a.da[place][2]-actnode->sol.a.da[place][5]*scal/sdc);
      }
   }
#endif
#else
   for (i=0; i<actfield->dis[0].numnp; i++)
   {
      actnode = &(actfield->dis[0].node[i]);
      switch (genprob.ndim)
      {
	case 3:
        fprintf(out," %6d %18.5E %18.5E %18.5E\n",actnode->Id+1,actnode->sol.a.da[place][0],actnode->sol.a.da[place][1],actnode->sol.a.da[place][2]);
	break;
	case 2:
        fprintf(out," %6d %18.5E %18.5E \n",actnode->Id+1,actnode->sol.a.da[place][0],actnode->sol.a.da[place][1]);
        break;
	default:
	  dserror("Unknown number of dimensions");
        break;
      }
   }
#endif
   fprintf(out,"END VALUES\n");
} /* end of (strncmp(string,"contact",stringlenght)==0) */
#endif
/*======================================================================*/
/*=============================================== result type is stress */
/*======================================================================*/
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
      fprintf(out,"# TIME %20.10f\n",totaltime);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cshell8_forces%c %cpcarat%c %d %s %s %c%s%c\n",sign,sign,sign,sign,step,resulttype,resultplace,sign,gpset,sign);
      fprintf(out,"COMPONENTNAMES %cForce-11%c,%cForce-22%c,%cForce-33%c,%cForce-12%c,%cForce-23%c,%cForce-13%c\n",sign,sign,sign,sign,sign,sign,sign,sign,sign,sign,sign,sign);
      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_shell8 || actele->numnp !=4) continue;
         forces = actele->e.s8->forces.a.d3[place];
         fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",actele->Id+1,
                             forces[0][gaussperm4[0]],
                             forces[1][gaussperm4[0]],
                             forces[9][gaussperm4[0]],
                             forces[2][gaussperm4[0]],
                             forces[4][gaussperm4[0]],
                             forces[3][gaussperm4[0]]);
         for (j=1; j<ngauss; j++)
         fprintf(out,"        %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                             forces[0][gaussperm4[j]],
                             forces[1][gaussperm4[j]],
                             forces[9][gaussperm4[j]],
                             forces[2][gaussperm4[j]],
                             forces[4][gaussperm4[j]],
                             forces[3][gaussperm4[j]]);
 
      }
      fprintf(out,"END VALUES\n");
      /*--- print the linear-in-thickness direction forces */
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT shell8_moments on FIELD %s\n",actgid->fieldname);
      fprintf(out,"# TIME %20.10f\n",totaltime);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cshell8_moments%c %cpcarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      fprintf(out,"COMPONENTNAMES %cMoment-11%c,%cMoment-22%c,%cMoment-33%c,%cMoment-12%c,%cMoment-23%c,%cMoment-13%c\n",sign,sign,sign,sign,sign,sign,sign,sign,sign,sign,sign,sign);
      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_shell8 || actele->numnp !=4) continue;
         forces = actele->e.s8->forces.a.d3[place];
         fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",actele->Id+1,
                             forces[5] [gaussperm4[0]],
                             forces[6] [gaussperm4[0]],
                             forces[15][gaussperm4[0]],
                             forces[7] [gaussperm4[0]],
                             forces[11][gaussperm4[0]],
                             forces[10][gaussperm4[0]]);
         for (j=1; j<ngauss; j++)
         fprintf(out,"        %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                             forces[5] [gaussperm4[j]],
                             forces[6] [gaussperm4[j]],
                             forces[15][gaussperm4[j]],
                             forces[7] [gaussperm4[j]],
                             forces[11][gaussperm4[j]],
                             forces[10][gaussperm4[j]]);
 
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
      fprintf(out,"# TIME %20.10f\n",totaltime);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cshell8_forces%c %cpcarat%c %d %s %s %c%s%c\n",sign,sign,sign,sign,step,resulttype,resultplace,sign,gpset,sign);
      fprintf(out,"COMPONENTNAMES %cForce-11%c,%cForce-22%c,%cForce-33%c,%cForce-12%c,%cForce-23%c,%cForce-13%c\n",sign,sign,sign,sign,sign,sign,sign,sign,sign,sign,sign,sign);
      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_shell8 || actele->numnp !=9) continue;
         forces = actele->e.s8->forces.a.d3[place];
         fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",actele->Id+1,
                             forces[0][gaussperm9[0]],
                             forces[1][gaussperm9[0]],
                             forces[9][gaussperm9[0]],
                             forces[2][gaussperm9[0]],
                             forces[4][gaussperm9[0]],
                             forces[3][gaussperm9[0]]);
         for (j=1; j<ngauss; j++)
         fprintf(out,"        %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                             forces[0][gaussperm9[j]],
                             forces[1][gaussperm9[j]],
                             forces[9][gaussperm9[j]],
                             forces[2][gaussperm9[j]],
                             forces[4][gaussperm9[j]],
                             forces[3][gaussperm9[j]]);
 
      }
      fprintf(out,"END VALUES\n");
      /*--- print the linear-in-thickness direction forces */
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# RESULT shell8_moments on FIELD %s\n",actgid->fieldname);
      fprintf(out,"# TIME %20.10f\n",totaltime);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cshell8_moments%c %cpcarat%c %d %s %s %c%s%c\n",sign,sign,sign,sign,step,resulttype,resultplace,sign,gpset,sign);
      fprintf(out,"COMPONENTNAMES %cMoment-11%c,%cMoment-22%c,%cMoment-33%c,%cMoment-12%c,%cMoment-23%c,%cMoment-13%c\n",sign,sign,sign,sign,sign,sign,sign,sign,sign,sign,sign,sign);
      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_shell8 || actele->numnp !=9) continue;
         forces = actele->e.s8->forces.a.d3[place];
         fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",actele->Id+1,
                             forces[5] [gaussperm9[0]],
                             forces[6] [gaussperm9[0]],
                             forces[15][gaussperm9[0]],
                             forces[7] [gaussperm9[0]],
                             forces[11][gaussperm9[0]],
                             forces[10][gaussperm9[0]]);
         for (j=1; j<ngauss; j++)
         fprintf(out,"        %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                             forces[5] [gaussperm9[j]],
                             forces[6] [gaussperm9[j]],
                             forces[15][gaussperm9[j]],
                             forces[7] [gaussperm9[j]],
                             forces[11][gaussperm9[j]],
                             forces[10][gaussperm9[j]]);
 
      }
      fprintf(out,"END VALUES\n");
   }
#endif /*D_SHELL8*/
#endif /* 0 */

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
      fprintf(out,"# TIME %20.10f\n",totaltime);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cwall1_forces%c %cpcarat%c %d %s %s %c%s%c\n",sign,sign,sign,sign,step,resulttype, resultplace,sign,gpset,sign);
      fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-xy%c,%cStress-zz%c,%cdummy%c,%cdummy%c\n",sign,sign,sign,sign,sign,sign,sign,sign,sign,sign,sign,sign);
      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
   /*---------------|| actele->numnp !=4--- */
         if (actele->eltyp != el_wall1 ) continue;
         stress=actele->e.w1->stress_GP.a.d3[place];
	 fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",actele->Id+1,
			     stress[0][gaussperm4[0]],
			     stress[1][gaussperm4[0]],
			     stress[2][gaussperm4[0]],
			     stress[3][gaussperm4[0]]);
         for (j=1; j<ngauss; j++)
         fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             stress[0][gaussperm4[j]],
			     stress[1][gaussperm4[j]],
			     stress[2][gaussperm4[j]],
			     stress[3][gaussperm4[j]]);
      }
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
      fprintf(out,"# TIME %20.10f\n",totaltime);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cwall1_forces%c %cpcarat%c %d %s %s %c%s%c\n",sign,sign,sign,sign,step,resulttype,resultplace,sign,gpset,sign);
      fprintf(out,"COMPONENTNAMES %cStress-xx%c,%cStress-yy%c,%cStress-xy%c,%cStress-zz%c,%cdummy%c,%cdummy%c\n",sign,sign,sign,sign,sign,sign,sign,sign,sign,sign,sign,sign);
      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
/*---------------|| actele->numnp < 8--- */
         if (actele->eltyp != el_wall1 ) continue;
         stress=actele->e.w1->stress_GP.a.d3[place];
	 fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",actele->Id+1,
			     stress[0][gaussperm9[0]],
			     stress[1][gaussperm9[0]],
			     stress[2][gaussperm9[0]],
			     stress[3][gaussperm9[0]]);
         for (j=1; j<ngauss; j++)
	 fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             stress[0][gaussperm9[j]],
			     stress[1][gaussperm9[j]],
			     stress[2][gaussperm9[j]],
			     stress[3][gaussperm9[j]]);
      }
      fprintf(out,"END VALUES\n");
/*----------------------------------------------------------------------*/      
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
      fprintf(out,"# TIME %20.10f\n",totaltime);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"RESULT %cinterface_stresses%c %cpcarat%c %d %s %s %c%s%c\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             step,
                                                             resulttype,
                                                             resultplace,
                                                             sign,gpset,sign
                                                             );
      actele = &(actfield->dis[0].element[0]);
      switch(actele->e.interf->stresstyp)
      {
      /*---------------------------------- local stresses are asked for ---*/
        case if_tn:
         fprintf(out,"COMPONENTNAMES %cstress-tang%c,%cstress-normal%c,%cdummy%c,%cdummy%c,%cdummy%c,%cdummy%c\n",
                 sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);
        break;
        /*--------------------------- transformation into global stresses ---*/
        case if_xy:
         fprintf(out,"COMPONENTNAMES %cstress-sxx%c,%cstress-syy%c,%cstress-sxy%c,%cdummy%c,%cdummy%c,%cdummy%c\n",
              sign,sign, sign,sign, sign,sign, sign,sign, sign,sign, sign,sign);
        break;
        default:
        break;
      }    
      fprintf(out,"VALUES\n");
      for (i=0; i<actfield->dis[0].numele; i++)
      {
         actele = &(actfield->dis[0].element[i]);
         if (actele->eltyp != el_interf || actele->numnp !=4) continue;

         /* ----------------------- gid's 1.and 4.GP get values of my 1.GP-- */
         /* ----------------------- gid's 2.and 3.GP get values of my 2.GP-- */
         stress=actele->e.interf->stress_GP.a.d3[place];
	 fprintf(out," %6d %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
                             actele->Id+1,
			     stress[0][0],
			     stress[1][0],
			     stress[2][0],
                             0.0,
                             0.0,
                             0.0
			     );
	 fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
			     stress[0][1],
			     stress[1][1],
			     stress[2][1],
                             0.0,
                             0.0,
                             0.0
			     );
	 fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
			     stress[0][1],
			     stress[1][1],
			     stress[2][1],
                             0.0,
                             0.0,
                             0.0
			     );
	 fprintf(out,"        %18.5E %18.5E %18.5E %18.5E %18.5E %18.5E \n",
			     stress[0][0],
			     stress[1][0],
			     stress[2][0],
                             0.0,
                             0.0,
                             0.0
			     );
         /* ------------------------------------------------- */
      }
      fprintf(out,"END VALUES\n");
   }
#endif
} /* end of (strncmp(string,"stress",stringlenght)==0) */
/*======================================================================*/
/*========================================= result type is velocity */
/*======================================================================*/
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
   fprintf(out,"# TIME %20.10f\n",totaltime);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %c%s%c %cpcarat%c %d %s %s\n",sign,string,sign,sign,sign,step,resulttype,resultplace);
   fprintf(out,"RESULTRANGESTABLE %c%s%c\n",sign,actgid->standardrangetable,sign);
   switch (genprob.ndim)
   {
     case 3:
       fprintf(out,"COMPONENTNAMES %c%s%c,%c%s%c,%c%s%c\n",sign,componentnames[0],sign,sign,componentnames[1],sign,sign,componentnames[2],sign);
     break;
     case 2:
       fprintf(out,"COMPONENTNAMES %c%s%c,%c%s%c\n",sign,componentnames[0],sign,sign,componentnames[1],sign);
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
        fprintf(out," %6d %18.5E %18.5E %18.5E\n",actnode->Id+1,
                                                  actnode->sol.a.da[place][0],
                                                  actnode->sol.a.da[place][1],
                                                  actnode->sol.a.da[place][2]);
	break;
	case 2:
        fprintf(out," %6d %18.5E %18.5E \n",actnode->Id+1,
                                            actnode->sol.a.da[place][0],
                                            actnode->sol.a.da[place][1]);
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
} /* end of out_gid_soldyn */




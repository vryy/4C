#include "../headers/standardtypes.h"

#ifdef D_SHELL8
  #include "../shell8/shell8.h"
#endif /*D_SHELL8*/
#ifdef D_SHELL9
  #include "../shell9/shell9.h"
#endif /*D_SHELL9*/
#ifdef D_WALL1
  #include "../wall1/wall1.h"
#endif /*D_WALL1*/
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
struct _IO_FLAGS        ioflags;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 *----------------------------------------------------------------------*/
extern struct _DYNAMIC *dyn;   
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA         *alldyn;   
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
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
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
 | pointer to allocate static variables if needed                       |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _STATIC_VAR  *statvar;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;                      

/*----------------------------------------------------------------------*
 |  print out general information about problem              m.gee 12/01|
 *----------------------------------------------------------------------*/
void out_general()
{
int        i,j,k,l;
int        is_shell9;
FILE      *out = allfiles.out_out;
FIELD     *actfield;
INTRA     *actintra;
ELEMENT   *actele;
NODE      *actnode;
int        myrank;
int        nprocs;
int        imyrank;
int        inprocs;

#ifdef DEBUG 
dstrc_enter("out_general");
#endif

/*----------------------------------------------------------------------*/
myrank = par.myrank;
nprocs = par.nprocs;
/*----------------------------------------------------------------------*/
if (myrank==0)
{
/*-------------------------------------------------------- print header */
fprintf(out,"________________________________________________________________________________\n\n");
fprintf(out,"p_carat outputfile\n");
fprintf(out,"________________________________________________________________________________\n\n");
for (i=0; i<5; i++)
fprintf(out,"%s\n",allfiles.title[i]);
fprintf(out,"________________________________________________________________________________\n");
fprintf(out,"\n");
/*------------------------------------------ print general problem data */
fprintf(out,"________________________________________________________________________________\n\n");
fprintf(out,"Total number of Fields    : %d\n",genprob.numfld);
fprintf(out,"Total number of Elements  : %d\n",genprob.nele);
fprintf(out,"Total number of Nodes     : %d\n",genprob.nnode);
fprintf(out,"Total number of Materials : %d\n",genprob.nmat);
switch(genprob.probtyp)
{
case prb_fsi:
fprintf(out,"Type of Problem           : Fluid-Structure-Interaction\n");
break;
case prb_structure:
fprintf(out,"Type of Problem           : Structural\n");
break;
case prb_fluid:
fprintf(out,"Type of Problem           : Fluid\n");
break;
case prb_opt:
fprintf(out,"Type of Problem           : Optimization\n");
break;
case prb_ale:
fprintf(out,"Type of Problem           : Ale\n");
break;
default:
dserror("Cannot print problem type");
break;
}
switch(genprob.timetyp)
{
case time_static:
fprintf(out,"Type of Time              : Static\n");
fprintf(out,"Total Number of Steps     : %d\n",statvar->nstep);
break;
case time_dynamic:
fprintf(out,"Type of Time              : Dynamic\n");
/*fprintf(out,"Total Number of Steps     : %d\n",dyn->nstep);*/
break;
default:
dserror("Cannot print time type");
break;
}
fprintf(out,"________________________________________________________________________________\n");
fprintf(out,"\n");
/*----------------------------------------------------------------------*/
for (i=0; i<genprob.numfld; i++)
{
   actfield = &(field[i]);
#ifdef PARALLEL 
   actintra = &(par.intra[i]);
#else
   actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
   if (!actintra) dserror("Allocation of INTRA failed");
   actintra->intra_rank     = 0;
   actintra->intra_nprocs   = 1;
#endif
fprintf(out,"================================================================================\n");
switch(actfield->fieldtyp)
{
case fluid:
fprintf(out,"FIELD: fluid\n");
break;
case ale:
fprintf(out,"FIELD: ale\n");
break;
case structure:
fprintf(out,"FIELD: structure\n");
break;
default:
dserror("Cannot print fieldtype");
break;
}
fprintf(out,"================================================================================\n");
fprintf(out,"________________________________________________________________________________\n\n");
fprintf(out,"Number of Elements  in this field : %d\n",actfield->dis[0].numele);
fprintf(out,"Number of Nodes     in this field : %d\n",actfield->dis[0].numnp);
fprintf(out,"Number of Dofs      in this field : %d\n",actfield->dis[0].numdf);
fprintf(out,"Number of Equations in this field : %d\n",actfield->dis[0].numeq);
fprintf(out,"________________________________________________________________________________\n\n");
fprintf(out,"Element connectivity in global Ids:\n");
for (j=0; j<actfield->dis[0].numele; j++)
{
actele = &(actfield->dis[0].element[j]);
fprintf(out,"glob_Id %6d Nnodes %2d Nodes: ",actele->Id,actele->numnp);
for (k=0; k<actele->numnp; k++) fprintf(out,"%6d ",actele->node[k]->Id);
fprintf(out,"\n");
}
fprintf(out,"Element connectivity in field-local Ids:\n");
for (j=0; j<actfield->dis[0].numele; j++)
{
actele = &(actfield->dis[0].element[j]);
fprintf(out,"loc_Id %6d Nnodes %2d Nodes: ",actele->Id_loc,actele->numnp);
for (k=0; k<actele->numnp; k++) fprintf(out,"%6d ",actele->node[k]->Id_loc);
fprintf(out,"\n");
}
fprintf(out,"________________________________________________________________________________\n\n");
fprintf(out,"Element types:\n");
for (j=0; j<actfield->dis[0].numele; j++)
{
actele = &(actfield->dis[0].element[j]);
switch(actele->eltyp)
{
case el_shell8:
fprintf(out,"ELE glob_Id %6d loc_Id %6d SHELL8\n",actele->Id,actele->Id_loc);
break;
case el_shell9:
fprintf(out,"ELE glob_Id %6d loc_Id %6d SHELL9\n",actele->Id,actele->Id_loc);
break;
case el_brick1:
fprintf(out,"ELE glob_Id %6d loc_Id %6d BRICK1\n",actele->Id,actele->Id_loc);
break;
case el_wall1:
fprintf(out,"ELE glob_Id %6d loc_Id %6d WALL1\n",actele->Id,actele->Id_loc);
break;
case el_fluid3:
fprintf(out,"ELE glob_Id %6d loc_Id %6d FLUID3\n",actele->Id,actele->Id_loc);
break;
case el_fluid2:
fprintf(out,"ELE glob_Id %6d loc_Id %6d FLUID2\n",actele->Id,actele->Id_loc);
break;
case el_ale3:
fprintf(out,"ELE glob_Id %6d loc_Id %6d ALE3\n",actele->Id,actele->Id_loc);
break;
case el_ale2:
fprintf(out,"ELE glob_Id %6d loc_Id %6d ALE2\n",actele->Id,actele->Id_loc);
break;
default:
dserror("Cannot print elementtype");
break;
}
}
fprintf(out,"________________________________________________________________________________\n\n");
fprintf(out,"Nodal Coordinates:\n");
for (j=0; j<actfield->dis[0].numnp; j++)
{
actnode = &(actfield->dis[0].node[j]);
fprintf(out,"NODE glob_Id %6d loc_Id %6d    %-18.5f %-18.5f %-18.5f \n",
        actnode->Id,actnode->Id_loc,actnode->x[0],actnode->x[1],actnode->x[2]);
}
fprintf(out,"________________________________________________________________________________\n\n");
#ifdef DEBUG
fprintf(out,"Degrees of Freedom:\n");
for (j=0; j<actfield->dis[0].numnp; j++)
{
actnode = &(actfield->dis[0].node[j]);

/* check if actnode belongs to a shell9 element */
#ifdef D_SHELL9
   is_shell9 = 0;
   for (k=0; k<actnode->numele; k++)
   {
     switch(actnode->element[k]->eltyp)
     {
     case el_shell9:
       is_shell9 = 1;    
     break;
     }
   }
   /* write dofs for shell9 */
   if (is_shell9 == 1) 
   {
     s9out_dof(actnode,j);
     continue;
   }
#endif /*D_SHELL9*/

switch (actnode->numdf)
{
case 2:
fprintf(out,"NODE glob_Id %6d loc_Id %6d    %6d    %6d\n",
        actnode->Id,actnode->Id_loc,actnode->dof[0],actnode->dof[1]);
break;
case 3:
fprintf(out,"NODE glob_Id %6d loc_Id %6d    %6d    %6d   %6d\n",
        actnode->Id,actnode->Id_loc,actnode->dof[0],actnode->dof[1],actnode->dof[2]);
break;
case 4:
fprintf(out,"NODE glob_Id %6d loc_Id %6d    %6d    %6d   %6d   %6d\n",
        actnode->Id,actnode->Id_loc,actnode->dof[0],actnode->dof[1],actnode->dof[2],actnode->dof[3]);
break;
case 5:
fprintf(out,"NODE glob_Id %6d loc_Id %6d    %6d    %6d   %6d   %6d   %6d\n",
        actnode->Id,actnode->Id_loc,actnode->dof[0],actnode->dof[1],actnode->dof[2],actnode->dof[3],actnode->dof[4]);
break;
case 6:
fprintf(out,"NODE glob_Id %6d loc_Id %6d    %6d    %6d   %6d   %6d   %6d   %6d\n",
        actnode->Id,actnode->Id_loc,actnode->dof[0],actnode->dof[1],actnode->dof[2],actnode->dof[3],actnode->dof[4],actnode->dof[5]);
break;
default:
break;
}
}
fprintf(out,"________________________________________________________________________________\n\n");
#endif
/* .... other stuff */




} /* end of (i=0; i<genprob.numfld; i++) */
/*----------------------------------------------------------------------*/
} /* end of if (myrank==0) */
/*----------------------------------------------------------------------*/
fflush(out);
#ifndef PARALLEL 
CCAFREE(actintra);
#endif
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_general */



/*----------------------------------------------------------------------*
 |  print out solution of a certain step                     m.gee 12/01|
 *----------------------------------------------------------------------*/
void out_sol(FIELD *actfield, PARTITION *actpart, INTRA *actintra, 
             int step, int place)
{
int        i,j,k,l;
int        is_shell9;    /*->shell9*/
int        num_klay,kl;  /*->shell9*/
FILE      *out = allfiles.out_out;
NODE      *actnode;
ELEMENT   *actele;
int        myrank;
int        nprocs;
int        imyrank;
int        inprocs;
int        ngauss;
int	     numnp;
int	     numff;
double     visc;
static FLUID_DYN_CALC  *dynvar;
#ifdef DEBUG 
dstrc_enter("out_sol");
#endif
/*----------------------------------------------------------------------*/
myrank = par.myrank;
nprocs = par.nprocs;
imyrank= actintra->intra_rank;
inprocs= actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
if (imyrank==0 && myrank==0)
{
/*-------------------------------------------------------- print header */
/*----------------------------------------------------------------------*/
fprintf(out,"================================================================================\n");
switch(actfield->fieldtyp)
{
case fluid:
fprintf(out,"FIELD: fluid\n");
break;
case ale:
fprintf(out,"FIELD: ale\n");
break;
case structure:
fprintf(out,"FIELD: structure\n");
break;
default:
dserror("Cannot print fieldtype");
break;
}
fprintf(out,"================================================================================\n");
fprintf(out,"Converged Solution in step %d\n",step); 
fprintf(out,"================================================================================\n");
/*-------------------------------------------------- print nodal values */
switch(actfield->fieldtyp)
{
case structure:
   if (ioflags.struct_disp_file==1)
   for (j=0; j<actfield->dis[0].numnp; j++)
   {
      actnode = &(actfield->dis[0].node[j]);

      /* check if actnode belongs to a shell9 element */
      #ifdef D_SHELL9
         is_shell9 = 0;
         for (k=0; k<actnode->numele; k++)
         {
           switch(actnode->element[k]->eltyp)
           {
           case el_shell9:
             is_shell9 = 1;
           break;
           }
         }
         /* print nodal values for shell9 */
         if (is_shell9 == 1)
         {
           if(actnode->numdf == 6) /*only one kinematic layer*/
           {
             fprintf(out,"NODE glob_Id %6d loc_Id %6d    ",actnode->Id,actnode->Id_loc);
             for (k=0; k<actnode->numdf; k++)
             {
                if (place >= actnode->sol.fdim) dserror("Cannot print solution step");
                fprintf(out,"%20.7E ",actnode->sol.a.da[place][k]);
             }
             fprintf(out,"\n");
           }
           else
             s9out_nodal_dis(actnode,place);
             continue;
         }
      #endif /*D_SHELL9*/

      fprintf(out,"NODE glob_Id %6d loc_Id %6d    ",actnode->Id,actnode->Id_loc);
      for (k=0; k<actnode->numdf; k++) 
      { 
         if (place >= actnode->sol.fdim) dserror("Cannot print solution step");
         fprintf(out,"%20.7#E ",actnode->sol.a.da[place][k]);
      }
      fprintf(out,"\n");     
   }
   fprintf(out,"________________________________________________________________________________\n\n");  
break;
case fluid:
   if (ioflags.fluid_sol_file==1)
   for (j=0; j<actfield->dis[0].numnp; j++)
   {
      actnode = &(actfield->dis[0].node[j]);
      fprintf(out,"NODE glob_Id %6d loc_Id %6d    ",actnode->Id,actnode->Id_loc);
      for (k=0; k<actnode->numdf; k++) 
      {
   	 if (place >= actnode->sol.fdim) dserror("Cannot print solution step");
   	 fprintf(out,"%20.7#E ",actnode->sol.a.da[place][k]);
      }
      fprintf(out,"\n");
   }
   fprintf(out,"________________________________________________________________________________\n\n");
break;
case ale:
   if (ioflags.ale_disp_file==1)
   for (j=0; j<actfield->dis[0].numnp; j++)
   {
      actnode = &(actfield->dis[0].node[j]);
      fprintf(out,"NODE glob_Id %6d loc_Id %6d    ",actnode->Id,actnode->Id_loc);
      for (k=0; k<actnode->numdf; k++) 
      { 
   	 if (place >= actnode->sol.fdim) dserror("Cannot print solution step");
   	 fprintf(out,"%20.7#E ",actnode->sol.a.da[place][k]);
      }
      fprintf(out,"\n");
   }
   fprintf(out,"________________________________________________________________________________\n\n");
break;
default: dserror("Cannot print fieldtype");
}
/*------------------------------------------------ print element values */
switch(actfield->fieldtyp)
{
case structure:
if (ioflags.struct_stress_file==1)
for (j=0; j<actfield->dis[0].numele; j++)
{
   actele = &(actfield->dis[0].element[j]);
   switch(actele->eltyp)
   {
   case el_shell8:
#ifdef D_SHELL8
       ngauss = actele->e.s8->nGP[0] * actele->e.s8->nGP[1];
       fprintf(out,"________________________________________________________________________________\n");
       fprintf(out,"Element glob_Id %d loc_Id %d                SHELL8\n",actele->Id,actele->Id_loc);
       switch(actele->e.s8->forcetyp)
       {
       case s8_xyz:
       fprintf(out,"Gaussian     Force-xx     Force-xy     Force-yx      Force-yy     Force-xz     Force-zx      Force-yz    Force-zy     Force-zz\n");
       break;
       case s8_rst:
       fprintf(out,"Gaussian     Force-rr     Force-rs     Force-sr      Force-ss     Force-rt     Force-tr      Force-st    Force-ts     Force-tt\n");
       break;
       case s8_rst_ortho:
       fprintf(out,"Gaussian     Force-rr     Force-rs     Force-sr      Force-ss     Force-rt     Force-tr      Force-st    Force-ts     Force-tt\n");
       break;
       default:
          dserror("Unknown type of element stresses");
       }
       for (i=0; i<ngauss; i++)
       {
       fprintf(out,"Gauss %d   %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
       i,
       actele->e.s8->forces.a.d3[place][0][i],
       actele->e.s8->forces.a.d3[place][2][i],
       actele->e.s8->forces.a.d3[place][8][i],
       actele->e.s8->forces.a.d3[place][1][i],
       actele->e.s8->forces.a.d3[place][3][i],
       actele->e.s8->forces.a.d3[place][16][i],
       actele->e.s8->forces.a.d3[place][4][i],
       actele->e.s8->forces.a.d3[place][17][i],
       actele->e.s8->forces.a.d3[place][9][i]
       );
       }
       switch(actele->e.s8->forcetyp)
       {
       case s8_xyz:
       fprintf(out,"Gaussian     Moment-xx    Moment-xy    Moment-yx     Moment-yy    Moment-xz    Moment-zx     Moment-yz    Moment-zy    Moment-zz\n"); 
       break;
       case s8_rst:
       fprintf(out,"Gaussian     Moment-rr    Moment-rs    Moment-sr     Moment-ss    Moment-rt    Moment-tr     Moment-st    Moment-ts    Moment-tt\n"); 
       break;
       case s8_rst_ortho:
       fprintf(out,"Gaussian     Moment-rr    Moment-rs    Moment-sr     Moment-ss    Moment-rt    Moment-tr     Moment-st    Moment-ts    Moment-tt\n"); 
       break;
       default:
          dserror("Unknown type of element stresses");
       }
       for (i=0; i<ngauss; i++)
       {
       fprintf(out,"Gauss %d   %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
       i,
       actele->e.s8->forces.a.d3[place][5][i],
       actele->e.s8->forces.a.d3[place][7][i],
       actele->e.s8->forces.a.d3[place][14][i],
       actele->e.s8->forces.a.d3[place][6][i],
       actele->e.s8->forces.a.d3[place][10][i],
       actele->e.s8->forces.a.d3[place][12][i],
       actele->e.s8->forces.a.d3[place][11][i],
       actele->e.s8->forces.a.d3[place][13][i],
       actele->e.s8->forces.a.d3[place][15][i]
       );
       }
#endif /*D_SHELL8*/
   break;
/*---------------------------------------------------------sh 9/02-------*/  
   case el_shell9:
#ifdef D_SHELL9
/*NOTE: It does not seem to be very interesting to write the forces of each kinematic layer; it is more
        usefull to look at the stresses. These are writen to the flavia.res-File, so that they could be
        visualized with gid.   sh 02/03 */
        
/*       ngauss   = actele->e.s9->nGP[0] * actele->e.s9->nGP[1];
       num_klay = actele->e.s9->num_klay;
       fprintf(out,"________________________________________________________________________________\n");
       fprintf(out,"Element glob_Id %d loc_Id %d                SHELL9\n",actele->Id,actele->Id_loc);
       switch(actele->e.s9->forcetyp)
       {
       case s9_xyz:
       fprintf(out,"Gaussian   Layer      Force-xx     Force-xy     Force-yx     Force-yy     Force-xz     Force-zx     Force-yz     Force-zy     Force-zz\n");
       break;
       case s9_rst:
       fprintf(out,"Gaussian   Layer      Force-rr     Force-rs     Force-sr     Force-ss     Force-rt     Force-tr     Force-st     Force-ts     Force-tt\n");
       break;
       case s9_rst_ortho:
       fprintf(out,"Gaussian   Layer      Force-rr     Force-rs     Force-sr     Force-ss     Force-rt     Force-tr     Force-st     Force-ts     Force-tt\n");
       break;
       default:
          dserror("Unknown type of element stresses");
       }
       for (i=0; i<ngauss; i++)
       {
          for (kl=0; kl<num_klay; kl++)  /*write forces for each kinematic layer*/
/*          {
          fprintf(out,"Gauss %d    Layer %d %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E \n",
          i,
          kl,
          actele->e.s9->forces.a.d3[place][0][(i*num_klay)+kl],
          actele->e.s9->forces.a.d3[place][2][(i*num_klay)+kl],
          actele->e.s9->forces.a.d3[place][8][(i*num_klay)+kl],
          actele->e.s9->forces.a.d3[place][1][(i*num_klay)+kl],
          actele->e.s9->forces.a.d3[place][3][(i*num_klay)+kl],
          actele->e.s9->forces.a.d3[place][16][(i*num_klay)+kl],
          actele->e.s9->forces.a.d3[place][4][(i*num_klay)+kl],
          actele->e.s9->forces.a.d3[place][17][(i*num_klay)+kl],
          actele->e.s9->forces.a.d3[place][9][(i*num_klay)+kl]
          );
          }
       }
       switch(actele->e.s9->forcetyp)
       {
       case s9_xyz:
       fprintf(out,"Gaussian   Layer      Moment-xx    Moment-xy    Moment-yx    Moment-yy    Moment-xz    Moment-zx    Moment-yz    Moment-zy    Moment-zz\n"); 
       break;
       case s9_rst:
       fprintf(out,"Gaussian   Layer      Moment-rr    Moment-rs    Moment-sr    Moment-ss    Moment-rt    Moment-tr    Moment-st    Moment-ts    Moment-tt\n"); 
       break;
       case s9_rst_ortho:
       fprintf(out,"Gaussian   Layer      Moment-rr    Moment-rs    Moment-sr    Moment-ss    Moment-rt    Moment-tr    Moment-st    Moment-ts    Moment-tt\n"); 
       break;
       default:
          dserror("Unknown type of element stresses");
       }
       for (i=0; i<ngauss; i++)
       {
          for (kl=0; kl<num_klay; kl++)  /*write forces for each kinematic layer*/
/*          {
          fprintf(out,"Gauss %d    Layer %d %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E \n",
          i,
          kl,
          actele->e.s9->forces.a.d3[place][5][(i*num_klay)+kl],
          actele->e.s9->forces.a.d3[place][7][(i*num_klay)+kl],
          actele->e.s9->forces.a.d3[place][14][(i*num_klay)+kl],
          actele->e.s9->forces.a.d3[place][6][(i*num_klay)+kl],
          actele->e.s9->forces.a.d3[place][10][(i*num_klay)+kl],
          actele->e.s9->forces.a.d3[place][12][(i*num_klay)+kl],
          actele->e.s9->forces.a.d3[place][11][(i*num_klay)+kl],
          actele->e.s9->forces.a.d3[place][13][(i*num_klay)+kl],
          actele->e.s9->forces.a.d3[place][15][(i*num_klay)+kl]
          );
          }
       }
*/#endif /*D_SHELL9*/
   break;
/*---------------------------------------------------------fh 03/02-------*/  
   case el_wall1:
#ifdef D_WALL1
       ngauss = actele->e.w1->nGP[0] * actele->e.w1->nGP[1];
       fprintf(out,"________________________________________________________________________________\n");
       fprintf(out,"Element glob_Id %d loc_Id %d                WALL1\n",actele->Id,actele->Id_loc);
       fprintf(out,"\n");
/* check wether stresses at Gauss Points are presented in global xy- or local rs-coordinate system and write stress type */
       switch(actele->e.w1->stresstyp)
       {
       case w1_xy:
       fprintf(out,"Gaussian     Stress-xx    Stress-yy    Stress-xy    Stress-zz    Max. P.S.    Min. P.S.    Angle\n");
       break;
       case w1_rs:
       fprintf(out,"Gaussian     Stress-rr    Stress-ss    Stress-rs    Stress-zz    Max. P.S.    Min. P.S.    Angle\n");
       break;
       default:
       fprintf(out,"Gaussian     Stress-xx    Stress-yy    Stress-xy    Stress-zz    Max. P.S.    Min. P.S.    Angle\n");
       }
       for (i=0; i<ngauss; i++)
       {
       fprintf(out,"Gauss %d   %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
       i,
       actele->e.w1->stress_GP.a.d3[place][0][i],
       actele->e.w1->stress_GP.a.d3[place][1][i],
       actele->e.w1->stress_GP.a.d3[place][2][i],
       actele->e.w1->stress_GP.a.d3[place][3][i],
       actele->e.w1->stress_GP.a.d3[place][4][i],
       actele->e.w1->stress_GP.a.d3[place][5][i],
       actele->e.w1->stress_GP.a.d3[place][6][i]
       );
       }           
#endif /*D_WALL1*/   
   break;
   
   default:
      dserror("unknown type of element");
   break;
   }
}
break;
case fluid:
break;
case ale:
break;
} /* end switch(actfield->fieldtyp) */

fprintf(out,"\n");
fprintf(out,"\n");

#if 0

for (j=0; j<actfield->dis[0].numele; j++)
{
   actele = &(actfield->dis[0].element[j]);
   switch(actele->eltyp)
   {
   case el_wall1:
   #ifdef D_WALL1
       
       numnp = actele->numnp;
       fprintf(out,"________________________________________________________________________________\n");
       fprintf(out,"Element glob_Id %d loc_Id %d                WALL1\n",actele->Id,actele->Id_loc);
       fprintf(out,"\n");
/* check wether stresses at FE Nodes are presented in global xy- or local rs-coordinate system and write stress type */
       switch(actele->e.w1->stresstyp)
       {
       case w1_xy:
       fprintf(out,"Nodal            Stress-xx    Stress-yy    Stress-xy    Stress-zz    Max. P.S.    Min. P.S.    Angle\n");
       break;
       case w1_rs:
       fprintf(out,"Nodal            Stress-rr    Stress-ss    Stress-rs    Stress-zz    Max. P.S.    Min. P.S.    Angle\n");
       break;
       default:
       fprintf(out,"Nodal            Stress-xx    Stress-yy    Stress-xy    Stress-zz    Max. P.S.    Min. P.S.    Angle\n");
       break;
       }
       for (i=0; i<numnp; i++)
       {
       fprintf(out,"Node %d  IEL %d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
       actele->node[i]->Id,
       i,
       actele->e.w1->stress_ND.a.d3[place][0][i],
       actele->e.w1->stress_ND.a.d3[place][1][i],
       actele->e.w1->stress_ND.a.d3[place][2][i],
       actele->e.w1->stress_ND.a.d3[place][3][i],
       actele->e.w1->stress_ND.a.d3[place][4][i],
       actele->e.w1->stress_ND.a.d3[place][5][i],
       actele->e.w1->stress_ND.a.d3[place][6][i]
       );
       }
   #endif /*D_WALL1*/                       
   break;
   
   }
}

#endif
/*----------------------------------------------------------------------*/  



/*----------------------------------------------------------------------*/
} /* end of if (myrank==0 && imyrank==0) */
/*----------------------------------------------------------------------*/
if (myrank==0 && imyrank==0) fflush(out);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_general */

/*!---------------------------------------------------------------------
\brief print fsi coupling informations

<pre>                                                         genk 01/03
		     
</pre>

\param *fluidfield    FIELD   (i)    

\return void                                                                             

------------------------------------------------------------------------*/
void out_fsi(FIELD *fluidfield)
{
int        i,j,k,l;
FILE      *out = allfiles.out_out;
NODE      *actfnode, *actanode, *actsnode;
GNODE     *actfgnode;
int        myrank;
int        nprocs;
int        imyrank;
int        inprocs;
int	   numnp;
int        numaf,numsf;

#ifdef DEBUG 
dstrc_enter("out_fsi");
#endif

#ifdef D_FSI
/*----------------------------------------------------------------------*/
myrank = par.myrank;
nprocs = par.nprocs;
numaf  = genprob.numaf;
numsf  = genprob.numsf;
/*----------------------------------------------------------------------*/
if (myrank==0)
{
fprintf(out,"================================================================================\n");
fprintf(out,"FSI node connectivity global Ids:\n");
fprintf(out,"================================================================================\n");
fprintf(out,"\n");
fprintf(out,"FLUID          ALE          STRUCTURE\n");
for (i=0;i<fluidfield->dis[0].numnp;i++)
{
   actfnode  = &(fluidfield->dis[0].node[i]);
   actfgnode = actfnode->gnode;
   actsnode  = actfgnode->mfcpnode[numsf];
   actanode  = actfgnode->mfcpnode[numaf];
   if (actsnode==NULL && actanode!=NULL)
   fprintf(out,"%-6d         %-6d       ------\n",actfnode->Id,actanode->Id);
   else if (actanode==NULL && actsnode!=NULL)
   fprintf(out,"%-6d         ------       %-6d\n",actfnode->Id,actsnode->Id);
   else if (actanode!=NULL && actsnode!=NULL)
   fprintf(out,"%-6d         %-6d       %-6d\n",actfnode->Id,actanode->Id,actsnode->Id); 
   else
   fprintf(out,"%-6d         ------       ------\n",actfnode->Id);    

}
fprintf(out,"________________________________________________________________________________\n\n");
fprintf(out,"================================================================================\n");
fprintf(out,"FSI node connectivity local Ids:\n");
fprintf(out,"================================================================================\n");
fprintf(out,"\n");
fprintf(out,"FLUID          ALE          STRUCTURE\n");
for (i=0;i<fluidfield->dis[0].numnp;i++)
{
   actfnode  = &(fluidfield->dis[0].node[i]);
   actfgnode = actfnode->gnode;
   actsnode  = actfgnode->mfcpnode[numsf];
   actanode  = actfgnode->mfcpnode[numaf];
   if (actsnode==NULL && actanode!=NULL)
   fprintf(out,"%-6d         %-6d       ------\n",actfnode->Id_loc,actanode->Id_loc);
   else if (actanode==NULL && actsnode!=NULL)
   fprintf(out,"%-6d         ------       %-6d\n",actfnode->Id_loc,actsnode->Id_loc);
   else if (actanode!=NULL && actsnode!=NULL)
   fprintf(out,"%-6d         %-6d       %-6d\n",actfnode->Id_loc,actanode->Id,actsnode->Id_loc); 
   else
   fprintf(out,"%-6d         ------       ------\n",actfnode->Id_loc);    

}
fprintf(out,"________________________________________________________________________________\n\n");
/*----------------------------------------------------------------------*/
} /* end of if (myrank==0 && imyrank==0) */
/*----------------------------------------------------------------------*/
if (myrank==0) fflush(out);
/*----------------------------------------------------------------------*/
#else
dserror("FSI functions not compiled in\n");
#endif
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_fsi */

/*!---------------------------------------------------------------------
\brief  print fluid multifield coupling informations

<pre>                                                         genk 01/03
		     
</pre>

\param *fluidfield    FIELD   (i)    

\return void                                                                             

------------------------------------------------------------------------*/
void out_fluidmf(FIELD *fluidfield)
{
int        i,j,k,l;
FILE      *out = allfiles.out_out;
NODE      *actfnode, *actanode;
GNODE     *actfgnode;
int        myrank;
int        nprocs;
int        imyrank;
int        inprocs;
int	   numnp;
int        numff;
int        numaf;

#ifdef DEBUG 
dstrc_enter("out_fluidmf");
#endif

#ifdef D_FSI
/*----------------------------------------------------------------------*/
myrank = par.myrank;
nprocs = par.nprocs;
numff  = genprob.numff;
numaf  = genprob.numaf;
/*----------------------------------------------------------------------*/
if (myrank==0)
{
fprintf(out,"================================================================================\n");
fprintf(out,"Fluid multiefield node connectivity global Ids:\n");
fprintf(out,"================================================================================\n");
fprintf(out,"\n");
fprintf(out,"FLUID          ALE          \n");
for (i=0;i<fluidfield->dis[0].numnp;i++)
{
   actfnode  = &(fluidfield->dis[0].node[i]);
   actfgnode = actfnode->gnode;
   actanode  = actfgnode->mfcpnode[numaf];
   if (actanode!=NULL)
   fprintf(out,"%-6d         %-6d       \n",actfnode->Id,actanode->Id); 
   else
   fprintf(out,"%-6d         ------       \n",actfnode->Id);    

}
fprintf(out,"________________________________________________________________________________\n\n");
fprintf(out,"================================================================================\n");
fprintf(out,"Fluid multiefield node connectivity global Ids:\n");
fprintf(out,"================================================================================\n");
fprintf(out,"\n");
fprintf(out,"FLUID          ALE          \n");
for (i=0;i<fluidfield->dis[0].numnp;i++)
{
   actfnode  = &(fluidfield->dis[0].node[i]);
   actfgnode = actfnode->gnode;
   actanode  = actfgnode->mfcpnode[numaf];
   if (actanode!=NULL)
   fprintf(out,"%-6d         %-6d       \n",actfnode->Id_loc,actanode->Id_loc);
   else
   fprintf(out,"%-6d         ------       \n",actfnode->Id_loc);    

}
fprintf(out,"________________________________________________________________________________\n\n");
/*----------------------------------------------------------------------*/
} /* end of if (myrank==0 && imyrank==0) */
/*----------------------------------------------------------------------*/
if (myrank==0) fflush(out);
/*----------------------------------------------------------------------*/
#else
dserror("FSI functions not compiled in\n");
#endif
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_fluidmf */

/*----------------------------------------------------------------------*
 |  print out solution of turbulence of a certain step      he    05/03 |
 *----------------------------------------------------------------------*/
void out_fluidtu(FIELD *actfield, INTRA *actintra, int step, int place)
{
int        i,j,k,l;
FILE      *out = allfiles.out_tur;
NODE      *actnode;
NODE      *actnode2;
GNODE     *actgnode2;	            
ELEMENT   *actele;
int        myrank;
int        nprocs;
int        imyrank;
int        inprocs;
int	     numnp;
int	     numff;
double     visc;
static FLUID_DYN_CALC  *dynvar;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("out_fluidmtu");
#endif
/*----------------------------------------------------------------------*/

for (numff=0;numff<genprob.numfld;numff++)
{
 actfield=&(field[numff]);
 if (actfield->fieldtyp==fluid)
 break;
} /* end loop over numff */
dynvar = &(alldyn[numff].fdyn->dynvar);

for(i=0; i<genprob.nmat; i++)
{
 if (mat[i].mattyp == m_fluid) visc=mat->m.fluid->viscosity;
}
/*----------------------------------------------------------------------*/
myrank = par.myrank;
nprocs = par.nprocs;
imyrank= actintra->intra_rank;
inprocs= actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
if (imyrank==0 && myrank==0)
{
/*-------------------------------------------------------- print header */
fprintf(out,"================================================================================\n");
fprintf(out,"Converged Solution in step %d\n",step); 
fprintf(out,"================================================================================\n");
fprintf(out,"________________________________________________________________________________\n\n");
for (j=0; j<actfield->dis[1].numnp; j++)
{
   actnode = &(actfield->dis[1].node[j]);
   fprintf(out,"NODE glob_Id %6d loc_Id %6d    ",actnode->Id,actnode->Id_loc);
   for (k=0; k<actnode->numdf; k++) 
   {
      if (place >= actnode->sol.fdim) dserror("Cannot print solution step");
      fprintf(out,"%20.7E %20.7E %20.7E",actnode->sol.a.da[place][k],actnode->sol.a.da[place][k+2],actnode->sol.a.da[place][k+1]);
   }
   fprintf(out,"\n");
}
fprintf(out,"________________________________________________________________________________\n\n");
for (j=0; j<actfield->dis[1].numnp; j++)
{
   actnode  = &(actfield->dis[1].node[j]);
   actnode2 = &(actfield->dis[0].node[j]);
   actgnode2 = actnode2->gnode;      
   if (FABS(actnode2->x[0]-dynvar->coord_scale[0])<EPS7)
   {
    fprintf(out,"COORD X %20.3E  Y+ %20.7E ",actnode2->x[0],dynvar->washvel*actnode->x[1]/visc);
    for (k=0; k<actnode->numdf; k++) 
    {
     if (place >= actnode->sol.fdim) dserror("Cannot print solution step");
     fprintf(out,"%20.7E %20.7E %20.7E",actnode->sol.a.da[place][k],actnode->sol.a.da[place][k+2],actnode2->sol.a.da[place][k]/dynvar->washvel);
    }
     fprintf(out,"\n");
   }
}
fprintf(out,"________________________________________________________________________________\n\n");
fprintf(out,"________________________________________________________________________________\n\n");
for (j=0; j<actfield->dis[0].numnp; j++)
{
   actnode2 = &(actfield->dis[0].node[j]);
   actgnode2 = actnode2->gnode;      
   if(actgnode2->dirich->dirich_onoff.a.iv[0]==1 && actgnode2->dirich->dirich_onoff.a.iv[1]==1)
   {
    if(actnode2->sol_increment.a.da[3][0] == 0.0 && actnode2->sol_increment.a.da[3][1] == 0.0)
    {
     fprintf(out,"COORD_X %20.3E C_F %20.7E",actnode2->x[0],actnode2->fluid_varia->c_f_shear);
     actnode2->fluid_varia->c_f_shear = ZERO;
     fprintf(out,"\n");
    }
   }
}
fprintf(out,"________________________________________________________________________________\n\n");

/*----------------------------------------------------------------------*/
} /* end of if (myrank==0 && imyrank==0) */
/*----------------------------------------------------------------------*/
if (myrank==0 && imyrank==0) fflush(out);
/*----------------------------------------------------------------------*/

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_fluidmtu */

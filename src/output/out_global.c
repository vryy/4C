#include "../headers/standardtypes.h"
#include "../shell8/shell8.h"
#include "../wall1/wall1.h"
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
   actintra    = (INTRA*)CALLOC(1,sizeof(INTRA));
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
fprintf(out,"Degrees of Freedom:\n");
for (j=0; j<actfield->dis[0].numnp; j++)
{
actnode = &(actfield->dis[0].node[j]);
switch (actnode->numdf)
{
case 2:
fprintf(out,"NODE glob_Id %6d loc_Id %6d    %6d    %6d\n",
        actnode->Id,actnode->Id_loc,actnode->dof[0],actnode->dof[1]);
case 3:
fprintf(out,"NODE glob_Id %6d loc_Id %6d    %6d    %6d   %6d\n",
        actnode->Id,actnode->Id_loc,actnode->dof[0],actnode->dof[1],actnode->dof[2]);
break;
case 4:
fprintf(out,"NODE glob_Id %6d loc_Id %6d    %6d    %6d   %6d   %6d\n",
        actnode->Id,actnode->Id_loc,actnode->dof[0],actnode->dof[1],actnode->dof[2],actnode->dof[3]);
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

/* .... other stuff */




} /* end of (i=0; i<genprob.numfld; i++) */
/*----------------------------------------------------------------------*/
} /* end of if (myrank==0) */
/*----------------------------------------------------------------------*/
fflush(out);
#ifndef PARALLEL 
FREE(actintra);
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
FILE      *out = allfiles.out_out;
NODE      *actnode;
ELEMENT   *actele;
int        myrank;
int        nprocs;
int        imyrank;
int        inprocs;
int        ngauss;
int	   numnp;
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
if (ioflags.struct_disp_file==1)
for (j=0; j<actfield->dis[0].numnp; j++)
{
   actnode = &(actfield->dis[0].node[j]);
   fprintf(out,"NODE glob_Id %6d loc_Id %6d    ",actnode->Id,actnode->Id_loc);
   for (k=0; k<actnode->numdf; k++) 
   {
      if (place >= actnode->sol.fdim) dserror("Cannot print solution step");
      fprintf(out,"%20.7E ",actnode->sol.a.da[place][k]);
   }
   fprintf(out,"\n");
}
fprintf(out,"________________________________________________________________________________\n\n");
if (ioflags.fluid_sol_file==1)
for (j=0; j<actfield->dis[0].numnp; j++)
{
   actnode = &(actfield->dis[0].node[j]);
   fprintf(out,"NODE glob_Id %6d loc_Id %6d    ",actnode->Id,actnode->Id_loc);
   for (k=0; k<actnode->numdf; k++) 
   {
      if (place >= actnode->sol.fdim) dserror("Cannot print solution step");
      fprintf(out,"%20.7E ",actnode->sol.a.da[place][k]);
   }
   fprintf(out,"\n");
}
fprintf(out,"________________________________________________________________________________\n\n");
/*------------------------------------------------ print element values */
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
#endif
   break;
/*---------------------------------------------------------fh 03/02-------*/  
   case el_wall1:
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
   
   break;
   
   default:
      dserror("unknown type of element");
   break;
   }
}
fprintf(out,"\n");
fprintf(out,"\n");

for (j=0; j<actfield->dis[0].numele; j++)
{
   actele = &(actfield->dis[0].element[j]);
   switch(actele->eltyp)
   {
   case el_wall1:
       
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
       
             
   break;
   
   }
}

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

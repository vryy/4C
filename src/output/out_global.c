#include "../headers/standardtypes.h"
#include "../shell8/shell8.h"
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
fprintf(out,"Total Number of Steps     : %d\n",dyn->nstep);
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
   actintra    = (INTRA*)calloc(1,sizeof(INTRA));
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
fprintf(out,"Number of Elements  in this field : %d\n",actfield->numele);
fprintf(out,"Number of Nodes     in this field : %d\n",actfield->numnp);
fprintf(out,"Number of Dofs      in this field : %d\n",actfield->numdf);
fprintf(out,"Number of Equations in this field : %d\n",actfield->numeq);
fprintf(out,"________________________________________________________________________________\n\n");
fprintf(out,"Element connectivity in global Ids:\n");
for (j=0; j<actfield->numele; j++)
{
actele = &(actfield->element[j]);
fprintf(out,"glob_Id %6d Nnodes %2d Nodes: ",actele->Id,actele->numnp);
for (k=0; k<actele->numnp; k++) fprintf(out,"%6d ",actele->node[k]->Id);
fprintf(out,"\n");
}
fprintf(out,"Element connectivity in field-local Ids:\n");
for (j=0; j<actfield->numele; j++)
{
actele = &(actfield->element[j]);
fprintf(out,"loc_Id %6d Nnodes %2d Nodes: ",actele->Id_loc,actele->numnp);
for (k=0; k<actele->numnp; k++) fprintf(out,"%6d ",actele->node[k]->Id_loc);
fprintf(out,"\n");
}
fprintf(out,"________________________________________________________________________________\n\n");
fprintf(out,"Element types:\n");
for (j=0; j<actfield->numele; j++)
{
actele = &(actfield->element[j]);
switch(actele->eltyp)
{
case el_shell8:
fprintf(out,"ELE glob_Id %6d loc_Id %6d SHELL8\n",actele->Id,actele->Id_loc);
break;
case el_brick1:
fprintf(out,"ELE glob_Id %6d loc_Id %6d BRICK1\n",actele->Id,actele->Id_loc);
break;
case el_fluid3:
fprintf(out,"ELE glob_Id %6d loc_Id %6d FLUID3\n",actele->Id,actele->Id_loc);
break;
case el_ale:
fprintf(out,"ELE glob_Id %6d loc_Id %6d ALE\n",actele->Id,actele->Id_loc);
break;
default:
dserror("Cannot print elementtype");
break;
}
}
fprintf(out,"________________________________________________________________________________\n\n");
fprintf(out,"Nodal Coordinates:\n");
for (j=0; j<actfield->numnp; j++)
{
actnode = &(actfield->node[j]);
fprintf(out,"NODE glob_Id %6d loc_Id %6d    %-18.5#f %-18.5#f %-18.5#f \n",
        actnode->Id,actnode->Id_loc,actnode->x[0],actnode->x[1],actnode->x[2]);
}
fprintf(out,"________________________________________________________________________________\n\n");

/* .... other stuff */




} /* end of (i=0; i<genprob.numfld; i++) */
/*----------------------------------------------------------------------*/
} /* end of if (myrank==0) */
/*----------------------------------------------------------------------*/
fflush(out);
#ifndef PARALLEL 
free(actintra);
#endif
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_general */



/*----------------------------------------------------------------------*
 |  print out solution of a certain step                     m.gee 12/01|
 *----------------------------------------------------------------------*/
void out_sol(FIELD *actfield, PARTITION *actpart, INTRA *actintra, int step)
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
for (j=0; j<actfield->numnp; j++)
{
   actnode = &(actfield->node[j]);
   fprintf(out,"NODE glob_Id %6d loc_Id %6d    ",actnode->Id,actnode->Id_loc);
   for (k=0; k<actnode->numdf; k++) 
   {
      if (step >= actnode->sol.fdim) dserror("Cannot print solution step");
      fprintf(out,"%20.7#E ",actnode->sol.a.da[step][k]);
   }
   fprintf(out,"\n");
}
fprintf(out,"________________________________________________________________________________\n\n");
/*------------------------------------------------ print element values */
if (ioflags.struct_stress_file==1)
for (j=0; j<actfield->numele; j++)
{
   actele = &(actfield->element[j]);
   switch(actele->eltyp)
   {
   case el_shell8:
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
       fprintf(out,"Gauss %d   %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E \n",
       i,
       actele->e.s8->forces.a.d3[step][0][i],
       actele->e.s8->forces.a.d3[step][2][i],
       actele->e.s8->forces.a.d3[step][8][i],
       actele->e.s8->forces.a.d3[step][1][i],
       actele->e.s8->forces.a.d3[step][3][i],
       actele->e.s8->forces.a.d3[step][16][i],
       actele->e.s8->forces.a.d3[step][4][i],
       actele->e.s8->forces.a.d3[step][17][i],
       actele->e.s8->forces.a.d3[step][9][i]
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
       fprintf(out,"Gauss %d   %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E \n",
       i,
       actele->e.s8->forces.a.d3[step][5][i],
       actele->e.s8->forces.a.d3[step][7][i],
       actele->e.s8->forces.a.d3[step][14][i],
       actele->e.s8->forces.a.d3[step][6][i],
       actele->e.s8->forces.a.d3[step][10][i],
       actele->e.s8->forces.a.d3[step][12][i],
       actele->e.s8->forces.a.d3[step][11][i],
       actele->e.s8->forces.a.d3[step][13][i],
       actele->e.s8->forces.a.d3[step][15][i]
       );
       }
       
   break;
   default:
      dserror("unknown type of element");
   break;
   }
}
/*----------------------------------------------------------------------*/
} /* end of if (myrank==0 && imyrank==0) */
/*----------------------------------------------------------------------*/
if (myrank==0 && imyrank==0) fflush(out);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_general */

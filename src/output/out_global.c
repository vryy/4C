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

#ifdef D_SHELL8
  #include "../shell8/shell8.h"
#endif /*D_SHELL8*/
#ifdef D_SHELL9
  #include "../shell9/shell9.h"
#endif /*D_SHELL9*/
#ifdef D_WALL1
  #include "../wall1/wall1.h"
#endif /*D_WALL1*/
#ifdef D_BEAM3
  #include "../beam3/beam3.h"
#endif /*D_BEAM3*/
#ifdef D_BRICK1
  #include "../brick1/brick1.h"
#endif /*D_BRICK1*/
#ifdef D_INTERF
  #include "../interf/interf.h"
#endif /*D_INTERF*/
#ifdef D_WALLGE
  #include "../wallge/wallge.h"
#endif /*D_WALLGE*/
#ifdef D_THERM2
  #include "../therm2/therm2.h"
#endif
#ifdef D_THERM3
  #include "../therm3/therm3.h"
#endif
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

#ifndef NO_TEXT_OUTPUT


  INT        i;

#ifdef DEBUG
  INT        j,k;
  ELEMENT   *actele;
  NODE      *actnode;
#endif

#ifdef D_SHELL9
  INT        is_shell9;
#endif

  FILE      *out = allfiles.out_out;
  FIELD     *actfield;
  INTRA     *actintra = NULL;
  INT        myrank;
  INT        nprocs;

#ifdef DEBUG
  dstrc_enter("out_general");
#endif

  /*----------------------------------------------------------------------*/
  myrank = par.myrank;
  nprocs = par.nprocs;

  if (myrank==0)
  {

    /* print header */
    fprintf(out,"CCARAT outputfile\n");
    fprintf(out,"________________________________________________________________________________\n\n");


    /* print problem title */
    for (i=0; i<5; i++)
      fprintf(out,"%s\n",allfiles.title[i]);
    fprintf(out,"________________________________________________________________________________\n\n");



    /* print general problem data */
    fprintf(out,"Total number of Fields    : %d\n",genprob.numfld);
    fprintf(out,"Total number of Elements  : %d\n",genprob.nele);
    fprintf(out,"Total number of Nodes     : %d\n",genprob.nnode);
    fprintf(out,"Total number of Materials : %d\n",genprob.nmat);

    switch(genprob.probtyp)
    {
      case prb_fsi:
        fprintf(out,"Type of Problem           : Fluid-Structure-Interaction\n");
        break;
      case prb_ssi:
        fprintf(out,"Type of Problem           : Structure-Structure-Interaction\n");
        break;
      case prb_structure:
        fprintf(out,"Type of Problem           : Structural\n");
        break;
      case prb_fluid:
        fprintf(out,"Type of Problem           : Fluid\n");
        break;
      case prb_fluid_pm:
        fprintf(out,"Type of Problem           : Fluid Projection\n");
        break;
      case prb_opt:
        fprintf(out,"Type of Problem           : Optimization\n");
        break;
      case prb_ale:
        fprintf(out,"Type of Problem           : Ale\n");
        break;
#ifdef D_TSI
      case prb_tsi:
        fprintf(out,"Type of Problem           : Thermal-Structure-Interaction\n");
        break;
#endif
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
        break;
      default:
        dserror("Cannot print time type");
        break;
    }
    fprintf(out,"________________________________________________________________________________\n\n");


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

      switch(actfield->fieldtyp)
      {
        case fluid:
          fprintf(out,"FIELD: fluid\n");
          fprintf(out,"============\n");
          break;
        case ale:
          fprintf(out,"FIELD: ale\n");
          fprintf(out,"==========\n");
          break;
        case structure:
          fprintf(out,"FIELD: structure\n");
          fprintf(out,"================\n");
          break;
#ifdef D_TSI
        case thermal:
          fprintf(out,"FIELD: thermal\n");
          fprintf(out,"================\n");
          break;
#endif
        default:
          dserror("Cannot print fieldtype");
          break;
      }

      fprintf(out,"   Number of Elements : %d\n",actfield->dis[0].numele);
      fprintf(out,"   Number of Nodes    : %d\n",actfield->dis[0].numnp);
      fprintf(out,"   Number of Dofs     : %d\n",actfield->dis[0].numdf);
      fprintf(out,"   Number of Equations: %d\n",actfield->dis[0].numeq);



#ifdef DEBUG

      fprintf(out,"________________________________________________________________________________\n\n");

      fprintf(out,"Element connectivity in global Ids:\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
        actele = &(actfield->dis[0].element[j]);
        fprintf(out,"glob_Id %6d Nnodes %2d Nodes: ",actele->Id,actele->numnp);
        for (k=0; k<actele->numnp; k++)
          fprintf(out,"%6d ",actele->node[k]->Id);
        fprintf(out,"\n");
      }

      fprintf(out,"Element connectivity in field-local Ids:\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
        actele = &(actfield->dis[0].element[j]);
        fprintf(out,"loc_Id %6d Nnodes %2d Nodes: ",actele->Id_loc,actele->numnp);
        for (k=0; k<actele->numnp; k++)
          fprintf(out,"%6d ",actele->node[k]->Id_loc);
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
          case el_fluid3_fast:
            fprintf(out,"ELE glob_Id %6d loc_Id %6d FLUID3_FAST\n",actele->Id,actele->Id_loc);
            break;
          case el_fluid2:
            fprintf(out,"ELE glob_Id %6d loc_Id %6d FLUID2\n",actele->Id,actele->Id_loc);
            break;
          case el_fluid2_pro:
            fprintf(out,"ELE glob_Id %6d loc_Id %6d FLUID2_PRO\n",actele->Id,actele->Id_loc);
            break;
          case el_fluid3_pro:
            fprintf(out,"ELE glob_Id %6d loc_Id %6d FLUID3_PRO\n",actele->Id,actele->Id_loc);
            break;
          case el_ale3:
            fprintf(out,"ELE glob_Id %6d loc_Id %6d ALE3\n",actele->Id,actele->Id_loc);
            break;
          case el_ale2:
            fprintf(out,"ELE glob_Id %6d loc_Id %6d ALE2\n",actele->Id,actele->Id_loc);
            break;
          case el_beam3:
            fprintf(out,"ELE glob_Id %6d loc_Id %6d BEAM3\n",actele->Id,actele->Id_loc);
            break;
          case el_axishell:
            fprintf(out,"ELE glob_Id %6d loc_Id %6d AXISHELL\n",actele->Id,actele->Id_loc);
            break;
          case el_interf:
            fprintf(out,"ELE glob_Id %6d loc_Id %6d INTERFACE\n",actele->Id,actele->Id_loc);
            break;
          case el_wallge:
            fprintf(out,"ELE glob_Id %6d loc_Id %6d WALLGE\n",actele->Id,actele->Id_loc);
            break;
#ifdef D_THERM2
          case el_therm2:
            fprintf(out,"ELE glob_Id %6d loc_Id %6d THERM2\n",actele->Id,actele->Id_loc);
            break;
#endif
#ifdef D_THERM3
          case el_therm3:
            fprintf(out,"ELE glob_Id %6d loc_Id %6d THERM3\n",actele->Id,actele->Id_loc);
            break;
#endif
          default:
            dserror("Cannot print elementtype");
            break;
        }
      }  /* for (j=0; j<actfield->dis[0].numele; j++) */
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

        /* check if actnode belongs to a shell9 element */
#ifdef D_SHELL9
        is_shell9 = 0;
        for (k=0; k<actnode->numele; k++)
        {
          if (actnode->element[k]->eltyp == el_shell9) is_shell9 = 1;
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
#ifdef D_TSI
          case 1:
            fprintf(out,"NODE glob_Id %6d loc_Id %6d    %6d\n",
                actnode->Id,actnode->Id_loc,actnode->dof[0]);
            break;
#endif
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
          case 7:
            fprintf(out,"NODE glob_Id %6d loc_Id %6d    %6d    %6d   %6d   %6d   %6d   %6d   %6d\n",
                actnode->Id,actnode->Id_loc,actnode->dof[0],actnode->dof[1],actnode->dof[2],actnode->dof[3],actnode->dof[4],actnode->dof[5],actnode->dof[6]);
            break;
          default:
            fprintf(out,"NODE glob_Id %6d loc_Id %6d    No output for actual numdf = %6d \n",
                actnode->Id,actnode->Id_loc,actnode->numdf);
            break;
        }
      }  /* for (j=0; j<actfield->dis[0].numnp; j++) */
      fprintf(out,"________________________________________________________________________________\n\n");

#endif /*ifdef DEBUG */


    fprintf(out,"\n");


    } /* end of (i=0; i<genprob.numfld; i++) */

    fprintf(out,"________________________________________________________________________________\n\n");
    fprintf(out,"\n\n");

  } /* end of if (myrank==0) */


  fflush(out);

#ifndef PARALLEL
  CCAFREE(actintra);
#endif

#ifdef DEBUG
  dstrc_exit();
#endif

  return;

#endif /* NO_TEXT_OUTPUT */

} /* end of out_general */



/*----------------------------------------------------------------------*
 |  print out solution of a certain step                     m.gee 12/01|
 *----------------------------------------------------------------------*/
void out_sol(
    FIELD              *actfield,
    PARTITION          *actpart,
    INT                 disnum,
    INTRA              *actintra,
    INT                 step,
    INT                 place
    )
{
#ifndef NO_TEXT_OUTPUT
INT        i,j,k;
#ifdef D_SHELL9
INT        is_shell9;    /*->shell9*/
#endif
FILE      *out = allfiles.out_out;
NODE      *actnode;
ELEMENT   *actele;
INT        myrank;
INT        nprocs;
INT        imyrank;
INT        inprocs;
INT        ngauss;
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
#ifdef D_TSI
case thermal:
fprintf(out,"FIELD: thermal\n");
break;
#endif
default:
dserror("Cannot print fieldtype");
break;
}
/*-------------------------------------------------- print nodal values */
switch(actfield->fieldtyp)
{
case structure:
   if (ioflags.struct_disp==1)
   {
   fprintf(out,"================================================================================\n");
   fprintf(out,"Converged Solution of Discretisation %d in step %d \n",0,step);
   fprintf(out,"================================================================================\n");
   for (j=0; j<actfield->dis[disnum].numnp; j++)
   {
      actnode = &(actfield->dis[disnum].node[j]);

      /* check if actnode belongs to a shell9 element */
      #ifdef D_SHELL9
         is_shell9 = 0;
         for (k=0; k<actnode->numele; k++)
         {
            if (actnode->element[k]->eltyp == el_shell9) is_shell9 = 1;
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
         fprintf(out,"%20.7E ",actnode->sol.a.da[place][k]);
      }
      fprintf(out,"\n");
   }
   fprintf(out,"________________________________________________________________________________\n\n");
   }
   break;


case fluid:
  if (ioflags.output_out==1 && ioflags.fluid_sol==1)
  {
    fprintf(out,"================================================================================\n");
    fprintf(out,"Converged Solution of Discretisation %d in step %d \n",disnum,step);
    fprintf(out,"================================================================================\n");
    for (j=0; j<actfield->dis[disnum].numnp; j++)
    {
      actnode = &(actfield->dis[disnum].node[j]);
      fprintf(out,"NODE glob_Id %6d loc_Id %6d    ",actnode->Id,actnode->Id_loc);
      for (k=0; k<actnode->sol.sdim; k++)
      {
        if (place >= actnode->sol.fdim) dserror("Cannot print solution step");
        fprintf(out,"%20.7E ",actnode->sol.a.da[place][k]);
      }
      fprintf(out,"\n");
    }
    fprintf(out,"________________________________________________________________________________\n\n");
  }
  break;


case ale:
   if (ioflags.ale_disp==1)
   {
   fprintf(out,"================================================================================\n");
   fprintf(out,"Converged Solution of Discretisation %d in step %d \n",0,step);
   fprintf(out,"================================================================================\n");
   for (j=0; j<actfield->dis[disnum].numnp; j++)
   {
      actnode = &(actfield->dis[disnum].node[j]);
      fprintf(out,"NODE glob_Id %6d loc_Id %6d    ",actnode->Id,actnode->Id_loc);
      for (k=0; k<actnode->numdf; k++)
      {
   	 if (place >= actnode->sol.fdim) dserror("Cannot print solution step");
   	 fprintf(out,"%20.7E ",actnode->sol.a.da[place][k]);
      }
      fprintf(out,"\n");
   }
   fprintf(out,"________________________________________________________________________________\n\n");
   }
break;


case thermal:
  if (ioflags.output_out==1 && ioflags.therm_temper==1)
  {
    fprintf(out,"================================================================================\n");
    fprintf(out,"Converged Solution of Discretisation %d in step %d \n",disnum,step);
    fprintf(out,"================================================================================\n");
    for (j=0; j<actfield->dis[disnum].numnp; j++)
    {
      actnode = &(actfield->dis[disnum].node[j]);
      fprintf(out,"NODE glob_Id %6d loc_Id %6d    ",actnode->Id,actnode->Id_loc);
      for (k=0; k<actnode->numdf; k++)
      {
        if (place >= actnode->sol.fdim) dserror("Cannot print solution step");
        fprintf(out,"%20.7E ",actnode->sol.a.da[place][k]);
      }
      fprintf(out,"\n");
    }
    fprintf(out,"________________________________________________________________________________\n\n");
  }
  break;
default: dserror("Cannot print fieldtype");
}
/*------------------------------------------------ print element values */
switch(actfield->fieldtyp)
{
case structure:
if (ioflags.struct_stress==1)
for (j=0; j<actfield->dis[disnum].numele; j++)
{
   actele = &(actfield->dis[disnum].element[j]);
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

#if 0
       ngauss   = actele->e.s9->nGP[0] * actele->e.s9->nGP[1];
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
          {
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
          {
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
#endif

#endif /*D_SHELL9*/
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
   case el_brick1:
#ifdef D_BRICK1
       ngauss =  actele->e.c1->nGP[0]
               * actele->e.c1->nGP[1]
               * actele->e.c1->nGP[2];
       fprintf(out,"________________________________________________________________________________\n");
       fprintf(out,"Element glob_Id %d loc_Id %d                BRICK1\n",actele->Id,actele->Id_loc);
       fprintf(out,"\n");
       switch(actele->e.c1->stresstyp)
       {
       case c1_gpxyz:
       fprintf(out,"INT.point   x-coord.     y-coord.     z-coord.     stress-xx    stress-yy    stress-zz    stress-xy    stress-xz    stress-yz\n");
       for (i=0; i<ngauss; i++)
       {
          fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
          i,
          actele->e.c1->stress_GP.a.d3[place][24][i],
          actele->e.c1->stress_GP.a.d3[place][25][i],
          actele->e.c1->stress_GP.a.d3[place][26][i],
          actele->e.c1->stress_GP.a.d3[place][ 6][i],
          actele->e.c1->stress_GP.a.d3[place][ 7][i],
          actele->e.c1->stress_GP.a.d3[place][ 8][i],
          actele->e.c1->stress_GP.a.d3[place][ 9][i],
          actele->e.c1->stress_GP.a.d3[place][10][i],
          actele->e.c1->stress_GP.a.d3[place][11][i]
          );
       }
       break;
       case c1_gprst:
       fprintf(out,"r,s,t    ---> local system on element level \n");
       fprintf(out,"rr,ss,tt ---> normal-stresses               \n");
       fprintf(out,"rs,st,tr ---> shear -stresses               \n\n");
       fprintf(out,"INT.point   x-coord.     y-coord.     z-coord.     stress-rr    stress-ss    stress-tt    stress-rs    stress-st    stress-tr\n");
       for (i=0; i<ngauss; i++)
       {
          fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
          i,
          actele->e.c1->stress_GP.a.d3[place][24][i],
          actele->e.c1->stress_GP.a.d3[place][25][i],
          actele->e.c1->stress_GP.a.d3[place][26][i],
          actele->e.c1->stress_GP.a.d3[place][0][i],
          actele->e.c1->stress_GP.a.d3[place][1][i],
          actele->e.c1->stress_GP.a.d3[place][2][i],
          actele->e.c1->stress_GP.a.d3[place][3][i],
          actele->e.c1->stress_GP.a.d3[place][4][i],
          actele->e.c1->stress_GP.a.d3[place][5][i]
          );
       }
       break;
       case c1_gp123:
       fprintf(out,"11,22,33 ---> principal-stresses                       \n");
       fprintf(out,"r1,s1,t1 ---> angles to the first  principal direction \n");
       fprintf(out,"r2,s2,t2 ---> angles to the second principal direction \n");
       fprintf(out,"r3,s3,t3 ---> angles to the third  principal direction \n\n");
       fprintf(out,"INT.point   stress-11    stress-22    stress-33  ang-r1  ang-s1   ang-t1    ang-r2   ang-s2   ang-t2   ang-r3   ang-s3   ang-t3\n");
       for (i=0; i<ngauss; i++)
       {
          fprintf(out,"  %-6d %12.3E %12.3E %12.3E %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f \n",
          i,
          actele->e.c1->stress_GP.a.d3[place][12][i],
          actele->e.c1->stress_GP.a.d3[place][13][i],
          actele->e.c1->stress_GP.a.d3[place][14][i] ,
          actele->e.c1->stress_GP.a.d3[place][15][i],
          actele->e.c1->stress_GP.a.d3[place][16][i],
          actele->e.c1->stress_GP.a.d3[place][17][i],
          actele->e.c1->stress_GP.a.d3[place][18][i],
          actele->e.c1->stress_GP.a.d3[place][19][i],
          actele->e.c1->stress_GP.a.d3[place][20][i],
          actele->e.c1->stress_GP.a.d3[place][21][i],
          actele->e.c1->stress_GP.a.d3[place][22][i],
          actele->e.c1->stress_GP.a.d3[place][23][i]
          );
       }
       break;
       case c1_nprst:
       fprintf(out,"elenode     stress-rr    stress-ss    stress-tt    stress-rs    stress-st    stress-tr\n");
       for (i=0; i<actele->numnp; i++)
       {
          fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
          i,
          actele->e.c1->stress_ND.a.d3[place][0][i],
          actele->e.c1->stress_ND.a.d3[place][1][i],
          actele->e.c1->stress_ND.a.d3[place][2][i] ,
          actele->e.c1->stress_ND.a.d3[place][3][i],
          actele->e.c1->stress_ND.a.d3[place][4][i],
          actele->e.c1->stress_ND.a.d3[place][5][i]
          );
       }
       break;
       case c1_np123:
       fprintf(out,"elenode     stress-11    stress-22    stress-33  ang-r1  ang-s1   ang-t1    ang-r2   ang-s2   ang-t2   ang-r3   ang-s3   ang-t3\n");
       for (i=0; i<actele->numnp; i++)
       {
          fprintf(out,"  %-6d %12.3E %12.3E %12.3E %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f \n",
          i,
          actele->e.c1->stress_ND.a.d3[place][12][i],
          actele->e.c1->stress_ND.a.d3[place][13][i],
          actele->e.c1->stress_ND.a.d3[place][14][i] ,
          actele->e.c1->stress_ND.a.d3[place][15][i],
          actele->e.c1->stress_ND.a.d3[place][16][i],
          actele->e.c1->stress_ND.a.d3[place][17][i],
          actele->e.c1->stress_ND.a.d3[place][18][i],
          actele->e.c1->stress_ND.a.d3[place][19][i],
          actele->e.c1->stress_ND.a.d3[place][20][i],
          actele->e.c1->stress_ND.a.d3[place][21][i],
          actele->e.c1->stress_ND.a.d3[place][22][i],
          actele->e.c1->stress_ND.a.d3[place][23][i]
          );
       }
       break;
       case c1_npxyz:
       fprintf(out,"elenode     stress-xx    stress-yy    stress-zz    stress-xy    stress-yz    stress-xz\n");
       for (i=0; i<actele->numnp; i++)
       {
          fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
          i,
          actele->e.c1->stress_ND.a.d3[place][ 6][i],
          actele->e.c1->stress_ND.a.d3[place][ 7][i],
          actele->e.c1->stress_ND.a.d3[place][ 8][i] ,
          actele->e.c1->stress_ND.a.d3[place][ 9][i],
          actele->e.c1->stress_ND.a.d3[place][10][i],
          actele->e.c1->stress_ND.a.d3[place][11][i]
          );
       }
       break;
       default:
       fprintf(out,"no stresses available\n");
       }
#endif /*D_BRICK1*/
   break;

   case el_beam3:
#ifdef D_BEAM3
       ngauss = actele->e.b3->nGP[0];
       fprintf(out,"________________________________________________________________________________\n");
       fprintf(out,"Element glob_Id %d loc_Id %d                BEAM3\n",actele->Id,actele->Id_loc);
       fprintf(out,"\n");
       fprintf(out,"Gaussian         Nx           Vy           Vz           Mx           My           Mz\n");
       for (i=0; i<ngauss; i++)
       {
       fprintf(out,"Gauss %d       %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
       i,
       actele->e.b3->force_GP.a.d3[place][0][i],
       actele->e.b3->force_GP.a.d3[place][1][i],
       actele->e.b3->force_GP.a.d3[place][2][i],
       actele->e.b3->force_GP.a.d3[place][3][i],
       actele->e.b3->force_GP.a.d3[place][4][i],
       actele->e.b3->force_GP.a.d3[place][5][i]

       );
       }
#endif /*D_BEAM3*/
   break;

   case el_interf:
#ifdef D_INTERF
       ngauss = 2;
       fprintf(out,"________________________________________________________________________________\n");
       fprintf(out,"Element glob_Id %d loc_Id %d                INTERF\n",actele->Id,actele->Id_loc);

       switch(actele->e.interf->stresstyp)
       {
       case if_xy:
       fprintf(out,"Gaussian     stresses-xx     stresses-yy    stresses-xy\n");
       break;
       case if_tn:
       fprintf(out,"Gaussian     stresses-tangential     stresses-normal\n");
       break;
       default:
          dserror("Unknown type of element stresses");
       }
       for (i=0; i<ngauss; i++)
       {
       fprintf(out,"Gauss %d   %12.3E %12.3E %12.3E\n",
       i,
       actele->e.interf->stress_GP.a.d3[place][0][i],
       actele->e.interf->stress_GP.a.d3[place][1][i],
       actele->e.interf->stress_GP.a.d3[place][2][i] );
       }
#endif /*D_INTERF*/
   break;

   case el_wallge:
#ifdef D_WALLGE
       ngauss = actele->e.wallge->nGP[0] * actele->e.wallge->nGP[1];
       fprintf(out,"________________________________________________________________________________\n");
       fprintf(out,"Element glob_Id %d loc_Id %d                WALL1\n",actele->Id,actele->Id_loc);
       fprintf(out,"\n");
/* check wether stresses at Gauss Points are presented in global xy- or local rs-coordinate system and write stress type */
       switch(actele->e.wallge->stresstyp)
       {
       case wge_xy:
       fprintf(out,"Gaussian     Stress-xx    Stress-yy    Stress-xy    Stress-zz    Max. P.S.    Min. P.S.    Angle\n");
       break;
       case wge_rs:
       fprintf(out,"Gaussian     Stress-rr    Stress-ss    Stress-rs    Stress-zz    Max. P.S.    Min. P.S.    Angle\n");
       break;
       default:
       fprintf(out,"Gaussian     Stress-xx    Stress-yy    Stress-xy    Stress-zz    Max. P.S.    Min. P.S.    Angle\n");
       }
       for (i=0; i<ngauss; i++)
       {
       fprintf(out,"Gauss %d   %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
       i,
       actele->e.wallge->stress_GP.a.d3[place][0][i],
       actele->e.wallge->stress_GP.a.d3[place][1][i],
       actele->e.wallge->stress_GP.a.d3[place][2][i],
       actele->e.wallge->stress_GP.a.d3[place][3][i],
       actele->e.wallge->stress_GP.a.d3[place][4][i],
       actele->e.wallge->stress_GP.a.d3[place][5][i],
       actele->e.wallge->stress_GP.a.d3[place][6][i]
       );
       }
#endif /*D_WALLGE*/
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
#ifdef D_TSI
case thermal:
  if (ioflags.therm_heatflux == 1)
  {
    for (j=0; j<actfield->dis[disnum].numele; j++)
    {
      actele = &(actfield->dis[disnum].element[j]);
      switch (actele->eltyp)
      {
        case el_therm2:
#ifdef D_THERM2
          /* print title */
          fprintf(out, "_________________________________________________"
                  "_______________________________\n");
          fprintf(out, "Element glob_Id %d loc_Id %d"
                  "                THERM2\n", actele->Id, actele->Id_loc);
          fprintf(out, "\n");
          /* check heat flux type */
          switch (actele->e.th2->hfluxtype)
          {
            /* none */
            case th2_hflux_none:
              break;
            /* at Gauss points */
            case th2_hflux_gpxy:
            case th2_hflux_gp12:
            case th2_hflux_gpxy12:
              fprintf(out, "Gaussian   "
                      "  heatflux-x   heatflux-y   heatflux-z"
                      "   abs.hflux        angle   heatflux-z\n");
              /* loop all Gauss points */
              for (i=0; i<actele->e.th2->gptotal; i++)
              {
                fprintf(out, "Gauss %2d   "
                        "%12.3E %12.3E %12.3E"
                        "%12.3E %12.3E %12.3E \n",
                        i,
                        actele->e.th2->hflux_gp.a.d3[place][0][i],
                        actele->e.th2->hflux_gp.a.d3[place][1][i],
                        actele->e.th2->hflux_gp.a.d3[place][2][i],
                        actele->e.th2->hflux_gp.a.d3[place][3][i],
                        actele->e.th2->hflux_gp.a.d3[place][4][i],
                        actele->e.th2->hflux_gp.a.d3[place][5][i]);
              }
              break;
            /* at element nodes */
            case th2_hflux_ndxy:
            case th2_hflux_nd12:
            case th2_hflux_ndxy12:
              fprintf(out, "Nodal      "
                      "  heatflux-x   heatflux-y   heatflux-z"
                      "   abs.hflux        angle   heatflux-z\n");
              /* loop all element nodes */
              for (i=0; i<actele->numnp; i++)
              {
                fprintf(out, "Node  %2d   "
                        "%12.3E %12.3E %12.3E"
                        "%12.3E %12.3E %12.3E \n",
                        i,
                        actele->e.th2->hflux_nd.a.d3[place][0][i],
                        actele->e.th2->hflux_nd.a.d3[place][1][i],
                        actele->e.th2->hflux_nd.a.d3[place][2][i],
                        actele->e.th2->hflux_nd.a.d3[place][3][i],
                        actele->e.th2->hflux_nd.a.d3[place][4][i],
                        actele->e.th2->hflux_nd.a.d3[place][5][i]);
              }
              break;
            default:
              dserror("Heat flux type is not available for printing!");
          }  /* end of */
#endif  /* end of #ifdef D_THERM2 */
          break;
        default:
          dserror("Unknown element type!");
          break;
      }  /* end of switch (actele->eltyp) */
    }  /* end of for (j=0; j<actfield->dis[disnum].numele; j++) */
  }  /* end of if (ioflags.therm_heatflux == 1) */
#endif  /* end of #ifdef D_TSI */
break;
default:
break;
} /* end switch(actfield->fieldtyp) */

fprintf(out,"\n");
fprintf(out,"\n");

#if 0

for (j=0; j<actfield->dis[disnum].numele; j++)
{
   actele = &(actfield->dis[disnum].element[j]);
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

   case el_beam3:
#ifdef D_BEAM3
       numnp = actele->numnp;
       fprintf(out,"________________________________________________________________________________\n");
       fprintf(out,"Element glob_Id %d loc_Id %d                BEAM3\n",actele->Id,actele->Id_loc);
       fprintf(out,"\n");
       fprintf(out,"Nodal            Nx           Vy           Vz           Mx           My           Mz\n");
       for (i=0; i<numnp; i++)
       {
       fprintf(out,"Node %d  IEL %d %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E %12.3#E \n",
       actele->node[i]->Id,
       i,
       actele->e.b3->force_ND.a.d3[place][0][i],
       actele->e.b3->force_ND.a.d3[place][1][i],
       actele->e.b3->force_ND.a.d3[place][2][i],
       actele->e.b3->force_ND.a.d3[place][3][i],
       actele->e.b3->force_ND.a.d3[place][4][i],
       actele->e.b3->force_ND.a.d3[place][5][i]
       );
       }
#endif /*D_BEAM3*/
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
#endif /* NO_TEXT_OUTPUT */
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

#if !defined(NO_TEXT_OUTPUT) && defined(DEBUG)

#ifdef D_FSI
  INT        i;
  FILE      *out = allfiles.out_out;
  NODE      *actfnode, *actanode, *actsnode;
  GNODE     *actfgnode;
  INT        myrank;
  INT        nprocs;
  INT        numaf,numsf;

#ifdef DEBUG
  dstrc_enter("out_fsi");
#endif

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
#endif /* NO_TEXT_OUTPUT */
} /* end of out_fsi */





/*!---------------------------------------------------------------------
\brief print ssi coupling informations

<pre>                                                         genk 01/03

</pre>

\param *masterfield    FIELD   (i)

\return void

------------------------------------------------------------------------*/
void out_ssi(FIELD *masterfield)
{
#ifndef NO_TEXT_OUTPUT

#ifdef D_SSI
INT        i;
FILE      *out = allfiles.out_out;
NODE      *actmnode, *actsnode;
GNODE     *actmgnode;
INT        myrank;
INT        nprocs;
INT        numsf;

#ifdef DEBUG
dstrc_enter("out_ssi");
#endif

/*----------------------------------------------------------------------*/
myrank = par.myrank;
nprocs = par.nprocs;
numsf  = 1;
/*----------------------------------------------------------------------*/
if (myrank==0)
{
fprintf(out,"================================================================================\n");
fprintf(out,"SSI node connectivity global Ids:\n");
fprintf(out,"================================================================================\n");
fprintf(out,"\n");
fprintf(out,"MASTER          SLAVE\n");
for (i=0;i<masterfield->dis[0].numnp;i++)
{
   actmnode  = &(masterfield->dis[0].node[i]);
   actmgnode = actmnode->gnode;
   actsnode  = actmgnode->mfcpnode[numsf];
   if (actsnode!=NULL)
   fprintf(out,"%-6d         %-6d\n",actmnode->Id,actsnode->Id);
   else
   fprintf(out,"%-6d         ------\n",actmnode->Id);

}
fprintf(out,"________________________________________________________________________________\n\n");
fprintf(out,"================================================================================\n");
fprintf(out,"SSI node connectivity local Ids:\n");
fprintf(out,"================================================================================\n");
fprintf(out,"\n");
fprintf(out,"MASTER          SLAVE\n");
for (i=0;i<masterfield->dis[0].numnp;i++)
{
   actmnode  = &(masterfield->dis[0].node[i]);
   actmgnode = actmnode->gnode;
   actsnode  = actmgnode->mfcpnode[numsf];
   if (actsnode!=NULL)
   fprintf(out,"%-6d         %-6d\n",actmnode->Id_loc,actsnode->Id_loc);
   else
   fprintf(out,"%-6d         ------\n",actmnode->Id_loc);

}

fprintf(out,"________________________________________________________________________________\n\n");
/*----------------------------------------------------------------------*/
} /* end of if (myrank==0 && imyrank==0) */
/*----------------------------------------------------------------------*/
if (myrank==0) fflush(out);
/*----------------------------------------------------------------------*/
#else
dserror("SSI functions not compiled in\n");
#endif
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of out_ssi */


/*----------------------------------------------------------------------*/
/*!
\brief Print TSI coupling informations

\param *structfield    FIELD   (i)
\return void

\author bborn
\date 03/06

*/
/*----------------------------------------------------------------------*/
#ifdef D_TSI
void out_tsi(FIELD *structfield)
{

#if !defined(NO_TEXT_OUTPUT) && defined(DEBUG)

  INT        i;
  INT        idis;  /* number of discretisations in structure field */
  FILE      *out = allfiles.out_out;
  NODE      *actsnode, *acttnode;
  GNODE     *actsgnode;
  ELEMENT   *actsele, *acttele;
  INT        myrank;
  INT        nprocs;
  INT        numtf;

#ifdef DEBUG
  dstrc_enter("out_tsi");
#endif

  /*--------------------------------------------------------------------*/
  myrank = par.myrank;
  nprocs = par.nprocs;
  numtf  = genprob.numtf;
  /*--------------------------------------------------------------------*/
  if (myrank==0)
  {
    /*------------------------------------------------------------------*/
    /* check wether multiple discretisations */
    if (structfield->ndis > 1)
    {
      dserror("Cannot print more than 1 structure discretisation!");
    }
    else
    {
      idis = 0;  /* hello again, this is C, therefore, 0-based index */
    }

    /*------------------------------------------------------------------*/
    /* connectivity of nodes */

    /*------------------------------------------------------------------*/
    /* global connectivity */

    /* print title */
    fprintf(out, "======================================================="
            "=========================\n");
    fprintf(out, "TSI node connectivity global Ids:\n");
    fprintf(out, "======================================================="
            "=========================\n");
    fprintf(out, "\n");
    fprintf(out, "STRUCTURE          THERMAL\n");

    /* loop over structure nodes */
    for (i=0; i<structfield->dis[idis].numnp; i++)
    {
      actsnode  = &(structfield->dis[idis].node[i]);
      actsgnode = actsnode->gnode;
      acttnode  = actsgnode->mfcpnode[numtf];

      if (acttnode != NULL)
      {
        fprintf(out, "%-6d             %-6d\n",
                actsnode->Id, acttnode->Id);
      }
      else
      {
        fprintf(out, "%-6d             ------\n", actsnode->Id);
      }
    }  /* end of for (i=0; i<fluidfield->dis[idis].numnp; i++) */


    /*------------------------------------------------------------------*/
    /* local connectivity */

    /* print title */
    fprintf(out, "_______________________________________________________"
            "_________________________\n\n");
    fprintf(out, "======================================================="
            "=========================\n");
    fprintf(out, "TSI node connectivity local Ids:\n");
    fprintf(out, "======================================================="
            "=========================\n");
    fprintf(out, "\n");
    fprintf(out,"STRUCTURE          THERMAL\n");

    /* loop all structure nodes */
    for (i=0; i<structfield->dis[idis].numnp; i++)
    {
      actsnode = &(structfield->dis[idis].node[i]);
      actsgnode = actsnode->gnode;
      acttnode = actsgnode->mfcpnode[numtf];
      if (acttnode != NULL)
      {
        fprintf(out, "%-6d             %-6d\n",
                actsnode->Id_loc, acttnode->Id_loc);
      }
      else
      {
        fprintf(out, "%-6d             ------\n",
                actsnode->Id_loc);
      }

    }  /* end of for (i=0; i<structfield->dis[idis].numnp; i++) */

    /* print footer */
    fprintf(out, "_______________________________________________________"
            "_________________________\n\n");

    /*------------------------------------------------------------------*/
    /* connectivity of elements */

    /* print title */
    fprintf(out, "======================================================="
            "=========================\n");
    fprintf(out, "TSI element connectivity global/local Ids:\n");
    fprintf(out, "======================================================="
            "=========================\n");
    fprintf(out, "\n");
    fprintf(out, "STRUCTURE glob     "
                 "STRUCTURE loc      "
                 "THERMAL glob       "
                 "THERMAL loc\n");

    /* loop over structure elements */
    for (i=0; i<structfield->dis[idis].numele; i++)
    {
      actsele = &(structfield->dis[idis].element[i]);
      acttele = NULL;
      switch (actsele->eltyp)
      {
#ifdef D_WALL1
          case el_wall1:
            acttele = actsele->e.w1->therm_ele;
            break;
#endif
#ifdef D_BRICK1
          case el_brick1:
            acttele = actsele->e.c1->therm_ele;
            break;
#endif
          default:
            continue;
      }  /* end of switch */

      if (acttele != NULL)
      {
        fprintf(out, "%-6d             "
                     "%-6d             "
                     "%-6d             "
                     "%-6d\n",
                actsele->Id, actsele->Id_loc, acttele->Id, acttele->Id_loc);
      }
      else
      {
        fprintf(out, "%-6d             ------\n", actsele->Id);
      }
    }  /* end of for (i=0; i<fluidfield->dis[idis].numele; i++) */

    /* print footer */
    fprintf(out, "_______________________________________________________"
            "_________________________\n\n");


  } /* end of if (myrank==0 && imyrank==0) */
  /*----------------------------------------------------------------------*/
  if (myrank == 0)
  {
    fflush(out);
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return;

#endif /* end of #ifdef !(NO_TEXT_OUTPUT) && defined(DEBUG) */

} /* end of out_tsi */
#endif /* end of #ifdef D_TSI */



/*!---------------------------------------------------------------------
\brief  print fluid multifield coupling informations

<pre>                                                         genk 01/03

</pre>

\param *fluidfield    FIELD   (i)

\return void

------------------------------------------------------------------------*/
void out_fluidmf(FIELD *fluidfield)
{
#ifndef NO_TEXT_OUTPUT

#ifdef D_FSI
INT        i;
FILE      *out = allfiles.out_out;
NODE      *actfnode, *actanode;
GNODE     *actfgnode;
INT        myrank;
INT        nprocs;
INT        numff;
INT        numaf;

#ifdef DEBUG
dstrc_enter("out_fluidmf");
#endif

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
fprintf(out,"Fluid multiefield node connectivity local Ids:\n");
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
#endif /* NO_TEXT_OUTPUT */
} /* end of out_fluidmf */

/*----------------------------------------------------------------------*
 |  print out solution of turbulence of a certain step      he    05/03 |
 *----------------------------------------------------------------------*/
void out_fluidtu(FIELD *actfield, INTRA *actintra, INT step, INT place)
{
#ifndef NO_TEXT_OUTPUT
INT        i,j,k;
FILE      *out = allfiles.out_tur;
NODE      *actnode;
NODE      *actnode2;
GNODE     *actgnode2;
INT        myrank;
INT        nprocs;
INT        imyrank;
INT        inprocs;
INT	     numff;
DOUBLE     visc = 0.0;
FLUID_DYNAMIC *fdyn;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("out_fluidmtu");
#endif
/*----------------------------------------------------------------------*/
fdyn   = alldyn[genprob.numff].fdyn;

for (numff=0;numff<genprob.numfld;numff++)
{
 actfield=&(field[numff]);
 if (actfield->fieldtyp==fluid)
 break;
} /* end loop over numff */

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
   if (FABS(actnode2->x[0]-fdyn->coord_scale[0])<EPS7)
   {
    fprintf(out,"COORD X %20.3E  Y+ %20.7E ",actnode2->x[0],fdyn->washvel*actnode->x[1]/visc);
    for (k=0; k<actnode->numdf; k++)
    {
     if (place >= actnode->sol.fdim) dserror("Cannot print solution step");
     fprintf(out,"%20.7E %20.7E %20.7E",actnode->sol.a.da[place][k],actnode->sol.a.da[place][k+2],actnode2->sol.a.da[place][k]/fdyn->washvel);
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
   if (actgnode2->dirich==NULL) continue;
   if(actgnode2->dirich->dirich_onoff.a.iv[0]==1 && actgnode2->dirich->dirich_onoff.a.iv[1]==1)
   {
    if(actnode2->sol_increment.a.da[3][0] == 0.0 && actnode2->sol_increment.a.da[3][1] == 0.0)
    {
#ifdef D_FLUID
     fprintf(out,"COORD_X %20.3E C_F %20.7E",actnode2->x[0],actnode2->fluid_varia->c_f_shear);
     actnode2->fluid_varia->c_f_shear = ZERO;
     fprintf(out,"\n");
#endif
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
#endif /* NO_TEXT_OUTPUT */
} /* end of out_fluidmtu */

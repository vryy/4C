/*!----------------------------------------------------------------------
  \file
  \brief controls output of fsi for gid

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

  ------------------------------------------------------------------------*/
#ifndef CCADISCRET
/*!
  \addtogroup FSI
  *//*! @{ (documentation module open)*/
#ifdef D_FSI
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "gid.h"
#include "../shell8/shell8.h"
#include "../shell9/shell9.h"
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
  |                                                       m.gee 06/01    |
  | pointer to allocate dynamic variables if needed                      |
  | dedfined in global_control.c                                         |
  | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;
/*----------------------------------------------------------------------*
  |                                                          mn 06/02    |
  | structure of flags to control output                                 |
  | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;
/*!----------------------------------------------------------------------
  \brief ranks and communicators

  <pre>                                                         m.gee 8/00
  This structure struct _PAR par; is defined in main_ccarat.c
  and the type is in partition.h
  </pre>

 *----------------------------------------------------------------------*/
extern struct _PAR   par;

extern struct _CURVE *curve;
/*!---------------------------------------------------------------------
  \brief routine to control output og fsi to gid

  <pre>                                                           mn 05/03

  This function writes the output of the current timestep into the
  output-file for gid.  The output includes:
  - fluid displacements (mesh deformation)
  - fluid velocities
  - fluid pressure
  - structure displacements
  - structure stresses
  - structure velocities    (not yet implemented)
  - structure accelerations (not yet implemented)
  </pre>

  \param *fluidfield    FIELD   (i)     the fluid field to output
  \param *structfield   FIELD   (i)     the structure field to output

  \return void

  ------------------------------------------------------------------------*/
void out_gid_sol_fsi(
    FIELD       *fluidfield,
    FIELD       *structfield,
    INT          disnumf,
    INT          disnums
    )
{

#ifndef NO_TEXT_OUTPUT

#ifdef SPLIT_HEX20
  INT           j;
  ELEMENT      *actele;
  INT           place;
#endif

  INT           i;

  static FLUID_DYNAMIC  *fdyn;
  static FSI_DYNAMIC    *fsidyn;
  static STRUCT_DYNAMIC *sdyn;

  INTRA        *actintras = NULL;
  INTRA        *actintraf;

  FILE         *out = allfiles.gidres;

  NODE         *actnode;
  NODE         *actanode;


  char         *resulttype;
  char         *resultplace;
  char         *rangetable;
  INT           ncomponent;
  char         *componentnames[18];

  char          sign='"';


#ifdef DEBUG
  dstrc_enter("out_gid_sol_fsi");
#endif


  if (structfield!=NULL)
    sdyn= alldyn[genprob.numsf].sdyn;
  fdyn= alldyn[genprob.numff].fdyn;
  fsidyn= alldyn[genprob.numfld].fsidyn;


#ifdef PARALLEL
  actintraf    = &(par.intra[genprob.numff]);
  if (structfield!=NULL)
  actintras    = &(par.intra[genprob.numsf]);
#else
  actintraf   = (INTRA*)CCACALLOC(1,sizeof(INTRA));
  if (!actintraf) dserror("Allocation of INTRA failed");
  actintraf->intra_fieldtyp = fluid;
  actintraf->intra_rank     = 0;
  actintraf->intra_nprocs   = 1;

  if (structfield!=NULL)
  {
    actintras   = (INTRA*)CCACALLOC(1,sizeof(INTRA));
    if (!actintras) dserror("Allocation of INTRA failed");
    actintras->intra_fieldtyp = structure;
    actintras->intra_rank     = 0;
    actintras->intra_nprocs   = 1;
  }
#endif

  if (structfield!=NULL)
    if ( sdyn->step != fdyn->step)
      dserror("Something is wrong!!");

  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"# Converged Solution of timestep %d\n",fsidyn->step);
  fprintf(out,"#-------------------------------------------------------------------------------\n");



  /* write header for displacements of structure & ale fields */
  /*----------------------------------------------------------*/
  if (ioflags.struct_disp == 1 || ioflags.ale_disp == 1)
  {
    resulttype        = "VECTOR";
    resultplace       = "ONNODES";
    rangetable        = "standard_structure";
    ncomponent        = genprob.ndim;
    componentnames[0] = "x-displ";
    componentnames[1] = "y-displ";
    componentnames[2] = "z-displ";
    /*-------------------------------------------------------------------*/
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"# RESULT DISPLACEMENTS on FIELD FSI, DIS %1i\n",disnumf);
    fprintf(out,"# TIME %18.5E \n",fsidyn->time);
    fprintf(out,"# STEP %6d    \n",fsidyn->step);
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"RESULT %cdisplacement%c %cccarat%c %d %s %s\n",
        sign,sign,
        sign,sign,
        fsidyn->step,
        resulttype,
        resultplace
        );
    fprintf(out,"RESULTRANGESTABLE %c%s%c\n",
        sign,rangetable,sign
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
  }


  /* write values for displ. of structure field */
  /*--------------------------------------------*/
  if (ioflags.struct_disp==1 && structfield!=NULL)
  {
    for (i=0; i<structfield->dis[0].numnp; i++)
    {
      actnode = &(structfield->dis[0].node[i]);
      switch (genprob.ndim)
      {
        case 3:
          fprintf(out," %6d %18.5E %18.5E %18.5E\n",
              actnode->Id+1,
              actnode->sol.a.da[0][0],
              actnode->sol.a.da[0][1],
              actnode->sol.a.da[0][2]
              );
          break;
        case 2:
          fprintf(out," %6d %18.5E %18.5E \n",
              actnode->Id+1,
              actnode->sol.a.da[0][0],
              actnode->sol.a.da[0][1]
              );
          break;
        default:
          dserror("Unknown number of dimensions");
          break;
      }
    }
  }



  /* write values for displacements of fluid (ale) field */
  /*-----------------------------------------------------*/
  if (ioflags.ale_disp == 1)
  {
    for (i=0; i<fluidfield->dis[0].numnp; i++)
    {
      actnode = &(fluidfield->dis[0].node[i]);
      actanode = (fluidfield->dis[0].node[i].gnode->mfcpnode[genprob.numaf]);
      /* ALE Gebiet */
      if (actanode != NULL)
      {
        switch (genprob.ndim)
        {
          case 3:
            fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                actnode->Id+1,
                    actanode->sol.a.da[0 /*fsidyn->actpos*/][0],
                    actanode->sol.a.da[0 /*fsidyn->actpos*/][1],
                    actanode->sol.a.da[0 /*fsidyn->actpos*/][2]
                );
            break;
          case 2:
            fprintf(out," %6d %18.5E %18.5E \n",
                actnode->Id+1,
                    actanode->sol.a.da[0 /*fsidyn->actpos*/][0],
                    actanode->sol.a.da[0 /*fsidyn->actpos*/][1]
                );
            break;
          default:
            dserror("Unknown number of dimensions");
            break;
        }
      }
      /* EULER GEBIET */
      else
      {
        switch (genprob.ndim)
        {
          case 3:
            fprintf(out," %6d %18.5E %18.5E %18.5E\n", actnode->Id+1, 0.0, 0.0, 0.0);
            break;
          case 2:
            fprintf(out," %6d %18.5E %18.5E \n", actnode->Id+1, 0.0, 0.0);
            break;
          default:
            dserror("Unknown number of dimensions");
            break;
        }
      }
    }


#ifdef SPLIT_HEX20
    /* write solutions for additional nodes */
    for (j=0; j<fluidfield->dis[0].numele; j++)
    {
      INT l;
      INT eleid,nodebase;
      DOUBLE x[3];
      actele = &(fluidfield->dis[0].element[j]);
      place  = fsidyn->actpos;

      /*
      ->gnode->mfcpnode[genprob.numaf]
      */

      if ( !(actele->eltyp == el_ale3 || actele->eltyp == el_fluid3_fast)
          || actele->numnp !=20) continue;


      eleid = actele->Id+1;
      nodebase = eleid * 1000;

      /* node 01 */
      for (l=0;l<3;l++)
        x[l] = 0.25*(actele->node[ 8]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
            + actele->node[ 9]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
            + actele->node[10]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
            + actele->node[11]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l] );
      fprintf(out," %6d %18.5e %18.5e %18.5e\n",
          nodebase + 1,
          x[0],
          x[1],
          x[2]);

      /* node 02 */
      for (l=0;l<3;l++)
        x[l] = 0.25*(actele->node[16]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
            + actele->node[17]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
            + actele->node[18]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
            + actele->node[19]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l] );
      fprintf(out," %6d %18.5e %18.5e %18.5e\n",
          nodebase + 2,
          x[0],
          x[1],
          x[2]);

      /* node 03 */
      for (l=0;l<3;l++)
        x[l] = 0.25*(actele->node[ 8]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
            + actele->node[13]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
            + actele->node[16]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
            + actele->node[12]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l] );
      fprintf(out," %6d %18.5e %18.5e %18.5e\n",
          nodebase + 3,
          x[0],
          x[1],
          x[2]);

      /* node 04 */
      for (l=0;l<3;l++)
        x[l] = 0.25*(actele->node[ 9]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
            + actele->node[14]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
            + actele->node[17]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
            + actele->node[13]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l] );
      fprintf(out," %6d %18.5e %18.5e %18.5e\n",
          nodebase + 4,
          x[0],
          x[1],
          x[2]);

      /* node 05 */
      for (l=0;l<3;l++)
        x[l] = 0.25*(actele->node[10]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
            + actele->node[14]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
            + actele->node[18]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
            + actele->node[15]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l] );
      fprintf(out," %6d %18.5e %18.5e %18.5e\n",
          nodebase + 5,
          x[0],
          x[1],
          x[2]);

      /* node 06 */
      for (l=0;l<3;l++)
        x[l] = 0.25*(actele->node[11]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
            + actele->node[15]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
            + actele->node[19]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
            + actele->node[12]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l] );
      fprintf(out," %6d %18.5e %18.5e %18.5e\n",
          nodebase + 6,
          x[0],
          x[1],
          x[2]);

      /* node 07 */
      for (l=0;l<3;l++)
        x[l] = 0.0833333333
          *(actele->node[ 8]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
              + actele->node[ 9]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
              + actele->node[10]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
              + actele->node[11]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
              + actele->node[12]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
              + actele->node[13]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
              + actele->node[14]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
              + actele->node[15]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
              + actele->node[16]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
              + actele->node[17]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
              + actele->node[18]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l]
              + actele->node[19]->gnode->mfcpnode[genprob.numaf]->sol.a.da[place][l] );
      fprintf(out," %6d %18.5e %18.5e %18.5e\n",
          nodebase + 7,
          x[0],
          x[1],
          x[2]);

    }  /* for (j=0; j<actfield->dis[0].numele; j++) */
#endif






  }


  if (ioflags.struct_disp == 1 || ioflags.ale_disp == 1)
  {
    fprintf(out,"END VALUES\n");
  }




  /* write fsi loads of structure field */
  /*------------------------------------*/
  /*   THIS DOES NOT WORK UP TO NOW!!!  */
#if 0
  if (ioflags.struct_disp==1||ioflags.fluid_sol==1)
  {
    resulttype        = "VECTOR";
    resultplace       = "ONNODES";
    rangetable        = "standard_structure";
    ncomponent        = genprob.ndim;
    componentnames[0] = "x-load";
    componentnames[1] = "y-load";
    componentnames[2] = "z-load";
    /*-------------------------------------------------------------------*/
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"# RESULT FSI LOADS on FIELD STRUCTURE, DIS %1i\n",disnums);
    fprintf(out,"# TIME %18.5E \n",fsidyn->time);
    fprintf(out,"# STEP %6d    \n",fsidyn->step);
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"RESULT %cfsiload%c %cccarat%c %d %s %s\n",
        sign,sign,
        sign,sign,
        fsidyn->step,
        resulttype,
        resultplace
        );
    fprintf(out,"RESULTRANGESTABLE %c%s%c\n",
        sign,rangetable,sign
        );
    /*-------------------------------------------------------------------*/
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
    /*-------------------------------------------------------------------*/
    fprintf(out,"VALUES\n");
  }
  if (ioflags.struct_disp==1 && structfield!=NULL)
  {
    for (i=0; i<structfield->dis[0].numnp; i++)
    {
      actnode = &(structfield->dis[0].node[i]);
      switch (genprob.ndim)
      {
        case 3:
          fprintf(out," %6d %18.5E %18.5E %18.5E\n",
              actnode->Id+1,
              fsidyn->fsiload[actnode->dof[0]],
              fsidyn->fsiload[actnode->dof[1]],
              fsidyn->fsiload[actnode->dof[2]]
              );
          break;
        case 2:
          fprintf(out," %6d %18.5E %18.5E \n",
              actnode->Id+1,
              fsidyn->fsiload[actnode->dof[0]],
              fsidyn->fsiload[actnode->dof[1]]
              );
          break;
        default:
          dserror("Unknown number of dimensions");
          break;
      }
    }
  }
  if (ioflags.struct_disp==1||ioflags.fluid_sol==1)
  {
    fprintf(out,"END VALUES\n");
  }
#endif



  /* write velocities and pressure of fluid field */
  /*----------------------------------------------*/
  if (ioflags.fluid_sol==1)
  {
    out_gid_sol("velocity",fluidfield,disnumf,actintraf,fsidyn->step,fsidyn->actpos,fsidyn->time);
    switch (genprob.probtyp)
    {
    case prb_fluid:
    case prb_fsi:
      out_gid_sol("pressure",fluidfield,disnumf,actintraf,fsidyn->step,fsidyn->actpos,fsidyn->time);
      break;
#ifdef D_FLUID_PM
    case prb_pfsi:
      switch (fdyn->dyntyp)
      {
      case dyntyp_pm_cont:
      case dyntyp_pm_cont_laplace:
        out_gid_sol("projected_pressure",fluidfield,disnumf+1,actintraf,fsidyn->step,fsidyn->actpos,fsidyn->time);
        break;
      case dyntyp_pm_discont:
	out_gid_sol("average_pressure",fluidfield,disnumf,actintraf,fsidyn->step,fsidyn->actpos,fsidyn->time);
        break;
      default:
        dserror("unsupported dynamic type %d for fluid projection", fdyn->dyntyp);
      }
      break;
#endif
    default:
      dserror("unsupported problem type %d for fluid in fsi", genprob.probtyp);
    }
  }




  /* write velocities, accelerations and stresses of structure field */
  /*-----------------------------------------------------------------*/
  if (ioflags.struct_disp==1 && structfield!=NULL)
  {
    out_gid_sol("velocities",structfield,disnums,actintras,fsidyn->step,1,fsidyn->time);
    out_gid_sol("accelerations",structfield,disnums,actintras,fsidyn->step,2,fsidyn->time);
  }
  if (ioflags.struct_stress==1 && structfield!=NULL)
    out_gid_sol("stress"      ,structfield,disnums,actintras,fsidyn->step,0,fsidyn->time);



  fflush(out);



#ifdef DEBUG
  dstrc_exit();
#endif


  return;


#endif /* NO_TEXT_OUTPUT */


} /* end of out_gid_sol_fsi */

#endif

/*! @} (documentation module close)*/
#endif

/*!----------------------------------------------------------------------
  \file
  \brief controls output of fsi for gid

  ------------------------------------------------------------------------*/
/*! 
  \addtogroup FSI
  *//*! @{ (documentation module open)*/
#ifdef D_FSI
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
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
void out_gid_sol_fsi(FIELD *fluidfield, FIELD *structfield)
{
  INT           i;

  static FLUID_DYNAMIC  *fdyn;
  static FSI_DYNAMIC    *fsidyn;
  static STRUCT_DYNAMIC *sdyn;

  INTRA        *actintras;
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
  /*----------------------------------------------------------------------*/
  if (structfield!=NULL)
  sdyn= alldyn[genprob.numsf].sdyn;
  fdyn= alldyn[genprob.numff].fdyn;
  fsidyn= alldyn[genprob.numfld].fsidyn;
  /*----------------------------------------------------------------------*/
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
  fprintf(out,"# Converged Solution of timestep %d\n",fdyn->step);
  fprintf(out,"#-------------------------------------------------------------------------------\n");

  /* ------------- write displacements of structure field */
  /*----------------------------------------------------------------------*/
  if (ioflags.struct_disp_gid==1||ioflags.fluid_sol_gid==1)
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
    fprintf(out,"# RESULT DISPLACEMENTS on FIELD FSI\n");
    fprintf(out,"# TIME %18.5E \n",fsidyn->time);   
    fprintf(out,"#-------------------------------------------------------------------------------\n");
    fprintf(out,"RESULT %cdisplacement%c %cccarat%c %d %s %s\n",
        sign,sign,
        sign,sign,
        fdyn->step,
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
  if (ioflags.struct_disp_gid==1 && structfield!=NULL)
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


    /* write displacements of fluid field */
    /*-------------------------------------------------------------------*/
  }
  if (ioflags.fluid_sol_gid==1)
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
                actanode->sol.a.da[fsidyn->actpos][0],
                actanode->sol.a.da[fsidyn->actpos][1],
                actanode->sol.a.da[fsidyn->actpos][2]
                );
            break;
          case 2:
            fprintf(out," %6d %18.5E %18.5E \n",
                actnode->Id+1,
                actanode->sol.a.da[fsidyn->actpos][0],
                actanode->sol.a.da[fsidyn->actpos][1]
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
  }
  if (ioflags.struct_disp_gid==1||ioflags.fluid_sol_gid==1)
  {
    fprintf(out,"END VALUES\n");
  }
  /*-------------------------------------------------------------------*/



  /* write velocities and pressure of fluid field */
  if (ioflags.fluid_sol_gid==1) 
  {
    out_gid_sol("velocity",fluidfield,actintraf,fdyn->step,fsidyn->actpos,fsidyn->time);
    out_gid_sol("pressure",fluidfield,actintraf,fdyn->step,fsidyn->actpos,fsidyn->time);
  }


  /* write velocities, accelerations and stresses of structure field */
  if (ioflags.struct_disp_gid==1 && structfield!=NULL)
  {
    out_gid_sol("velocities",structfield,actintras,sdyn->step,1,fsidyn->time);
    out_gid_sol("accelerations",structfield,actintras,sdyn->step,2,fsidyn->time);
  }
  if (ioflags.struct_stress_gid==1 && structfield!=NULL)
    out_gid_sol("stress"      ,structfield,actintras,sdyn->step,0,fsidyn->time);


  /*----------------------------------------------------------------------*/
  fflush(out);
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
} /* end of out_gid_sol_fsi */

#endif

/*! @} (documentation module close)*/

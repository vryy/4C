/*!----------------------------------------------------------------------
\file
\brief controls output of ssi for gid

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*!
  \addtogroup SSI
  *//*! @{ (documentation module open)*/
#ifdef D_SSI
#include "../headers/standardtypes.h"
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
\brief routine to control output og ssi to gid

<pre>                                               mfirl / genk 02/04

This function writes the output of the current timestep into the
output-file for gid.  The output includes:
- slave displacements (mesh deformation)
- slave stresses
- structure displacements
- structure stresses
- structure velocities    (not yet implemented)
- structure accelerations (not yet implemented)
</pre>

\param *slavefield    FIELD   (i)     the slave field to output
\param *masterfield   FIELD   (i)     the master field to output

\return void

------------------------------------------------------------------------*/
void out_gid_sol_ssi(FIELD *slavefield, FIELD *masterfield)
{
#ifndef NO_TEXT_OUTPUT
INT           i;

static STRUCT_DYNAMIC *sdyn_s;
static SSI_DYNAMIC    *ssidyn;
static STRUCT_DYNAMIC *sdyn_m;

INTRA        *actintram;
INTRA        *actintras;

FILE         *out = allfiles.gidres;

NODE         *actnode;

char         *resulttype;
char         *resultplace;
char         *rangetable;
INT           ncomponent;
char         *componentnames[18];

char          sign='"';

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("out_gid_sol_ssi");
#endif

/*----------------------------------------------------------------------*/
if (masterfield!=NULL)
sdyn_m= alldyn[1].sdyn;
sdyn_s= alldyn[0].sdyn;
ssidyn= alldyn[genprob.numfld].ssidyn;

/*----------------------------------------------------------------------*/
#ifdef PARALLEL
actintras    = &(par.intra[genprob.numff]);
if (masterfield!=NULL)
actintram    = &(par.intra[genprob.numsf]);
#else
actintras   = (INTRA*)CCACALLOC(1,sizeof(INTRA));
if (!actintras) dserror("Allocation of INTRA failed");
actintras->intra_fieldtyp = structure;
actintras->intra_rank     = 0;
actintras->intra_nprocs   = 1;
if (masterfield!=NULL)
{
actintram   = (INTRA*)CCACALLOC(1,sizeof(INTRA));
if (!actintram) dserror("Allocation of INTRA failed");
actintram->intra_fieldtyp = structure;
actintram->intra_rank     = 0;
actintram->intra_nprocs   = 1;
}
else
actintram = NULL; /* neither parallel nor master field! */
#endif

if (masterfield!=NULL)
if ( sdyn_m->step != sdyn_s->step)
  dserror("Something is wrong!!");

fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"# Converged Solution of timestep %d\n",ssidyn->step);
fprintf(out,"#-------------------------------------------------------------------------------\n");

/* ------------- write displacements of master field */
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
  fprintf(out,"# RESULT DISPLACEMENTS on FIELD SSI\n");
  fprintf(out,"# TIME %20.10f \n",ssidyn->time);
  fprintf(out,"#-------------------------------------------------------------------------------\n");
  fprintf(out,"RESULT %cdisplacement%c %cccarat%c %d %s %s\n",
      sign,sign,
      sign,sign,
      ssidyn->step,
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
if (ioflags.struct_disp_gid==1 && masterfield!=NULL)
{
  for (i=0; i<masterfield->dis[0].numnp; i++)
  {
    actnode = &(masterfield->dis[0].node[i]);
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


  /* write displacements of slave field */
  /*-------------------------------------------------------------------*/
}
if (ioflags.struct_disp_gid==1)
{
  for (i=0; i<slavefield->dis[0].numnp; i++)
  {
    actnode = &(slavefield->dis[0].node[i]);
    if (actnode != NULL)
    {
      switch (genprob.ndim)
      {
        case 3:
          fprintf(out," %6d %18.5E %18.5E %18.5E\n",
              actnode->Id+1,
              actnode->sol.a.da[ssidyn->actpos][0],
              actnode->sol.a.da[ssidyn->actpos][1],
              actnode->sol.a.da[ssidyn->actpos][2]
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
}
if (ioflags.struct_disp_gid==1||ioflags.fluid_sol_gid==1)
{
  fprintf(out,"END VALUES\n");
}
/*---------------------------------------------------- some output ---*/
if (ioflags.struct_stress_gid==1 && masterfield!=NULL)
  out_gid_sol("stress"      ,masterfield,actintram,ssidyn->step,0,ssidyn->time);
if (ioflags.struct_stress_gid==1 && slavefield!=NULL)
  out_gid_sol("stress"      ,slavefield,actintram,ssidyn->step,0,ssidyn->time);

fflush(out);

/*--------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
#endif /* NO_TEXT_OUTPUT */
} /* end of out_gid_sol_ssi */

#endif

/*! @} (documentation module close)*/

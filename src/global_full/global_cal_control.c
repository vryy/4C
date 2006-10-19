/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../fluid_full/fluid_prototypes.h"

#ifdef D_WALLGE
#include "../wallge/wallge.h"
#include "../wallge/wallge_prototypes.h"
#endif

#ifdef D_FSI
#include "../fsi_full/fsi_prototypes.h"
#endif

#ifdef D_SSI
#include "../ssi_full/ssi_prototypes.h"
#endif

#ifdef D_TSI
#include "../tsi_full/tsi_prototypes.h"
#endif

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB      genprob;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

#ifdef D_MLSTRUCT
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld submeshFIELDs, defined in global_control.c          |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *sm_field;
#endif /* D_MLSTRUCT */

/*----------------------------------------------------------------------*
 |  routine to control execution phase                   m.gee 6/01     |
 *----------------------------------------------------------------------*/
void ntacal()
{
INT i;
FIELD *actfield;
#ifdef D_MLSTRUCT
FIELD *actsmfield;
#endif /* D_MLSTRUCT */

#ifdef DEBUG
dstrc_enter("ntacal");
#endif
/*----------------------------------------------------------------------*/
/*------------------------do initial partitioning of nodes and elements */
#ifdef PERF
  perf_begin(12);
#endif
  part_fields();
#ifdef PERF
  perf_end(12);
#endif

#ifdef D_FLUID
/*---------------------------------- set dofs for implicit free surface */
if (genprob.numff>=0) fluid_freesurf_setdofs();
/*----------------------------- modify coordinates for special problems */
if (genprob.numff>=0) fluid_modcoor();
#endif
/*-------------------------------- set dofs for gradient enhanced model */
#ifdef D_WALLGE
if (genprob.graderw>0) wge_setdof();
#endif

/*------------------------------------------------ assign dofs to nodes */
#ifdef PERF
  perf_begin(13);
#endif
  for(i=0; i<genprob.numfld; i++)
  {
    actfield = &(field[i]);
    if (actfield->ndis==1) assign_dof(actfield);
    if (actfield->ndis>1) assign_dof_ndis(actfield);
  }
#ifdef PERF
  perf_end(13);
#endif
/*--------------------------------- assign dofs to nodes in the submesh */
#ifdef D_MLSTRUCT
if (genprob.multisc_struct == 1)
{
   actsmfield = &(sm_field[0]);
   assign_dof(actsmfield);
   if (actsmfield->dis[0].numeq >= 3000)
      dserror(">3000 submesh-DOF's->size of fields ele.e.w1.fint_mi,..for static cond.have to be enlarged\n");
}
#endif /* D_MLSTRUCT */
/*--------make the procs know their own nodes and elements a bit better */
#ifdef PERF
  perf_begin(14);
#endif
part_assignfield();
#ifdef PERF
  perf_end(14);
#endif


  /* divide fast elements into vectors */
#ifdef D_FLUID3_F
  divide_fast();
#endif


  /* check the values of the defines for MAXNODE etc. */
#ifdef CHECK_MAX
  if (par.myrank==0)
    check_max_sizes();
#endif


/*-------------------calculate system matrices parallel storage formats */
#ifdef PERF
  perf_begin(15);
#endif

  mask_global_matrices();

#ifdef PERF
  perf_end(15);
#endif
#ifdef D_MLSTRUCT
if (genprob.multisc_struct == 1)
{
  mask_submesh_matrices();
}
#endif /* D_MLSTRUCT */

/*------------------- inherit local co-ordinate systems to the elements */
locsys_inherit_to_node();

/*------------------------------------------------ write general output */
out_general();
/*--------------------------------------------------- write mesh to gid */
if (par.myrank==0 && ioflags.output_gid)
{
  out_gid_sol_init();
  if (genprob.probtyp != prb_structure && genprob.probtyp != prb_fsi)
    out_gid_msh();
}
/*------------------------ program to control execution of optimization */
/*------------------ call control programs of static or dynamic control */

switch (genprob.probtyp) {
case prb_structure:

  switch (genprob.timetyp) {
  case time_static:
    calsta();
    break;
  case time_dynamic:
    caldyn();
    break;
  default:
    dserror("Unspecified time handling");
  }

  break;

#ifdef D_FLUID
case prb_fluid:
case prb_fluid_pm:
  dyn_fluid();
  break;
#endif

#ifdef D_FSI
case prb_fsi:
  dyn_fsi(0);
  break;
#endif

#ifdef D_SSI
case prb_ssi:
  dyn_ssi();
  break;
#endif

#ifdef D_ALE
case prb_ale:
  dyn_ale();
  break;
#endif

#ifdef D_OPTIM
case prb_opt:
  caloptmain();
  break;
#endif

#ifdef D_TSI
case prb_tsi:
  tsi_dyn();
  break;
#endif

default:
  dserror("solution of unknown problemtyp requested");
break;
}

/*------------------------------------------------------- check results */
#ifdef RESULTTEST
global_result_test();
#endif

/*--------------------------------------------------- write warnings ---*/
dswarning(2,0);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of ntacal */

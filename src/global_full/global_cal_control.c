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

#ifdef D_LS
#include "../ls/ls_prototypes.h"
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
extern struct _GENPROB     genprob;
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
/*----------------------------------------------------------------------*
 |  routine to control execution phase                   m.gee 6/01     |
 *----------------------------------------------------------------------*/
void ntacal()
{
INT i;
FIELD *actfield;

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
/*--------make the procs know their own nodes and elements a bit better */
#ifdef PERF
  perf_begin(14);
#endif
part_assignfield();
#ifdef PERF
  perf_end(14);
#endif


  /* check the values of the defines for MAXNODE etc. */
#ifdef CHECK_MAX
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

/*------------------- inherit local co-ordinate systems to the elements */
locsys_inherit_to_node();

/*------------------------------------------------ write general output */
out_general();
/*--------------------------------------------------- write mesh to gid */
if (par.myrank==0)
{
   if (ioflags.struct_disp_gid||ioflags.struct_stress_gid||ioflags.fluid_sol_gid
       ||ioflags.ale_disp_gid||ioflags.fluid_stress_gid)
   {
      out_gid_sol_init();
      if (genprob.probtyp != prb_structure && genprob.probtyp != prb_fsi)
         out_gid_msh();
   }
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
  dyn_fluid();
  break;
#endif

#ifdef D_FSI
case prb_fsi:
  dyn_fsi(0);
  break;
#endif

#ifdef D_ALE
case prb_ale:
  dyn_ale();
  break;
#endif

#ifdef D_LS
case prb_twophase:
case prb_levelset:
  ls_dyn();
  break;
#endif

#ifdef D_OPTIM
case prb_opt:
  caloptmain();
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

/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
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

/* header for DRT style input */
#ifdef CCADISCRET
#include "../drt_lib/global_inp_control2.H"
#include "../drt_fluid/fluid_dyn_nln_drt.H"
#include "../drt_fluid/xfluid_dyn_nln_drt.H"
#include "../drt_fluid/condif_drt.H"
#include "../drt_fsi/ale_dyn.H"
#include "../drt_fsi/fsi_dyn.H"
#include "../drt_elch/elch_dyn.H"
#endif

/*----------------------------------------------------------------------*
 |  routine to control execution phase                   m.gee 6/01     |
 *----------------------------------------------------------------------*/
void ntacal()
{
#ifndef CCADISCRET
INT i;
FIELD *actfield;
#endif
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

#ifndef CCADISCRET /* the 'old' ccarat style discretization management */
  part_fields();
#endif
#ifdef PERF
  perf_end(12);
#endif

#ifdef D_FLUID
#ifndef CCADISCRET
/*---------------------------------- set dofs for implicit free surface */
if (genprob.numff>=0) fluid_freesurf_setdofs();
/*----------------------------- modify coordinates for special problems */
if (genprob.numff>=0) fluid_modcoor();
#endif
#endif
/*-------------------------------- set dofs for gradient enhanced model */
#ifdef D_WALLGE
#ifndef CCADISCRET
if (genprob.graderw>0) wge_setdof();
#endif
#endif

/*------------------------------------------------ assign dofs to nodes */
#ifdef PERF
  perf_begin(13);
#endif
#ifndef CCADISCRET
  for(i=0; i<genprob.numfld; i++)
  {
    actfield = &(field[i]);
    if (actfield->ndis==1) assign_dof(actfield);
    if (actfield->ndis>1) assign_dof_ndis(actfield);
  }
#endif
#ifdef PERF
  perf_end(13);
#endif

/*--------------------------------- assign dofs to nodes in the submesh */
#ifndef CCADISCRET
#ifdef D_MLSTRUCT
if (genprob.multisc_struct == 1)
{
   actsmfield = &(sm_field[0]);
   assign_dof(actsmfield);
   if (actsmfield->dis[0].numeq >= 3000)
      dserror(">3000 submesh-DOF's->size of fields ele.e.w1.fint_mi,..for static cond.have to be enlarged\n");
}
#endif /* D_MLSTRUCT */
#endif
/*--------make the procs know their own nodes and elements a bit better */
#ifndef CCADISCRET
#ifdef PERF
  perf_begin(14);
#endif
part_assignfield();
#ifdef PERF
  perf_end(14);
#endif
#endif

  /* divide fast elements into vectors */
#ifndef CCADISCRET
#ifdef D_FLUID3_F
  divide_fast();
#endif
#endif

  /* check the values of the defines for MAXNODE etc. */
#ifndef CCADISCRET
  if (par.myrank==0)
    check_max_sizes();
#endif


/*-------------------calculate system matrices parallel storage formats */
#ifndef CCADISCRET
#ifdef PERF
  perf_begin(15);
#endif
  mask_global_matrices();
#ifdef PERF
  perf_end(15);
#endif
#endif

#ifndef CCADISCRET
#ifdef D_MLSTRUCT
if (genprob.multisc_struct == 1)
{
  mask_submesh_matrices();
}
#endif /* D_MLSTRUCT */
#endif

/*------------------- inherit local co-ordinate systems to the elements */
#ifndef CCADISCRET
locsys_inherit_to_node();
#endif

#ifndef CCADISCRET
/*------------------------------------------------ write general output */
out_general();
#endif
/*--------------------------------------------------- write mesh to gid */
#ifndef CCADISCRET
if (par.myrank==0 && ioflags.output_gid)
{
  out_gid_sol_init();
  if (genprob.probtyp != prb_structure && genprob.probtyp != prb_fsi && genprob.probtyp != prb_pfsi)
    out_gid_msh();
}
#else
#endif
/*------------------------ program to control execution of optimization */
/*------------------ call control programs of static or dynamic control */

switch (genprob.probtyp) {
case prb_structure:
case prb_struct_multi:

  switch (genprob.timetyp) {
  case time_static:
    calsta();
/*#ifndef CCADISCRET
#else
    dserror("calsta with DRT not yet impl.");
#endif*/
    break;
  case time_dynamic:
#ifndef CCADISCRET
  caldyn();
#else
  caldyn_drt();
#endif
    break;
  default:
    dserror("Unspecified time handling");
  }

  break;

#ifdef D_FLUID
case prb_fluid:
case prb_fluid_pm:
#ifndef CCADISCRET
  dyn_fluid();
#else
  dyn_fluid_drt();
#endif
  break;
case prb_condif:
#ifndef CCADISCRET
  dserror("dyn_condif without DRT not implemented");
#else
  dyn_condif_drt();
#endif
  break;
case prb_fluid_xfem:
#ifdef CCADISCRET
  xdyn_fluid_drt();
#else
  dserror("prb_fluid_xfem without DRT not implemented");
#endif
#endif
  break;

#ifdef D_FSI
case prb_fsi:
case prb_pfsi:
#ifndef CCADISCRET
  dyn_fsi(0);
#else
  fsi_ale_drt();
#endif
  break;
case prb_fsi_xfem:
#ifndef CCADISCRET
  dserror("prb_xfsi without DRT not impl.");;
#else
  xfsi_drt();
#endif
  break;
#endif

#ifdef D_SSI
case prb_ssi:
  dyn_ssi();
  break;
#endif

#ifdef D_ALE
case prb_ale:
#ifndef CCADISCRET
  dyn_ale();
#else
  dyn_ale_drt();
#endif
  break;
#endif

#ifdef D_OPTIM
case prb_opt:
  caloptmain();
  break;
#endif

#ifdef D_TSI
case prb_tsi:
#ifndef CCADISCRET
  tsi_dyn();
#else
  dserror("tsi_dyn() with DRT not yet impl.");
#endif
  break;
#endif

case prb_elch:
#ifndef CCADISCRET
  dserror("Electrochemistry module not available in CCARAT");
#else
  elch_dyn();
#endif
  break;

default:
  dserror("solution of unknown problemtyp requested");
break;
}

#ifdef CCADISCRET
drt_problem_done();
#endif

/*------------------------------------------------------- check results */
#ifndef CCADISCRET
#ifdef RESULTTEST
global_result_test();
#endif
#endif

/*--------------------------------------------------- write warnings ---*/
dswarning(2,0);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of ntacal */

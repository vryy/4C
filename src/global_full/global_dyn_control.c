/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../fluid_full/fluid_prototypes.h"
#include "../fsi_full/fsi_prototypes.h"
#ifdef CCADISCRET
#include "../drt_structure/stru_dyn_nln_drt.H"
#endif
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;
/*----------------------------------------------------------------------*
 |  routine to control dynamic execution                 m.gee 5/01     |
 *----------------------------------------------------------------------*/
void caldyn()
{
#ifdef DEBUG
  dstrc_enter("caldyn");
#endif
/*------------------------------ switch into different time integrators */
  switch (alldyn[0].sdyn->Typ)
  {
    /* central differences time integration */
    case centr_diff:
#ifndef CCADISCRET
      dyn_nln_stru_expl();
#else
      dserror("no central differences in DRT");
#endif
      break;
    /* generalized alfa time integration */
    case gen_alfa:
#ifndef CCADISCRET
      dyn_nln_structural();
#else
      /* stru_genalpha_drt(); */
      dyn_nlnstructural_drt();
#endif
      break;
    /* Generalized Energy-Momentum time integration */
    case Gen_EMM:
#ifdef GEMM
      dyn_nln_gemm();
#else
      dserror("GEMM not supported");
#endif
      break;
    default:
      dserror("Untimely unknown time integration scheme");
      break;
  }


/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
return;
} /* end of caldyn */

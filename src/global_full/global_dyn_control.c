/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0771 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../fluid_full/fluid_prototypes.h"
#include "../fsi_full/fsi_prototypes.h"
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
/*----------------------------------------------------------------------*/

if (alldyn[0].sdyn->Typ == centr_diff) {
  /* central difference time integration */
  dyn_nln_stru_expl();
}
else if (alldyn[0].sdyn->Typ == gen_alfa) {
  /* generalized alfa time integration */
  dyn_nln_structural();
}
else if (alldyn[0].sdyn->Typ == Gen_EMM) {
#ifdef GEMM
  dyn_nln_gemm(); /* Generalized Energy-Momentum time integration */
#else
  dserror("GEMM not supported");
#endif
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of caldyn */

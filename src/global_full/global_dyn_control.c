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
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
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
switch (genprob.probtyp)
{
case prb_structure:
    if (alldyn[0].sdyn->Typ == centr_diff)
    {
       dyn_nln_stru_expl();/* central difference time integration */
    }
    if (alldyn[0].sdyn->Typ == gen_alfa)
    {
       dyn_nln_structural();/* generalized alfa time integration */
    }
    if (alldyn[0].sdyn->Typ == Gen_EMM)
    {
#ifdef GEMM
       dyn_nln_gemm(); /* Generalized Energy-Momentum time integration */
#endif 
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
default:
    dserror("Dynamic solution of unknown Problemtyp requested");
break;
}



/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of caldyn */

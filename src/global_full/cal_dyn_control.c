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
       dyn_nln_stru_expl();/* generalized alfa time integration */
    }
    if (alldyn[0].sdyn->Typ == gen_alfa)
    {
       dyn_nln_structural();/* central difference time integration */
    }
break;
case prb_fluid:
    dyn_fluid();
break;
case prb_fsi:
    dyn_fsi(0);
break;    
case prb_ale:
    dyn_ale();
break;
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

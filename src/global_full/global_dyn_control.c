#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
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
    dyn_nln_structural();
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

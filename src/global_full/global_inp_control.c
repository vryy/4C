#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 | input of control, element and load information         m.gee 8/00    |
 *----------------------------------------------------------------------*/
void ntainp()
{
#ifdef DEBUG 
trace.actroutine = trace.actroutine->next;
strncpy(trace.actroutine->name,"ntainp",49);
trace.actroutine->dsroutcontrol=dsin;
trace.deepness++;
#endif

/*--------------------------------------------- input of tracing option */
#ifdef DEBUG 
inptrace();
#endif
/*----------------------- input of not mesh or time based problem data  */
inpctr();
/*----------------------------------- input of design data if necessary */
if (genprob.design==1) inpdesign();
/*------------------------------------------------------input of meshes */
inpfield();
/*---------------------------------------------------design-fe topology */
if (genprob.design==1) inpdesign_topology();
/*------------------------------------------------- input of conditions */
inp_conditions();
/*-------------------------------------------------- input of materials */
inp_material();
/*--------------------------------------------- all reading is over here*/

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ntainp */

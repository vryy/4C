#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | tracing variables                                                    |
 | defined in pss_ds.c                                                  |
 *----------------------------------------------------------------------*/
#ifdef DEBUG
extern struct _TRACE         trace;
#endif
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 | input of control, element and load information         m.gee 8/00    |
 *----------------------------------------------------------------------*/
void ntainp()
{
/* 
   the input of the tracing option has not been done yet, so
   we have to make the dstrc_enter 'by hand' 
*/
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
/*=========================== tracing is active and working from now on */
/*----------------------- input of not mesh or time based problem data  */
inpctr();
/*----------------------------------- input of design data if necessary */
if (genprob.design==1) inpdesign();
/*------------------------------------------------------input of meshes */
inpfield();
/*--------------------------------------- input of general dynamic data */
if (genprob.timetyp==time_dynamic) inpctrdyn();
/*---------------------------------------- input of general static data */
else inpctrstat();
/*-----------------------------------design-design & design-fe topology */
if (genprob.design==1) 
{
   inpdesign_topology_design();
   inpdesign_topology_fe();
}
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

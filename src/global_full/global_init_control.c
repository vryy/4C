#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 | Initialize program service systems                     m.gee 8/00    |
 *----------------------------------------------------------------------*/
void ntaini(int argc, char *argv[])
{
/*---------------------------------------------------initialize tracing */
#ifdef DEBUG 
dsinit();
#endif
/*-------------------------------------------------------initialize I/O */
ntadev(argc,argv);
/*------------------------------------------initialize free-field-input */
frinit();

#ifdef DEBUG 
trace.actroutine->dsroutcontrol=dsout;
trace.actroutine = trace.actroutine->prev;
trace.deepness--;
#endif
return;
} /* end of ntaini */

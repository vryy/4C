#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
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
int i,j;
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
/*------------------------------------------------ input of design data */
/* read description of design volumes, surfaces, lines and nodes, allocate
   and fill the structures DVOL, DSURF, DLINE and DNODE 
*/
inpdesign();
/* build the topology among the design elements, allocate and create the
   pointers connecting DVOL, DSURF, DLINE and DNODE 
*/
inpdesign_topology_design();
/*------------------------------------------------------input of meshes */
/* read the fe-nodes in NODE, the fe-elements in ELEMENT and build the 
   topology among ELEMENTs and NODEs. Assign the disctretization to
   the different fields in a multifield problem dependent on type
   of element
*/    
inpfield();
/*------------------ built the detailed topology of the discretizations */
/* for each existing field and discretization build a detailed topology 
   of GNODEs, GLINEs, GSURFs and GVOLs and connect them to the topology
   of NODEs and ELEMENTs
*/   
for (i=0; i<genprob.numfld; i++)
for (j=0; j<field[i].ndis; j++)
inp_detailed_topology(&(field[i].dis[j]));
/*---------------------------------------------------design-fe topology */
/* Read which node is on which design object from file.
   For each field and discretization build the pointers among design
   and FE-objects DVOL<->GVOL,DSURF<->GSURF,DLINE<->GLINE,DNODE<->GNODE.
*/   
inpdesign_topology_fe();
/*--------------------------------------- input of general dynamic data */
if (genprob.timetyp==time_dynamic) inpctrdyn();
/*---------------------------------------- input of general static data */
else inpctrstat();
/*------------------------------------------------- input of conditions */
/* dirichlet/coupling/neumann conditions are read from file to the
   DVOLS/DSURFS/DLINES/DNODES
*/   
inp_conditions();
/*----- inherit the dirichlet and coupling conditions inside the design */
/* conditions are inherited 'downwards': DVOL->DSURF->DLINE->DNODE      */
/* BUT: if a 'lower object already has its own condition, it does NOT   */
/* inherit from above, it does inherit its own condition further down   */
inherit_dirich_coup_indesign();
/* set pointers in the discretization to the design dirichlet conditions*/
for (i=0; i<genprob.numfld; i++)
for (j=0; j<field[i].ndis; j++)
inherit_design_dis_dirichlet(&(field[i].dis[j]));
/*--- set pointers in the discretization to the design couple conditions*/
for (i=0; i<genprob.numfld; i++)
for (j=0; j<field[i].ndis; j++)
inherit_design_dis_couple(&(field[i].dis[j]));
/*-------------------------------------------------- input of materials */
inp_material();
/*------------------------------------- input of initial data for fluid */
inp_fluid_start_data();
/*--------------------------------------------- all reading is over here*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ntainp */

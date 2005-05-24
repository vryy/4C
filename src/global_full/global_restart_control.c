/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../headers/restart.h"
#include "../shell8/shell8.h"
#include "../shell9/shell9.h"
#include "../wall1/wall1.h"
#include "../beam3/beam3.h"
#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3.h"
#include "../fluid2/fluid2_prototypes.h"
#include "../fluid3/fluid3_prototypes.h"

#ifdef D_FLUID3_F
#include "../fluid3_fast/f3f_prototypes.h"
#endif

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | struct _GENPROB       genprob; defined in global_control.c           |
 *----------------------------------------------------------------------*/
extern struct _GENPROB       genprob;
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 | write restart                                         m.gee 05/02    |
 | of nln-structural dynamics                                           |
 | How to restart:                                                      |
 | 1.) run calculation with RESTART 0 and NUMSTEP 10 in the input file  |
 | 2.) cp xxx.pss xxx_0.pss                                             |
 | 3.) run calculation with RESTART 9 and NUMSTEP 20 in the input file  |
 |     by starting with                                                 |
 |     cca_seq_10.20.exe input.dat xxx restart                          |
 | 4.) cp xxx.respss xxx_1.pss                                          |
 | 5.) cp xxx.respss xxx.pss                                            |
 | 6.) run calculation with RESTART 19 and NUMSTEP 30 in the input file |
 |     by starting with                                                 |
 |     cca_seq_10.20.exe input.dat xxx restart                          |
 | 7.) cp xxx.respss xxx_2.pss                                          |
 | 8.) cp xxx.respss xxx.pss                                            |
 | generally: restart is read from .pss and writes to .respss           |
 |            a restarted calculation runs from step                    |
 |            RESTART i+1 to NUMSTEP j-1                                |
 *----------------------------------------------------------------------*/
void restart_write_nlnstructdyn(STRUCT_DYNAMIC  *sdyn,
                                STRUCT_DYN_CALC *dynvar,
                                FIELD           *actfield,
                                PARTITION       *actpart,
                                INTRA           *actintra,
                                CALC_ACTION     *action,
                                INT nrhs,  DIST_VECTOR *rhs,
                                INT nsol,  DIST_VECTOR *sol,
                                INT ndis,  DIST_VECTOR *dispi,
                                INT nvel,  DIST_VECTOR *vel,
                                INT nacc,  DIST_VECTOR *acc,
                                INT nfie,  DIST_VECTOR *fie,
                                INT nwork, DIST_VECTOR *work,
                                ARRAY *intforce_a,
                                ARRAY *dirich_a,
                                CONTAINER    *container)     /*!< contains variables defined in container.h */
{
#ifndef NO_TEXT_OUTPUT
INT                  i;
INT                  ierr;
INT                  numnp;
long int           **node_handles;
INT                  numele;
long int           **ele_handles;
char                 resname[100];
long int             longdummy;
FILE                *out;
RESTART_DYNSTRUCT    res;
DIST_VECTOR         *distwrite;
NODE                *actnode;
ELEMENT             *actele;
#ifdef DEBUG
dstrc_enter("restart_write_nlnstructdyn");
#endif
#ifdef PERF
  perf_begin(74);
#endif
/*----------------------------------------------------------------------*/
out = allfiles.out_pss;
/*-------- check the step we are in and create the name "res<step>" */
res.step = sdyn->step;
sprintf(resname,"sres%d",res.step);
/*------------------- write the structure sdyn to the restart structure */
res.sdyn = *sdyn;
/*----------------- write the structure dynvar to the restart structure */
res.dynvar = *dynvar;
/* write all distributed vectors and store the handles associated with them */
res.dist_vec_rhs[0]   = nrhs;
res.dist_vec_sol[0]   = nsol;
res.dist_vec_dispi[0] = ndis;
res.dist_vec_vel[0]   = nvel;
res.dist_vec_acc[0]   = nacc;
res.dist_vec_fie[0]   = nfie;
res.dist_vec_work[0]  = nwork;
/*------------------------------------------------------ write the *rhs */
 for (i=0; i<nrhs; i++)
 {
    distwrite = &(rhs[i]);
    pss_write_array(&(distwrite->vec),&(res.dist_vec_rhs[i+1]),out,&ierr);
    if (ierr != 1) dserror("Error writing restart data");
 }
/*------------------------------------------------------ write the *sol */
 for (i=0; i<nsol; i++)
 {
    distwrite = &(sol[i]);
    pss_write_array(&(distwrite->vec),&(res.dist_vec_sol[i+1]),out,&ierr);
    if (ierr != 1) dserror("Error writing restart data");
 }
/*------------------------------------------------------ write the *dispi */
 for (i=0; i<ndis; i++)
 {
    distwrite = &(dispi[i]);
    pss_write_array(&(distwrite->vec),&(res.dist_vec_dispi[i+1]),out,&ierr);
    if (ierr != 1) dserror("Error writing restart data");
 }
/*------------------------------------------------------ write the *vel */
 for (i=0; i<nvel; i++)
 {
    distwrite = &(vel[i]);
    pss_write_array(&(distwrite->vec),&(res.dist_vec_vel[i+1]),out,&ierr);
    if (ierr != 1) dserror("Error writing restart data");
 }
/*------------------------------------------------------ write the *acc */
 for (i=0; i<nacc; i++)
 {
    distwrite = &(acc[i]);
    pss_write_array(&(distwrite->vec),&(res.dist_vec_acc[i+1]),out,&ierr);
    if (ierr != 1) dserror("Error writing restart data");
 }
/*------------------------------------------------------ write the *fie */
 for (i=0; i<nfie; i++)
 {
    distwrite = &(fie[i]);
    pss_write_array(&(distwrite->vec),&(res.dist_vec_fie[i+1]),out,&ierr);
    if (ierr != 1) dserror("Error writing restart data");
 }
/*------------------------------------------------------ write the *work */
 for (i=0; i<nwork; i++)
 {
    distwrite = &(work[i]);
    pss_write_array(&(distwrite->vec),&(res.dist_vec_work[i+1]),out,&ierr);
    if (ierr != 1) dserror("Error writing restart data");
 }
/*----------------------------------------- write the ARRAY *intforce_a */
pss_write_array(intforce_a,&(res.intforce),out,&ierr);
/*------------------------------------------- write the ARRAY *dirich_a */
pss_write_array(dirich_a,&(res.dirich),out,&ierr);
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   now write the data, that is stored in the nodes of this processor's
   partition

   a structure node as the following ARRAYS:
   actnode->sol
   actnode->sol_increment
   actnode->sol_residual
   in case of a multifield problem there exists
   actnode->sol_mf

   so we need to store an array of numnp x 4 handles
*/
/*----------------------------------------------------------------------*/
numnp = actpart->pdis[0].numnp;
res.node_handles = (long int**)CCAMALLOC(numnp*sizeof(long int*));
node_handles = res.node_handles;
res.node_handles[0] = (long int*)CCAMALLOC(4*numnp*sizeof(long int));
for (i=1; i<numnp; i++)
node_handles[i] = &(node_handles[0][i*4]);
/*----------------------------------------------------------------------*/
/* now we loop the nodes on the partition and each node writes his ARRAYs */
for (i=0; i<numnp; i++)
{
   actnode = actpart->pdis[0].node[i];
   pss_write_array(&(actnode->sol),&(node_handles[i][0]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   pss_write_array(&(actnode->sol_increment),&(node_handles[i][1]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   pss_write_array(&(actnode->sol_residual),&(node_handles[i][2]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   if (actnode->sol_mf.Typ == cca_DA)
   {
   pss_write_array(&(actnode->sol_mf),&(node_handles[i][3]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   }

}
/*
   all nodal data is written, so we now write the node_handles and store
   the handle to this in res.handle_of_node_handles
*/
pss_write("nod_hand",numnp,4,sizeof(long int),node_handles[0],&(res.handle_of_node_handles),out,&ierr);
if (ierr != 1) dserror("Error writing restart data");
/*----------------- delete the res.node_handles but keep the dimensions */
res.node_fdim = numnp;
res.node_sdim = 4;
CCAFREE(res.node_handles[0]);
CCAFREE(res.node_handles);
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   we come to the most difficult part, the writing of the element data.
   One element is assumed not to have more then MAXRECORDPERELE different
   records to write of arbitary type

   The element is called with the element action action = write_restart
   Again, all elements to the processors partition are written
*/
/*----------------------------------------------------------------------*/
numele = actpart->pdis[0].numele;
res.ele_handles = (long int**)CCAMALLOC(numele*sizeof(long int*));
ele_handles = res.ele_handles;
res.ele_handles[0] = (long int*)CCAMALLOC(MAXRECORDPERELE*numele*sizeof(long int));
for (i=1; i<numele; i++)
ele_handles[i] = &(ele_handles[0][i*MAXRECORDPERELE]);
/*--------------------- now loop element and switch for type of element */
*action = write_restart;
for (i=0; i<actpart->pdis[0].numele; i++)
{
   actele = actpart->pdis[0].element[i];
   switch(actele->eltyp)/*======================= call element routines */
   {
#ifdef D_SHELL8
   case el_shell8:
      container->kstep    = 0;
      container->handsize = MAXRECORDPERELE;
      container->handles  = ele_handles[i];
      shell8(actfield,actpart,actintra,actele,
             NULL,NULL,NULL,action,container);
   break;
#endif
   case el_shell9:
       dserror("Restart 'nlnstructdyn' for shell9 not yet impl.");
   break;
   case el_brick1:

   break;
#ifdef D_WALL1
   case el_wall1:
      container->handsize = MAXRECORDPERELE;
      container->handles  = ele_handles[i];
      wall1(actpart,actintra,actele,NULL,NULL,NULL,
            action,container);
   break;
#endif
   case el_fluid2:
       dserror("Restart for fluid2 not yet impl.");
   break;
   case el_fluid3:
       dserror("Restart for fluid3 not yet impl.");
   break;
   case el_ale3:
       dserror("Restart for ale not yet impl.");
   break;
#ifdef D_BEAM3
   case el_beam3:
      container->handsize = MAXRECORDPERELE;
      container->handles  = ele_handles[i];
      beam3(actfield,actpart,actintra,actele,
            NULL,NULL,NULL,action,container);
   break;
#endif
   case el_interf:
   break;
   case el_wallge:
   break;
   case el_none:
      dserror("Typ of element unknown");
   break;
   default:
      dserror("Typ of element unknown");
   }/* end of calling elements */
}
/*----------------------------------------------------------------------*/
/*
   all ele data is written, so write the ele_handles and store the handle to
   it in res.handle_of_ele_handles
*/
res.ele_fdim = numele;
res.ele_sdim = MAXRECORDPERELE;
pss_write("ele_hand",numele,MAXRECORDPERELE,sizeof(long int),ele_handles[0],&(res.handle_of_ele_handles),out,&ierr);
if (ierr != 1) dserror("Error writing restart data");
/*------------------ delete the res.ele_handles but keep the dimensions */
CCAFREE(res.ele_handles[0]);
CCAFREE(res.ele_handles);
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   the only thing to do is now write the RESTART_DYNSTRUCT itself with
   its unique name resname = "res<step>"
   NOTE:
   names are limited to 9 characters, so a step larger then sres99999
   can not be restarted at the moment !!!!
*/
pss_write(resname,1,1,sizeof(RESTART_DYNSTRUCT),&res,&longdummy,out,&ierr);
if (ierr != 1) dserror("Error writing restart data");
/*----------------------------------------------------------------------*/
/*
#ifdef DEBUG
pss_status_to_err();
#endif
*/
/*----------------------------------------------------------------------*/
/*
   now write a notice to the err file, that this retsrat step was
   successfully written
*/
if (res.step != 0)
{
   fprintf(allfiles.out_err,"===========================================\n");
   fprintf(allfiles.out_err,"RESTART MESSAGE\n");
   fprintf(allfiles.out_err,"In step %d restart data was written\n",res.step);
   fprintf(allfiles.out_err,"Calculation can be restarted using \n");
   fprintf(allfiles.out_err,"RESTART      %d \n",res.step);
   fprintf(allfiles.out_err,"in the input file\n");
   fprintf(allfiles.out_err,"===========================================\n");
   fflush(allfiles.out_err);
}
/*----------------------------------------------------------------------*/
#ifdef PERF
  perf_end(74);
#endif
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of restart_write_nlnstructdyn */




/*----------------------------------------------------------------------*
 | read restart                                          m.gee 05/02    |
 | of nln-structural dynamics                                           |
 *----------------------------------------------------------------------*/
void restart_read_nlnstructdyn(INT restart,
                               STRUCT_DYNAMIC  *sdyn,
                               STRUCT_DYN_CALC *dynvar,
                               FIELD           *actfield,
                               PARTITION       *actpart,
                               INTRA           *actintra,
                               CALC_ACTION     *action,
                               INT nrhs,  DIST_VECTOR *rhs,
                               INT nsol,  DIST_VECTOR *sol,
                               INT ndis,  DIST_VECTOR *dispi,
                               INT nvel,  DIST_VECTOR *vel,
                               INT nacc,  DIST_VECTOR *acc,
                               INT nfie,  DIST_VECTOR *fie,
                               INT nwork, DIST_VECTOR *work,
                               ARRAY *intforce_a,
                               ARRAY *dirich_a,
                               CONTAINER    *container)     /*!< contains variables defined in container.h */
{
#ifndef NO_TEXT_OUTPUT
INT                  i;
INT                  ierr;
long int             reshandle;
INT                  byte;
INT                  dims[3];
INT                  numnp;
long int           **node_handles;
INT                  numele;
long int           **ele_handles;
char                 resname[100];
FILE                *in;
RESTART_DYNSTRUCT    res;
DIST_VECTOR         *distread;
NODE                *actnode;
ELEMENT             *actele;

#ifdef PARALLEL
INT                  fdim;
INT                  j;
INT                  sender;
#endif

#ifdef DEBUG
dstrc_enter("restart_read_nlnstructdyn");
#endif
/*----------------------------------------------------------------------*/
in = allfiles.in_pss;
/*----------------------------------- check the step that shall be read */
sprintf(resname,"sres%d",restart);
pss_chck(resname,&reshandle,in,&ierr);
if (ierr != 1) dserror("Cannot restart, step doesn't exist in pss-file");
/*----------------------------- the structure res exists, so we read it */
pss_read_name_handle(resname,&i,&i,&byte,&res,&reshandle,in,&ierr);
if (ierr != 1) dserror("Restart structure exists, but cannot read it");
/*------------------------------- write the structure res.sdyn to *sdyn */
*sdyn   = res.sdyn;
/*--------------------------- write the structure res.dynvar to *dynvar */
*dynvar = res.dynvar;
/*---------------------------------------------------- make some checks */
if (res.dist_vec_rhs[0]   != nrhs) dserror("Mismatch in dimensions in restart");
if (res.dist_vec_sol[0]   != nsol) dserror("Mismatch in dimensions in restart");
if (res.dist_vec_dispi[0] != ndis) dserror("Mismatch in dimensions in restart");
if (res.dist_vec_vel[0]   != nvel) dserror("Mismatch in dimensions in restart");
if (res.dist_vec_acc[0]   != nacc) dserror("Mismatch in dimensions in restart");
if (res.dist_vec_fie[0]   != nfie) dserror("Mismatch in dimensions in restart");
if (res.dist_vec_work[0]  != nwork) dserror("Mismatch in dimensions in restart");
/*---------------------------------------- read the distributed vectors */
/*------------------------------------------------------------ read rhs */
for (i=0; i<nrhs; i++)
{
   distread = &(rhs[i]);
   pss_read_array_name_handle(distread->vec.name,&(distread->vec),&(res.dist_vec_rhs[i+1]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
}
/*------------------------------------------------------------ read sol */
for (i=0; i<nsol; i++)
{
   distread = &(sol[i]);
   pss_read_array_name_handle(distread->vec.name,&(distread->vec),&(res.dist_vec_sol[i+1]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
}
/*------------------------------------------------------------ read dispi */
for (i=0; i<ndis; i++)
{
   distread = &(dispi[i]);
   pss_read_array_name_handle(distread->vec.name,&(distread->vec),&(res.dist_vec_dispi[i+1]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
}
/*------------------------------------------------------------ read vel */
for (i=0; i<nvel; i++)
{
   distread = &(vel[i]);
   pss_read_array_name_handle(distread->vec.name,&(distread->vec),&(res.dist_vec_vel[i+1]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
}
/*------------------------------------------------------------ read acc */
for (i=0; i<nacc; i++)
{
   distread = &(acc[i]);
   pss_read_array_name_handle(distread->vec.name,&(distread->vec),&(res.dist_vec_acc[i+1]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
}
/*------------------------------------------------------------ read fie */
for (i=0; i<nfie; i++)
{
   distread = &(fie[i]);
   pss_read_array_name_handle(distread->vec.name,&(distread->vec),&(res.dist_vec_fie[i+1]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
}
/*------------------------------------------------------------ read work */
for (i=0; i<nwork; i++)
{
   distread = &(work[i]);
   pss_read_array_name_handle(distread->vec.name,&(distread->vec),&(res.dist_vec_work[i+1]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
}
/*------------------------------------------- read the ARRAY *intforce_a */
pss_read_array_name_handle(intforce_a->name,intforce_a,&(res.intforce),in,&ierr);
if (ierr != 1) dserror("Cannot read restart data");
/*--------------------------------------------- read the ARRAY *dirich_a */
pss_read_array_name_handle(dirich_a->name,dirich_a,&(res.dirich),in,&ierr);
if (ierr != 1) dserror("Cannot read restart data");
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   now read the data that is stored in the nodes of this processors
   partition

   a structure node as the following ARRAYS:
   actnode->sol
   actnode->sol_increment
   actnode->sol_residual
   in case of a multifield problem there exists
   actnode->sol_mf

   so we need to read an array of numnp x 4 handles first
*/
/*----------------------------------------------------------------------*/
numnp = actpart->pdis[0].numnp;
if (numnp != res.node_fdim || 4 != res.node_sdim)
    dserror("Mismatch in number of nodes on reading restart");
/*----------------------------------------- define the array of handles */
numnp = actpart->pdis[0].numnp;
res.node_handles = (long int**)CCAMALLOC(numnp*sizeof(long int*));
node_handles = res.node_handles;
res.node_handles[0] = (long int*)CCAMALLOC(4*numnp*sizeof(long int));
for (i=1; i<numnp; i++)
node_handles[i] = &(node_handles[0][i*4]);
/*------------------------------------------- read the array of handles */
pss_read_name_handle("nod_hand",&(res.node_fdim),&(res.node_sdim),&i,
                     node_handles[0],&res.handle_of_node_handles,in,&ierr);
if (ierr != 1) dserror("Cannot read restart data");
/*---------------- now we loop the nodes and each node reads his ARRAYs */
for (i=0; i<numnp; i++)
{
   actnode = actpart->pdis[0].node[i];
   /*------------------------- check for the dimensions of actnode->sol */
   pss_getdims_name_handle(actnode->sol.name,&dims[0],&dims[1],&dims[2],&(node_handles[i][0]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if ((unsigned)dims[2] != sizeof(DOUBLE)) dserror("Cannot read restart data");
   /*------------------------ redefine it, if dimension mismatch occurs */
   if (dims[0] != actnode->sol.fdim ||
       dims[1] != actnode->sol.sdim)
   amredef(&(actnode->sol),dims[0],dims[1],"DA");
   /*---------------------------------------------------------- read it */
   pss_read_array_name_handle(actnode->sol.name,&(actnode->sol),&(node_handles[i][0]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   /*-------------------------------------------------------------------*/
   /*--------------- check for the dimensions of actnode->sol_increment */
   pss_getdims_name_handle(actnode->sol_increment.name,&dims[0],&dims[1],&dims[2],&(node_handles[i][1]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if ((unsigned)dims[2] != sizeof(DOUBLE)) dserror("Cannot read restart data");
   /*------------------------ redefine it, if dimension mismatch occurs */
   if (dims[0] != actnode->sol_increment.fdim ||
       dims[1] != actnode->sol_increment.sdim)
   amredef(&(actnode->sol_increment),dims[0],dims[1],"DA");
   /*---------------------------------------------------------- read it */
   pss_read_array_name_handle(actnode->sol_increment.name,&(actnode->sol_increment),&(node_handles[i][1]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   /*-------------------------------------------------------------------*/
   /*---------------- check for the dimensions of actnode->sol_residual */
   pss_getdims_name_handle(actnode->sol_residual.name,&dims[0],&dims[1],&dims[2],&(node_handles[i][2]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if ((unsigned)dims[2] != sizeof(DOUBLE)) dserror("Cannot read restart data");
   /*------------------------ redefine it, if dimension mismatch occurs */
   if (dims[0] != actnode->sol_residual.fdim ||
       dims[1] != actnode->sol_residual.sdim)
   amredef(&(actnode->sol_residual),dims[0],dims[1],"DA");
   /*---------------------------------------------------------- read it */
   pss_read_array_name_handle(actnode->sol_residual.name,&(actnode->sol_residual),&(node_handles[i][2]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   /*---------------------- check for the dimensions of actnode->sol_mf */
   if (actnode->sol_mf.Typ == cca_DA)
   {
   pss_getdims_name_handle(actnode->sol_mf.name,&dims[0],&dims[1],&dims[2],&(node_handles[i][3]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if ((unsigned)dims[2] != sizeof(DOUBLE)) dserror("Cannot read restart data");
   /*------------------------ redefine it, if dimension mismatch occurs */
   if (dims[0] != actnode->sol_mf.fdim ||
       dims[1] != actnode->sol_mf.sdim)
   amredef(&(actnode->sol_mf),dims[0],dims[1],"DA");
   /*---------------------------------------------------------- read it */
   pss_read_array_name_handle(actnode->sol_mf.name,&(actnode->sol_mf),&(node_handles[i][3]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   }
   /*-------------------------------------------------------------------*/

   /*-------------------------------------------------------------------*/
} /* end of for (i=0; i<numnp; i++) */
/*------------------------------- delete the handle array for the nodes */
/*amdel(&(res.node_handles));*/
CCAFREE(res.node_handles[0]);
CCAFREE(res.node_handles);
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   now start reading the element data
*/
/*----------------------------------------------------------------------*/
numele = actpart->pdis[0].numele;
if (numele != res.ele_fdim || MAXRECORDPERELE != res.ele_sdim)
    dserror("Mismatch in number of elements on reading restart");
/*----------------------------------------- define the array of handles */
res.ele_handles = (long int**)CCAMALLOC(numele*sizeof(long int*));
ele_handles = res.ele_handles;
res.ele_handles[0] = (long int*)CCAMALLOC(MAXRECORDPERELE*numele*sizeof(long int));
for (i=1; i<numele; i++)
ele_handles[i] = &(ele_handles[0][i*MAXRECORDPERELE]);
/*------------------------------------------- read the array of handles */
pss_read_name_handle("ele_hand",&res.ele_fdim,&res.ele_sdim,&i,
                     ele_handles[0],&res.handle_of_ele_handles,in,&ierr);
if (ierr != 1) dserror("Cannot read restart data");
/*--------------------- now loop element and switch for type of element */
*action = read_restart;
for (i=0; i<actpart->pdis[0].numele; i++)
{
   actele = actpart->pdis[0].element[i];
   switch(actele->eltyp)/*======================= call element routines */
   {
#ifdef D_SHELL8
   case el_shell8:
      container->kstep    = 0;
      container->handsize = MAXRECORDPERELE;
      container->handles  = ele_handles[i];
      shell8(actfield,actpart,actintra,actele,
             NULL,NULL,NULL,action,container);
   break;
#endif
   case el_shell9:
      dserror("Restart 'nlnstructdyn' for shell9 not yet impl.");
   break;
   case el_brick1:
       dserror("Restart for brick not yet impl.");
   break;
#ifdef D_WALL1
   case el_wall1:
       container->handsize = MAXRECORDPERELE;
       container->handles  = ele_handles[i];
       wall1(actpart,actintra,actele,NULL,NULL,NULL,action,container);
   break;
#endif
   case el_fluid2:
       dserror("Restart for fluid2 not yet impl.");
   break;
   case el_fluid3:
       dserror("Restart for fluid3 not yet impl.");
   break;
   case el_ale3:
       dserror("Restart for ale not yet impl.");
   break;
   case el_none:
      dserror("Typ of element unknown");
   break;
   default:
      dserror("Typ of element unknown");
   }/* end of calling elements */
}
/*----------------------------------------------------------------------*/
/*----------------------------- delete the handle array of the elements */
CCAFREE(res.ele_handles[0]);
CCAFREE(res.ele_handles);
/*----------------------------------------------------------------------*/
/*
   now we have to make the arrays node->sol, node->sol_increment,
   node->sol_residual redundant for the whole field
*/
#ifdef PARALLEL
numnp = actfield->dis[0].numnp;
for (i=0; i<numnp; i++)
{
   actnode = &(actfield->dis[0].node[i]);
   sender  = actnode->proc;
   /* now we make sol */
   fdim = actnode->sol.fdim;
   MPI_Bcast(&fdim,1,MPI_INT,sender,actintra->MPI_INTRA_COMM);
   if (fdim != actnode->sol.fdim)
   amredef(&(actnode->sol),fdim,actnode->sol.sdim,"DA");
   j = actnode->sol.fdim * actnode->sol.sdim;
   MPI_Bcast(actnode->sol.a.da[0],j,MPI_DOUBLE,sender,actintra->MPI_INTRA_COMM);
   /* now we make sol_increment */
   fdim = actnode->sol_increment.fdim;
   MPI_Bcast(&fdim,1,MPI_INT,sender,actintra->MPI_INTRA_COMM);
   if (fdim != actnode->sol_increment.fdim)
   amredef(&(actnode->sol_increment),fdim,actnode->sol_increment.sdim,"DA");
   j = actnode->sol_increment.fdim * actnode->sol_increment.sdim;
   MPI_Bcast(actnode->sol_increment.a.da[0],j,MPI_DOUBLE,sender,actintra->MPI_INTRA_COMM);
   /* now we make sol_residual */
   fdim = actnode->sol_residual.fdim;
   MPI_Bcast(&fdim,1,MPI_INT,sender,actintra->MPI_INTRA_COMM);
   if (fdim != actnode->sol_residual.fdim)
   amredef(&(actnode->sol_residual),fdim,actnode->sol_residual.sdim,"DA");
   j = actnode->sol_residual.fdim * actnode->sol_residual.sdim;
   MPI_Bcast(actnode->sol_residual.a.da[0],j,MPI_DOUBLE,sender,actintra->MPI_INTRA_COMM);
   /* now we make sol_mf */
   if (actnode->sol_mf.Typ == cca_DA)
   {
   fdim = actnode->sol_mf.fdim;
   MPI_Bcast(&fdim,1,MPI_INT,sender,actintra->MPI_INTRA_COMM);
   if (fdim != actnode->sol_mf.fdim)
   amredef(&(actnode->sol_mf),fdim,actnode->sol_mf.sdim,"DA");
   j = actnode->sol_mf.fdim * actnode->sol_mf.sdim;
   MPI_Bcast(actnode->sol_mf.a.da[0],j,MPI_DOUBLE,sender,actintra->MPI_INTRA_COMM);
   }
}
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of restart_read_nlnstructdyn */


/*----------------------------------------------------------------------------------------*
 | write restart                                         ah    08/02                      |
 | of nonlinear - structural statics                                                      |
 | How to restart:                                                                        |
 | 1.) run calculation with RESTART 0 and NUMSTEP 10 in the input file                    |
 | 2.) cp xxx.pss xxx_0.pss                                                               |
 | 3.) run calculation with RESTART 9 and NUMSTEP 20 and maybe differ-                    |
 |     ent STEPSIZE, MAXITER and RESTARTEVRY in the input file                            |
 |     by starting with                                                                   |
 |     cca_seq_10.20.exe i/input.dat o/xxx restart                                        |
 | 4.) cp xxx.respss xxx_1.pss                                                            |
 | 5.) cp xxx.respss xxx.pss                                                              |
 | 6.) run calculation with RESTART 19 and NUMSTEP 30 in the input file                   |
 |     by starting with                                                                   |
 |     cca_seq_10.20.exe input.dat xxx restart                                            |
 | 7.) cp xxx.respss xxx_2.pss                                                            |
 | 8.) cp xxx.respss xxx.pss                                                              |
 | generally: restart is read from .pss and writes to .respss                             |
 |            a restarted calculation runs from step                                      |
 |            RESTART i+1 to NUMSTEP j-1                                                  |
 *----------------------------------------------------------------------------------------*/
void restart_write_nlnstructstat(STATIC_VAR     *statvar, /*-------------- static input --*/
                STANLN          *nln_data,  /*-- control variables for global NR-Iterat --*/
                FIELD           *actfield,  /*---------------------------- actual field --*/
                PARTITION       *actpart,   /*------------------------ actual partition --*/
                INTRA           *actintra,  /*---------------- actual intra comunicator --*/
                CALC_ACTION     *action,    /*---------- element action = write-restart --*/
                INT kstep,                  /*------------------------ actual load step --*/
                INT nrhs,  DIST_VECTOR *rhs,/*-- Fext processorpart of actual load step --*/
                INT nsol,  DIST_VECTOR *sol,/* solution processorpart     --"--         --*/
                INT ndis,  DIST_VECTOR *dispi,/*- displacement processorpart  --"--     --*/
                CONTAINER    *container)     /*!< contains variables defined in container.h */
{
#ifndef NO_TEXT_OUTPUT
INT                  i;                      /*-------------------------------- counter --*/
INT                  ierr;                   /*----------------------------- error-flag --*/
INT                  numnp;                  /*----------------- number of nodal points --*/
long int           **node_handles;           /*---------- handles for node- information --*/
INT                  numele;                 /*--------------------- number of elements --*/
long int           **ele_handles;            /*---------- handles for node- information --*/
char                 resname[100];           /*------- restartname = restart+stepnumber --*/
long int             longdummy;
FILE                *out;                    /*--------------- file to write restart in --*/
RESTART_STATSTRUCT   res;                    /*------ res has acces to all restart info --*/
DIST_VECTOR         *distwrite;              /*--------- array, which has to be written --*/
NODE                *actnode;                /*---------------------------- actual node --*/
ELEMENT             *actele;                 /*------------------------- actual element --*/
#ifdef DEBUG
dstrc_enter("restart_write_nlnstructstat");
#endif
#ifdef PERF
  perf_begin(74);
#endif
/*----------------------------------------------------------------------*/
out = allfiles.out_pss;
/*--------- check the step we are in and create the name "res<step>" ---*/
res.step = kstep;
sprintf(resname,"sres%d",res.step);
/*------------- write the structure statvar to the restart structure ---*/
res.statvar = *statvar;
/*------------ write the structure nln_data to the restart structure ---*/
res.nln_data = *nln_data;
/* write all distributed vectors and store the handles associated with them */
res.dist_vec_rhs[0]   = nrhs;
res.dist_vec_sol[0]   = nsol;
res.dist_vec_dispi[0] = ndis;

/*------------------------------------------------------ write the *rhs */
 for (i=0; i<nrhs; i++)
 {
    distwrite = &(rhs[i]);
    pss_write_array(&(distwrite->vec),&(res.dist_vec_rhs[i+1]),out,&ierr);
    if (ierr != 1) dserror("Error writing restart data");
 }
 /*------------------------------------------------------ write the *sol */
 for (i=0; i<nsol; i++)
 {
    distwrite = &(sol[i]);
    pss_write_array(&(distwrite->vec),&(res.dist_vec_sol[i+1]),out,&ierr);
    if (ierr != 1) dserror("Error writing restart data");
 }
/*------------------------------------------------------ write the *displ */
 for (i=0; i<ndis; i++)
 {
    distwrite = &(dispi[i]);
    pss_write_array(&(distwrite->vec),&(res.dist_vec_dispi[i+1]),out,&ierr);
    if (ierr != 1) dserror("Error writing restart data");
 }
/*----------------------------------------------------------------------*/
/*
   now write the data, that is stored in the nodes of this processor's
   partition

   a node has 3 ARRAYS:
   actnode->sol
   actnode->sol_increment
   actnode->sol_residual

   so we need to store an array of numnp x 3 handles
*/
/*----------------------------------------------------------------------*/
numnp = actpart->pdis[0].numnp;
res.node_handles = (long int**)CCAMALLOC(numnp*sizeof(long int*));
node_handles = res.node_handles;
res.node_handles[0] = (long int*)CCAMALLOC(3*numnp*sizeof(long int));
for (i=1; i<numnp; i++)
node_handles[i] = &(node_handles[0][i*3]);
/*----------------------------------------------------------------------*/
/* now we loop the nodes on the partition and each node writes his ARRAYs */
for (i=0; i<numnp; i++)
{
   actnode = actpart->pdis[0].node[i];
   pss_write_array(&(actnode->sol),&(node_handles[i][0]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   pss_write_array(&(actnode->sol_increment),&(node_handles[i][1]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   pss_write_array(&(actnode->sol_residual),&(node_handles[i][2]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
}
/*
   all nodal data is written, so we now write the node_handles and store
   the handle to this in res.handle_of_node_handles
*/
pss_write("nod_hand",numnp,3,sizeof(long int),node_handles[0],&(res.handle_of_node_handles),out,&ierr);
if (ierr != 1) dserror("Error writing restart data");
/*----------------- delete the res.node_handles but keep the dimensions */
res.node_fdim = numnp;
res.node_sdim = 3;
CCAFREE(res.node_handles[0]);
CCAFREE(res.node_handles);
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   we come to the most difficult part, the writing of the element data.
   One element is assumed not to have more then MAXRECORDPERELE different
   records to write of arbitary type

   The element is called with the element action action = write_restart
   Again, all elements to the processors partition are written
*/
/*----------------------------------------------------------------------*/
numele = actpart->pdis[0].numele;
res.ele_handles = (long int**)CCAMALLOC(numele*sizeof(long int*));
ele_handles = res.ele_handles;
res.ele_handles[0] = (long int*)CCAMALLOC(MAXRECORDPERELE*numele*sizeof(long int));
for (i=1; i<numele; i++)
ele_handles[i] = &(ele_handles[0][i*MAXRECORDPERELE]);
/*--------------------- now loop element and switch for type of element */
*action = write_restart;
for (i=0; i<actpart->pdis[0].numele; i++)
{
   actele = actpart->pdis[0].element[i];
   switch(actele->eltyp)/*======================= call element routines */
   {
#ifdef D_SHELL8
   case el_shell8:
      container->kstep    = 0;
      container->handsize = MAXRECORDPERELE;
      container->handles  = ele_handles[i];
      shell8(actfield,actpart,actintra,actele,
             NULL,NULL,NULL,action,container);
   break;
#endif

#ifdef D_SHELL9
   case el_shell9:
      container->kstep    = 0;
      container->handsize = MAXRECORDPERELE;
      container->handles  = ele_handles[i];
      shell9(actfield,actintra,actele,
             NULL,NULL,NULL,action,container);
   break;
#endif
   case el_brick1:
   break;
#ifdef D_WALL1
   case el_wall1:
      container->handsize = MAXRECORDPERELE;
      container->handles  = ele_handles[i];
      wall1(actpart,actintra,actele,
            NULL,NULL,NULL,action,container);
   break;
#endif
   case el_fluid2:
       dserror("Restart for fluid2 not yet impl.");
   break;
   case el_fluid3:
       dserror("Restart for fluid3 not yet impl.");
   break;
   case el_ale3:
       dserror("Restart for ale not yet impl.");
   break;
#ifdef D_BEAM3
   case el_beam3:
       container->handsize = MAXRECORDPERELE;
       container->handles = ele_handles[i];
       beam3(actfield,actpart,actintra,actele,
             NULL,NULL,NULL,action,container);
   break;
#endif
   case el_interf:

   break;
   case el_wallge:

   break;
   case el_none:
      dserror("Typ of element unknown");
   break;
   default:
      dserror("Typ of element unknown");
   }/* end of calling elements */
}
/*----------------------------------------------------------------------*/
/*
   all ele data is written, so write the ele_handles and store the handle to
   it in res.handle_of_ele_handles
*/
res.ele_fdim = numele;
res.ele_sdim = MAXRECORDPERELE;
pss_write("ele_hand",numele,MAXRECORDPERELE,sizeof(long int),ele_handles[0],&(res.handle_of_ele_handles),out,&ierr);
if (ierr != 1) dserror("Error writing restart data");
/*------------------ delete the res.ele_handles but keep the dimensions */
CCAFREE(res.ele_handles[0]);
CCAFREE(res.ele_handles);
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   the only thing to do is now write the RESTART_STRUCT itself with
   its unique name resname = "res<step>"
   NOTE:
   names are limited to 9 characters, so a step larger then sres99999
   can not be restarted at the moment !!!!
*/
pss_write(resname,1,1,sizeof(RESTART_STATSTRUCT),&res,&longdummy,out,&ierr);
if (ierr != 1) dserror("Error writing restart data");
/*----------------------------------------------------------------------*/
/*
#ifdef DEBUG
pss_status_to_err();
#endif
*/
/*----------------------------------------------------------------------*/
/*
   now write a notice to the err file, that this retsrat step was
   successfully written
*/
if (res.step != 0)
{
   fprintf(allfiles.out_err,"===========================================\n");
   fprintf(allfiles.out_err,"RESTART MESSAGE\n");
   fprintf(allfiles.out_err,"In step %d structure restart data was written\n",res.step);
   fprintf(allfiles.out_err,"Calculation can be restarted using \n");
   fprintf(allfiles.out_err,"RESTART      %d \n",res.step);
   fprintf(allfiles.out_err,"in the input file\n");
   fprintf(allfiles.out_err,"===========================================\n");
   fflush(allfiles.out_err);
}


/*----------------------------------------------------------------------*/
#ifdef PERF
  perf_end(74);
#endif
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of restart_write_nlnstructstat*/




/*-------------------------------------------------------------------------------------------*
 | read restart                                                                     ah 08/02 |
 | of nln-structural statics                                                                 |
 *-------------------------------------------------------------------------------------------*/
void restart_read_nlnstructstat(INT restart,   /* restart step */
                STATIC_VAR     *statvar,       /* static input */
                STANLN          *nln_data,     /* control variables for global NR-Iterat */
                FIELD           *actfield,     /* actual field */
                PARTITION       *actpart,      /* actual partition */
                INTRA           *actintra,     /* actual intra comunicator */
                CALC_ACTION     *action,       /* element action = write-restart */
                INT nrhs,  DIST_VECTOR *rhs,   /* Fext processorpart of actual load step */
                INT nsol,  DIST_VECTOR *sol,   /* solution processorpart     --"--       */
                INT ndis,  DIST_VECTOR *dispi, /* displacement processorpart  --"--      */
                CONTAINER    *container)       /*!< contains variables defined in container.h */
{
#ifndef NO_TEXT_OUTPUT
INT                  i;                      /* counters */
INT                  ierr;                     /* error-flag */
long int             reshandle;                /* restart handle */
INT                  byte;                     /* length of binary file */
INT                  dims[3];                  /* dimensions of actnode->sol */
INT                  numnp;                    /* number of nodal points */
long int           **node_handles;             /* handles for node- information */
INT                  numele;                   /* number of elements */
long int           **ele_handles;              /* handles for element - information */
char                 resname[100];             /* restartname = restart+stepnumber */
FILE                *in;                       /* file to be read in  */
RESTART_STATSTRUCT  res;                       /* res has acces to all restart info */
DIST_VECTOR         *distread;                 /* array, which has to be read */
NODE                *actnode;                  /* actual node */
ELEMENT             *actele;                   /* actual element */

#ifdef PARALLEL
INT                  j;                      /* counters */
INT                  fdim;                     /* fdim = actfield->actnode->sol.fdim */
INT                  sender;                   /* sender in parallel comunication */
#endif

#ifdef DEBUG
dstrc_enter("restart_read_nlnstructstat");
#endif
/*----------------------------------------------------------------------*/
in = allfiles.in_pss;
/*----------------------------------- check the step that shall be read */
sprintf(resname,"sres%d",restart);
pss_chck(resname,&reshandle,in,&ierr);
if (ierr != 1) dserror("Cannot restart, step doesn't exist in pss-file");
/*----------------------------- the structure res exists, so we read it */
pss_read_name_handle(resname,&i,&i,&byte,&res,&reshandle,in,&ierr);
if (ierr != 1) dserror("Restart structure exists, but cannot read it");
/*------------------------------- write the structure res.sdyn to *sdyn */
*statvar   = res.statvar;
/*--------------------------- write the structure res.dynvar to *dynvar */
*nln_data = res.nln_data;
/*---------------------------------------------------- make some checks */
if (res.dist_vec_rhs[0]   != nrhs) dserror("Mismatch in dimensions in restart");
if (res.dist_vec_sol[0]   != nsol) dserror("Mismatch in dimensions in restart");
if (res.dist_vec_dispi[0] != ndis) dserror("Mismatch in dimensions in restart");
/*---------------------------------------- read the distributed vectors */
/*------------------------------------------------------------ read rhs */
for (i=0; i<nrhs; i++)
{
   distread = &(rhs[i]);
   pss_read_array_name_handle(distread->vec.name,&(distread->vec),&(res.dist_vec_rhs[i+1]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
}
/*------------------------------------------------------------ read sol */
for (i=0; i<nsol; i++)
{
   distread = &(sol[i]);
   pss_read_array_name_handle(distread->vec.name,&(distread->vec),&(res.dist_vec_sol[i+1]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
}
/*------------------------------------------------------------ read dispi */
for (i=0; i<ndis; i++)
{
   distread = &(dispi[i]);
   pss_read_array_name_handle(distread->vec.name,&(distread->vec),&(res.dist_vec_dispi[i+1]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
}
/*----------------------------------------------------------------------*/
/*
   now read the data that is stored in the nodes of this processors
   partition

   a node has 3 ARRAYS:
   actnode->sol
   actnode->sol_increment
   actnode->sol_residual

   so we need to read an array of numnp x 3 handles first
*/
/*----------------------------------------------------------------------*/
numnp = actpart->pdis[0].numnp;
if (numnp != res.node_fdim || 3 != res.node_sdim)
    dserror("Mismatch in number of nodes on reading restart");
/*----------------------------------------- define the array of handles */
numnp = actpart->pdis[0].numnp;
res.node_handles = (long int**)CCAMALLOC(numnp*sizeof(long int*));
node_handles = res.node_handles;
res.node_handles[0] = (long int*)CCAMALLOC(3*numnp*sizeof(long int));
for (i=1; i<numnp; i++)
node_handles[i] = &(node_handles[0][i*3]);
/*------------------------------------------- read the array of handles */
pss_read_name_handle("nod_hand",&(res.node_fdim),&(res.node_sdim),&i,
                     node_handles[0],&res.handle_of_node_handles,in,&ierr);
if (ierr != 1) dserror("Cannot read restart data");
/*---------------- now we loop the nodes and each node reads his ARRAYs */
for (i=0; i<numnp; i++)
{
   actnode = actpart->pdis[0].node[i];
   /*------------------------- check for the dimensions of actnode->sol */
   pss_getdims_name_handle(actnode->sol.name,&dims[0],&dims[1],&dims[2],&(node_handles[i][0]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if ((unsigned)dims[2] != sizeof(DOUBLE)) dserror("Cannot read restart data");
   /*------------------------ redefine it, if dimension mismatch occurs */
   if (dims[0] != actnode->sol.fdim ||
       dims[1] != actnode->sol.sdim)
   amredef(&(actnode->sol),dims[0],dims[1],"DA");
   /*---------------------------------------------------------- read it */
   pss_read_array_name_handle(actnode->sol.name,&(actnode->sol),&(node_handles[i][0]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   /*-------------------------------------------------------------------*/
   /*--------------- check for the dimensions of actnode->sol_increment */
   pss_getdims_name_handle(actnode->sol_increment.name,&dims[0],&dims[1],&dims[2],&(node_handles[i][1]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if ((unsigned)dims[2] != sizeof(DOUBLE)) dserror("Cannot read restart data");
   /*------------------------ redefine it, if dimension mismatch occurs */
   if (dims[0] != actnode->sol_increment.fdim ||
       dims[1] != actnode->sol_increment.sdim)
   amredef(&(actnode->sol_increment),dims[0],dims[1],"DA");
   /*---------------------------------------------------------- read it */
   pss_read_array_name_handle(actnode->sol_increment.name,&(actnode->sol_increment),&(node_handles[i][1]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   /*-------------------------------------------------------------------*/
   /*---------------- check for the dimensions of actnode->sol_residual */
   pss_getdims_name_handle(actnode->sol_residual.name,&dims[0],&dims[1],&dims[2],&(node_handles[i][2]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if ((unsigned)dims[2] != sizeof(DOUBLE)) dserror("Cannot read restart data");
   /*------------------------ redefine it, if dimension mismatch occurs */
   if (dims[0] != actnode->sol_residual.fdim ||
       dims[1] != actnode->sol_residual.sdim)
   amredef(&(actnode->sol_residual),dims[0],dims[1],"DA");
   /*---------------------------------------------------------- read it */
   pss_read_array_name_handle(actnode->sol_residual.name,&(actnode->sol_residual),&(node_handles[i][2]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   /*-------------------------------------------------------------------*/
} /* end of for (i=0; i<numnp; i++) */
/*------------------------------- delete the handle array for the nodes */
/*amdel(&(res.node_handles));*/
CCAFREE(res.node_handles[0]);
CCAFREE(res.node_handles);
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   now start reading the element data
*/
/*----------------------------------------------------------------------*/
numele = actpart->pdis[0].numele;
if (numele != res.ele_fdim || MAXRECORDPERELE != res.ele_sdim)
    dserror("Mismatch in number of elements on reading restart");
/*----------------------------------------- define the array of handles */
res.ele_handles = (long int**)CCAMALLOC(numele*sizeof(long int*));
ele_handles = res.ele_handles;
res.ele_handles[0] = (long int*)CCAMALLOC(MAXRECORDPERELE*numele*sizeof(long int));
for (i=1; i<numele; i++)
ele_handles[i] = &(ele_handles[0][i*MAXRECORDPERELE]);
/*------------------------------------------- read the array of handles */
pss_read_name_handle("ele_hand",&res.ele_fdim,&res.ele_sdim,&i,
                     ele_handles[0],&res.handle_of_ele_handles,in,&ierr);
if (ierr != 1) dserror("Cannot read restart data");
/*--------------------- now loop element and switch for type of element */
*action = read_restart;
for (i=0; i<actpart->pdis[0].numele; i++)
{
   actele = actpart->pdis[0].element[i];
   switch(actele->eltyp)/*======================= call element routines */
   {
#ifdef D_SHELL8
   case el_shell8:
      container->kstep    = 0;
      container->handsize = MAXRECORDPERELE;
      container->handles  = ele_handles[i];
      shell8(actfield,actpart,actintra,actele,
             NULL,NULL,NULL,action,container);
   break;
#endif

#ifdef D_SHELL9
   case el_shell9:
      container->kstep    = 0;
      container->handsize = MAXRECORDPERELE;
      container->handles  = ele_handles[i];
      shell9(actfield,actintra,actele,
             NULL,NULL,NULL,action,container);
   break;
#endif
   case el_brick1:
   break;
#ifdef D_WALL1
   case el_wall1:
       container->handsize = MAXRECORDPERELE;
       container->handles  = ele_handles[i];
       wall1(actpart,actintra,actele,NULL,NULL,NULL,action,container);
   break;
#endif
   case el_fluid2:
       dserror("Restart for fluid2 not yet impl.");
   break;
   case el_fluid3:
       dserror("Restart for fluid3 not yet impl.");
   break;
   case el_ale3:
       dserror("Restart for ale not yet impl.");
   break;
#ifdef D_BEAM3
   case el_beam3:
       container->kstep = 0;
       container->handsize = MAXRECORDPERELE;
       container->handles = ele_handles[i];
       beam3(actfield,actpart,actintra,actele,
             NULL,NULL,NULL,action,container);
   break;
#endif
   case el_none:
      dserror("Typ of element unknown");
   break;
   default:
      dserror("Typ of element unknown");
   }/* end of calling elements */
}
/*----------------------------------------------------------------------*/
/*----------------------------- delete the handle array of the elements */
CCAFREE(res.ele_handles[0]);
CCAFREE(res.ele_handles);
/*----------------------------------------------------------------------*/
/*
   now we have to make the arrays node->sol, node->sol_increment,
   node->sol_residual redundant for the whole field
*/
#ifdef PARALLEL
numnp = actfield->dis[0].numnp;
for (i=0; i<numnp; i++)
{
   actnode = &(actfield->dis[0].node[i]);
   sender  = actnode->proc;
   /* now we make sol */
   fdim = actnode->sol.fdim;
   MPI_Bcast(&fdim,1,MPI_INT,sender,actintra->MPI_INTRA_COMM);
   if (fdim != actnode->sol.fdim)
   amredef(&(actnode->sol),fdim,actnode->sol.sdim,"DA");
   j = actnode->sol.fdim * actnode->sol.sdim;
   MPI_Bcast(actnode->sol.a.da[0],j,MPI_DOUBLE,sender,actintra->MPI_INTRA_COMM);
   /* now we make sol_increment */
   fdim = actnode->sol_increment.fdim;
   MPI_Bcast(&fdim,1,MPI_INT,sender,actintra->MPI_INTRA_COMM);
   if (fdim != actnode->sol_increment.fdim)
   amredef(&(actnode->sol_increment),fdim,actnode->sol_increment.sdim,"DA");
   j = actnode->sol_increment.fdim * actnode->sol_increment.sdim;
   MPI_Bcast(actnode->sol_increment.a.da[0],j,MPI_DOUBLE,sender,actintra->MPI_INTRA_COMM);
   /* now we make sol_residual */
   fdim = actnode->sol_residual.fdim;
   MPI_Bcast(&fdim,1,MPI_INT,sender,actintra->MPI_INTRA_COMM);
   if (fdim != actnode->sol_residual.fdim)
   amredef(&(actnode->sol_residual),fdim,actnode->sol_residual.sdim,"DA");
   j = actnode->sol_residual.fdim * actnode->sol_residual.sdim;
   MPI_Bcast(actnode->sol_residual.a.da[0],j,MPI_DOUBLE,sender,actintra->MPI_INTRA_COMM);
}
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /*end of restart_read_nlnstructstat */
/*!---------------------------------------------------------------------
\brief write restart data of fluid field

<pre>                                                         genk 09/03

write restart data of FLUID field from pss-file

</pre>
\param   *fdyn	     FLUID_DYNAMIC    (i)
\param   *actfield   FIELD            (i)   actual field
\param   *actpart    PARTITION        (i)   actual partition
\param   *actintra   INTRA            (i)   actual intra-communicator
\param   *action     CALC_ACTION      (i)
\param   *container  CONTAINER        (i)

\return void

------------------------------------------------------------------------*/
void restart_write_fluiddyn(FLUID_DYNAMIC   *fdyn,
                            FIELD	    *actfield,
                            PARTITION	    *actpart,
                            INTRA	    *actintra,
			    CALC_ACTION     *action,
			    CONTAINER       *container)
{
#ifndef NO_TEXT_OUTPUT
INT                  i;
INT                  ierr;
INT                  numnp;
long int           **node_handles;
INT                  numele;
long int           **ele_handles;
char                 resname[100];
char                 iposname[100];
long int             longdummy;
FILE                *out;
RESTART_DYNFLUID     res;
NODE                *actnode;
ELEMENT             *actele;

#ifdef DEBUG
dstrc_enter("restart_write_fluiddyn");
#endif
#ifdef PERF
  perf_begin(74);
#endif
/*----------------------------------------------------------------------*/
out = allfiles.out_pss;
/*-------- check the step we are in and create the name "res<step>" */
res.step = fdyn->step;
res.time = fdyn->acttime;
sprintf(resname,"fres%d",res.step);
sprintf(iposname,"ipos%d",res.step);

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   now write the data, that is stored in the nodes of this processor's
   partition

   for a fluid node two arrays are important:
   actnode->sol
   actnode->sol_increment
   additionally there may be a multifield array
   actnode->sol_mf;

   so we need to store an array of numnp x 3 handles
*/
/*----------------------------------------------------------------------*/
numnp = actpart->pdis[0].numnp;
res.node_handles = (long int**)CCAMALLOC(numnp*sizeof(long int*));
if (!res.node_handles) dserror("Allocation of memory failed");
node_handles = res.node_handles;
res.node_handles[0] = (long int*)CCAMALLOC(3*numnp*sizeof(long int));
if (!res.node_handles[0]) dserror("Allocation of memory failed");
for (i=1; i<numnp; i++)
node_handles[i] = &(node_handles[0][i*3]);
/*----------------------------------------------------------------------*/
/* now we loop the nodes on the partition and each node writes his ARRAYs */
for (i=0; i<numnp; i++)
{
   actnode = actpart->pdis[0].node[i];
   pss_write_array(&(actnode->sol),&(node_handles[i][0]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   pss_write_array(&(actnode->sol_increment),&(node_handles[i][1]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   if (actnode->sol_mf.Typ == cca_DA)
   {
   pss_write_array(&(actnode->sol_mf),&(node_handles[i][2]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   }
}
/*
   all nodal data is written, so we now write the node_handles and store
   the handle to this in res.handle_of_node_handles
*/
pss_write("nod_hand",numnp,3,sizeof(long int),node_handles[0],&(res.handle_of_node_handles),out,&ierr);
if (ierr != 1) dserror("Error writing restart data");
/*----------------- delete the res.node_handles but keep the dimensions */
res.node_fdim = numnp;
res.node_sdim = 3;
CCAFREE(res.node_handles[0]);
CCAFREE(res.node_handles);

/*----------------------------------------------------------------------*/
/*
   we come to the most difficult part, the writing of the element data.
   One element is assumed not to have more then 5 different records to
   write of arbitary type

   The element is called with the element action action = write_restart
   Again, all elements to the processors partition are written

   At the moment this is only necessary for a free surface with surface
   tension effects. Then we have to save the nodal curvature.
*/
/*----------------------------------------------------------------------*/
numele = actpart->pdis[0].numele;
res.ele_handles = (long int**)CCAMALLOC(numele*sizeof(long int*));
ele_handles = res.ele_handles;
res.ele_handles[0] = (long int*)CCAMALLOC(5*numele*sizeof(long int));
for (i=1; i<numele; i++)
ele_handles[i] = &(ele_handles[0][i*5]);
/*------------------- now loop element and switch for type of element */
*action = write_restart;
for (i=0; i<actpart->pdis[0].numele; i++)
{
   actele = actpart->pdis[0].element[i];
   switch(actele->eltyp)/*===================== call element routines */
   {
#ifdef D_FLUID2
   case el_fluid2:
      container->handsize = 5;
      container->handles  = ele_handles[i];
      fluid2(actpart,actintra,actele,NULL,
             NULL,NULL,NULL,NULL,
             action,NULL,NULL,container);
   break;
#endif

#ifdef D_FLUID3
   case el_fluid3:
      container->handsize = 5;
      container->handles  = ele_handles[i];
      fluid3(actpart,actintra,actele,
             NULL,NULL,NULL,NULL,
             action,NULL,NULL,container);
   break;
#endif

#ifdef D_FLUID3_F
   case el_fluid3_fast:
      container->handsize = 5;
      container->handles  = ele_handles[i];
      fluid3_fast(actpart,actintra,&actele,
             NULL,NULL,
             NULL,NULL,
             action,NULL,NULL,container,1);
   break;
#endif
   case el_none:
      dserror("Typ of element unknown");
   break;

   default:
      dserror("Typ of element unknown");
   }/* end of calling elements */
}
/*----------------------------------------------------------------------*/
/*
   all ele data is written, so write the ele_handles and store the handle to
   it in res.handle_of_ele_handles
*/
res.ele_fdim = numele;
res.ele_sdim = 5;
pss_write("ele_hand",numele,5,sizeof(long int),ele_handles[0],&(res.handle_of_ele_handles),out,&ierr);
if (ierr != 1) dserror("Error writing restart data");
/*---------------- delete the res.ele_handles but keep the dimensions */
CCAFREE(res.ele_handles[0]);
CCAFREE(res.ele_handles);

/*----------------------------------------------------------------------*/
/*
   the only thing to do is now write the RESTART_DYNFLUID itself with
   its unique name resname = "fres<step>"
   NOTE:
   names are limited to 9 characters, so a step larger then fres99999
   can not be restarted at the moment !!!!
*/
pss_write(iposname,1,1,sizeof(ARRAY_POSITION),&actfield->dis[0].ipos,&longdummy,out,&ierr);
if (ierr != 1) dserror("Error writing restart data");
pss_write(resname,1,1,sizeof(RESTART_DYNFLUID),&res,&longdummy,out,&ierr);
if (ierr != 1) dserror("Error writing restart data");
/*----------------------------------------------------------------------*/
/*
#ifdef DEBUG
pss_status_to_err();
#endif
*/
/*----------------------------------------------------------------------*/
/*
   now write a notice to the err file, that this retsrat step was
   successfully written
*/
if (res.step != 0)
{
   fprintf(allfiles.out_err,"===========================================\n");
   fprintf(allfiles.out_err,"RESTART MESSAGE\n");
   fprintf(allfiles.out_err,"In step %d fluid restart data were written\n",res.step);
   fprintf(allfiles.out_err,"Calculation can be restarted using \n");
   fprintf(allfiles.out_err,"RESTART      %d \n",res.step);
   fprintf(allfiles.out_err,"in the input file\n");
   fprintf(allfiles.out_err,"===========================================\n");
   fflush(allfiles.out_err);
}
/*----------------------------------------------------------------------*/
#ifdef PERF
  perf_end(74);
#endif
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of restart_write_fluiddyn */

/*!---------------------------------------------------------------------
\brief read restart data of ale field

<pre>                                                         genk 09/03

read restart data of FLUID field from pss file

</pre>
\param    restart    INT              (i)   restart step
\param   *fdyn	     FLUID_DYNAMIC    (i)
\param   *actfield   FIELD            (i)   actual field
\param   *actpart    PARTITION        (i)   actual partition
\param   *actintra   INTRA            (i)   actual intra-communicator
\param   *action     CALC_ACTION      (i)
\param   *container  CONTAINER        (i)

\return void

------------------------------------------------------------------------*/
void restart_read_fluiddyn(INT restart,
                           FLUID_DYNAMIC   *fdyn,
                           FIELD	   *actfield,
                           PARTITION	   *actpart,
                           INTRA	   *actintra,
			   CALC_ACTION     *action,
			   CONTAINER       *container)
{
#ifndef NO_TEXT_OUTPUT
INT                  i;
INT                  ierr;
long int             reshandle;
INT                  byte;
INT                  dims[3];
INT                  numnp;
INT                  numele;
long int           **ele_handles;
long int           **node_handles;
char                 resname[100];
char                 iposname[100];
FILE                *in;
RESTART_DYNFLUID     res;
NODE                *actnode;
ELEMENT             *actele;

#ifdef PARALLEL
INT                  j;
INT                  fdim;
INT                  sender;
#endif

#ifdef DEBUG
dstrc_enter("restart_read_fluiddyn");
#endif
/*----------------------------------------------------------------------*/
in = allfiles.in_pss;
/*----------------------------------- check the step that shall be read */
sprintf(resname,"fres%d",restart);
sprintf(iposname,"ipos%d",restart);
pss_chck(resname,&reshandle,in,&ierr);
if (ierr != 1) dserror("Cannot restart, step doesn't exist in pss-file");
/*----------------------------- the structure res exists, so we read it */
pss_read_name_handle(resname,&i,&i,&byte,&res,&reshandle,in,&ierr);
if (ierr != 1) dserror("Restart structure exists, but cannot read it");
fdyn->step=res.step;
fdyn->acttime=res.time;
/*------------------------------------------------------- get positions */
pss_chck(resname,&reshandle,in,&ierr);
if (ierr != 1) dserror("Cannot restart, step doesn't exist in pss-file");
pss_chck(iposname,&reshandle,in,&ierr);
/*---------------------------- the structure ipos exists, so we read it */
pss_read_name_handle(iposname,&i,&i,&byte,&actfield->dis[0].ipos,&reshandle,in,&ierr);
if (ierr != 1) dserror("Restart structure exists, but cannot read it");


/*----------------------------------------------------------------------*/
/*
   now read the data that is stored in the nodes of this processors
   partition

   for a fluid node two arrays are important:
   actnode->sol
   actnode->sol_increment
   additionally there may be a multifield array
   actnode->sol_mf;

   so we need to read an array of numnp x 3 handles first
*/
/*----------------------------------------------------------------------*/
numnp = actpart->pdis[0].numnp;
if (numnp != res.node_fdim || 3 != res.node_sdim)
    dserror("Mismatch in number of nodes on reading restart");
/*----------------------------------------- define the array of handles */
numnp = actpart->pdis[0].numnp;
res.node_handles = (long int**)CCAMALLOC(numnp*sizeof(long int*));
node_handles = res.node_handles;
res.node_handles[0] = (long int*)CCAMALLOC(3*numnp*sizeof(long int));
for (i=1; i<numnp; i++)
node_handles[i] = &(node_handles[0][i*3]);
/*------------------------------------------- read the array of handles */
pss_read_name_handle("nod_hand",&(res.node_fdim),&(res.node_sdim),&i,
                     node_handles[0],&res.handle_of_node_handles,in,&ierr);
if (ierr != 1) dserror("Cannot read restart data");
/*---------------- now we loop the nodes and each node reads his ARRAYs */
for (i=0; i<numnp; i++)
{
   actnode = actpart->pdis[0].node[i];
   /*------------------------- check for the dimensions of actnode->sol */
   pss_getdims_name_handle(actnode->sol.name,&dims[0],&dims[1],&dims[2],&(node_handles[i][0]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if ((unsigned)dims[2] != sizeof(DOUBLE)) dserror("Cannot read restart data");
   /*------------------------ redefine it, if dimension mismatch occurs */
   if (dims[0] != actnode->sol.fdim ||
       dims[1] != actnode->sol.sdim)
   amredef(&(actnode->sol),dims[0],dims[1],"DA");
   /*---------------------------------------------------------- read it */
   pss_read_array_name_handle(actnode->sol.name,&(actnode->sol),&(node_handles[i][0]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   /*-------------------------------------------------------------------*/
   /*--------------- check for the dimensions of actnode->sol_increment */
   pss_getdims_name_handle(actnode->sol_increment.name,&dims[0],&dims[1],&dims[2],&(node_handles[i][1]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if ((unsigned)dims[2] != sizeof(DOUBLE)) dserror("Cannot read restart data");
   /*------------------------ redefine it, if dimension mismatch occurs */
   if (dims[0] != actnode->sol_increment.fdim ||
       dims[1] != actnode->sol_increment.sdim)
   amredef(&(actnode->sol_increment),dims[0],dims[1],"DA");
   /*---------------------------------------------------------- read it */
   pss_read_array_name_handle(actnode->sol_increment.name,&(actnode->sol_increment),&(node_handles[i][1]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   /*-------------------------------------------------------------------*/
   /*---------------- check for the dimensions of actnode->sol_mf */
   if (actnode->sol_mf.Typ == cca_DA)
   {
   pss_getdims_name_handle(actnode->sol_mf.name,&dims[0],&dims[1],&dims[2],&(node_handles[i][2]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if ((unsigned)dims[2] != sizeof(DOUBLE)) dserror("Cannot read restart data");
   /*------------------------ redefine it, if dimension mismatch occurs */
   if (dims[0] != actnode->sol_mf.fdim ||
       dims[1] != actnode->sol_mf.sdim)
   amredef(&(actnode->sol_mf),dims[0],dims[1],"DA");
   /*---------------------------------------------------------- read it */
   pss_read_array_name_handle(actnode->sol_mf.name,&(actnode->sol_mf),&(node_handles[i][2]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   }
   /*-------------------------------------------------------------------*/
} /* end of for (i=0; i<numnp; i++) */
/*------------------------------- delete the handle array for the nodes */
/*amdel(&(res.node_handles));*/
CCAFREE(res.node_handles[0]);
CCAFREE(res.node_handles);
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   now start reading the element data
*/
/*----------------------------------------------------------------------*/
numele = actpart->pdis[0].numele;
if (numele != res.ele_fdim || 5 != res.ele_sdim)
    dserror("Mismatch in number of elements on reading restart");
/*----------------------------------------- define the array of handles */
res.ele_handles = (long int**)CCAMALLOC(numele*sizeof(long int*));
ele_handles = res.ele_handles;
res.ele_handles[0] = (long int*)CCAMALLOC(5*numele*sizeof(long int));
for (i=1; i<numele; i++)
ele_handles[i] = &(ele_handles[0][i*5]);
/*------------------------------------------- read the array of handles */
pss_read_name_handle("ele_hand",&res.ele_fdim,&res.ele_sdim,&i,
                     ele_handles[0],&res.handle_of_ele_handles,in,&ierr);
if (ierr != 1) dserror("Cannot read restart data");
/*--------------------- now loop element and switch for type of element */
*action = read_restart;
for (i=0; i<actpart->pdis[0].numele; i++)
{
   actele = actpart->pdis[0].element[i];
   switch(actele->eltyp)/*======================= call element routines */
   {
   case el_fluid2:
#ifdef D_FLUID2
      container->handsize = 5;
      container->handles  = ele_handles[i];
      fluid2(actpart,actintra,actele,NULL,
             NULL,NULL,NULL,NULL,
             action,NULL,NULL,container);
#endif
   break;
   case el_fluid3:
#ifdef D_FLUID3
      container->handsize = 5;
      container->handles  = ele_handles[i];
      fluid3(actpart,actintra,actele,
             NULL,NULL,NULL,NULL,
             action,NULL,NULL,container);
#endif
   break;
   case el_none:
      dserror("Typ of element unknown");
   break;
   default:
      dserror("Typ of element unknown");
   }/* end of calling elements */
}
/*----------------------------------------------------------------------*/
/*----------------------------- delete the handle array of the elements */
CCAFREE(res.ele_handles[0]);
CCAFREE(res.ele_handles);
/*----------------------------------------------------------------------*/
/*
   now we have to make the arrays node->sol, node->sol_increment,
   node->sol_mf redundant for the whole field
*/
#ifdef PARALLEL
numnp = actfield->dis[0].numnp;
for (i=0; i<numnp; i++)
{
   actnode = &(actfield->dis[0].node[i]);
   sender  = actnode->proc;
   /* now we make sol */
   fdim = actnode->sol.fdim;
   MPI_Bcast(&fdim,1,MPI_INT,sender,actintra->MPI_INTRA_COMM);
   if (fdim != actnode->sol.fdim)
   amredef(&(actnode->sol),fdim,actnode->sol.sdim,"DA");
   j = actnode->sol.fdim * actnode->sol.sdim;
   MPI_Bcast(actnode->sol.a.da[0],j,MPI_DOUBLE,sender,actintra->MPI_INTRA_COMM);
   /* now we make sol_increment */
   fdim = actnode->sol_increment.fdim;
   MPI_Bcast(&fdim,1,MPI_INT,sender,actintra->MPI_INTRA_COMM);
   if (fdim != actnode->sol_increment.fdim)
   amredef(&(actnode->sol_increment),fdim,actnode->sol_increment.sdim,"DA");
   j = actnode->sol_increment.fdim * actnode->sol_increment.sdim;
   MPI_Bcast(actnode->sol_increment.a.da[0],j,MPI_DOUBLE,sender,actintra->MPI_INTRA_COMM);
   /* now we make sol_mf */
   if (actnode->sol_mf.Typ == cca_DA)
   {
   fdim = actnode->sol_mf.fdim;
   MPI_Bcast(&fdim,1,MPI_INT,sender,actintra->MPI_INTRA_COMM);
   if (fdim != actnode->sol_mf.fdim)
   amredef(&(actnode->sol_mf),fdim,actnode->sol_mf.sdim,"DA");
   j = actnode->sol_mf.fdim * actnode->sol_mf.sdim;
   MPI_Bcast(actnode->sol_mf.a.da[0],j,MPI_DOUBLE,sender,actintra->MPI_INTRA_COMM);
   }
}
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of restart_read_fluiddyn */

/*!---------------------------------------------------------------------
\brief write restart data of ale field

<pre>                                                         genk 09/03

write restart data of ALE field from pss file

</pre>
\param   *adyn	     ALE_DYNAMIC      (i)
\param   *actfield   FIELD            (i)   actual field
\param   *actpart    PARTITION        (i)   actual partition
\param   *actintra   INTRA            (i)   actual intra-communicator

\return void

------------------------------------------------------------------------*/
void restart_write_aledyn(ALE_DYNAMIC       *adyn,
                            FIELD	    *actfield,
                            PARTITION	    *actpart,
                            INTRA	    *actintra)
{
#ifndef NO_TEXT_OUTPUT
INT                  i;
INT                  ierr;
INT                  numnp;
long int           **node_handles;
char                 resname[100];
long int             longdummy;
FILE                *out;
RESTART_DYNALE       res;
NODE                *actnode;

#ifdef DEBUG
dstrc_enter("restart_write_aledyn");
#endif
#ifdef PERF
  perf_begin(74);
#endif
/*----------------------------------------------------------------------*/
out = allfiles.out_pss;
/*-------- check the step we are in and create the name "res<step>" */
res.step = adyn->step;
res.time = adyn->time;
sprintf(resname,"ares%d",res.step);

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   now write the data, that is stored in the nodes of this processor's
   partition

   for a ale node two arrays are important:
   actnode->sol
   actnode->sol_increment
   additionally there may be a multifield array
   actnode->sol_mf;

   so we need to store an array of numnp x 3 handles
*/
/*----------------------------------------------------------------------*/
numnp = actpart->pdis[0].numnp;
res.node_handles = (long int**)CCAMALLOC(numnp*sizeof(long int*));
node_handles = res.node_handles;
res.node_handles[0] = (long int*)CCAMALLOC(3*numnp*sizeof(long int));
for (i=1; i<numnp; i++)
node_handles[i] = &(node_handles[0][i*3]);
/*----------------------------------------------------------------------*/
/* now we loop the nodes on the partition and each node writes his ARRAYs */
for (i=0; i<numnp; i++)
{
   actnode = actpart->pdis[0].node[i];
   pss_write_array(&(actnode->sol),&(node_handles[i][0]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   pss_write_array(&(actnode->sol_increment),&(node_handles[i][1]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   if (actnode->sol_mf.Typ == cca_DA)
   {
   pss_write_array(&(actnode->sol_mf),&(node_handles[i][2]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   }
}
/*
   all nodal data is written, so we now write the node_handles and store
   the handle to this in res.handle_of_node_handles
*/
pss_write("nod_hand",numnp,3,sizeof(long int),node_handles[0],&(res.handle_of_node_handles),out,&ierr);
if (ierr != 1) dserror("Error writing restart data");
/*----------------- delete the res.node_handles but keep the dimensions */
res.node_fdim = numnp;
res.node_sdim = 3;
CCAFREE(res.node_handles[0]);
CCAFREE(res.node_handles);

/*----------------------------------------------------------------------*/
/*
   the only thing to do is now write the RESTART_DYNSTRUCT itself with
   its unique name resname = "res<step>"
   NOTE:
   names are limited to 9 characters, so a step larger then ares99999
   can not be restarted at the moment !!!!
*/
pss_write(resname,1,1,sizeof(RESTART_DYNALE),&res,&longdummy,out,&ierr);
if (ierr != 1) dserror("Error writing restart data");
/*----------------------------------------------------------------------*/
/*
#ifdef DEBUG
pss_status_to_err();
#endif
*/
/*----------------------------------------------------------------------*/
/*
   now write a notice to the err file, that this retsrat step was
   successfully written
*/
if (res.step != 0)
{
   fprintf(allfiles.out_err,"===========================================\n");
   fprintf(allfiles.out_err,"RESTART MESSAGE\n");
   fprintf(allfiles.out_err,"In step %d ale restart data were written\n",res.step);
   fprintf(allfiles.out_err,"Calculation can be restarted using \n");
   fprintf(allfiles.out_err,"RESTART      %d \n",res.step);
   fprintf(allfiles.out_err,"in the input file\n");
   fprintf(allfiles.out_err,"===========================================\n");
   fflush(allfiles.out_err);
}
/*----------------------------------------------------------------------*/
#ifdef PERF
  perf_end(74);
#endif
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of restart_write_fluiddyn */

/*!---------------------------------------------------------------------
\brief read restart data of ale field

<pre>                                                         genk 09/03

read restart data of ALE field from pss file

</pre>
 \param   *adyn	     ALE_DYNAMIC      (i)
\param   *actfield   FIELD            (i)   actual field
\param   *actpart    PARTITION        (i)   actual partition
\param   *actintra   INTRA            (i)   actual intra-communicator

\return void

------------------------------------------------------------------------*/
void restart_read_aledyn(INT          restart,
                         ALE_DYNAMIC *adyn,
                         FIELD	     *actfield,
                         PARTITION   *actpart,
                         INTRA	     *actintra)
{
#ifndef NO_TEXT_OUTPUT
INT                  i;
INT                  ierr;
long int             reshandle;
INT                  byte;
INT                  dims[3];
INT                  numnp;
long int           **node_handles;
char                 resname[100];
FILE                *in;
RESTART_DYNALE       res;
NODE                *actnode;

#ifdef PARALLEL
INT                  j;
INT                  fdim;
INT                  sender;
#endif

#ifdef DEBUG
dstrc_enter("restart_read_aledyn");
#endif
/*----------------------------------------------------------------------*/
in = allfiles.in_pss;
/*----------------------------------- check the step that shall be read */
sprintf(resname,"ares%d",restart);
pss_chck(resname,&reshandle,in,&ierr);
if (ierr != 1) dserror("Cannot restart, step doesn't exist in pss-file");
/*----------------------------- the structure res exists, so we read it */
pss_read_name_handle(resname,&i,&i,&byte,&res,&reshandle,in,&ierr);
if (ierr != 1) dserror("Restart structure exists, but cannot read it");
adyn->step=res.step;
adyn->time=res.time;

/*----------------------------------------------------------------------*/
/*
   now read the data that is stored in the nodes of this processors
   partition

   a node has 3 ARRAYS:
   actnode->sol
   actnode->sol_increment
   actnode->sol_residual

   so we need to read an array of numnp x 3 handles first
*/
/*----------------------------------------------------------------------*/
numnp = actpart->pdis[0].numnp;
if (numnp != res.node_fdim || 3 != res.node_sdim)
    dserror("Mismatch in number of nodes on reading restart");
/*----------------------------------------- define the array of handles */
numnp = actpart->pdis[0].numnp;
res.node_handles = (long int**)CCAMALLOC(numnp*sizeof(long int*));
node_handles = res.node_handles;
res.node_handles[0] = (long int*)CCAMALLOC(3*numnp*sizeof(long int));
for (i=1; i<numnp; i++)
node_handles[i] = &(node_handles[0][i*3]);
/*------------------------------------------- read the array of handles */
pss_read_name_handle("nod_hand",&(res.node_fdim),&(res.node_sdim),&i,
                     node_handles[0],&res.handle_of_node_handles,in,&ierr);
if (ierr != 1) dserror("Cannot read restart data");
/*---------------- now we loop the nodes and each node reads his ARRAYs */
for (i=0; i<numnp; i++)
{
   actnode = actpart->pdis[0].node[i];
   /*------------------------- check for the dimensions of actnode->sol */
   pss_getdims_name_handle(actnode->sol.name,&dims[0],&dims[1],&dims[2],&(node_handles[i][0]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if ((unsigned)dims[2] != sizeof(DOUBLE)) dserror("Cannot read restart data");
   /*------------------------ redefine it, if dimension mismatch occurs */
   if (dims[0] != actnode->sol.fdim ||
       dims[1] != actnode->sol.sdim)
   amredef(&(actnode->sol),dims[0],dims[1],"DA");
   /*---------------------------------------------------------- read it */
   pss_read_array_name_handle(actnode->sol.name,&(actnode->sol),&(node_handles[i][0]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   /*-------------------------------------------------------------------*/
   /*--------------- check for the dimensions of actnode->sol_increment */
   pss_getdims_name_handle(actnode->sol_increment.name,&dims[0],&dims[1],&dims[2],&(node_handles[i][1]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if ((unsigned)dims[2] != sizeof(DOUBLE)) dserror("Cannot read restart data");
   /*------------------------ redefine it, if dimension mismatch occurs */
   if (dims[0] != actnode->sol_increment.fdim ||
       dims[1] != actnode->sol_increment.sdim)
   amredef(&(actnode->sol_increment),dims[0],dims[1],"DA");
   /*---------------------------------------------------------- read it */
   pss_read_array_name_handle(actnode->sol_increment.name,&(actnode->sol_increment),&(node_handles[i][1]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   /*-------------------------------------------------------------------*/
   /*---------------- check for the dimensions of actnode->sol_mf */
   if (actnode->sol_mf.Typ == cca_DA)
   {
   pss_getdims_name_handle(actnode->sol_mf.name,&dims[0],&dims[1],&dims[2],&(node_handles[i][2]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if ((unsigned)dims[2] != sizeof(DOUBLE)) dserror("Cannot read restart data");
   /*------------------------ redefine it, if dimension mismatch occurs */
   if (dims[0] != actnode->sol_mf.fdim ||
       dims[1] != actnode->sol_mf.sdim)
   amredef(&(actnode->sol_mf),dims[0],dims[1],"DA");
   /*---------------------------------------------------------- read it */
   pss_read_array_name_handle(actnode->sol_mf.name,&(actnode->sol_mf),&(node_handles[i][2]),in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   }
   /*-------------------------------------------------------------------*/
} /* end of for (i=0; i<numnp; i++) */
/*------------------------------- delete the handle array for the nodes */
/*amdel(&(res.node_handles));*/
CCAFREE(res.node_handles[0]);
CCAFREE(res.node_handles);
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/*
   now we have to make the arrays node->sol, node->sol_increment,
   node->sol_residual redundant for the whole field
*/
#ifdef PARALLEL
numnp = actfield->dis[0].numnp;
for (i=0; i<numnp; i++)
{
   actnode = &(actfield->dis[0].node[i]);
   sender  = actnode->proc;
   /* now we make sol */
   fdim = actnode->sol.fdim;
   MPI_Bcast(&fdim,1,MPI_INT,sender,actintra->MPI_INTRA_COMM);
   if (fdim != actnode->sol.fdim)
   amredef(&(actnode->sol),fdim,actnode->sol.sdim,"DA");
   j = actnode->sol.fdim * actnode->sol.sdim;
   MPI_Bcast(actnode->sol.a.da[0],j,MPI_DOUBLE,sender,actintra->MPI_INTRA_COMM);
   /* now we make sol_increment */
   fdim = actnode->sol_increment.fdim;
   MPI_Bcast(&fdim,1,MPI_INT,sender,actintra->MPI_INTRA_COMM);
   if (fdim != actnode->sol_increment.fdim)
   amredef(&(actnode->sol_increment),fdim,actnode->sol_increment.sdim,"DA");
   j = actnode->sol_increment.fdim * actnode->sol_increment.sdim;
   MPI_Bcast(actnode->sol_increment.a.da[0],j,MPI_DOUBLE,sender,actintra->MPI_INTRA_COMM);
   /* now we make sol_residual */
   if (actnode->sol_mf.Typ == cca_DA)
   {
   fdim = actnode->sol_mf.fdim;
   MPI_Bcast(&fdim,1,MPI_INT,sender,actintra->MPI_INTRA_COMM);
   if (fdim != actnode->sol_mf.fdim)
   amredef(&(actnode->sol_mf),fdim,actnode->sol_mf.sdim,"DA");
   j = actnode->sol_mf.fdim * actnode->sol_mf.sdim;
   MPI_Bcast(actnode->sol_mf.a.da[0],j,MPI_DOUBLE,sender,actintra->MPI_INTRA_COMM);
   }
}
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of restart_read_aledyn */

/*!---------------------------------------------------------------------
\brief write restart data for FSI-problems

<pre>                                                         genk 09/03

write some specific parameters for restart of FSI-problems

</pre>
\param *fsidyn	  FSI_DYNAMIC (i)

\return void

------------------------------------------------------------------------*/
void restart_write_fsidyn(FSI_DYNAMIC       *fsidyn)
{
#ifndef NO_TEXT_OUTPUT
INT                  ierr;
long int             longdummy;
char                 resname[100];
FILE                *out;
RESTART_DYNFSI       res;


#ifdef DEBUG
dstrc_enter("restart_write_fsidyn");
#endif
#ifdef PERF
  perf_begin(74);
#endif
/*----------------------------------------------------------------------*/
out = allfiles.out_pss;
/*-------- check the step we are in and create the name "res<step>" */
res.step  = fsidyn->step;
res.time  = fsidyn->time;
res.relax = fsidyn->relax;
sprintf(resname,"ires%d",res.step);

/*----------------------------------------------------------------------*/
/*
   the only thing to do is now write the RESTART_DYNSTRUCT itself with
   its unique name resname = "res<step>"
   NOTE:
   names are limited to 9 characters, so a step larger then ires99999
   can not be restarted at the moment !!!!
*/
pss_write(resname,1,1,sizeof(RESTART_DYNFSI),&res,&longdummy,out,&ierr);
if (ierr != 1) dserror("Error writing restart data");
/*----------------------------------------------------------------------*/
/*
#ifdef DEBUG
pss_status_to_err();
#endif
*/
/*----------------------------------------------------------------------*/
/*
   now write a notice to the err file, that this retsrat step was
   successfully written
*/
if (res.step != 0)
{
   fprintf(allfiles.out_err,"===========================================\n");
   fprintf(allfiles.out_err,"RESTART MESSAGE\n");
   fprintf(allfiles.out_err,"In step %d fsi restart data were written\n",res.step);
   fprintf(allfiles.out_err,"Calculation can be restarted using \n");
   fprintf(allfiles.out_err,"RESTART      %d \n",res.step);
   fprintf(allfiles.out_err,"in the input file\n");
   fprintf(allfiles.out_err,"===========================================\n");
   fflush(allfiles.out_err);
}
/*----------------------------------------------------------------------*/
#ifdef PERF
  perf_end(74);
#endif
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of restart_write_fsidyn */

/*!---------------------------------------------------------------------
\brief read restart data for FSI problems

<pre>                                                         genk 09/03

read some specific parameters for restart of FSI-problems

</pre>
\param *fsidyn	 FSI_DYNAMIC (i)

\return void

------------------------------------------------------------------------*/
void restart_read_fsidyn(INT restart,
                         FSI_DYNAMIC   *fsidyn)
{
#ifndef NO_TEXT_OUTPUT
INT                  i,ierr;
long int             reshandle;
INT                  byte;
char                 resname[100];
FILE                *in;
RESTART_DYNFSI       res;

#ifdef DEBUG
dstrc_enter("restart_read_fsidyn");
#endif
/*----------------------------------------------------------------------*/
in = allfiles.in_pss;
/*----------------------------------- check the step that shall be read */
sprintf(resname,"ires%d",restart);
pss_chck(resname,&reshandle,in,&ierr);
if (ierr != 1) dserror("Cannot restart, step doesn't exist in pss-file");
/*----------------------------- the structure res exists, so we read it */
pss_read_name_handle(resname,&i,&i,&byte,&res,&reshandle,in,&ierr);
if (ierr != 1) dserror("Restart structure exists, but cannot read it");
fsidyn->step=res.step;
fsidyn->time=res.time;
fsidyn->relax=res.relax;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of restart_read_fsidyn */

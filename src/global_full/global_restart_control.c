#include "../headers/standardtypes.h"
#include "../headers/solution.h"
#include "../headers/restart.h"
#include "../shell8/shell8.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | struct _GENPROB       genprob; defined in global_control.c           |
 *----------------------------------------------------------------------*/
extern struct _GENPROB       genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 | write restart                                         m.gee 05/02    |
 | of nln-structural dynamics                                           |
 | NOTE:                                                                |
 | restart can only be performed with the last written restart step     |
 | noticed in the err-file                                              |
 | restarting in an earlier step works, but it leads to                 |
 | duplicate restart entries in the pss-file. This is o.k. as long as   |
 | the run goes further in number of steps as the one before            |
 | and writes a restart step which is outside the old step range        |
 | (bissle raetselhaft oder?)                                           |
 | Klartext:                                                            |
 | Wenn man nicht mit dem letzten restart-Schritt restartet koennen     |
 | lustige Effekte auftreten.                                           |
 *----------------------------------------------------------------------*/
void restart_write_nlnstructdyn(STRUCT_DYNAMIC  *sdyn,
                                STRUCT_DYN_CALC *dynvar,
                                FIELD           *actfield,
                                PARTITION       *actpart,
                                INTRA           *actintra,
                                CALC_ACTION     *action,
                                int nrhs,  DIST_VECTOR *rhs,
                                int nsol,  DIST_VECTOR *sol,
                                int ndis,  DIST_VECTOR *dispi,
                                int nvel,  DIST_VECTOR *vel,
                                int nacc,  DIST_VECTOR *acc,
                                int nfie,  DIST_VECTOR *fie,
                                int nwork, DIST_VECTOR *work,
                                ARRAY *intforce_a,
                                ARRAY *dirich_a)
{
int                  i;
int                  ierr;
int                  numnp;
int                **node_handles;
int                  numele;
int                **ele_handles;
char                 resname[100];
RESTART_DYNSTRUCT    res;
DIST_VECTOR         *distwrite;
NODE                *actnode;
ELEMENT             *actele;
#ifdef DEBUG 
dstrc_enter("restart_write_nlnstructdyn");
#endif
/*----------------------------------------------------------------------*/
/*-------- check the step we are in and create the name "res<step>" */
res.step = sdyn->step;
sprintf(resname,"res%d",res.step);
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
    pss_write_array(&(distwrite->vec),&(res.dist_vec_rhs[i+1]),&ierr);
    if (ierr != 1) dserror("Error writing restart data");
 }
/*------------------------------------------------------ write the *sol */
 for (i=0; i<nsol; i++)
 {
    distwrite = &(sol[i]);
    pss_write_array(&(distwrite->vec),&(res.dist_vec_sol[i+1]),&ierr);
    if (ierr != 1) dserror("Error writing restart data");
 }
/*------------------------------------------------------ write the *dispi */
 for (i=0; i<ndis; i++)
 {
    distwrite = &(dispi[i]);
    pss_write_array(&(distwrite->vec),&(res.dist_vec_dispi[i+1]),&ierr);
    if (ierr != 1) dserror("Error writing restart data");
 }
/*------------------------------------------------------ write the *vel */
 for (i=0; i<nvel; i++)
 {
    distwrite = &(vel[i]);
    pss_write_array(&(distwrite->vec),&(res.dist_vec_vel[i+1]),&ierr);
    if (ierr != 1) dserror("Error writing restart data");
 }
/*------------------------------------------------------ write the *acc */
 for (i=0; i<nacc; i++)
 {
    distwrite = &(acc[i]);
    pss_write_array(&(distwrite->vec),&(res.dist_vec_acc[i+1]),&ierr);
    if (ierr != 1) dserror("Error writing restart data");
 }
/*------------------------------------------------------ write the *fie */
 for (i=0; i<nfie; i++)
 {
    distwrite = &(fie[i]);
    pss_write_array(&(distwrite->vec),&(res.dist_vec_fie[i+1]),&ierr);
    if (ierr != 1) dserror("Error writing restart data");
 }
/*------------------------------------------------------ write the *work */
 for (i=0; i<nwork; i++)
 {
    distwrite = &(work[i]);
    pss_write_array(&(distwrite->vec),&(res.dist_vec_work[i+1]),&ierr);
    if (ierr != 1) dserror("Error writing restart data");
 }
/*----------------------------------------- write the ARRAY *intforce_a */
pss_write_array(intforce_a,&(res.intforce),&ierr);
/*------------------------------------------- write the ARRAY *dirich_a */
pss_write_array(dirich_a,&(res.dirich),&ierr);
/*----------------------------------------------------------------------*/
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
node_handles = amdef("nodehand",&(res.node_handles),numnp,3,"IA");
/*----------------------------------------------------------------------*/
/* now we loop the nodes on the partition and each node writes his ARRAYs */
for (i=0; i<numnp; i++)
{
   actnode = actpart->pdis[0].node[i];
   pss_write_array(&(actnode->sol),&(node_handles[i][0]),&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   pss_write_array(&(actnode->sol_increment),&(node_handles[i][1]),&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   pss_write_array(&(actnode->sol_residual),&(node_handles[i][2]),&ierr);
   if (ierr != 1) dserror("Error writing restart data");
}
/* 
   all nodal data is written, so we now write the node_handles and store
   the handle to this in res.handle_of_node_handles
*/
pss_write_array(&(res.node_handles),&(res.handle_of_node_handles),&ierr);   
if (ierr != 1) dserror("Error writing restart data");
/*----------------- delete the res.node_handles but keep the dimensions */
amdel(&(res.node_handles));
res.node_handles.fdim = numnp;
res.node_handles.sdim = 3;
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/* 
   we come to the most difficult part, the writing of the element data.
   One element is assumed not to have more then 5 different records to
   write of arbitary type
   
   The element is called with the element action action = write_restart
   Again, all elements to the processors partition are written
*/
/*----------------------------------------------------------------------*/
numele = actpart->pdis[0].numele;
ele_handles = amdef("elehandl",&(res.ele_handles),numele,5,"IA");
/*--------------------- now loop element and switch for type of element */
*action = write_restart;
for (i=0; i<actpart->pdis[0].numele; i++)
{
   actele = actpart->pdis[0].element[i];
   switch(actele->eltyp)/*======================= call element routines */
   {
   case el_shell8:
      shell8(actfield,actpart,actintra,actele,
             NULL,NULL,NULL,
             0,res.ele_handles.sdim,ele_handles[i],action);
   break;
   case el_brick1:
       dserror("Restart for brick not yet impl.");
   break;
   case el_wall1:
       dserror("Restart for wall not yet impl.");
   break;
   case el_fluid2: 
       dserror("Restart for fluid2 not yet impl.");
   break;
   case el_fluid3: 
       dserror("Restart for fluid3 not yet impl.");
   break;
   case el_ale:
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
/* 
   all ele data is written, so write the ele_handles and store the handle to
   it in res.handle_of_ele_handles
*/   
pss_write_array(&(res.ele_handles),&(res.handle_of_ele_handles),&ierr);   
if (ierr != 1) dserror("Error writing restart data");
/*------------------ delete the res.ele_handles but keep the dimensions */
amdel(&(res.ele_handles));
res.ele_handles.fdim = numele;
res.ele_handles.sdim = 5;
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   the only thing to do is now write the RESTART_DYNSTRUCT itself with 
   its unique name resname = "res<step>"
   NOTE:
   names are limited to 9 characters, so a step larger then res999999
   can not be restarted at the moment !!!!
*/   
pss_write(resname,1,1,sizeof(RESTART_DYNSTRUCT),&res,&i,&ierr);
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
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of restart_write_nlnstructdyn */




/*----------------------------------------------------------------------*
 | read restart                                          m.gee 05/02    |
 | of nln-structural dynamics                                           |
 *----------------------------------------------------------------------*/
void restart_read_nlnstructdyn(int restart,
                               STRUCT_DYNAMIC  *sdyn,
                               STRUCT_DYN_CALC *dynvar,
                               FIELD           *actfield,
                               PARTITION       *actpart,
                               INTRA           *actintra,
                               CALC_ACTION     *action,
                               int nrhs,  DIST_VECTOR *rhs,
                               int nsol,  DIST_VECTOR *sol,
                               int ndis,  DIST_VECTOR *dispi,
                               int nvel,  DIST_VECTOR *vel,
                               int nacc,  DIST_VECTOR *acc,
                               int nfie,  DIST_VECTOR *fie,
                               int nwork, DIST_VECTOR *work,
                               ARRAY *intforce_a,
                               ARRAY *dirich_a)
{
int                  i;
int                  ierr;
int                  reshandle;
int                  byte;
int                  dims[3];
int                  numnp;
int                **node_handles;
int                  numele;
int                **ele_handles;
char                 resname[100];
RESTART_DYNSTRUCT    res;
DIST_VECTOR         *distread;
NODE                *actnode;
ELEMENT             *actele;
#ifdef DEBUG 
dstrc_enter("restart_read_nlnstructdyn");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------- check the step that shall be read */
sprintf(resname,"res%d",restart);
pss_chck(resname,&reshandle,&ierr);
if (ierr != 1) dserror("Cannot restart, step doesn't exist in pss-file");
/*----------------------------- the structure res exists, so we read it */
pss_read_name_handle(resname,&i,&i,&byte,&res,&reshandle,&ierr);
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
   pss_read_array_name_handle(distread->vec.name,&(distread->vec),&(res.dist_vec_rhs[i+1]),&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
}
/*------------------------------------------------------------ read sol */
for (i=0; i<nsol; i++)
{
   distread = &(sol[i]);
   pss_read_array_name_handle(distread->vec.name,&(distread->vec),&(res.dist_vec_sol[i+1]),&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
}
/*------------------------------------------------------------ read dispi */
for (i=0; i<ndis; i++)
{
   distread = &(dispi[i]);
   pss_read_array_name_handle(distread->vec.name,&(distread->vec),&(res.dist_vec_dispi[i+1]),&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
}
/*------------------------------------------------------------ read vel */
for (i=0; i<nvel; i++)
{
   distread = &(vel[i]);
   pss_read_array_name_handle(distread->vec.name,&(distread->vec),&(res.dist_vec_vel[i+1]),&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
}
/*------------------------------------------------------------ read acc */
for (i=0; i<nacc; i++)
{
   distread = &(acc[i]);
   pss_read_array_name_handle(distread->vec.name,&(distread->vec),&(res.dist_vec_acc[i+1]),&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
}
/*------------------------------------------------------------ read fie */
for (i=0; i<nfie; i++)
{
   distread = &(fie[i]);
   pss_read_array_name_handle(distread->vec.name,&(distread->vec),&(res.dist_vec_fie[i+1]),&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
}
/*------------------------------------------------------------ read work */
for (i=0; i<nwork; i++)
{
   distread = &(work[i]);
   pss_read_array_name_handle(distread->vec.name,&(distread->vec),&(res.dist_vec_work[i+1]),&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
}
/*------------------------------------------- read the ARRAY *intforce_a */
pss_read_array_name_handle(intforce_a->name,intforce_a,&(res.intforce),&ierr);
if (ierr != 1) dserror("Cannot read restart data");
/*--------------------------------------------- read the ARRAY *dirich_a */
pss_read_array_name_handle(dirich_a->name,dirich_a,&(res.dirich),&ierr);
if (ierr != 1) dserror("Cannot read restart data");
/*----------------------------------------------------------------------*/
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
if (numnp != res.node_handles.fdim || 3 != res.node_handles.sdim)
    dserror("Mismatch in number of nodes on reading restart");
/*----------------------------------------- define the array of handles */
node_handles = amdef("nodehand",&(res.node_handles),numnp,3,"IA");
/*------------------------------------------- read the array of handles */
pss_read_array_name_handle(res.node_handles.name,&(res.node_handles),&(res.handle_of_node_handles),&ierr);
if (ierr != 1) dserror("Cannot read restart data");
/*---------------- now we loop the nodes and each node reads his ARRAYs */
for (i=0; i<numnp; i++)
{
   actnode = actpart->pdis[0].node[i];
   /*------------------------- check for the dimensions of actnode->sol */
   pss_getdims_name_handle(actnode->sol.name,&dims[0],&dims[1],&dims[2],&(node_handles[i][0]),&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (dims[2] != sizeof(double)) dserror("Cannot read restart data");
   /*------------------------ redefine it, if dimension mismatch occurs */
   if (dims[0] != actnode->sol.fdim ||
       dims[1] != actnode->sol.sdim)
   amredef(&(actnode->sol),dims[0],dims[1],"DA");
   /*---------------------------------------------------------- read it */
   pss_read_array_name_handle(actnode->sol.name,&(actnode->sol),&(node_handles[i][0]),&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   /*-------------------------------------------------------------------*/
   /*--------------- check for the dimensions of actnode->sol_increment */
   pss_getdims_name_handle(actnode->sol_increment.name,&dims[0],&dims[1],&dims[2],&(node_handles[i][1]),&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (dims[2] != sizeof(double)) dserror("Cannot read restart data");
   /*------------------------ redefine it, if dimension mismatch occurs */
   if (dims[0] != actnode->sol_increment.fdim ||
       dims[1] != actnode->sol_increment.sdim)
   amredef(&(actnode->sol_increment),dims[0],dims[1],"DA");
   /*---------------------------------------------------------- read it */
   pss_read_array_name_handle(actnode->sol_increment.name,&(actnode->sol_increment),&(node_handles[i][1]),&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   /*-------------------------------------------------------------------*/
   /*---------------- check for the dimensions of actnode->sol_residual */
   pss_getdims_name_handle(actnode->sol_residual.name,&dims[0],&dims[1],&dims[2],&(node_handles[i][2]),&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (dims[2] != sizeof(double)) dserror("Cannot read restart data");
   /*------------------------ redefine it, if dimension mismatch occurs */
   if (dims[0] != actnode->sol_residual.fdim ||
       dims[1] != actnode->sol_residual.sdim)
   amredef(&(actnode->sol_residual),dims[0],dims[1],"DA");
   /*---------------------------------------------------------- read it */
   pss_read_array_name_handle(actnode->sol_residual.name,&(actnode->sol_residual),&(node_handles[i][2]),&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   /*-------------------------------------------------------------------*/
} /* end of for (i=0; i<numnp; i++) */
/*------------------------------- delete the handle array for the nodes */
amdel(&(res.node_handles));
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   now start reading the element data 
*/
/*----------------------------------------------------------------------*/
numele = actpart->pdis[0].numele;
if (numele != res.ele_handles.fdim || 5 != res.ele_handles.sdim)
    dserror("Mismatch in number of elements on reading restart");
/*----------------------------------------- define the array of handles */
ele_handles = amdef("elehandl",&(res.ele_handles),numele,5,"IA");
/*------------------------------------------- read the array of handles */
pss_read_array_name_handle(res.ele_handles.name,&(res.ele_handles),&(res.handle_of_ele_handles),&ierr);
if (ierr != 1) dserror("Cannot read restart data");
/*--------------------- now loop element and switch for type of element */
*action = read_restart;
for (i=0; i<actpart->pdis[0].numele; i++)
{
   actele = actpart->pdis[0].element[i];
   switch(actele->eltyp)/*======================= call element routines */
   {
   case el_shell8:
      shell8(actfield,actpart,actintra,actele,
             NULL,NULL,NULL,
             0,res.ele_handles.sdim,ele_handles[i],action);
   break;
   case el_brick1:
       dserror("Restart for brick not yet impl.");
   break;
   case el_wall1:
       dserror("Restart for wall not yet impl.");
   break;
   case el_fluid2: 
       dserror("Restart for fluid2 not yet impl.");
   break;
   case el_fluid3: 
       dserror("Restart for fluid3 not yet impl.");
   break;
   case el_ale:
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
amdel(&(res.ele_handles));
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of restart_read_nlnstructdyn */

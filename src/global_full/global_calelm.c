#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "../shell8/shell8.h"
#include "../wall1/wall1.h"
#include "../brick1/brick1.h"
#include "../fluid3/fluid3.h"
#include "../ale3/ale3.h"
#include "../ale2/ale2.h"
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
enum _CALC_ACTION calc_action[MAXFIELD];
/*----------------------------------------------------------------------*
 | global dense matrices for element routines             m.gee 7/01    |
 *----------------------------------------------------------------------*/
struct _ARRAY estif_global;    /* element stiffness matrix              */
struct _ARRAY emass_global;    /* element mass matrix                   */  
struct _ARRAY etforce_global;  /* element Time RHS                      */
struct _ARRAY eiforce_global;  /* element Iteration RHS                 */
struct _ARRAY edforce_global;  /* element dirichlet RHS                 */
struct _ARRAY intforce_global;
/*----------------------------------------------------------------------*
 |  routine to call elements                             m.gee 6/01     |
 *----------------------------------------------------------------------*/
void calelm(FIELD        *actfield,     /* active field */        
            SOLVAR       *actsolv,      /* active SOLVAR */
            PARTITION    *actpart,      /* my partition of this field */
            INTRA        *actintra,     /* my intra-communicator */
            int           sysarray1,    /* number of first sparse system matrix */
            int           sysarray2,    /* number of secnd system matrix, if present, else -1 */
            CONTAINER    *container,    /*!< contains variables defined in container.h */
            CALC_ACTION  *action)       /* calculation option passed to element routines */     
/*----------------------------------------------------------------------*/
{
int               i;
int               hasdirich=0;      /* flag                             */
int               hasext=0;         /* flag                             */
ELEMENT          *actele;
SPARSE_TYP        sysarray1_typ;
SPARSE_TYP        sysarray2_typ;
ASSEMBLE_ACTION   assemble_action;

#ifdef DEBUG 
dstrc_enter("calelm");
#endif
/*----------------------------------------------------------------------*/
/*-------------- zero the parallel coupling exchange buffers if present */  
#ifdef PARALLEL 
/*------------------------ check the send & recv buffers from sysarray1 */
if (sysarray1 != -1)
{
   switch(actsolv->sysarray_typ[sysarray1])
   {
   case msr:
      if (actsolv->sysarray[sysarray1].msr->couple_d_send)
         amzero(actsolv->sysarray[sysarray1].msr->couple_d_send);
      if (actsolv->sysarray[sysarray1].msr->couple_d_recv)
         amzero(actsolv->sysarray[sysarray1].msr->couple_d_recv);
   break;
   case parcsr:
      if (actsolv->sysarray[sysarray1].parcsr->couple_d_send)
         amzero(actsolv->sysarray[sysarray1].parcsr->couple_d_send);
      if (actsolv->sysarray[sysarray1].parcsr->couple_d_recv)
         amzero(actsolv->sysarray[sysarray1].parcsr->couple_d_recv);
   break;
   case ucchb:
      if (actsolv->sysarray[sysarray1].ucchb->couple_d_send)
         amzero(actsolv->sysarray[sysarray1].ucchb->couple_d_send);
      if (actsolv->sysarray[sysarray1].ucchb->couple_d_recv)
         amzero(actsolv->sysarray[sysarray1].ucchb->couple_d_recv);
   break;
   case dense:
      if (actsolv->sysarray[sysarray1].dense->couple_d_send)
         amzero(actsolv->sysarray[sysarray1].dense->couple_d_send);
      if (actsolv->sysarray[sysarray1].dense->couple_d_recv)
         amzero(actsolv->sysarray[sysarray1].dense->couple_d_recv);
   break;
   case rc_ptr:
      if (actsolv->sysarray[sysarray1].rc_ptr->couple_d_send)
         amzero(actsolv->sysarray[sysarray1].rc_ptr->couple_d_send);
      if (actsolv->sysarray[sysarray1].rc_ptr->couple_d_recv)
         amzero(actsolv->sysarray[sysarray1].rc_ptr->couple_d_recv);
   break;
   case ccf:
      if (actsolv->sysarray[sysarray1].ccf->couple_d_send)
         amzero(actsolv->sysarray[sysarray1].ccf->couple_d_send);
      if (actsolv->sysarray[sysarray1].ccf->couple_d_recv)
         amzero(actsolv->sysarray[sysarray1].ccf->couple_d_recv);
   break;
   case skymatrix:
      if (actsolv->sysarray[sysarray1].sky->couple_d_send)
         amzero(actsolv->sysarray[sysarray1].sky->couple_d_send);
      if (actsolv->sysarray[sysarray1].sky->couple_d_recv)
         amzero(actsolv->sysarray[sysarray1].sky->couple_d_recv);
   break;
   case spoolmatrix:
      if (actsolv->sysarray[sysarray1].spo->couple_d_send)
         amzero(actsolv->sysarray[sysarray1].spo->couple_d_send);
      if (actsolv->sysarray[sysarray1].spo->couple_d_recv)
         amzero(actsolv->sysarray[sysarray1].spo->couple_d_recv);
   break;
   case bdcsr:;
   break;
   default:
      dserror("Unknown typ of system matrix");
   break;
   }
}
/*------------------------ check the send & recv buffers from sysarray2 */
if (sysarray2 != -1)
{
   switch(actsolv->sysarray_typ[sysarray2])
   {
   case msr:
      if (actsolv->sysarray[sysarray2].msr->couple_d_send)
         amzero(actsolv->sysarray[sysarray2].msr->couple_d_send);
      if (actsolv->sysarray[sysarray2].msr->couple_d_recv)
         amzero(actsolv->sysarray[sysarray2].msr->couple_d_send);
   break;
   case parcsr:
      if (actsolv->sysarray[sysarray2].parcsr->couple_d_send)
         amzero(actsolv->sysarray[sysarray2].parcsr->couple_d_send);
      if (actsolv->sysarray[sysarray2].parcsr->couple_d_recv)
         amzero(actsolv->sysarray[sysarray2].parcsr->couple_d_send);
   break;
   case ucchb:
      if (actsolv->sysarray[sysarray2].ucchb->couple_d_send)
         amzero(actsolv->sysarray[sysarray2].ucchb->couple_d_send);
      if (actsolv->sysarray[sysarray2].ucchb->couple_d_recv)
         amzero(actsolv->sysarray[sysarray2].ucchb->couple_d_send);
   break;
   case dense:
      if (actsolv->sysarray[sysarray2].dense->couple_d_send)
         amzero(actsolv->sysarray[sysarray2].dense->couple_d_send);
      if (actsolv->sysarray[sysarray2].dense->couple_d_recv)
         amzero(actsolv->sysarray[sysarray2].dense->couple_d_send);
   break;
   case rc_ptr:
      if (actsolv->sysarray[sysarray2].rc_ptr->couple_d_send)
         amzero(actsolv->sysarray[sysarray2].rc_ptr->couple_d_send);
      if (actsolv->sysarray[sysarray2].rc_ptr->couple_d_recv)
         amzero(actsolv->sysarray[sysarray2].rc_ptr->couple_d_recv);
   break;
   case ccf:
      if (actsolv->sysarray[sysarray2].ccf->couple_d_send)
         amzero(actsolv->sysarray[sysarray2].ccf->couple_d_send);
      if (actsolv->sysarray[sysarray2].ccf->couple_d_recv)
         amzero(actsolv->sysarray[sysarray2].ccf->couple_d_recv);
   break;
   case skymatrix:
      if (actsolv->sysarray[sysarray2].sky->couple_d_send)
         amzero(actsolv->sysarray[sysarray2].sky->couple_d_send);
      if (actsolv->sysarray[sysarray2].sky->couple_d_recv)
         amzero(actsolv->sysarray[sysarray2].sky->couple_d_recv);
   break;
   case spoolmatrix:
      if (actsolv->sysarray[sysarray2].spo->couple_d_send)
         amzero(actsolv->sysarray[sysarray2].spo->couple_d_send);
      if (actsolv->sysarray[sysarray2].spo->couple_d_recv)
         amzero(actsolv->sysarray[sysarray2].spo->couple_d_recv);
   break;
   case bdcsr:;
   break;
   default:
      dserror("Unknown typ of system matrix");
   break;
   }
}
#endif
/* =======================================================call elements */
/*---------------------------------------------- loop over all elements */
for (i=0; i<actpart->pdis[0].numele; i++)
{
   /*------------------------------------ set pointer to active element */
   actele = actpart->pdis[0].element[i];
   /* if present, init the element vectors intforce_global and dirich_global */
   if (container->dvec) amzero(&intforce_global);
   switch(actele->eltyp)/*======================= call element routines */
   {
   case el_shell8:
      container->handsize = 0;
      container->handles  = NULL;
      shell8(actfield,actpart,actintra,actele,
             &estif_global,&emass_global,&intforce_global,
             action,container);
   break;
   case el_brick1:
      brick1(actpart,actintra,actele,
             &estif_global,&emass_global,&intforce_global,
             action,container);
   break;
   case el_wall1:
      container->handsize = 0;
      container->handles  = NULL;
      wall1(actpart,actintra,actele,
            &estif_global,&emass_global,&intforce_global,
            action, container);
   break;
   case el_fluid2: 
      fluid2(actpart,actintra,actele,
             &estif_global,&emass_global,
             &etforce_global,&eiforce_global,&edforce_global,
	     action,&hasdirich,&hasext,container);
   break;
   case el_fluid3: 
      fluid3(actpart,actintra,actele,
             &estif_global,&emass_global,
             &etforce_global,&eiforce_global,&edforce_global,
	     action,&hasdirich,&hasext,container); 
   break;
   case el_ale3:
	ale3(actpart,actintra,actele,
        &estif_global,
        action,container);
   break;
   case el_ale2:
	ale2(actpart,actintra,actele,
        &estif_global,
        action,container);
   break;
   case el_none:
      dserror("Typ of element unknown");
   break;
   default:
      dserror("Typ of element unknown");
   }/* end of calling elements */


   switch(*action)/*=== call assembly dependent on calculation-flag */
   {
   case calc_struct_linstiff     : assemble_action = assemble_one_matrix; break;
   case calc_struct_nlnstiff     : assemble_action = assemble_one_matrix; break;
   case calc_struct_nlnstiffmass : assemble_action = assemble_two_matrix; break;
   case calc_struct_internalforce: assemble_action = assemble_do_nothing; break;
   case calc_struct_eleload      : assemble_action = assemble_do_nothing; break;
   case calc_struct_stress       : assemble_action = assemble_do_nothing; break;
   case calc_struct_update_istep : assemble_action = assemble_do_nothing; break;
   case calc_ale_stiff           : assemble_action = assemble_one_matrix; break;
   case calc_ale_rhs             : assemble_action = assemble_do_nothing; break;
   case calc_fluid               : assemble_action = assemble_one_matrix; break;
   default: dserror("Unknown type of assembly"); break;
   }
   /*--------------------------- assemble one or two system matrices */
   assemble(sysarray1,
            &estif_global,
            sysarray2,
            &emass_global,
            actpart,
            actsolv,
            actintra,
            actele,
            assemble_action,
            container);
   /*---------------------------- assemble the vector intforce_global */
   switch(container->fieldtyp)
   {
   case structure:
      if (container->dvec)
      assemble_intforce(actele,&intforce_global,container);
   
      /*------ assemble the rhs vector of condensed dirichlet conditions */
      if (container->dirich && container->isdyn==0)
      assemble_dirich(actele,&estif_global,container);
      if (container->dirich && container->isdyn==1)
      assemble_dirich_dyn(actele,&estif_global,&emass_global,container);
   break;
   case fluid:
      if (container->nif!=0)
      {
         container->dvec = container->ftimerhs;
         assemble_intforce(actele,&etforce_global,container);
      }
   /*-------------- assemble the vector eiforce_global to iteration rhs */
   if (container->nii+hasext!=0)
      {   
         container->dvec = container->fiterhs;
         assemble_intforce(actele,&eiforce_global,container); 
      }
   /*-------------- assemble the vector edforce_global to iteration rhs */
   if (hasdirich!=0)
      {
         container->dvec = container->fiterhs;
         assemble_intforce(actele,&edforce_global,container);
      }   
   break;
   case ale:
   break;
   default:
      dserror("fieldtyp unknown!");
   }
}/* end of loop over elements */
/*----------------------------------------------------------------------*/
/*                    in parallel coupled dofs have to be exchanged now */
/*             (if there are any inter-proc couplings, which is tested) */
/*----------------------------------------------------------------------*/
#ifdef PARALLEL 
switch(*action)
{
case calc_struct_linstiff      : assemble_action = assemble_one_exchange; break;
case calc_struct_nlnstiff      : assemble_action = assemble_one_exchange; break;
case calc_struct_internalforce : assemble_action = assemble_do_nothing; break;
case calc_struct_nlnstiffmass  : assemble_action = assemble_two_exchange; break;
case calc_struct_eleload       : assemble_action = assemble_do_nothing; break;
case calc_struct_stress        : assemble_action = assemble_do_nothing; break;
case calc_struct_update_istep  : assemble_action = assemble_do_nothing; break;
case calc_ale_stiff            : assemble_action = assemble_one_exchange; break;
case calc_ale_rhs              : assemble_action = assemble_do_nothing; break;
case calc_fluid                : assemble_action = assemble_one_exchange; break;
default: dserror("Unknown type of assembly"); break;
}
/*------------------------------ exchange coupled dofs, if there are any */
assemble(sysarray1,
         NULL,
         sysarray2,
         NULL,
         actpart,
         actsolv,
         actintra,
         actele,
         assemble_action,
         container);
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of calelm */





/*----------------------------------------------------------------------*
 |  routine to call elements to init                     m.gee 7/01     |
 *----------------------------------------------------------------------*/
void calinit(FIELD       *actfield,   /* the active physical field */ 
             PARTITION   *actpart,    /* my partition of this field */
             CALC_ACTION *action,
             CONTAINER   *container)  /*!< contains variables defined in container.h */
{
int i;                        /* a counter */
int is_shell8=0;              /* flags to check for presents of certain element types */
int is_brick1=0;
int is_wall1 =0;
int is_fluid2=0;
int is_fluid3=0;
int is_ale3=0;
int is_ale2=0;

ELEMENT *actele;              /* active element */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("calinit");
#endif
/*-------------------------- define dense element matrices for assembly */
amdef("estif",&estif_global,(MAXNOD*MAXDOFPERNODE),(MAXNOD*MAXDOFPERNODE),"DA");
amdef("emass",&emass_global,(MAXNOD*MAXDOFPERNODE),(MAXNOD*MAXDOFPERNODE),"DA");
amdef("etforce",&etforce_global,(MAXNOD*MAXDOFPERNODE),1,"DV");
amdef("eiforce",&eiforce_global,(MAXNOD*MAXDOFPERNODE),1,"DV");
amdef("edforce",&edforce_global,(MAXNOD*MAXDOFPERNODE),1,"DV");
amdef("inforce",&intforce_global,(MAXNOD*MAXDOFPERNODE),1,"DV");
/*--------------------what kind of elements are there in this example ? */
for (i=0; i<actfield->dis[0].numele; i++)
{
   actele = &(actfield->dis[0].element[i]);
   switch(actele->eltyp)
   {
   case el_shell8:
      is_shell8=1;
   break;
   case el_brick1:
      is_brick1=1;
   break;
   case el_wall1:
      is_wall1=1;
   break;
   case el_fluid2:
      is_fluid2=1;
   break;
   case el_fluid3:
      is_fluid3=1;
   break;
   case el_ale3:
      is_ale3=1;
   break;
   case el_ale2:
      is_ale2=1;
   break;
   default:
      dserror("Unknown typ of element");
   break;   
   }
}/* end of loop over all elements */
/*--------------------- init the element routines for all present types */
container->kstep = 0;  
/*------------------------------- init all kind of routines for shell8  */
if (is_shell8==1)
{
   container->handsize = 0;
   container->handles  = NULL;
   shell8(actfield,actpart,NULL,NULL,&estif_global,&emass_global,&intforce_global,
          action,container);
}
/*-------------------------------- init all kind of routines for brick1 */
if (is_brick1==1)
{
   brick1(actpart,NULL,NULL,&estif_global,&emass_global,NULL,action,container);
}
/*-------------------------------- init all kind of routines for wall1  */
if (is_wall1==1)
{
   container->handsize = 0;
   container->handles  = NULL;
   wall1(actpart,NULL,NULL,&estif_global,&emass_global,&intforce_global,
         action,container);
}
/*-------------------------------- init all kind of routines for fluid2 */
if (is_fluid2==1)
{
   fluid2(actpart,NULL,NULL,
          &estif_global,&emass_global,
          &etforce_global,&eiforce_global,&edforce_global,
          action,NULL,NULL,container);
}
/*-------------------------------- init all kind of routines for fluid3 */
if (is_fluid3==1)
{
   fluid3(actpart,NULL,NULL,
          &estif_global,&emass_global,
          &etforce_global,&eiforce_global,&edforce_global,
          action,NULL,NULL,container);
}
/*----------------------------------- init all kind of routines for ale */
if (is_ale3==1)
{
   ale3(actpart,NULL,NULL,&estif_global,action,container);
}
/*----------------------------------- init all kind of routines for ale */
if (is_ale2==1)
{
   ale2(actpart,NULL,NULL,&estif_global,action,container);
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of calinit */







/*----------------------------------------------------------------------*
 |  in here the element's results are made redundant     m.gee 12/01    |
 *----------------------------------------------------------------------*/
void calreduce(FIELD       *actfield, /* the active field */
               PARTITION   *actpart,  /* my partition of this field */
               INTRA       *actintra, /* the field's intra-communicator */
               CALC_ACTION *action,   /* action for element routines */
               CONTAINER   *container)/* contains variables defined in container.h */
{
int i;
int is_shell8=0;
int is_brick1=0;
int is_wall1 =0;
int is_fluid1=0;
int is_fluid3=0;
int is_ale3=0;
ELEMENT *actele;
#ifdef DEBUG 
dstrc_enter("calreduce");
#endif
/*----------------------------------------------------------------------*/
/*--------------------what kind of elements are there in this example ? */
for (i=0; i<actfield->dis[0].numele; i++)
{
   actele = &(actfield->dis[0].element[i]);
   switch(actele->eltyp)
   {
   case el_shell8:
      is_shell8=1;
   break;
   case el_brick1:
      is_brick1=1;
   break;
   case el_wall1:
      is_wall1=1;
   break;
   case el_fluid2:
      is_fluid1=1;
   break;
   case el_fluid3:
      is_fluid3=1;
   break;
   case el_ale3:
      is_ale3=1;
   break;
   default:
      dserror("Unknown typ of element");
   break;   
   }
}/* end of loop over all elements */
/*-------------------------------------------reduce results for shell8  */
if (is_shell8==1)
{
   container->handsize = 0;
   container->handles  = NULL;
   shell8(actfield,actpart,actintra,NULL,NULL,NULL,NULL,action,container);
}
/*--------------------------------------------reduce results for brick1 */
if (is_brick1==1)
{
}
/*---------------------------------------------reduce results for wall1 */
if (is_wall1==1)
{
}
/*--------------------------------------------reduce results for fluid1 */
if (is_fluid1==1)
{
}
/*--------------------------------------------reduce results for fluid3 */
if (is_fluid3==1)
{
}
/*-----------------------------------------------reduce results for ale */
if (is_ale3==1)
{
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of calreduce */

#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 | global dense matrices for element routines             m.gee 7/01    |
 *----------------------------------------------------------------------*/
struct _ARRAY estif_global;
struct _ARRAY emass_global;
/*----------------------------------------------------------------------*
 |  routine to call elements                             m.gee 6/01     |
 *----------------------------------------------------------------------*/
void calelm(FIELD      *actfield,     /* active field */        
            SOLVAR     *actsolv,      /* active SOLVAR */
            PARTITION  *actpart,      /* my partition of this field */
            INTRA      *actintra,     /* my intra-communicator */
            int         sysarray1,    /* number of first sparse system matrix */
            int         sysarray2,    /* number of secnd system matrix, if present, else -1 */
            double     *dvec,         /* global redundant vector passed to elements */
            int         global_numeq, /* size of dvec */
            int         kstep,        /* time in increment step we are in */
            int         calc_option)  /* calculation option passed to element routines */
/*----------------------------------------------------------------------*/
/*
  calc_options which are already in use:
    
    calc_option=0 init the element routines, do no calculation
    calc_option=1 calc linear structural stiffness matrix 
    calc_option=2 calc nonlinear structural stiffness matrix and internal forces
    calc_option=5 calc structural forces and stresses
    calc_option=6 calc load vector of element loads
    calc_option=7 Allreduce stress results to make results redundant on all procs
    
    Structural Finite elements should use these calc_options.
    Fluid      Finite elements can use these calc_options with similar meaning,
               there won't be a conflict
    Fluid      Finite elements can use these calc_options with other meaning,
               there won't be a conflict
    Fluid      Finite elements can use other calc_options and add them to,
               this list
               
    Some of these calc_options call an assembly to a sparse matrix or a dist. vector.!!!!               
               
                                                                  m.gee 01/02
*/
/*----------------------------------------------------------------------*/
{
int               i;
ELEMENT          *actele;
SPARSE_TYP        sysarray1_typ;
SPARSE_TYP        sysarray2_typ;


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
   default:
      dserror("Unknown typ of system matrix");
   break;
   }
}
#endif
/* =======================================================call elements */
/*---------------------------------------------- loop over all elements */
for (i=0; i<actpart->numele; i++)
{
   actele = actpart->element[i];
   switch(actele->eltyp)/*======================= call element routines */
   {
   case el_shell8:
      shell8(actfield,actpart,actintra,actele,&estif_global,&emass_global,dvec,global_numeq,kstep,calc_option);
   break;
   case el_brick1:
   break;
   case el_fluid1: 
   break;
   case el_fluid3:
   break;
   case el_ale:
   break;
   case el_none:
      dserror("Typ of element unknown");
   break;
   default:
      dserror("Typ of element unknown");
   }/* end of calling elements */


   switch(calc_option)/*=== call assembly dependent on calculation-flag */
   {
   case 1:/* in structural field, calculate linear stiffness only, static analysis */
      assemble(sysarray1,&estif_global,sysarray2,NULL,actpart,actsolv,actintra,actele,0);
   break;/* end of assembly */
   case 2:/* in structural field, calculate nonlinear stiffness only, static analysis */
      assemble(sysarray1,&estif_global,sysarray2,NULL,actpart,actsolv,actintra,actele,0);
   break;/* end of assembly */
   }
}/* end of loop over elements */
/*----------------------------------------------------------------------*/
/*                    in parallel coupled dofs have to be exchanged now */
/*             (if there are any inter-proc couplings, which is tested) */
/*----------------------------------------------------------------------*/
#ifdef PARALLEL 
switch(calc_option)
{
case 1: /* in structural field, do linear stiffness only, static analysis */
   assemble(sysarray1,NULL,sysarray2,NULL,actpart,actsolv,actintra,actele,1);
break;
case 2: /* in structural field, do nonlinear stiffness only, static analysis */
   assemble(sysarray1,NULL,sysarray2,NULL,actpart,actsolv,actintra,actele,1);
break;
}
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
void calinit(FIELD      *actfield,   /* the actove physical field */ 
             PARTITION  *actpart)    /* my partition of this field */
{
int i;                        /* a counter */
int is_shell8=0;              /* flags to check for presents of certain element types */
int is_brick1=0;
int is_fluid1=0;
int is_fluid3=0;
int is_ale=0;
ELEMENT *actele;              /* active element */
#ifdef DEBUG 
dstrc_enter("calinit");
#endif
/*----------------------------------------------------------------------*/
/*-------------------------- define dense element matrices for assembly */
amdef("estif",&estif_global,(MAXNOD*MAXDOFPERNODE),(MAXNOD*MAXDOFPERNODE),"DA");
amdef("emass",&emass_global,(MAXNOD*MAXDOFPERNODE),(MAXNOD*MAXDOFPERNODE),"DA");
/*--------------------what kind of elements are there in this example ? */
for (i=0; i<actfield->numele; i++)
{
   actele = &(actfield->element[i]);
   switch(actele->eltyp)
   {
   case el_shell8:
      is_shell8=1;
   break;
   case el_brick1:
      is_brick1=1;
   break;
   case el_fluid1:
      is_fluid1=1;
   break;
   case el_fluid3:
      is_fluid3=1;
   break;
   case el_ale:
      is_ale=1;
   break;
   default:
      dserror("Unknown typ of element");
   break;   
   }
}/* end of loop over all elements */
/*--------------------- init the element routines for all present types */
/*------------------------------- init all kind of routines for shell8  */
if (is_shell8==1)
{
   shell8(actfield,actpart,NULL,NULL,&estif_global,&emass_global,NULL,0,0,0);
}
/*-------------------------------- init all kind of routines for brick1 */
if (is_brick1==1)
{
}
/*-------------------------------- init all kind of routines for fluid1 */
if (is_fluid1==1)
{
}
/*-------------------------------- init all kind of routines for fluid3 */
if (is_fluid3==1)
{
}
/*----------------------------------- init all kind of routines for ale */
if (is_ale==1)
{
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
void calreduce(FIELD      *actfield, /* the active field */
               PARTITION  *actpart,  /* my partition of this field */
               INTRA      *actintra, /* the field's intra-communicator */
               int         kstep)    /* the actual time or incremental step */
{
int i;
int is_shell8=0;
int is_brick1=0;
int is_fluid1=0;
int is_fluid3=0;
int is_ale=0;
ELEMENT *actele;
#ifdef DEBUG 
dstrc_enter("calreduce");
#endif
/*----------------------------------------------------------------------*/
/*--------------------what kind of elements are there in this example ? */
for (i=0; i<actfield->numele; i++)
{
   actele = &(actfield->element[i]);
   switch(actele->eltyp)
   {
   case el_shell8:
      is_shell8=1;
   break;
   case el_brick1:
      is_brick1=1;
   break;
   case el_fluid1:
      is_fluid1=1;
   break;
   case el_fluid3:
      is_fluid3=1;
   break;
   case el_ale:
      is_ale=1;
   break;
   default:
      dserror("Unknown typ of element");
   break;   
   }
}/* end of loop over all elements */
/*--------------------- init the element routines for all present types */
/*------------------------------- init all kind of routines for shell8  */
if (is_shell8==1)
{
   shell8(actfield,actpart,actintra,NULL,NULL,NULL,NULL,0,kstep,7);
}
/*-------------------------------- init all kind of routines for brick1 */
if (is_brick1==1)
{
}
/*-------------------------------- init all kind of routines for fluid1 */
if (is_fluid1==1)
{
}
/*-------------------------------- init all kind of routines for fluid3 */
if (is_fluid3==1)
{
}
/*----------------------------------- init all kind of routines for ale */
if (is_ale==1)
{
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of calreduce */










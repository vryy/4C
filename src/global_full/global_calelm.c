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
void calelm(FIELD      *actfield, 
               SOLVAR     *actsolv, 
               PARTITION  *actpart, 
               INTRA      *actintra,
               int         sysarray1,
               int         sysarray2,
               double     *dvec,
               int         global_numeq,
               int         kstep,
               int         calc_option)
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
   }/* end of element routines */


   switch(calc_option)/*=== call assembly dependent on calculation-flag */
   {
   case 1:/*-----------calculate linear stiffness only, static analysis */
      assemble(sysarray1,&estif_global,sysarray2,NULL,actpart,actsolv,actintra,actele,0);
   break;/* end of assembly */
   case 2:/*--------calculate nonlinear stiffness only, static analysis */
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
case 1: /*--------------------do linear stiffness only, static analysis */
   assemble(sysarray1,NULL,sysarray2,NULL,actpart,actsolv,actintra,actele,1);
break;
case 2: /*-----------------do nonlinear stiffness only, static analysis */
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
void calinit(FIELD      *actfield, 
             PARTITION  *actpart)
{
int i;
int is_shell8=0;
int is_brick1=0;
int is_fluid1=0;
int is_fluid3=0;
int is_ale=0;
ELEMENT *actele;
#ifdef DEBUG 
dstrc_enter("calinit");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------------------------- init the local matrices */
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
void calreduce(FIELD      *actfield, 
               PARTITION  *actpart,
               INTRA      *actintra,
               int         kstep)
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










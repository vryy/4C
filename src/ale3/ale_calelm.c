/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale_calelm' which calls the element routines 

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "../shell8/shell8.h"
#include "../wall1/wall1.h"
#include "../brick1/brick1.h"
#include "../fluid3/fluid3.h"
#include "ale3.h"
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
struct _ARRAY estif_global;
struct _ARRAY emass_global;
struct _ARRAY intforce_global;

/*! 
\addtogroup Ale 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief  routine to call elements                                            

<pre>                                                              mn 06/02 
This routine checks the ale-element-type and calls the appropriate
element-routine

</pre>
\param *actfield  FIELD       (i)  active field
\param *actsolv   SOLVAR      (i)  active solvar
\param *actpart   PARTITION   (i)  my partition of the field
\param *actintra  INTRA       (i)  my intra-communicator
\param sysarray1  INT         (i)  number of first sparse system matrix
\param sysarray2  INT         (i)  number of second sparse system matrix
\param *action    CALC_ACTION (i)  calculation option passed to elements

\warning There is nothing special to this routine
\return void                                               
\sa calling: ale2(), ale3(); called by: dyn_ale()

*----------------------------------------------------------------------*/

void ale_calelm(
    FIELD        *actfield,     /* active field */        
    SOLVAR       *actsolv,      /* active SOLVAR */
    PARTITION    *actpart,      /* my partition of this field */
    INTRA        *actintra,     /* my intra-communicator */
    INT           sysarray1,    /* number of first sparse system matrix */
    INT           sysarray2,    /* number of secnd system matrix,*/
    CALC_ACTION  *action)       /* calculation option passed to element routines */
/*----------------------------------------------------------------------*/
{
INT               i;
ELEMENT          *actele;
SPARSE_TYP        sysarray1_typ;
SPARSE_TYP        sysarray2_typ;
ASSEMBLE_ACTION   assemble_action;

#ifdef DEBUG 
dstrc_enter("ale_calelm");
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
   switch(actele->eltyp)/*======================= call element routines */
   {
   case el_ale3:
	ale3(actpart,actintra,actele,
        &estif_global,
        action);
   break;
   case el_ale2:
	ale2(actpart,actintra,actele,
        &estif_global,
        action);
   break;
   case el_none:
      dserror("Typ of element unknown");
   break;
   default:
      dserror("Typ of element unknown");
   }/* end of calling elements */


   switch(*action)/*=== call assembly dependent on calculation-flag */
   {
   case calc_ale_stiff           : assemble_action = assemble_one_matrix; break;
   case calc_ale_rhs             : assemble_action = assemble_do_nothing; break;
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
            assemble_action);
}/* end of loop over elements */
/*----------------------------------------------------------------------*/
/*                    in parallel coupled dofs have to be exchanged now */
/*             (if there are any inter-proc couplings, which is tested) */
/*----------------------------------------------------------------------*/
#ifdef PARALLEL 
switch(*action)
{
case calc_ale_stiff            : assemble_action = assemble_one_exchange; break;
case calc_ale_rhs              : assemble_action = assemble_do_nothing; break;
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
         assemble_action);
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ale_calelm */
#endif
/*! @} (documentation module close)*/

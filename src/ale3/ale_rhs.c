/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale_rhs' which calculates the rhs for
dirichlet boundary conditions for ale elements

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

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
 | global dense matrices for element routines             m.gee 7/01    |
 *----------------------------------------------------------------------*/
struct _ARRAY estif_global;
struct _ARRAY emass_global;
struct _ARRAY intforce_global;

/*! 
\addtogroup Ale 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief controls the calculation of the rhs for ale elements

<pre>                                                              mn 06/02 
This routine controls the calculation of the rhs for ale elements.

</pre>
\param *actfield    FIELD       (i)  active field
\param *actsolv     SOLVAR      (i)  active solvar
\param *actpart     PARTITION   (i)  my partition of the field
\param *actintra    INTRA       (i)  my intra-communicator
\param sysarray1    INT         (i)  number of first sparse system matrix
\param sysarray2    INT         (i)  number of second sparse system matrix
\param *dirich      DOUBLE      (i)  global redundant vector of dirichlet forces
\param global_numeq INT         (i)  size of dvec
\param kstep        INT         (i)  time in increment step we are in
\param *action      CALC_ACTION (i)  calculation option passed to elements

\warning There is nothing special to this routine
\return void                                               
\sa calling: ale2(), ale3(); ale_caldirich(); called by: dyn_ale

*----------------------------------------------------------------------*/
void ale_rhs(FIELD        *actfield,     /* active field */        
            SOLVAR       *actsolv,      /* active SOLVAR */
            PARTITION    *actpart,      /* my partition of this field */
            INTRA        *actintra,     /* my intra-communicator */
            INT           sysarray1,    /* number of first sparse system matrix */
            INT           sysarray2,    /* number of secnd system matrix, if present, else -1 */
            DOUBLE       *dirich,       /* global redundant vector of dirichlet forces */
            INT           global_numeq, /* size of dvec */
            INT           kstep,        /* time in increment step we are in */
            CONTAINER    *container,
            CALC_ACTION  *action)       /* calculation option passed to element routines */
/*----------------------------------------------------------------------*/
{
INT               i,j,k;
INT               hasdirich;
ELEMENT          *actele;
ASSEMBLE_ACTION   assemble_action;
GNODE            *actgnode;


#ifdef DEBUG 
dstrc_enter("ale_rhs");
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
   case oll:
      if (actsolv->sysarray[sysarray1].oll->couple_d_send)
         amzero(actsolv->sysarray[sysarray1].oll->couple_d_send);
      if (actsolv->sysarray[sysarray1].oll->couple_d_recv)
         amzero(actsolv->sysarray[sysarray1].oll->couple_d_recv);
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
   case oll:
      if (actsolv->sysarray[sysarray2].oll->couple_d_send)
         amzero(actsolv->sysarray[sysarray2].oll->couple_d_send);
      if (actsolv->sysarray[sysarray2].oll->couple_d_recv)
         amzero(actsolv->sysarray[sysarray2].oll->couple_d_recv);
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
   amzero(&intforce_global);
   hasdirich = 0;

   /* check if there are inhomogeneous dirich conditions present for this element */
   for (j=0; j<actele->numnp; j++)
   {
      actgnode = actele->node[j]->gnode;   
      if (actgnode->dirich==NULL) 
       continue;
      else
      {
	if (actgnode->dirich->dirich_type==dirich_freesurf)
	{
	   hasdirich=1;
	   goto out;
	}
	for(k=0; k<actele->node[j]->numdf; k++)
	{
	  if (actgnode->dirich->dirich_val.a.dv[k]!=0.0)
	  {
             hasdirich=1;
	     goto out;
          }
	  else if (actgnode->dirich->dirich_type==dirich_FSI)
	  {
             hasdirich=1;
	     goto out;
          }	  
        }
      }
   }   					  
   out:
   *action = calc_ale_stiff;
   /*----------------------------------------------------------------------*/
   if (hasdirich==1) /* --> nodes with DBC for this element */
   {
     switch(actele->eltyp)/*======================= call element routines */
     {
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

     /*------ assemble the rhs vector of condensed dirichlet conditions */
     ale_caldirich(actele,dirich,global_numeq,&estif_global);
   }
   /*----------------------------------------------------------------------*/
}/* end of loop over elements */
/*----------------------------------------------------------------------*/
/*                    in parallel coupled dofs have to be exchanged now */
/*             (if there are any inter-proc couplings, which is tested) */
/*----------------------------------------------------------------------*/
#ifdef PARALLEL 
*action = calc_ale_rhs;
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
         assemble_action,
         container);
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ale_rhs */
#endif
/*! @} (documentation module close)*/

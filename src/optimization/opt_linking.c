/*!----------------------------------------------------------------------
\file
\brief subroutine to encode and decode grdobj, grdcon and var when linking is used   AS 03/03

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0711 - 685-6575
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET

/*----------------------------------------------------------------------*/
#ifdef D_OPTIM                   /* include optimization code to ccarat */
/*----------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#include "../headers/optimization.h"
#include "opt_prototypes.h"
/*!
\addtogroup OPTIMIZATION
*//*! @{ (documentation module open)*/



/*!----------------------------------------------------------------------
\brief the optimization main structure
<pre>                                                            al 06/01
defined in opt_cal_main.c
</pre>
*----------------------------------------------------------------------*/
 struct _OPTI *opt;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | global variable *partition, vector of lenght numfld of structures    |
 | PARTITION is defined in global_control.c                             |
 *----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
extern enum _CALC_ACTION calc_action[MAXFIELD];
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate static variables if needed                       |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _STATIC_VAR  *statvar;
/*----------------------------------------------------------------------*
 |                                                         al 06/02     |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


void linking (INT phase, DOUBLE *grdobj, DOUBLE *grdcon, DOUBLE *var,
             DOUBLE *grdobj_lin, DOUBLE *grdcon_lin, DOUBLE *var_lin)
{
/* -------------------------------------------------------------------------- */
  INT i, numvarlin;              /* counter                          */
  INT n, o;                      /* new position, old position       */
  ELEMENT *actele;              /* active element              */
  PARTITION    *actpart;        /* pointer to the fields PARTITION structure */
  FIELD        *actfield;       /* pointer to the structural FIELD */
  CALC_ACTION  *action;         /* pointer to the structures cal_action enum */
/*----------------------------------------------------------------------*/
CONTAINER     container;        /* contains variables defined in container.h */
/*--------------------------------------------------- set some pointers */
  actfield    = &(field[0]);
  actpart     = &(partition[0]);
  action      = &(calc_action[0]);
  container.fieldtyp  = actfield->fieldtyp;
/* ---------------------------------------------------------------------------------- */
  numvarlin = opt->strat.fsd->numvar_lin;
/* --------------------------------------- write linked grdobj and grdcon vectors --- */
  if(phase == 1)
  {
   for(i=0; i<actfield->dis[0].numele; i++)
   {
    actele     = &(actfield->dis[0].element[i]);
    if(actele->optdata==NULL) continue; /* element does not take part in opt. */
    if(actele->optdata[0]==0) continue; /* position in variable vector        */
    o = actele->optdata[0];	      /* old unlinked position  	    */
    n = actele->optdata[2];	      /* new linked position		    */
    grdobj_lin[n-1] += grdobj[o-1] * actele->mylinweight;
    grdcon_lin[n-1] += grdcon[o-1];
    /* grdcon additive, because: c=(1-m/m_{ref}) = 1-1/m_{ref}* integral\rho dV
                                   \approx 1-1/m_{ref}\sum \rho_{i} V_{i}
     \partial c / \partial \rho_{i}=-1/m_{ref} V_{i}   */
   }

   for(i=0; i< numvarlin; i++)
   {
    var_lin[i] = var[i];
   }
  }

/* ------------------------------ rewrite var from OC algo to original var length --- */
  if(phase == 2)
  {
   for(i=0; i<actfield->dis[0].numele; i++)
   {
    actele     = &(actfield->dis[0].element[i]);
    if(actele->optdata==NULL) continue; /* element does not take part in opt. */
    if(actele->optdata[0]==0) continue; /* position in variable vector        */
    o = actele->optdata[0];           /* old unlinked position              */
    n = actele->optdata[2];           /* new linked position                */
    var[o-1] = var_lin[n-1];
   }
  printf("variables are linked\n");
  }

}
/*----------------------------------------------------------------------*/
#endif /* stop including optimization code to ccarat :*/
/*----------------------------------------------------------------------*/

/*! @} (documentation module close)*/

/* linking prerequisites:
   ----------------------
   - equal numbers of elements in materials to be linked
   - elements to be linked are numbered in the same sense
   - OV-LIN followed by number of linking rules
   - linking rules: link material 2 with 1
                    link material 3 with 2
		    link material 4 with 3 and so on */

#endif

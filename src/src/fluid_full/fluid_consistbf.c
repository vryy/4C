/*!---------------------------------------------------------------------
\file
\brief controlling calculation of consistent boundary forces

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6572
</pre>

---------------------------------------------------------------------*/
#ifndef CCADISCRET
/*!
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fluid_prototypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid2/fluid2_prototypes.h"
#include "../fluid2_pro/fluid2pro.h"
#include "../fluid2_pro/fluid2pro_prototypes.h"
#include "../fluid2_is/fluid2_is.h"
#include "../fluid3/fluid3.h"
#ifdef D_ALE
#include "../ale2/ale2.h"
#endif

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PAR   par;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 | global dense matrices for element routines             m.gee 7/01    |
 *----------------------------------------------------------------------*/
extern struct _ARRAY eforce_global;   /* element RHS                    */



/*!---------------------------------------------------------------------
\brief controls the evaluation of consistent boundary forces for l&d

<pre>                                                        chfoe 09/04
This routine controls the evaluation of consistent boundary forces along
dirichlet boundary lines for lift & drag values.
Within an init phase (init=0) all elements which participate at a liftdrag
line are searched and labled by setting the flag
actele->e.f2->force_on on the respective liftdrag evaluation flag.
In the main part these elements are visited again, its interesting nodes
(the ones which are situated on the lift&drag line) are found and the
elemental residual vector is calculated.
Lift&drag forces are the respective sum of the elemental residual entities.

</pre>
\param   *actpdis        PARTDISCRET    (i)
\param   *container      CONTAINER      (i)
\param   *ipos                          (i)   node array positions
\param    init           INT            (i)   init flag

\warning At the moment only elements with the same number of dofs at all
         nodes can be treated.

\return void

------------------------------------------------------------------------*/
void fluid_cbf(PARTDISCRET *actpdis,
	       CONTAINER   *container,
               ARRAY_POSITION *ipos,
               INT          init)
{

#ifdef D_FLUID2
INT i, j, k;
INT ld_id;                            /* Id of the actual liftdrag line */
INT numdf;
INT hasldline;
INT hasdirich, hasext;
INT force_on_node[MAXNOD];

DOUBLE   rho;

DOUBLE   xforce, yforce;
DOUBLE   center[2];
DOUBLE   xy[2];                  /* relative coordinates of actual node */
/*static DOUBLE **estif;*/           /* pointer to global ele-stif          */
/*static DOUBLE  *edforce;*/

NODE    *actnode;
GNODE   *actgnode;
ELEMENT *actele;
GLINE   *actgline;
DLINE   *actdline;

#endif


static INT num_ldele; /* number of elements along lift&drag line */
static DOUBLE  *eforce;
static DOUBLE  *liftdrag;        /* pointer to liftdrag field           */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("fluid_calc_cforce");
#endif

/*--------------------------------------------------- initialisation ---*/
if (init==0)
{
   num_ldele = 0;

   /*--- preselect elements that have lift&drag line or fsi coupling ---*/
   switch (genprob.ndim)
   {
   case 2:
#ifdef D_FLUID2
      for(i=0; i<actpdis->numele; i++)
      {
        INT found_ld = -1;
         actele = actpdis->element[i];
         if (par.myrank!=actele->proc) continue; /* only my elments */
         hasldline = 0;
         for(j=0; j<actele->g.gsurf->ngline; j++)
         {
            actgline = actele->g.gsurf->gline[j];
            if (actgline->dline == NULL) continue; /* not on boundary line */
            actdline = actgline->dline;
            if (actdline->liftdrag == NULL) continue; /* not on l&d line   */
            hasldline++; /* actele has lift&drag line */
            if(found_ld != 0 && found_ld != found_ld)
              dswarning(1,13);
            found_ld = actdline->liftdrag->liftdrag;
         }
         /*--- care for elements which have only nodes touching ld line ---*/
         if(!hasldline)
         {
            for(j=0; j<actele->numnp; j++)
            {
               actnode  = actele->node[j];
               actgnode = actnode->gnode;
               for(k=0; k<actgnode->ngline; k++)
               {
                  actgline = actgnode->gline[k];
                  if (actgline->dline == NULL) continue; /* not on boundary line */
                  actdline = actgline->dline;
                  if (actdline->liftdrag == NULL) continue; /* not on l&d line   */
                  hasldline++;
                  if(found_ld != 0 && found_ld != found_ld)
                    dswarning(1,13);
                  found_ld = actdline->liftdrag->liftdrag;
               }
            }
         }
         if (hasldline)
         {
           num_ldele++;         /* count number of elements along ld-line */

           switch (actele->eltyp)
           {
           case el_fluid2:
             if(actele->e.f2->force_on != 0 &&
                actele->e.f2->force_on != found_ld)
               dswarning(1,13);
             actele->e.f2->force_on = found_ld;
             break;
#ifdef D_FLUID2_PRO
           case el_fluid2_pro:
             if(actele->e.f2pro->force_on != 0 &&
                actele->e.f2pro->force_on != found_ld)
               dswarning(1,13);
             actele->e.f2pro->force_on = found_ld;
             break;
#endif
#ifdef D_FLUID2_IS
           case el_fluid2_is:
             if(actele->e.f2is->force_on != 0 &&
                actele->e.f2is->force_on != found_ld)
               dswarning(1,13);
             actele->e.f2is->force_on = found_ld;
             break;
#endif
           default:
             dserror("element type %d unknown", actele->eltyp);
           }
         }
      }
#endif
   break;
   case 3:
#ifdef D_FLUID3
      dserror("liftdrag nodeforce for 3D not yet implemented! (go for liftdrag stress instead)");
#endif
   break;
   }
   goto end;
}

liftdrag = container->liftdrag;

eforce  = eforce_global.a.dv;

/*--------------- evaluate lift & drag forces depending on dimension ---*/
switch (genprob.ndim)
{
case 2: /* problem is two-dimensional */
#ifdef D_FLUID2
   rho = mat[actpdis->element[0]->mat-1].m.fluid->density;
   numdf = 3;
   for(i=0; i<actpdis->numele; i++)
   {
      actele = actpdis->element[i];

     /*------------------------------------------------------------------*/
      switch (actele->eltyp)
      {
      case el_fluid2:
        if (actele->e.f2->force_on==0)
          continue; /* element not of interest*/
        /*------------------------------ get liftdrag Id of element ---*/
        ld_id = actele->e.f2->force_on;
        break;
#ifdef D_FLUID2_PRO
      case el_fluid2_pro:
        if (actele->e.f2pro->force_on==0)
          continue; /* element not of interest*/
        /*------------------------------ get liftdrag Id of element ---*/
        ld_id = actele->e.f2pro->force_on;
        break;
#endif
#ifdef D_FLUID2_IS
      case el_fluid2_is:
        if (actele->e.f2is->force_on==0)
          continue; /* element not of interest*/
        /*------------------------------ get liftdrag Id of element ---*/
        ld_id = actele->e.f2is->force_on;
        break;
#endif
      default:
        dserror("element type!");
      }

      for(j=0; j<MAXNOD; j++)
        force_on_node[j] = -1;
      hasldline = 0;
      for(j=0; j<actele->g.gsurf->ngline; j++)
      {
         actgline = actele->g.gsurf->gline[j];
         if (actgline->dline == NULL) continue; /* not on boundary line */
         actdline = actgline->dline;
         if (actdline->liftdrag == NULL) continue; /* not on l&d line   */
         hasldline++; /* actele has lift&drag line */
         /* get center of liftdrag momentum (only one per node!) */
         center[0] = actdline->liftdrag->ld_center[0];
         center[1] = actdline->liftdrag->ld_center[1];

         /*--------- get list of Ids of interesting nodes at element ---*/
         for(k=0; k<actgline->ngnode; k++)
         {
           INT l;
           for (l=0; l<actele->numnp; l++)
           {
             if (actele->node[l]==actgline->gnode[k]->node)
             {
               force_on_node[l] = actgline->gnode[k]->node->Id;
               break;
             }
           }
         }
      }
      /*--- care for elements which have only nodes touching ld line ---*/
      if (!hasldline)           /* one node per element only!!! */
      {
         for(j=0; j<actele->numnp; j++)
         {
            actnode  = actele->node[j];
            actgnode = actnode->gnode;
            for(k=0; k<actgnode->ngline; k++)
            {
               actgline = actgnode->gline[k];
               if (actgline->dline == NULL) continue; /* not on boundary line */
               actdline = actgline->dline;
               if (actdline->liftdrag == NULL) continue; /* not on l&d line   */
               hasldline++;
               force_on_node[j] = actnode->Id;
               /* center of liftdrag angular momentum (only one per node!) */
               center[0] = actdline->liftdrag->ld_center[0];
               center[1] = actdline->liftdrag->ld_center[1];
            }
         }
      }
      if (hasldline)
      {
        INT dof;

         /*--- get force vector ---*/

         switch (actele->eltyp)
         {
         case el_fluid2:
           f2_caleleres(actele,&eforce_global,ipos,&hasdirich,&hasext);
           break;
#ifdef D_FLUID2_PRO
         case el_fluid2_pro:
           f2pro_caleleres(actele,&eforce_global,ipos,&hasdirich,&hasext);
           break;
#endif
#ifdef D_FLUID2_IS
         case el_fluid2_is:
           f2is_caleleres(actele,&eforce_global,ipos,&hasdirich,&hasext);
           break;
#endif
         default:
           dserror("element type!");
         }

         /*--------------------- perform matrix-vector multiplication...
           ------------------------------------ ...on selected lines ---*/
         dof = 0;
         for(j=0; j<actele->numnp; j++)
         {
            actnode = actele->node[j];

            if (force_on_node[j]!=-1)
            {
              /* warning: this is not yet for the ale case!!! */
              xy[0] = actnode->x[0] - center[0];
              xy[1] = actnode->x[1] - center[1];

              xforce = eforce[dof  ]*rho;
              yforce = eforce[dof+1]*rho;

              /* write nodal result from this ele to total ld sum */
              liftdrag[(ld_id-1)*6+0] += xforce;
              liftdrag[(ld_id-1)*6+1] += yforce;
              liftdrag[(ld_id-1)*6+5] += (xforce * xy[1] - yforce*xy[0]);
            }
            dof += actnode->numdf;
         }
      }
   }
#endif
break;
case 3: /* problem is three-dimensional */
#ifdef D_FLUID3
   dserror("consistent nodal fluid forces for 3D not yet implemented");
#endif
break;
default:
   dserror("genprob->ndim not valid");
}

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
}

#endif
#endif

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
/*!
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../solver/solver.h" 
#include "fluid_prototypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid2/fluid2_prototypes.h"
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
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;

static FLUID_DYNAMIC *fdyn;
/*!----------------------------------------------------------------------
\brief positions of physical values in node arrays

<pre>                                                        chfoe 11/04

This structure contains the positions of the various fluid solutions 
within the nodal array of sol_increment.a.da[ipos][dim].

extern variable defined in fluid_service.c
</pre>

------------------------------------------------------------------------*/
extern struct _FLUID_POSITION ipos;

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
\param    init           INT            (i)   init flag

\warning At the moment only elements with the same number of dofs at all
         nodes can be treated.

\return void

------------------------------------------------------------------------*/
void fluid_cbf(PARTDISCRET *actpdis,
	       CONTAINER   *container,
               INT          init)
{
INT i, j, k, l;
INT line;
INT ld_id;                            /* Id of the actual liftdrag line */
INT numdf;
INT foundsome, foundit;
INT hasdirich, hasext;
INT force_on_node[MAXNOD];
INT nfnode;  /* number of nodes of actele where forces are searched for */
static INT num_ldele; /* number of elements along lift&drag line */

DOUBLE   timefac;
DOUBLE   xforce, yforce;
DOUBLE   center[2];
DOUBLE   xy[2];                  /* relative coordinates of actual node */
static DOUBLE **estif;           /* pointer to global ele-stif          */
static DOUBLE  *eforce;
static DOUBLE  *edforce;
static DOUBLE  *liftdrag;        /* pointer to liftdrag field           */

NODE    *actnode;
GNODE   *actgnode;
ELEMENT *actele;
GLINE   *actgline;
DLINE   *actdline;

DOUBLE testforce;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("fluid_calc_cforce");
#endif

/*--------------------------------------------------- initialisation ---*/
if (init==0)
{
   fdyn      = alldyn[genprob.numff].fdyn;
   num_ldele = 0;

   /*--- preselect elements that have lift&drag line or fsi coupling ---*/
   switch (genprob.ndim)
   {
   case 2:
#ifdef D_FLUID2
      for(i=0; i<actpdis->numele; i++)
      {
         actele = actpdis->element[i];
         if (par.myrank!=actele->proc) continue; /* only my elments */
         foundsome = 0;
         for(j=0; j<actele->g.gsurf->ngline; j++)
         {
            actgline = actele->g.gsurf->gline[j];
            if (actgline->dline == NULL) continue; /* not on boundary line */
            actdline = actgline->dline;
            if (actdline->liftdrag == NULL) continue; /* not on l&d line   */
            foundsome++; /* actele has lift&drag line */
            if(actele->e.f2->force_on != 0 && 
               actele->e.f2->force_on != actdline->liftdrag->liftdrag)
               dswarning(1,13);
            actele->e.f2->force_on = actdline->liftdrag->liftdrag;
         }
         /*--- care for elements which have only nodes touching ld line ---*/
         if(!foundsome)
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
                  foundsome++;
                  if(actele->e.f2->force_on != 0 && 
                     actele->e.f2->force_on != actdline->liftdrag->liftdrag)
                     dswarning(1,13);
                  actele->e.f2->force_on = actdline->liftdrag->liftdrag;
               }
            }
         }
         if (foundsome) num_ldele++; /* count number of elements along ld-line */
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
timefac  = 1.0/fdyn->thsl;

eforce  = eforce_global.a.dv;

/*--------------- evaluate lift & drag forces depending on dimension ---*/
switch (genprob.ndim)
{
case 2: /* problem is two-dimensional */
#ifdef D_FLUID2
   numdf = 3;
   for(i=0; i<actpdis->numele; i++)
   { 
      
      actele = actpdis->element[i];

     /*------------------------------------------------------------------*/
      if (actele->e.f2->force_on==0) continue; /* element not of interest*/
      nfnode = 0;
      for(j=0; j<MAXNOD; j++)
        force_on_node[j] = -1;
      foundsome = 0;
      for(j=0; j<actele->g.gsurf->ngline; j++)
      {
         actgline = actele->g.gsurf->gline[j];
         if (actgline->dline == NULL) continue; /* not on boundary line */
         actdline = actgline->dline;
         if (actdline->liftdrag == NULL) continue; /* not on l&d line   */
         foundsome++; /* actele has lift&drag line */
         /* get center of liftdrag momentum (only one per node!) */
         center[0] = actdline->liftdrag->ld_center[0];
         center[1] = actdline->liftdrag->ld_center[1];
         
         /*--------- get list of Ids of interesting nodes at element ---*/
         for(k=0; k<actgline->ngnode; k++)
         {
            l=0;
            foundit=0;
            if((actgline->gnode[k]->node->x[0] == 0.0 &&
                actgline->gnode[k]->node->x[1] == 1.0) ||
               (actgline->gnode[k]->node->x[0] == 0.0 &&
                actgline->gnode[k]->node->x[1] == 1.0)) continue;
            while(force_on_node[l]!=-1 && l<MAXNOD)
            {
               if(force_on_node[l]==actgline->gnode[k]->node->Id)
                  foundit++;
               l++;
            }
            if(foundit) continue;
            force_on_node[l] = actgline->gnode[k]->node->Id;
            nfnode++;
         }
      }
      /*--- care for elements which have only nodes touching ld line ---*/
      if(!foundsome)  /* one node per element only!!! */
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
               foundsome++;
               force_on_node[0] = actnode->Id;
               if(foundsome==1) nfnode++;
               /* center of liftdrag angular momentum (only one per node!) */
               center[0] = actdline->liftdrag->ld_center[0];
               center[1] = actdline->liftdrag->ld_center[1];
            }
         }
      }
      if(foundsome)
      {
         /*------------------------------ get liftdrag Id of element ---*/
         ld_id = actele->e.f2->force_on;
         /*------------------ which nodes of actele are of interest? ---*/
         j=0;
         while(force_on_node[j]!=-1 && j<MAXNOD)
         {
            foundit=0;
            l=0;
            while(!foundit && l<actele->numnp)
            {
               if(force_on_node[j]==actele->node[l]->Id)
               {
                  force_on_node[j]=l;
                  foundit++;
               }
               l++;
            }
            j++;
            dsassert(foundit,"something is wrong!");
         }
         /*--- get force vector ---*/
         f2_caleleres(actele,&eforce_global,&hasdirich,&hasext);

         /*--------------------- perform matrix-vector multiplication...
           ------------------------------------ ...on selected lines ---*/
         for(j=0; j<nfnode; j++)
         {
            actnode = actele->node[force_on_node[j]];
            /*warning: this is not yet for the ale case!!! */
            xy[0] = actnode->x[0] - center[0];
            xy[1] = actnode->x[1] - center[1];
            xforce = yforce = 0.0;
            testforce = 0.0;
            line = force_on_node[j] * 3;

            xforce += eforce[line]   * timefac;
            yforce += eforce[line+1] * timefac;
            
            testforce += eforce[24] * timefac;
            
            /* write nodal result from this ele to total ld sum */
            liftdrag[(ld_id-1)*6+0] += xforce;
            liftdrag[(ld_id-1)*6+1] += yforce;
            liftdrag[(ld_id-1)*6+5] += (xforce * xy[1] - yforce*xy[0]);
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

/*!----------------------------------------------------------------------
\file
\brief setting dirichlet conditions for fluid

<pre>
Maintainer: Thomas Hettich
            hettich@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hettich/
            0771 - 685-6575
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fluid_prototypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;

static FLUID_DYNAMIC *fdyn;
/*!---------------------------------------------------------------------
\brief routine to set dirichlet boundary conditions at time <time>

<pre>                                                         he  12/02

in this routine the dirichlet boundary conditions for fluid2 and fluid3
elements are set at time <T=fdyn->time>.
the actual dirichlet values are written to the solution history of the
nodes:
    'actnode->sol_increment.a.da[3][j]
                                 |
                            time (n+1)

</pre>
\param *actfield    FIELD         (i)  actual field (fluid)
\param *lower_limit_kappa  DOUBLE (o) lower limit for kappa
\param *lower_limit_kappa  DOUBLE (o) lower limit for epsilon
\return void
------------------------------------------------------------------------*/
void fluid_setdirich_tu(
                       FIELD  *actfield,
                       DOUBLE   *lower_limit_kappa,
                       DOUBLE   *lower_limit_eps
                       )
{
INT        i,j;
INT        numnp_total;              /* total number of fluid nodes     */
INT        numele_total;             /* total number of fluid elements  */
INT        numdf;	                   /* number of fluid dofs       	*/
DOUBLE     k_2,int_lenght;
GNODE     *actgnode;	             /* actual GNODE		            */
NODE      *actnode;	             /* actual NODE		            */
NODE      *actnode2;	             /* NODE from RANS                  */

#ifdef DEBUG
dstrc_enter("fluid_setdirich_tu");
#endif

/*----------------------------------------------------- set some values */
fdyn = alldyn[genprob.numff].fdyn;

numnp_total  = actfield->dis[1].numnp;
numele_total = actfield->dis[1].numele;
numdf        = 1;
int_lenght   = fdyn->lenght;

*lower_limit_kappa=0.0;
*lower_limit_eps  =0.0;

/*-------------------- loop all nodes and set actual dirichlet condition */
for (i=0;i<numnp_total;i++)
{
   actnode  = &(actfield->dis[1].node[i]);
   actnode2 = &(actfield->dis[0].node[i]);
   actgnode = actnode->gnode;
   if (actgnode->dirich==NULL)
   continue;
   for (j=0;j<numdf;j++) /* loop dofs */
   {
      if (actgnode->dirich->dirich_onoff.a.iv[j]==0 || actgnode->dirich->dirich_onoff.a.iv[j+1]==0)
         continue;
      else
      k_2 = 0.0000375*(pow(actnode2->sol_increment.a.da[3][j],2)+pow(actnode2->sol_increment.a.da[3][j+1],2));

      actnode->sol_increment.a.da[3][j]   = k_2 ;
      actnode->sol_increment.a.da[1][j]   = actnode->sol_increment.a.da[3][j];
      actnode->sol_increment.a.da[3][j+1] = 0.005*int_lenght*sqrt(actnode->sol_increment.a.da[3][j]);
      actnode->sol_increment.a.da[2][j+1] = actnode->sol_increment.a.da[3][j+1];
      actnode->sol_increment.a.da[3][j+2] = 0.09 * pow(actnode->sol_increment.a.da[3][j],1.5) / (0.005*int_lenght);
      actnode->sol_increment.a.da[1][j+2] = actnode->sol_increment.a.da[3][j+2];
      *lower_limit_kappa = DMAX(*lower_limit_kappa,actnode->sol_increment.a.da[3][j]*0.0001);
      *lower_limit_eps   = DMAX(*lower_limit_eps,actnode->sol_increment.a.da[3][j+2]*0.0001);

   } /*end loop over nodes */
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_settdirich*/
/*!---------------------------------------------------------------------
\brief routine to calculate the element dirichlet load vector

<pre>                                                         he 12/02

in this routine the element load vector due to dirichlet conditions
is calcluated. The prescribed values are taken from the node solution
history at (n+1) 'dirich[j] = actnode->sol_increment.a.da[3][j]'.
the element load vector 'dforce' is calculated by eveluating
</pre>
\code
      dforces[i] -= estif[i][j] * dirich[j];
\endcode

\param  *actele    ELEMENT   (i)   actual element
\param  *dforces   DOUBLE    (o)   dirichlet force vector
\param **estif     DOUBLE    (i)   element stiffness matrix
\param  *hasdirich INT       (o)   flag if s.th. was written to dforces

\return void

------------------------------------------------------------------------*/
void fluid_caldirich_tu(
                    ELEMENT   *actele,
		        DOUBLE    *dforces,
                    DOUBLE   **estif,
		        INT       *hasdirich
		       )
{

INT         i,j;
INT         numdf;                      /* number of fluid dofs         */
INT         nd=0;
DOUBLE      dirich[MAXDOFPERELE];       /* dirichlet values of act. ele */
INT         dirich_onoff[MAXDOFPERELE]; /* dirichlet flags of act. ele  */
GNODE      *actgnode;	                /* actual GNODE                 */
NODE       *actnode;	                /* actual NODE                  */

#ifdef DEBUG
dstrc_enter("fluid_caldirich_tu");
#endif

/*------------------------- check if there are any dirichlet conditions *
                                          for the nodes of this element */
for (i=0; i<actele->numnp; i++)
{
   actgnode = actele->node[i]->gnode;
   if (actgnode->dirich==NULL)
      continue;
   else
      *hasdirich=1;
      break;
} /* end loop over nodes */

if (*hasdirich==0) /* --> no nodes with DBC for this element */
   goto end;

/*---------------------------------- set number of dofs on this element */
for (i=0; i<actele->numnp; i++) nd += actele->node[i]->numdf;

/*---------------------------- init the vectors dirich and dirich_onoff */
for (i=0; i<nd; i++)
{
   dirich[i] = 0.0;
   dirich_onoff[i] = 0;
}

/*-------------------------------- fill vectors dirich and dirich_onoff */
/*                               dirichlet values at (n+1) were already */
/*                           written to the nodes (sol_increment[3][j]) */
for (i=0; i<actele->numnp; i++) /* loop nodes */
{
   numdf    = actele->node[i]->numdf;
   actnode  = actele->node[i];
   actgnode = actnode->gnode;
   for (j=0; j<numdf; j++) /* loop dofs */
   {
      if (actgnode->dirich==NULL || actgnode->dirich->dirich_onoff.a.iv[0]==0 ||
          actgnode->dirich->dirich_onoff.a.iv[1]==0) continue;
      dirich_onoff[i*numdf+j] = actgnode->dirich->dirich_onoff.a.iv[j];
      if(fdyn->kapeps_flag==0)
      dirich[i*numdf+j] = actnode->sol_increment.a.da[3][j];
      if(fdyn->kapeps_flag==1)
      dirich[i*numdf+j] = actnode->sol_increment.a.da[3][j+2];
   } /* end loop over dofs */
} /* end loop over nodes */
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*------------------------------------- do nothing for supported row */
   if (dirich_onoff[i]!=0) continue;
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      /*---------------------------- do nothing for unsupported columns */
      if (dirich_onoff[j]==0) continue;
      dforces[i] -= estif[i][j] * dirich[j];
   }/* loop j over columns */
}/* loop i over rows */

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_caldirich*/

#endif

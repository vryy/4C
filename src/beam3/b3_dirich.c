/*!----------------------------------------------------------------------
\file
\brief contains the routine 'b3_setdirich', which puts all dirichlet
values to the nodal solution

<pre>
Maintainer: Frank Huber
            huber@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/huber/
            0771 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BEAM3
#include "../headers/standardtypes.h"
#include "beam3.h"

/*!
\addtogroup BEAM3
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief puts all dirichlet values to the nodal solution

<pre>                                                              fh 12/02
This routine puts all values for dirichlet conditions a priori to the
nodal displacement solution vector

</pre>
\param *actfield FIELD    (i/o)  actual field


\warning There is nothing special in this routine
\return void
\sa calling:   ---;
    called by: beam3()

*----------------------------------------------------------------------*/
void b3_setdirich(FIELD     *actfield)
{
GNODE                *actgnode;     /* actual gnode */
NODE                 *actnode;      /* actual node */
INT                   i,j;          /* some loopers */
INT                   numnp_total;  /* sum of all nodes of the field */
DOUBLE                initval;

#ifdef DEBUG
dstrc_enter("b3_setdirich");
#endif

numnp_total  = actfield->dis[0].numnp;

/*------------------------------------------------- loop over all nodes */
for (i=0;i<numnp_total;i++)
{
   actnode  = &(actfield->dis[0].node[i]);
   actgnode = actnode->gnode;
   if (actgnode->dirich==NULL)
         continue;
   for (j=0;j<actnode->numdf;j++)
   {
      if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
         continue;
      initval  = actgnode->dirich->dirich_val.a.dv[j];
      actnode->sol.a.da[0][j] = initval;
   }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of b3_setdirich*/
#endif
/*! @} (documentation module close)*/

/*!----------------------------------------------------------------------
\file
\brief init ale field

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
/*!
\addtogroup Ale
*//*! @{ (documentation module open)*/
/*!----------------------------------------------------------------------
\brief

<pre>                                                         genk 03/03

</pre>

\return void

*----------------------------------------------------------------------*/
void ale_init(
    FIELD *actfield
    )
{

  INT i;
  INT numnp_total;
  NODE *actnode;

#ifdef DEBUG
  dstrc_enter("ale_init");
#endif

  /* set some values */
  numnp_total  = actfield->dis[0].numnp;

  /* allocate space for solution history */
  for (i=0;i<numnp_total;i++)
  {
    actnode=&(actfield->dis[0].node[i]);
    amredef(&(actnode->sol_increment),2,actnode->sol_increment.sdim,"DA");
    amzero(&(actnode->sol_increment));
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ale_init*/


#endif
/*! @} (documentation module close)*/

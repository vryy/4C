/*!----------------------------------------------------------------------
\file
\brief assign dofs

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID2_PRO
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2_PRO
#include "../headers/standardtypes.h"
/*!---------------------------------------------------------------------
\brief assign dofs to nodes for the Q2Q1 fluid element

<pre>                                                         genk 10/03


</pre>
\param  *actnode	         NODE     (i)
\param  *counter	         INT      (i)
\return void

------------------------------------------------------------------------*/

void f2pro_ass_dof_q2q1(NODE *actnode, INT *counter)
{
INT dirich;
INT couple;
INT geocouple;

#ifdef DEBUG
dstrc_enter("f2pro_ass_dof_q2q1");
#endif

/* this element has only one dof for the second discretisation */

dirich=0;
couple=0;
geocouple=0;
/*-------------------------- dof has dirichlet condition */
if (actnode->gnode->dirich!=NULL && actnode->gnode->dirich->dirich_onoff.a.iv[2]!=0)
dirich=1;
if (actnode->gnode->couple != NULL)
{
   dserror("coupling conditions for FLUID2_PRO Q2Q1 not implemented yet!\n");
}
/*-------------------------------------------------------*/
if (couple==1 && geocouple==1)
   dserror("geostationary dof coupling conflicts");
if (dirich==0 && couple==0 && geocouple==0)
{
   actnode->dof[0] = *counter; (*counter)++;
}
else if (dirich==1 && couple==0 && geocouple==0)
{
   actnode->dof[0]=-1;
}
else if (dirich==0 && couple==1 && geocouple==0)
{
   dserror("coupling conditions for FLUID2_PRO Q2Q1 not implemented yet!\n");
}
else if (dirich==0 && couple==0 && geocouple==1)
{
   dserror("coupling conditions for FLUID2_PRO Q2Q1 not implemented yet!\n");
}
else
   dserror("this is not OK! - you better go home now!\n");

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2pro_ass_dof_q2q1 */
#endif
/*! @} (documentation module close)*/

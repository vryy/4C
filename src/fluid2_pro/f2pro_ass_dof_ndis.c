/*!----------------------------------------------------------------------
\file
\brief assign dofs 

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
   /*---------- dof has geostationary coupling condition */
   /* if (actnode->gnode->couple->couple.a.ia[2][0]!=0) geocouple=1;
   /*------------------------ dof has coupling condition */
   /* if (actnode->gnode->couple->couple.a.ia[2][1]!=0)    couple=1; */
}
/*-------------------------------------------------------*/
if (couple==1 && geocouple==1) dserror("geostationary dof coupling conflicts");
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
/*  coupleID = actnode->gnode->couple->couple.a.ia[2][1];
   find_assign_coupset(actfield,coupleID,&counter); */
}
else if (dirich==0 && couple==0 && geocouple==1)
{
   dserror("coupling conditions for FLUID2_PRO Q2Q1 not implemented yet!\n");
   /* the coupling Id */
   /*coupleID = actnode->gnode->couple->couple.a.ia[2][0];
   /* find a geometrically compatibel node */
   /*iscouple_find_node_comp(
			actnode,
			actfield,
			&partnernode,
			coupleID,
			l);
   if (partnernode==NULL) dserror("Cannot do geostationary coupling"); 
   /* check wheher there already has been a dof assigned to partnernode */
   /*dof = partnernode->dof[l];
   if (dof==-2 && actnode->dof[0]==-2)
   {
      actnode->dof[l] = counter;
      partnernode->dof[l] = counter;
      counter++;
   }
   else if (dof!=-2 && actnode->dof[0]==-2)
   {
      actnode->dof[0] = dof;
   }
   else if (dof==-2 && actnode->dof[0]!=-2)
   {
      partnernode->dof[0] = actnode->dof[0];
   } */
}
 /* else if (dirich==1 && couple==0 && geocouple==1) 
 dserror("Case dirichlet condition in geocoupleset not yet implemented");*/

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2pro_ass_dof_q2q1 */
#endif
/*! @} (documentation module close)*/

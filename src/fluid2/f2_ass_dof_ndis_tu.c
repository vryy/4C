/*!----------------------------------------------------------------------
\file
\brief assign dofs to nodes for one fluid2 element for turbulence

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"

/*----------------------------------------------------------------------*
 |  put dofs to nodes                                   he    12/02     |
 |  for turbulence-models                                               |
 *----------------------------------------------------------------------*/
void f2tu_ass_dof(NODE *actnode, int *counter)
{
int dirich;
int couple;
int geocouple;

#ifdef DEBUG 
dstrc_enter("f2tu_ass_dof");
#endif

/* this element has only one dof for the second discretisation */

dirich=0;
couple=0;
geocouple=0;
/*-------------------------- dof has dirichlet condition */
if (actnode->gnode->dirich!=NULL && actnode->gnode->dirich->dirich_onoff.a.iv[0]!=0
&& actnode->gnode->dirich->dirich_onoff.a.iv[1]!=0) 
dirich=1;

if (actnode->gnode->couple != NULL)
{
   dserror("coupling conditions for FLUID2_TU not implemented yet!\n");
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
   dserror("coupling conditions for FLUID2_TU not implemented yet!\n");
/*  coupleID = actnode->gnode->couple->couple.a.ia[2][1];
   find_assign_coupset(actfield,coupleID,&counter); */
}
else if (dirich==0 && couple==0 && geocouple==1)
{
   dserror("coupling conditions for FLUID2_TU not implemented yet!\n");
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2tu_ass_dof */

#endif

#include "../headers/standardtypes.h"
#ifdef PARALLEL
#ifdef SUN
#include "../../../lib_sun/metis-4.0/Lib/metis.h"
#endif
#ifdef AZUSA
#include "../../../../lib_ita1/metis/metis.h"
#endif
#ifdef HPUX10
#include "/bau/stat33/users/statik/lib/METIS/metis.h"
#endif
#ifdef HPUX11
#include "/bau/stat33/users/statik/lib/METIS/metis.h"
#endif
#endif
/*----------------------------------------------------------------------*
 | find a partner coupling compatible node                m.gee 8/00    |
 *----------------------------------------------------------------------*/
void iscouple_find_node_comp(NODE  *actnode, 
                             FIELD *searchfield, 
                             NODE **partnernode,
                             int    coupleID,
                             int    dof)
{
int i,j,l;
int ierr;
double tol = EPS8;
#ifdef DEBUG 
dstrc_enter("iscouple_find_node_comp");
#endif
/*----------------------------------------------------------------------*/
*partnernode=NULL;
   for (i=0; i<searchfield->dis[0].numnp; i++)
   {
      /* no coupling conditions on this node */
      if (searchfield->dis[0].node[i].gnode->couple==NULL) continue;
      /* I do not find myself */
      if (searchfield->dis[0].node[i].Id_loc == actnode->Id_loc) continue;
         /* check for the right coupling set */
         /*
         Note: At the moment only the given dof is checked here */
         if (searchfield->dis[0].node[i].gnode->couple->couple.a.ia[dof][0]==coupleID)
         {
            /* check geometrical distance */
            cheque_distance(actnode->x,searchfield->dis[0].node[i].x,tol,&ierr);
            /* distance is too large ( > 1.0E-08 ) */
            if (ierr==0) continue;
            /* I found the correct node */
            else
            {
               *partnernode = &(searchfield->dis[0].node[i]);
               goto finish;
            }
         }
   }
finish:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of iscouple_find_node_comp */



/*----------------------------------------------------------------------*
 | check distance between nodes                           m.gee 8/00    |
 *----------------------------------------------------------------------*/
void cheque_distance(double *x1, double *x2, double tol, int *ierr)
{
int i,j;
double v[3];
double lenght;
#ifdef DEBUG 
dstrc_enter("cheque_distance");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<3; i++) v[i] = x2[i] - x1[i];
lenght = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
lenght = sqrt(lenght);
if (lenght <= tol) *ierr=1;
else *ierr=0;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of cheque_distance */





/*----------------------------------------------------------------------*
 | assign a dofcoupling set                               m.gee 5/01    |
 *----------------------------------------------------------------------*/
void find_assign_coupset(FIELD *actfield, 
                            int    coupleID, 
                            int   *counter)
{
int i,j,k,l;
int dof;
#ifdef DEBUG 
dstrc_enter("find_assign_coupset");
#endif
/*----------------------------------------------------------------------*/
/*----- check in couple set for already given dof / dirichlet condition */
/*
   if there is one dof in the coupling set which is dirichlet-conditioned, then the
   whole coupling set has to be dirichlet-conditioned 
*/
dof = -2;
/*------------------------ check for a dirichlet condition in coupleset */
for (i=0; i<actfield->dis[0].numnp; i++)
{
   if (actfield->dis[0].node[i].gnode->couple==NULL &&
       actfield->dis[0].node[i].gnode->dirich==NULL) continue;
   if (actfield->dis[0].node[i].gnode->couple==NULL) continue;
   for (l=0; l<actfield->dis[0].node[i].numdf; l++)
   {
      if (actfield->dis[0].node[i].gnode->couple->couple.a.ia[l][1]!=coupleID) continue;
      if (actfield->dis[0].node[i].gnode->dirich!=NULL)
      {
         if (actfield->dis[0].node[i].gnode->dirich->dirich_onoff.a.iv[l]!=0)
         {
            dof = -1;
            goto assign;
         }
      }
   }
}   
for (i=0; i<actfield->dis[0].numnp; i++)
{
   if (actfield->dis[0].node[i].gnode->couple==NULL &&
       actfield->dis[0].node[i].gnode->dirich==NULL) continue;
   if (actfield->dis[0].node[i].gnode->couple==NULL) continue;
   for (l=0; l<actfield->dis[0].node[i].numdf; l++)
   {
      if (actfield->dis[0].node[i].gnode->couple->couple.a.ia[l][1]!=coupleID) continue;
      if (actfield->dis[0].node[i].dof[l]!=-2)
      {
         dof = actfield->dis[0].node[i].dof[l];
         goto assign;
      }
   }
}
/*----------------------------------------------------------------------*/
assign:
if (dof != -2)
for (i=0; i<actfield->dis[0].numnp; i++)
{
   if (actfield->dis[0].node[i].gnode->couple==NULL &&
       actfield->dis[0].node[i].gnode->dirich==NULL) continue;
   if (actfield->dis[0].node[i].gnode->couple==NULL) continue;
   for (l=0; l<actfield->dis[0].node[i].numdf; l++)
   {
      if (actfield->dis[0].node[i].gnode->couple->couple.a.ia[l][1]!=coupleID) continue;
      
      if (actfield->dis[0].node[i].dof[l]==dof) goto finish;
      actfield->dis[0].node[i].dof[l]=dof;
   }
}
else
{
for (i=0; i<actfield->dis[0].numnp; i++)
{
   if (actfield->dis[0].node[i].gnode->couple==NULL &&
       actfield->dis[0].node[i].gnode->dirich==NULL) continue;
   if (actfield->dis[0].node[i].gnode->couple==NULL) continue;
   for (l=0; l<actfield->dis[0].node[i].numdf; l++)
   {
      if (actfield->dis[0].node[i].gnode->couple->couple.a.ia[l][1]!=coupleID) continue;
      
      actfield->dis[0].node[i].dof[l]=*counter;
   }
}
(*counter)=(*counter)+1;
}
finish:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of find_assign_coupset */

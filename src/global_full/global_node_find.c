#include "../headers/standardtypes.h"
#include "/bau/stat33/users/statik/lib/METIS/metis.h"

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
   for (i=0; i<searchfield->numnp; i++)
   {
      /* no conditions on this node */
      if (searchfield->node[i].c==NULL) continue;
      /* no coupling conditions on this node */
      if (searchfield->node[i].c->iscoupled==0) continue;
      /* I do not find myself */
      if (searchfield->node[i].Id_loc == actnode->Id_loc) continue;
         /* check for the right coupling set */
         /*
         Note: At the moment only the given dof is checked here */
         */
         if (searchfield->node[i].c->couple.a.ia[dof][0]==coupleID)
         {
            /* check geometrical distance */
            cheque_distance(actnode->x,searchfield->node[i].x,tol,&ierr);
            /* distance is too large ( > 1.0E-08 ) */
            if (ierr==0) continue;
            /* I found the correct node */
            else
            {
               *partnernode = &(searchfield->node[i]);
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
for (i=0; i<actfield->numnp; i++)
{
   if (actfield->node[i].c==NULL) continue;
   if (actfield->node[i].c->iscoupled==0) continue;
   for (l=0; l<actfield->node[i].numdf; l++)
   {
      if (actfield->node[i].c->couple.a.ia[l][1]!=coupleID) continue;
      if (actfield->node[i].c->isdirich==1)
      {
         if (actfield->node[i].c->dirich_onoff.a.iv[l]!=0)
         {
            dof = -1;
            goto assign;
         }
      }
   }
}   
for (i=0; i<actfield->numnp; i++)
{
   if (actfield->node[i].c==NULL) continue;
   if (actfield->node[i].c->iscoupled==0) continue;
   for (l=0; l<actfield->node[i].numdf; l++)
   {
      if (actfield->node[i].c->couple.a.ia[l][1]!=coupleID) continue;
      if (actfield->node[i].dof[l]!=-2)
      {
         dof = actfield->node[i].dof[l];
         goto assign;
      }
   }
}
/*----------------------------------------------------------------------*/
assign:
if (dof != -2)
for (i=0; i<actfield->numnp; i++)
{
   if (actfield->node[i].c==NULL) continue;
   if (actfield->node[i].c->iscoupled==0) continue;
   for (l=0; l<actfield->node[i].numdf; l++)
   {
      if (actfield->node[i].c->couple.a.ia[l][1]!=coupleID) continue;
      
      if (actfield->node[i].dof[l]==dof) goto finish;
      actfield->node[i].dof[l]=dof;
   }
}
else
{
for (i=0; i<actfield->numnp; i++)
{
   if (actfield->node[i].c==NULL) continue;
   if (actfield->node[i].c->iscoupled==0) continue;
   for (l=0; l<actfield->node[i].numdf; l++)
   {
      if (actfield->node[i].c->couple.a.ia[l][1]!=coupleID) continue;
      
      actfield->node[i].dof[l]=*counter;
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

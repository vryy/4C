#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 |  put dofs to nodes                                    m.gee 5/01     |
 |  the numbering of the dofs is done the same on all procs             |
 |  because it seems practical that all                                 |
 |  procs know the dof numbering of all fields.                         |
 |  (there is also no communication at all in this routine)             |
 *----------------------------------------------------------------------*/
void assign_dof()
{
int i,j,k,l;                       /* some counters */
int counter;                       
int coupleID;                      /* Id of a coupling set from gid */
int dof;                           /* dof in progress */
int couple,geocouple,dirich;       /* flags for conditions */
FIELD   *actfield;                 /* the actual field */
ELEMENT *actele;                   /* the element in progress */
NODE    *actnode, *partnernode;    /* node and coupling partner in progress */

#ifdef DEBUG 
dstrc_enter("assign_dof");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------------------- loop the meshes */
for (i=0; i<genprob.numfld; i++)
{
   actfield = &(field[i]);
/*--------------------------------------------- find new node numbering */
/* notice: bandwith opimization has been done by gid, so just do a
   field-local renumbering of this */
   for (j=0; j<actfield->dis[0].numnp; j++) actfield->dis[0].node[j].Id_loc = j;

/* loop all elements and put numdf to their nodes depending on kind of
   element. Notice that nodes which couple different kind of elements have
   the higher number of dofs */
   for (j=0; j<actfield->dis[0].numele; j++)
   {
      actele = &(actfield->dis[0].element[j]);
      switch(actele->eltyp)
      {
      case el_shell8:
         for (k=0; k<actele->numnp; k++)
         {
            if (actele->node[k]->numdf < NUMDOF_SHELL8) actele->node[k]->numdf=NUMDOF_SHELL8;
         }
         break;
      case el_brick1:
         for (k=0; k<actele->numnp; k++)
         {
            if (actele->node[k]->numdf < 3) actele->node[k]->numdf=3;
         }
         break;
      case el_wall1:
         for (k=0; k<actele->numnp; k++)
         {
            if (actele->node[k]->numdf < 2) actele->node[k]->numdf=2;
         }
         break;
      case el_fluid2:
         for (k=0; k<actele->numnp; k++)
         {
            if (actele->node[k]->numdf < 3) actele->node[k]->numdf=3;
         }
         break;
      case el_fluid3:
         for (k=0; k<actele->numnp; k++)
         {
            if (actele->node[k]->numdf < 4) actele->node[k]->numdf=4;
         }
         break;
      case el_ale3:
         for (k=0; k<actele->numnp; k++)
         {
            if (actele->node[k]->numdf < 3) actele->node[k]->numdf=3;
         }
         break;
      case el_ale2:
         for (k=0; k<actele->numnp; k++)
         {
            if (actele->node[k]->numdf < 2) actele->node[k]->numdf=2;
         }
         break;
      default:
         dserror("Unknown type of element, cannot assign number of dofs");
         break;
      }
   }
/*---------------------------------------- assign the dofs to the nodes */
   counter=0;
   for (j=0; j<actfield->dis[0].numnp; j++)
   {
      actfield->dis[0].node[j].dof  = (int*)CALLOC(actfield->dis[0].node[j].numdf,sizeof(int));
      if (!(actfield->dis[0].node[j].dof)) 
         dserror("Allocation of dof in NODE failed");
      /*------------------------- allocate the arrays to hold solutions */
      amdef("sol",&(actfield->dis[0].node[j].sol),1,actfield->dis[0].node[j].numdf,"DA");
      amzero(&(actfield->dis[0].node[j].sol));

      amdef("sol_incr",&(actfield->dis[0].node[j].sol_increment),1,actfield->dis[0].node[j].numdf,"DA");
      amzero(&(actfield->dis[0].node[j].sol_increment));

      amdef("sol_res",&(actfield->dis[0].node[j].sol_residual),1,actfield->dis[0].node[j].numdf,"DA");
      amzero(&(actfield->dis[0].node[j].sol_residual));
      /*------------------------------------------- init all dofs to -2 */
      for (l=0; l<actfield->dis[0].node[j].numdf; l++) actfield->dis[0].node[j].dof[l]=-2;
   }   
/*------- eliminate geostationary coupling conditions that conflict with 
   dofcoupling sets by putting the geostat coupling to the coupling set */
   coupleID=0;
   for (j=0; j<actfield->dis[0].numnp; j++)
   {
      actnode = &(actfield->dis[0].node[j]);
      if (actnode->gnode->couple==NULL) continue;
      for (l=0; l<actnode->numdf; l++)
      {
         if (actnode->gnode->couple->couple.a.ia[l][1]!=0)
         {
            coupleID = actnode->gnode->couple->couple.a.ia[l][1];
            /* found a conflict with geostationary coupling and delete it */
            if (actnode->gnode->couple->couple.a.ia[l][0]!=0)
            {
                iscouple_find_node_comp(
                                           actnode,
                                           actfield,
                                           &partnernode,
                                           actnode->gnode->couple->couple.a.ia[l][0],
                                           l);
                if (partnernode==NULL) dserror("Cannot do geostationary coupling");                                         
                partnernode->gnode->couple->couple.a.ia[l][1] =  coupleID;
                partnernode->gnode->couple->couple.a.ia[l][0] =  0;
                actnode->gnode->couple->couple.a.ia[l][0]     =  0;                                 
            }
         }
      }
   }

/*--------------------------------------------------------- assign dofs */
   for (j=0; j<actfield->dis[0].numnp; j++)
   {
      actnode = &(actfield->dis[0].node[j]);
/*------------------------------- the node does not have any conditions */
      if (actnode->gnode->couple==NULL && actnode->gnode->dirich==NULL)
      {
         for (l=0; l<actnode->numdf; l++)
         {
            actnode->dof[l] = counter;
            counter++;
         }
      }
/*--------------------------------------- the node does have conditions */
      else
      {
         /*----- the node does not have dirichlet or coupling condition */
         if (actnode->gnode->couple==NULL && actnode->gnode->dirich==NULL)
         {
            for (l=0; l<actnode->numdf; l++)
            {
               actnode->dof[l] = counter;
               counter++;
            }
         }
         else/* the node does have dirichlet and/or coupling conditions */
         {
            for (l=0; l<actnode->numdf; l++)
            {
               dirich=0;
               couple=0;
               geocouple=0;
               /*-------------------------- dof has dirichlet condition */
               if (actnode->gnode->dirich!=NULL && actnode->gnode->dirich->dirich_onoff.a.iv[l]!=0) 
               dirich=1;
               if (actnode->gnode->couple != NULL)
               {
                  /*---------- dof has geostationary coupling condition */
                  if (actnode->gnode->couple->couple.a.ia[l][0]!=0) geocouple=1;
                  /*------------------------ dof has coupling condition */
                  if (actnode->gnode->couple->couple.a.ia[l][1]!=0)    couple=1;
               }
               /*-------------------------------------------------------*/
               if (couple==1 && geocouple==1) dserror("geostationary dof coupling conflicts");
               if (dirich==0 && couple==0 && geocouple==0)
               {
                  actnode->dof[l] = counter; counter++;
               }
               else if (dirich==1 && couple==0 && geocouple==0)
               {
                  actnode->dof[l]=-1;
               }
               else if (dirich==0 && couple==1 && geocouple==0)
               {
                  coupleID = actnode->gnode->couple->couple.a.ia[l][1];
                  find_assign_coupset(actfield,coupleID,&counter);
               }
               else if (dirich==0 && couple==0 && geocouple==1)
               {
                  /* the coupling Id */
                  coupleID = actnode->gnode->couple->couple.a.ia[l][0];
                  /* find a geometrically compatibel node */
                  iscouple_find_node_comp(
                                             actnode,
                                             actfield,
                                             &partnernode,
                                             coupleID,
                                             l);
                  if (partnernode==NULL) dserror("Cannot do geostationary coupling"); 
                  /* check wheher there already has been a dof assigned to partnernode */
                  dof = partnernode->dof[l];
                  if (dof==-2 && actnode->dof[l]==-2)
                  {
                     actnode->dof[l] = counter;
                     partnernode->dof[l] = counter;
                     counter++;
                  }
                  else if (dof!=-2 && actnode->dof[l]==-2)
                  {
                     actnode->dof[l] = dof;
                  }
                  else if (dof==-2 && actnode->dof[l]!=-2)
                  {
                     partnernode->dof[l] = actnode->dof[l];
                  }
               }
             /*  else if (dirich==1 && couple==0 && geocouple==1) 
               dserror("Case dirichlet condition in geocoupleset not yet implemented");*/
            }/* end of loops over dofs */
         }/* end of has dirich and/or coupling condition */
      }
   } /* end of loop over nodes */
   /* Now all free dofs are numbered, so now number the dirichlet conditioned
      dofs from here on */
   actfield->dis[0].numeq = counter;
   for (j=0; j<actfield->dis[0].numnp; j++)
   {
      actnode = &(actfield->dis[0].node[j]);
      if (actnode->gnode->couple==NULL && actnode->gnode->dirich==NULL) continue;
      for (l=0; l<actnode->numdf; l++)
      {
         if (actnode->dof[l]==-1)
         {
            actnode->dof[l] = counter;
            counter++;
         }
      }
   }
   actfield->dis[0].numdf = counter;   
} /* end of loop over fields */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of assign_dof */




#include "../headers/standardtypes.h"
#include "../fluid2_pro/fluid2pro.h"
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
 |  put dofs to nodes                                    genk 08/02     |
 |  the numbering of the dofs is done the same on all procs             |
 |  because it seems practical that all                                 |
 |  procs know the dof numbering of all fields.                         |
 |  (there is also no communication at all in this routine)             |
 |  for severel discretisations this seems to get a little bit more     |
 |  more complicated and has to be specified for the different cases    |
 |  the vector nodelfag[] is introduced                                 |
 |  this flag tells us how we have to assign the dofs to the nodes      |
 |    nodeflag[j]=11: node j belongs to a FLUID2_PRO element            |
 |                    actual discretisation: 1                          |
 |   continue with                                                      |
 |     nodeflag[j]=20+kk: node j belongs to a ??? element               |
 |                    actual discretisation: kk                         | 
 |                                                                      |
 *----------------------------------------------------------------------*/
void assign_dof_ndis(FIELD *actfield)
{
int i,j,k,l,kk;                    /* some counters */
int counter;                       
int cpro;
int coupleID;                      /* Id of a coupling set from gid */
int dof;                           /* dof in progress */
int couple,geocouple,dirich;       /* flags for conditions */
ELEMENT *actele;                   /* the element in progress */
NODE    *actnode, *partnernode;    /* node and coupling partner in progress */
int max;
ARRAY   nodeflag_a;
int *nodeflag;
#ifdef DEBUG 
dstrc_enter("assign_dof_ndis");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------ loop discretisations */
for (kk=0;kk<actfield->ndis;kk++)
{
   nodeflag = amdef("nodeflag",&nodeflag_a,actfield->dis[kk].numnp,1,"IV");
   amzero(&nodeflag_a);
   cpro=0;
/*--------------------------------------------- find new node numbering */
/* notice: bandwith opimization has been done by gid, so just do a
   field-local renumbering of this */
   for (j=0; j<actfield->dis[kk].numnp; j++) actfield->dis[kk].node[j].Id_loc = j;

/* loop all elements and put numdf to their nodes depending on kind of
   element. Notice that nodes which couple different kind of elements have
   the higher number of dofs */
   for (j=0; j<actfield->dis[kk].numele; j++)
   {
      actele = &(actfield->dis[kk].element[j]);
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
      case el_fluid2_pro:
        if(actele->e.f2pro->dm==dm_q2q1)
	 {
	    cpro++;
            for (k=0; k<actele->numnp; k++)
            {           
	       nodeflag[actele->node[k]->Id_loc]=10+kk;
	       if (kk==0)
	       {
                  if (actele->node[k]->numdf < 2) actele->node[k]->numdf=2;
               }
	       if (kk==1)
	       {
	          if (actele->node[k]->numdf < 1) actele->node[k]->numdf=1;
	       }
	    } /* end of loop over numnp */
         }
	 else	 
	    dserror("assign_dof for DISMODE of FLUID2_PRO not implemented yet!");	
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
      } /* end switch (actele->eltyp) */
   } /* end of loop over elements */
/*------------------------------------------------------ do some checke */
   if (cpro!=actfield->dis[kk].numele)
      dserror("FLUID2_PRO elements can not be mixed with other fluid elements!\n");
/*---------------------------------------- assign the dofs to the nodes */
   counter=0;
   for (j=0; j<actfield->dis[kk].numnp; j++)
   {
      actfield->dis[kk].node[j].dof  = (int*)CALLOC(actfield->dis[kk].node[j].numdf,sizeof(int));
      if (!(actfield->dis[kk].node[j].dof)) 
         dserror("Allocation of dof in NODE failed");
      /*------------------------- allocate the arrays to hold solutions */
      max = IMAX(3,actfield->dis[kk].node[j].numdf);
      
      amdef("sol",&(actfield->dis[kk].node[j].sol),1,max,"DA");
      amzero(&(actfield->dis[kk].node[j].sol));

      amdef("sol_incr",&(actfield->dis[kk].node[j].sol_increment),1,max,"DA");
      amzero(&(actfield->dis[kk].node[j].sol_increment));

      amdef("sol_res",&(actfield->dis[kk].node[j].sol_residual),1,max,"DA");
      amzero(&(actfield->dis[kk].node[j].sol_residual));
      /*------------------------------------------- init all dofs to -2 */
      for (l=0; l<actfield->dis[kk].node[j].numdf; l++) actfield->dis[kk].node[j].dof[l]=-2;
   }   
/*------- eliminate geostationary coupling conditions that conflict with 
   dofcoupling sets by putting the geostat coupling to the coupling set */
   coupleID=0;
   for (j=0; j<actfield->dis[kk].numnp; j++)
   {
      actnode = &(actfield->dis[kk].node[j]);
      if (actnode->gnode->couple==NULL) continue;
      if (nodeflag[actnode->Id_loc]==10 || nodeflag[actnode->Id_loc]==11)
	   dserror("geostationary coupling not implemented yet for FLUID2_PRO\n");
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
/* this is getting a little bit more complicated for several 
   discretisations. The problem is that the dirichlet conditions from the
   input file may not fit any more to the actual nodes of this discretisation.
   e.g. element FLUID2_PRO:
      - the nodes of the second discretisation will have only one dof (pressure)
      - however we still need the dirichlet information from the input-file
      - but the only dof of the second discretisation corresponds to the third
        dof in gnode->dirich
*/      
   for (j=0; j<actfield->dis[kk].numnp; j++)
   {
      actnode = &(actfield->dis[kk].node[j]);
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
	    switch(nodeflag[actnode->Id_loc])
	    {
	    case 11:
	       f2pro_ass_dof_q2q1(actnode,&counter,0);
	    break;
	    default:
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
	    }/* end switch(nodeflag[actnode->Id]) */
         }/* end of has dirich and/or coupling condition */
      }
   } /* end of loop over nodes */
   /* Now all free dofs are numbered, so now number the dirichlet conditioned
      dofs from here on */
   actfield->dis[kk].numeq = counter;
   for (j=0; j<actfield->dis[kk].numnp; j++)
   {
      actnode = &(actfield->dis[kk].node[j]);
      if (actnode->gnode->couple==NULL && actnode->gnode->dirich==NULL) continue;
      for (l=0; l<actnode->numdf; l++)
      {
         if (actnode->dof[l]==-1)
         {
            actnode->dof[l] = counter;
            counter++;
         }
      } /* end of loop over dofs */
   } /* end of loop over nodes */
   actfield->dis[kk].numdf = counter;   
   amdel(&nodeflag_a);
} /* end of loop over fields */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of assign_dof */




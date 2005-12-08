/*!----------------------------------------------------------------------
\file
\brief coupling conditions for multifield fluid problems

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#ifdef D_FSI
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fluid_prototypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*!---------------------------------------------------------------------
\brief create fluid multifield coupling conditions

<pre>                                                         genk 01/03

In this routine the fluid multifield coupling conditions are created.
Problems with freesurfaces have to be solved in an ALE-framework.
 create mulitfield solution history and set pointers to the corresponding
 nodes of the other fields

</pre>

\param  *fluidfield  FIELD    (i)
\param  *alefield    FIELD    (i)

\return void

------------------------------------------------------------------------*/
void fluid_initmfcoupling(
                           FIELD         *fluidfield,
			   FIELD         *alefield
		         )
{
INT     numfnp,  numanp;                     /* number of nodes         */
INT     numfld;
INT     i,j;                                 /* simply some counters    */
INT     ierr;                                /* flag                    */
INT     afound;                              /* flag                    */
INT     numff,numaf;
DOUBLE  tol=EPS8;                            /* tolerance for node dist */
NODE   *actfnode=NULL, *actanode=NULL;       /* actual nodes            */
GNODE  *actfgnode=NULL,*actagnode=NULL;      /* actual gnodes           */
ARRAY   aindex_a;
INT    *aindex;
DLINE  *actdline=NULL;
DNODE  *actdnode;

#ifdef DEBUG
dstrc_enter("fluid_initmfcoupling");
#endif

/*--------------------------------------------- actual number of fields */
numfld=genprob.numfld;
numff =genprob.numff;
numaf =genprob.numaf;

/*---------------------------- find number of nodes in different fields */
numfnp  = fluidfield->dis[0].numnp;
numanp  = alefield->dis[0].numnp;

/*---------------------------------------------------- loop fluid nodes */
/* multifield solution history of fluid nodes:
   actfnode->sol_mf.a.da[0][i]: velocity transfered to ale
  ----------------------------------------------------------------------*/
for (i=0;i<numfnp;i++)
{
   actfnode  = &(fluidfield->dis[0].node[i]);
   actfgnode = actfnode->gnode;
   amdef("sol_mf",&actfnode->sol_mf,1,actfnode->numdf,"DA");
   amzero(&(actfnode->sol_mf));
   actfgnode->mfcpnode=(NODE**)CCACALLOC(numfld,sizeof(NODE*));
   for (j=0;j<numfld;j++) actfgnode->mfcpnode[j]=NULL;
} /* end of loop over fluid nodes */

/*------------------------------------------------------ loop ale nodes */
/* multifield solution history of ale nodes:
   actanode->sol_mf.a.da[0][i]: displacements at (n)
   actanode->sol_mf.a.da[1][i]: displacements at (n+1)
  ----------------------------------------------------------------------*/
for (i=0;i<numanp;i++)
{
   actanode  = &(alefield->dis[0].node[i]);
   actagnode = actanode->gnode;
   amdef("sol_mf",&actanode->sol_mf,3,actanode->numdf,"DA");
   amzero(&(actanode->sol_mf));
   actagnode->mfcpnode=(NODE**)CCACALLOC(numfld,sizeof(NODE*));
   for (j=0;j<numfld;j++) actagnode->mfcpnode[j]=NULL;
} /* end of loop over ale nodes */

/*------------------------------- create and initialsise index arrays */
aindex = amdef("aindex",&aindex_a,numanp,1,"IV");
for (i=0;i<numanp;i++) aindex[i]=1;

/* find multifield coupled nodes and set ptrs to the corresponding
   nodes of  the other fields ------------------------------------------*/
/*--------------------------------------------------- loop fluid nodes  */
for (i=0;i<numfnp;i++)
{
   actfnode  = &(fluidfield->dis[0].node[i]);
   actfgnode = actfnode->gnode;
   afound=0;
   /*----------------------- loop ale nodes and find corresponding node */
   for (j=0;j<numanp;j++)
   {
      if (aindex[j]==0) continue;
      actanode  = &(alefield->dis[0].node[j]);
      cheque_distance(&(actfnode->x[0]),&(actanode->x[0]),tol,&ierr);
      if(ierr==0) continue;
      afound++;
      actagnode = actanode->gnode;
      aindex[j]=0;
      break;
   } /* end of loop over ale nodes */

   /*------------------------------ set pointers to corresponding nodes */
   actfgnode->mfcpnode[numff]=actfnode;
   if(afound>0)
   {
      actfgnode->mfcpnode[numaf]=actanode;
      actagnode->mfcpnode[numff]=actfnode;
      actagnode->mfcpnode[numaf]=actanode;
   }
}/* end of loop over fluidnodes */

/*------------------------------------------------- plausibility checks */
for (i=0;i<numanp;i++)
{
   actanode  = &(alefield->dis[0].node[i]);
   actagnode = actanode->gnode;
   dsassert(actagnode->mfcpnode!=NULL,"No fluid node for ale node!\n");
   if (actagnode->freesurf==NULL) continue;
   dsassert(actagnode->dirich!=NULL,"No dirich condition for freesurf ale node!\n");
   dsassert(actagnode->dirich->dirich_type==dirich_freesurf,
   "wrong dirch_type at freesurface for ale node!\n");
   if (actanode->locsysId>0)
   {
      actfnode = actagnode->mfcpnode[genprob.numff];
      dsassert(actfnode!=NULL,"cannot read from NULL pointer\n");
      if (actfnode->locsysId != actanode->locsysId)
         dserror("locsysId at free surface not the same for ale and fluid!\n");
   }
}

/*-------------------------------------------------- print out coupling */
if (genprob.visual==0)
out_fluidmf(fluidfield);

/*--------------------------------------------- init the slip dirich BC */
for (i=0;i<numfnp;i++)
{
   actfnode  = &(fluidfield->dis[0].node[i]);
   actfgnode = actfnode->gnode;
   /*-------------------------------- we are looking for a design nodes */
   if (actfgnode->ondesigntyp==ondnode)
   {
      actdnode=actfgnode->d.dnode;
      /*-------------------------------------- find dline with slip DBC */
      ierr=0;
      for(j=0;j<actdnode->ndline;j++)
      {
         actdline=actdnode->dline[j];
         if(actdline->slipdirich!=NULL)
         {
            ierr++;
            break;
         }
      }
      if (ierr==0) continue;
      /*------ check if actual node is first or last node on this dline */
      actdnode=actdline->dnode[0];
      cheque_distance(&(actdnode->x[0]),&(actfnode->x[0]),tol,&ierr);
      if (ierr!=0)
      {
         if (actfgnode->freesurf!=NULL)
         {
            dsassert(actfgnode->slipdirich->firstnode==NULL,
               "lastnode already set for slipdirich!\n");
            actdline->slipdirich->firstnode=actfnode;
         }
         else
         {
            dsassert(actfgnode->slipdirich->lastnode==NULL,
               "lastnode already set for slipdirich!\n");
            actdline->slipdirich->lastnode=actfnode;
         }
         continue;
      }
      actdnode=actdline->dnode[1];
      cheque_distance(&(actdnode->x[0]),&(actfnode->x[0]),tol,&ierr);
      if (ierr!=0)
      {
         if (actfgnode->freesurf!=NULL)
         {
            dsassert(actfgnode->slipdirich->firstnode==NULL,
             "lastnode already set for slipdirich!\n");
            actdline->slipdirich->firstnode=actfnode;
         }
         else
         {
            dsassert(actfgnode->slipdirich->lastnode==NULL,
               "lastnode already set for slipdirich!\n");
            actdline->slipdirich->lastnode=actfnode;
         }
         continue;
      }
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_initmfcoupling*/
#endif
#endif
/*! @} (documentation module close)*/

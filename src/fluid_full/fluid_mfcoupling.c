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
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
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
NODE   *actfnode, *actanode;                 /* actual nodes            */
GNODE  *actfgnode,*actagnode;                /* actual gnodes           */

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
   if (actfgnode->mfcpnode==NULL) 
      dserror("Allocation of coupling node pointers failed");
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
   amdef("sol_mf",&actanode->sol_mf,2,actanode->numdf,"DA");
   amzero(&(actanode->sol_mf));
   actagnode->mfcpnode=(NODE**)CCACALLOC(numfld,sizeof(NODE*));
   if (actagnode->mfcpnode==NULL) 
      dserror("Allocation of coupling node pointers failed");
   for (j=0;j<numfld;j++) actagnode->mfcpnode[j]=NULL;
} /* end of loop over ale nodes */

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
      actanode  = &(alefield->dis[0].node[j]);
      cheque_distance(&(actfnode->x[0]),&(actanode->x[0]),tol,&ierr);
      if(ierr==0) continue;
      afound++;
      actagnode = actanode->gnode;
      break;          
   } /* end of loop over ale nodes */      
   
   /*------------------------------ set pointers to corresponding nodes */
                 actfgnode->mfcpnode[numff]=actfnode;
   if(afound>0)  actfgnode->mfcpnode[numaf]=actanode;   
   if(afound>0)  actagnode->mfcpnode[numff]=actfnode;
   if(afound>0)  actagnode->mfcpnode[numaf]=actanode;	    
}/* end of loop over fluidnodes */   
   
/*------------------------------------------------- plausibility checks */
for (i=0;i<numanp;i++)
{
   actanode  = &(alefield->dis[0].node[i]); 
   actagnode = actanode->gnode;
   if (actagnode->freesurf==NULL) continue; 
   dsassert(actagnode->dirich!=NULL,"No dirich condition for freesurf ale node!\n"); 
   dsassert(actagnode->dirich->dirich_type==dirich_freesurf,
   "wrong dirch_type at freesurface for ale node!\n");      
}   

/*-------------------------------------------------- print out coupling */
if (genprob.visual==0)
out_fluidmf(fluidfield);   

   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_initmfcoupling*/   
#endif   
#endif   	     
/*! @} (documentation module close)*/

/*!----------------------------------------------------------------------
\file
\brief coupling conditions for fsi problems

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*! 
\addtogroup FSI
*//*! @{ (documentation module open)*/
#ifdef D_FSI
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "fsi_prototypes.h"    
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;
/*!---------------------------------------------------------------------
\brief create fsi coupling conditions

<pre>                                                         genk 09/02

In this routine the fsi-coupling conditions from the input file are
evaluated and transformed into Neumann conditions for structure
and into Dirichlet conditions for fluid and ale.
			     
</pre>   

\return void 

------------------------------------------------------------------------*/
void fsi_createfsicoup()
{
INT i;                                    /* simply some counters       */
INT hasdirich,hascouple,hasfsi,hasneum;   /* different flags            */

DLINE    *actdline;
FIELDTYP  fieldtyp;

#ifdef DEBUG 
dstrc_enter("fsi_createfsicoup");
#endif

/*--------------------------------------------------------- loop dlines */
for (i=0; i<design->ndline; i++)
{
   hasdirich=0;
   hascouple=0;
   hasfsi   =0;
   hasneum  =0;
   actdline = &(design->dline[i]);
   /*--------------------------------------------- check for conditions */
   if (actdline->dirich!=NULL) hasdirich++;
   if (actdline->couple!=NULL) hascouple++;
   if (actdline->fsicouple!=NULL) hasfsi++;
   if (actdline->neum!=NULL) hasneum++;
   if (hasfsi==0) continue;
   fieldtyp=actdline->fsicouple->fieldtyp;
   switch (fieldtyp)
   {
   case structure:
      dsassert(hasneum==0,"neumann- and fsi-coupling condition defined on same DLINE\n");
      dsassert(hasdirich==0,"dirich- and fsi-coupling condition defined on same DLINE\n");
      dsassert(hascouple==0,"coupling- and fsi-coupling condition defined on same DLINE\n");      
      actdline->neum = (NEUM_CONDITION*)CCACALLOC(1,sizeof(NEUM_CONDITION));
      if (!actdline->neum) dserror("Allocation of memory failed");
      actdline->neum->neum_type = neum_FSI;
   break;
   case fluid:
      dsassert(hasdirich==0,"dirich- and fsi-coupling condition defined on same DLINE\n");
      dsassert(hascouple==0,"coupling- and fsi-coupling condition defined on same DLINE\n");
      dsassert(hasneum==0,"neumann- and fsi-coupling condition defined on same DLINE\n");
      /*----------- allocate space for a dirichlet condition in this dline */
      actdline->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
      if (!actdline->dirich) dserror("Allocation of memory failed");  
      amdef("onoff",&(actdline->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
      amzero(&(actdline->dirich->dirich_onoff));
      amdef("val",&(actdline->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
      amdef("curve",&(actdline->dirich->curve),MAXDOFPERNODE,1,"IV"); 
      amzero(&(actdline->dirich->dirich_val));
      amzero(&(actdline->dirich->curve));
      amdef("function",&(actdline->dirich->funct),MAXDOFPERNODE,1,"IV");
      amzero(&(actdline->dirich->funct));
      /*----------------------------------- initialise for fsi-coupling */
      actdline->dirich->dirich_onoff.a.iv[0] = 1;   
      actdline->dirich->dirich_onoff.a.iv[1] = 1;       
      actdline->dirich->dirich_type=dirich_FSI;
   break;   
   case ale:
       dsassert(hasdirich==0,"dirich- and fsi-coupling condition defined on same DLINE\n");
       dsassert(hascouple==0,"coupling- and fsi-coupling condition defined on same DLINE\n");
       dsassert(hasneum==0,"neumann- and fsi-coupling condition defined on same DLINE\n");
      /*-------- allocate space for a dirichlet condition in this dline */
      actdline->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
      if (!actdline->dirich) dserror("Allocation of memory failed");  
      amdef("onoff",&(actdline->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
      amzero(&(actdline->dirich->dirich_onoff));
      amdef("val",&(actdline->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
      amdef("curve",&(actdline->dirich->curve),MAXDOFPERNODE,1,"IV"); 
      amzero(&(actdline->dirich->dirich_val));
      amzero(&(actdline->dirich->curve));
      amdef("function",&(actdline->dirich->funct),MAXDOFPERNODE,1,"IV");
      amzero(&(actdline->dirich->funct));
      /*----------------------------------- initialise for fsi-coupling */
      actdline->dirich->dirich_onoff.a.iv[0] = 1;   
      actdline->dirich->dirich_onoff.a.iv[1] = 1;   
      actdline->dirich->dirich_type=dirich_FSI;
   break;
   default:
      dserror("fieldtyp unknown!\n");
   } /* end switch(fieldtyp) */
} /* end of loop over dlines */


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fsi_creatcoup*/

/*!--------------------------------------------------------------------- 
\brief initialise fsi coupling conditions

<pre>                                                         genk 09/02

 create mulitfield solution history and set pointers to the corresponding
 nodes of the other fields
			     
</pre>   

\param *structfield   FIELD         (i)      structure field 
\param *fluidfield    FIELD         (i)      fluid field     
\param *alefield      FIELD         (i)      ale field       

\return void 

------------------------------------------------------------------------*/
void fsi_initcoupling( 
                          FIELD       *structfield,
                          FIELD       *fluidfield, 
		          FIELD       *alefield
		      )
{
INT     numfnp, numsnp, numanp;              /* number of nodes         */  
INT     numdf;                               /* number of dofs          */
INT     numc;                                /* number of columns in mf */
INT     numaf,numsf,numff;
INT     i,j;                                 /* simply some counters    */
INT     ierr;                                /* flag                    */
INT     sfound,afound;                       /* flag                    */          
DOUBLE  tol=EPS8;                            /* tolerance for node dist */
NODE   *actfnode, *actsnode, *actanode;      /* actual nodes            */
GNODE  *actfgnode,*actsgnode,*actagnode;     /* actual gnodes           */

#ifdef DEBUG 
dstrc_enter("fsi_initcoupling");
#endif

/*---------------------------- find number of nodes in different fields */
numsnp  = structfield->dis[0].numnp;
numfnp  = fluidfield->dis[0].numnp;
numanp  = alefield->dis[0].numnp;
numaf   = genprob.numaf;
numff   = genprob.numff;
numsf   = genprob.numsf;

/*---------------------- allocate space for mulitfield solution history *
  the following data are necassary:
  ale needs displacements of the structural fsi coupling nodes
  fluid needs velocity of all ale nodes
  structure needs stresses of the fluid fsi coupling nodes
 *----------------------------------------------------------------------*/

/*---------------------------------------------------- loop fluid nodes */
/* multifield solution history of fluid nodes:
   actfnode->sol_mf.a.da[0][i]: velocity transfered to ale
   actfnode->sol_mf.a.da[1][i]: stresses transfered to structure
   2D: SIGMA11,SIGMA22,SIGMA12
   3D: SIGMA11,SIGMA22,SIGMA33,SIGMA12,SIGMA13,SIGMA23                      
  ----------------------------------------------------------------------*/   
for (i=0;i<numfnp;i++)
{
   actfnode  = &(fluidfield->dis[0].node[i]);
   actfgnode = actfnode->gnode;
   numdf = actfnode->numdf;
   if (numdf==3 || numdf==5) numc=numdf;
   else if (numdf==4) numc=IMAX(6,numdf);
   else  dserror("number of fluid dofs not possible!\n");
   amdef("sol_mf",&actfnode->sol_mf,2,numc,"DA");
   amzero(&(actfnode->sol_mf));
   actfgnode->mfcpnode=(NODE**)CCACALLOC(3,sizeof(NODE*));
   if (actfgnode->mfcpnode==NULL) 
      dserror("Allocation of coupling node pointers failed");
   for (j=0;j<3;j++) actfgnode->mfcpnode[j]=NULL;
} /* end of loop over fluid nodes */

/*------------------------------------------------ loop structure nodes */
/* multifield solution history of structural nodes:
   actsnode->sol_mf.a.da[0][i]: actual displacements transfered to ale
   actsnode->sol_mf.a.da[1][i]: displacements of old iteration step 
   actsnode->sol_mf.a.da[2][i]: displacements of old time step
   actsnode->sol_mf.a.da[3][i]: dispi
   actsnode->sol_mf.a.da[4][i]: couplingforces at the end of time step
   actsnode->sol_mf.a.da[5][i]: coupling forces at the beginning of time step
 -----------------------------------------------------------------------*/   
for (i=0;i<numsnp;i++)
{
   actsnode  = &(structfield->dis[0].node[i]);
   actsgnode = actsnode->gnode;
   numdf = actsnode->numdf;
   amdef("sol_mf",&actsnode->sol_mf,6,numdf,"DA");
   amzero(&(actsnode->sol_mf));
   actsgnode->mfcpnode=(NODE**)CCACALLOC(3,sizeof(NODE*));
   if (actsgnode->mfcpnode==NULL) 
      dserror("Allocation of coupling node pointers failed");
   for (j=0;j<3;j++) actsgnode->mfcpnode[j]=NULL;
} /* end of loop over struct nodes */

/*------------------------------------------------------ loop ale nodes */
/* multifield solution history of ale nodes:
   actanode->sol_mf.a.da[0][i]: displacements at (n)
   actanode->sol_mf.a.da[1][i]: displacements at (n+1)
  ----------------------------------------------------------------------*/
for (i=0;i<numanp;i++)
{
   actanode  = &(alefield->dis[0].node[i]);
   actagnode = actanode->gnode;
   numdf = actanode->numdf;
   amdef("sol_mf",&actanode->sol_mf,2,numdf,"DA");
   amzero(&(actanode->sol_mf));
   actagnode->mfcpnode=(NODE**)CCACALLOC(3,sizeof(NODE*));
   if (actagnode->mfcpnode==NULL) 
      dserror("Allocation of coupling node pointers failed");
   for (j=0;j<3;j++) actagnode->mfcpnode[j]=NULL;
} /* end of loop over ale nodes */

/* find fsi coupled nodes and set ptrs to the corresponding nodes of 
   the other fields --------------------------------------------------*/
/*------------------------------------------------- loop fluid nodes  */    
for (i=0;i<numfnp;i++)
{
   actfnode  = &(fluidfield->dis[0].node[i]);
   actfgnode = actfnode->gnode;
   sfound=0;
   afound=0;
   /*------------------ loop struct nodes and find corresponding node */
   for (j=0;j<numsnp;j++)
   {
      actsnode   = &(structfield->dis[0].node[j]);
      cheque_distance(&(actfnode->x[0]),&(actsnode->x[0]),tol,&ierr);
      if(ierr==0) continue;
      sfound++;
      actsgnode = actsnode->gnode;
      break;      
   } /* end of loop over structnodes */  
   /*--------------------- loop ale nodes and find corresponding node */
   for (j=0;j<numanp;j++)
   {
      actanode  = &(alefield->dis[0].node[j]);
      cheque_distance(&(actfnode->x[0]),&(actanode->x[0]),tol,&ierr);
      if(ierr==0) continue;
      afound++;
      actagnode = actanode->gnode;
      break;          
   } /* end of loop over ale nodes */      

   /*--------------------------- set pointers to corresponding nodes */
   if(sfound>0)   actfgnode->mfcpnode[numsf]=actsnode;
                  actfgnode->mfcpnode[numff]=actfnode;
   if(afound>0)   actfgnode->mfcpnode[numaf]=actanode; 
   if(sfound>0)   actsgnode->mfcpnode[numsf]=actsnode;
   if(sfound>0)   actsgnode->mfcpnode[numff]=actfnode;
   if(sfound>0)   actsgnode->mfcpnode[numaf]=actanode;    
   if(sfound>0)   actagnode->mfcpnode[numsf]=actsnode;
   if(afound>0)   actagnode->mfcpnode[numff]=actfnode;
   if(afound>0)   actagnode->mfcpnode[numaf]=actanode;  	   
}/* end of loop over fluidnodes */

/*------------------------------------------------ plausibility checks */
for (i=0;i<numfnp;i++)
{
   actfnode  = &(fluidfield->dis[0].node[i]);
   actfgnode = actfnode->gnode;
   if (actfgnode->fsicouple==NULL) continue;
   actsnode = actfgnode->mfcpnode[numsf];
   actanode = actfgnode->mfcpnode[numaf];
   actsgnode = actsnode->gnode;
   actagnode = actanode->gnode;
   dsassert(actsgnode->fsicouple!=NULL,"FSI Coupling Condition Fluid-Struct not the same!\n");
   dsassert(actagnode->fsicouple!=NULL,"FSI Coupling Condition Fluid-Ale not the same!\n");
   dsassert(actfgnode->fsicouple->fsi_coupleId==actsgnode->fsicouple->fsi_coupleId,
            "FSI Coupling Condition Fluid-Struct: wrong coupleId\n");
   dsassert(actfgnode->fsicouple->fsi_coupleId==actagnode->fsicouple->fsi_coupleId,
            "FSI Coupling Condition Fluid-Ale: wrong coupleId\n");	
   dsassert(actfgnode->fsicouple->fsi_mesh==actsgnode->fsicouple->fsi_mesh,
            "FSI Coupling Condition Fluid-Struct: wrong mesh\n"); 
   dsassert(actfgnode->fsicouple->fsi_mesh==actagnode->fsicouple->fsi_mesh,
            "FSI Coupling Condition Fluid-Ale: wrong mesh\n");
   dsassert(actfgnode->dirich!=NULL,"No dirich condition for fsi-coupling fluid node!\n");
   if(actfgnode->dirich->dirich_type!=dirich_FSI) continue;
   dsassert(actfgnode->dirich->dirich_onoff.a.iv[0]==1,"wrong onoff(0) at fluid node!\n");
   dsassert(actfgnode->dirich->dirich_onoff.a.iv[1]==1,"wrong onoff(1) at fluid node!\n");
   for (j=2;j<MAXDOFPERNODE;j++)
   dsassert(actfgnode->dirich->dirich_onoff.a.iv[j]==0,"wrong onoff() at fluid node!\n");
} /* end of loop over fluid nodes */

for (i=0;i<numanp;i++)
{
   actanode  = &(alefield->dis[0].node[i]);
   actagnode = actanode->gnode;
   if (actagnode->fsicouple==NULL) continue; 
   dsassert(actagnode->dirich!=NULL,"No dirich condition for fsi-coupling ale node!\n"); 
   if(actagnode->dirich->dirich_type!=dirich_FSI) continue;
   dsassert(actagnode->dirich->dirich_onoff.a.iv[0]==1,"wrong onoff(0) at ale node!\n");
   dsassert(actagnode->dirich->dirich_onoff.a.iv[1]==1,"wrong onoff(1) at ale node!\n");
   for (j=2;j<MAXDOFPERNODE;j++)
   dsassert(actagnode->dirich->dirich_onoff.a.iv[j]==0,"wrong onoff() at ale node!\n");       
}
/*------------------------------------------------ print out coupling */
if (genprob.visual==0)
out_fsi(fluidfield);

/*--------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fsi_initcoup*/

/*!---------------------------------------------------------------------
\brief determine structural fsi interface dofs

<pre>                                                         genk 01/03

 create array sid (structural interface dofs)
			     
</pre>   

\param *structfield   FIELD         (i)      structure field 
\param *fsidyn        FSI_DYNAMIC   (i)               
\return void 

------------------------------------------------------------------------*/
void fsi_struct_intdofs(
                          FIELD       *structfield, 
			  FSI_DYNAMIC *fsidyn
		       )
{
INT    i,j;                    /* some counters                         */
INT    numnp_total;            /* total number of structure dofs        */             
INT    dof;                    /* actual dof                            */
INT    counter=0;              
INT    numaf;
INT   *sid;                    /* structural interface dofs             */
NODE  *actsnode, *actanode;    /* actual nodes                          */
GNODE *actsgnode, *actagnode;  /* actual gnodes                         */

#ifdef DEBUG 
dstrc_enter("fsi_struct_intdofs");
#endif

numnp_total = structfield->dis[0].numnp;
numaf       = genprob.numaf;

sid = amdef("sid",&fsidyn->sid,structfield->dis[0].numdf,1,"IV");
amzero(&fsidyn->sid);


/*---------------------------------------------------------- loop nodes */
for (i=0;i<numnp_total;i++)
{
   actsnode  = &(structfield->dis[0].node[i]);
   actsgnode = actsnode->gnode;
   actanode  = actsgnode->mfcpnode[numaf];
   if (actanode==NULL) continue; 
   actagnode = actanode->gnode;
   /*----------------------------------------- check for coupling nodes */
   dsassert(actagnode->dirich!=NULL,"no dirich condition for coupled ALE node\n");
   if (actagnode->dirich->dirich_type!=dirich_FSI) continue;
   for (j=0;j<actanode->numdf;j++)
   {
      dof = actsnode->dof[j];      
      sid[dof]=1;
      counter++;
   }      
} /* end of loop over nodes */

fsidyn->numsid=counter;

#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fsi_struct_intdofs*/
#endif

/*! @} (documentation module close)*/

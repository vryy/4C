/*!----------------------------------------------------------------------
\file
\brief coupling conditions for ssi problems

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
#ifndef CCADISCRET
/*!
\addtogroup SSI
*//*! @{ (documentation module open)*/
#ifdef D_SSI
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "ssi_prototypes.h"
#include "../wall1/wall1.h"
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
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

/*!---------------------------------------------------------------------
\brief create ssi coupling conditions

<pre>                                                         genk 09/02

In this routine the ssi-coupling conditions from the input file are
evaluated and transformed into Neumann/Dirichlet conditions

ssi_couptyp=ssi_master -> Neumann condition
ssi_couptyp=ssi_slave  -> Dirichlet condition

</pre>

\return void

------------------------------------------------------------------------*/
void ssi_createssicoup()
{
INT i,j;                                  /* simply some counters       */
INT hasdirich,hascouple,hasssi,hasneum;   /* different flags            */
DLINE    *actdline;
DNODE    *actdnode;

SSI_COUPTYP ssi_line_couptyp;

#ifdef DEBUG
dstrc_enter("ssi_createssicoup");
#endif

/*---------------------------------------------------------- loop DLINE */
/*----------- inherit ssi coupling condition to all dnodes on dlines ---*/
for (i=0; i<design->ndline; i++)
{
   actdline = &(design->dline[i]);
   /*------------------------- do nothing if DLINE has no ssi condition */
   if (actdline->ssicouple == NULL) continue;
   if (actdline->ssicouple->ssi_couptyp==ssi_none) continue;
   /*------------------------------ loop the DNODEs related to actdline */
   for (j=0; j<actdline->ndnode; j++)
   {
      actdnode = actdline->dnode[j];
      /* check if there is a ssicouple condition at the actdnode, is not */
      /* allocate the memory for a ssicouple condition */
      if(actdnode->ssicouple == NULL)
      {
        actdnode->ssicouple = (SSI_COUPLE_CONDITION*)CCACALLOC(1,sizeof(SSI_COUPLE_CONDITION));
        if (!actdnode->ssicouple) dserror("Allocation of memory failed");
      }
      /*----- inherit the dirichlet condition from actdline to actdnode */
      actdnode->ssicouple->ssi_couptyp = actdline->ssicouple->ssi_couptyp;
      actdnode->ssicouple->ssi_coupleId = actdline->ssicouple->ssi_coupleId;
      actdnode->ssicouple->ssi_mesh = actdline->ssicouple->ssi_mesh;
   }/* loop j over dnodes */
}/* loop i over dlines */


/*--------------------------------------------------------- loop dlines */
for (i=0; i<design->ndline; i++)
{
   hasdirich=0;
   hascouple=0;
   hasssi   =0;
   hasneum  =0;
   actdline = &(design->dline[i]);
   /*--------------------------------------------- check for conditions */
   if (actdline->dirich!=NULL) hasdirich++;
   if (actdline->couple!=NULL) hascouple++;
/*   if (actdline->ssicouple->ssi_couptyp!=ssi_none) hasssi++; */
   if (actdline->ssicouple!=NULL) hasssi++;
   if (actdline->neum!=NULL) hasneum++;
   if (hasssi==0) continue;
   ssi_line_couptyp=actdline->ssicouple->ssi_couptyp;
   switch (ssi_line_couptyp)
   {
   case ssi_master:
      dsassert(hasneum==0,"neumann- and ssi-coupling condition defined on same DLINE\n");
      dsassert(hasdirich==0,"dirich- and ssi-coupling condition defined on same DLINE\n");
      dsassert(hascouple==0,"coupling- and ssi-coupling condition defined on same DLINE\n");
      actdline->neum = (NEUM_CONDITION*)CCACALLOC(1,sizeof(NEUM_CONDITION));
      if (!actdline->neum) dserror("Allocation of memory failed");
      actdline->neum->neum_type = neum_SSI;
   break;
   case ssi_slave:
      dsassert(hasdirich==0,"dirich- and ssi-coupling condition defined on same DLINE\n");
      dsassert(hascouple==0,"coupling- and ssi-coupling condition defined on same DLINE\n");
      dsassert(hasneum==0,"neumann- and ssi-coupling condition defined on same DLINE\n");
      /*----------- allocate space for a dirichlet condition in this dline */
      actdline->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
      if (!actdline->dirich) dserror("Allocation of memory failed");
      amdef("onoff",&(actdline->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
      amzero(&(actdline->dirich->dirich_onoff));
      amdef("val",&(actdline->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
      amdef("curve",&(actdline->dirich->curve),MAXDOFPERNODE,1,"IV");
      amdef("function",&(actdline->dirich->funct),MAXDOFPERNODE,1,"IV");
      amzero(&(actdline->dirich->dirich_val));
      amzero(&(actdline->dirich->curve));
      amzero(&(actdline->dirich->funct));
      /*----------------------------------- initialise for ssi-coupling */
      actdline->dirich->dirich_onoff.a.iv[0] = 1;
      actdline->dirich->dirich_onoff.a.iv[1] = 1;
      actdline->dirich->dirich_type=dirich_SSI;
   break;
   default:
      dserror("ssi_couptyp unknown!\n");
   } /* end switch(fieldtyp) */
} /* end of loop over dlines */


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of ssi_createssicoup */


/*!---------------------------------------------------------------------
\brief initialise fsi coupling conditions

<pre>                                                         genk 09/02

 create mulitfield solution history and set pointers to the corresponding
 nodes of the other fields

</pre>

\param *structfield   FIELD         (i)      structure field
\param *fluidfield    FIELD         (i)      fluid field

\return void

------------------------------------------------------------------------*/
void ssi_initcoupling(
                          FIELD       *masterfield,
                          FIELD       *slavefield,
                          INTERFACES  *int_faces
		      )
{
INT     nummnp, numsnp;             /* number of nodes                  */
INT     numdf;                      /* number of dofs                   */
INT     numsf,nummf;
INT     i,j,a;                      /* simply some counters             */
INT     ierr;                       /* flag                             */
INT     sfound;                     /* flag                             */
DOUBLE  tol=EPS8;                   /* tolerance for node dist          */
NODE   *actmnode, *actsnode=NULL;   /* actual nodes                     */
GNODE  *actmgnode,*actsgnode=NULL;  /* actual gnodes                    */
GLINE  *actgline;                   /* actual glines                    */
ARRAY   int_ids;                    /* a vector with the interface Ids  */
DOUBLE  *int_ids_a;
#ifdef DEBUG
dstrc_enter("ssi_initcoupling");
#endif

/*---------------------------- find number of nodes in different fields */
numsnp  = slavefield->dis[0].numnp;
nummnp  = masterfield->dis[0].numnp;
nummf   = 0;
numsf   = 1;

/*---------------------- allocate space for mulitfield solution history *
  the following data are necassary:
  master (NEUMANN) needs displacements transfered to slave (DIRICHLET)
  slave (DIRICHLET) needs forces transfered to master (NEUMANN)
 *----------------------------------------------------------------------*/

/*---------------------------------------------------- loop master nodes */
/* multifield solution history of master nodes:
   actfnode->sol_mf.a.da[0][i]: displacements transfered to slave
  ----------------------------------------------------------------------*/
for (i=0;i<nummnp;i++)
{
   actmnode  = &(masterfield->dis[0].node[i]);
   actmgnode = actmnode->gnode;
   numdf = actmnode->numdf;
   amdef("sol_mf",&actmnode->sol_mf,7,numdf,"DA");
   amzero(&(actmnode->sol_mf));
   amredef(&actmnode->sol,11,numdf,"DA");
   amzero(&(actmnode->sol));
   actmgnode->mfcpnode=(NODE**)CCACALLOC(2,sizeof(NODE*));
   actmgnode->mfcpnode[nummf]=actmnode;
   actmgnode->mfcpnode[numsf]=NULL;
} /* end of loop over fluid nodes */

/*--------------------------------------------------- loop slave nodes */
/* multifield solution history of slave nodes:
   actsnode->sol_mf.a.da[0][i]: couplingforces at the end of time step
   actsnode->sol_mf.a.da[1][i]: coupling forces at the beginning of time step
 -----------------------------------------------------------------------*/
for (i=0;i<numsnp;i++)
{
   actsnode  = &(slavefield->dis[0].node[i]);
   actsgnode = actsnode->gnode;
   numdf = actsnode->numdf;
   amdef("sol_mf",&actsnode->sol_mf,7,numdf,"DA");
   amzero(&(actsnode->sol_mf));
   amredef(&actsnode->sol,11,numdf,"DA");
   amzero(&(actsnode->sol));
   actsgnode->mfcpnode=(NODE**)CCACALLOC(2,sizeof(NODE*));
   actsgnode->mfcpnode[nummf]=NULL;
   actsgnode->mfcpnode[numsf]=NULL;
} /* end of loop over struct nodes */

/* find ssi coupled nodes and set ptrs to the corresponding nodes of
   the other fields --------------------------------------------------*/
/*------------------------------------------------- loop master nodes  */
for (i=0;i<nummnp;i++)
{
   actmnode  = &(masterfield->dis[0].node[i]);
   actmgnode = actmnode->gnode;
   sfound=0;
   /*------------------- loop slave nodes and find corresponding node */
   for (j=0;j<numsnp;j++)
   {
      actsnode   = &(slavefield->dis[0].node[j]);
      cheque_distance(&(actmnode->x[0]),&(actsnode->x[0]),tol,&ierr);
      if(ierr==0) continue;
      sfound++;
      actsgnode = actsnode->gnode;
      break;
   } /* end of loop over structnodes */

   /*--------------------------- set pointers to corresponding nodes */
   if(sfound>0)   actmgnode->mfcpnode[numsf]=actsnode;
   if(sfound>0)   actsgnode->mfcpnode[nummf]=actmnode;
}/* end of loop over fluidnodes */

/*------------------------------------------------ print out coupling */
if (genprob.visual==0)
out_ssi(masterfield);

/* detect the number of glines with a coupling condition */
/* loop glines of master field */
a = 0;
for(i=0; i<masterfield->dis->ngline; i++)
{
  actgline = &(masterfield->dis->gline[i]);
  if(actgline->ssicouple == NULL) continue;
  if(actgline->ssicouple->ssi_couptyp == ssi_master)
    a++;
}
/* store the couplingId's of the glines with a coupling condition  in */
/* the vector int_ids.a.iv */
int_ids_a = amdef("int_ids", &int_ids,a,1,"IV");
amzero(&int_ids);
a=0;
for(i=0; i<masterfield->dis->ngline; i++)
{
  actgline = &(masterfield->dis->gline[i]);
  if(actgline->ssicouple == NULL) continue;
  if(actgline->ssicouple->ssi_couptyp == ssi_master)
  {
    int_ids.a.iv[a] = actgline->ssicouple->ssi_coupleId;
    a++;
  }
}
/* sort the vector int_ids*/
a=0;
for(i=0; i<int_ids.fdim; i++)
{
  for(j=0; j<int_ids.fdim-1; j++)
  {
    if(int_ids.a.iv[j]>int_ids.a.iv[j+1])
    {
      a = int_ids.a.iv[j];
      int_ids.a.iv[j] = int_ids.a.iv[j+1];
      int_ids.a.iv[j+1] = a;
    }
  }
}
/* detect the number of interfaces */
a=1;
for(i=0; i<int_ids.fdim-1; i++)
{
  if(int_ids.a.iv[i] != int_ids.a.iv[i+1])
    a++;
}

/* store the number of interfaces in int_faces->numint */
int_faces->numint = a;

/* allocate memory for the vector of interface Id's */
int_faces->int_ids = (INT*)CCACALLOC(int_faces->numint,
                     sizeof(INT));

/* store the interface Id's in the vector int_faces->int_ids */
/* loop over the vector int_ids.a.iv which contains the interface Id's*/
a=1;
for(i=0; i<int_ids.fdim-1; i++)
{
  if(i==0)
    int_faces->int_ids[i] = int_ids.a.iv[i];
  if(int_ids.a.iv[i] != int_ids.a.iv[i+1] && int_ids.a.iv[i+1] != 0)
  {
    int_faces->int_ids[a] = int_ids.a.iv[i+1];
    a++;
  }
}

/*--------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of ssi_initcoupling*/

/*!---------------------------------------------------------------------
\brief determine structural fsi interface dofs

<pre>                                                         genk 01/03

 create array sid (structural interface dofs)

</pre>

\param *structfield   FIELD         (i)      structure field
\param *ssidyn        SSI_DYNAMIC   (i)
\return void

------------------------------------------------------------------------*/
void ssi_master_intdofs(
                          FIELD       *masterfield,
			  SSI_DYNAMIC *ssidyn
		       )
{
INT    i,j;                    /* some counters                         */
INT    numnp_total;            /* total number of structure dofs        */
INT    dof;                    /* actual dof                            */
INT    counter=0;
INT   *sid;                    /* structural interface dofs             */
NODE  *actmnode;               /* actual nodes                          */
GNODE *actmgnode;              /* actual gnodes                         */

#ifdef DEBUG
dstrc_enter("ssi_master_intdofs");
#endif

numnp_total = masterfield->dis[0].numnp;

sid = amdef("sid",&ssidyn->sid,masterfield->dis[0].numdf,1,"IV");
amzero(&ssidyn->sid);


/*---------------------------------------------------------- loop nodes */
for (i=0;i<numnp_total;i++)
{
   actmnode  = &(masterfield->dis[0].node[i]);
   actmgnode = actmnode->gnode;
   if (actmgnode->ssicouple == NULL) continue;
   if (actmgnode->ssicouple->ssi_couptyp==ssi_none) continue;
   for (j=0;j<actmnode->numdf;j++)
   {
      dof = actmnode->dof[j];
      sid[dof]=1;
      counter++;
   }
} /* end of loop over nodes */

ssidyn->numsid=counter;

#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of ssi_master_intdofs*/
#endif
/*! @} (documentation module close)*/
#endif

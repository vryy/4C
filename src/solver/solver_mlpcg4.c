#ifdef MLPCG

/*!---------------------------------------------------------------------
\file
\brief contains the multilevel preconditioner for shells

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0771 - 685-6572
</pre>

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/prototypes_mlpcg.h"
INT cmp_int(const void *a, const void *b );
DOUBLE cmp_double(const void *a, const void *b );
/*! 
\addtogroup MLPCG 
*//*! @{ (documentation module open)*/
/*!----------------------------------------------------------------------
\brief the multilevel preconditioner main structure

<pre>                                                         m.gee 09/02    
defined in solver_mlpcg.c
</pre>

*----------------------------------------------------------------------*/
extern struct _MLPRECOND mlprecond;



/*!---------------------------------------------------------------------
\brief create the tentative prolongator from actlev+1 to actlev                                              

<pre>                                                        m.gee 9/02 
create the tentative prolongator from actlev+1 to actlev  
from the aggregation done before , this routine only from 0 to 1
</pre>
\param actlev      MLLEVEL*     (i/o) the active level in the ml-precond.                   
\param actintra   INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_P(MLLEVEL  *actlev, INTRA *actintra)
{
INT          i,j,k,counter=0;
INT          myrank,nproc;
DBCSR       *actstiff;
DBCSR       *P;
AGG         *actagg;
MLLEVEL     *prevlev;
INT          nrow,ncol;
DOUBLE       aggblock[1000][500];
INT          rindex[1000],cindex[500];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_P");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
actstiff    = actlev->csr;
/*--------------------------- allocate prolongator, if it doesn't exist */
if (actlev->P == NULL)
{
   actlev->P = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
   P = actlev->P;
/*----------------------------------------------------------------------*/
/* 
the prolongator takes the row blocks from the actstiff  
*/
   am_alloc_copy(&(actstiff->blocks),&(P->blocks));   
/*----------------------------------------------------------------------*/
/* 
make a good guess of the size of the prolongator and init the csr matrix
*/
   counter = actlev->agg[0].numdf * actstiff->numeq;
/*----------------------------------------------------------------------*/
/* 
make the guess 80% too large 
*/
   counter = (INT)(counter*4.0);
/*----------------------------------------------------------------------*/
/* 
open the csr matrix to add to 
*/
   mlpcg_csr_open(P,
                  actstiff->update.a.iv[0],
                  actstiff->update.a.iv[actstiff->update.fdim-1],
                  actstiff->numeq_total,
                  counter,
                  actintra);
/*----------------------------------------------------------------------*/
/*
get first interproc row in P 
*/
   P->firstcoupledof=actstiff->firstcoupledof;
} /* end of if (actlev->P == NULL) */
else
   mlpcg_csr_zero(actlev->P,actintra);
P = actlev->P;
/*-------- get the pointer to the previous level, to get the R matrices */
prevlev = actlev-1;
/*----------------------------------------------------------------------*/
/* 
loop the aggregates and create the tentative prolongator  
*/                  
for (i=0; i<actlev->nagg; i++)
{
   actagg = &(actlev->agg[i]);
   /* create the tentative prolongator diagonal block of this aggregate */
   mlpcg_precond_oneP_vanek(actagg,aggblock,rindex,cindex,&nrow,&ncol,actstiff,prevlev);
   if (!(actagg->tentP))
   {
      actagg->tentP = (ARRAY*)CCACALLOC(1,sizeof(ARRAY));
      amdef("tentP",actagg->tentP,nrow,ncol,"DA");
   }
   else if (actagg->tentP->fdim != nrow) 
       dserror("Size mismatch in aggregate");
   else if (actagg->tentP_nrow != nrow)  
       dserror("Size mismatch in aggregate");
   if (!(actagg->tentP_rindex))
      actagg->tentP_rindex = (INT*)CCAMALLOC(nrow*sizeof(INT));
   /*------------------------------- put rindex to actagg->tentP_rindex */
   if (mlprecond.ncall==0)
      for (j=0; j<nrow; j++) actagg->tentP_rindex[j] = rindex[j];
   /*--------------------------- put the aggblock to actagg->tentP.a.da */
   for (j=0; j<ncol; j++)
   for (k=0; k<nrow; k++)
      actagg->tentP->a.da[k][j] = aggblock[k][j];
   /*------------------ set number of rows in this piece of Prolongator */
   actagg->tentP_nrow = nrow;
} /* end of for (i=0; i<actlev->nagg; i++) */
/*--------------- loop all aggregates again and fill the DBCSR matrix P */   
#if 1
for (i=0; i<actlev->nagg; i++)
{
   for (j=0; j<actlev->agg[i].numdf; j++) /* column loop */
   for (k=0; k<actlev->agg[i].tentP_nrow; k++) /* row loop */
   if (FABS(actlev->agg[i].tentP->a.da[k][j])>EPS15)
   mlpcg_csr_setentry(P,
                      actlev->agg[i].tentP->a.da[k][j],
                      actlev->agg[i].tentP_rindex[k],
                      actlev->agg[i].dof[j],
                      actintra);
}
#endif
#if 0
for (i=0; i<actlev->nagg; i++)
{
   for (j=0; j<actlev->agg[i].numdf; j++) /* column loop */
   for (k=0; k<actlev->agg[i].tentP_nrow; k++) /* row loop */
   aggblock[k][j] = actlev->agg[i].tentP->a.da[k][j];
   mlpcg_csr_setblock(P,
                      aggblock,
                      actlev->agg[i].tentP_rindex,
                      actlev->agg[i].dof,
                      actlev->agg[i].tentP_nrow,
                      actlev->agg[i].numdf,
                      actintra);
}
#endif
/*------------------------------- make the smoothing of the prolongator */
if (mlprecond.omega>0.0)
   mlpcg_smoothP(P,actlev->csr,actintra);
/*------------------------ tent. prolongator is ready, close the matrix */
mlpcg_csr_close(P);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_P */



/*!---------------------------------------------------------------------
\brief create the tentative prolongator from actlev+1 to actlev                                              

<pre>                                                        m.gee 9/02 
create the tentative prolongator from actlev+1 to actlev  
from the aggregation done before , this routine only from 0 to 1
</pre>
\param actlev      MLLEVEL*     (i/o) the active level in the ml-precond.                   
\param actintra   INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_P0(MLLEVEL  *actlev, INTRA *actintra)
{
INT          i,j,k,counter=0;
INT          myrank,nproc;
DBCSR       *actstiff;
DBCSR       *P;
AGG         *actagg;
INT          nrow,ncol;
DOUBLE       aggblock[1000][500];
INT          rindex[1000],cindex[500];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_P0");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
actstiff    = actlev->csr;
/*--------------------------- allocate prolongator, if it doesn't exist */
if (actlev->P == NULL)
   actlev->P = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
P = actlev->P;
/*----------------------------------------------------------------------*/
/* 
the prolongator takes the row blocks from the actstiff  
*/
if (mlprecond.ncall==0)
{
   am_alloc_copy(&(actstiff->blocks),&(P->blocks));   
/*----------------------------------------------------------------------*/
/* 
make a good guess of the size of the prolongator and init the csr matrix
*/
   counter = actlev->agg[0].numdf * actstiff->numeq;
/*----------------------------------------------------------------------*/
/* 
make the guess 80% too large 
*/
   counter = (INT)(counter*4.0);
/*----------------------------------------------------------------------*/
/* 
open the csr matrix to add to 
*/
   mlpcg_csr_open(P,
                  actstiff->update.a.iv[0],
                  actstiff->update.a.iv[actstiff->update.fdim-1],
                  actstiff->numeq_total,
                  counter,
                  actintra);
}
else if (mlprecond.mod==0)
   mlpcg_csr_zero(P,actintra);
/*----------------------------------------------------------------------*/
/*
get first interproc row in P 
*/
if (mlprecond.ncall==0)
   P->firstcoupledof=actstiff->firstcoupledof;
/*----------------------------------------------------------------------*/
/*
get the directors of the shell elements 
*/
if (mlprecond.ncall==0)
   mlpcg_precond_getdirs();
/*----------------------------------------------------------------------*/
/* 
loop the aggregates and create the tentative prolongator  
*/                  
for (i=0; i<actlev->nagg; i++)
{
   actagg = &(actlev->agg[i]);
   /* create the tentative prolongator diagonal block of this aggregate */
   mlpcg_precond_oneP0_vanek(actagg,aggblock,rindex,cindex,&nrow,&ncol,actstiff);
   if (!(actagg->tentP))
   {
      actagg->tentP = (ARRAY*)CCACALLOC(1,sizeof(ARRAY));
      amdef("tentP",actagg->tentP,nrow,ncol,"DA");
   }
   else if (actagg->tentP->fdim != nrow) 
       dserror("Size mismatch in aggregate");
   else if (actagg->tentP_nrow != nrow)  
       dserror("Size mismatch in aggregate");
   if (!(actagg->tentP_rindex))
      actagg->tentP_rindex = (INT*)CCAMALLOC(nrow*sizeof(INT));
   /*------------------------------- put rindex to actagg->tentP_rindex */
   if (mlprecond.ncall==0)
      for (j=0; j<nrow; j++) actagg->tentP_rindex[j] = rindex[j];
   /*--------------------------- put the aggblock to actagg->tentP.a.da */
   for (j=0; j<ncol; j++)
   for (k=0; k<nrow; k++)
      actagg->tentP->a.da[k][j] = aggblock[k][j];
   /*------------------ set number of rows in this piece of Prolongator */
   actagg->tentP_nrow = nrow;
} /* end of for (i=0; i<actlev->nagg; i++) */
/*--------------- loop all aggregates again and fill the DBCSR matrix P */   
#if 1
for (i=0; i<actlev->nagg; i++)
{
   for (j=0; j<actlev->agg[i].numdf; j++) /* column loop */
   for (k=0; k<actlev->agg[i].tentP_nrow; k++) /* row loop */
   if (FABS(actlev->agg[i].tentP->a.da[k][j])>EPS15)
   mlpcg_csr_setentry(P,
                      actlev->agg[i].tentP->a.da[k][j],
                      actlev->agg[i].tentP_rindex[k],
                      actlev->agg[i].dof[j],
                      actintra);
}
#endif
#if 0
for (i=0; i<actlev->nagg; i++)
{
   for (j=0; j<actlev->agg[i].numdf; j++) /* column loop */
   for (k=0; k<actlev->agg[i].tentP_nrow; k++) /* row loop */
   aggblock[k][j] = actlev->agg[i].tentP->a.da[k][j];
   mlpcg_csr_setblock(P,
                      aggblock,
                      actlev->agg[i].tentP_rindex,
                      actlev->agg[i].dof,
                      actlev->agg[i].tentP_nrow,
                      actlev->agg[i].numdf,
                      actintra);
}
#endif
/*------------------------------- make the smoothing of the prolongator */
if (mlprecond.omega>0.0)
   mlpcg_smoothP(P,actlev->csr,actintra);
/*------------------------ tent. prolongator is ready, close the matrix */
mlpcg_csr_close(P);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_P0 */



/*!---------------------------------------------------------------------
\brief create the tentative prolongator for one aggregate                                          

<pre>                                                        m.gee 11/02 

</pre>
\param actagg         AGG*    (i/o) the active aggregate
\param aggblock       DOUBLE[1000][500] (o) the aggregate's block prolongator
\param rindex         INT[1000]         (o) global indizes of aggblock
\param cindex         INT[500]          (o) global indizes of aggblock
\param nrow           INT*              (o) dimension of rindex
\param ncol           INT*              (o) dimension of cindex
\param actstiff       DBCSR*            (i) fine grid stiffness matrix
\param prevlevel      MLLEVEL*          (i) last level
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_oneP_vanek(AGG     *actagg,
                             DOUBLE   aggblock[][500],
                             INT      rindex[],
                             INT      cindex[],
                             INT     *nrow,
                             INT     *ncol,
                             DBCSR   *actstiff,
                             MLLEVEL *prevlevel)
{
INT           i,j,k,counter;
AGG          *prevagg[200];
AGG          *actprevagg;
INT          *actblock;
INT           dof;
DOUBLE        x0,y0,z0,x,y,z;
INT           index;
INT           shift,bins[5000];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_oneP_vanek");
#endif
/*----------------------------------------------------------------------*/
/*------------------------ get the size of the tentative prolong. block */
*nrow=0;
for (i=0; i<actagg->nblock; i++)
   (*nrow) += actagg->block[i][0];
   
*ncol = actagg->numdf;
if (*nrow >= 1000 || *ncol >= 500)
dserror("Local variable aggblock[1000][500] too small");
/*--------------------------------------- get the global column indizes */
for (i=0; i<actagg->numdf; i++)
   cindex[i] = actagg->dof[i];
/*------------------------------------------ get the global row indizes */   
counter=0;
for (i=0; i<actagg->nblock; i++)
   for (j=0; j<actagg->block[i][0]; j++)
   {
      rindex[counter] = actagg->block[i][j+1];
      counter++;
   }
dsassert(counter==*nrow,"Number of dofs in prolongator wrong");
qsort((INT*)rindex,*nrow,sizeof(INT),cmp_int);
/*------------------------------ get the nodal patch from the partition */
/*======================================================================*/
dsassert(actagg->nblock<=200,"Local variable prevagg[200] too small");
dsassert(*ncol==6,"number of rbm's has to be 6 in this case");
/* 
the index of the agg in prevlev corresponds to the index of the block in
blocks
*/
counter=0;
for (i=0; i<actagg->nblock; i++)
{
   actblock = actagg->block[i];
   dof      = actblock[1];
   for (j=0; j<prevlevel->nagg; j++)
   {
       actprevagg = &(prevlevel->agg[j]);
       for (k=0; k<actprevagg->numdf; k++)
       {
          if (dof == actprevagg->dof[k])
          {
             prevagg[counter] = actprevagg;
             counter++;
             goto nextblock;
          }
       }
   }
   nextblock:;
}
dsassert(counter==actagg->nblock,"Cannot find aggregates on lower level");
/*--------------------------- make the coordinates of the new aggregate */
/*
        transX   transY  transZ   rotX       rotY       rotZ
-------------------------------------------------------------
x     |    1       0       0       0          z-z0      -y+y0
y     |    0       1       0      -z+z0       0          x-x0
z     |    0       0       1       y-y0      -x+x0       0
rotx  |    0       0       0       1          0          0
roty  |    0       0       0       0          1          0
rotz  |    0       0       0       0          0          1
*/
x0=0.0;y0=0.0;z0=0.0;
for (i=0; i<actagg->nblock; i++)
{
   x0 += prevagg[i]->x[0];
   y0 += prevagg[i]->x[1];
   z0 += prevagg[i]->x[2];
}
x0 /= (actagg->nblock);
y0 /= (actagg->nblock);
z0 /= (actagg->nblock);
actagg->x[0] = x0;
actagg->x[1] = y0;
actagg->x[2] = z0;
/*----------------------- loop the nodes and put values to the aggblock */
if (*nrow > 18000) dserror("local bins too small");
init_quick_find(rindex,*nrow,&shift,bins);
for (i=0; i<actagg->nblock; i++)
{
   /* these are the aggregates on prevlev which are part of actagg */
   actprevagg = prevagg[i];
   x          = actprevagg->x[0];
   y          = actprevagg->x[1];
   z          = actprevagg->x[2];
   /* loop the dofs of the previous aggregate */
   for (j=0; j<actprevagg->numdf; j++)
   {
      dof = actprevagg->dof[j];
      /* find the dof in the row index */
      index = quick_find(dof,rindex,*nrow,shift,bins);
      if (index==-1) dserror("Cannot find aggregate-local dof");
      switch(j)
      {
      case 0: /* x mid-surface */
         aggblock[index][0] = 1.0;
         aggblock[index][1] = 0.0;
         aggblock[index][2] = 0.0;
         aggblock[index][3] = 0.0;
         aggblock[index][4] = z-z0;
         aggblock[index][5] = y0-y;
      break;
      case 1: /* y mid-surface */
         aggblock[index][0] = 0.0;
         aggblock[index][1] = 1.0;
         aggblock[index][2] = 0.0;
         aggblock[index][3] = z0-z;
         aggblock[index][4] = 0.0;
         aggblock[index][5] = x-x0;
      break;
      case 2: /* z mid-surface */
         aggblock[index][0] = 0.0;
         aggblock[index][1] = 0.0;
         aggblock[index][2] = 1.0;
         aggblock[index][3] = y-y0;
         aggblock[index][4] = x0-x;
         aggblock[index][5] = 0.0;
      break;
      case 3: /* rotx */
         aggblock[index][0] = 0.0;
         aggblock[index][1] = 0.0;
         aggblock[index][2] = 0.0;
         aggblock[index][3] = 1.0;
         aggblock[index][4] = 0.0;
         aggblock[index][5] = 0.0;
      break;
      case 4: /* roty */
         aggblock[index][0] = 0.0;
         aggblock[index][1] = 0.0;
         aggblock[index][2] = 0.0;
         aggblock[index][3] = 0.0;
         aggblock[index][4] = 1.0;
         aggblock[index][5] = 0.0;
      break;
      case 5: /* rotz */
         aggblock[index][0] = 0.0;
         aggblock[index][1] = 0.0;
         aggblock[index][2] = 0.0;
         aggblock[index][3] = 0.0;
         aggblock[index][4] = 0.0;
         aggblock[index][5] = 1.0;
      break;
      default:
         dserror("j out of range");
      break;
      }
   }/* end of for (j=0; j<actprevagg->numdf; j++) */
} /* end of for (i=0; i<actagg->nblock; i++) */
/*======================================================================*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_oneP_vanek */


/*!---------------------------------------------------------------------
\brief create the tentative prolongator for one aggregate                                          

<pre>                                                        m.gee 10/02 

</pre>
\param actagg         AGG*    (i/o) the active aggregate
\param aggblock       DOUBLE[1000][500] (o) the aggregate's block prolongator
\param rindex         INT[1000]         (o) global indizes of aggblock
\param cindex         INT[500]          (o) global indizes of aggblock
\param nrow           INT*              (o) dimension of rindex
\param ncol           INT*              (o) dimension of cindex
\param actstiff       DBCSR*            (i) fine grid stiffness matrix
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_oneP0_vanek(AGG     *actagg,
                               DOUBLE   aggblock[][500],
                               INT      rindex[],
                               INT      cindex[],
                               INT     *nrow,
                               INT     *ncol,
                               DBCSR   *actstiff)
{
INT           i,j,k,counter;
NODE         *node[200];
NODE         *actnode;
PARTDISCRET  *actpdis;
INT          *actblock;
INT           dof;
DOUBLE        x0,y0,z0,x,y,z,a1,a2,a3;
INT           index;
INT           shift,bins[5000];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_oneP0_vanek");
#endif
/*----------------------------------------------------------------------*/
/*------------------------ get the size of the tentative prolong. block */
*nrow=0;
for (i=0; i<actagg->nblock; i++)
   (*nrow) += actagg->block[i][0];
   
*ncol = actagg->numdf;
if (*nrow >= 1000 || *ncol >= 500)
dserror("Local variable aggblock[1000][500] too small");
/*--------------------------------------- get the global column indizes */
for (i=0; i<actagg->numdf; i++)
   cindex[i] = actagg->dof[i];
/*------------------------------------------ get the global row indizes */   
counter=0;
for (i=0; i<actagg->nblock; i++)
   for (j=0; j<actagg->block[i][0]; j++)
   {
      rindex[counter] = actagg->block[i][j+1];
      counter++;
   }
dsassert(counter==*nrow,"Number of dofs in prolongator wrong");
qsort((INT*)rindex,*nrow,sizeof(INT),cmp_int);
/*------------------------------ get the nodal patch from the partition */
/*======================================================================*/
/*======================================================================*/
#if 1 /*-- this is the calculation of the rigid body mode for SHELL8 !!!*/
/*======================================================================*/
/*======================================================================*/
dsassert(actagg->nblock<200,"Local variable node[200] too small");
dsassert(*ncol==6,"number of rbm's has to be 6 in this case");
/* 
the index of the node in node corresponds to the index of the block in
blocks
*/
counter=0;
actpdis = mlprecond.partdis;
for (i=0; i<actagg->nblock; i++)
{
   actblock = actagg->block[i];
   dof      = actblock[1];
   for (j=0; j<actpdis->numnp; j++)
   {
      actnode = actpdis->node[j];
      for (k=0; k<actnode->numdf; k++)
      {
         if (actnode->dof[k]>=actstiff->numeq_total) continue;
         if (dof == actnode->dof[k])
         {
            node[counter] = actnode;
            counter++;
            goto nextblock;
         }
      }
   }
   nextblock:;
}
dsassert(counter==actagg->nblock,"Cannot find nodes to matrix blocks");
/*----------------------- make the nodal coordinate center of the patch */
/* 
IMPORTANT NOTE:
in the nonlinear dynaic or static case, it could be advisable to use the
deformed configuration to create the rigid body modes !

the discrete representation of the rigid body of shell8 is:

      transX   transY  transZ   rotX       rotY       rotZ
-----------------------------------------------------------
x   |    1       0       0       0          z-z0      -y+y0
y   |    0       1       0      -z+z0       0          x-x0
z   |    0       0       1       y-y0      -x+x0       0
dx  |    0       0       0       0          a3        -a2
dy  |    0       0       0      -a3         0          a1
dz  |    0       0       0       a2        -a1         0

*/
x0=0.0;y0=0.0;z0=0.0;
for (i=0; i<actagg->nblock; i++)
{
   x0 += node[i]->x[0];
   y0 += node[i]->x[1];
   z0 += node[i]->x[2];
}
x0 /= (actagg->nblock);
y0 /= (actagg->nblock);
z0 /= (actagg->nblock);
actagg->x[0] = x0;
actagg->x[1] = y0;
actagg->x[2] = z0;
/*----------------------- loop the nodes and put values to the aggblock */
if (*nrow > 18000) dserror("local bins too small");
init_quick_find(rindex,*nrow,&shift,bins);
for (i=0; i<actagg->nblock; i++)
{
   actnode = node[i];
   for (j=0; j<mlprecond.director.fdim; j++)
   if (actnode == mlprecond.node[j]) break;
   dsassert(j<mlprecond.director.fdim,"Could not find node");
   x  = actnode->x[0];
   y  = actnode->x[1];
   z  = actnode->x[2];
   a1 = mlprecond.director.a.da[j][0];
   a2 = mlprecond.director.a.da[j][1];
   a3 = mlprecond.director.a.da[j][2];
   /*-------------------------------------------- make rigid body modes */
   /* loop the dofs of the nodes */
   for (j=0; j<actnode->numdf; j++)
   {
      dof = actnode->dof[j];
      if (dof >= actstiff->numeq_total) continue;
      /* find the dof in the row index */
      index = quick_find(dof,rindex,*nrow,shift,bins);
      if (index==-1) dserror("Cannot find aggregate-local dof");
      switch(j)
      {
      case 0:/* x mid-surface */
         aggblock[index][0] = 1.0;
         aggblock[index][1] = 0.0;
         aggblock[index][2] = 0.0;
         aggblock[index][3] = 0.0;
         aggblock[index][4] = z-z0;
         aggblock[index][5] = -y+y0;
      break;
      case 1:/* y mid-surface */
         aggblock[index][0] = 0.0;
         aggblock[index][1] = 1.0;
         aggblock[index][2] = 0.0;
         aggblock[index][3] = -z+z0;
         aggblock[index][4] = 0.0;
         aggblock[index][5] = x-x0;
      break;
      case 2:/* z mid-surface */
         aggblock[index][0] = 0.0;
         aggblock[index][1] = 0.0;
         aggblock[index][2] = 1.0;
         aggblock[index][3] = y-y0;
         aggblock[index][4] = -x+x0;
         aggblock[index][5] = 0.0;
      break;
      case 3:/* dx of difference vector */
         aggblock[index][0] = 0.0;
         aggblock[index][1] = 0.0;
         aggblock[index][2] = 0.0;
         aggblock[index][3] = 0.0;
         aggblock[index][4] = a3;
         aggblock[index][5] = -a2;
      break;
      case 4:/* dy of difference vector */
         aggblock[index][0] = 0.0;
         aggblock[index][1] = 0.0;
         aggblock[index][2] = 0.0;
         aggblock[index][3] = -a3;
         aggblock[index][4] = 0.0;
         aggblock[index][5] = a1;
      break;
      case 5:/* dz of difference vector */
         aggblock[index][0] = 0.0;
         aggblock[index][1] = 0.0;
         aggblock[index][2] = 0.0;
         aggblock[index][3] = a2;
         aggblock[index][4] = -a1;
         aggblock[index][5] = 0.0;
      break;
      default:
         dserror("j out of range");
      break;
      }/* end of switch(j) that is type of dof of the node */
   }/* end of for (j=0; j<actnode->numdf; j++) */
} /* end of for (i=0; i<actagg->nblock; i++) */
/*======================================================================*/
/*======================================================================*/
#endif /*- this is the calculation of the rigid body mode for SHELL8 !!!*/
/*======================================================================*/
/*======================================================================*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_oneP0_vanek */



/*!---------------------------------------------------------------------
\brief smoothes the columns of the prolongator P by damped Jacobi                                             

<pre>                                                        m.gee 10/02 

</pre>
\param P          DBCSR*       (i/o) the Prolongator
\param block      DOUBLE[][500](i)   working matrix
\param rindex     INT*         (i)   row indize of block
\param cindex     INT*         (i)   column indize of block
\param nrow       INT*         (i)   row dimension of block
\param ncol       INT*         (i)   column dimension of block
\param actstiff   DBCSR*       (i)   the stiffness matrix to be smoothe 
\param agg        AGG*         (i)   the aggregates the proongator belongs to
\param nagg       INT          (i)   the number of aggregates on this proc 
\param actintra   INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_smoothP(DBCSR *P, DBCSR *actstiff, INTRA *actintra)
{
INT           i,j,k,n,m,counter;
INT           myrank,nproc,tag;
DOUBLE        omega,fac;
INT           numeq,numeq_total;
INT           firstcol,lastcol,maxcol,mincol,intercol;
INT          *ia,*ja,*update;
DOUBLE       *a;
INT           colstart,colend,actrow,actcol,actcolP,diag;
INT           owner,index;
ARRAY         a_a;

INT           rcol[1000];
DOUBLE        col[1000];
INT           nr;

DOUBLE        sum;
INT           foundit;

INT           myinters[MAXPROC][2],myinterr[MAXPROC][2];
INT           needits[MAXPROC][MAXPROC],needitr[MAXPROC][MAXPROC];
INT           nsend,nrecv;
INT           nc;
ARRAY         isend_a,dsend_a;
INT         **isend;
DOUBLE      **dsend;
INT           sendsize[2];
INT           recvsize[2];
ARRAY         irecv_a,drecv_a;
INT         **irecv;
DOUBLE      **drecv;
INT          *rc;
DOUBLE       *c;
INT           shift;
INT           bins[5000];
INT           max,min;
#ifdef PARALLEL 
MPI_Status    status;  
MPI_Request  *request;   
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_smoothP");
#endif
/*----------------------------------------------------------------------*/
mlpcg_csr_close(P); /* this should be unneccessary, but for unknown reason, it's not */
/*----------------------------------------------------------------------*/
myrank      = actintra->intra_rank;
nproc       = actintra->intra_nprocs;
omega       = mlprecond.omega;
numeq       = actstiff->numeq;
numeq_total = actstiff->numeq_total;
update      = actstiff->update.a.iv;
ia          = actstiff->ia.a.iv;
ja          = actstiff->ja.a.iv;
/*----------------------------------------------------------------------*/
/*
for (i=0; i<numeq; i++)
{
   actrow = P->update.a.iv[i];
   for (j = P->ia.a.iv[i]; j<P->ia.a.iv[i+1]; j++)
   {
      if (P->ja.a.iv[j]==-1)
         continue;
      actcol = P->ja.a.iv[j];
      printf("actrow %d actcol %d val %E\n",actrow,actcol,P->a.a.dv[j]);
   }
}
*/
/*----------------------------------------------------------------------*/
mlpcg_matvec_init(actstiff,actintra);
/*--------------------------------- build the smoother stiffness matrix */
/* copy the values array */
am_alloc_copy(&(actstiff->a),&(a_a));
a = a_a.a.dv;
/* scale the rows by -omega/diag[actrow] */
for (i=0; i<numeq; i++)
{
   actrow = update[i];
   colstart = ia[i];
   colend   = ia[i+1];
   /* find the diagonal element */
   for (j=colstart; j<colend; j++)
   {
      if (actrow==ja[j])
      {
         if (FABS(a[j])<EPS14) dserror("Zero diagonal element detected");
         fac  = -omega/a[j];
         diag = j;
         break;
      }
   }
   /* scale the row */
   for (j=colstart; j<colend; j++)
      a[j] *= fac;
   /* add a one to the diagonal element */
   a[diag] += 1.0;
}
/*------------------------------------------ find the column range of P */
/*
   firtscol = first column on this processor
   lastcol  = last column on this processor
   mincol   = globally lowest column (=0)
   maxcol   = globally highest column (=lastcol on proc nproc-1)
*/

firstcol = VERYLARGEINT;
lastcol  = 0;
for (i=0; i<P->ia.a.iv[P->numeq]; i++)
{
   actcol = P->ja.a.iv[i];
   if (actcol==-1) continue;
   if (actcol<firstcol) firstcol = actcol;
   if (actcol>lastcol)  lastcol  = actcol;
}
i = lastcol;/* this is a dummy operation with lastcol necessary, because */
lastcol = i;/* otherwise the compiler strangely eliminates the value of lastcol (??) */
mincol = 0;
#ifdef PARALLEL
MPI_Allreduce(&lastcol,&maxcol,1,MPI_INT,MPI_MAX,actintra->MPI_INTRA_COMM);
#else
maxcol = lastcol;
#endif
/*-------------------------------- find the first interproc column of P */
/*
   intercol = the first local column, where rows occur, that have interproc 
              off-diagonal entries 
*/
intercol = VERYLARGEINT;
index    = mlpcg_getindex(P->firstcoupledof,P->update.a.iv,numeq);
if (index==-1) dserror("Cannot find local dof");
for (i=index; i<numeq; i++)
{
   actrow   = P->update.a.iv[i];
   colstart = P->ia.a.iv[i];
   colend   = P->ia.a.iv[i+1];
   for (j=colstart; j<colend; j++)
   {
      if (P->ja.a.iv[j]==-1) continue;
      if (P->ja.a.iv[j] < intercol)
         intercol = P->ja.a.iv[j];
      break;
   }
}
dsassert(intercol<VERYLARGEINT,"No intercol found");
/*======================================make interproc smoothing part I */
#ifdef PARALLEL 
if (nproc > 1) {
/* make myinterrows array */
for (n=0; n<nproc; n++)
for (m=0; m<2; m++) 
   myinters[n][m] = 0;
myinters[myrank][0] = P->firstcoupledof;
myinters[myrank][1] = P->update.a.iv[numeq-1];
MPI_Allreduce(&(myinters[0][0]),&(myinterr[0][0]),MAXPROC*2,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
/* loop my piece of matrix, and check, from who I need the interproc columns */
/* 
   I need the interproc columns from      processors n : needitr[myrank][n]==1
   I have to send my interproc columns to processors n : needitr[n][myrank]==1
*/
for (n=0; n<nproc; n++)
for (m=0; m<nproc; m++) 
   needits[n][m] = 0;
for (j=0; j<ia[numeq]; j++)
{
   actcol = ja[j];
   if (actcol==-1) continue;
   owner  = mlpcg_getowner(actcol,actstiff->owner,nproc);
   if (owner==myrank) 
      continue;
   for (n=0; n<nproc; n++)
   {
      if (n==myrank) 
         continue;
      if (actcol >= myinterr[n][0] && actcol <= myinterr[n][1])
         needits[myrank][n] = 1;
   }
}
MPI_Allreduce(&(needits[0][0]),&(needitr[0][0]),MAXPROC*MAXPROC,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
/*------------------- count number of sends and receives I have to make */
nsend=0;
nrecv=0;
for (n=0; n<nproc; n++)
{
   if (n==myrank) continue;
   if (needitr[myrank][n]==1) nrecv++;
   if (needitr[n][myrank]==1) nsend++;
}
/*----------------------- allocate sendbuffers for my interproc columns */
/* number of interproc columns */
nc    = lastcol - intercol + 1;
isend = amdef("tmp",&isend_a,nc,50,"IA");
dsend = amdef("tmp",&dsend_a,nc,50,"DA");
/* fill the sendbuffers */
counter=0;
for (i=intercol; i<=lastcol; i++)
{
   actcolP = i;
   mlpcg_extractcollocal(P,actcolP,col,rcol,&nr);
   if (nr > isend_a.sdim-1)
   {
      isend = amredef(&isend_a,nc,nr+1,"IA");
      dsend = amredef(&dsend_a,nc,nr+1,"DA");
   }
   dsassert(counter<isend_a.fdim,"buffer overflow");
   isend[counter][0] = nr;
   dsend[counter][0] = (DOUBLE)actcolP;
   for (j=0; j<nr; j++)
   {
      dsassert(j+1<isend_a.sdim,"buffer overflow");
      dsassert(j+1<dsend_a.sdim,"buffer overflow");
      isend[counter][j+1] = rcol[j];
      dsend[counter][j+1] = col[j];
   }
   counter++;
}
dsassert(nc==counter,"number of sends wrong");
/* allocate requests */
request = (MPI_Request*)CCAMALLOC(3*nsend*sizeof(MPI_Request));
/* loop all procs and make sends */
sendsize[0] = isend_a.fdim;
sendsize[1] = isend_a.sdim;
counter=0;
for (n=0; n<nproc; n++)
{
   if (n==myrank || needitr[n][myrank]!=1)
      continue;
   MPI_Isend(sendsize,2,MPI_INT,n,counter,actintra->MPI_INTRA_COMM,&(request[counter]));
   counter++;
   MPI_Isend(isend[0],sendsize[0]*sendsize[1],MPI_INT,n,counter,actintra->MPI_INTRA_COMM,&(request[counter]));
   counter++;
   MPI_Isend(dsend[0],sendsize[0]*sendsize[1],MPI_DOUBLE,n,counter,actintra->MPI_INTRA_COMM,&(request[counter]));
   counter++;
}
if (counter != nsend*3) dserror("Number of sends wrong");
/* end of if (nproc > 1) */
}
#endif
/*======================================================================*/






/*=================================================make local smoothing */
/*------------------------ loop all columns in P and do local smoothing */
for (i=firstcol; i<=lastcol; i++)
{
   /* extract the column from the prolongator P */
   actcolP = i;
   mlpcg_extractcollocal(P,actcolP,col,rcol,&nr);
   if (nr==0)
      continue;
   if (nr > 18000) dserror("static bins too small");
   init_quick_find(rcol,nr,&shift,bins); 
   min = rcol[0];
   max = rcol[nr-1]; 
   /* make the multiplication with the modfied stiffness */
   /* loop the row of the mod. stiffness matrix */
   for (j=0; j<numeq; j++)
   {
      actrow   = update[j];
      colstart = ia[j];
      colend   = ia[j+1];
      sum      = 0.0;
      foundit  = 0;
      if (ja[colstart]>max || ja[colend-1]<min)
         continue;
      /* loop the columns of the stiffness matrix */
      for (k=colstart; k<colend; k++)
      {
         if (ja[k] < min || ja[k] > max)
            continue;
         /* find the column in the prolongator column */
         index = quick_find(ja[k],rcol,nr,shift,bins);
         if (index==-1)
            continue;
         foundit=1;
         /* make the multiplication */
         sum += a[k] * col[index];
      } 
      /* put the value back to the prolongator */
      if (foundit && FABS(sum)>EPS14)
         mlpcg_csr_setentry(P,sum,actrow,actcolP,actintra);
   }      
}
/*======================================================================*/

/*
for (i=0; i<numeq; i++)
{
   actrow = P->update.a.iv[i];
   for (j = P->ia.a.iv[i]; j<P->ia.a.iv[i+1]; j++)
   {
      if (P->ja.a.iv[j]==-1)
         continue;
      actcol = P->ja.a.iv[j];
      printf("actrow %d actcol %d val %E\n",actrow,actcol,P->a.a.dv[j]);
   }
}
*/


/*=====================================make interproc smoothing part II */
#ifdef PARALLEL 
if (nproc > 1) {
/* make receive */
while (nrecv != 0)
{
   MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG,actintra->MPI_INTRA_COMM,&status);
   /* get the sender , size and tag */
   n   = status.MPI_SOURCE;
   tag = status.MPI_TAG;
   MPI_Get_count(&status,MPI_INT,&i);
   dsassert(i==2,"first message not of length 2");
   /* receive the size */
   MPI_Recv(recvsize,2,MPI_INT,n,tag,actintra->MPI_INTRA_COMM,&status);
   /* allocate the right buffer */
   irecv = amdef("tmp",&irecv_a,recvsize[0],recvsize[1],"IA");
   drecv = amdef("tmp",&drecv_a,recvsize[0],recvsize[1],"DA");
   /* receive the integer message */
   MPI_Recv(irecv[0],recvsize[0]*recvsize[1],MPI_INT,n,tag+1,actintra->MPI_INTRA_COMM,&status);
   MPI_Get_count(&status,MPI_INT,&i);
   dsassert(i==recvsize[0]*recvsize[1],"buffer overflow");
   /* receive the DOUBLE message */
   MPI_Recv(drecv[0],recvsize[0]*recvsize[1],MPI_DOUBLE,n,tag+2,actintra->MPI_INTRA_COMM,&status);
   MPI_Get_count(&status,MPI_DOUBLE,&i);
   dsassert(i==recvsize[0]*recvsize[1],"buffer overflow");
   /* decrease number of unreceived sets */
   nrecv--;
   /* get number of columns received */
   nc = recvsize[0];
   /* loop the received columns */
   for (i=0; i<nc; i++)
   {
      nr      =   irecv[i][0];
      rc      = &(irecv[i][1]);
      actcolP = (INT)drecv[i][0];
      c       = &(drecv[i][1]);
      if (nr==0)
         continue;
      if (nr > 18000) dserror("local bins too small");
      init_quick_find(rc,nr,&shift,bins);
      min = rc[0];
      max = rc[nr-1];
      /* make the multiplication with the modified stiffness */
      /* loop the processors rows */
      for (j=0; j<numeq; j++)
      {
         actrow   = update[j];
         colstart = ia[j];
         colend   = ia[j+1];
         if (ja[colstart]>max || ja[colend-1]<min)
            continue;
         sum      = 0.0;
         foundit  = 0;
         /* loop the columns in the row */
         for (k=colstart; k<colend; k++)
         {
            if (ja[k] < min || ja[k] > max) continue;
            /* find actcol in the received prolongator columnd */
            index = quick_find(ja[k],rc,nr,shift,bins);
            if (index==-1) 
               continue;
            foundit = 1;
            /* make the multiplication */
            sum += a[k] * c[index];
         }
         /* if found, then put value to P */
         if (foundit && FABS(sum)>EPS14)
            mlpcg_csr_setentry(P,sum,actrow,actcolP,actintra);
      }
   }
   amdel(&irecv_a);
   amdel(&drecv_a);
}
for (i=0; i<nsend*3; i++)
   MPI_Wait(&(request[i]),&status);
CCAFREE(request);   
/* end of if (nproc > 1) */
}
#endif
/*======================================================================*/





/*============================================================= tidy up */
amdel(&a_a);
#ifdef PARALLEL 
if (nproc > 1)
{
   amdel(&isend_a);
   amdel(&dsend_a);
}
#endif
/*======================================================================*/
/*
for (i=0; i<numeq; i++)
{
   actrow = P->update.a.iv[i];
   for (j = P->ia.a.iv[i]; j<P->ia.a.iv[i+1]; j++)
   {
      if (P->ja.a.iv[j]==-1)
         continue;
      actcol = P->ja.a.iv[j];
      printf("actrow %d actcol %d val %E\n",actrow,actcol,P->a.a.dv[j]);
   }
}
*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_smoothP */






#endif
/*! @} (documentation module close)*/

/*!---------------------------------------------------------------------
\file
\brief contains the multilevel preconditioner for shells

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/prototypes_mlpcg.h"
int cmp_int(const void *a, const void *b );
double cmp_double(const void *a, const void *b );
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
from the aggregation done before 
</pre>
\param actlev      MLLEVEL*     (i/o) the active level in the ml-precond.                   
\param actintra   INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_P0(MLLEVEL  *actlev, INTRA *actintra)
{
int          i,j,k,n,counter=0;
int          myrank,nproc;
DBCSR       *actstiff;
DBCSR       *P;
double     **R;
AGG         *actagg;
int          min;
int          nrow,ncol;
double       aggblock[1000][500];
int          rindex[1000],cindex[500];
int          ndofs[MAXPROC],ndofr[MAXPROC];
int          firstdof;
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
{
   actlev->P = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
   if (!actlev->P) dserror("Allocation of memory failed");
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
make the guess 20% too large 
*/
   counter = (int)(counter*1.2);
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
/*----------------------------------------------------------------------*/
/*
get the directors of the shell elements 
*/
mlpcg_precond_getdirs();
/*----------------------------------------------------------------------*/
/* 
loop the aggregates and create the tentative prolongator  
*/                  
for (i=0; i<actlev->nagg; i++)
{
   actagg = &(actlev->agg[i]);
   /* create the tentative prolongator diagonal block of this aggregate */
   mlpcg_precond_oneP0(actagg,aggblock,rindex,cindex,&nrow,&ncol,actstiff);
   if (!(actagg->tentP))
   {
      actagg->tentP = (ARRAY*)CCACALLOC(1,sizeof(ARRAY));
      if (!(actagg->tentP)) dserror("Allocation of memory failed");
      amdef("tentP",actagg->tentP,nrow,ncol,"DA");
   }
   if (!(actagg->tentP_rindex))
      actagg->tentP_rindex = (int*)CCAMALLOC(nrow*sizeof(int));
   if (!(actagg->tentP_rindex)) dserror("Allocation of memory failed");
   /*------------------------------- put rindex to actagg->tentP_rindex */
   for (j=0; j<nrow; j++) actagg->tentP_rindex[j] = rindex[j];
   /*--------------------------- put the aggblock to actagg->tentP.a.da */
   for (j=0; j<ncol; j++)
   for (k=0; k<nrow; k++)
   actagg->tentP->a.da[k][j] = aggblock[k][j];
   /*------------------ set number of rows in this piece of Prolongator */
   actagg->tentP_nrow = nrow;
   /*--------- make gram-schmidt orthogonalization of this block P = QR */
   /* do this later.....
   if (actagg->R==NULL)
   {
      actagg->R = (ARRAY*)MALLOC(sizeof(ARRAY));
      if (!(actagg->R)) dserror("Allocation of memory failed");
      R = amdef("R",actagg->R,actagg->numdf,actagg->numdf,"DA");
   }
   else
   {
      dsassert(actagg->R->Typ==cca_DA,"R not allocated in aggregate");
      R = actagg->R->a.da;
   }
   mlpcg_precond_gramschmidt(aggblock,R,nrow,ncol);*/
   /*------------------------ store the R part to use in the next level */
} /* end of for (i=0; i<actlev->nagg; i++) */
/*--------------- loop all aggregates again and fill the DBCSR matrix P */   
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
/*------------------------------- make the smoothing of the prolongator */
mlpcg_smoothP(P,aggblock,rindex,cindex,&nrow,&ncol,actlev->csr,
              actlev->agg,actlev->nagg,actintra);
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

<pre>                                                        m.gee 10/02 

</pre>
\param actagg         AGG*    (i/o) the active aggregate
\param aggblock       double[1000][500] (o) the aggregate's block prolongator
\param rindex         int[1000]         (o) global indizes of aggblock
\param cindex         int[500]          (o) global indizes of aggblock
\param nrow           int*              (o) dimension of rindex
\param ncol           int*              (o) dimension of cindex
\param actstiff       DBCSR*            (i) fine grid stiffness matrix
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_oneP0(AGG     *actagg,
                         double   aggblock[][500],
                         int      rindex[],
                         int      cindex[],
                         int     *nrow,
                         int     *ncol,
                         DBCSR   *actstiff)
{
int           i,j,k,l,counter;
NODE         *node[200];
NODE         *actnode;
PARTDISCRET  *actpdis;
int          *actblock;
int           dof;
double        x0,y0,z0,x,y,z,a1,a2,a3;
int           index;
int           foundit;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_oneP0");
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
qsort((int*)rindex,*nrow,sizeof(int),cmp_int);
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
/*----------------------- loop the nodes and put values to the aggblock */
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
      index = find_index(dof,rindex,*nrow);
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
/*--------------------------------- check the aggblock for zero columns */
#if 0 /* delete zero columns */
for (i=0; i<*ncol; i++)
{
   foundit=0;
   for (j=0; j<*nrow; j++)
      if (FABS(aggblock[j][i])>EPS14) 
      {
         foundit=1;
         break;
      }
   if (!foundit)/* detected an empty prolongator column ! */
   {
      printf("WARNING: zero prolongator column in aggregate in column %d\n",actagg->dof[i]);
      /* delete this dof from the aggregate , it's in place i*/
      if (i < actagg->numdf-1)
         for (j=i; j<actagg->numdf-1; j++)
         {
            actagg->dof[j] = actagg->dof[j+1];
            cindex[j] = cindex[j+1];
            for (k=0; k<*nrow; k++)
               aggblock[k][j] = aggblock[k][j+1];
         }
      actagg->numdf--;
      *ncol = *ncol - 1;
   }
}
#endif
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
} /* end of mlpcg_precond_oneP0 */



/*!---------------------------------------------------------------------
\brief smoothes the columns of the prolongator P by damped Jacobi                                             

<pre>                                                        m.gee 10/02 

</pre>
\param P          DBCSR*       (i/o) the Prolongator
\param block      double[][500](i)   working matrix
\param rindex     int*         (i)   row indize of block
\param cindex     int*         (i)   column indize of block
\param nrow       int*         (i)   row dimension of block
\param ncol       int*         (i)   column dimension of block
\param actstiff   DBCSR*       (i)   the stiffness matrix to be smoothe 
\param agg        AGG*         (i)   the aggregates the proongator belongs to
\param nagg       int          (i)   the number of aggregates on this proc 
\param actintra   INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_smoothP(DBCSR *P, double block[][500], int *rindex, int *cindex,
                  int *nrow, int *ncol, DBCSR *actstiff, 
                  AGG *agg, int nagg, INTRA *actintra)
{
int           i,j,k,l,n,m,counter;
int           myrank,nproc,flag,tag;
double        omega,fac;
int           numeq,numeq_total;
int           firstcol,lastcol,maxcol,mincol,intercol;
int          *ia,*ja,*update;
double       *a;
int           colstart,colend,actrow,actcol,actcolP,diag;
int           owner,index;
ARRAY         a_a;

int           rcol[1000];
double        col[1000];
int           nr;

double        sum;
int           foundit;

int           myinters[MAXPROC][2],myinterr[MAXPROC][2];
int           needits[MAXPROC][MAXPROC],needitr[MAXPROC][MAXPROC];
int           nsend,nrecv;
int           nc;
ARRAY         isend_a,dsend_a;
int         **isend;
double      **dsend;
int           sendsize[2];
int           recvsize[2];
ARRAY         irecv_a,drecv_a;
int         **irecv;
double      **drecv;
int          *rc;
double       *c;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_smoothP");
#endif
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
mincol = 0;
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
/*======================================make interproc smoothing part I */
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
   /* make the multiplication with the modfied stiffness */
   /* loop the row of the mod. stiffness matrix */
   for (j=0; j<numeq; j++)
   {
      actrow   = update[j];
      colstart = ia[j];
      colend   = ia[j+1];
      sum      = 0.0;
      foundit  = 0;
      /* loop the columns of the stiffness matrix */
      for (k=colstart; k<colend; k++)
      {
         actcol = ja[k];
         /* find the column in the prolongator column */
         index = find_index(actcol,rcol,nr);
         if (index==-1)
            continue;
         foundit=1;
         /* make the multiplication */
         sum += a[k] * col[index];
      } 
      /* put the value back to the prolongator */
      if (foundit)
         mlpcg_csr_setentry(P,sum,actrow,actcolP,actintra);
   }      
}
/*======================================================================*/





/*=====================================make interproc smoothing part II */
/*======================================================================*/





/*============================================================= tidy up */
amdel(&a_a);
/*======================================================================*/



/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_smoothP */






/*! @} (documentation module close)*/

/*!---------------------------------------------------------------------
\file
\brief contains the multilevel preconditioner for shells

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/prototypes_mlpcg.h"
#include "../shell8/shell8.h"
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
\brief smoothed aggregation algebraic multigrid preconditioner                                              

<pre>                                                        m.gee 10/02 
on input r is the residuum, it is not altered on output !
on input z is zero!         it is     altered on output
</pre>
\param actlevel     int          (i)   the active level in the grid hierachy
\param z_a          ARRAY*       (o)   the correction on this active level
\param r_a          ARRAY*       (i)   the residuum on this active level
\param actintra     INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_amg(int    level,
                       ARRAY *z_a,
                       ARRAY *r_a,
                       INTRA *actintra)
{
int      i;
DBCSR   *stiff;
MLLEVEL *actlevel,*nextlevel;
int      nlevel;
int      numeq;

double  *z,*r;
ARRAY    rwork_a;
double  *rwork;
ARRAY    zwork_a;
double  *zwork;

double  *zc,*rc;
ARRAY    rcwork_a;
double  *rcwork;
ARRAY    zcwork_a;
double  *zcwork;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_amg");
#endif
/*----------------------------------------------------------------------*/
nlevel   = mlprecond.numlev;
actlevel = &(mlprecond.level[level]);
stiff    = actlevel->csr;
numeq    = stiff->numeq;
z        = z_a->a.dv;
r        = r_a->a.dv;
/*-------------------------------------allocate working copy of r and z */
rwork = amdef("rwork",&rwork_a,numeq,1,"DV");
zwork = amdef("zwork",&zwork_a,numeq,1,"DV");
/*----------------------------------------------------- copy r to rwork */
mlpcg_updvec(rwork,r,1.0,1,numeq);

/*=================================== on coarsest level do coarse solve */
if (level == nlevel-1)
{
   mlpcg_precond_coarsesolv(z,r,actlevel,actintra);
   goto exit;
}
/*======================================================================*/


/*---------------------------------------------------tidy up this level */
exit:
amdel(&rwork_a);
amdel(&zwork_a);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_amg */
/*!---------------------------------------------------------------------
\brief make presmoothing                                              

<pre>                                                        m.gee 10/02 

</pre>
\param z            double*      (o)   the smoothed residuum, on input zero!
\param r            double*      (i)   the residuum
\param csr          DBCSR*       (i)   the matrix on the active level 
\param lev          MLLEVEL*     (i)   the structure of the active level
\param actintra     INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_presmo(double *z, double *r, DBCSR *csr, MLLEVEL *lev, INTRA *actintra)
{
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_presmo");
#endif
/*----------------------------------------------------------------------*/
switch(lev->presmoother)
{
case pre_none:
break;
case pre_fwdGS:
break;
case pre_Jacobi:
   mlpcg_precond_smoJacobi(z,r,csr,lev->presweep,actintra);
break;
case pre_ilu:
break;
default:
   dserror("Unknown type of presmoother");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_presmo */
/*!---------------------------------------------------------------------
\brief make poistsmoothing                                              

<pre>                                                        m.gee 10/02 

</pre>
\param z            double*      (o)   the smoothed residuum, on input zero!
\param r            double*      (i)   the residuum
\param csr          DBCSR*       (i)   the matrix on the active level 
\param lev          MLLEVEL*     (i)   the structure of the active level
\param actintra     INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_postsmo(double *z, double *r, DBCSR *csr, MLLEVEL *lev, INTRA *actintra)
{
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_postsmo");
#endif
/*----------------------------------------------------------------------*/
switch(lev->postsmoother)
{
case post_none:
break;
case post_bckGS:
break;
case post_Jacobi:
   mlpcg_precond_smoJacobi(z,r,csr,lev->postsweep,actintra);
break;
case post_ilu:
break;
default:
   dserror("Unknown type of postsmoother");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_postsmo */
/*!---------------------------------------------------------------------
\brief make coarsest solution                                              

<pre>                                                        m.gee 11/02 

</pre>
\param z            double*      (o)   the smoothed residuum, on input zero!
\param r            double*      (i)   the residuum
\param lev          MLLEVEL*     (i)   the structure of the active level
\param actintra     INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_coarsesolv(double *z, double *r, MLLEVEL *lev, INTRA *actintra)
{
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_coarsesolv");
#endif
/*----------------------------------------------------------------------*/
switch(lev->coarsesolv)
{
case co_none:
break;
case co_ilu:
break;
case co_lapack:
   mlpcg_precond_lapacksolve(z,r,lev->csr,actintra);
break;
case co_spooles:
break;
default:
   dserror("Unknown type of coarse solver");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_coarsesolv */

 

/*!---------------------------------------------------------------------
\brief create multilevel preconditioner                                              

<pre>                                                        m.gee 9/02 

</pre>
\param mlpcgvars    MLPCGVARS*   (i)   variables needed for mlpcg                   
\param bdcsr        DBCSR*       (i)   the distributed  csr matrix                   
\param actintra     INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_create(DBCSR     *bdcsr, 
                          MLPCGVARS *mlpcgvars,
                          INTRA     *actintra)
{
int        i,j;
MLLEVEL   *actlev;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_create");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------------------------- set field and partition */
mlprecond.fielddis = mlpcgvars->fielddis;
mlprecond.partdis  = mlpcgvars->partdis;
mlprecond.omega    = mlpcgvars->p_omega;
/*----------------------------------------------------- allocate levels */
mlprecond.numlev = mlpcgvars->numlev;
mlprecond.level = (MLLEVEL*)CCACALLOC(mlprecond.numlev,sizeof(MLLEVEL));
if (!mlprecond.level) dserror("Allocation of MLPRECOND failed");
/*-------------- loop levels and set smoothers and prolongation damping */
for (i=0; i<mlprecond.numlev-1; i++)
{
   actlev               = &(mlprecond.level[i]);
   actlev->coarsesolv   = co_none;
   actlev->co_ilu_n     = 0;
   actlev->presmoother  = mlpcgvars->presmoother;
   actlev->presweep     = mlpcgvars->presweep;
   actlev->postsmoother = mlpcgvars->postsmoother;
   actlev->postsweep    = mlpcgvars->postsweep;
}
   actlev               = &(mlprecond.level[mlprecond.numlev-1]);
   actlev->coarsesolv   = mlpcgvars->coarsesolv;
   actlev->co_ilu_n     = mlpcgvars->co_ilu_n;
   actlev->presmoother  = pre_none;
   actlev->presweep     = 0;
   actlev->postsmoother = post_none;
   actlev->postsweep    = 0;
/*-- loop levels again and allocate a csr matrix on all levels except 0 */
for (i=1; i<mlprecond.numlev; i++)
{
   actlev               = &(mlprecond.level[i]);
   actlev->csr = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
   if (!actlev->csr) dserror("Allocation of DBCSR failed");
}
/*--------------------- set pointer in lowest level to the bdcsr matrix */
mlprecond.level[0].csr = bdcsr;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_create */



/*!---------------------------------------------------------------------
\brief init multilevel preconditioner                                              

<pre>                                                        m.gee 9/02 

</pre>
\param bdcsr      DBCSR*       (i)   the distributed  csr matrix                   
\param actintra   INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_init(DBCSR  *bdcsr,MLPCGVARS *mlpcgvars, INTRA *actintra)
{
int        i,j;
MLLEVEL   *actlev;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_init");
#endif
/*-------------------------------- prepare the levels and the smoothers */
actlev = &(mlprecond.level[0]);
/*============================================== make level0 separately */
/*------------------------------ init the smoothers on the active level */

/* make aggregates on the lowest level, using the grid and the partition */
mlpcg_precond_agg(actlev,actintra);
/*------- set dof numbers in the domain decomposition of the aggregates */
mlpcg_precond_aggsetdofs(actlev,mlpcgvars->numdf,actintra);
/*---------------------------------- create prolongator for finest grid */
mlpcg_precond_P0(actlev,actintra);
/*------------------------------- restrict the matrix to the next level */
mlpcg_precond_restrictK(actlev,actlev+1,actintra);
/*------------------------------------------------------- make printout */
if (actintra->intra_rank==0)
printf("level 0: size %d nnz %d\n",mlprecond.level[0].csr->numeq_total,
                                   mlprecond.level[0].csr->ja.fdim);
/*-------------------------------------------- loop level 1 to numlev-2 */
for (i=1; i<mlprecond.numlev-1; i++)
{
   actlev = &(mlprecond.level[i]);
   /*--------------------------- init the smoothers on the active level */
   /*---------------------------------- make aggregates on active level */
   /*-------------------------------- set dof numbers of the aggregates */
   /*----------------------------------------------- create prolongator */
   /*-------------------------------- restrict matrix to the next level */
   /*---------------------------------------------------- make printout */
   if (actintra->intra_rank==0)
   printf("level %d: size %d nnz %d\n",i,mlprecond.level[i].csr->numeq_total,
                                         mlprecond.level[i].csr->ja.fdim);
}
/*------------------------ make level numlev-1 (coarse grid) separately */
/*------------------------------------------------------- make printout */
if (actintra->intra_rank==0)
printf("level %d: size %d nnz %d\n",i,mlprecond.level[i].csr->numeq_total,
                                      mlprecond.level[i].csr->ja.fdim);

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_init */


/*!---------------------------------------------------------------------
\brief prolonge z down the finer level                                              

<pre>                                                        m.gee 11/02 
z = P * zc
</pre>
\param zc         double*     (o)   residuum in the dimension of the next level                   
\param z          double*     (i)   residuum on the active level
\param P          DBCSR*      (i)   the prolongator matrix
\param coarsecsr  DBCSR*      (i)   the fine level stiffness matrix
\param actintra   INTRA*      (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_prolongz(double *zc, double *z, DBCSR *P, DBCSR *coarsecsr, 
                            INTRA *actintra)
{
int        i,j,n,m,counter;
int        actrow,actcol,colstart,colend,index,owner,flag,tag,length;
int        myrank,nproc;
int       *update,*ia,*ja,numeq;
double    *a;
int        ineeds[MAXPROC][MAXPROC],ineedr[MAXPROC][MAXPROC];
int        nsend=0,nrecv=0;
ARRAY      irecv_a,drecv_a;
int       *irecv;
double    *drecv;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_prolongz");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
numeq  = P->numeq;
update = P->update.a.iv;
ia     = P->ia.a.iv;
ja     = P->ja.a.iv;
a      = P->a.a.dv;
/*----------------------------------------------------------------------*/
for (i=0; i<P->numeq; i++) z[i] = 0.0;

/*================================ make interproc multiplication part I */
/*======================================================================*/



/*============================== make the local piece of multiplication */
for (i=0; i<numeq; i++)
{
   actrow   = update[i];
   colstart = ia[i];
   colend   = ia[i+1];
   for (j=colstart; j<colend; j++)
   {
      actcol = ja[j];
      owner = mlpcg_getowner(actcol,coarsecsr->owner,nproc);
      if (owner != myrank)
         continue;
      index = mlpcg_getindex(actcol,coarsecsr->update.a.iv,coarsecsr->numeq);
      if (index==-1) dserror("Cannot find local dof");
      z[i] += a[j] * zc[index];
   }
}
/*======================================================================*/

/*=============================== make interproc multiplication part II */
/*----------------------------- allocate a guess of the receive buffers */
/*======================================================================*/




/*============================================================= tidy up */
/*======================================================================*/

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_prolongz */



/*!---------------------------------------------------------------------
\brief restrict the residuum to the next level                                              

<pre>                                                        m.gee 11/02 
rc = Ptransposed * r
</pre>
\param rc         double*     (o)   residuum in the dimension of the next level                   
\param r          double*     (i)   residuum on the active level
\param P          DBCSR*      (i)   the prolongator matrix
\param coarsecsr  DBCSR*      (i)   the coarse level stiffness matrix
\param actintra   INTRA*      (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_restrictr(double *rc, double *r, DBCSR *P, DBCSR *coarsecsr, 
                            INTRA *actintra)
{
int        i,j,n,m,counter;
int        myrank,nproc;
int        index,flag;
int        sender,tag,length;
int        numeq,*update,*ia,*ja;
double    *a;
int        actrow,actcol,colstart,colend,owner,rowindex;
int        cindex;
int        sendtos[MAXPROC][MAXPROC],sendtor[MAXPROC][MAXPROC];
int        nsend=0,nrecv=0,sendsize;
ARRAY      send_rindex_a,send_val_a,recv_rindex_a,recv_val_a;
int       *send_rindex,*recv_rindex;
double    *send_val,*recv_val;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_restrictr");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*----------------- if the column version of P doesn't exist, create it */
numeq       = P->numeq;
update      = P->update.a.iv;
ia          = P->ia.a.iv;
ja          = P->ja.a.iv;
a           = P->a.a.dv;
/*--------------------------- make the matvec init of the coarse matrix */
mlpcg_matvec_init(coarsecsr,actintra);
/*----------------------------------------------------------------------*/
for (i=0; i<coarsecsr->numeq; i++) rc[i] = 0.0;
/*---------------------- loop my interproc columns, and check the owner */
/*--------------------------------------- allocate guess of sendbuffers */
/*=========================make the local multiplication in daxpy style */
for (i=0; i<numeq; i++) /* loop rows of P */
{
   actrow   = update[i];
   colstart = ia[i];
   colend   = ia[i+1];
   for (j=colstart; j<colend; j++)
   {
      actcol = ja[j];
      owner  = mlpcg_getowner(actcol,coarsecsr->owner,nproc);
      if (owner != myrank)
         continue;
      cindex = mlpcg_getindex(actcol,coarsecsr->update.a.iv,coarsecsr->numeq);
      if (cindex==-1) dserror("Cannot find local dof");
      rc[cindex] += a[j] * r[i];
   }
}
/*================= make the interprocessor part of the multiplication */
/*==============================================================tidy up */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_restrictr */






/*!---------------------------------------------------------------------
\brief restrict the system matrix to the next level                                              

<pre>                                                        m.gee 9/02 

</pre>
\param actlev      MLLEVEL*     (i)   the active level in the ml-precond.                   
\param nextlev     MLLEVEL*     (i/o) the next coarser level.
\param actintra   INTRA*        (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_restrictK(MLLEVEL  *actlev, MLLEVEL *nextlev, INTRA *actintra)
{
int        i,j,counter;
int        myrank,nproc;
int        startdof,enddof;
int        numeq_total,nnz_guess;
DBCSR     *work=NULL;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_restrictK");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*------------ get the row numbers of the nextlev->csr stiffness matrix */
/*-------- these correspond to the columns of the prolongator actlev->P */
/* the blocks of the prolongator actlev->P are rowblocks (correspond to */
/* the rows and columns of actlev->csr) so we have to compute the column */
/* blocks from the aggregates themselfs.As the BDCSR is continous on proc*/
/* it is o.k. to take first and last dof from aggregates                */
startdof = actlev->agg[0].dof[0];
enddof   = actlev->agg[actlev->nagg-1].dof[actlev->agg[actlev->nagg-1].numdf-1];
/*------------- now we have to collect the unknown total number of dofs */
if (nextlev->csr->numeq_total==0)
{
   counter=0;
   for (i=0; i<actlev->nagg; i++) 
      counter += actlev->agg[i].numdf;
   numeq_total=counter;
}
/*------- as guess for the nnz serves the bigger finer grid matrix size */
nnz_guess = actlev->csr->nnz;
/*---------------- set numeq and numeq_total in the nextlev->csr matrix */
/*------------------------------------------ now we can open the matrix */
/*--- note that the variable nextlev->csr->firscoupledof is not yet set */
mlpcg_csr_open(nextlev->csr,startdof,enddof,numeq_total,nnz_guess,actintra);
/*----------------------------------------------------------------------*/
/* 
   for the restriction, we need a working csr matrix of dimensions      
   nrow        = nrow of nextlev->csr
   ncol        = nrow of actlev->csr
   numeq_total = numeq_total of nextlev->csr
   The size can be estimated to nnz_guess from actlev->csr
   The splitting of the dofs on the processores correponds to nextlev->csr
*/
work = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
if (!work) dserror("Allocation of memory failed");
mlpcg_csr_open(work,startdof,enddof,numeq_total,nnz_guess,actintra);
/*------------------------------------- we can now make the restriction */
mlpcg_precond_PtKP(actlev->P,actlev->csr,nextlev->csr,work,actlev->agg,
                   actlev->nagg,actintra);
/*------------------------------------ destroy the temporary csr matrix */
mlpcg_csr_destroy(work);
work = CCAFREE(work);
/*----------------------------------------------------------------------*/
/* 
  due to boundary conditions there may be empty columns in the prolongator.
  This leads to zero rows and columns and zero diagonal element in the
  coarse grid stiffness matrices. Without influence on the physical behaviour,
  a ONE is set on the diagonal of such rows 
*/
mlpcg_precond_checkdirich(nextlev->csr,actintra);
/*--------------------------- close the stiffness matrix on coarse grid */
mlpcg_csr_close(nextlev->csr);
/*----------------------- check for the first coupled row in the matrix */
mlpcg_precond_check_fcd(nextlev->csr,actintra);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_restrictK */




/*!---------------------------------------------------------------------
\brief check for the first coupled row in the matrix                                             

<pre>                                                        m.gee 10/02 
build the owner array inside csr and check for the first coupled 
row in the matrix
</pre>
\param matrix     DBCSR* (i/o) the csr matrix to be ckecked
\param actintra   INTRA* (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_check_fcd(DBCSR *matrix, INTRA *actintra)
{
int           i,j,n,actrow,actcol,colstart,colend;
int           myrank,nproc;
int           foundit,owner;
int           numeq,*update,*ia,*ja;
double       *a;
double        one=1.0;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_check_fcd");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
numeq  = matrix->numeq;
update = matrix->update.a.iv;
ia     = matrix->ia.a.iv;
ja     = matrix->ja.a.iv;
a      = matrix->a.a.dv;
/*----------------------------------------------------------------------*/
matrix->owner[myrank][0] = update[0];
matrix->owner[myrank][1] = update[numeq-1];
/*----------------------------------------------------------------------*/
owner = myrank;
for (i=0; i<numeq; i++)
{
   actrow   = update[i];
   colstart = ia[i];
   colend   = ia[i+1];
   for (j=colstart; j<colend; j++)
   {
      actcol = ja[j];
      owner = mlpcg_getowner(actcol,matrix->owner,nproc);
      if (owner != myrank)
      {
         matrix->firstcoupledof = actrow;
         goto exit;
      }
   }
}
/*----------------------------------------------------------------------*/
exit:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_check_fcd */

/*!---------------------------------------------------------------------
\brief check for zero diagonal entries                                             

<pre>                                                        m.gee 10/02 
  due to boundary conditions there may be empty columns in the prolongator.
  This leads to zero rows and columns and zero diagonal element in the
  coarse grid stiffness matrices. Without influence on the physical behaviour,
  a ONE is set on the diagonal of such rows 
</pre>
\param matrix     DBCSR* (i/o) the csr matrix to be ckecked
\param actintra   INTRA* (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_checkdirich(DBCSR *matrix, INTRA *actintra)
{
int           i,j,actrow,actcol,colstart,colend;
int           foundit;
int           numeq,*update,*ia,*ja;
double       *a;
double        one=1.0;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_checkdirich");
#endif
/*----------------------------------------------------------------------*/
numeq  = matrix->numeq;
update = matrix->update.a.iv;
ia     = matrix->ia.a.iv;
ja     = matrix->ja.a.iv;
a      = matrix->a.a.dv;
/*----------------------------------------------------------------------*/
for (i=0; i<numeq; i++)
{
   foundit  = 0;
   actrow   = update[i];
   colstart = ia[i];
   colend   = ia[i+1];
   for (j=colstart; j<colend; j++)
   {
      actcol = ja[j];
      if (actcol == actrow)
      {
         foundit=1;
         if (FABS(a[j])<EPS14) 
         a[j] = 1.0;
         break;
      }
   }
   if (!foundit)
   {
      mlpcg_csr_addentry(matrix,one,actrow,actcol,actintra);
      numeq  = matrix->numeq;
      update = matrix->update.a.iv;
      ia     = matrix->ia.a.iv;
      ja     = matrix->ja.a.iv;
      a      = matrix->a.a.dv;
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_checkdirich */


/*!---------------------------------------------------------------------
\brief get directors from the shell8 element for rbm's                                              

<pre>                                                        m.gee 9/02 

</pre>
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_getdirs(void)
{
int           i,j,k,l,counter;
PARTDISCRET  *actpdis;
NODE         *actnode;
ELEMENT      *actele;
double      **director;
double        fac;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_getdirs");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------- check presence of shell8 elements */
#ifndef D_SHELL8
dserror("SHELL8 not compiled in! MLPCG is a SHELL8 solver!");
#endif
/*-------------------------------------------- get pointer to partition */   
actpdis = mlprecond.partdis;
/*------------------------- allocate pointers to directors in mlprecond */
if (mlprecond.director.Typ==cca_XX)
{
   director = amdef("direct",&(mlprecond.director),actpdis->numnp,3,"DA");
   mlprecond.node      = (NODE**)CCAMALLOC((actpdis->numnp)*sizeof(NODE*));
}   
if (!(mlprecond.node)) 
dserror("Allocation of memory failed");
/*------------------- loop all nodes in actpdis and get there directors */
/* 
IMPORTANT NOTE:
in the nonlinear dynamic or static case, it could be advisable to use the
deformed configuration (total lagrange!) to create the rigid body modes !
*/
counter=0;
for (i=0; i<actpdis->numnp; i++)
{
   actnode = actpdis->node[i];
   mlprecond.node[counter] = actnode;
   /* get element connected to this node */
   actele = actnode->element[0];
   /* loop nodes on element to find right one */
   for (j=0; j<actele->numnp; j++)
      if (actele->node[j] == actnode) break;
   dsassert(j!=actele->numnp,"Couldn't find node in element");
#ifdef D_SHELL8
   /* h/2 */
   fac = actele->e.s8->thick_node.a.dv[j] * actele->e.s8->sdc / 2.0 ;
   /* averaged director of unit length bischoff style */
   director[counter][0] = actele->e.s8->a3ref.a.da[0][j]*fac;
   director[counter][1] = actele->e.s8->a3ref.a.da[1][j]*fac;
   director[counter][2] = actele->e.s8->a3ref.a.da[2][j]*fac;
#endif
   counter++;
}
dsassert(counter==actpdis->numnp,"Number of nodes wrong");
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_getdirs */

/*! @} (documentation module close)*/

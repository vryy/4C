#include "../headers/standardtypes.h"
#include "../headers/solution.h"

/*----------------------------------------------------------------------*
 |                                                           m.gee 03/02|
 | prototypes of functions only visible in this file                    |
 *----------------------------------------------------------------------*/
static void solserv_cp_rc_ptrmask(INTRA *actintra, RC_PTR *from, RC_PTR *to);
static void solserv_cp_ucchbmask(UCCHB *from, UCCHB *to);
static void solserv_cp_skymask(SKYMATRIX *from, SKYMATRIX *to);
static void solserv_cp_msrmask(AZ_ARRAY_MSR *from, AZ_ARRAY_MSR *to);
static void solserv_cp_densemask(DENSE *from, DENSE *to);
static void solserv_matvec_msr(INTRA *actintra,AZ_ARRAY_MSR *msr,double *work1,double *work2);
static void solserv_matvec_sky(INTRA *actintra,SKYMATRIX *sky,double *work1,double *work2);
static void solserv_matvec_dense(INTRA *actintra,DENSE *dense,double *work1,double *work2);

/*----------------------------------------------------------------------*
 |  make A = A * factor                                      m.gee 02/02|
 | SPARSE_TYP *Atyp (i)        type of sparse matrix                    |
 | SPARSE_ARRAY *A  (i/o)      sparse matrix                            |
 | double factor    (i)        factor                                   |
 *----------------------------------------------------------------------*/
void solserv_scal_mat(SPARSE_TYP *Atyp,SPARSE_ARRAY *A,double factor)
{
int      i;

#ifdef DEBUG 
dstrc_enter("solserv_scal_mat");
#endif
/*----------------------------------------------------------------------*/
switch (*Atyp)
{
case mds:
   dserror("not implemented for MLIB yet");
break;
case msr:
   amscal(&(A->msr->val),&factor);
break;
case parcsr:
   dserror("not implemented for HYPRE yet");
break;
case ucchb:
   amscal(&(A->ucchb->a),&factor);
break;
case dense:
   amscal(&(A->dense->A),&factor);
break;
case rc_ptr:
   amscal(&(A->rc_ptr->A_loc),&factor);
break;
case skymatrix:
   amscal(&(A->sky->A),&factor);
break;
case sparse_none:
   dserror("Unknown typ of sparse distributed system matrix");
break;
default:
   dserror("Unknown typ of sparse distributed system matrix");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_scal_mat */


/*----------------------------------------------------------------------*
 |  make A = A + B * factor                                  m.gee 02/02|
 | INTRA *actintra  (i)     intra-communicator the matrices live on     |
 | SPARSE_TYP *Atyp (i)     type of sparse matrix                       |
 | SPARSE_ARRAY *A  (i/o)   sparse matrix                               |
 | SPARSE_TYP *Btyp (i)     type of sparse matrix                       |
 | SPARSE_ARRAY *B  (i)     sparse matrix                               |
 | double factor    (i)     factor                                      |
 *----------------------------------------------------------------------*/
void solserv_add_mat(INTRA *actintra,
                     SPARSE_TYP *Atyp,
                     SPARSE_ARRAY *A,
                     SPARSE_TYP *Btyp,
                     SPARSE_ARRAY *B,
                     double factor)
{
int      i;

#ifdef DEBUG 
dstrc_enter("solserv_add_mat");
#endif
/*----------------------------------------------------------------------*/
if (*Atyp != *Btyp) dserror("Incompatible types of sparse matrices");
switch (*Atyp)
{
case mds:
   dserror("not implemented for MLIB yet");
break;
case msr:
   amadd(&(A->msr->val),&(B->msr->val),factor,0);
break;
case parcsr:
   dserror("not implemented for HYPRE yet");
break;
case ucchb:
   amadd(&(A->ucchb->a),&(B->ucchb->a),factor,0);
break;
case dense:
   amadd(&(A->dense->A),&(B->dense->A),factor,0);
break;
case rc_ptr:
   amadd(&(A->rc_ptr->A_loc),&(B->rc_ptr->A_loc),factor,0);
break;
case skymatrix:
   amadd(&(A->sky->A),&(B->sky->A),factor,0);
break;
case sparse_none:
   dserror("Unknown typ of sparse distributed system matrix");
break;
default:
   dserror("Unknown typ of sparse distributed system matrix");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_add_mat */


/*----------------------------------------------------------------------*
 |  get dimensions of sparse matrix                          m.gee 02/02|
 | SPARSE_ARRAY mat  (i)     sparse matrix                              |
 | SPARSE_TYP mattyp (i)     type of sparse matrix                      |
 | int *numeq        (o)     proc-local dimension of sparse matrix      |
 | int *numeq_total  (o)     global dimension of sparse matrix          |
 *----------------------------------------------------------------------*/
void solserv_getmatdims(SPARSE_ARRAY mat,SPARSE_TYP mattyp,
                        int *numeq, int *numeq_total)
{
#ifdef DEBUG 
dstrc_enter("solserv_getmatdims");
#endif
/*----------------------------------------------------------------------*/
switch (mattyp)
{
case mds:
   *numeq       = mat.mds->numeq;
   *numeq_total = *numeq;
break;
case msr:
   *numeq       = mat.msr->numeq;
   *numeq_total = mat.msr->numeq_total;
break;
case parcsr:
   *numeq       = mat.parcsr->numeq;
   *numeq_total = mat.parcsr->numeq_total;
break;
case ucchb:
   *numeq       = mat.ucchb->numeq;
   *numeq_total = mat.ucchb->numeq_total;
break;
case dense:
   *numeq       = mat.dense->numeq;
   *numeq_total = mat.dense->numeq_total;
break;
case rc_ptr:
   *numeq       = mat.rc_ptr->numeq;
   *numeq_total = mat.rc_ptr->numeq_total;
break;
case skymatrix:
   *numeq       = mat.sky->numeq;
   *numeq_total = mat.sky->numeq_total;
break;
case sparse_none:
   dserror("Unknown typ of sparse distributed system matrix");
break;
default:
   dserror("Unknown typ of sparse distributed system matrix");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_getmatdims */



/*----------------------------------------------------------------------*
 |  init a distributed matrix to zero - collective call !    m.gee 10/01|
 | INTRA *actintra  (i)     intra-communicator the matrices live on     |
 | SPARSE_ARRAY *mat  (i/o)    sparse matrix                            |
 | SPARSE_TYP *mattyp (i)    type of sparse matrix                      |
 *----------------------------------------------------------------------*/
void solserv_zero_mat(INTRA *actintra, SPARSE_ARRAY *mat,SPARSE_TYP *mattyp)
{
int                  imyrank;
int                  inprocs;
int                  ilower,iupper,jlower,jupper;
int                  err=1;
#ifdef DEBUG 
dstrc_enter("solserv_zero_mat");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
switch (*mattyp)
{
case mds:
   solver_mlib(NULL,NULL,mat->mds,NULL,NULL,2);
break;
case msr:
   amzero(&(mat->msr->val));
   mat->msr->is_factored=0;
break;
case parcsr:/*---- this stupid package does not have a zero function!!!!*/
#ifdef HYPRE_PACKAGE
   err=HYPRE_IJMatrixDestroy(mat->parcsr->ij_matrix);
   if (err) dserror("Cannot destroy PARCSR matrix");
   ilower = jlower = mat->parcsr->perm.a.ia[imyrank][0];
   iupper = jupper = mat->parcsr->perm.a.ia[imyrank][mat->parcsr->perm_sizes.a.iv[imyrank]-1];
   err=HYPRE_IJMatrixCreate(actintra->MPI_INTRA_COMM,ilower,iupper,jlower,jupper,&(mat->parcsr->ij_matrix));
   if (err) dserror("Cannot Create PARCSR matrix");
   err=HYPRE_IJMatrixSetObjectType(mat->parcsr->ij_matrix,HYPRE_PARCSR);
   if (err) dserror("Cannot Set PARCSR matrix object type");
   err=HYPRE_IJMatrixInitialize(mat->parcsr->ij_matrix);
   if (err) dserror("Cannot init PARCSR matrix");
#endif
   mat->parcsr->is_factored=0;
break;
case ucchb:
   amzero(&(mat->ucchb->a));
   mat->ucchb->is_factored=0;
break;
case dense:
   amzero(&(mat->dense->A));
   mat->dense->is_factored=0;
break;
case rc_ptr:
   amzero(&(mat->rc_ptr->A_loc));
   mat->rc_ptr->is_factored=0;
break;
case skymatrix:
   amzero(&(mat->sky->A));
   mat->sky->is_factored=0;
break;
case sparse_none:
   dserror("Unknown typ of sparse distributed system matrix");
break;
default:
   dserror("Unknown typ of sparse distributed system matrix");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_zero_mat */



/*----------------------------------------------------------------------*
 |  copies the sparsity mask of a system matrix              m.gee 02/02|
 |  for the new sparsity mask, a suitable structure is allocated        |
 | INTRA        *actintra (i)  intra-communicator the matrices live on  |
 | SPARSE_TYP   *typfrom  (i)  type of sparsity mask to be copied       |
 | SPARSE_ARRAY *matfrom  (i)  sparsity mask to be copied               |
 | SPARSE_TYP   *typto    (o)  type of sparsity mask to be copied to    |
 | SPARSE_ARRAY *matto    (o)  sparsity mask to be allocated and copied |
 *----------------------------------------------------------------------*/
void solserv_alloc_cp_sparsemask(INTRA        *actintra, 
                                 SPARSE_TYP   *typfrom,
                                 SPARSE_ARRAY *matfrom,
                                 SPARSE_TYP   *typto,
                                 SPARSE_ARRAY *matto)
{
#ifdef DEBUG 
dstrc_enter("solserv_alloc_cp_sparsemask");
#endif
/*----------------------------------------------------------------------*/
   *typto = *typfrom;
/*----------------------------------------------------------------------*/
switch (*typfrom)
{
case mds:
   matto->mds = (ML_ARRAY_MDS*)CALLOC(1,sizeof(ML_ARRAY_MDS));
   if (!matto->mds) dserror("Allocation of memory failed");
   dserror("Copy of msd matrix not yet implemented ");
break;
case msr:
   matto->msr = (AZ_ARRAY_MSR*)CALLOC(1,sizeof(AZ_ARRAY_MSR));
   if (!matto->msr) dserror("Allocation of memory failed");
   solserv_cp_msrmask(matfrom->msr,matto->msr);
break;
case parcsr:
   matto->parcsr = (H_PARCSR*)CALLOC(1,sizeof(H_PARCSR));
   if (!matto->parcsr) dserror("Allocation of memory failed");
   dserror("Copy of parcsr matrix not yet implemented ");
break;
case ucchb:
   matto->ucchb = (UCCHB*)CALLOC(1,sizeof(UCCHB));
   if (!matto->ucchb) dserror("Allocation of memory failed");
   solserv_cp_ucchbmask(matfrom->ucchb,matto->ucchb);
break;
case dense:
   matto->dense = (DENSE*)CALLOC(1,sizeof(DENSE));
   if (!matto->dense) dserror("Allocation of memory failed");
   solserv_cp_densemask(matfrom->dense,matto->dense);
break;
case rc_ptr:
   matto->rc_ptr = (RC_PTR*)CALLOC(1,sizeof(RC_PTR));
   if (!matto->rc_ptr) dserror("Allocation of memory failed");
   solserv_cp_rc_ptrmask(actintra,matfrom->rc_ptr,matto->rc_ptr);
break;
case skymatrix:
   matto->sky = (SKYMATRIX*)CALLOC(1,sizeof(SKYMATRIX));
   if (!matto->sky) dserror("Allocation of memory failed");
   solserv_cp_skymask(matfrom->sky,matto->sky);
break;
case sparse_none:
   dserror("Unknown typ of sparse distributed system matrix");
break;
default:
   dserror("Unknown typ of sparse distributed system matrix");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_alloc_cp_sparsemask */


/*----------------------------------------------------------------------*
 |  copies the sparsity mask of a rc_ptr matrix              m.gee 02/02|
 |  for the new sparsity mask, all memory is allocated                  |
 |  called by solserv_alloc_cp_sparsemask only!                         |
 *----------------------------------------------------------------------*/
static void solserv_cp_rc_ptrmask(INTRA *actintra, RC_PTR *from, RC_PTR *to)
{
int i;
#ifdef DEBUG 
dstrc_enter("solserv_cp_rc_ptrmask");
#endif
/*----------------------------------------------------------------------*/
/* copy all information, which is directly included in the structure */
*to = *from;

/* alloccopy update */
am_alloc_copy(&(from->update),&(to->update));

/* alloccopy irn_loc */
am_alloc_copy(&(from->irn_loc),&(to->irn_loc));

/* alloccopy jcn_loc */
am_alloc_copy(&(from->jcn_loc),&(to->jcn_loc));

/* alloccopy A_loc */
am_alloc_copy(&(from->A_loc),&(to->A_loc));

/* alloccopy rowptr */
am_alloc_copy(&(from->rowptr),&(to->rowptr));

/* alloccopy irn_glob and jcn_glob on imyrank=0 */
if (actintra->intra_rank==0)
{
   am_alloc_copy(&(from->irn_glob),&(to->irn_glob));
   am_alloc_copy(&(from->jcn_glob),&(to->jcn_glob));
}

/* backups are allocated in solver_mumps */

/* numcoupsend,numcouprecv,couple_d_send,couple_i_send,couple_d_recv,couple_i_recv
   are allocated if necessary in the assembly init phase */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_cp_rc_ptrmask */


/*----------------------------------------------------------------------*
 |  copies the sparsity mask of a ucchb matrix               m.gee 02/02|
 |  for the new sparsity mask, all memory is allocated                  |
 |  called by solserv_alloc_cp_sparsemask only!                         |
 *----------------------------------------------------------------------*/
static void solserv_cp_ucchbmask(UCCHB *from, UCCHB *to)
{
int i;
#ifdef DEBUG 
dstrc_enter("solserv_cp_ucchbmask");
#endif
/*----------------------------------------------------------------------*/
/* copy all information, which is directly included in the structure */
*to = *from;

/* alloccopy update */
am_alloc_copy(&(from->update),&(to->update));

/* alloccopy a */
am_alloc_copy(&(from->a),&(to->a));

/* alloccopy asub */
am_alloc_copy(&(from->asub),&(to->asub));

/* alloccopy xa */
am_alloc_copy(&(from->xa),&(to->xa));


/* backups are allocated in solver_superlu */

/* numcoupsend,numcouprecv,couple_d_send,couple_i_send,couple_d_recv,couple_i_recv
   are allocated if necessary in the assembly init phase */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_cp_ucchbmask */


/*----------------------------------------------------------------------*
 |  copies the sparsity mask of a skyline matrix             m.gee 02/02|
 |  for the new sparsity mask, all memory is allocated                  |
 |  called by solserv_alloc_cp_sparsemask only!                         |
 *----------------------------------------------------------------------*/
static void solserv_cp_skymask(SKYMATRIX *from, SKYMATRIX *to)
{
int i;
#ifdef DEBUG 
dstrc_enter("solserv_cp_skymask");
#endif
/*----------------------------------------------------------------------*/
/* copy all information, which is directly included in the structure */
*to = *from;

/* alloccopy update */
am_alloc_copy(&(from->update),&(to->update));

/* alloccopy maxa */
am_alloc_copy(&(from->maxa),&(to->maxa));

/* maxaf is allocated in solver_colsol init phase */

/* alloccopy A */
am_alloc_copy(&(from->A),&(to->A));

/* val_backup is allocated in solver_aztec */

/* numcoupsend,numcouprecv,couple_d_send,couple_i_send,couple_d_recv,couple_i_recv
   are allocated if necessary in the assembly init phase */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_cp_skymask */


/*----------------------------------------------------------------------*
 |  copies the sparsity mask of a msr matrix                 m.gee 02/02|
 |  for the new sparsity mask, all memory is allocated                  |
 |  called by solserv_alloc_cp_sparsemask only!                         |
 *----------------------------------------------------------------------*/
static void solserv_cp_msrmask(AZ_ARRAY_MSR *from, AZ_ARRAY_MSR *to)
{
int i;
#ifdef DEBUG 
dstrc_enter("solserv_cp_msrmask");
#endif
/*----------------------------------------------------------------------*/
/* copy all information, which is directly included in the structure */
*to = *from;

/* alloccopy update */
am_alloc_copy(&(from->update),&(to->update));

/* shift and bins are allocated by assembly routines */

/* alloccopy bindx */
am_alloc_copy(&(from->bindx),&(to->bindx));

/* bindx_backup is allocated in solver_aztec */

/* alloccopy val */
am_alloc_copy(&(from->val),&(to->val));

/* val_backup is allocated in solver_aztec */

/* numcoupsend,numcouprecv,couple_d_send,couple_i_send,couple_d_recv,couple_i_recv
   are allocated if necessary in the assembly init phase */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_cp_msrmask */


/*----------------------------------------------------------------------*
 |  copies the sparsity mask of a dense matrix               m.gee 02/02|
 |  for the new sparsity mask, all memory is allocated                  |
 |  called by solserv_alloc_cp_sparsemask only!                         |
 *----------------------------------------------------------------------*/
static void solserv_cp_densemask(DENSE *from, DENSE *to)
{
#ifdef DEBUG 
dstrc_enter("solserv_cp_densemask");
#endif
/*----------------------------------------------------------------------*/
/* copy all information, which is directly included in the structure */
*to = *from;

/* alloccopy update */
am_alloc_copy(&(from->update),&(to->update));

/* make A */
amdef(from->A.name,&(to->A),from->A.fdim,from->A.sdim,"DA");

/* ipiv, lwork and work are allocated in solver_lapack in the init phase */

/* couple_d_send,couple_i_send,couple_d_recv,couple_i_recv
   are allocated if necessary in the assembly init phase */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_cp_densemask */



/*----------------------------------------------------------------------*
 |  make matrix vector multiplication                        m.gee 02/02|
 | INTRA        *actintra (i)  intra-communicator the matrices live on  |
 | DIST_VECTOR  *result   (o)  result = mat * vec                       |
 | SPARSE_ARRAY *mat      (i)  sparse matrix                            |
 | SPARSE_TYP   *mattyp   (i)  type of sparse matrix                    |
 | DIST_VECTOR  *vec      (i)  vector to be multiplied with             |
 *----------------------------------------------------------------------*/
void solserv_sparsematvec(INTRA        *actintra,
                          DIST_VECTOR  *result,
                          SPARSE_ARRAY *mat,
                          SPARSE_TYP   *mattyp,
                          DIST_VECTOR  *vec)
{
static ARRAY   work1_a;
static double *work1;
static ARRAY   work2_a;
static double *work2;

#ifdef DEBUG 
dstrc_enter("solserv_sparsematvec");
#endif
/*----------------------------------------------------------------------*/
if (work1_a.Typ != DV) work1 = amdef("work1",&work1_a,vec->numeq_total,1,"DV");
if (work2_a.Typ != DV) work2 = amdef("work2",&work2_a,vec->numeq_total,1,"DV");
if (work1_a.fdim < vec->numeq_total)
{
           amdel(&work1_a);
   work1 = amdef("work1",&work1_a,vec->numeq_total,1,"DV");
           amdel(&work2_a);
   work2 = amdef("work2",&work2_a,vec->numeq_total,1,"DV");
}
/*----------------------------------------------------------------------*/
switch (*mattyp)
{
case mds:
   dserror("Matrix-Vector Product for MLIB not implemented");
break;
case msr:
   solserv_reddistvec(vec,mat,mattyp,work1,vec->numeq_total,actintra);
   solserv_matvec_msr(actintra,mat->msr,work1,work2);
   solserv_distribdistvec(result,mat,mattyp,work2,vec->numeq_total,actintra);
break;
case parcsr:
   solserv_reddistvec(vec,mat,mattyp,work1,vec->numeq_total,actintra);
   solserv_distribdistvec(result,mat,mattyp,work2,vec->numeq_total,actintra);
   dserror("Matrix-Vector Product for HYPRE not implemented");
break;
case ucchb:
   solserv_reddistvec(vec,mat,mattyp,work1,vec->numeq_total,actintra);
   solserv_distribdistvec(result,mat,mattyp,work2,vec->numeq_total,actintra);
   dserror("Matrix-Vector Product for SuperLU not implemented");
break;
case dense:
   solserv_reddistvec(vec,mat,mattyp,work1,vec->numeq_total,actintra);
   solserv_matvec_dense(actintra,mat->dense,work1,work2);
   solserv_distribdistvec(result,mat,mattyp,work2,vec->numeq_total,actintra);
break;
case rc_ptr:
   solserv_reddistvec(vec,mat,mattyp,work1,vec->numeq_total,actintra);
   solserv_distribdistvec(result,mat,mattyp,work2,vec->numeq_total,actintra);
   dserror("Matrix-Vector Product for MUMPS not implemented");
break;
case skymatrix:
   solserv_reddistvec(vec,mat,mattyp,work1,vec->numeq_total,actintra);
   solserv_matvec_sky(actintra,mat->sky,work1,work2);
   solserv_distribdistvec(result,mat,mattyp,work2,vec->numeq_total,actintra);
break;
case sparse_none:
   dserror("Unknown typ of sparse distributed system matrix");
break;
default:
   dserror("Unknown typ of sparse distributed system matrix");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_sparsematvec */




/*----------------------------------------------------------------------*
 |  make matrix vector multiplication with dense matrix      m.gee 02/02|
 |  called by solserv_sparsematvec only !                               |
 *----------------------------------------------------------------------*/
static void solserv_matvec_msr(INTRA        *actintra,
                        AZ_ARRAY_MSR *msr,
                        double       *work1,
                        double       *work2)
{
#ifdef AZTEC_PACKAGE
int         i,j,dof;
int         start,end,lenght,j_index;
int         myrank;
int         nprocs;
int         numeq;
int         numeq_total;
int        *update;
int        *bindx;
double     *val;

#ifdef DEBUG 
dstrc_enter("solserv_matvec_msr");
#endif
/*----------------------------------------------------------------------*/
myrank     = actintra->intra_rank;
nprocs     = actintra->intra_nprocs;
numeq_total= msr->numeq_total;
numeq      = msr->numeq;
update     = msr->update.a.iv;
bindx      = msr->bindx.a.iv;
val        = msr->val.a.dv;
/*------------------------------------------- now loop my own equations */
for (i=0; i<numeq; i++)
{
   dof = update[i];
   work2[dof]=0.0;
   /*----------- the dof number is dof, the place in val and bindx is i */
   /*---------------------------------------- the place in work1 is dof */
   /*------------------------------------- multiply main diagonal entry */
   work2[dof] += val[i] * work1[dof];
   /*------------------------------------ multiply off-diagonal entries */
   start  = bindx[i];
   end    = bindx[i+1];
   lenght = end - start;
   for (j=start; j<end; j++)
   {
      j_index = bindx[j];
      work2[dof] += val[j] * work1[j_index];
   }
}
/* every proc added his complete part to the full-sized vector work2.
   There is no need to alreduce this vector here, because it will anyway
   bre distributed to a DIS_VECTOR */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
# endif /* end of #ifdef AZTEC_PACKAGE */
return;
} /* end of solserv_matvec_msr */





/*----------------------------------------------------------------------*
 |  make matrix vector multiplication with skyline matrix    m.gee 02/02|
C     ******************************************************************
C     *  INSTITUT FUER BAUSTATIK  *  CREATED / CHANGED                 *
C     *  UNIVERSITAET  STUTTGART  *  04.05.92/           O. PETERSEN   *
C     ******************************************************************
         ported to C                                     m.gee
 |  called by solserv_sparsematvec only !                               |
 *----------------------------------------------------------------------*/
static void solserv_matvec_sky(INTRA        *actintra,
                        SKYMATRIX    *sky,
                        double       *work1,
                        double       *work2)
{
int     i,j,kl,ku,kdiff;
int     nrn,neq1,neq;
double *A;
int    *maxa;

#ifdef DEBUG 
dstrc_enter("solserv_matvec_sky");
#endif
/*---------------------- the dense matrix is redundant on all processes */
nrn  = sky->A.fdim;
neq  = sky->numeq_total;
neq1 = neq + 1;
A    = sky->A.a.dv;
maxa = sky->maxa.a.iv;
/*----------------------------------------------------------------------*/
for (i=0; i<neq; i++)
{
   work2[i]=0.0;
   kl    = maxa[i];
   ku    = maxa[i+1];
   kdiff = ku-kl;
   work2[i] += A[kl] * work1[i];
   for (j=1; j<kdiff; j++)
   {
      work2[i-j] += A[kl+j] * work1[i];
      work2[i]   += A[kl+j] * work1[i-j];
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_matvec_sky */



/*----------------------------------------------------------------------*
 |  make matrix vector multiplication with dense matrix      m.gee 02/02|
 |  called by solserv_sparsematvec only !                               |
 *----------------------------------------------------------------------*/
static void solserv_matvec_dense(INTRA        *actintra,
                          DENSE        *dense,
                          double       *work1,
                          double       *work2)
{
int      i,j,k;
int      I,J;
double   sum;
double **A;
#ifdef DEBUG 
dstrc_enter("solserv_matvec_dense");
#endif
/*---------------------- the dense matrix is redundant on all processes */
I = dense->numeq_total;
J = I;
A = dense->A.a.da;
/*----------------------------------------------------------------------*/
for (i=0; i<I; i++)
{
   sum = 0.0;
   for (j=0; j<J; j++)
   {
         sum += A[i][j] * work1[j];
   }
   work2[i] = sum;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_matvec_dense */

/*!----------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"

/*----------------------------------------------------------------------*
 |                                                           m.gee 03/02|
 | prototypes of functions only visible in this file                    |
 *----------------------------------------------------------------------*/
static void solserv_cp_rc_ptrmask(INTRA *actintra, RC_PTR *from, RC_PTR *to);
static void solserv_cp_ccfmask(INTRA *actintra, CCF *from, CCF *to);
static void solserv_cp_ucchbmask(UCCHB *from, UCCHB *to);
static void solserv_cp_skymask(SKYMATRIX *from, SKYMATRIX *to);
static void solserv_cp_msrmask(AZ_ARRAY_MSR *from, AZ_ARRAY_MSR *to);
static void solserv_cp_densemask(DENSE *from, DENSE *to);
static void solserv_cp_spomask(SPOOLMAT *from, SPOOLMAT *to);
static void solserv_cp_bdcsrmask(DBCSR *from, DBCSR *to);
static void solserv_matvec_rc_ptr(INTRA  *actintra,RC_PTR *rcptr,DOUBLE *work1,DOUBLE *work2);
static void solserv_matvec_spo(INTRA *actintra,SPOOLMAT *spo,DOUBLE *work1,DOUBLE *work2);
static void solserv_matvec_ccf(INTRA  *actintra,CCF *ccf,DOUBLE *work1,DOUBLE *work2);
#if 0
static void solserv_matvec_msr(INTRA *actintra,AZ_ARRAY_MSR *msr,DOUBLE *work1,DOUBLE *work2);
#endif
static void solserv_matvec_sky(INTRA *actintra,SKYMATRIX *sky,DOUBLE *work1,DOUBLE *work2);
static void solserv_matvec_dense(INTRA *actintra,DENSE *dense,DOUBLE *work1,DOUBLE *work2);
static void oll_matvec(INTRA *actintra, OLL *oll, DOUBLE *work1, DOUBLE *work2);
/*----------------------------------------------------------------------*
 |  make A = A * factor                                      m.gee 02/02|
 | SPARSE_TYP *Atyp (i)        type of sparse matrix                    |
 | SPARSE_ARRAY *A  (i/o)      sparse matrix                            |
 | DOUBLE factor    (i)        factor                                   |
 *----------------------------------------------------------------------*/
void solserv_scal_mat(SPARSE_TYP *Atyp,SPARSE_ARRAY *A,DOUBLE factor)
{

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
case spoolmatrix:
   amscal(&(A->spo->A_loc),&factor);
break;
case bdcsr:
   amscal(&(A->bdcsr->a),&factor);
break;
case ccf:
   amscal(&(A->ccf->Ax),&factor);
break;
case skymatrix:
   amscal(&(A->sky->A),&factor);
break;
case oll:
   oll_scal(A->oll,factor);
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
 | DOUBLE factor    (i)     factor                                      |
 *----------------------------------------------------------------------*/
void solserv_add_mat(INTRA *actintra,
                     SPARSE_TYP *Atyp,
                     SPARSE_ARRAY *A,
                     SPARSE_TYP *Btyp,
                     SPARSE_ARRAY *B,
                     DOUBLE factor)
{

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
case spoolmatrix:
#ifdef D_CONTACT
   add_spooles_matrix(A->spo,B->spo,factor,0,actintra);
#else
   amadd(&(A->spo->A_loc),&(B->spo->A_loc),factor,0);
#endif
break;
case bdcsr:
   amadd(&(A->bdcsr->a),&(B->bdcsr->a),factor,0);
break;
case ccf:
   amadd(&(A->ccf->Ax),&(B->ccf->Ax),factor,0);
break;
case skymatrix:
   amadd(&(A->sky->A),&(B->sky->A),factor,0);
   break;
case oll:
   if(B->oll->is_masked==0) break;
   if(A->oll->is_masked==0) 
   {
     oll_copy(B->oll,A->oll);
     oll_scal(B->oll,factor);
   }
   else
   {
     oll_add(A->oll,B->oll,factor);
   }
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
 | INT *numeq        (o)     proc-local dimension of sparse matrix      |
 | INT *numeq_total  (o)     global dimension of sparse matrix          |
 *----------------------------------------------------------------------*/
void solserv_getmatdims(SPARSE_ARRAY* mat,SPARSE_TYP mattyp,
                        INT *numeq, INT *numeq_total)
{
#ifdef DEBUG 
dstrc_enter("solserv_getmatdims");
#endif
/*----------------------------------------------------------------------*/
switch (mattyp)
{
case mds:
   *numeq       = mat->mds->numeq;
   *numeq_total = *numeq;
break;
case msr:
   *numeq       = mat->msr->numeq;
   *numeq_total = mat->msr->numeq_total;
break;
case parcsr:
   *numeq       = mat->parcsr->numeq;
   *numeq_total = mat->parcsr->numeq_total;
break;
case ucchb:
   *numeq       = mat->ucchb->numeq;
   *numeq_total = mat->ucchb->numeq_total;
break;
case dense:
   *numeq       = mat->dense->numeq;
   *numeq_total = mat->dense->numeq_total;
break;
case rc_ptr:
   *numeq       = mat->rc_ptr->numeq;
   *numeq_total = mat->rc_ptr->numeq_total;
break;
case ccf:
   *numeq       = mat->ccf->numeq;
   *numeq_total = mat->ccf->numeq_total;
break;
case skymatrix:
   *numeq       = mat->sky->numeq;
   *numeq_total = mat->sky->numeq_total;
break;
case spoolmatrix:
   *numeq       = mat->spo->numeq;
   *numeq_total = mat->spo->numeq_total;
break;
case bdcsr:
   *numeq       = mat->bdcsr->numeq;
   *numeq_total = mat->bdcsr->numeq_total;
break;
case oll:
   *numeq       = mat->oll->numeq;
   *numeq_total = mat->oll->numeq_total;
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

INT                  imyrank;
INT                  inprocs;

#ifdef HYPRE_PACKAGE
INT                  ilower,iupper,jlower,jupper;
INT                  err=1;
#endif

#ifdef D_CONTACT
INT                  minusone=-1;
#endif

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
   mat->msr->is_transformed=0;
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
case ccf:
   amzero(&(mat->ccf->Ax));
   mat->ccf->is_factored=0;
break;
case skymatrix:
   amzero(&(mat->sky->A));
   mat->sky->is_factored=0;
break;
case bdcsr:
   amzero(&(mat->bdcsr->a));
   mat->bdcsr->is_factored=0;
break;
case spoolmatrix:
   amzero(&(mat->spo->A_loc));
#ifdef D_CONTACT
   aminit(&(mat->spo->irn_loc),&minusone);
   aminit(&(mat->spo->jcn_loc),&minusone);
#endif
   mat->spo->is_factored=0;
#ifdef SPOOLES_PACKAGE
   if (mat->spo->ncall > 0)
   {
      FrontMtx_free(mat->spo->frontmtx);
      InpMtx_free(mat->spo->newA);
      DenseMtx_free(mat->spo->newY);
/*      DenseMtx_free(mat->spo->mtxX);*/
      ETree_free(mat->spo->frontETree);
      SubMtxManager_free(mat->spo->mtxmanager);
      IV_free(mat->spo->newToOldIV);
      IV_free(mat->spo->oldToNewIV);
      IV_free(mat->spo->ownersIV);
      IV_free(mat->spo->vtxmapIV);
      IV_free(mat->spo->ownedColumnsIV);
      SolveMap_free(mat->spo->solvemap);
      IVL_free(mat->spo->symbfacIVL);
   }
#endif
break;
case oll:
   oll_zero(mat->oll);
   /*solserv_zero_mat(actintra,&(mat->oll->sysarray[0]),&(mat->oll->sysarray_typ[0]));*/
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
   matto->mds = (ML_ARRAY_MDS*)CCACALLOC(1,sizeof(ML_ARRAY_MDS));
   dserror("Copy of msd matrix not yet implemented ");
break;
case msr:
   matto->msr = (AZ_ARRAY_MSR*)CCACALLOC(1,sizeof(AZ_ARRAY_MSR));
   solserv_cp_msrmask(matfrom->msr,matto->msr);
break;
case parcsr:
   matto->parcsr = (H_PARCSR*)CCACALLOC(1,sizeof(H_PARCSR));
   dserror("Copy of parcsr matrix not yet implemented ");
break;
case ucchb:
   matto->ucchb = (UCCHB*)CCACALLOC(1,sizeof(UCCHB));
   solserv_cp_ucchbmask(matfrom->ucchb,matto->ucchb);
break;
case dense:
   matto->dense = (DENSE*)CCACALLOC(1,sizeof(DENSE));
   solserv_cp_densemask(matfrom->dense,matto->dense);
break;
case rc_ptr:
   matto->rc_ptr = (RC_PTR*)CCACALLOC(1,sizeof(RC_PTR));
   solserv_cp_rc_ptrmask(actintra,matfrom->rc_ptr,matto->rc_ptr);
break;
case ccf:
   matto->ccf = (CCF*)CCACALLOC(1,sizeof(CCF));
   solserv_cp_ccfmask(actintra,matfrom->ccf,matto->ccf);
break;
case skymatrix:
   matto->sky = (SKYMATRIX*)CCACALLOC(1,sizeof(SKYMATRIX));
   solserv_cp_skymask(matfrom->sky,matto->sky);
break;
case spoolmatrix:
   matto->spo = (SPOOLMAT*)CCACALLOC(1,sizeof(SPOOLMAT));
   solserv_cp_spomask(matfrom->spo,matto->spo);
break;
case bdcsr:
   matto->bdcsr = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
   solserv_cp_bdcsrmask(matfrom->bdcsr,matto->bdcsr);
break;
case oll:
   matto->oll = (OLL*)CCACALLOC(1,sizeof(OLL));
   oll_cp_mask(matfrom->oll,matto->oll);
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
 |  copies the sparsity mask of a mlpcg matrix               m.gee 01/03|
 |  for the new sparsity mask, all memory is allocated                  |
 |  called by solserv_alloc_cp_sparsemask only!                         |
 *----------------------------------------------------------------------*/
static void solserv_cp_bdcsrmask(DBCSR *from, DBCSR *to)
{
#ifdef DEBUG 
dstrc_enter("solserv_cp_spomask");
#endif
/*----------------------------------------------------------------------*/
/* copy all information, which is directly included in the structure */
*to = *from;

/* alloccopy update */
am_alloc_copy(&(from->update),&(to->update));
/* alloccopy a */
am_alloc_copy(&(from->a),&(to->a));
/* alloccopy ja */
am_alloc_copy(&(from->ja),&(to->ja));
/* alloccopy ia */
am_alloc_copy(&(from->ia),&(to->ia));
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_cp_spomask */


/*----------------------------------------------------------------------*
 |  copies the sparsity mask of a spoole matrix              m.gee 02/02|
 |  for the new sparsity mask, all memory is allocated                  |
 |  called by solserv_alloc_cp_sparsemask only!                         |
 *----------------------------------------------------------------------*/
static void solserv_cp_spomask(SPOOLMAT *from, SPOOLMAT *to)
{
#ifdef DEBUG 
dstrc_enter("solserv_cp_spomask");
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
/* alloccopy rowptr */
am_alloc_copy(&(from->rowptr),&(to->rowptr));
/* alloccopy A_loc */
am_alloc_copy(&(from->A_loc),&(to->A_loc));
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_cp_spomask */


/*----------------------------------------------------------------------*
 |  copies the sparsity mask of a rc_ptr matrix              m.gee 02/02|
 |  for the new sparsity mask, all memory is allocated                  |
 |  called by solserv_alloc_cp_sparsemask only!                         |
 *----------------------------------------------------------------------*/
static void solserv_cp_rc_ptrmask(INTRA *actintra, RC_PTR *from, RC_PTR *to)
{
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


static void solserv_cp_ccfmask(INTRA *actintra, CCF *from, CCF *to)
{
#ifdef DEBUG 
dstrc_enter("solserv_cp_ccfmask");
#endif
/*----------------------------------------------------------------------*/
/* copy all information, which is directly included in the structure */
*to = *from;

/* alloccopy update */
am_alloc_copy(&(from->update),&(to->update));

/* alloccopy irn_loc */
am_alloc_copy(&(from->Ap),&(to->Ap));

/* alloccopy jcn_loc */
am_alloc_copy(&(from->Ai),&(to->Ai));

/* alloccopy A_loc */
am_alloc_copy(&(from->Ax),&(to->Ax));

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_cp_ccfmask */


/*----------------------------------------------------------------------*
 |  copies the sparsity mask of a ucchb matrix               m.gee 02/02|
 |  for the new sparsity mask, all memory is allocated                  |
 |  called by solserv_alloc_cp_sparsemask only!                         |
 *----------------------------------------------------------------------*/
static void solserv_cp_ucchbmask(UCCHB *from, UCCHB *to)
{
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
static DOUBLE *work1;
static ARRAY   work2_a;
static DOUBLE *work2;

#ifdef DEBUG 
dstrc_enter("solserv_sparsematvec");
#endif
/*----------------------------------------------------------------------*/
if (*mattyp != msr) {
  if (work1_a.Typ != cca_DV) work1 = amdef("work1",&work1_a,vec->numeq_total,1,"DV");
  if (work2_a.Typ != cca_DV) work2 = amdef("work2",&work2_a,vec->numeq_total,1,"DV");
  if (work1_a.fdim < vec->numeq_total) {
    amdel(&work1_a); work1 = amdef("work1",&work1_a,vec->numeq_total,1,"DV");
    amdel(&work2_a); work2 = amdef("work2",&work2_a,vec->numeq_total,1,"DV");
  }
}
/*----------------------------------------------------------------------*/
switch (*mattyp)
{
case mds:
   dserror("Matrix-Vector Product for MLIB not implemented");
break;
#ifdef AZTEC_PACKAGE
case msr:
  /*
   solserv_reddistvec(vec,mat,mattyp,work1,vec->numeq_total,actintra);
   solserv_matvec_msr(actintra,mat->msr,work1,work2);
   solserv_distribdistvec(result,mat,mattyp,work2,vec->numeq_total,actintra);
  */
{
  struct _AZ_ARRAY_MSR* msr_array = mat->msr;

  INT i;
  
  DOUBLE     *tmpsol;
  ARRAY       tmpsol_a;

  DOUBLE     *tmprhs;
  ARRAY       tmprhs_a;

/*----------------------- transform matrix to processor local numbering */
    if (msr_array->is_transformed==0) {
      msr_array->is_transformed = 1;
      
      if (msr_array->Amat != NULL) {
        AZ_matrix_destroy(&(msr_array->Amat)); msr_array->Amat        =NULL;
      }
      if (msr_array->external != NULL) {
        free(msr_array->external);             msr_array->external      =NULL;
      }
      if (msr_array->update_index != NULL) {
        free(msr_array->update_index);         msr_array->update_index  =NULL;
      }
      if (msr_array->extern_index != NULL) {
        free(msr_array->extern_index);         msr_array->extern_index  =NULL;
      }
      if (msr_array->data_org != NULL) {
        free(msr_array->data_org);             msr_array->data_org      =NULL;
      }
      
      /* Make backup copy of bindx, as it is permuted in
       * solution. This has to be done on demand as we need the same
       * thing for solving linear systems.
       *
       * In a sense we abuse bindx_backup because it no longer
       * contains backup data. Instead all communication with aztec
       * relys on bindx_backup (transformed) and the outside world ---
       * that is the assembling --- used the original bindx. */
      if (msr_array->bindx_backup.Typ == cca_XX) {
        am_alloc_copy(&(msr_array->bindx),&(msr_array->bindx_backup));
      }
      else {
        amcopy(&(msr_array->bindx),&(msr_array->bindx_backup));
      }

      AZ_transform(msr_array->proc_config,
                   &(msr_array->external),
                   msr_array->bindx_backup.a.iv,
                   msr_array->val.a.dv,
                   msr_array->update.a.iv,
                   &(msr_array->update_index),
                   &(msr_array->extern_index),
                   &(msr_array->data_org),
                   msr_array->numeq,
                   NULL,
                   NULL,
                   NULL,
                   NULL,
                   AZ_MSR_MATRIX);

      /* create Aztec structure AZ_MATRIX */
      msr_array->Amat = AZ_matrix_create(msr_array->data_org[AZ_N_internal]+
                                         msr_array->data_org[AZ_N_border]);

      /* attach dmsr-matrix to this structure */
      AZ_set_MSR(msr_array->Amat, 
                 msr_array->bindx_backup.a.iv, 
                 msr_array->val.a.dv, 
                 msr_array->data_org, 
                 0, 
                 NULL, 
                 AZ_LOCAL);
    }
    

  /* reorder rhs-vector */
  tmprhs = amdef("tmprhs",&tmprhs_a,msr_array->numeq+msr_array->N_external,1,"DV");
  for (i=0; i<msr_array->numeq; ++i) {
    tmprhs[i] = vec->vec.a.dv[i];
  }
  AZ_reorder_vec(tmprhs,
                 msr_array->data_org,
                 msr_array->update_index,
                 NULL);
  
  /* allocate temporary solution vector */
  tmpsol = amdef("tmpsol",&tmpsol_a,msr_array->numeq,1,"DV");
  amzero(&tmpsol_a);

  /* multiply */
  /*msr_array->Amat->matvec(tmprhs, tmpsol, msr_array->Amat, msr_array->proc_config);*/
  AZ_MSR_matvec_mult(tmprhs, tmpsol, msr_array->Amat, msr_array->proc_config);

  /* invorder solv vector */
  AZ_invorder_vec(tmpsol,
                  msr_array->data_org,
                  msr_array->update_index,
                  NULL,
                  result->vec.a.dv);
  
  /* delete temporary rhs */
  amdel(&tmprhs_a);

  /* delete temporary solution vector */
  amdel(&tmpsol_a);
}
break;
#endif
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
   solserv_matvec_rc_ptr(actintra,mat->rc_ptr,work1,work2);
   solserv_distribdistvec(result,mat,mattyp,work2,vec->numeq_total,actintra);
break;
case spoolmatrix:
   solserv_reddistvec(vec,mat,mattyp,work1,vec->numeq_total,actintra);
   solserv_matvec_spo(actintra,mat->spo,work1,work2);
   solserv_distribdistvec(result,mat,mattyp,work2,vec->numeq_total,actintra);
break;
case bdcsr:
   mlpcg_matvec(result->vec.a.dv,mat->bdcsr,vec->vec.a.dv,1.0,1,actintra);
break;
case ccf:
   solserv_reddistvec(vec,mat,mattyp,work1,vec->numeq_total,actintra);
   solserv_matvec_ccf(actintra,mat->ccf,work1,work2);
   solserv_distribdistvec(result,mat,mattyp,work2,vec->numeq_total,actintra);
break;
case skymatrix:
   solserv_reddistvec(vec,mat,mattyp,work1,vec->numeq_total,actintra);
   solserv_matvec_sky(actintra,mat->sky,work1,work2);
   solserv_distribdistvec(result,mat,mattyp,work2,vec->numeq_total,actintra);
break;
case oll:
   solserv_reddistvec(vec,mat,mattyp,work1,vec->numeq_total,actintra);
   oll_matvec(actintra,mat->oll,work1,work2);
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
 |  make matrix vector multiplication with spooles matrix    m.gee 04/02|
 |  called by solserv_sparsematvec only !                               |
 *----------------------------------------------------------------------*/
static void solserv_matvec_spo(INTRA        *actintra,
                                  SPOOLMAT       *spo,
                                  DOUBLE       *work1,
                                  DOUBLE       *work2)/* work2 is the result */
{
#ifdef SPOOLES_PACKAGE
INT         i,j,dof;
INT         start,end,lenght,j_index;
INT         myrank;
INT         nprocs;
INT         numeq;
INT         numeq_total;
INT        *update;
INT        *row,*irn,*jcn;
DOUBLE     *A;
INT         rowstart,rowend;
INT         col;

#ifdef DEBUG 
dstrc_enter("solserv_matvec_spo");
#endif
/*----------------------------------------------------------------------*/
myrank     = actintra->intra_rank;
nprocs     = actintra->intra_nprocs;
numeq_total= spo->numeq_total;
numeq      = spo->numeq;
update     = spo->update.a.iv;
row        = spo->rowptr.a.iv;
irn        = spo->irn_loc.a.iv;
jcn        = spo->jcn_loc.a.iv;
A          = spo->A_loc.a.dv;
/*------------------------------------------- now loop my own equations */
for (i=0; i<numeq; i++)
{
   dof = update[i];
   work2[dof]=0.0;
   rowstart  = row[i];
   rowend    = row[i+1];
   for (j=rowstart; j<rowend;j++)
   {
      col = jcn[j];
      work2[dof] += A[j] * work1[col];
   }
}
/* every proc added his complete part to the full-sized vector work2.
   There is no need to alreduce this vector here, because it will anyway
   bre distributed to a DIS_VECTOR */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
# endif /* end of #ifdef SPOOLES_PACKAGE */
return;
} /* end of solserv_matvec_spo */


/*----------------------------------------------------------------------*
 |  make matrix vector multiplication with rc_ptr matrix     m.gee 04/02|
 |  called by solserv_sparsematvec only !                               |
 *----------------------------------------------------------------------*/
static void solserv_matvec_rc_ptr(INTRA        *actintra,
                                  RC_PTR       *rcptr,
                                  DOUBLE       *work1,
                                  DOUBLE       *work2)/* work2 is the result */
{
#ifdef MUMPS_PACKAGE
INT         i,j,dof;
INT         start,end,lenght,j_index;
INT         myrank;
INT         nprocs;
INT         numeq;
INT         numeq_total;
INT        *update;
INT        *row,*irn,*jcn;
DOUBLE     *A;
INT         rowstart,rowend;
INT         col;

#ifdef DEBUG 
dstrc_enter("solserv_matvec_rc_ptr");
#endif
/*----------------------------------------------------------------------*/
myrank     = actintra->intra_rank;
nprocs     = actintra->intra_nprocs;
numeq_total= rcptr->numeq_total;
numeq      = rcptr->numeq;
update     = rcptr->update.a.iv;
row        = rcptr->rowptr.a.iv;
irn        = rcptr->irn_loc.a.iv;
jcn        = rcptr->jcn_loc.a.iv;
A          = rcptr->A_loc.a.dv;
/*------------------------------------------- now loop my own equations */
for (i=0; i<numeq; i++)
{
   dof = update[i];
   work2[dof]=0.0;
   rowstart  = row[i];
   rowend    = row[i+1];
   for (j=rowstart; j<rowend;j++)
   {
      col = jcn[j];
      work2[dof] += A[j] * work1[col];
   }
}
/* every proc added his complete part to the full-sized vector work2.
   There is no need to alreduce this vector here, because it will anyway
   bre distributed to a DIS_VECTOR */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
# endif /* end of #ifdef MUMPS_PACKAGE */
return;
} /* end of solserv_matvec_rc_ptr */


/*----------------------------------------------------------------------*
 |  make matrix vector multiplication with ccf matrix  s.offermanns 05/02|
 |  called by solserv_sparsematvec only !                               |
 *----------------------------------------------------------------------*/
static void solserv_matvec_ccf(INTRA        *actintra,
                                  CCF          *ccf,
                                  DOUBLE       *work1,
                                  DOUBLE       *work2)/* work2 is the result */
{
#ifdef UMFPACK
INT         i,j,dof;
INT         myrank;
INT         nprocs;
INT         numeq_total;
INT        *update;
INT        *Ap,*Ai;
DOUBLE     *Ax;
INT         start,end;
INT         col;

#ifdef DEBUG 
dstrc_enter("solserv_matvec_ccf");
#endif
/*----------------------------------------------------------------------*/
myrank      = actintra->intra_rank;
nprocs      = actintra->intra_nprocs;
numeq_total = ccf->numeq_total;
update      = ccf->update.a.iv;
Ap          = ccf->Ap.a.iv;
Ai          = ccf->Ai.a.iv;
Ax          = ccf->Ax.a.dv;
/*------------------------------------------- now loop my own equations */
for (i=0; i<numeq_total; i++)
{
   dof = i;
   work2[dof]=0.0;
   start     = Ap[i];
   end       = Ap[i+1];
   for (j=start; j<end;j++)
   {
      col = Ai[j];
      work2[dof] += Ax[j] * work1[col];
   }
}
/* every proc added his complete part to the full-sized vector work2.
   There is no need to alreduce this vector here, because it will anyway
   bre distributed to a DIS_VECTOR */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif /* end of ifdef UMFPACK */
return;
} /* end of solserv_matvec_ccf */



#if 0
/* No longer needed. We use aztec's function now. */
/*----------------------------------------------------------------------*
 |  make matrix vector multiplication with msr matrix        m.gee 02/02|
 |  called by solserv_sparsematvec only !                               |
 *----------------------------------------------------------------------*/
static void solserv_matvec_msr(
    INTRA        *actintra,
    AZ_ARRAY_MSR *msr,
    DOUBLE       *work1,
    DOUBLE       *work2)
{
#ifdef AZTEC_PACKAGE
INT         i,j,dof;
INT         start,end,lenght,j_index;
INT         myrank;
INT         nprocs;
INT         numeq;
INT         numeq_total;
INT        *update;
INT        *bindx;
DOUBLE     *val;

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
#endif





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
                        DOUBLE       *work1,
                        DOUBLE       *work2)
{
INT     i,j,kl,ku,kdiff;
INT     nrn,neq1,neq;
DOUBLE *A;
INT    *maxa;

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
                          DOUBLE       *work1,
                          DOUBLE       *work2)
{
INT      i,j;
INT      I,J;
DOUBLE   sum;
DOUBLE **A;
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


/*----------------------------------------------------------------------*
 |  make matrix vector multiplication with oll     matrix      m.n 04/03|
 |  called by solserv_sparsematvec only !                               |
 *----------------------------------------------------------------------*/
static void oll_matvec(
    INTRA        *actintra,
    OLL          *oll,
    DOUBLE       *work1,
    DOUBLE       *work2)        /* work2 is the result */
{
  INT         i,dof;
  INT         numeq;
  INT         numeq_total;
  INT        *update;
  MATENTRY  **row;
  MATENTRY   *actentry;

#ifdef DEBUG 
  dstrc_enter("oll_matvec");
#endif
  /*----------------------------------------------------------------------*/
  numeq_total= oll->numeq_total;
  numeq      = oll->numeq;
  update     = oll->update.a.iv;
  row = oll->row;

  /*------------------------------------------- now loop my own equations */
  for (i=0; i<oll->rdim; i++)
  {
    dof = update[i];
    work2[dof]=0.0;
    actentry = row[i];
    while(actentry != NULL )
    {
      work2[dof] += actentry->val * work1[actentry->c];
      actentry = actentry->rnext;
    } /* end while row i */
  } /* end for all rows */

  /* every proc added his complete part to the full-sized vector work2.
     There is no need to alreduce this vector here, because it will anyway
     be distributed to a DIS_VECTOR */
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
} /* end of oll_matvec */




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
#ifndef CCADISCRET
#include "../headers/standardtypes.h"
#include "../solver/solver.h"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | struct _GENPROB       genprob; defined in global_control.c           |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 |                                                           m.gee 03/02|
 | prototypes of functions only visible in this file                    |
 *----------------------------------------------------------------------*/
#ifdef MUMPS_PACKAGE
static void solserv_cp_rc_ptrmask(INTRA *actintra, RC_PTR *from, RC_PTR *to);
static void solserv_matvec_rc_ptr(INTRA  *actintra,RC_PTR *rcptr,DOUBLE *work1,DOUBLE *work2);
#endif

#ifdef UMFPACK
static void solserv_cp_ccfmask(INTRA *actintra, CCF *from, CCF *to);
static void solserv_matvec_ccf(INTRA  *actintra,CCF *ccf,DOUBLE *work1,DOUBLE *work2);
#endif

#ifdef PARSUPERLU_PACKAGE
static void solserv_cp_ucchbmask(UCCHB *from, UCCHB *to);
#endif

static void solserv_cp_skymask(SKYMATRIX *from, SKYMATRIX *to);
static void solserv_matvec_sky(INTRA *actintra,SKYMATRIX *sky,DOUBLE *work1,DOUBLE *work2);

#ifdef AZTEC_PACKAGE
static void solserv_cp_msrmask(AZ_ARRAY_MSR *from, AZ_ARRAY_MSR *to);
static void solserv_matvec_msr(INTRA *actintra,AZ_ARRAY_MSR *msr,DOUBLE *work1,DOUBLE *work2);
#endif

static void solserv_cp_densemask(DENSE *from, DENSE *to);
static void solserv_matvec_dense(INTRA *actintra,DENSE *dense,DOUBLE *work1,DOUBLE *work2);

#ifdef SPOOLES_PACKAGE
static void solserv_cp_spomask(SPOOLMAT *from, SPOOLMAT *to);
static void solserv_matvec_spo(INTRA *actintra,SPOOLMAT *spo,DOUBLE *work1,DOUBLE *work2);
#endif

#ifdef MLPCG
static void solserv_cp_bdcsrmask(DBCSR *from, DBCSR *to);
#endif

static void oll_matvec(INTRA *actintra, OLL *oll, DOUBLE *work1, DOUBLE *work2);

/*----------------------------------------------------------------------*
 |                                                           m.gee 09/06|
 | prototypes of functions from solver_trilinos_service.cpp             |
 *----------------------------------------------------------------------*/
#ifdef TRILINOS_PACKAGE
extern void trilinos_cp_matrixmask(TRILINOSMATRIX  *from, TRILINOSMATRIX* to);
extern void trilinos_zero_matrix(TRILINOSMATRIX *tri);
extern void add_trilinos_matrix(TRILINOSMATRIX* from, TRILINOSMATRIX* to, double factor);
extern void close_trilinos_matrix(struct _TRILINOSMATRIX *tri);
extern void matvec_trilinos(DIST_VECTOR* y, DIST_VECTOR* x, TRILINOSMATRIX* A);
extern void matvec_trilinos_trans(DIST_VECTOR* y, DIST_VECTOR* x, TRILINOSMATRIX* A);
extern void scale_trilinos_matrix(TRILINOSMATRIX* A, double factor);
#endif

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

#ifdef MLIB_PACKAGE
case mds:
   dserror("not implemented for MLIB yet");
break;
#endif

#ifdef TRILINOS_PACKAGE
case trilinos:
   scale_trilinos_matrix(A->trilinos,factor);
break;
#endif

#ifdef AZTEC_PACKAGE
case msr:
   amscal(&(A->msr->val),&factor);
break;
#endif

#ifdef HYPRE_PARCSR
case parcsr:
   dserror("not implemented for HYPRE yet");
break;
#endif

#ifdef PARSUPERLU_PACKAGE
case ucchb:
   amscal(&(A->ucchb->a),&factor);
break;
#endif

case dense:
   amscal(&(A->dense->A),&factor);
break;

#ifdef MUMPS_PACKAGE
case rc_ptr:
   amscal(&(A->rc_ptr->A_loc),&factor);
break;
#endif

#ifdef SPOOLES_PACKAGE
case spoolmatrix:
   amscal(&(A->spo->A_loc),&factor);
break;
#endif

#ifdef MLPCG
case bdcsr:
   amscal(&(A->bdcsr->a),&factor);
break;
#endif

case ccf:
#ifdef UMFPACK
   amscal(&(A->ccf->Ax),&factor);
break;
#endif

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
 |                                                           m.gee 09/06|
 | close a system matrix                                                |
 | this might be necessary for certain dynamic system matrices          |
 | such as spooles, oll, trilinos                                       |
 *----------------------------------------------------------------------*/
void solserv_close_mat(INTRA *actintra, SPARSE_TYP* Atyp,SPARSE_ARRAY* A)
{

#ifdef DEBUG
dstrc_enter("solserv_close_mat");
#endif
/*----------------------------------------------------------------------*/
switch (*Atyp)
{

#ifdef MLIB_PACKAGE
case mds:
break;
#endif

#ifdef TRILINOS_PACKAGE
case trilinos:
  close_trilinos_matrix(A->trilinos);
break;
#endif

#ifdef AZTEC_PACKAGE
case msr:
break;
#endif

#ifdef HYPRE_PARCSR
case parcsr:
break;
#endif

#ifdef PARSUPERLU_PACKAGE
case ucchb:
break;
#endif

case dense:
break;

#ifdef MUMPS_PACKAGE
case rc_ptr:
break;
#endif

#ifdef SPOOLES_PACKAGE
case spoolmatrix:
   close_spooles_matrix(A->spo,actintra);
break;
#endif

#ifdef MLPCG
case bdcsr:
break;
#endif

case ccf:
#ifdef UMFPACK
break;
#endif

case skymatrix:
break;

case oll:
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
} /* end of solserv_close_mat */


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

#ifdef MLIB_PACKAGE
case mds:
   dserror("not implemented for MLIB yet");
break;
#endif

#ifdef TRILINOS_PACKAGE
case trilinos:
  add_trilinos_matrix(B->trilinos,A->trilinos,factor);
break;
#endif

#ifdef AZTEC_PACKAGE
case msr:
   amadd(&(A->msr->val),&(B->msr->val),factor,0);
break;
#endif

#ifdef HYPRE_PACKAGE
case parcsr:
   dserror("not implemented for HYPRE yet");
break;
#endif

#ifdef PARSUPERLU_PACKAGE
case ucchb:
   amadd(&(A->ucchb->a),&(B->ucchb->a),factor,0);
break;
#endif

case dense:
   amadd(&(A->dense->A),&(B->dense->A),factor,0);
break;

#ifdef MUMPS_PACKAGE
case rc_ptr:
   amadd(&(A->rc_ptr->A_loc),&(B->rc_ptr->A_loc),factor,0);
break;
#endif

#ifdef SPOOLES_PACKAGE
case spoolmatrix:
#ifdef D_CONTACT
   add_spooles_matrix(A->spo,B->spo,factor,0,actintra);
#else
   amadd(&(A->spo->A_loc),&(B->spo->A_loc),factor,0);
#endif
break;
#endif

#ifdef MLPCG
case bdcsr:
   amadd(&(A->bdcsr->a),&(B->bdcsr->a),factor,0);
break;
#endif

#ifdef UMFPACK
case ccf:
   amadd(&(A->ccf->Ax),&(B->ccf->Ax),factor,0);
break;
#endif

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

#ifdef MLIB
case mds:
   *numeq       = mat->mds->numeq;
   *numeq_total = *numeq;
break;
#endif

#ifdef AZTEC_PACKAGE
case msr:
   *numeq       = mat->msr->numeq;
   *numeq_total = mat->msr->numeq_total;
break;
#endif

#ifdef HYPRE_PACKAGE
case parcsr:
   *numeq       = mat->parcsr->numeq;
   *numeq_total = mat->parcsr->numeq_total;
break;
#endif

#ifdef PARSUPERLU_PACKAGE
case ucchb:
   *numeq       = mat->ucchb->numeq;
   *numeq_total = mat->ucchb->numeq_total;
break;
#endif

case dense:
   *numeq       = mat->dense->numeq;
   *numeq_total = mat->dense->numeq_total;
break;

#ifdef MUMPS_PACKAGE
case rc_ptr:
   *numeq       = mat->rc_ptr->numeq;
   *numeq_total = mat->rc_ptr->numeq_total;
break;
#endif

#ifdef UMFPACK
case ccf:
   *numeq       = mat->ccf->numeq;
   *numeq_total = mat->ccf->numeq_total;
break;
#endif

case skymatrix:
   *numeq       = mat->sky->numeq;
   *numeq_total = mat->sky->numeq_total;
break;

#ifdef SPOOLES_PACKAGE
case spoolmatrix:
   *numeq       = mat->spo->numeq;
   *numeq_total = mat->spo->numeq_total;
break;
#endif

#ifdef MLPCG
case bdcsr:
   *numeq       = mat->bdcsr->numeq;
   *numeq_total = mat->bdcsr->numeq_total;
break;
#endif

case oll:
   *numeq       = mat->oll->numeq;
   *numeq_total = mat->oll->numeq_total;
break;

#ifdef TRILINOS_PACKAGE
case trilinos:
   *numeq       = mat->trilinos->numeq;
   *numeq_total = mat->trilinos->numeq_total;
break;
#endif

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
#ifdef MLIB_PACKAGE
case mds:
   solver_mlib(NULL,NULL,mat->mds,NULL,NULL,2);
break;
#endif

#ifdef AZTEC_PACKAGE
case msr:
   amzero(&(mat->msr->val));
   mat->msr->is_factored=0;
   mat->msr->is_transformed=0;
break;
#endif

#ifdef HYPRE_PACKAGE
case parcsr:/*---- this stupid package does not have a zero function!!!!*/
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
   mat->parcsr->is_factored=0;
break;
#endif

#ifdef PARSUPERLU_PACKAGE
case ucchb:
   amzero(&(mat->ucchb->a));
   mat->ucchb->is_factored=0;
break;
#endif

case dense:
   amzero(&(mat->dense->A));
   mat->dense->is_factored=0;
break;

#ifdef MUMPS_PACKAGE
case rc_ptr:
   amzero(&(mat->rc_ptr->A_loc));
   mat->rc_ptr->is_factored=0;
break;
#endif

#ifdef UMFPACK
case ccf:
   amzero(&(mat->ccf->Ax));
  /* mat->ccf->is_factored=0; is never used but for struct2_ml,
                              and there I need it not to be set to zero */
break;
#endif

case skymatrix:
   amzero(&(mat->sky->A));
   mat->sky->is_factored=0;
break;

#ifdef MLPCG
case bdcsr:
   amzero(&(mat->bdcsr->a));
   mat->bdcsr->is_factored=0;
break;
#endif

#ifdef SPOOLES_PACKAGE
case spoolmatrix:
   amzero(&(mat->spo->A_loc));
#ifdef D_CONTACT
   aminit(&(mat->spo->irn_loc),&minusone);
   aminit(&(mat->spo->jcn_loc),&minusone);
#endif
   mat->spo->is_factored=0;
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
break;
#endif

case oll:
   oll_zero(mat->oll);
   /*solserv_zero_mat(actintra,&(mat->oll->sysarray[0]),&(mat->oll->sysarray_typ[0]));*/
break;

#ifdef TRILINOS_PACKAGE
case trilinos:
   trilinos_zero_matrix(mat->trilinos);
break;
#endif

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

#ifdef MLIB_PACKAGE
case mds:
   matto->mds = (ML_ARRAY_MDS*)CCACALLOC(1,sizeof(ML_ARRAY_MDS));
   dserror("Copy of msd matrix not yet implemented ");
break;
#endif

#ifdef AZTEC_PACKAGE
case msr:
   matto->msr = (AZ_ARRAY_MSR*)CCACALLOC(1,sizeof(AZ_ARRAY_MSR));
   solserv_cp_msrmask(matfrom->msr,matto->msr);
break;
#endif

#ifdef HYPRE_PACKAGE
case parcsr:
   matto->parcsr = (H_PARCSR*)CCACALLOC(1,sizeof(H_PARCSR));
   dserror("Copy of parcsr matrix not yet implemented ");
break;
#endif

#ifdef PARSUPERLU_PACKAGE
case ucchb:
   matto->ucchb = (UCCHB*)CCACALLOC(1,sizeof(UCCHB));
   solserv_cp_ucchbmask(matfrom->ucchb,matto->ucchb);
break;
#endif

case dense:
   matto->dense = (DENSE*)CCACALLOC(1,sizeof(DENSE));
   solserv_cp_densemask(matfrom->dense,matto->dense);
break;

#ifdef MUMPS_PACKAGE
case rc_ptr:
   matto->rc_ptr = (RC_PTR*)CCACALLOC(1,sizeof(RC_PTR));
   solserv_cp_rc_ptrmask(actintra,matfrom->rc_ptr,matto->rc_ptr);
break;
#endif

#ifdef UMFPACK
case ccf:
   matto->ccf = (CCF*)CCACALLOC(1,sizeof(CCF));
   solserv_cp_ccfmask(actintra,matfrom->ccf,matto->ccf);
break;
#endif

case skymatrix:
   matto->sky = (SKYMATRIX*)CCACALLOC(1,sizeof(SKYMATRIX));
   solserv_cp_skymask(matfrom->sky,matto->sky);
break;

#ifdef SPOOLES_PACKAGE
case spoolmatrix:
   matto->spo = (SPOOLMAT*)CCACALLOC(1,sizeof(SPOOLMAT));
   solserv_cp_spomask(matfrom->spo,matto->spo);
break;
#endif

#ifdef MLPCG
case bdcsr:
   matto->bdcsr = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
   solserv_cp_bdcsrmask(matfrom->bdcsr,matto->bdcsr);
break;
#endif

case oll:
   matto->oll = (OLL*)CCACALLOC(1,sizeof(OLL));
   oll_cp_mask(matfrom->oll,matto->oll);
break;

#ifdef TRILINOS_PACKAGE
case trilinos:
   matto->trilinos = (TRILINOSMATRIX*)CCACALLOC(1,sizeof(TRILINOSMATRIX));
   trilinos_cp_matrixmask(matfrom->trilinos,matto->trilinos);
break;
#endif

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


#ifdef MLPCG
/*----------------------------------------------------------------------*
 |  copies the sparsity mask of a mlpcg matrix               m.gee 01/03|
 |  for the new sparsity mask, all memory is allocated                  |
 |  called by solserv_alloc_cp_sparsemask only!                         |
 *----------------------------------------------------------------------*/
static void solserv_cp_bdcsrmask(DBCSR *from, DBCSR *to)
{
#ifdef DEBUG
dstrc_enter("solserv_cp_bdcsrmask");
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
} /* end of solserv_cp_bdcsrmask */
#endif


#ifdef SPOOLES_PACKAGE
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
#endif


#ifdef MUMPS_PACKAGE
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
#endif



#ifdef UMFPACK
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
#endif




#ifdef PARSUPERLU_PACKAGE
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
#endif




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



#ifdef AZTEC_PACKAGE
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
#endif



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
#ifdef TRILINOS_PACKAGE
if (genprob.usetrilinosalgebra)
if (*mattyp==trilinos)
{
  matvec_trilinos(result,vec,mat->trilinos);
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}
#endif
/*----------------------------------------------------------------------*/
/* these are redundant working vectors needed by some matvec versions */
if (work1_a.Typ != cca_DV) work1 = amdef("work1",&work1_a,vec->numeq_total,1,"DV");
if (work2_a.Typ != cca_DV) work2 = amdef("work2",&work2_a,vec->numeq_total,1,"DV");
if (work1_a.fdim < vec->numeq_total)
{
  amdel(&work1_a); work1 = amdef("work1",&work1_a,vec->numeq_total,1,"DV");
  amdel(&work2_a); work2 = amdef("work2",&work2_a,vec->numeq_total,1,"DV");
}

/*----------------------------------------------------------------------*/
switch (*mattyp)
{
#ifdef MLIB_PACKAGE
case mds:
   dserror("Matrix-Vector Product for MLIB not implemented");
break;
#endif

#ifdef TRILINOS_PACKAGE
case trilinos:
  dserror("Put 'ALGEBRA Trilinos' in the '---PROBLEM TYP' block of your input file to use trilinos algebra");
break;
#endif

#ifdef AZTEC_PACKAGE
case msr:
   solserv_reddistvec(vec,mat,mattyp,work1,vec->numeq_total,actintra);
   solserv_matvec_msr(actintra,mat->msr,work1,work2);
   solserv_distribdistvec(result,mat,mattyp,work2,vec->numeq_total,actintra);
break;
#endif

#ifdef HYPRE_PACKAGE
case parcsr:
   solserv_reddistvec(vec,mat,mattyp,work1,vec->numeq_total,actintra);
   solserv_distribdistvec(result,mat,mattyp,work2,vec->numeq_total,actintra);
   dserror("Matrix-Vector Product for HYPRE not implemented");
break;
#endif

#ifdef PARSUPERLU_PACKAGE
case ucchb:
   solserv_reddistvec(vec,mat,mattyp,work1,vec->numeq_total,actintra);
   solserv_distribdistvec(result,mat,mattyp,work2,vec->numeq_total,actintra);
   dserror("Matrix-Vector Product for SuperLU not implemented");
break;
#endif

case dense:
   solserv_reddistvec(vec,mat,mattyp,work1,vec->numeq_total,actintra);
   solserv_matvec_dense(actintra,mat->dense,work1,work2);
   solserv_distribdistvec(result,mat,mattyp,work2,vec->numeq_total,actintra);
break;

#ifdef MUMPS_PACKAGE
case rc_ptr:
   solserv_reddistvec(vec,mat,mattyp,work1,vec->numeq_total,actintra);
   solserv_matvec_rc_ptr(actintra,mat->rc_ptr,work1,work2);
   solserv_distribdistvec(result,mat,mattyp,work2,vec->numeq_total,actintra);
break;
#endif

#ifdef SPOOLES_PACKAGE
case spoolmatrix:
   solserv_reddistvec(vec,mat,mattyp,work1,vec->numeq_total,actintra);
   solserv_matvec_spo(actintra,mat->spo,work1,work2);
   solserv_distribdistvec(result,mat,mattyp,work2,vec->numeq_total,actintra);
break;
#endif

#ifdef MLPCG
case bdcsr:
   mlpcg_matvec(result->vec.a.dv,mat->bdcsr,vec->vec.a.dv,1.0,1,actintra);
break;
#endif

#ifdef UMFPACK
case ccf:
   solserv_reddistvec(vec,mat,mattyp,work1,vec->numeq_total,actintra);
   solserv_matvec_ccf(actintra,mat->ccf,work1,work2);
   solserv_distribdistvec(result,mat,mattyp,work2,vec->numeq_total,actintra);
break;
#endif

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



#ifdef SPOOLES_PACKAGE
/*----------------------------------------------------------------------*
 |  make matrix vector multiplication with spooles matrix    m.gee 04/02|
 |  called by solserv_sparsematvec only !                               |
 *----------------------------------------------------------------------*/
static void solserv_matvec_spo(INTRA        *actintra,
                                  SPOOLMAT       *spo,
                                  DOUBLE       *work1,
                                  DOUBLE       *work2)/* work2 is the result */
{
INT         i,j,dof;
/*INT         start,end,lenght,j_index;*/
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
return;
} /* end of solserv_matvec_spo */
# endif /* ifdef SPOOLES_PACKAGE */




#ifdef MUMPS_PACKAGE
/*----------------------------------------------------------------------*
 |  make matrix vector multiplication with rc_ptr matrix     m.gee 04/02|
 |  called by solserv_sparsematvec only !                               |
 *----------------------------------------------------------------------*/
static void solserv_matvec_rc_ptr(INTRA        *actintra,
                                  RC_PTR       *rcptr,
                                  DOUBLE       *work1,
                                  DOUBLE       *work2)/* work2 is the result */
{
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
return;
} /* end of solserv_matvec_rc_ptr */
# endif /* ifdef MUMPS_PACKAGE */





#ifdef UMFPACK
/*----------------------------------------------------------------------*
 |  make matrix vector multiplication with ccf matrix  s.offermanns 05/02|
 |  called by solserv_sparsematvec only !                               |
 *----------------------------------------------------------------------*/
static void solserv_matvec_ccf(INTRA        *actintra,
                                  CCF          *ccf,
                                  DOUBLE       *work1,
                                  DOUBLE       *work2)/* work2 is the result */
{
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
return;
} /* end of solserv_matvec_ccf */
#endif /* ifdef UMFPACK */



#ifdef AZTEC_PACKAGE
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
return;
} /* end of solserv_matvec_msr */
# endif /* ifdef AZTEC_PACKAGE */





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



#endif

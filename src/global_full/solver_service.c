#include "../headers/standardtypes.h"
#include "../headers/solution.h"

/*----------------------------------------------------------------------*
 |  get dimensions of sparse matrix                          m.gee 02/02|
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
 *----------------------------------------------------------------------*/
int solserv_cp_rc_ptrmask(INTRA *actintra, RC_PTR *from, RC_PTR *to)
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
 *----------------------------------------------------------------------*/
int solserv_cp_ucchbmask(UCCHB *from, UCCHB *to)
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
 *----------------------------------------------------------------------*/
int solserv_cp_skymask(SKYMATRIX *from, SKYMATRIX *to)
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
 *----------------------------------------------------------------------*/
int solserv_cp_msrmask(AZ_ARRAY_MSR *from, AZ_ARRAY_MSR *to)
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
 *----------------------------------------------------------------------*/
int solserv_cp_densemask(DENSE *from, DENSE *to)
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

